from __future__ import annotations

from collections import OrderedDict, defaultdict
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import TypedDict

from .gff import child_ids, parse_attributes
from .gff_feature_sync import (
    synchronize_transcript_children,
    transcript_boundary_rows,
    validate_parent_child_containment,
)
from .step_logging import write_id_list, write_step_log, write_step_metrics
from .transcript_models import group_direct_cds
from .utils import MSSPackError, atomic_text_writer


@dataclass
class _MrnaModel:
    line: list[str]
    exons: list[list[str]] = field(default_factory=list)
    cdss: list[list[str]] = field(default_factory=list)
    children: list[list[str]] = field(default_factory=list)


@dataclass
class _GeneModel:
    line: list[str]
    mrnas: OrderedDict[str, _MrnaModel] = field(default_factory=OrderedDict)


class InframeFixSummary(TypedDict):
    updated_gene_models: int
    unchanged_gene_models: int
    removed_features: int
    synchronized_features: int
    updated_gene_ids: list[str]


def _safe_phase(value: str) -> int:
    return int(value) if value.isdigit() else 0


def _compute_mrna_boundaries(features: list[list[str]]) -> tuple[int | None, int | None]:
    if not features:
        return None, None
    starts = [int(feature[3]) for feature in features]
    ends = [int(feature[4]) for feature in features]
    return min(starts), max(ends)


def _compute_total_cds_length(cdss: list[list[str]]) -> int:
    return sum(max(0, int(cds[4]) - int(cds[3]) + 1) for cds in cdss)


def _find_terminal_cds(cdss: list[list[str]], strand: str, *, first: bool) -> list[str] | None:
    if not cdss:
        return None
    if strand == "+":
        chooser = min if first else max
    else:
        chooser = max if first else min
    return chooser(cdss, key=lambda fields: int(fields[3]))


def _update_matching_exons(
    exons: list[list[str]],
    *,
    old_start: int,
    old_end: int,
    new_start: int,
    new_end: int,
    strand: str,
) -> None:
    for exon_line in exons:
        if (
            int(exon_line[3]) == old_start
            and int(exon_line[4]) == old_end
            and exon_line[6] == strand
        ):
            exon_line[3] = str(new_start)
            exon_line[4] = str(new_end)


def _adjust_first_cds_frame(
    cdss: list[list[str]],
    strand: str,
    exons: list[list[str]],
) -> int:
    first_cds = _find_terminal_cds(cdss, strand, first=True)
    if first_cds is None:
        return 0
    old_start = int(first_cds[3])
    old_end = int(first_cds[4])
    old_phase = _safe_phase(first_cds[7])
    if old_phase == 0:
        return 0

    if strand == "+":
        new_start = old_start + old_phase
        first_cds[3] = str(new_start)
        _update_matching_exons(
            exons,
            old_start=old_start,
            old_end=old_end,
            new_start=new_start,
            new_end=old_end,
            strand=strand,
        )
    else:
        new_end = old_end - old_phase
        first_cds[4] = str(new_end)
        _update_matching_exons(
            exons,
            old_start=old_start,
            old_end=old_end,
            new_start=old_start,
            new_end=new_end,
            strand=strand,
        )
    first_cds[7] = "0"
    return old_phase


def _truncate_last_cds_to_multiple_of_three(
    cdss: list[list[str]],
    strand: str,
    exons: list[list[str]],
    removed_row_ids: set[int],
) -> int:
    remainder = _compute_total_cds_length(cdss) % 3
    remaining = remainder
    while remaining and cdss:
        last_cds = _find_terminal_cds(cdss, strand, first=False)
        assert last_cds is not None
        old_start, old_end = int(last_cds[3]), int(last_cds[4])
        length = old_end - old_start + 1
        if length <= remaining:
            # Remove matching exons before mutating the CDS coordinates.
            _remove_cds_and_matching_exons(cdss, exons, last_cds, removed_row_ids)
            remaining -= length
            continue
        new_start = old_start + remaining if strand == "-" else old_start
        new_end = old_end - remaining if strand == "+" else old_end
        last_cds[3], last_cds[4] = str(new_start), str(new_end)
        _update_matching_exons(
            exons, old_start=old_start, old_end=old_end,
            new_start=new_start, new_end=new_end, strand=strand,
        )
        remaining = 0
    return remainder


def _remove_cds_and_matching_exons(
    cdss: list[list[str]],
    exons: list[list[str]],
    cds_line: list[str],
    removed_row_ids: set[int],
) -> None:
    if cds_line in cdss:
        cdss.remove(cds_line)
    removed_row_ids.add(id(cds_line))

    cds_start = int(cds_line[3])
    cds_end = int(cds_line[4])
    cds_strand = cds_line[6]
    matching_exons = [
        exon
        for exon in exons
        if int(exon[3]) == cds_start
        and int(exon[4]) == cds_end
        and exon[6] == cds_strand
    ]
    for exon in matching_exons:
        exons.remove(exon)
        removed_row_ids.add(id(exon))


def fix_gff_to_inframe(
    *,
    input_path: str | Path,
    output_path: str | Path,
    log_path: str | Path,
    updated_gene_ids_path: str | Path | None = None,
    metrics_path: str | Path | None = None,
) -> InframeFixSummary:
    started_at = datetime.now()
    all_lines: list[str | list[str]] = []
    feature_lines: list[list[str]] = []
    children_of: dict[str, list[list[str]]] = defaultdict(list)
    typed_lines: list[tuple[list[str], str, str]] = []

    with Path(input_path).open("r", encoding="utf-8") as handle:
        in_fasta = False
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if in_fasta:
                all_lines.append(line)
                continue
            if line == "##FASTA":
                in_fasta = True
                all_lines.append(line)
                continue
            if not line or line.startswith("#"):
                all_lines.append(line)
                continue
            cols = line.split("\t")
            if len(cols) < 9:
                all_lines.append(line)
                continue
            all_lines.append(cols)
            feature_lines.append(cols)

    for cols in feature_lines:
        attrs = parse_attributes(cols[8])
        feature_id = attrs.get("ID", "")
        typed_lines.append((cols, cols[2], feature_id))
        for parent_id in child_ids(attrs.get("Parent")):
            children_of[parent_id].append(cols)

    gene_lines_by_id: dict[str, list[str]] = {}
    gene_ids_in_order: list[str] = []
    for cols, feature_type, feature_id in typed_lines:
        if feature_type == "gene" and feature_id:
            gene_ids_in_order.append(feature_id)
            gene_lines_by_id[feature_id] = cols

    genes: OrderedDict[str, _GeneModel] = OrderedDict()
    for gene_id in gene_ids_in_order:
        genes[gene_id] = _GeneModel(line=gene_lines_by_id[gene_id])

    for gene_id, gene_data in genes.items():
        mrna_map = gene_data.mrnas
        for mline in (
            child
            for child in children_of.get(gene_id, [])
            if child[2] in ("mRNA", "transcript")
        ):
            mrna_id = parse_attributes(mline[8]).get("ID")
            if mrna_id:
                mrna_map[mrna_id] = _MrnaModel(line=mline)

    # Root transcripts and gene-direct CDS are valid models too. Virtual containers
    # organize them for adjustment; they are never emitted as extra GFF records.
    owned = {transcript_id for gene in genes.values() for transcript_id in gene.mrnas}
    for row, feature_type, feature_id in typed_lines:
        if feature_type in ("mRNA", "transcript") and feature_id and feature_id not in owned:
            parent_ids = child_ids(parse_attributes(row[8]).get("Parent")) or [feature_id]
            for gene_id in parent_ids:
                gene = genes.setdefault(gene_id, _GeneModel(line=row.copy()))
                gene.mrnas[feature_id] = _MrnaModel(line=row)
    for gene_id, gene in genes.items():
        direct_children = children_of.get(gene_id, [])
        groups = group_direct_cds(gene_id, (
            (parse_attributes(child[8]).get("ID", ""), child)
            for child in direct_children if child[2] == "CDS"
        ))
        for model_id, cdss in groups.items():
            children = (
                [child for child in direct_children if child[2] not in {"mRNA", "transcript"}]
                if len(groups) == 1 else [
                    child for child in direct_children
                    if any(child is cds for cds in cdss)
                    or (child[2] == "exon" and any(child[3:5] == cds[3:5] for cds in cdss))
                ]
            )
            line = gene.line.copy()
            line[3] = str(min(int(child[3]) for child in children))
            line[4] = str(max(int(child[4]) for child in children))
            gene.mrnas[model_id] = _MrnaModel(
                line=line, cdss=cdss, children=children,
                exons=[child for child in children if child[2] == "exon"],
            )

    for gene_data in genes.values():
        mrna_map = gene_data.mrnas
        for mrna_id, mrna_data in mrna_map.items():
            if mrna_data.children:
                continue
            for child in children_of.get(mrna_id, []):
                mrna_data.children.append(child)
                if child[2] == "exon":
                    mrna_data.exons.append(child)
                elif child[2] == "CDS":
                    mrna_data.cdss.append(child)

    removed_row_ids: set[int] = set()
    num_updated = 0
    num_unchanged = 0
    synchronized_features = 0
    updated_gene_ids: list[str] = []
    adjusted_parent_ids: set[str] = set()

    for gene_id, gene_data in genes.items():
        gene_changed = False
        mrna_map = gene_data.mrnas
        for mrna_data in mrna_map.values():
            transcript_changed = False
            exons = mrna_data.exons
            cdss = mrna_data.cdss
            mrna_line = mrna_data.line
            if not cdss:
                continue
            if _compute_total_cds_length(cdss) < 3:
                continue
            strand = mrna_line[6]

            while cdss:
                adjusted_bases = _adjust_first_cds_frame(cdss, strand, exons)
                if not adjusted_bases:
                    break
                transcript_changed = True
                invalid_first_cdss = [cds for cds in cdss if int(cds[3]) > int(cds[4])]
                if not invalid_first_cdss:
                    break
                for invalid_cds in invalid_first_cdss:
                    _remove_cds_and_matching_exons(
                        cdss,
                        exons,
                        invalid_cds,
                        removed_row_ids,
                    )

            remainder = _truncate_last_cds_to_multiple_of_three(
                cdss, strand, exons, removed_row_ids,
            )
            if remainder:
                transcript_changed = True

            new_start, new_end = _compute_mrna_boundaries(transcript_boundary_rows(
                mrna_data.children, removed_row_ids=removed_row_ids,
            ))
            if new_start is None or new_end is None:
                continue
            if new_start != int(mrna_line[3]) or new_end != int(mrna_line[4]):
                transcript_changed = True
            mrna_line[3] = str(new_start)
            mrna_line[4] = str(new_end)
            if transcript_changed:
                transcript_id = parse_attributes(mrna_line[8]).get("ID", "")
                if transcript_id:
                    adjusted_parent_ids.add(transcript_id)
                synchronized_features += synchronize_transcript_children(
                    transcript_row=mrna_line,
                    child_rows=mrna_data.children,
                    removed_row_ids=removed_row_ids,
                )
                gene_changed = True

        if gene_changed:
            transcript_starts = [int(model.line[3]) for model in mrna_map.values()]
            transcript_ends = [int(model.line[4]) for model in mrna_map.values()]
            direct_children = [
                row for row in children_of.get(gene_id, []) if id(row) not in removed_row_ids
            ]
            transcript_starts.extend(int(row[3]) for row in direct_children)
            transcript_ends.extend(int(row[4]) for row in direct_children)
            if transcript_starts and transcript_ends:
                gene_line = gene_data.line
                gene_line[3] = str(min(transcript_starts))
                gene_line[4] = str(max(transcript_ends))
                adjusted_parent_ids.add(gene_id)
            num_updated += 1
            updated_gene_ids.append(gene_id)
        else:
            num_unchanged += 1

    hierarchy_issues = validate_parent_child_containment(
        feature_lines,
        scope_parent_ids=adjusted_parent_ids,
        removed_row_ids=removed_row_ids,
    )
    if hierarchy_issues:
        issue_text = "; ".join(issue.message for issue in hierarchy_issues[:5])
        raise MSSPackError(f"Coordinate adjustment produced an invalid GFF hierarchy: {issue_text}")

    with atomic_text_writer(Path(output_path)) as handle:
        for item in all_lines:
            if isinstance(item, str):
                handle.write(item + "\n")
                continue
            if id(item) in removed_row_ids:
                continue
            handle.write("\t".join(item) + "\n")

    write_step_log(
        log_path=Path(log_path),
        command=f"msspack internal update-gff-to-inframe --input {input_path} --output {output_path}",
        step="update-gff-to-inframe",
        started_at=started_at,
        count_unit="genes",
        input_total=len(genes),
        changed_total=num_updated,
        output_total=len(genes),
        details=[
            f"Number of unchanged gene models: {num_unchanged:,}",
            f"Removed features: {len(removed_row_ids):,}",
            f"Synchronized dependent features: {synchronized_features:,}",
        ],
    )
    if updated_gene_ids_path is not None:
        write_id_list(Path(updated_gene_ids_path), updated_gene_ids)
    if metrics_path is not None:
        write_step_metrics(
            metrics_path=Path(metrics_path),
            step="update-gff-to-inframe",
            count_unit="genes",
            input_total=len(genes),
            changed_total=num_updated,
            output_total=len(genes),
            details={
                "unchanged_gene_models": num_unchanged,
                "removed_features": len(removed_row_ids),
                "synchronized_features": synchronized_features,
                "updated_gene_ids_path": str(updated_gene_ids_path) if updated_gene_ids_path else "",
            },
        )

    return {
        "updated_gene_models": num_updated,
        "unchanged_gene_models": num_unchanged,
        "removed_features": len(removed_row_ids),
        "synchronized_features": synchronized_features,
        "updated_gene_ids": updated_gene_ids,
    }
