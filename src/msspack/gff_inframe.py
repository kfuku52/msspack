from __future__ import annotations

from collections import OrderedDict, defaultdict
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Optional, TypedDict

from .gff import child_ids, parse_attributes
from .step_logging import write_id_list, write_step_log, write_step_metrics


@dataclass
class _MrnaModel:
    line: list[str]
    exons: list[list[str]] = field(default_factory=list)
    cdss: list[list[str]] = field(default_factory=list)


@dataclass
class _GeneModel:
    line: list[str]
    mrnas: OrderedDict[str, _MrnaModel] = field(default_factory=OrderedDict)


class InframeFixSummary(TypedDict):
    updated_gene_models: int
    unchanged_gene_models: int
    removed_features: int
    updated_gene_ids: list[str]


def _safe_phase(value: str) -> int:
    return int(value) if value.isdigit() else 0


def _compute_mrna_boundaries(exons: list[list[str]]) -> tuple[Optional[int], Optional[int]]:
    if not exons:
        return None, None
    starts = [int(exon[3]) for exon in exons]
    ends = [int(exon[4]) for exon in exons]
    return min(starts), max(ends)


def _compute_total_cds_length(cdss: list[list[str]]) -> int:
    return sum((int(cds[4]) - int(cds[3]) + 1) for cds in cdss)


def _find_terminal_cds(cdss: list[list[str]], strand: str, *, first: bool) -> Optional[list[str]]:
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
) -> tuple[int, Optional[list[str]]]:
    remainder = _compute_total_cds_length(cdss) % 3
    if remainder == 0:
        return 0, None

    last_cds = _find_terminal_cds(cdss, strand, first=False)
    if last_cds is None:
        return 0, None
    old_start = int(last_cds[3])
    old_end = int(last_cds[4])

    if strand == "+":
        new_end = old_end - remainder
        last_cds[4] = str(new_end)
        if new_end >= old_start:
            _update_matching_exons(
                exons,
                old_start=old_start,
                old_end=old_end,
                new_start=old_start,
                new_end=new_end,
                strand=strand,
            )
    else:
        new_start = old_start + remainder
        last_cds[3] = str(new_start)
        if new_start <= old_end:
            _update_matching_exons(
                exons,
                old_start=old_start,
                old_end=old_end,
                new_start=new_start,
                new_end=old_end,
                strand=strand,
            )
    return remainder, last_cds


def _remove_cds_and_matching_exons(
    cdss: list[list[str]],
    exons: list[list[str]],
    cds_line: list[str],
    removed_ids: set[str],
) -> None:
    if cds_line in cdss:
        cdss.remove(cds_line)
    cds_id = parse_attributes(cds_line[8]).get("ID")
    if cds_id:
        removed_ids.add(cds_id)

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
        exon_id = parse_attributes(exon[8]).get("ID")
        if exon_id:
            removed_ids.add(exon_id)


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
        for raw_line in handle:
            line = raw_line.rstrip("\n")
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
        for mline in (child for child in children_of.get(gene_id, []) if child[2] == "mRNA"):
            mrna_id = parse_attributes(mline[8]).get("ID")
            if mrna_id:
                mrna_map[mrna_id] = _MrnaModel(line=mline)

    for gene_data in genes.values():
        mrna_map = gene_data.mrnas
        for mrna_id, mrna_data in mrna_map.items():
            for child in children_of.get(mrna_id, []):
                if child[2] == "exon":
                    mrna_data.exons.append(child)
                elif child[2] == "CDS":
                    mrna_data.cdss.append(child)

    removed_ids: set[str] = set()
    num_updated = 0
    num_unchanged = 0
    updated_gene_ids: list[str] = []

    for gene_id, gene_data in genes.items():
        gene_changed = False
        mrna_map = gene_data.mrnas
        for mrna_data in mrna_map.values():
            exons = mrna_data.exons
            cdss = mrna_data.cdss
            mrna_line = mrna_data.line
            if not cdss:
                continue
            strand = mrna_line[6]

            if _adjust_first_cds_frame(cdss, strand, exons):
                gene_changed = True

            remainder, last_cds = _truncate_last_cds_to_multiple_of_three(cdss, strand, exons)
            if remainder:
                gene_changed = True
                if last_cds is not None and int(last_cds[3]) > int(last_cds[4]):
                    _remove_cds_and_matching_exons(cdss, exons, last_cds, removed_ids)

            new_start, new_end = _compute_mrna_boundaries(exons)
            if new_start is None or new_end is None:
                continue
            if new_start != int(mrna_line[3]) or new_end != int(mrna_line[4]):
                gene_changed = True
            mrna_line[3] = str(new_start)
            mrna_line[4] = str(new_end)

        if gene_changed:
            exon_starts: list[int] = []
            exon_ends: list[int] = []
            for mrna_data in mrna_map.values():
                for exon_line in mrna_data.exons:
                    exon_starts.append(int(exon_line[3]))
                    exon_ends.append(int(exon_line[4]))
            if exon_starts and exon_ends:
                gene_line = gene_data.line
                gene_line[3] = str(min(exon_starts))
                gene_line[4] = str(max(exon_ends))
            num_updated += 1
            updated_gene_ids.append(gene_id)
        else:
            num_unchanged += 1

    with Path(output_path).open("w", encoding="utf-8") as handle:
        for item in all_lines:
            if isinstance(item, str):
                handle.write(item + "\n")
                continue
            feature_identifier = parse_attributes(item[8]).get("ID")
            if feature_identifier and feature_identifier in removed_ids:
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
            f"Removed features: {len(removed_ids):,}",
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
                "removed_features": len(removed_ids),
                "updated_gene_ids_path": str(updated_gene_ids_path) if updated_gene_ids_path else "",
            },
        )

    return {
        "updated_gene_models": num_updated,
        "unchanged_gene_models": num_unchanged,
        "removed_features": len(removed_ids),
        "updated_gene_ids": updated_gene_ids,
    }
