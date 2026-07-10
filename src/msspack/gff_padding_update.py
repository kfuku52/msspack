from __future__ import annotations

from collections import OrderedDict
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import TypedDict

from .gff import child_ids, parse_attributes
from .utils import write_text


def _safe_phase(value: str) -> int:
    return int(value) if value.isdigit() else 0


def _parse_padding_log(log_path: str | Path) -> dict[str, dict[str, int]]:
    entries: dict[str, dict[str, int]] = {}
    with Path(log_path).open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if "original_seqlen=" not in line or "head_padding=" not in line:
                continue
            parts = [chunk.strip() for chunk in line.split(",")]
            record_id = parts[0]
            values: dict[str, int] = {}
            for part in parts[1:]:
                key, sep, value = part.partition("=")
                if not sep:
                    continue
                try:
                    values[key.strip()] = int(value.strip())
                except ValueError:
                    continue
            entries[record_id] = {
                "head_padding": values.get("head_padding", 0),
                "tail_padding": values.get("tail_padding", 0),
                "original_num_stop": values.get("original_num_stop", 0),
                "new_num_stop": values.get("new_num_stop", 0),
            }
    return entries


def _group_gene_records(
    gff_path: str | Path,
) -> tuple[list[str], OrderedDict[str, list[list[str]]], dict[str, str], list[str]]:
    header_lines: list[str] = []
    fasta_lines: list[str] = []
    feature_rows: list[tuple[list[str], dict[str, str]]] = []
    gene_dict: OrderedDict[str, list[list[str]]] = OrderedDict()
    transcript_to_gene: dict[str, str] = {}
    gene_order: list[str] = []

    with Path(gff_path).open("r", encoding="utf-8") as handle:
        in_fasta = False
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if in_fasta:
                fasta_lines.append(line)
                continue
            if line == "##FASTA":
                in_fasta = True
                fasta_lines.append(line)
                continue
            if not line or line.startswith("#"):
                header_lines.append(line)
                continue
            fields = line.split("\t")
            if len(fields) < 9:
                continue
            attrs = parse_attributes(fields[8])
            feature_rows.append((fields, attrs))

    for fields, attrs in feature_rows:
        feature_type = fields[2]
        if feature_type == "gene":
            gene_id = attrs.get("ID")
            if gene_id and gene_id not in gene_dict:
                gene_dict[gene_id] = []
                gene_order.append(gene_id)
        elif feature_type in ("mRNA", "transcript"):
            transcript_id = attrs.get("ID")
            for gene_id in child_ids(attrs.get("Parent")):
                if transcript_id and gene_id and transcript_id not in transcript_to_gene:
                    transcript_to_gene[transcript_id] = gene_id
                if gene_id and gene_id not in gene_dict:
                    gene_dict[gene_id] = []
                    gene_order.append(gene_id)

    for fields, attrs in feature_rows:
        feature_type = fields[2]
        matched_gene_ids: list[str] = []
        if feature_type == "gene":
            gene_id = attrs.get("ID")
            if gene_id:
                matched_gene_ids = [gene_id]
        elif feature_type in ("mRNA", "transcript"):
            matched_gene_ids = child_ids(attrs.get("Parent"))
        else:
            seen_gene_ids: set[str] = set()
            for parent_id in child_ids(attrs.get("Parent")):
                gene_id = transcript_to_gene.get(parent_id, parent_id)
                if gene_id and gene_id not in seen_gene_ids:
                    matched_gene_ids.append(gene_id)
                    seen_gene_ids.add(gene_id)
        for gene_id in matched_gene_ids:
            gene_dict.setdefault(gene_id, []).append(fields)

    ordered_gene_dict = OrderedDict((gene_id, gene_dict[gene_id]) for gene_id in gene_order)
    for gene_id, records in gene_dict.items():
        if gene_id not in ordered_gene_dict:
            ordered_gene_dict[gene_id] = records
    return header_lines, ordered_gene_dict, transcript_to_gene, fasta_lines


@dataclass
class _FeatureSpan:
    rec: list[str]
    start: int
    end: int
    phase: str


class PaddingUpdateSummary(TypedDict):
    genes_with_stops: list[str]
    updated_genes: list[str]
    forced_first_cds_warnings: list[str]


def _length_of_feature(start: int, end: int) -> int:
    return max(0, end - start + 1)


def _clip_from_5prime_plus(features: list[_FeatureSpan], leftover: int) -> int:
    index = 0
    while index < len(features) and leftover > 0:
        feature = features[index]
        feature_len = _length_of_feature(feature.start, feature.end)
        if leftover >= feature_len:
            leftover -= feature_len
            features.pop(index)
        else:
            feature.start += leftover
            leftover = 0
            index += 1
    return leftover


def _clip_from_5prime_minus(features: list[_FeatureSpan], leftover: int) -> int:
    index = 0
    while index < len(features) and leftover > 0:
        feature = features[index]
        feature_len = _length_of_feature(feature.start, feature.end)
        if leftover >= feature_len:
            leftover -= feature_len
            features.pop(index)
        else:
            feature.end -= leftover
            leftover = 0
            index += 1
    return leftover


def _clip_from_3prime_plus(features: list[_FeatureSpan], leftover: int) -> int:
    index = len(features) - 1
    while index >= 0 and leftover > 0:
        feature = features[index]
        feature_len = _length_of_feature(feature.start, feature.end)
        if leftover >= feature_len:
            leftover -= feature_len
            features.pop(index)
            index -= 1
        else:
            feature.end -= leftover
            leftover = 0
            index -= 1
    return leftover


def _clip_from_3prime_minus(features: list[_FeatureSpan], leftover: int) -> int:
    index = len(features) - 1
    while index >= 0 and leftover > 0:
        feature = features[index]
        feature_len = _length_of_feature(feature.start, feature.end)
        if leftover >= feature_len:
            leftover -= feature_len
            features.pop(index)
            index -= 1
        else:
            feature.start += leftover
            leftover = 0
            index -= 1
    return leftover


def _remove_zero_length(features: Iterable[_FeatureSpan]) -> list[_FeatureSpan]:
    return [feature for feature in features if feature.end >= feature.start]


def _update_phase(old_start: int, new_start: int, old_phase: int) -> int:
    return (old_phase + (old_start - new_start)) % 3


def apply_padding_to_gff(
    *,
    gff_path: str | Path,
    padding_log_path: str | Path,
    output_path: str | Path,
    genes_with_stops_path: str | Path,
    updated_genes_path: str | Path,
) -> PaddingUpdateSummary:
    padding_entries = _parse_padding_log(padding_log_path)
    header_lines, gene_dict, transcript_to_gene, fasta_lines = _group_gene_records(gff_path)
    gene_info: dict[str, dict[str, int]] = {}
    for record_id, padding_info in padding_entries.items():
        gene_id = transcript_to_gene.get(record_id, record_id)
        gene_info[gene_id] = padding_info

    genes_with_stops: list[str] = []
    updated_genes: list[str] = []
    final_lines: list[str] = list(header_lines)
    forced_first_cds_warnings: list[str] = []

    for gene_id, raw_records in gene_dict.items():
        info = gene_info.get(gene_id)
        if info is None:
            final_lines.extend("\t".join(rec) for rec in raw_records)
            continue

        if info["new_num_stop"] > 0:
            genes_with_stops.append(gene_id)
            final_lines.extend("\t".join(rec) for rec in raw_records)
            continue

        strand = "+"
        gene_lines = [rec for rec in raw_records if rec[2] == "gene"]
        if gene_lines:
            strand = gene_lines[0][6]
        else:
            mrna_lines = [rec for rec in raw_records if rec[2] in ("mRNA", "transcript")]
            if mrna_lines:
                strand = mrna_lines[0][6]

        exons: list[_FeatureSpan] = []
        cdss: list[_FeatureSpan] = []
        gene_structs: list[_FeatureSpan] = []
        mrna_structs: list[_FeatureSpan] = []
        others: list[_FeatureSpan] = []
        for rec in raw_records:
            item = _FeatureSpan(
                rec=rec,
                start=int(rec[3]),
                end=int(rec[4]),
                phase=rec[7],
            )
            if rec[2] == "gene":
                gene_structs.append(item)
            elif rec[2] in ("mRNA", "transcript"):
                mrna_structs.append(item)
            elif rec[2] == "exon":
                exons.append(item)
            elif rec[2] == "CDS":
                cdss.append(item)
            else:
                others.append(item)

        def sorter(item: _FeatureSpan) -> tuple[int, int]:
            return item.start, item.end

        exons.sort(key=sorter, reverse=(strand == "-"))
        cdss.sort(key=sorter, reverse=(strand == "-"))

        head_trim = (3 - info["head_padding"]) % 3
        tail_trim = (3 - info["tail_padding"]) % 3

        if head_trim > 0:
            if strand == "+":
                _clip_from_5prime_plus(exons, head_trim)
                _clip_from_5prime_plus(cdss, head_trim)
            else:
                _clip_from_5prime_minus(exons, head_trim)
                _clip_from_5prime_minus(cdss, head_trim)
        if tail_trim > 0:
            if strand == "+":
                _clip_from_3prime_plus(exons, tail_trim)
                _clip_from_3prime_plus(cdss, tail_trim)
            else:
                _clip_from_3prime_minus(exons, tail_trim)
                _clip_from_3prime_minus(cdss, tail_trim)

        exons = _remove_zero_length(exons)
        cdss = _remove_zero_length(cdss)
        if not cdss:
            forced_first_cds_warnings.append(
                f"Gene {gene_id}: padding adjustment would remove every CDS; model was left unchanged."
            )
            final_lines.extend("\t".join(rec) for rec in raw_records)
            continue
        updated_genes.append(gene_id)
        exons.sort(key=sorter)
        cdss.sort(key=sorter)

        if cdss:
            if len(cdss) == 1:
                cds = cdss[0]
                new_phase = _update_phase(int(cds.rec[3]), cds.start, _safe_phase(cds.rec[7]))
                if new_phase != 0:
                    forced_first_cds_warnings.append(
                        f"Gene {gene_id}: after 5' trim, single CDS forced phase from {new_phase} to 0."
                    )
                    new_phase = 0
                cds.phase = str(new_phase)
            else:
                first_cds = cdss[0]
                first_phase = _update_phase(
                    int(first_cds.rec[3]),
                    first_cds.start,
                    _safe_phase(first_cds.rec[7]),
                )
                if first_phase != 0:
                    forced_first_cds_warnings.append(
                        f"Gene {gene_id}: after 5' trim, first CDS forced phase from {first_phase} to 0."
                    )
                    first_phase = 0
                first_cds.phase = str(first_phase)

                last_cds = cdss[-1]
                last_cds.phase = str(
                    _update_phase(
                        int(last_cds.rec[3]),
                        last_cds.start,
                        _safe_phase(last_cds.rec[7]),
                    )
                )

        boundary_features = exons or cdss
        min_start = min(feature.start for feature in boundary_features)
        max_end = max(feature.end for feature in boundary_features)
        for feature in gene_structs + mrna_structs:
            feature.start = min_start
            feature.end = max_end
            feature.rec[3] = str(min_start)
            feature.rec[4] = str(max_end)
        for feature in exons:
            feature.rec[3] = str(feature.start)
            feature.rec[4] = str(feature.end)
        for feature in cdss:
            feature.rec[3] = str(feature.start)
            feature.rec[4] = str(feature.end)
            feature.rec[7] = feature.phase

        updated_records = [*gene_structs, *mrna_structs, *exons, *cdss, *others]
        rank_map = {"gene": 0, "mRNA": 1, "transcript": 1, "exon": 2, "CDS": 3}
        updated_records.sort(
            key=lambda feature: (
                rank_map.get(feature.rec[2], 4),
                int(feature.rec[3]),
            )
        )
        final_lines.extend("\t".join(feature.rec) for feature in updated_records)

    final_lines.extend(fasta_lines)

    write_text(Path(output_path), "\n".join(final_lines) + "\n")
    write_text(
        Path(genes_with_stops_path),
        "\n".join(genes_with_stops) + ("\n" if genes_with_stops else ""),
    )
    write_text(
        Path(updated_genes_path),
        "\n".join(updated_genes) + ("\n" if updated_genes else ""),
    )
    return {
        "genes_with_stops": genes_with_stops,
        "updated_genes": updated_genes,
        "forced_first_cds_warnings": forced_first_cds_warnings,
    }
