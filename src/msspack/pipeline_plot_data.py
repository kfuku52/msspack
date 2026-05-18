from __future__ import annotations

import json
import re
from itertools import combinations
from pathlib import Path
from typing import Optional

from .pipeline_plot_models import (
    GENE_SET_SPECS,
    GeneOverlapRow,
    ParsedStepRecord,
    PipelineGeneSet,
    PipelinePlotDataBundle,
    PipelinePlotMetrics,
)
from .utils import MSSPackError

COUNT_LINE_RE = re.compile(r"^(?P<label>[^|][^:]+):\s*(?P<value>.+?)\s*$")
STANDARD_LOG_FIELDS = {
    "Step",
    "Status",
    "Started at",
    "Ended at",
    "Duration seconds",
    "Count unit",
    "Input total",
    "Changed total",
    "Output total",
}


def _parse_count(value: object | None) -> Optional[int]:
    if value is None:
        return None
    if isinstance(value, int):
        return value
    if isinstance(value, float):
        return int(value)
    if not isinstance(value, str):
        raise MSSPackError(f"Unsupported count value type: {type(value).__name__}")
    stripped = value.strip()
    if stripped.lower() == "n/a":
        return None
    return int(stripped.replace(",", ""))


def _parse_step_metrics(path: Path) -> ParsedStepRecord:
    if not path.exists():
        raise MSSPackError(f"Required step metrics are missing: {path}")
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise MSSPackError(f"Invalid step metrics payload: {path}")
    details = payload.get("details")
    if not isinstance(details, dict):
        details = {}
    step = str(payload.get("step", "")).strip()
    if not step:
        raise MSSPackError(f"Missing step name in {path}")
    return ParsedStepRecord(
        path=path,
        step=step,
        count_unit=str(payload.get("count_unit", "n/a")),
        input_total=_parse_count(payload.get("input_total")),
        changed_total=_parse_count(payload.get("changed_total")),
        output_total=_parse_count(payload.get("output_total")),
        details=details,
    )


def parse_step_log(path: Path) -> ParsedStepRecord:
    if not path.exists():
        raise MSSPackError(f"Required step log is missing: {path}")
    text = path.read_text(encoding="utf-8")
    step = ""
    count_unit = ""
    input_total: Optional[int] = None
    changed_total: Optional[int] = None
    output_total: Optional[int] = None
    details: dict[str, object] = {}
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line or line.startswith("$ ") or line.startswith("|"):
            continue
        match = COUNT_LINE_RE.match(line)
        if match is None:
            continue
        label = match.group("label").strip()
        value = match.group("value").strip()
        if label == "Step":
            step = value
            continue
        if label == "Count unit":
            count_unit = value
            continue
        if label == "Input total":
            input_total = _parse_count(value)
            continue
        if label == "Changed total":
            changed_total = _parse_count(value)
            continue
        if label == "Output total":
            output_total = _parse_count(value)
            continue
        if label in STANDARD_LOG_FIELDS or label == "Summary table":
            continue
        details[label] = value
    if not step:
        raise MSSPackError(f"Could not parse step name from {path}")
    return ParsedStepRecord(
        path=path,
        step=step,
        count_unit=count_unit or "n/a",
        input_total=input_total,
        changed_total=changed_total,
        output_total=output_total,
        details=details,
    )


def _load_step_record(log_dir: Path, stem: str) -> ParsedStepRecord:
    metrics_path = log_dir / f"{stem}.metrics.json"
    if metrics_path.exists():
        return _parse_step_metrics(metrics_path)
    return parse_step_log(log_dir / f"{stem}.log")


def _required_count(value: Optional[int], label: str, path: Path) -> int:
    if value is None:
        raise MSSPackError(f"Missing {label} in {path}")
    return value


def _detail_value(record: ParsedStepRecord, *keys: str, prefix: str | None = None) -> object | None:
    for key in keys:
        if key in record.details:
            return record.details[key]
    if prefix is not None:
        for key, value in record.details.items():
            if key.startswith(prefix):
                return value
    return None


def _detail_count(
    record: ParsedStepRecord,
    *,
    keys: tuple[str, ...] = (),
    prefix: str | None = None,
) -> int:
    value = _detail_value(record, *keys, prefix=prefix)
    if value is None:
        if keys:
            joined = ", ".join(keys)
            raise MSSPackError(f"Missing detail '{joined}' in {record.path}")
        raise MSSPackError(f"Missing detail starting with '{prefix}' in {record.path}")
    parsed = _parse_count(value)
    if parsed is None:
        raise MSSPackError(f"Could not parse count detail from {record.path}")
    return parsed


def _detail_path(record: ParsedStepRecord, *keys: str) -> Optional[Path]:
    value = _detail_value(record, *keys)
    if value is None:
        return None
    text = str(value).strip()
    return Path(text) if text else None


def _parse_step_records(log_dir: Path) -> dict[str, ParsedStepRecord]:
    return {
        "drop_duplicate_coordinate_gene": _load_step_record(log_dir, "06.drop-duplicate-coordinate-gene"),
        "select_one_mrna": _load_step_record(log_dir, "07.select-one-mrna"),
        "update_gff_to_inframe": _load_step_record(log_dir, "09.update-gff-to-inframe"),
        "update_gff_with_padding": _load_step_record(log_dir, "11.update-gff-with-padding"),
        "mss_cds_to_misc": _load_step_record(log_dir, "16.mss-cds-to-misc"),
    }


def _build_metrics_from_records(records: dict[str, ParsedStepRecord]) -> PipelinePlotMetrics:
    dedup = records["drop_duplicate_coordinate_gene"]
    select_mrna = records["select_one_mrna"]
    inframe = records["update_gff_to_inframe"]
    padding = records["update_gff_with_padding"]
    cds_to_misc = records["mss_cds_to_misc"]

    initial_genes = _required_count(dedup.input_total, "Input total", dedup.path)
    duplicate_removed_genes = _required_count(dedup.changed_total, "Changed total", dedup.path)
    genes_after_dedup = _required_count(dedup.output_total, "Output total", dedup.path)

    transcript_changed_genes = _required_count(select_mrna.changed_total, "Changed total", select_mrna.path)
    genes_after_single_mrna = _required_count(select_mrna.output_total, "Output total", select_mrna.path)
    transcript_unchanged_genes = genes_after_single_mrna - transcript_changed_genes
    if transcript_unchanged_genes < 0:
        raise MSSPackError(
            f"Transcript-selection changed genes exceed total output in {select_mrna.path}"
        )
    removed_mrnas = _detail_count(select_mrna, keys=("removed_mrnas", "Removed mRNAs"))

    inframe_updated_genes = _required_count(inframe.changed_total, "Changed total", inframe.path)
    genes_after_inframe = _required_count(inframe.output_total, "Output total", inframe.path)
    inframe_unchanged_genes = _detail_count(
        inframe,
        keys=("unchanged_gene_models", "Number of unchanged gene models"),
    )

    genes_after_padding = _required_count(padding.output_total, "Output total", padding.path)
    padding_updated_genes = _detail_count(
        padding,
        keys=("updated_genes", "Number of updated genes (new_num_stop=0)"),
    )
    genes_with_stops = _detail_count(
        padding,
        keys=("genes_with_stops", "Number of genes with stops (new_num_stop>0)"),
    )
    padding_unchanged_genes = genes_after_padding - padding_updated_genes - genes_with_stops
    if padding_unchanged_genes < 0:
        raise MSSPackError(f"Padding unchanged genes became negative for {padding.path}")

    converted_to_misc_genes = _required_count(cds_to_misc.changed_total, "Changed total", cds_to_misc.path)
    total_cds_input = _detail_count(cds_to_misc, keys=("cds_input", "Total number of CDS in input"))
    total_cds_output = _detail_count(cds_to_misc, keys=("cds_output", "Total number of CDS in output"))
    misc_feature_output = _detail_count(
        cds_to_misc,
        keys=("misc_feature_output", "Total number of misc_feature in output"),
    )
    final_cds_genes = genes_after_padding - converted_to_misc_genes
    if final_cds_genes < 0:
        raise MSSPackError(
            f"Converted-to-misc genes exceed available genes in {cds_to_misc.path}"
        )

    return PipelinePlotMetrics(
        initial_genes=initial_genes,
        duplicate_removed_genes=duplicate_removed_genes,
        genes_after_dedup=genes_after_dedup,
        transcript_changed_genes=transcript_changed_genes,
        transcript_unchanged_genes=transcript_unchanged_genes,
        removed_mrnas=removed_mrnas,
        genes_after_single_mrna=genes_after_single_mrna,
        inframe_updated_genes=inframe_updated_genes,
        inframe_unchanged_genes=inframe_unchanged_genes,
        genes_after_inframe=genes_after_inframe,
        padding_updated_genes=padding_updated_genes,
        genes_with_stops=genes_with_stops,
        padding_unchanged_genes=padding_unchanged_genes,
        genes_after_padding=genes_after_padding,
        converted_to_misc_genes=converted_to_misc_genes,
        final_cds_genes=final_cds_genes,
        total_cds_input=total_cds_input,
        total_cds_output=total_cds_output,
        misc_feature_output=misc_feature_output,
        sources={key: str(value.path) for key, value in records.items()},
    )


def _read_id_list(path: Path) -> tuple[str, ...]:
    if not path.exists():
        return ()
    seen: set[str] = set()
    identifiers: list[str] = []
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        identifier = raw_line.strip()
        if identifier and identifier not in seen:
            identifiers.append(identifier)
            seen.add(identifier)
    return tuple(identifiers)


def _build_gene_sets(output_root: Path, records: dict[str, ParsedStepRecord]) -> tuple[PipelineGeneSet, ...]:
    log_dir = output_root / "logs"
    intermediate_dir = output_root / "intermediate"
    record_for_key = {
        "duplicate_removed_genes": records["drop_duplicate_coordinate_gene"],
        "transcript_changed_genes": records["select_one_mrna"],
        "inframe_updated_genes": records["update_gff_to_inframe"],
        "padding_updated_genes": records["update_gff_with_padding"],
        "genes_with_stops": records["update_gff_with_padding"],
        "converted_to_misc_genes": records["mss_cds_to_misc"],
    }
    detail_keys = {
        "duplicate_removed_genes": ("removed_gene_ids_path",),
        "transcript_changed_genes": ("changed_gene_ids_path",),
        "inframe_updated_genes": ("updated_gene_ids_path",),
        "padding_updated_genes": ("updated_genes_path",),
        "genes_with_stops": ("genes_with_stops_path",),
        "converted_to_misc_genes": ("converted_gene_ids_path",),
    }
    groups: list[PipelineGeneSet] = []
    for key, label, color, fallback_name in GENE_SET_SPECS:
        record = record_for_key[key]
        path = _detail_path(record, *detail_keys[key])
        if path is None:
            base_dir = intermediate_dir if fallback_name.startswith("11.") else log_dir
            path = base_dir / fallback_name
        groups.append(
            PipelineGeneSet(
                key=key,
                label=label,
                color=color,
                path=path,
                gene_ids=_read_id_list(path),
            )
        )
    return tuple(groups)


def _build_overlap_rows(gene_sets: tuple[PipelineGeneSet, ...]) -> tuple[GeneOverlapRow, ...]:
    nonempty = [gene_set for gene_set in gene_sets if gene_set.count > 0]
    if not nonempty:
        return ()
    id_sets = {gene_set.key: set(gene_set.gene_ids) for gene_set in nonempty}
    rows: list[GeneOverlapRow] = []
    for size in range(1, len(nonempty) + 1):
        for combo in combinations(nonempty, size):
            member_keys = tuple(gene_set.key for gene_set in combo)
            member_labels = tuple(gene_set.label for gene_set in combo)
            shared = set(id_sets[member_keys[0]])
            for gene_set in combo[1:]:
                shared &= id_sets[gene_set.key]
            if not shared:
                continue
            outside = set()
            for gene_set in nonempty:
                if gene_set.key in member_keys:
                    continue
                outside |= id_sets[gene_set.key]
            exclusive = tuple(sorted(shared - outside))
            if not exclusive:
                continue
            rows.append(
                GeneOverlapRow(
                    member_keys=member_keys,
                    member_labels=member_labels,
                    count=len(exclusive),
                    gene_ids=exclusive,
                )
            )
    rows.sort(key=lambda row: (-row.count, len(row.member_keys), row.label))
    return tuple(rows)


def _build_summary_payload(
    metrics: PipelinePlotMetrics,
    gene_sets: tuple[PipelineGeneSet, ...],
    overlap_rows: tuple[GeneOverlapRow, ...],
) -> dict[str, object]:
    return {
        "metrics": metrics.to_dict(),
        "sources": metrics.sources,
        "gene_sets": {gene_set.key: gene_set.to_dict() for gene_set in gene_sets},
        "overlap": {
            "row_count": len(overlap_rows),
            "rows": [row.to_dict() for row in overlap_rows],
        },
    }


def collect_pipeline_plot_data(output_root: Path, log_dir: Path) -> PipelinePlotDataBundle:
    records = _parse_step_records(log_dir)
    metrics = _build_metrics_from_records(records)
    gene_sets = _build_gene_sets(output_root, records)
    overlap_rows = _build_overlap_rows(gene_sets)
    return PipelinePlotDataBundle(
        records=records,
        metrics=metrics,
        gene_sets=gene_sets,
        overlap_rows=overlap_rows,
        summary_payload=_build_summary_payload(metrics, gene_sets, overlap_rows),
    )


def parse_pipeline_plot_metrics(log_dir: Path) -> PipelinePlotMetrics:
    return _build_metrics_from_records(_parse_step_records(log_dir))
