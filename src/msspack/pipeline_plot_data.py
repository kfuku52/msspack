from __future__ import annotations

import csv
import json
import re
from pathlib import Path

from .gff import parse_attributes
from .pipeline_plot_models import (
    GENE_SET_SPECS,
    SANKEY_COLORS,
    AnnotationConsistencyGroup,
    AnnotationConsistencySummary,
    FunctionalAnnotationGroup,
    FunctionalAnnotationSummary,
    ParsedStepRecord,
    PipelineGeneSet,
    PipelinePlotDataBundle,
    PipelinePlotMetrics,
)
from .utils import MSSPackError
from .validation import ValidationSummary, load_validation_summary

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

CONSISTENCY_THRESHOLD_RE = re.compile(
    r"identity>=(?P<identity>[0-9]+(?:\.[0-9]+)?)%;\s*"
    r"mutual coverage>=(?P<coverage>[0-9]+(?:\.[0-9]+)?)%"
)


def _parse_count(value: object | None) -> int | None:
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
    input_total: int | None = None
    changed_total: int | None = None
    output_total: int | None = None
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


def _required_count(value: int | None, label: str, path: Path) -> int:
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


def _detail_path(record: ParsedStepRecord, *keys: str) -> Path | None:
    value = _detail_value(record, *keys)
    if value is None:
        return None
    text = str(value).strip()
    return Path(text) if text else None


def _parse_step_records(log_dir: Path) -> dict[str, ParsedStepRecord]:
    return {
        "drop_duplicate_coordinate_gene": _load_step_record(
            log_dir, "06.drop-duplicate-coordinate-gene"
        ),
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

    transcript_changed_genes = _required_count(
        select_mrna.changed_total, "Changed total", select_mrna.path
    )
    genes_after_single_mrna = _required_count(
        select_mrna.output_total, "Output total", select_mrna.path
    )
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

    converted_to_misc_genes = _required_count(
        cds_to_misc.changed_total, "Changed total", cds_to_misc.path
    )
    total_cds_input = _detail_count(cds_to_misc, keys=("cds_input", "Total number of CDS in input"))
    total_cds_output = _detail_count(
        cds_to_misc, keys=("cds_output", "Total number of CDS in output")
    )
    misc_feature_output = _detail_count(
        cds_to_misc,
        keys=("misc_feature_output", "Total number of misc_feature in output"),
    )
    final_cds_genes = genes_after_padding - converted_to_misc_genes
    if final_cds_genes < 0:
        raise MSSPackError(f"Converted-to-misc genes exceed available genes in {cds_to_misc.path}")

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


def _build_gene_sets(
    output_root: Path, records: dict[str, ParsedStepRecord]
) -> tuple[PipelineGeneSet, ...]:
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


def _build_summary_payload(
    metrics: PipelinePlotMetrics,
    gene_sets: tuple[PipelineGeneSet, ...],
    functional_annotation: FunctionalAnnotationSummary | None,
    annotation_consistency: AnnotationConsistencySummary | None,
    validation_summary: ValidationSummary | None,
) -> dict[str, object]:
    payload: dict[str, object] = {
        "metrics": metrics.to_dict(),
        "sources": metrics.sources,
        "gene_sets": {gene_set.key: gene_set.to_dict() for gene_set in gene_sets},
    }
    if functional_annotation is not None:
        payload["functional_annotation"] = functional_annotation.to_dict()
    if annotation_consistency is not None:
        payload["annotation_consistency"] = annotation_consistency.to_dict()
    if validation_summary is not None:
        payload["ddbj_validation"] = validation_summary.to_dict()
    return payload


def _annotation_source_style(source: str) -> tuple[str, str]:
    normalized = source.casefold()
    if normalized == "swissprot":
        return "Swiss-Prot", SANKEY_COLORS["annotation_similarity"]
    if normalized == "uniref90":
        return "UniRef90", SANKEY_COLORS["annotation_uniref"]
    if normalized == "pfam":
        return "Pfam", SANKEY_COLORS["annotation_domain"]
    if normalized == "cdd":
        return "CDD", SANKEY_COLORS["annotation_cdd"]
    if normalized == "existing":
        return "Existing product", SANKEY_COLORS["annotation_existing"]
    if normalized == "none":
        return "Unannotated", SANKEY_COLORS["annotation_missing"]
    return source, SANKEY_COLORS["annotation_similarity"]


def _annotation_source_sort_key(source: str) -> tuple[int, str]:
    normalized = source.casefold()
    if normalized == "existing":
        return 1, normalized
    if normalized == "none":
        return 2, normalized
    return 0, normalized


def _load_plot_gene_id_map(output_root: Path) -> dict[str, str]:
    gff_path = output_root / "intermediate" / "12.gff.final-sorted.gff"
    if not gff_path.is_file():
        return {}

    gene_ids: set[str] = set()
    parents_by_id: dict[str, str] = {}
    with gff_path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise MSSPackError(
                    f"Expected 9 GFF columns at {gff_path}:{line_number}, found {len(fields)}"
                )
            try:
                attributes = parse_attributes(fields[8])
            except ValueError as exc:
                raise MSSPackError(
                    f"Invalid GFF attributes at {gff_path}:{line_number}: {exc}"
                ) from exc
            identifier = attributes.get("ID", "").strip()
            if not identifier:
                continue
            if fields[2].casefold() in {"gene", "pseudogene"}:
                gene_ids.add(identifier)
            parent = attributes.get("Parent", "").split(",", 1)[0].strip()
            if parent:
                parents_by_id[identifier] = parent

    gene_id_map = {identifier: identifier for identifier in gene_ids}
    for identifier in parents_by_id:
        current = identifier
        visited: set[str] = set()
        while current not in gene_ids and current in parents_by_id and current not in visited:
            visited.add(current)
            current = parents_by_id[current]
        if current in gene_ids:
            gene_id_map[identifier] = current
    return gene_id_map


def _plot_gene_identifier(row: dict[str, str], gene_id_map: dict[str, str]) -> str:
    feature_id = row["ID"].strip()
    if feature_id in gene_id_map:
        return gene_id_map[feature_id]
    return row["Locus_tag"].strip() or feature_id


def load_functional_annotation_summary(
    output_root: Path,
) -> FunctionalAnnotationSummary | None:
    evidence_path = output_root / "final" / "functional-annotation.tsv"
    if not evidence_path.is_file():
        return None
    gene_id_map = _load_plot_gene_id_map(output_root)
    grouped_ids: dict[str, list[str]] = {}
    seen_ids: set[str] = set()
    with evidence_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"ID", "Locus_tag", "source"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise MSSPackError(
                "Functional annotation evidence must contain ID, Locus_tag, and source columns: "
                f"{evidence_path}"
            )
        for row in reader:
            identifier = _plot_gene_identifier(row, gene_id_map)
            source = row["source"].strip() or "none"
            if not identifier:
                raise MSSPackError(
                    f"Functional annotation evidence contains an empty identifier: {evidence_path}"
                )
            if identifier in seen_ids:
                raise MSSPackError(
                    f"Functional annotation evidence contains duplicate locus tag {identifier}: "
                    f"{evidence_path}"
                )
            seen_ids.add(identifier)
            grouped_ids.setdefault(source, []).append(identifier)

    groups: list[FunctionalAnnotationGroup] = []
    used_keys: set[str] = set()
    for source in sorted(grouped_ids, key=_annotation_source_sort_key):
        base_key = re.sub(r"[^a-z0-9]+", "_", source.casefold()).strip("_") or "none"
        key = base_key
        suffix = 2
        while key in used_keys:
            key = f"{base_key}_{suffix}"
            suffix += 1
        used_keys.add(key)
        label, color = _annotation_source_style(source)
        groups.append(
            FunctionalAnnotationGroup(
                key=key,
                source=source,
                label=label,
                color=color,
                locus_tags=tuple(grouped_ids[source]),
            )
        )
    return FunctionalAnnotationSummary(path=evidence_path, groups=tuple(groups))


def load_annotation_consistency_summary(
    output_root: Path,
) -> AnnotationConsistencySummary | None:
    audit_path = output_root / "final" / "functional-annotation-consistency.tsv"
    summary_path = output_root / "final" / "functional-annotation-consistency-summary.tsv"
    source_pair_path = output_root / "final" / "functional-annotation-source-pairs.tsv"
    if not audit_path.is_file():
        return None
    missing = [path for path in (summary_path, source_pair_path) if not path.is_file()]
    if missing:
        raise MSSPackError(
            "Functional annotation consistency outputs are incomplete; missing: "
            + ", ".join(str(path) for path in missing)
        )
    gene_id_map = _load_plot_gene_id_map(output_root)
    styles = {
        "consistent": (
            "Consistent",
            SANKEY_COLORS["consistency_consistent"],
        ),
        "resolved": (
            "Auto-resolved family variation",
            SANKEY_COLORS["consistency_review"],
        ),
        "review": (
            "Needs name review",
            SANKEY_COLORS["consistency_review"],
        ),
        "no_close_family_peer": (
            "No annotated close-family peer",
            SANKEY_COLORS["consistency_no_peer"],
        ),
        "unannotated": (
            "Unannotated",
            SANKEY_COLORS["consistency_unannotated"],
        ),
    }
    grouped_ids: dict[str, list[str]] = {key: [] for key in styles}
    seen: set[str] = set()
    with audit_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"ID", "Locus_tag", "name_consistency"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise MSSPackError(
                "Functional annotation consistency evidence must contain ID, Locus_tag, "
                f"and name_consistency columns: {audit_path}"
            )
        for row in reader:
            identifier = _plot_gene_identifier(row, gene_id_map)
            status = row["name_consistency"].strip()
            if status in {"no_comparable_family", "no_near_identical_peer"}:
                # Backward-compatible reading of audit tables written before
                # the gene-level comparison tier was relaxed to 70/80.
                status = "no_close_family_peer"
            if status not in styles:
                raise MSSPackError(
                    f"Unknown functional annotation consistency status {status!r}: {audit_path}"
                )
            if not identifier or identifier in seen:
                raise MSSPackError(
                    "Functional annotation consistency evidence contains an empty or "
                    f"duplicate locus tag: {audit_path}"
                )
            seen.add(identifier)
            grouped_ids[status].append(identifier)
    groups = tuple(
        AnnotationConsistencyGroup(
            key=key,
            label=styles[key][0],
            color=styles[key][1],
            locus_tags=tuple(grouped_ids[key]),
        )
        for key in styles
        if grouped_ids[key]
    )
    identity_threshold: float | None = None
    coverage_threshold: float | None = None
    with summary_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or not {"tier", "threshold"}.issubset(reader.fieldnames):
            raise MSSPackError(
                "Functional annotation consistency summary must contain tier and threshold "
                f"columns: {summary_path}"
            )
        for row in reader:
            if row["tier"].strip() != "family":
                continue
            match = CONSISTENCY_THRESHOLD_RE.fullmatch(row["threshold"].strip())
            if match is None:
                raise MSSPackError(
                    "Invalid close-family consistency threshold in "
                    f"{summary_path}: {row['threshold']!r}"
                )
            identity_threshold = float(match.group("identity"))
            coverage_threshold = float(match.group("coverage"))
            break
    if identity_threshold is None or coverage_threshold is None:
        raise MSSPackError(
            f"Functional annotation consistency summary has no family tier: {summary_path}"
        )
    return AnnotationConsistencySummary(
        path=audit_path,
        summary_path=summary_path,
        source_pair_path=source_pair_path,
        groups=groups,
        comparison_tier="family",
        identity_threshold=identity_threshold,
        coverage_threshold=coverage_threshold,
    )


def collect_pipeline_plot_data(output_root: Path, log_dir: Path) -> PipelinePlotDataBundle:
    records = _parse_step_records(log_dir)
    metrics = _build_metrics_from_records(records)
    gene_sets = _build_gene_sets(output_root, records)
    functional_annotation = load_functional_annotation_summary(output_root)
    if (
        functional_annotation is not None
        and functional_annotation.total != metrics.genes_after_padding
    ):
        raise MSSPackError(
            "Functional annotation evidence contains "
            f"{functional_annotation.total:,} genes; expected {metrics.genes_after_padding:,}"
        )
    annotation_consistency = load_annotation_consistency_summary(output_root)
    if (
        annotation_consistency is not None
        and annotation_consistency.total != metrics.genes_after_padding
    ):
        raise MSSPackError(
            "Functional annotation consistency evidence contains "
            f"{annotation_consistency.total:,} genes; expected {metrics.genes_after_padding:,}"
        )
    validation_summary = load_pipeline_validation_summary(output_root)
    return PipelinePlotDataBundle(
        records=records,
        metrics=metrics,
        gene_sets=gene_sets,
        functional_annotation=functional_annotation,
        annotation_consistency=annotation_consistency,
        validation_summary=validation_summary,
        summary_payload=_build_summary_payload(
            metrics,
            gene_sets,
            functional_annotation,
            annotation_consistency,
            validation_summary,
        ),
    )


def load_pipeline_validation_summary(output_root: Path) -> ValidationSummary | None:
    """Load the newest validation summary produced by pack or explicit validate."""
    final_dir = output_root / "final"
    candidates = (
        final_dir / "ddbj-validation-summary.json",
        final_dir / "validation" / "ddbj-validation-summary.json",
    )
    existing = [path for path in candidates if path.is_file()]
    if not existing:
        return None
    newest = max(existing, key=lambda path: path.stat().st_mtime_ns)
    return load_validation_summary(newest)


def parse_pipeline_plot_metrics(log_dir: Path) -> PipelinePlotMetrics:
    return _build_metrics_from_records(_parse_step_records(log_dir))
