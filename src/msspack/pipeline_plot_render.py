from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import cast
from xml.sax.saxutils import escape

from .chart_primitives import (
    CHART_FONT_SIZE_PT,
    GRID_RGB,
    MUTED_RGB,
    SVG_FONT_SIZE,
    TEXT_RGB,
    write_single_page_pdf,
)
from .chart_primitives import hex_to_rgb as _hex_to_rgb
from .chart_primitives import pdf_helvetica_text_width as _pdf_helvetica_text_width
from .chart_primitives import pdf_text_command as _pdf_text_command
from .chart_primitives import pdf_top_to_bottom as _pdf_top_to_bottom
from .pipeline_plot_models import (
    SANKEY_COLORS,
    AnnotationConsistencySummary,
    EventCount,
    FunctionalAnnotationSummary,
    PipelineGeneSet,
    PipelinePlotArtifacts,
    PipelinePlotMetrics,
    SankeyBuscoSummary,
    SankeyLink,
    SankeyNode,
)
from .step_logging import count_fasta_records
from .utils import MSSPackError, ensure_dir, write_text

PDF_POINTS_PER_INCH = 72.0
SANKEY_WIDTH_IN = 7.2
SANKEY_WIDTH = SANKEY_WIDTH_IN * PDF_POINTS_PER_INCH
SANKEY_HEIGHT = 288.0
SANKEY_BUSCO_HEIGHT = 420.0
SANKEY_CONSISTENCY_BAND_HEIGHT = 112.0
SANKEY_LOWER_BAND_TOP_GAP = 8.0
SANKEY_BUSCO_BAND_HEIGHT = SANKEY_BUSCO_HEIGHT - SANKEY_HEIGHT - SANKEY_LOWER_BAND_TOP_GAP
SANKEY_SUMMARY_ROW_HEIGHT = 136.0
SANKEY_LINK_OPACITY = 0.72
SANKEY_BUSCO_COLORS = {
    "single_copy": "#2ca25f",
    "duplicated": "#3b82f6",
    "fragmented": "#f59e0b",
    "missing": "#ef4444",
}
SANKEY_BUSCO_LABELS = {
    "single_copy": "Single-copy",
    "duplicated": "Duplicated",
    "fragmented": "Fragmented",
    "missing": "Missing",
}
SANKEY_STAGE_LABELS = {
    "Input": "Input\nGFF",
    "Duplicate removal": "Coordinate\nduplicate\nremoval",
    "Transcript choice": "mRNA\nselection",
    "Frame correction": "Frame\ncorrection",
    "Padding analysis": "CDS\nboundary\nadjustment",
    "Functional annotation": "Functional\nannotation",
    "Final feature fate": "Output\nann.txt",
}
SANKEY_NODE_LABELS = {
    "start": "Input",
    "after_dedup": "Kept",
    "duplicate_removed": "Removed",
    "transcript_changed": "Reduced to\none mRNA",
    "transcript_unchanged": "Already one\nmRNA per gene",
    "inframe_updated": "Frame fixed",
    "inframe_unchanged": "Frame OK",
    "padding_updated": "Adjusted",
    "genes_with_stops": "Stops",
    "padding_unchanged": "No adjustment",
    "final_cds": "Final CDS",
    "final_misc": "Final misc",
}
_SANKEY_NODE_PRIORITY = {
    # Similarity databases precede domain fallbacks, preserved products, and misses.
    "annotation_swissprot": 0,
    "annotation_uniref90": 2,
    "annotation_pfam": 3,
    "annotation_cdd": 4,
    "annotation_existing": 5,
    "annotation_none": 6,
}
SANKEY_GENE_SET_METRICS = {
    "duplicate_removed_genes": "duplicate_removed_genes",
    "transcript_changed_genes": "transcript_changed_genes",
    "inframe_updated_genes": "inframe_updated_genes",
    "padding_updated_genes": "padding_updated_genes",
    "genes_with_stops": "genes_with_stops",
    "converted_to_misc_genes": "converted_to_misc_genes",
}


def _sankey_node_sort_key(node: SankeyNode) -> tuple[int, int, str, str]:
    priority = _SANKEY_NODE_PRIORITY.get(node.id)
    if priority is None and node.id.startswith("annotation_"):
        # A configured close-reference database sits between Swiss-Prot and
        # UniRef90. Its user-defined name cannot be enumerated in advance.
        priority = 1
    if priority is None:
        priority = 50
    return node.stage, priority, node.label.casefold(), node.id


def build_plot_artifacts(output_root: Path) -> PipelinePlotArtifacts:
    root = ensure_dir(output_root / "plots")
    return PipelinePlotArtifacts(
        root=root,
        summary_json=root / "pipeline-flow-summary.json",
        summary_tsv=root / "pipeline-flow-summary.tsv",
        gene_flow_tsv=root / "pipeline-gene-flow.tsv",
        gene_flow_svg=root / "pipeline-gene-flow.sankey.svg",
        gene_flow_pdf=root / "pipeline-gene-flow.sankey.pdf",
        event_counts_tsv=root / "pipeline-event-counts.tsv",
        event_counts_svg=root / "pipeline-event-counts.svg",
        event_counts_pdf=root / "pipeline-event-counts.pdf",
        coordinate_duplicates_tsv=root / "coordinate-duplicate-gene-models.tsv",
        coordinate_duplicates_svg=root / "coordinate-duplicate-gene-models.svg",
        coordinate_duplicates_pdf=root / "coordinate-duplicate-gene-models.pdf",
        name_consistency_tsv=root / "functional-annotation-name-consistency.tsv",
        name_consistency_svg=root / "functional-annotation-name-consistency.svg",
        name_consistency_pdf=root / "functional-annotation-name-consistency.pdf",
        source_consistency_tsv=root / "functional-annotation-source-consistency.tsv",
        source_consistency_svg=root / "functional-annotation-source-consistency.svg",
        source_consistency_pdf=root / "functional-annotation-source-consistency.pdf",
    )


def _required_busco_count(
    counts: dict[str, object],
    *,
    key: str,
    summary_label: str,
    comparison_path: Path,
) -> int:
    value = counts.get(key)
    if not isinstance(value, int) or isinstance(value, bool) or value < 0:
        raise MSSPackError(
            f"BUSCO CDS summary has an invalid {key} count for {summary_label}: {comparison_path}"
        )
    return value


def _busco_input_sequence_count(
    raw_summary: dict[str, object],
    *,
    summary_label: str,
    comparison_path: Path,
) -> int:
    value = raw_summary.get("input_sequence_count")
    if isinstance(value, int) and not isinstance(value, bool) and value > 0:
        return value
    input_fasta = raw_summary.get("input_fasta")
    if isinstance(input_fasta, str):
        fasta_path = Path(input_fasta)
        if fasta_path.is_file():
            sequence_count = count_fasta_records(fasta_path)
            if sequence_count > 0:
                return sequence_count
    raise MSSPackError(
        "BUSCO CDS summary has no usable input sequence count for "
        f"{summary_label}: {comparison_path}"
    )


def load_sankey_busco_summaries(output_root: Path) -> tuple[SankeyBuscoSummary, ...]:
    comparison_path = output_root / "busco" / "cds" / "comparison.json"
    if not comparison_path.exists():
        return ()
    try:
        payload = json.loads(comparison_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise MSSPackError(f"Could not read BUSCO CDS comparison: {comparison_path}") from exc
    if not isinstance(payload, dict) or not isinstance(payload.get("summaries"), list):
        raise MSSPackError(f"BUSCO CDS comparison has no summaries list: {comparison_path}")

    stage_specs = {
        "input": ("Input CDS", 0),
        "processed": ("Boundary-adjusted CDS", 4),
    }
    summaries: dict[str, SankeyBuscoSummary] = {}
    for raw_summary in payload["summaries"]:
        if not isinstance(raw_summary, dict):
            raise MSSPackError(
                f"BUSCO CDS comparison contains an invalid summary: {comparison_path}"
            )
        raw_label = raw_summary.get("label")
        if not isinstance(raw_label, str) or raw_label not in stage_specs:
            continue
        counts = raw_summary.get("counts")
        lineage_dataset = raw_summary.get("lineage_dataset")
        if not isinstance(counts, dict) or not isinstance(lineage_dataset, str):
            raise MSSPackError(
                f"BUSCO CDS summary is incomplete for {raw_label}: {comparison_path}"
            )
        typed_counts = cast(dict[str, object], counts)
        label, stage = stage_specs[raw_label]
        summary = SankeyBuscoSummary(
            label=label,
            stage=stage,
            lineage_dataset=lineage_dataset,
            input_sequences=_busco_input_sequence_count(
                cast(dict[str, object], raw_summary),
                summary_label=raw_label,
                comparison_path=comparison_path,
            ),
            total_buscos=_required_busco_count(
                typed_counts,
                key="total_buscos",
                summary_label=raw_label,
                comparison_path=comparison_path,
            ),
            single_copy=_required_busco_count(
                typed_counts,
                key="single_copy",
                summary_label=raw_label,
                comparison_path=comparison_path,
            ),
            duplicated=_required_busco_count(
                typed_counts,
                key="duplicated",
                summary_label=raw_label,
                comparison_path=comparison_path,
            ),
            fragmented=_required_busco_count(
                typed_counts,
                key="fragmented",
                summary_label=raw_label,
                comparison_path=comparison_path,
            ),
            missing=_required_busco_count(
                typed_counts,
                key="missing",
                summary_label=raw_label,
                comparison_path=comparison_path,
            ),
        )
        if summary.total_buscos <= 0:
            raise MSSPackError(f"BUSCO CDS total is zero for {raw_label}: {comparison_path}")
        segment_total = sum(count for _, count in summary.segment_counts())
        if segment_total != summary.total_buscos:
            raise MSSPackError(
                f"BUSCO CDS categories sum to {segment_total:,}, not {summary.total_buscos:,}, "
                f"for {raw_label}: {comparison_path}"
            )
        summaries[raw_label] = summary

    missing_labels = [label for label in stage_specs if label not in summaries]
    if missing_labels:
        raise MSSPackError(
            f"BUSCO CDS comparison is missing {', '.join(missing_labels)}: {comparison_path}"
        )
    return tuple(summaries[label] for label in stage_specs)


def _summary_rows(metrics: PipelinePlotMetrics) -> list[tuple[str, int, str, str]]:
    return [
        (
            "initial_genes",
            metrics.initial_genes,
            "genes",
            "Input gene models before duplicate removal.",
        ),
        (
            "duplicate_removed_genes",
            metrics.duplicate_removed_genes,
            "genes",
            "Genes removed because they shared identical coordinates.",
        ),
        (
            "genes_after_dedup",
            metrics.genes_after_dedup,
            "genes",
            "Genes retained after duplicate removal.",
        ),
        (
            "transcript_changed_genes",
            metrics.transcript_changed_genes,
            "genes",
            "Genes whose transcript model changed during single-mRNA selection.",
        ),
        (
            "transcript_unchanged_genes",
            metrics.transcript_unchanged_genes,
            "genes",
            "Genes unchanged by single-mRNA selection.",
        ),
        (
            "removed_mrnas",
            metrics.removed_mrnas,
            "transcripts",
            "Transcript models removed during single-mRNA selection.",
        ),
        (
            "genes_after_single_mrna",
            metrics.genes_after_single_mrna,
            "genes",
            "Genes remaining after collapsing to one mRNA per gene.",
        ),
        (
            "inframe_updated_genes",
            metrics.inframe_updated_genes,
            "genes",
            "Genes whose CDS coordinates changed during in-frame correction.",
        ),
        (
            "inframe_unchanged_genes",
            metrics.inframe_unchanged_genes,
            "genes",
            "Genes left unchanged by in-frame correction.",
        ),
        (
            "genes_after_inframe",
            metrics.genes_after_inframe,
            "genes",
            "Genes remaining after frame correction.",
        ),
        (
            "padding_updated_genes",
            metrics.padding_updated_genes,
            "genes",
            "Genes updated by padding analysis without introducing stops.",
        ),
        (
            "genes_with_stops",
            metrics.genes_with_stops,
            "genes",
            "Genes flagged as containing stops after padding analysis.",
        ),
        (
            "padding_unchanged_genes",
            metrics.padding_unchanged_genes,
            "genes",
            "Genes left unchanged by padding analysis.",
        ),
        (
            "genes_after_padding",
            metrics.genes_after_padding,
            "genes",
            "Genes carried into the MSS rendering stage.",
        ),
        (
            "converted_to_misc_genes",
            metrics.converted_to_misc_genes,
            "genes",
            "Genes whose CDS blocks were rewritten to misc_feature in the final MSS output.",
        ),
        (
            "final_cds_genes",
            metrics.final_cds_genes,
            "genes",
            "Genes that remain as CDS blocks in the final MSS output.",
        ),
        (
            "total_cds_input",
            metrics.total_cds_input,
            "features",
            "CDS feature count before CDS-to-misc conversion.",
        ),
        (
            "total_cds_output",
            metrics.total_cds_output,
            "features",
            "CDS feature count after CDS-to-misc conversion.",
        ),
        (
            "misc_feature_output",
            metrics.misc_feature_output,
            "features",
            "misc_feature count after CDS-to-misc conversion.",
        ),
    ]


def _validated_sankey_gene_id_sets(
    metrics: PipelinePlotMetrics,
    gene_sets: tuple[PipelineGeneSet, ...],
) -> dict[str, set[str]]:
    id_sets = {gene_set.key: set(gene_set.gene_ids) for gene_set in gene_sets}
    for key, metric_name in SANKEY_GENE_SET_METRICS.items():
        expected = int(getattr(metrics, metric_name))
        identifiers = id_sets.get(key)
        if identifiers is None:
            raise MSSPackError(f"Gene-ID set required for Sankey transitions is missing: {key}")
        if len(identifiers) != expected:
            raise MSSPackError(
                f"Gene-ID set '{key}' has {len(identifiers):,} IDs; expected {expected:,}"
            )

    active_total = metrics.genes_after_dedup
    stage_totals = {
        "representative transcript": (
            metrics.transcript_changed_genes + metrics.transcript_unchanged_genes
        ),
        "frame correction": metrics.inframe_updated_genes + metrics.inframe_unchanged_genes,
        "padding analysis": (
            metrics.padding_updated_genes
            + metrics.genes_with_stops
            + metrics.padding_unchanged_genes
        ),
        "final feature fate": metrics.converted_to_misc_genes + metrics.final_cds_genes,
    }
    for stage, total in stage_totals.items():
        if total != active_total:
            raise MSSPackError(
                f"Sankey stage '{stage}' contains {total:,} genes; expected {active_total:,}"
            )

    duplicate_ids = id_sets["duplicate_removed_genes"]
    active_keys = tuple(key for key in SANKEY_GENE_SET_METRICS if key != "duplicate_removed_genes")
    active_observed_ids = set().union(*(id_sets[key] for key in active_keys))
    duplicate_overlap = duplicate_ids & active_observed_ids
    if duplicate_overlap:
        raise MSSPackError(
            f"{len(duplicate_overlap):,} duplicate-removed genes also appear in later Sankey stages"
        )
    padding_overlap = id_sets["padding_updated_genes"] & id_sets["genes_with_stops"]
    if padding_overlap:
        raise MSSPackError(
            f"{len(padding_overlap):,} genes occur in both padding-updated and stop-containing sets"
        )
    if len(active_observed_ids) > active_total:
        raise MSSPackError(
            f"Sankey gene-ID sets contain {len(active_observed_ids):,} active genes; "
            f"only {active_total:,} are available"
        )
    return id_sets


def build_sankey(
    metrics: PipelinePlotMetrics,
    gene_sets: tuple[PipelineGeneSet, ...],
    functional_annotation: FunctionalAnnotationSummary | None = None,
    annotation_consistency: AnnotationConsistencySummary | None = None,
) -> tuple[list[str], list[SankeyNode], list[SankeyLink]]:
    stage_labels = [
        "Input",
        "Duplicate removal",
        "Transcript choice",
        "Frame correction",
        "Padding analysis",
    ]
    if functional_annotation is not None:
        stage_labels.append("Functional annotation")
    if annotation_consistency is not None and functional_annotation is None:
        raise MSSPackError("Annotation consistency pie data requires functional annotation data")
    stage_labels.append("Final feature fate")
    final_stage = len(stage_labels) - 1
    nodes: list[SankeyNode] = [
        SankeyNode("start", "Input genes", 0, metrics.initial_genes, SANKEY_COLORS["start"]),
        SankeyNode(
            "after_dedup",
            "After duplicate removal",
            1,
            metrics.genes_after_dedup,
            SANKEY_COLORS["kept"],
        ),
        SankeyNode(
            "final_cds",
            "Final CDS genes",
            final_stage,
            metrics.final_cds_genes,
            SANKEY_COLORS["final_cds"],
        ),
    ]
    if metrics.duplicate_removed_genes > 0:
        nodes.append(
            SankeyNode(
                "duplicate_removed",
                "Duplicate genes removed",
                1,
                metrics.duplicate_removed_genes,
                SANKEY_COLORS["removed"],
            )
        )
    if metrics.transcript_changed_genes > 0:
        nodes.append(
            SankeyNode(
                "transcript_changed",
                "Reduced to one mRNA",
                2,
                metrics.transcript_changed_genes,
                SANKEY_COLORS["transcript_changed"],
            )
        )
    if metrics.transcript_unchanged_genes > 0:
        nodes.append(
            SankeyNode(
                "transcript_unchanged",
                "Already one mRNA per gene",
                2,
                metrics.transcript_unchanged_genes,
                SANKEY_COLORS["transcript_unchanged"],
            )
        )
    if metrics.inframe_updated_genes > 0:
        nodes.append(
            SankeyNode(
                "inframe_updated",
                "Frame fixed",
                3,
                metrics.inframe_updated_genes,
                SANKEY_COLORS["inframe_updated"],
            )
        )
    if metrics.inframe_unchanged_genes > 0:
        nodes.append(
            SankeyNode(
                "inframe_unchanged",
                "Frame OK",
                3,
                metrics.inframe_unchanged_genes,
                SANKEY_COLORS["inframe_unchanged"],
            )
        )
    if metrics.padding_updated_genes > 0:
        nodes.append(
            SankeyNode(
                "padding_updated",
                "Adjusted",
                4,
                metrics.padding_updated_genes,
                SANKEY_COLORS["padding_updated"],
            )
        )
    if metrics.genes_with_stops > 0:
        nodes.append(
            SankeyNode(
                "genes_with_stops",
                "Stops",
                4,
                metrics.genes_with_stops,
                SANKEY_COLORS["genes_with_stops"],
            )
        )
    if metrics.padding_unchanged_genes > 0:
        nodes.append(
            SankeyNode(
                "padding_unchanged",
                "No adjustment",
                4,
                metrics.padding_unchanged_genes,
                SANKEY_COLORS["padding_unchanged"],
            )
        )
    if metrics.converted_to_misc_genes > 0:
        nodes.append(
            SankeyNode(
                "final_misc",
                "Final misc_feature genes",
                final_stage,
                metrics.converted_to_misc_genes,
                SANKEY_COLORS["final_misc"],
            )
        )
    if functional_annotation is not None:
        annotation_stage = stage_labels.index("Functional annotation")
        for group in functional_annotation.groups:
            if group.count > 0:
                nodes.append(
                    SankeyNode(
                        f"annotation_{group.key}",
                        group.label,
                        annotation_stage,
                        group.count,
                        group.color,
                    )
                )
    links: list[SankeyLink] = [
        SankeyLink("start", "after_dedup", metrics.genes_after_dedup, SANKEY_COLORS["kept"])
    ]
    if metrics.duplicate_removed_genes > 0:
        links.append(
            SankeyLink(
                "start",
                "duplicate_removed",
                metrics.duplicate_removed_genes,
                SANKEY_COLORS["removed"],
            )
        )
    if metrics.transcript_changed_genes > 0:
        links.append(
            SankeyLink(
                "after_dedup",
                "transcript_changed",
                metrics.transcript_changed_genes,
                SANKEY_COLORS["transcript_changed"],
            )
        )
    if metrics.transcript_unchanged_genes > 0:
        links.append(
            SankeyLink(
                "after_dedup",
                "transcript_unchanged",
                metrics.transcript_unchanged_genes,
                SANKEY_COLORS["transcript_unchanged"],
            )
        )

    id_sets = _validated_sankey_gene_id_sets(metrics, gene_sets)
    observed_ids = set().union(
        id_sets["transcript_changed_genes"],
        id_sets["inframe_updated_genes"],
        id_sets["padding_updated_genes"],
        id_sets["genes_with_stops"],
        id_sets["converted_to_misc_genes"],
    )
    residual_count = metrics.genes_after_dedup - len(observed_ids)
    active_ids = observed_ids
    if functional_annotation is not None:
        active_ids = set().union(*(set(group.locus_tags) for group in functional_annotation.groups))
        if len(active_ids) != metrics.genes_after_dedup:
            raise MSSPackError(
                "Functional annotation Sankey groups contain "
                f"{len(active_ids):,} genes; expected {metrics.genes_after_dedup:,}"
            )
        missing_event_ids = observed_ids - active_ids
        if missing_event_ids:
            raise MSSPackError(
                f"{len(missing_event_ids):,} pipeline-event genes are absent from functional "
                "annotation evidence"
            )
        residual_count = 0
    if annotation_consistency is not None:
        consistency_ids = set().union(
            *(set(group.locus_tags) for group in annotation_consistency.groups)
        )
        if consistency_ids != active_ids:
            raise MSSPackError(
                "Annotation consistency Sankey groups do not contain the same genes as "
                "functional annotation evidence"
            )
    transcript_members = {
        "transcript_unchanged": active_ids - id_sets["transcript_changed_genes"],
        "transcript_changed": id_sets["transcript_changed_genes"],
    }
    frame_members = {
        "inframe_unchanged": active_ids - id_sets["inframe_updated_genes"],
        "inframe_updated": id_sets["inframe_updated_genes"],
    }
    padding_event_ids = id_sets["padding_updated_genes"] | id_sets["genes_with_stops"]
    padding_members = {
        "padding_unchanged": active_ids - padding_event_ids,
        "padding_updated": id_sets["padding_updated_genes"],
        "genes_with_stops": id_sets["genes_with_stops"],
    }
    final_members = {
        "final_cds": active_ids - id_sets["converted_to_misc_genes"],
        "final_misc": id_sets["converted_to_misc_genes"],
    }
    node_colors = {node.id: node.color for node in nodes}

    def append_transitions(
        source_members: dict[str, set[str]],
        target_members: dict[str, set[str]],
        *,
        default_source: str,
        default_target: str,
    ) -> None:
        transition_counts: dict[tuple[str, str], int] = {}
        for source_id, source_ids in source_members.items():
            if source_id not in node_colors:
                continue
            for target_id, target_ids in target_members.items():
                if target_id not in node_colors:
                    continue
                count = len(source_ids & target_ids)
                if count > 0:
                    transition_counts[(source_id, target_id)] = count
        if residual_count > 0:
            key = (default_source, default_target)
            transition_counts[key] = transition_counts.get(key, 0) + residual_count
        for (source_id, target_id), count in transition_counts.items():
            links.append(SankeyLink(source_id, target_id, count, node_colors[target_id]))

    append_transitions(
        transcript_members,
        frame_members,
        default_source="transcript_unchanged",
        default_target="inframe_unchanged",
    )
    append_transitions(
        frame_members,
        padding_members,
        default_source="inframe_unchanged",
        default_target="padding_unchanged",
    )
    if functional_annotation is None:
        append_transitions(
            padding_members,
            final_members,
            default_source="padding_unchanged",
            default_target="final_cds",
        )
    else:
        annotation_members = {
            f"annotation_{group.key}": set(group.locus_tags)
            for group in functional_annotation.groups
        }
        default_annotation = next(iter(annotation_members))
        append_transitions(
            padding_members,
            annotation_members,
            default_source="padding_unchanged",
            default_target=default_annotation,
        )
        append_transitions(
            annotation_members,
            final_members,
            default_source=default_annotation,
            default_target="final_cds",
        )
    return stage_labels, sorted(nodes, key=_sankey_node_sort_key), links


def build_event_counts(metrics: PipelinePlotMetrics) -> list[EventCount]:
    return [
        EventCount(
            "duplicate_removed_genes",
            "Coordinate duplicates removed",
            metrics.duplicate_removed_genes,
            "genes",
            SANKEY_COLORS["removed"],
        ),
        EventCount(
            "transcript_changed_genes",
            "Genes changed by mRNA selection",
            metrics.transcript_changed_genes,
            "genes",
            SANKEY_COLORS["transcript_changed"],
        ),
        EventCount(
            "removed_mrnas", "Removed mRNAs", metrics.removed_mrnas, "transcripts", "#7c2d12"
        ),
        EventCount(
            "inframe_updated_genes",
            "Genes changed by frame correction",
            metrics.inframe_updated_genes,
            "genes",
            SANKEY_COLORS["inframe_updated"],
        ),
        EventCount(
            "padding_updated_genes",
            "CDS boundary-adjusted genes",
            metrics.padding_updated_genes,
            "genes",
            SANKEY_COLORS["padding_updated"],
        ),
        EventCount(
            "genes_with_stops",
            "Genes with stops",
            metrics.genes_with_stops,
            "genes",
            SANKEY_COLORS["genes_with_stops"],
        ),
        EventCount(
            "converted_to_misc_genes",
            "Genes converted to misc_feature",
            metrics.converted_to_misc_genes,
            "genes",
            SANKEY_COLORS["final_misc"],
        ),
    ]


def write_summary_json(payload: dict[str, object], output_path: Path) -> Path:
    return write_text(output_path, json.dumps(payload, indent=2, sort_keys=True) + "\n")


def write_summary_tsv(
    metrics: PipelinePlotMetrics,
    output_path: Path,
    functional_annotation: FunctionalAnnotationSummary | None = None,
    annotation_consistency: AnnotationConsistencySummary | None = None,
) -> Path:
    lines = ["metric\tvalue\tunit\tdescription"]
    for metric, value, unit, description in _summary_rows(metrics):
        lines.append(f"{metric}\t{value}\t{unit}\t{description}")
    if functional_annotation is not None:
        for group in functional_annotation.groups:
            lines.append(
                f"functional_annotation_{group.key}\t{group.count}\tgenes\t"
                f"Genes in the {group.label} functional-annotation outcome."
            )
    if annotation_consistency is not None:
        for consistency_group in annotation_consistency.groups:
            lines.append(
                f"annotation_consistency_{consistency_group.key}\t"
                f"{consistency_group.count}\tgenes\t"
                f"Genes in the {consistency_group.label} name-consistency outcome."
            )
    return write_text(output_path, "\n".join(lines) + "\n")


def write_gene_flow_tsv(
    stage_labels: list[str],
    nodes: list[SankeyNode],
    links: list[SankeyLink],
    output_path: Path,
) -> Path:
    node_map = {node.id: node for node in nodes}
    lines = ["source\ttarget\tcount\tsource_label\ttarget_label\tsource_stage\ttarget_stage"]
    for link in links:
        source = node_map[link.source]
        target = node_map[link.target]
        lines.append(
            "\t".join(
                [
                    source.id,
                    target.id,
                    str(link.count),
                    source.label,
                    target.label,
                    stage_labels[source.stage],
                    stage_labels[target.stage],
                ]
            )
        )
    return write_text(output_path, "\n".join(lines) + "\n")


def write_event_counts_tsv(events: list[EventCount], output_path: Path) -> Path:
    lines = ["metric\tlabel\tcount\tunit"]
    for event in events:
        lines.append(f"{event.key}\t{event.label}\t{event.count}\t{event.unit}")
    return write_text(output_path, "\n".join(lines) + "\n")


@dataclass
class _LaidOutNode:
    node: SankeyNode
    x: float
    y: float
    width: float
    height: float


@dataclass
class _LaidOutLink:
    link: SankeyLink
    source: _LaidOutNode
    target: _LaidOutNode
    source_y: float
    target_y: float
    height: float


def _sankey_stage_label(label: str) -> str:
    return SANKEY_STAGE_LABELS.get(label, label)


def _sankey_node_label(node: SankeyNode) -> str:
    return SANKEY_NODE_LABELS.get(node.id, node.label)


def _sankey_label_lines(label: str) -> tuple[str, ...]:
    return tuple(label.splitlines()) or (label,)


def _format_gene_count(count: int) -> str:
    return f"{count:,}"


def _sankey_layout(
    stage_labels: list[str],
    nodes: list[SankeyNode],
    links: list[SankeyLink],
) -> tuple[dict[str, _LaidOutNode], list[_LaidOutLink], dict[str, float]]:
    width = SANKEY_WIDTH
    chart_height = SANKEY_HEIGHT
    left = 28.0
    right = 56.0 if len(stage_labels) >= 7 else left
    top = 78.0
    bottom = 24.0
    node_width = 9.0
    stage_gap = (width - left - right - node_width) / max(1, len(stage_labels) - 1)
    content_height = chart_height - top - bottom
    node_gap = 18.0
    flow_min_height = 3.5
    node_min_height = 7.0
    totals_by_stage: dict[int, int] = {}
    for node in nodes:
        totals_by_stage[node.stage] = totals_by_stage.get(node.stage, 0) + node.count
    nodes_by_stage: dict[int, list[SankeyNode]] = {}
    for node in nodes:
        nodes_by_stage.setdefault(node.stage, []).append(node)
    max_stage_total = max(totals_by_stage.values()) if totals_by_stage else 0

    def scaled_geometry(
        scale: float,
    ) -> tuple[
        dict[tuple[str, str], float],
        dict[str, float],
        dict[str, float],
        dict[str, float],
    ]:
        display_flow_heights = {
            (link.source, link.target): (
                max(link.count * scale, flow_min_height) if link.count > 0 else 0.0
            )
            for link in links
        }
        outgoing_totals = {node.id: 0.0 for node in nodes}
        incoming_totals = {node.id: 0.0 for node in nodes}
        for link in links:
            display_height = display_flow_heights[(link.source, link.target)]
            outgoing_totals[link.source] += display_height
            incoming_totals[link.target] += display_height
        node_heights = {
            node.id: max(
                node.count * scale,
                outgoing_totals[node.id],
                incoming_totals[node.id],
                node_min_height,
            )
            for node in nodes
        }
        return display_flow_heights, outgoing_totals, incoming_totals, node_heights

    def fits(scale: float) -> bool:
        *_, node_heights = scaled_geometry(scale)
        return all(
            sum(node_heights[node.id] for node in stage_nodes)
            + max(0, len(stage_nodes) - 1) * node_gap
            <= content_height
            for stage_nodes in nodes_by_stage.values()
        )

    scale_low = 0.0
    scale_high = content_height / max_stage_total if max_stage_total > 0 else 0.0
    if fits(scale_high):
        scale = scale_high
    else:
        for _ in range(60):
            scale_mid = (scale_low + scale_high) / 2.0
            if fits(scale_mid):
                scale_low = scale_mid
            else:
                scale_high = scale_mid
        scale = scale_low

    display_flow_heights, outgoing_totals, incoming_totals, node_heights = scaled_geometry(scale)

    laid_out_nodes: dict[str, _LaidOutNode] = {}
    for stage_index, stage_nodes in nodes_by_stage.items():
        ordered_nodes = sorted(stage_nodes, key=_sankey_node_sort_key)
        heights = [node_heights[node.id] for node in ordered_nodes]
        total_height = sum(heights) + max(0, len(heights) - 1) * node_gap
        cursor = top + (content_height - total_height) / 2.0
        x = left + stage_index * stage_gap
        for node, node_height in zip(ordered_nodes, heights, strict=True):
            laid_out_nodes[node.id] = _LaidOutNode(
                node=node,
                x=x,
                y=cursor,
                width=node_width,
                height=node_height,
            )
            cursor += node_height + node_gap

    source_offsets = {node_id: 0.0 for node_id in laid_out_nodes}
    target_offsets = {node_id: 0.0 for node_id in laid_out_nodes}
    node_order = {node.id: index for index, node in enumerate(nodes)}
    ordered_links = sorted(
        links,
        key=lambda link: (
            laid_out_nodes[link.source].node.stage,
            node_order[link.source],
            node_order[link.target],
        ),
    )
    laid_out_links: list[_LaidOutLink] = []
    for link in ordered_links:
        source = laid_out_nodes[link.source]
        target = laid_out_nodes[link.target]
        display_height = display_flow_heights[(link.source, link.target)]
        source_y = source.y + source_offsets[link.source]
        target_y = target.y + target_offsets[link.target]
        source_offsets[link.source] += display_height
        target_offsets[link.target] += display_height
        laid_out_links.append(
            _LaidOutLink(
                link=link,
                source=source,
                target=target,
                source_y=source_y,
                target_y=target_y,
                height=display_height,
            )
        )

    return (
        laid_out_nodes,
        laid_out_links,
        {
            "width": width,
            "height": chart_height,
            "left": left,
            "top": top,
            "stage_gap": stage_gap,
            "node_width": node_width,
        },
    )


def _sankey_link_path(link: _LaidOutLink) -> str:
    x1 = link.source.x + link.source.width
    x2 = link.target.x
    y1 = link.source_y
    y2 = link.target_y
    height = link.height
    curve = max(10.0, (x2 - x1) * 0.35)
    return (
        f"M {x1:.2f} {y1:.2f} "
        f"C {x1 + curve:.2f} {y1:.2f}, {x2 - curve:.2f} {y2:.2f}, {x2:.2f} {y2:.2f} "
        f"L {x2:.2f} {y2 + height:.2f} "
        f"C {x2 - curve:.2f} {y2 + height:.2f}, {x1 + curve:.2f} {y1 + height:.2f}, {x1:.2f} {y1 + height:.2f} Z"
    )


def _sankey_label_anchor(node: _LaidOutNode, total_stages: int) -> tuple[float, str]:
    if total_stages >= 7:
        # With annotation evidence present, opposing labels compete for the
        # final inter-stage space. Keep every label on the downstream side and
        # reserve a wider right margin for the output node.
        return node.x + node.width + 4.0, "start"
    if node.node.stage >= total_stages - 2:
        return node.x - 4.0, "end"
    return node.x + node.width + 4.0, "start"


def _sankey_label_y(node: _LaidOutNode, total_stages: int) -> float:
    if total_stages > 6 and node.node.id == "padding_unchanged" and node.height >= 60.0:
        return node.y + node.height * 0.35
    return node.y + min(max(node.height / 2.0, 9.0), node.height - 2.0)


def _svg_multiline_text(
    *,
    x: float,
    y: float,
    anchor: str,
    css_class: str,
    lines: tuple[str, ...],
    line_height: float,
) -> str:
    tspans = [
        f'<tspan x="{x:.2f}" dy="{0.0 if index == 0 else line_height:g}">{escape(line)}</tspan>'
        for index, line in enumerate(lines)
    ]
    return (
        f'<text x="{x:.2f}" y="{y:.2f}" text-anchor="{anchor}" class="{css_class}">'
        + "".join(tspans)
        + "</text>"
    )


def _sankey_stage_x(meta: dict[str, float], stage: int) -> float:
    return meta["left"] + stage * meta["stage_gap"] + meta["node_width"] / 2.0


def _busco_band_title(summaries: tuple[SankeyBuscoSummary, ...]) -> str:
    lineages = {summary.lineage_dataset for summary in summaries}
    totals = {summary.total_buscos for summary in summaries}
    lineage = next(iter(lineages)) if len(lineages) == 1 else "mixed lineages"
    total = f"n={next(iter(totals)):,}" if len(totals) == 1 else "different n"
    return f"BUSCO results (CDS; {lineage}; {total})"


def _svg_pie_wedge(
    *,
    cx: float,
    cy: float,
    radius: float,
    start_angle: float,
    end_angle: float,
    color: str,
) -> str:
    start_x = cx + radius * math.cos(start_angle)
    start_y = cy + radius * math.sin(start_angle)
    end_x = cx + radius * math.cos(end_angle)
    end_y = cy + radius * math.sin(end_angle)
    large_arc = 1 if end_angle - start_angle > math.pi else 0
    return (
        f'<path d="M {cx:.2f} {cy:.2f} L {start_x:.2f} {start_y:.2f} '
        f'A {radius:.2f} {radius:.2f} 0 {large_arc} 1 {end_x:.2f} {end_y:.2f} Z" '
        f'fill="{color}" stroke="white" stroke-width="0.6"/>'
    )


_CONSISTENCY_PIE_ORDER = (
    "consistency_consistent",
    "consistency_resolved",
    "consistency_review",
    "consistency_no_close_family_peer",
    "consistency_unannotated",
)


def _consistency_pie_nodes(
    summary: AnnotationConsistencySummary | None,
) -> tuple[SankeyNode, ...]:
    if summary is None:
        return ()
    nodes = {
        f"consistency_{group.key}": SankeyNode(
            id=f"consistency_{group.key}",
            label=group.label,
            stage=-1,
            count=group.count,
            color=group.color,
        )
        for group in summary.groups
        if group.count > 0
    }
    return tuple(nodes[node_id] for node_id in _CONSISTENCY_PIE_ORDER if node_id in nodes)


def _sankey_band_geometry(
    *,
    has_consistency: bool,
    busco_count: int,
) -> tuple[float, float | None, float | None]:
    next_top = SANKEY_HEIGHT + SANKEY_LOWER_BAND_TOP_GAP
    consistency_top: float | None = None
    if has_consistency and busco_count == 2:
        return next_top + SANKEY_SUMMARY_ROW_HEIGHT, next_top, next_top
    if has_consistency:
        consistency_top = next_top
        next_top += SANKEY_CONSISTENCY_BAND_HEIGHT
    busco_top = next_top if busco_count else None
    if busco_top is not None:
        height = busco_top + SANKEY_BUSCO_BAND_HEIGHT
    elif consistency_top is not None:
        height = consistency_top + SANKEY_CONSISTENCY_BAND_HEIGHT
    else:
        height = SANKEY_HEIGHT
    return height, consistency_top, busco_top


def _consistency_pie_label(node: SankeyNode, total: int) -> str:
    percentage = 100.0 * node.count / total if total else 0.0
    return f"{node.label} {node.count:,} ({percentage:.1f}%)"


def _sankey_summary_guide_geometry(
    stage_labels: list[str],
    laid_out_nodes: dict[str, _LaidOutNode],
    meta: dict[str, float],
) -> tuple[tuple[float, float], ...]:
    guides: list[tuple[float, float]] = []
    for stage_label in ("Input", "Padding analysis", "Functional annotation"):
        if stage_label not in stage_labels:
            continue
        stage = stage_labels.index(stage_label)
        stage_nodes = [node for node in laid_out_nodes.values() if node.node.stage == stage]
        if not stage_nodes:
            continue
        guides.append(
            (
                _sankey_stage_x(meta, stage),
                max(node.y + node.height for node in stage_nodes),
            )
        )
    return tuple(guides)


def _append_summary_guides_svg(
    parts: list[str],
    stage_labels: list[str],
    laid_out_nodes: dict[str, _LaidOutNode],
    meta: dict[str, float],
    *,
    top_y: float,
) -> None:
    for x, start_y in _sankey_summary_guide_geometry(
        stage_labels,
        laid_out_nodes,
        meta,
    ):
        parts.append(
            f'<line x1="{x:.2f}" y1="{start_y:.2f}" x2="{x:.2f}" y2="{top_y:.2f}" '
            'stroke="#64748b" stroke-width="0.8" stroke-dasharray="3 2"/>'
        )


_CONSISTENCY_PIE_SHORT_LABELS = {
    "consistency_consistent": "Consistent",
    "consistency_resolved": "Auto-resolved variation",
    "consistency_review": "Needs review",
    "consistency_no_close_family_peer": "No close-family peer",
    "consistency_unannotated": "Unannotated",
}


def _summary_consistency_legend_lines(
    node: SankeyNode,
    percentage: float,
) -> tuple[str, ...]:
    percentage_text = f"{percentage:.1f}%"
    if node.id == "consistency_resolved":
        return ("Auto-resolved family", f"variation {percentage_text}")
    if node.id == "consistency_no_close_family_peer":
        return ("No annotated", f"close-family peer {percentage_text}")
    label = _CONSISTENCY_PIE_SHORT_LABELS.get(node.id, node.label)
    return (f"{label} {percentage_text}",)


def _consistency_tier_label(summary: AnnotationConsistencySummary) -> str:
    return {
        "near_identical": "Near-identical",
        "family": "Close family peer",
        "broad": "Broad homolog",
    }.get(summary.comparison_tier, summary.comparison_tier.replace("_", " ").title())


def _compact_percentage(value: float) -> str:
    return f"{value:g}"


def _append_summary_pie_row_svg(
    parts: list[str],
    summaries: tuple[SankeyBuscoSummary, ...],
    consistency_summary: AnnotationConsistencySummary,
    consistency_nodes: tuple[SankeyNode, ...],
    meta: dict[str, float],
    *,
    top_y: float,
) -> None:
    width = meta["width"]
    panel_width = width / 3.0
    center_y = top_y + 90.0
    radius = 17.0
    for panel_index in range(3):
        panel_x = panel_index * panel_width + 3.0
        parts.append(
            f'<rect x="{panel_x:.2f}" y="{top_y:.2f}" width="{panel_width - 6.0:.2f}" '
            f'height="{SANKEY_SUMMARY_ROW_HEIGHT - 4.0:.2f}" rx="4" fill="none" '
            'stroke="#cbd5e1" stroke-width="0.7"/>'
        )

    for panel_index, summary in enumerate(summaries):
        panel_left = panel_index * panel_width
        title_x = panel_left + panel_width / 2.0
        pie_x = panel_left + 21.0
        legend_x = pie_x + radius + 5.0
        parts.append(
            f'<text x="{title_x:.2f}" y="{top_y + 14.0:.2f}" text-anchor="middle" '
            f'class="summary-pie-title">{escape(summary.label)} BUSCO</text>'
        )
        subtitle = f"C={summary.complete_pct:.1f}%; BUSCO genes n={summary.total_buscos:,}"
        parts.append(
            f'<text x="{title_x:.2f}" y="{top_y + 28.0:.2f}" text-anchor="middle" '
            f'class="summary-pie-note">{escape(subtitle)}</text>'
        )
        parts.append(
            f'<text x="{title_x:.2f}" y="{top_y + 40.0:.2f}" text-anchor="middle" '
            f'class="summary-pie-note">CDS input n={summary.input_sequences:,}</text>'
        )
        parts.append(
            f'<text x="{title_x:.2f}" y="{top_y + 52.0:.2f}" text-anchor="middle" '
            f'class="summary-pie-note">{escape(summary.lineage_dataset)}</text>'
        )
        angle = -math.pi / 2.0
        for key, count in summary.segment_counts():
            if count <= 0:
                continue
            next_angle = angle + 2.0 * math.pi * count / summary.total_buscos
            parts.append(
                _svg_pie_wedge(
                    cx=pie_x,
                    cy=center_y,
                    radius=radius,
                    start_angle=angle,
                    end_angle=next_angle,
                    color=SANKEY_BUSCO_COLORS[key],
                )
            )
            angle = next_angle
        for legend_index, (key, count) in enumerate(summary.segment_counts()):
            legend_y = top_y + 65.0 + legend_index * 13.0
            percentage = 100.0 * count / summary.total_buscos
            parts.append(
                f'<rect x="{legend_x:.2f}" y="{legend_y - 6.0:.2f}" width="6" height="6" '
                f'rx="1" fill="{SANKEY_BUSCO_COLORS[key]}"/>'
            )
            parts.append(
                f'<text x="{legend_x + 9.0:.2f}" y="{legend_y:.2f}" '
                f'class="summary-pie-legend">{SANKEY_BUSCO_LABELS[key]} {percentage:.1f}%</text>'
            )

    total = sum(node.count for node in consistency_nodes)
    if total <= 0:
        return
    panel_left = panel_width * 2.0
    title_x = panel_left + panel_width / 2.0
    pie_x = panel_left + 21.0
    legend_x = pie_x + radius + 5.0
    parts.append(
        f'<text x="{title_x:.2f}" y="{top_y + 14.0:.2f}" text-anchor="middle" '
        f'class="summary-pie-title">Name consistency (n={total:,})</text>'
    )
    parts.append(
        f'<text x="{title_x:.2f}" y="{top_y + 28.0:.2f}" text-anchor="middle" '
        f'class="summary-pie-note">{escape(_consistency_tier_label(consistency_summary))}</text>'
    )
    parts.append(
        f'<text x="{title_x:.2f}" y="{top_y + 40.0:.2f}" text-anchor="middle" '
        f'class="summary-pie-note">id&gt;='
        f'{_compact_percentage(consistency_summary.identity_threshold)}%, cov&gt;='
        f'{_compact_percentage(consistency_summary.coverage_threshold)}%</text>'
    )
    angle = -math.pi / 2.0
    for node in consistency_nodes:
        next_angle = angle + 2.0 * math.pi * node.count / total
        parts.append(
            _svg_pie_wedge(
                cx=pie_x,
                cy=center_y,
                radius=radius,
                start_angle=angle,
                end_angle=next_angle,
                color=node.color,
            )
        )
        angle = next_angle
    legend_y = top_y + 52.0
    line_height = 10.7
    item_gap = 1.3
    for node in consistency_nodes:
        percentage = 100.0 * node.count / total
        lines = _summary_consistency_legend_lines(node, percentage)
        parts.append(
            f'<rect x="{legend_x:.2f}" y="{legend_y - 6.0:.2f}" width="6" height="6" '
            f'rx="1" fill="{node.color}"/>'
        )
        parts.append(
            _svg_multiline_text(
                x=legend_x + 9.0,
                y=legend_y,
                anchor="start",
                css_class="summary-pie-legend",
                lines=lines,
                line_height=line_height,
            )
        )
        legend_y += len(lines) * line_height + item_gap


def _append_consistency_pie_svg(
    parts: list[str],
    summary: AnnotationConsistencySummary,
    nodes: tuple[SankeyNode, ...],
    meta: dict[str, float],
    *,
    top_y: float,
) -> None:
    total = sum(node.count for node in nodes)
    if total <= 0:
        return
    center_x = meta["width"] / 2.0
    center_y = top_y + 75.0
    radius = 21.0
    width = meta["width"]
    parts.append(
        f'<line x1="16" y1="{top_y:.2f}" x2="{width - 16.0:.2f}" y2="{top_y:.2f}" '
        'stroke="#cbd5e1" stroke-width="0.7"/>'
    )
    parts.append(
        f'<text x="{center_x:.2f}" y="{top_y + 13.0:.2f}" text-anchor="middle" '
        'class="consistency-title">Name consistency (genes)</text>'
    )
    parts.append(
        f'<text x="{center_x:.2f}" y="{top_y + 25.0:.2f}" text-anchor="middle" '
        f'class="consistency-note">{escape(_consistency_tier_label(summary))} threshold</text>'
    )
    parts.append(
        f'<text x="{center_x:.2f}" y="{top_y + 37.0:.2f}" text-anchor="middle" '
        f'class="consistency-note">&gt;={_compact_percentage(summary.identity_threshold)}% '
        f'identity / &gt;={_compact_percentage(summary.coverage_threshold)}% mutual coverage</text>'
    )
    angle = -math.pi / 2.0
    for node in nodes:
        next_angle = angle + 2.0 * math.pi * node.count / total
        parts.append(
            _svg_pie_wedge(
                cx=center_x,
                cy=center_y,
                radius=radius,
                start_angle=angle,
                end_angle=next_angle,
                color=node.color,
            )
        )
        angle = next_angle
    legend_x = center_x - 212.0
    for index, node in enumerate(nodes):
        legend_y = top_y + 55.0 + index * 13.0
        parts.append(
            f'<rect x="{legend_x:.2f}" y="{legend_y - 6.0:.2f}" width="6" height="6" '
            f'rx="1" fill="{node.color}" stroke="#334155" stroke-width="0.4"/>'
        )
        parts.append(
            f'<text x="{legend_x + 9.0:.2f}" y="{legend_y:.2f}" '
            f'class="consistency-legend">{escape(_consistency_pie_label(node, total))}</text>'
        )


def _append_busco_svg(
    parts: list[str],
    summaries: tuple[SankeyBuscoSummary, ...],
    meta: dict[str, float],
    *,
    top_y: float,
) -> None:
    center_y = top_y + 54.0
    radius = 25.0
    width = meta["width"]
    parts.append(
        f'<line x1="16" y1="{top_y:.2f}" x2="{width - 16.0:.2f}" y2="{top_y:.2f}" '
        'stroke="#cbd5e1" stroke-width="0.7"/>'
    )
    parts.append(
        f'<text x="{width / 2.0:.2f}" y="{top_y + 13.0:.2f}" text-anchor="middle" class="busco-title">'
        f"{escape(_busco_band_title(summaries))}</text>"
    )
    for summary in summaries:
        center_x = _sankey_stage_x(meta, summary.stage)
        parts.append(
            f'<line x1="{center_x:.2f}" y1="{top_y - 8.0:.2f}" '
            f'x2="{center_x:.2f}" y2="{top_y + 23.0:.2f}" '
            'stroke="#94a3b8" stroke-width="0.7" stroke-dasharray="2 2"/>'
        )
        angle = -math.pi / 2.0
        for key, count in summary.segment_counts():
            if count <= 0:
                continue
            next_angle = angle + 2.0 * math.pi * count / summary.total_buscos
            parts.append(
                _svg_pie_wedge(
                    cx=center_x,
                    cy=center_y,
                    radius=radius,
                    start_angle=angle,
                    end_angle=next_angle,
                    color=SANKEY_BUSCO_COLORS[key],
                )
            )
            angle = next_angle
        legend_x = center_x + radius + 7.0
        for index, (key, count) in enumerate(summary.segment_counts()):
            legend_y = top_y + 31.0 + index * 14.0
            percentage = 100.0 * count / summary.total_buscos
            parts.append(
                f'<rect x="{legend_x:.2f}" y="{legend_y - 6.0:.2f}" width="6" height="6" '
                f'rx="1" fill="{SANKEY_BUSCO_COLORS[key]}"/>'
            )
            parts.append(
                f'<text x="{legend_x + 9.0:.2f}" y="{legend_y:.2f}" class="busco-legend">'
                f"{SANKEY_BUSCO_LABELS[key]} {percentage:.1f}%</text>"
            )
        summary_x = 8.0 if summary.stage == 0 else center_x
        summary_anchor = "start" if summary.stage == 0 else "middle"
        parts.append(
            f'<text x="{summary_x:.2f}" y="{top_y + 89.0:.2f}" '
            f'text-anchor="{summary_anchor}" class="busco-label">'
            f"{escape(summary.label)} (CDS input n={summary.input_sequences:,})</text>"
        )
        parts.append(
            f'<text x="{summary_x:.2f}" y="{top_y + 103.0:.2f}" '
            f'text-anchor="{summary_anchor}" class="busco-value">'
            f"Complete {summary.complete_pct:.1f}%</text>"
        )


def write_sankey_svg(
    stage_labels: list[str],
    nodes: list[SankeyNode],
    links: list[SankeyLink],
    output_path: Path,
    *,
    busco_summaries: tuple[SankeyBuscoSummary, ...] = (),
    annotation_consistency: AnnotationConsistencySummary | None = None,
) -> Path:
    laid_out_nodes, laid_out_links, meta = _sankey_layout(stage_labels, nodes, links)
    consistency_nodes = _consistency_pie_nodes(annotation_consistency)
    width = float(meta["width"])
    height, consistency_top, busco_top = _sankey_band_geometry(
        has_consistency=bool(consistency_nodes),
        busco_count=len(busco_summaries),
    )
    svg_width = f"{width / PDF_POINTS_PER_INCH:g}in"
    svg_height = f"{height / PDF_POINTS_PER_INCH:g}in"
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{svg_width}" height="{svg_height}" viewBox="0 0 {width:.2f} {height:.2f}">',
        f"<style>text{{font-family:Helvetica,Arial,sans-serif;font-size:{SVG_FONT_SIZE};fill:#111827}} .title{{font-weight:700}} .subtitle{{fill:#4b5563}} .stage{{font-weight:700;fill:#334155}} .count{{fill:#475569}} .busco-title,.busco-label,.consistency-title,.summary-pie-title{{font-weight:700;fill:#334155}} .busco-legend,.busco-value,.consistency-note,.summary-pie-note,.summary-pie-legend{{fill:#475569}}</style>",
        '<rect width="100%" height="100%" fill="white"/>',
        '<text x="16" y="16" class="title">Stage-wise pipeline gene flow</text>',
        '<text x="16" y="31" class="subtitle">Ribbon widths are proportional to gene counts; node values are genes.</text>',
    ]
    total_stages = len(stage_labels)
    for index, label in enumerate(stage_labels):
        stage_x = _sankey_stage_x(meta, index)
        stage_label = _sankey_stage_label(label)
        if not stage_label:
            continue
        stage_lines = _sankey_label_lines(stage_label)
        parts.append(
            _svg_multiline_text(
                x=stage_x,
                y=58.0 - (len(stage_lines) - 1) * 5.0,
                anchor="middle",
                css_class="stage",
                lines=stage_lines,
                line_height=10.0,
            )
        )
    for link in laid_out_links:
        parts.append(
            f'<path d="{_sankey_link_path(link)}" fill="{link.link.color}" '
            f'fill-opacity="{SANKEY_LINK_OPACITY:.2f}" stroke="none"/>'
        )
    for node in laid_out_nodes.values():
        parts.append(
            f'<rect x="{node.x:.2f}" y="{node.y:.2f}" width="{node.width:.2f}" height="{node.height:.2f}" rx="4" fill="{node.node.color}" stroke="#0f172a" stroke-width="0.8"/>'
        )
    for node in laid_out_nodes.values():
        label_x, anchor = _sankey_label_anchor(node, total_stages)
        label_y = _sankey_label_y(node, total_stages)
        label_lines = _sankey_label_lines(_sankey_node_label(node.node))
        label_line_height = 10.0 if len(label_lines) > 1 else 12.0
        label_start_y = label_y - 4.0 - (len(label_lines) - 1) * label_line_height / 2.0
        parts.append(
            _svg_multiline_text(
                x=label_x,
                y=label_start_y,
                anchor=anchor,
                css_class="label",
                lines=label_lines,
                line_height=label_line_height,
            )
        )
        parts.append(
            f'<text x="{label_x:.2f}" y="{label_start_y + len(label_lines) * label_line_height:.2f}" text-anchor="{anchor}" class="count">{_format_gene_count(node.node.count)}</text>'
        )
    combined_summary_row = (
        consistency_top is not None and busco_top == consistency_top and len(busco_summaries) == 2
    )
    if combined_summary_row:
        assert consistency_top is not None
        assert annotation_consistency is not None
        _append_summary_guides_svg(
            parts,
            stage_labels,
            laid_out_nodes,
            meta,
            top_y=consistency_top,
        )
        _append_summary_pie_row_svg(
            parts,
            busco_summaries,
            annotation_consistency,
            consistency_nodes,
            meta,
            top_y=consistency_top,
        )
    elif consistency_top is not None:
        assert annotation_consistency is not None
        _append_consistency_pie_svg(
            parts,
            annotation_consistency,
            consistency_nodes,
            meta,
            top_y=consistency_top,
        )
    if not combined_summary_row and busco_summaries and busco_top is not None:
        _append_busco_svg(
            parts,
            busco_summaries,
            meta,
            top_y=busco_top,
        )
    parts.append("</svg>")
    return write_text(output_path, "\n".join(parts) + "\n")


def _sankey_pdf_path_commands(link: _LaidOutLink, page_height: float) -> str:
    x1 = link.source.x + link.source.width
    x2 = link.target.x
    y1 = link.source_y
    y2 = link.target_y
    h = link.height
    curve = max(10.0, (x2 - x1) * 0.35)
    py1 = _pdf_top_to_bottom(page_height, y1)
    py2 = _pdf_top_to_bottom(page_height, y2)
    py1b = _pdf_top_to_bottom(page_height, y1 + h)
    py2b = _pdf_top_to_bottom(page_height, y2 + h)
    return (
        f"{x1:.2f} {py1:.2f} m "
        f"{x1 + curve:.2f} {py1:.2f} {x2 - curve:.2f} {py2:.2f} {x2:.2f} {py2:.2f} c "
        f"{x2:.2f} {py2b:.2f} l "
        f"{x2 - curve:.2f} {py2b:.2f} {x1 + curve:.2f} {py1b:.2f} {x1:.2f} {py1b:.2f} c h f"
    )


def _sankey_link_rgb(color: str) -> tuple[float, float, float]:
    red, green, blue = _hex_to_rgb(color)
    opacity = SANKEY_LINK_OPACITY
    return (
        1.0 - (1.0 - red) * opacity,
        1.0 - (1.0 - green) * opacity,
        1.0 - (1.0 - blue) * opacity,
    )


def _pdf_pie_wedge_commands(
    *,
    page_height: float,
    cx: float,
    cy: float,
    radius: float,
    start_angle: float,
    end_angle: float,
) -> str:
    start_x = cx + radius * math.cos(start_angle)
    start_y = cy + radius * math.sin(start_angle)
    commands = [
        f"{cx:.2f} {_pdf_top_to_bottom(page_height, cy):.2f} m",
        f"{start_x:.2f} {_pdf_top_to_bottom(page_height, start_y):.2f} l",
    ]
    angle = start_angle
    while angle < end_angle - 1e-12:
        next_angle = min(end_angle, angle + math.pi / 2.0)
        delta = next_angle - angle
        tangent_scale = 4.0 / 3.0 * math.tan(delta / 4.0)
        x1 = cx + radius * math.cos(angle)
        y1 = cy + radius * math.sin(angle)
        x2 = cx + radius * math.cos(next_angle)
        y2 = cy + radius * math.sin(next_angle)
        control1_x = x1 - tangent_scale * radius * math.sin(angle)
        control1_y = y1 + tangent_scale * radius * math.cos(angle)
        control2_x = x2 + tangent_scale * radius * math.sin(next_angle)
        control2_y = y2 - tangent_scale * radius * math.cos(next_angle)
        commands.append(
            f"{control1_x:.2f} {_pdf_top_to_bottom(page_height, control1_y):.2f} "
            f"{control2_x:.2f} {_pdf_top_to_bottom(page_height, control2_y):.2f} "
            f"{x2:.2f} {_pdf_top_to_bottom(page_height, y2):.2f} c"
        )
        angle = next_angle
    commands.append("h f")
    return " ".join(commands)


def _centered_pdf_text_command(
    *,
    page_height: float,
    center_x: float,
    y_top: float,
    text: str,
    font: str,
    size: int,
    color: tuple[float, float, float],
    bold: bool = False,
) -> str:
    return _pdf_text_command(
        page_height=page_height,
        x=center_x - _pdf_helvetica_text_width(text, size=size, bold=bold) / 2.0,
        y_top=y_top,
        text=text,
        font=font,
        size=size,
        color=color,
    )


def _append_consistency_pie_pdf(
    commands: list[str],
    summary: AnnotationConsistencySummary,
    nodes: tuple[SankeyNode, ...],
    meta: dict[str, float],
    page_height: float,
    *,
    top_y: float,
) -> None:
    total = sum(node.count for node in nodes)
    if total <= 0:
        return
    center_x = meta["width"] / 2.0
    center_y = top_y + 75.0
    radius = 21.0
    width = meta["width"]
    separator_y = _pdf_top_to_bottom(page_height, top_y)
    commands.append(
        f"{GRID_RGB[0]:.3f} {GRID_RGB[1]:.3f} {GRID_RGB[2]:.3f} RG 0.7 w "
        f"16 {separator_y:.2f} m {width - 16.0:.2f} {separator_y:.2f} l S"
    )
    commands.append(
        _centered_pdf_text_command(
            page_height=page_height,
            center_x=center_x,
            y_top=top_y + 13.0,
            text="Name consistency (genes)",
            font="F2",
            size=CHART_FONT_SIZE_PT,
            color=MUTED_RGB,
            bold=True,
        )
    )
    commands.append(
        _centered_pdf_text_command(
            page_height=page_height,
            center_x=center_x,
            y_top=top_y + 25.0,
            text=f"{_consistency_tier_label(summary)} threshold",
            font="F1",
            size=CHART_FONT_SIZE_PT,
            color=MUTED_RGB,
        )
    )
    commands.append(
        _centered_pdf_text_command(
            page_height=page_height,
            center_x=center_x,
            y_top=top_y + 37.0,
            text=(
                f">={_compact_percentage(summary.identity_threshold)}% identity / "
                f">={_compact_percentage(summary.coverage_threshold)}% mutual coverage"
            ),
            font="F1",
            size=CHART_FONT_SIZE_PT,
            color=MUTED_RGB,
        )
    )
    angle = -math.pi / 2.0
    for node in nodes:
        next_angle = angle + 2.0 * math.pi * node.count / total
        red, green, blue = _hex_to_rgb(node.color)
        commands.append(f"{red:.3f} {green:.3f} {blue:.3f} rg")
        commands.append(
            _pdf_pie_wedge_commands(
                page_height=page_height,
                cx=center_x,
                cy=center_y,
                radius=radius,
                start_angle=angle,
                end_angle=next_angle,
            )
        )
        angle = next_angle
    legend_x = center_x - 212.0
    for index, node in enumerate(nodes):
        legend_y = top_y + 55.0 + index * 13.0
        red, green, blue = _hex_to_rgb(node.color)
        rect_y = _pdf_top_to_bottom(page_height, legend_y - 6.0, 6.0)
        commands.append(f"{red:.3f} {green:.3f} {blue:.3f} rg {legend_x:.2f} {rect_y:.2f} 6 6 re f")
        commands.append(
            _pdf_text_command(
                page_height=page_height,
                x=legend_x + 9.0,
                y_top=legend_y,
                text=_consistency_pie_label(node, total),
                font="F1",
                size=CHART_FONT_SIZE_PT,
                color=MUTED_RGB,
            )
        )


def _append_summary_guides_pdf(
    commands: list[str],
    stage_labels: list[str],
    laid_out_nodes: dict[str, _LaidOutNode],
    meta: dict[str, float],
    page_height: float,
    *,
    top_y: float,
) -> None:
    for x, start_y in _sankey_summary_guide_geometry(
        stage_labels,
        laid_out_nodes,
        meta,
    ):
        pdf_start = _pdf_top_to_bottom(page_height, start_y)
        pdf_end = _pdf_top_to_bottom(page_height, top_y)
        commands.append(
            f"0.392 0.455 0.545 RG 0.8 w [3 2] 0 d {x:.2f} {pdf_start:.2f} m "
            f"{x:.2f} {pdf_end:.2f} l S [] 0 d"
        )


def _append_summary_pie_row_pdf(
    commands: list[str],
    summaries: tuple[SankeyBuscoSummary, ...],
    consistency_summary: AnnotationConsistencySummary,
    consistency_nodes: tuple[SankeyNode, ...],
    meta: dict[str, float],
    page_height: float,
    *,
    top_y: float,
) -> None:
    width = meta["width"]
    panel_width = width / 3.0
    center_y = top_y + 90.0
    radius = 17.0
    box_height = SANKEY_SUMMARY_ROW_HEIGHT - 4.0
    box_y = _pdf_top_to_bottom(page_height, top_y, box_height)
    for panel_index in range(3):
        panel_x = panel_index * panel_width + 3.0
        commands.append(
            f"{GRID_RGB[0]:.3f} {GRID_RGB[1]:.3f} {GRID_RGB[2]:.3f} RG 0.7 w "
            f"{panel_x:.2f} {box_y:.2f} {panel_width - 6.0:.2f} {box_height:.2f} re S"
        )

    for panel_index, summary in enumerate(summaries):
        panel_left = panel_index * panel_width
        title_x = panel_left + panel_width / 2.0
        pie_x = panel_left + 21.0
        legend_x = pie_x + radius + 5.0
        commands.append(
            _centered_pdf_text_command(
                page_height=page_height,
                center_x=title_x,
                y_top=top_y + 14.0,
                text=f"{summary.label} BUSCO",
                font="F2",
                size=CHART_FONT_SIZE_PT,
                color=MUTED_RGB,
                bold=True,
            )
        )
        commands.append(
            _centered_pdf_text_command(
                page_height=page_height,
                center_x=title_x,
                y_top=top_y + 28.0,
                text=(
                    f"C={summary.complete_pct:.1f}%; "
                    f"BUSCO genes n={summary.total_buscos:,}"
                ),
                font="F1",
                size=CHART_FONT_SIZE_PT,
                color=MUTED_RGB,
            )
        )
        commands.append(
            _centered_pdf_text_command(
                page_height=page_height,
                center_x=title_x,
                y_top=top_y + 40.0,
                text=f"CDS input n={summary.input_sequences:,}",
                font="F1",
                size=CHART_FONT_SIZE_PT,
                color=MUTED_RGB,
            )
        )
        commands.append(
            _centered_pdf_text_command(
                page_height=page_height,
                center_x=title_x,
                y_top=top_y + 52.0,
                text=summary.lineage_dataset,
                font="F1",
                size=CHART_FONT_SIZE_PT,
                color=MUTED_RGB,
            )
        )
        angle = -math.pi / 2.0
        for key, count in summary.segment_counts():
            if count <= 0:
                continue
            next_angle = angle + 2.0 * math.pi * count / summary.total_buscos
            red, green, blue = _hex_to_rgb(SANKEY_BUSCO_COLORS[key])
            commands.append(f"{red:.3f} {green:.3f} {blue:.3f} rg")
            commands.append(
                _pdf_pie_wedge_commands(
                    page_height=page_height,
                    cx=pie_x,
                    cy=center_y,
                    radius=radius,
                    start_angle=angle,
                    end_angle=next_angle,
                )
            )
            angle = next_angle
        for legend_index, (key, count) in enumerate(summary.segment_counts()):
            legend_y = top_y + 65.0 + legend_index * 13.0
            percentage = 100.0 * count / summary.total_buscos
            red, green, blue = _hex_to_rgb(SANKEY_BUSCO_COLORS[key])
            rect_y = _pdf_top_to_bottom(page_height, legend_y - 6.0, 6.0)
            commands.append(
                f"{red:.3f} {green:.3f} {blue:.3f} rg {legend_x:.2f} {rect_y:.2f} 6 6 re f"
            )
            commands.append(
                _pdf_text_command(
                    page_height=page_height,
                    x=legend_x + 9.0,
                    y_top=legend_y,
                    text=f"{SANKEY_BUSCO_LABELS[key]} {percentage:.1f}%",
                    font="F1",
                    size=CHART_FONT_SIZE_PT,
                    color=MUTED_RGB,
                )
            )

    total = sum(node.count for node in consistency_nodes)
    if total <= 0:
        return
    panel_left = panel_width * 2.0
    title_x = panel_left + panel_width / 2.0
    pie_x = panel_left + 21.0
    legend_x = pie_x + radius + 5.0
    commands.append(
        _centered_pdf_text_command(
            page_height=page_height,
            center_x=title_x,
            y_top=top_y + 14.0,
            text=f"Name consistency (n={total:,})",
            font="F2",
            size=CHART_FONT_SIZE_PT,
            color=MUTED_RGB,
            bold=True,
        )
    )
    commands.append(
        _centered_pdf_text_command(
            page_height=page_height,
            center_x=title_x,
            y_top=top_y + 28.0,
            text=_consistency_tier_label(consistency_summary),
            font="F1",
            size=CHART_FONT_SIZE_PT,
            color=MUTED_RGB,
        )
    )
    commands.append(
        _centered_pdf_text_command(
            page_height=page_height,
            center_x=title_x,
            y_top=top_y + 40.0,
            text=(
                f"id>={_compact_percentage(consistency_summary.identity_threshold)}%, "
                f"cov>={_compact_percentage(consistency_summary.coverage_threshold)}%"
            ),
            font="F1",
            size=CHART_FONT_SIZE_PT,
            color=MUTED_RGB,
        )
    )
    angle = -math.pi / 2.0
    for node in consistency_nodes:
        next_angle = angle + 2.0 * math.pi * node.count / total
        red, green, blue = _hex_to_rgb(node.color)
        commands.append(f"{red:.3f} {green:.3f} {blue:.3f} rg")
        commands.append(
            _pdf_pie_wedge_commands(
                page_height=page_height,
                cx=pie_x,
                cy=center_y,
                radius=radius,
                start_angle=angle,
                end_angle=next_angle,
            )
        )
        angle = next_angle
    legend_y = top_y + 52.0
    line_height = 10.7
    item_gap = 1.3
    for node in consistency_nodes:
        percentage = 100.0 * node.count / total
        lines = _summary_consistency_legend_lines(node, percentage)
        red, green, blue = _hex_to_rgb(node.color)
        rect_y = _pdf_top_to_bottom(page_height, legend_y - 6.0, 6.0)
        commands.append(f"{red:.3f} {green:.3f} {blue:.3f} rg {legend_x:.2f} {rect_y:.2f} 6 6 re f")
        for line_index, line in enumerate(lines):
            commands.append(
                _pdf_text_command(
                    page_height=page_height,
                    x=legend_x + 9.0,
                    y_top=legend_y + line_index * line_height,
                    text=line,
                    font="F1",
                    size=CHART_FONT_SIZE_PT,
                    color=MUTED_RGB,
                )
            )
        legend_y += len(lines) * line_height + item_gap


def _append_busco_pdf(
    commands: list[str],
    summaries: tuple[SankeyBuscoSummary, ...],
    meta: dict[str, float],
    page_height: float,
    *,
    top_y: float,
) -> None:
    center_y = top_y + 54.0
    radius = 25.0
    width = meta["width"]
    separator_y = _pdf_top_to_bottom(page_height, top_y)
    commands.append(
        f"{GRID_RGB[0]:.3f} {GRID_RGB[1]:.3f} {GRID_RGB[2]:.3f} RG 0.7 w "
        f"16 {separator_y:.2f} m {width - 16.0:.2f} {separator_y:.2f} l S"
    )
    commands.append(
        _centered_pdf_text_command(
            page_height=page_height,
            center_x=width / 2.0,
            y_top=top_y + 13.0,
            text=_busco_band_title(summaries),
            font="F2",
            size=CHART_FONT_SIZE_PT,
            color=MUTED_RGB,
            bold=True,
        )
    )
    for summary in summaries:
        center_x = _sankey_stage_x(meta, summary.stage)
        connector_top = _pdf_top_to_bottom(page_height, top_y - 8.0)
        connector_bottom = _pdf_top_to_bottom(page_height, top_y + 23.0)
        commands.append(
            f"0.58 0.64 0.72 RG 0.7 w [2 2] 0 d {center_x:.2f} {connector_top:.2f} m "
            f"{center_x:.2f} {connector_bottom:.2f} l S [] 0 d"
        )
        angle = -math.pi / 2.0
        for key, count in summary.segment_counts():
            if count <= 0:
                continue
            next_angle = angle + 2.0 * math.pi * count / summary.total_buscos
            red, green, blue = _hex_to_rgb(SANKEY_BUSCO_COLORS[key])
            commands.append(f"{red:.3f} {green:.3f} {blue:.3f} rg")
            commands.append(
                _pdf_pie_wedge_commands(
                    page_height=page_height,
                    cx=center_x,
                    cy=center_y,
                    radius=radius,
                    start_angle=angle,
                    end_angle=next_angle,
                )
            )
            angle = next_angle
        legend_x = center_x + radius + 7.0
        for index, (key, count) in enumerate(summary.segment_counts()):
            legend_y = top_y + 31.0 + index * 14.0
            percentage = 100.0 * count / summary.total_buscos
            red, green, blue = _hex_to_rgb(SANKEY_BUSCO_COLORS[key])
            rect_y = _pdf_top_to_bottom(page_height, legend_y - 6.0, 6.0)
            commands.append(
                f"{red:.3f} {green:.3f} {blue:.3f} rg {legend_x:.2f} {rect_y:.2f} 6 6 re f"
            )
            commands.append(
                _pdf_text_command(
                    page_height=page_height,
                    x=legend_x + 9.0,
                    y_top=legend_y,
                    text=f"{SANKEY_BUSCO_LABELS[key]} {percentage:.1f}%",
                    font="F1",
                    size=CHART_FONT_SIZE_PT,
                    color=MUTED_RGB,
                )
            )
        if summary.stage == 0:
            commands.append(
                _pdf_text_command(
                    page_height=page_height,
                    x=8.0,
                    y_top=top_y + 89.0,
                    text=f"{summary.label} (CDS input n={summary.input_sequences:,})",
                    font="F2",
                    size=CHART_FONT_SIZE_PT,
                    color=MUTED_RGB,
                )
            )
            commands.append(
                _pdf_text_command(
                    page_height=page_height,
                    x=8.0,
                    y_top=top_y + 103.0,
                    text=f"Complete {summary.complete_pct:.1f}%",
                    font="F1",
                    size=CHART_FONT_SIZE_PT,
                    color=MUTED_RGB,
                )
            )
        else:
            commands.append(
                _centered_pdf_text_command(
                    page_height=page_height,
                    center_x=center_x,
                    y_top=top_y + 89.0,
                    text=f"{summary.label} (CDS input n={summary.input_sequences:,})",
                    font="F2",
                    size=CHART_FONT_SIZE_PT,
                    color=MUTED_RGB,
                    bold=True,
                )
            )
            commands.append(
                _centered_pdf_text_command(
                    page_height=page_height,
                    center_x=center_x,
                    y_top=top_y + 103.0,
                    text=f"Complete {summary.complete_pct:.1f}%",
                    font="F1",
                    size=CHART_FONT_SIZE_PT,
                    color=MUTED_RGB,
                )
            )


def write_sankey_pdf(
    stage_labels: list[str],
    nodes: list[SankeyNode],
    links: list[SankeyLink],
    output_path: Path,
    *,
    busco_summaries: tuple[SankeyBuscoSummary, ...] = (),
    annotation_consistency: AnnotationConsistencySummary | None = None,
) -> Path:
    laid_out_nodes, laid_out_links, meta = _sankey_layout(stage_labels, nodes, links)
    consistency_nodes = _consistency_pie_nodes(annotation_consistency)
    width = float(meta["width"])
    height, consistency_top, busco_top = _sankey_band_geometry(
        has_consistency=bool(consistency_nodes),
        busco_count=len(busco_summaries),
    )
    commands = [
        f"1 1 1 rg 0 0 {width:.2f} {height:.2f} re f",
        _pdf_text_command(
            page_height=height,
            x=16,
            y_top=16,
            text="Stage-wise pipeline gene flow",
            font="F2",
            size=CHART_FONT_SIZE_PT,
            color=TEXT_RGB,
        ),
        _pdf_text_command(
            page_height=height,
            x=16,
            y_top=31,
            text="Ribbon widths are proportional to gene counts; node values are genes.",
            font="F1",
            size=CHART_FONT_SIZE_PT,
            color=MUTED_RGB,
        ),
    ]
    total_stages = len(stage_labels)
    for index, label in enumerate(stage_labels):
        stage_label = _sankey_stage_label(label)
        if not stage_label:
            continue
        stage_lines = _sankey_label_lines(stage_label)
        stage_x = _sankey_stage_x(meta, index)
        stage_start_y = 58.0 - (len(stage_lines) - 1) * 5.0
        for line_index, stage_line in enumerate(stage_lines):
            commands.append(
                _pdf_text_command(
                    page_height=height,
                    x=stage_x
                    - _pdf_helvetica_text_width(
                        stage_line,
                        size=CHART_FONT_SIZE_PT,
                        bold=True,
                    )
                    / 2.0,
                    y_top=stage_start_y + line_index * 10.0,
                    text=stage_line,
                    font="F2",
                    size=CHART_FONT_SIZE_PT,
                    color=MUTED_RGB,
                )
            )
    for link in laid_out_links:
        r, g, b = _sankey_link_rgb(link.link.color)
        commands.append(f"{r:.3f} {g:.3f} {b:.3f} rg")
        commands.append(_sankey_pdf_path_commands(link, height))
    for node in laid_out_nodes.values():
        r, g, b = _hex_to_rgb(node.node.color)
        rect_y = _pdf_top_to_bottom(height, node.y, node.height)
        commands.append(
            f"{r:.3f} {g:.3f} {b:.3f} rg {node.x:.2f} {rect_y:.2f} {node.width:.2f} {node.height:.2f} re f"
        )
        commands.append(
            f"{TEXT_RGB[0]:.3f} {TEXT_RGB[1]:.3f} {TEXT_RGB[2]:.3f} RG 0.8 w {node.x:.2f} {rect_y:.2f} {node.width:.2f} {node.height:.2f} re S"
        )
    for node in laid_out_nodes.values():
        label_x, anchor = _sankey_label_anchor(node, total_stages)
        label_y = _sankey_label_y(node, total_stages)
        label_lines = _sankey_label_lines(_sankey_node_label(node.node))
        label_line_height = 10.0 if len(label_lines) > 1 else 12.0
        label_start_y = label_y - 4.0 - (len(label_lines) - 1) * label_line_height / 2.0
        count = _format_gene_count(node.node.count)
        for line_index, label_line in enumerate(label_lines):
            line_x = label_x
            if anchor == "end":
                line_x = max(
                    10.0,
                    label_x
                    - _pdf_helvetica_text_width(
                        label_line,
                        size=CHART_FONT_SIZE_PT,
                    ),
                )
            commands.append(
                _pdf_text_command(
                    page_height=height,
                    x=line_x,
                    y_top=label_start_y + line_index * label_line_height,
                    text=label_line,
                    font="F1",
                    size=CHART_FONT_SIZE_PT,
                    color=TEXT_RGB,
                )
            )
        count_x = label_x
        if anchor == "end":
            count_x = max(
                10.0,
                label_x - _pdf_helvetica_text_width(count, size=CHART_FONT_SIZE_PT),
            )
        commands.append(
            _pdf_text_command(
                page_height=height,
                x=count_x,
                y_top=label_start_y + len(label_lines) * label_line_height,
                text=count,
                font="F1",
                size=CHART_FONT_SIZE_PT,
                color=MUTED_RGB,
            )
        )
    combined_summary_row = (
        consistency_top is not None and busco_top == consistency_top and len(busco_summaries) == 2
    )
    if combined_summary_row:
        assert consistency_top is not None
        assert annotation_consistency is not None
        _append_summary_guides_pdf(
            commands,
            stage_labels,
            laid_out_nodes,
            meta,
            height,
            top_y=consistency_top,
        )
        _append_summary_pie_row_pdf(
            commands,
            busco_summaries,
            annotation_consistency,
            consistency_nodes,
            meta,
            height,
            top_y=consistency_top,
        )
    elif consistency_top is not None:
        assert annotation_consistency is not None
        _append_consistency_pie_pdf(
            commands,
            annotation_consistency,
            consistency_nodes,
            meta,
            height,
            top_y=consistency_top,
        )
    if not combined_summary_row and busco_summaries and busco_top is not None:
        _append_busco_pdf(
            commands,
            busco_summaries,
            meta,
            height,
            top_y=busco_top,
        )
    return write_single_page_pdf(
        width=width, height=height, commands=commands, output_path=output_path
    )


def _nice_axis_max(max_count: int) -> int:
    if max_count <= 0:
        return 1
    step = 10 ** max(0, math.floor(math.log10(max_count / 4.0)))
    return int(math.ceil(max_count / step) * step)


def _event_chart_geometry(events: list[EventCount]) -> tuple[float, ...]:
    width = SANKEY_WIDTH
    left = 190.0
    top = 78.0
    bar_height = 16.0
    row_gap = 31.0
    bar_width = 273.0
    height = top + max(0, len(events) - 1) * row_gap + bar_height + 18.0
    return width, height, left, top, bar_height, row_gap, bar_width


def _event_ticks(axis_max: int) -> tuple[float, ...]:
    return tuple(axis_max * fraction / 4.0 for fraction in range(5))


def _format_axis_tick(value: float) -> str:
    return f"{int(value):,}" if value.is_integer() else f"{value:,.1f}"


def write_event_counts_svg(events: list[EventCount], output_path: Path) -> Path:
    width, height, left, top, bar_height, row_gap, bar_width = _event_chart_geometry(events)
    max_count = max((event.count for event in events), default=0)
    axis_max = _nice_axis_max(max_count)
    grid_top = top - 8.0
    grid_bottom = height - 18.0
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width / PDF_POINTS_PER_INCH:g}in" height="{height / PDF_POINTS_PER_INCH:g}in" viewBox="0 0 {width:.2f} {height:.2f}">',
        f"<style>text{{font-family:Helvetica,Arial,sans-serif;fill:#111827}} .title{{font-size:{SVG_FONT_SIZE};font-weight:700}} .subtitle{{font-size:{SVG_FONT_SIZE};fill:#4b5563}} .label{{font-size:{SVG_FONT_SIZE};font-weight:700}} .unit{{font-size:{SVG_FONT_SIZE};fill:#64748b}} .value{{font-size:{SVG_FONT_SIZE};fill:#334155}} .tick{{font-size:{SVG_FONT_SIZE};fill:#6b7280}}</style>",
        '<rect width="100%" height="100%" fill="white"/>',
        '<text x="16" y="16" class="title">Pipeline event counts</text>',
        '<text x="16" y="31" class="subtitle">Step-level counts from packaging logs. Removed mRNAs are transcript counts; the other bars are gene counts.</text>',
    ]
    for tick in _event_ticks(axis_max):
        x = left + bar_width * tick / axis_max
        parts.append(
            f'<line x1="{x:.2f}" y1="{grid_top:.2f}" x2="{x:.2f}" y2="{grid_bottom:.2f}" '
            'stroke="#e2e8f0" stroke-width="0.7"/>'
        )
        parts.append(
            f'<text x="{x:.2f}" y="62" text-anchor="middle" class="tick">'
            f"{_format_axis_tick(tick)}</text>"
        )
    for index, event in enumerate(events):
        y = top + index * row_gap
        bar_len = bar_width * event.count / axis_max
        parts.append(f'<text x="16" y="{y + 7.0:.2f}" class="label">{escape(event.label)}</text>')
        parts.append(f'<text x="16" y="{y + 20.0:.2f}" class="unit">{escape(event.unit)}</text>')
        parts.append(
            f'<rect x="{left:.2f}" y="{y:.2f}" width="{bar_width:.2f}" height="{bar_height:.2f}" '
            'rx="2" fill="#f8fafc" stroke="#cbd5e1" stroke-width="0.7"/>'
        )
        if event.count > 0:
            if bar_len >= 1.2:
                parts.append(
                    f'<rect x="{left:.2f}" y="{y:.2f}" width="{bar_len:.2f}" '
                    f'height="{bar_height:.2f}" rx="1" fill="{event.color}"/>'
                )
            else:
                marker_x = left + bar_len
                parts.append(
                    f'<line x1="{marker_x:.2f}" y1="{y:.2f}" x2="{marker_x:.2f}" '
                    f'y2="{y + bar_height:.2f}" stroke="{event.color}" stroke-width="1.2"/>'
                )
        parts.append(
            f'<text x="{left + bar_width + 8.0:.2f}" y="{y + 11.0:.2f}" class="value">'
            f"{event.count:,}</text>"
        )
    parts.append("</svg>")
    return write_text(output_path, "\n".join(parts) + "\n")


def write_event_counts_pdf(events: list[EventCount], output_path: Path) -> Path:
    width, height, left, top, bar_height, row_gap, bar_width = _event_chart_geometry(events)
    max_count = max((event.count for event in events), default=0)
    axis_max = _nice_axis_max(max_count)
    grid_top = top - 8.0
    grid_bottom = height - 18.0
    commands = [
        f"1 1 1 rg 0 0 {width:.2f} {height:.2f} re f",
        _pdf_text_command(
            page_height=height,
            x=16,
            y_top=16,
            text="Pipeline event counts",
            font="F2",
            size=CHART_FONT_SIZE_PT,
            color=TEXT_RGB,
        ),
        _pdf_text_command(
            page_height=height,
            x=16,
            y_top=31,
            text="Step-level counts from packaging logs. Removed mRNAs are transcript counts; the other bars are gene counts.",
            font="F1",
            size=CHART_FONT_SIZE_PT,
            color=MUTED_RGB,
        ),
    ]
    for tick in _event_ticks(axis_max):
        x = left + bar_width * tick / axis_max
        y1 = _pdf_top_to_bottom(height, grid_top)
        y2 = _pdf_top_to_bottom(height, grid_bottom)
        commands.append(
            f"{GRID_RGB[0]:.3f} {GRID_RGB[1]:.3f} {GRID_RGB[2]:.3f} RG 0.7 w {x:.2f} {y1:.2f} m {x:.2f} {y2:.2f} l S"
        )
        tick_text = _format_axis_tick(tick)
        commands.append(
            _centered_pdf_text_command(
                page_height=height,
                center_x=x,
                y_top=62.0,
                text=tick_text,
                font="F1",
                size=CHART_FONT_SIZE_PT,
                color=MUTED_RGB,
            )
        )
    for index, event in enumerate(events):
        y = top + index * row_gap
        bar_len = bar_width * event.count / axis_max
        bg_y = _pdf_top_to_bottom(height, y, bar_height)
        commands.append(
            _pdf_text_command(
                page_height=height,
                x=16,
                y_top=y + 7.0,
                text=event.label,
                font="F2",
                size=CHART_FONT_SIZE_PT,
                color=TEXT_RGB,
            )
        )
        commands.append(
            _pdf_text_command(
                page_height=height,
                x=16,
                y_top=y + 20.0,
                text=event.unit,
                font="F1",
                size=CHART_FONT_SIZE_PT,
                color=MUTED_RGB,
            )
        )
        commands.append(
            f"0.973 0.980 0.988 rg {left:.2f} {bg_y:.2f} {bar_width:.2f} {bar_height:.2f} re f"
        )
        commands.append(
            f"{GRID_RGB[0]:.3f} {GRID_RGB[1]:.3f} {GRID_RGB[2]:.3f} RG 0.7 w {left:.2f} {bg_y:.2f} {bar_width:.2f} {bar_height:.2f} re S"
        )
        if event.count > 0:
            red, green, blue = _hex_to_rgb(event.color)
            if bar_len >= 1.2:
                commands.append(
                    f"{red:.3f} {green:.3f} {blue:.3f} rg {left:.2f} {bg_y:.2f} {bar_len:.2f} {bar_height:.2f} re f"
                )
            else:
                marker_x = left + bar_len
                commands.append(
                    f"{red:.3f} {green:.3f} {blue:.3f} RG 1.2 w {marker_x:.2f} {bg_y:.2f} m {marker_x:.2f} {bg_y + bar_height:.2f} l S"
                )
        commands.append(
            _pdf_text_command(
                page_height=height,
                x=left + bar_width + 8.0,
                y_top=y + 11.0,
                text=f"{event.count:,}",
                font="F1",
                size=CHART_FONT_SIZE_PT,
                color=MUTED_RGB,
            )
        )
    return write_single_page_pdf(
        width=width, height=height, commands=commands, output_path=output_path
    )


def update_plot_manifest(
    manifest_path: Path,
    *,
    artifacts: PipelinePlotArtifacts,
    metrics: PipelinePlotMetrics,
    gene_sets: tuple[PipelineGeneSet, ...],
    coordinate_duplicate_summary: dict[str, object],
    annotation_consistency: AnnotationConsistencySummary | None = None,
) -> None:
    if manifest_path.exists():
        payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    else:
        payload = {}
    existing_plots = payload.get("plots")
    plots = dict(existing_plots) if isinstance(existing_plots, dict) else {}
    pipeline_payload: dict[str, object] = {
        "enabled": True,
        "stage_wise": True,
        "summary_json": str(artifacts.summary_json),
        "summary_tsv": str(artifacts.summary_tsv),
        "gene_flow_tsv": str(artifacts.gene_flow_tsv),
        "gene_flow_svg": str(artifacts.gene_flow_svg),
        "gene_flow_pdf": str(artifacts.gene_flow_pdf),
        "event_counts_tsv": str(artifacts.event_counts_tsv),
        "event_counts_svg": str(artifacts.event_counts_svg),
        "event_counts_pdf": str(artifacts.event_counts_pdf),
        "coordinate_duplicates": {
            "tsv": str(artifacts.coordinate_duplicates_tsv),
            "svg": str(artifacts.coordinate_duplicates_svg),
            "pdf": str(artifacts.coordinate_duplicates_pdf),
            **coordinate_duplicate_summary,
        },
        "metrics": metrics.to_dict(),
        "sources": metrics.sources,
        "gene_sets": {gene_set.key: gene_set.to_dict() for gene_set in gene_sets},
    }
    if annotation_consistency is not None:
        pipeline_payload["annotation_consistency"] = {
            "name_consistency_tsv": str(artifacts.name_consistency_tsv),
            "name_consistency_svg": str(artifacts.name_consistency_svg),
            "name_consistency_pdf": str(artifacts.name_consistency_pdf),
            "source_consistency_tsv": str(artifacts.source_consistency_tsv),
            "source_consistency_svg": str(artifacts.source_consistency_svg),
            "source_consistency_pdf": str(artifacts.source_consistency_pdf),
            "audit_tsv": str(annotation_consistency.path),
            "summary_tsv": str(annotation_consistency.summary_path),
            "source_pair_tsv": str(annotation_consistency.source_pair_path),
        }
    plots["pipeline"] = pipeline_payload
    payload["plots"] = plots
    write_text(manifest_path, json.dumps(payload, indent=2, sort_keys=True) + "\n")


def summarize_pipeline_plots(artifacts: PipelinePlotArtifacts) -> list[str]:
    payload = json.loads(artifacts.summary_json.read_text(encoding="utf-8"))
    metrics = payload["metrics"]
    return [
        (
            "pipeline-flow"
            f"\tgenes={metrics['initial_genes']}"
            f"\tdedup_removed={metrics['duplicate_removed_genes']}"
            f"\ttranscript_changed={metrics['transcript_changed_genes']}"
            f"\tremoved_mrnas={metrics['removed_mrnas']}"
            f"\tinframe_updated={metrics['inframe_updated_genes']}"
            f"\tpadding_updated={metrics['padding_updated_genes']}"
            f"\tgenes_with_stops={metrics['genes_with_stops']}"
            f"\tmisc_feature_genes={metrics['converted_to_misc_genes']}"
        )
    ]
