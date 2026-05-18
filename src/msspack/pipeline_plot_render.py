from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path
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
from .chart_primitives import pdf_text_command as _pdf_text_command
from .chart_primitives import pdf_top_to_bottom as _pdf_top_to_bottom
from .pipeline_plot_models import (
    SANKEY_COLORS,
    EventCount,
    GeneOverlapRow,
    PipelineGeneSet,
    PipelinePlotArtifacts,
    PipelinePlotMetrics,
    SankeyLink,
    SankeyNode,
)
from .utils import ensure_dir, write_text

PDF_POINTS_PER_INCH = 72.0
SANKEY_WIDTH_IN = 7.2
SANKEY_WIDTH = SANKEY_WIDTH_IN * PDF_POINTS_PER_INCH
SANKEY_HEIGHT = 288.0
SANKEY_STAGE_LABELS = {
    "Duplicate removal": "Dedup",
    "Transcript choice": "Transcript",
    "After transcript choice": "Selected",
    "Frame correction": "Frame",
    "After frame correction": "Framed",
    "Padding analysis": "Padding",
    "After padding": "Padded",
    "Final feature fate": "Final",
}
SANKEY_NODE_LABELS = {
    "start": "Input",
    "after_dedup": "Kept",
    "duplicate_removed": "Removed",
    "transcript_changed": "Pruned",
    "transcript_unchanged": "Unchanged",
    "after_transcript": "Selected",
    "inframe_updated": "Frame fixed",
    "inframe_unchanged": "Frame OK",
    "after_inframe": "Framed",
    "padding_updated": "Padded",
    "genes_with_stops": "Stops",
    "padding_unchanged": "No padding",
    "after_padding": "Padded",
    "final_cds": "Final CDS",
    "final_misc": "Final misc",
}
SANKEY_UNLABELED_NODES = {"after_transcript", "after_inframe", "after_padding"}


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
        overlap_tsv=root / "pipeline-gene-overlap.tsv",
        overlap_svg=root / "pipeline-gene-overlap.svg",
        overlap_pdf=root / "pipeline-gene-overlap.pdf",
    )


def _summary_rows(metrics: PipelinePlotMetrics) -> list[tuple[str, int, str, str]]:
    return [
        ("initial_genes", metrics.initial_genes, "genes", "Input gene models before duplicate removal."),
        (
            "duplicate_removed_genes",
            metrics.duplicate_removed_genes,
            "genes",
            "Genes removed because they shared identical coordinates.",
        ),
        ("genes_after_dedup", metrics.genes_after_dedup, "genes", "Genes retained after duplicate removal."),
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
        ("removed_mrnas", metrics.removed_mrnas, "transcripts", "Transcript models removed during single-mRNA selection."),
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
        ("genes_after_inframe", metrics.genes_after_inframe, "genes", "Genes remaining after frame correction."),
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
        ("genes_after_padding", metrics.genes_after_padding, "genes", "Genes carried into the MSS rendering stage."),
        (
            "converted_to_misc_genes",
            metrics.converted_to_misc_genes,
            "genes",
            "Genes whose CDS blocks were rewritten to misc_feature in the final MSS output.",
        ),
        ("final_cds_genes", metrics.final_cds_genes, "genes", "Genes that remain as CDS blocks in the final MSS output."),
        ("total_cds_input", metrics.total_cds_input, "features", "CDS feature count before CDS-to-misc conversion."),
        ("total_cds_output", metrics.total_cds_output, "features", "CDS feature count after CDS-to-misc conversion."),
        ("misc_feature_output", metrics.misc_feature_output, "features", "misc_feature count after CDS-to-misc conversion."),
    ]


def build_sankey(metrics: PipelinePlotMetrics) -> tuple[list[str], list[SankeyNode], list[SankeyLink]]:
    stage_labels = [
        "Input",
        "Duplicate removal",
        "Transcript choice",
        "After transcript choice",
        "Frame correction",
        "After frame correction",
        "Padding analysis",
        "After padding",
        "Final feature fate",
    ]
    nodes: list[SankeyNode] = [
        SankeyNode("start", "Input genes", 0, metrics.initial_genes, SANKEY_COLORS["start"]),
        SankeyNode("after_dedup", "After duplicate removal", 1, metrics.genes_after_dedup, SANKEY_COLORS["kept"]),
        SankeyNode("after_transcript", "After single-mRNA selection", 3, metrics.genes_after_single_mrna, SANKEY_COLORS["merge"]),
        SankeyNode("after_inframe", "After frame correction", 5, metrics.genes_after_inframe, SANKEY_COLORS["merge"]),
        SankeyNode("after_padding", "After padding analysis", 7, metrics.genes_after_padding, SANKEY_COLORS["merge"]),
        SankeyNode("final_cds", "Final CDS genes", 8, metrics.final_cds_genes, SANKEY_COLORS["final_cds"]),
    ]
    if metrics.duplicate_removed_genes > 0:
        nodes.append(SankeyNode("duplicate_removed", "Duplicate genes removed", 1, metrics.duplicate_removed_genes, SANKEY_COLORS["removed"]))
    if metrics.transcript_changed_genes > 0:
        nodes.append(SankeyNode("transcript_changed", "Genes with transcript pruning", 2, metrics.transcript_changed_genes, SANKEY_COLORS["transcript_changed"]))
    if metrics.transcript_unchanged_genes > 0:
        nodes.append(SankeyNode("transcript_unchanged", "Genes unchanged by transcript choice", 2, metrics.transcript_unchanged_genes, SANKEY_COLORS["transcript_unchanged"]))
    if metrics.inframe_updated_genes > 0:
        nodes.append(SankeyNode("inframe_updated", "Genes updated to restore frame", 4, metrics.inframe_updated_genes, SANKEY_COLORS["inframe_updated"]))
    if metrics.inframe_unchanged_genes > 0:
        nodes.append(SankeyNode("inframe_unchanged", "Genes unchanged by frame correction", 4, metrics.inframe_unchanged_genes, SANKEY_COLORS["inframe_unchanged"]))
    if metrics.padding_updated_genes > 0:
        nodes.append(SankeyNode("padding_updated", "Genes updated by padding", 6, metrics.padding_updated_genes, SANKEY_COLORS["padding_updated"]))
    if metrics.genes_with_stops > 0:
        nodes.append(SankeyNode("genes_with_stops", "Genes with stops after padding", 6, metrics.genes_with_stops, SANKEY_COLORS["genes_with_stops"]))
    if metrics.padding_unchanged_genes > 0:
        nodes.append(SankeyNode("padding_unchanged", "Genes unchanged by padding", 6, metrics.padding_unchanged_genes, SANKEY_COLORS["padding_unchanged"]))
    if metrics.converted_to_misc_genes > 0:
        nodes.append(SankeyNode("final_misc", "Final misc_feature genes", 8, metrics.converted_to_misc_genes, SANKEY_COLORS["final_misc"]))

    links: list[SankeyLink] = [SankeyLink("start", "after_dedup", metrics.genes_after_dedup, SANKEY_COLORS["kept"])]
    if metrics.duplicate_removed_genes > 0:
        links.append(SankeyLink("start", "duplicate_removed", metrics.duplicate_removed_genes, SANKEY_COLORS["removed"]))
    if metrics.transcript_changed_genes > 0:
        links.extend(
            [
                SankeyLink("after_dedup", "transcript_changed", metrics.transcript_changed_genes, SANKEY_COLORS["transcript_changed"]),
                SankeyLink("transcript_changed", "after_transcript", metrics.transcript_changed_genes, SANKEY_COLORS["transcript_changed"]),
            ]
        )
    if metrics.transcript_unchanged_genes > 0:
        links.extend(
            [
                SankeyLink("after_dedup", "transcript_unchanged", metrics.transcript_unchanged_genes, SANKEY_COLORS["transcript_unchanged"]),
                SankeyLink("transcript_unchanged", "after_transcript", metrics.transcript_unchanged_genes, SANKEY_COLORS["transcript_unchanged"]),
            ]
        )
    if metrics.inframe_updated_genes > 0:
        links.extend(
            [
                SankeyLink("after_transcript", "inframe_updated", metrics.inframe_updated_genes, SANKEY_COLORS["inframe_updated"]),
                SankeyLink("inframe_updated", "after_inframe", metrics.inframe_updated_genes, SANKEY_COLORS["inframe_updated"]),
            ]
        )
    if metrics.inframe_unchanged_genes > 0:
        links.extend(
            [
                SankeyLink("after_transcript", "inframe_unchanged", metrics.inframe_unchanged_genes, SANKEY_COLORS["inframe_unchanged"]),
                SankeyLink("inframe_unchanged", "after_inframe", metrics.inframe_unchanged_genes, SANKEY_COLORS["inframe_unchanged"]),
            ]
        )
    if metrics.padding_updated_genes > 0:
        links.extend(
            [
                SankeyLink("after_inframe", "padding_updated", metrics.padding_updated_genes, SANKEY_COLORS["padding_updated"]),
                SankeyLink("padding_updated", "after_padding", metrics.padding_updated_genes, SANKEY_COLORS["padding_updated"]),
            ]
        )
    if metrics.genes_with_stops > 0:
        links.extend(
            [
                SankeyLink("after_inframe", "genes_with_stops", metrics.genes_with_stops, SANKEY_COLORS["genes_with_stops"]),
                SankeyLink("genes_with_stops", "after_padding", metrics.genes_with_stops, SANKEY_COLORS["genes_with_stops"]),
            ]
        )
    if metrics.padding_unchanged_genes > 0:
        links.extend(
            [
                SankeyLink("after_inframe", "padding_unchanged", metrics.padding_unchanged_genes, SANKEY_COLORS["padding_unchanged"]),
                SankeyLink("padding_unchanged", "after_padding", metrics.padding_unchanged_genes, SANKEY_COLORS["padding_unchanged"]),
            ]
        )
    links.append(SankeyLink("after_padding", "final_cds", metrics.final_cds_genes, SANKEY_COLORS["final_cds"]))
    if metrics.converted_to_misc_genes > 0:
        links.append(SankeyLink("after_padding", "final_misc", metrics.converted_to_misc_genes, SANKEY_COLORS["final_misc"]))
    return stage_labels, sorted(nodes, key=lambda node: (node.stage, node.label)), links


def build_event_counts(metrics: PipelinePlotMetrics) -> list[EventCount]:
    return [
        EventCount("duplicate_removed_genes", "Duplicate genes removed", metrics.duplicate_removed_genes, "genes", SANKEY_COLORS["removed"]),
        EventCount("transcript_changed_genes", "Genes with transcript pruning", metrics.transcript_changed_genes, "genes", SANKEY_COLORS["transcript_changed"]),
        EventCount("removed_mrnas", "Removed mRNAs", metrics.removed_mrnas, "transcripts", "#7c2d12"),
        EventCount("inframe_updated_genes", "Genes updated to restore frame", metrics.inframe_updated_genes, "genes", SANKEY_COLORS["inframe_updated"]),
        EventCount("padding_updated_genes", "Genes updated by padding", metrics.padding_updated_genes, "genes", SANKEY_COLORS["padding_updated"]),
        EventCount("genes_with_stops", "Genes with stops after padding", metrics.genes_with_stops, "genes", SANKEY_COLORS["genes_with_stops"]),
        EventCount("converted_to_misc_genes", "Genes converted to misc_feature", metrics.converted_to_misc_genes, "genes", SANKEY_COLORS["final_misc"]),
    ]


def write_summary_json(payload: dict[str, object], output_path: Path) -> Path:
    return write_text(output_path, json.dumps(payload, indent=2, sort_keys=True) + "\n")


def write_summary_tsv(metrics: PipelinePlotMetrics, output_path: Path) -> Path:
    lines = ["metric\tvalue\tunit\tdescription"]
    for metric, value, unit, description in _summary_rows(metrics):
        lines.append(f"{metric}\t{value}\t{unit}\t{description}")
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


def write_overlap_tsv(
    gene_sets: tuple[PipelineGeneSet, ...],
    overlap_rows: tuple[GeneOverlapRow, ...],
    output_path: Path,
) -> Path:
    lines = ["rank\texclusive_gene_count\tmember_keys\tmember_labels"]
    for index, row in enumerate(overlap_rows, start=1):
        lines.append(
            "\t".join(
                [
                    str(index),
                    str(row.count),
                    ",".join(row.member_keys),
                    " + ".join(row.member_labels),
                ]
            )
        )
    if not overlap_rows and gene_sets:
        lines.append("0\t0\t\t")
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


def _sankey_node_is_labeled(node: SankeyNode) -> bool:
    return node.id not in SANKEY_UNLABELED_NODES


def _format_gene_count(count: int) -> str:
    suffix = "gene" if count == 1 else "genes"
    return f"{count:,} {suffix}"


def _sankey_layout(
    stage_labels: list[str],
    nodes: list[SankeyNode],
    links: list[SankeyLink],
) -> tuple[dict[str, _LaidOutNode], list[_LaidOutLink], dict[str, float]]:
    width = SANKEY_WIDTH
    chart_height = SANKEY_HEIGHT
    left = 28.0
    top = 78.0
    bottom = 24.0
    node_width = 9.0
    stage_gap = (width - left * 2 - node_width) / max(1, len(stage_labels) - 1)
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
    max_stage_total = max(totals_by_stage.values()) if totals_by_stage else 1
    max_gap_total = max(
        ((len(stage_nodes) - 1) * node_gap for stage_nodes in nodes_by_stage.values()),
        default=0.0,
    )
    scale = max(0.1, (content_height - max_gap_total) / max_stage_total)

    display_flow_heights: dict[tuple[str, str], float] = {}
    for link in links:
        display_flow_heights[(link.source, link.target)] = (
            max(link.count * scale, flow_min_height) if link.count > 0 else 0.0
        )

    outgoing_totals = {node.id: 0.0 for node in nodes}
    incoming_totals = {node.id: 0.0 for node in nodes}
    for link in links:
        display_height = display_flow_heights[(link.source, link.target)]
        outgoing_totals[link.source] += display_height
        incoming_totals[link.target] += display_height

    laid_out_nodes: dict[str, _LaidOutNode] = {}
    for stage_index, stage_nodes in nodes_by_stage.items():
        ordered_nodes = sorted(stage_nodes, key=lambda node: node.label)
        heights = [
            max(node.count * scale, outgoing_totals[node.id], incoming_totals[node.id], node_min_height)
            for node in ordered_nodes
        ]
        total_height = sum(heights) + max(0, len(heights) - 1) * node_gap
        cursor = top + (content_height - total_height) / 2.0
        x = left + stage_index * stage_gap
        for node, node_height in zip(ordered_nodes, heights):
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
        key=lambda link: (laid_out_nodes[link.source].node.stage, node_order[link.source], node_order[link.target]),
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
        laid_out_links.append(_LaidOutLink(link=link, source=source, target=target, source_y=source_y, target_y=target_y, height=display_height))

    return laid_out_nodes, laid_out_links, {
        "width": width,
        "height": chart_height,
        "left": left,
        "top": top,
        "stage_gap": stage_gap,
        "node_width": node_width,
    }


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
    if node.node.stage >= total_stages - 2:
        return node.x - 4.0, "end"
    return node.x + node.width + 4.0, "start"


def write_sankey_svg(
    stage_labels: list[str],
    nodes: list[SankeyNode],
    links: list[SankeyLink],
    output_path: Path,
) -> Path:
    laid_out_nodes, laid_out_links, meta = _sankey_layout(stage_labels, nodes, links)
    width = float(meta["width"])
    height = float(meta["height"])
    svg_width = f"{width / PDF_POINTS_PER_INCH:g}in"
    svg_height = f"{height / PDF_POINTS_PER_INCH:g}in"
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{svg_width}" height="{svg_height}" viewBox="0 0 {width:.2f} {height:.2f}">',
        f'<style>text{{font-family:Helvetica,Arial,sans-serif;fill:#111827}} .title{{font-size:{SVG_FONT_SIZE};font-weight:700}} .subtitle{{font-size:{SVG_FONT_SIZE};fill:#4b5563}} .stage{{font-size:{SVG_FONT_SIZE};font-weight:700;fill:#334155}} .label{{font-size:{SVG_FONT_SIZE}}} .count{{font-size:{SVG_FONT_SIZE};fill:#475569}}</style>',
        '<rect width="100%" height="100%" fill="white"/>',
        '<text x="16" y="16" class="title">Stage-wise pipeline gene flow</text>',
        '<text x="16" y="31" class="subtitle">Synthetic example; ribbons summarize stage counts.</text>',
    ]
    total_stages = len(stage_labels)
    for index, label in enumerate(stage_labels):
        stage_x = meta["left"] + index * meta["stage_gap"] + meta["node_width"] / 2.0
        parts.append(
            f'<text x="{stage_x:.2f}" y="58" text-anchor="middle" class="stage">{escape(_sankey_stage_label(label))}</text>'
        )
    for link in laid_out_links:
        parts.append(f'<path d="{_sankey_link_path(link)}" fill="{link.link.color}" fill-opacity="0.55" stroke="none"/>')
    for node in laid_out_nodes.values():
        parts.append(f'<rect x="{node.x:.2f}" y="{node.y:.2f}" width="{node.width:.2f}" height="{node.height:.2f}" rx="4" fill="{node.node.color}" stroke="#0f172a" stroke-width="0.8"/>')
    for node in laid_out_nodes.values():
        if not _sankey_node_is_labeled(node.node):
            continue
        label_x, anchor = _sankey_label_anchor(node, total_stages)
        label_y = node.y + min(max(node.height / 2.0, 9.0), node.height - 2.0)
        parts.append(
            f'<text x="{label_x:.2f}" y="{label_y - 4:.2f}" text-anchor="{anchor}" class="label">{escape(_sankey_node_label(node.node))}</text>'
        )
        parts.append(
            f'<text x="{label_x:.2f}" y="{label_y + 8:.2f}" text-anchor="{anchor}" class="count">{_format_gene_count(node.node.count)}</text>'
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


def write_sankey_pdf(
    stage_labels: list[str],
    nodes: list[SankeyNode],
    links: list[SankeyLink],
    output_path: Path,
) -> Path:
    laid_out_nodes, laid_out_links, meta = _sankey_layout(stage_labels, nodes, links)
    width = float(meta["width"])
    height = float(meta["height"])
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
            text="Synthetic example; ribbons summarize stage counts.",
            font="F1",
            size=CHART_FONT_SIZE_PT,
            color=MUTED_RGB,
        ),
    ]
    total_stages = len(stage_labels)
    for index, label in enumerate(stage_labels):
        stage_label = _sankey_stage_label(label)
        stage_x = meta["left"] + index * meta["stage_gap"] + meta["node_width"] / 2.0
        commands.append(
            _pdf_text_command(
                page_height=height,
                x=stage_x - 2.1 * len(stage_label),
                y_top=58,
                text=stage_label,
                font="F2",
                size=CHART_FONT_SIZE_PT,
                color=MUTED_RGB,
            )
        )
    for link in laid_out_links:
        r, g, b = _hex_to_rgb(link.link.color)
        commands.append(f"{r:.3f} {g:.3f} {b:.3f} rg")
        commands.append(_sankey_pdf_path_commands(link, height))
    for node in laid_out_nodes.values():
        r, g, b = _hex_to_rgb(node.node.color)
        rect_y = _pdf_top_to_bottom(height, node.y, node.height)
        commands.append(f"{r:.3f} {g:.3f} {b:.3f} rg {node.x:.2f} {rect_y:.2f} {node.width:.2f} {node.height:.2f} re f")
        commands.append(f"{TEXT_RGB[0]:.3f} {TEXT_RGB[1]:.3f} {TEXT_RGB[2]:.3f} RG 0.8 w {node.x:.2f} {rect_y:.2f} {node.width:.2f} {node.height:.2f} re S")
    for node in laid_out_nodes.values():
        if not _sankey_node_is_labeled(node.node):
            continue
        label_x, anchor = _sankey_label_anchor(node, total_stages)
        label_y = node.y + min(max(node.height / 2.0, 9.0), node.height - 2.0)
        label = _sankey_node_label(node.node)
        count = _format_gene_count(node.node.count)
        if anchor == "end":
            label_x = max(10.0, label_x - 3.4 * len(label))
        commands.append(
            _pdf_text_command(
                page_height=height,
                x=label_x,
                y_top=label_y - 4.0,
                text=label,
                font="F1",
                size=CHART_FONT_SIZE_PT,
                color=TEXT_RGB,
            )
        )
        count_x = label_x
        if anchor == "end":
            count_x = max(10.0, label_x - 3.4 * max(0, len(count) - len(label)))
        commands.append(
            _pdf_text_command(
                page_height=height,
                x=count_x,
                y_top=label_y + 8.0,
                text=count,
                font="F1",
                size=CHART_FONT_SIZE_PT,
                color=MUTED_RGB,
            )
        )
    return write_single_page_pdf(width=width, height=height, commands=commands, output_path=output_path)


def write_event_counts_svg(events: list[EventCount], output_path: Path) -> Path:
    width = 1200
    height = 140 + len(events) * 72
    left = 360
    top = 110
    bar_height = 26
    bar_gap = 72
    bar_width = 740
    max_count = max((event.count for event in events), default=1)
    axis_max = max(1, int(math.ceil(max_count / 5.0) * 5))
    ticks = [0, axis_max / 4.0, axis_max / 2.0, axis_max * 3.0 / 4.0, axis_max]
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        f'<style>text{{font-family:Helvetica,Arial,sans-serif;fill:#111827}} .title{{font-size:{SVG_FONT_SIZE};font-weight:700}} .subtitle{{font-size:{SVG_FONT_SIZE};fill:#4b5563}} .label{{font-size:{SVG_FONT_SIZE};font-weight:700}} .unit{{font-size:{SVG_FONT_SIZE};fill:#64748b}} .value{{font-size:{SVG_FONT_SIZE};fill:#334155}} .tick{{font-size:{SVG_FONT_SIZE};fill:#6b7280}}</style>',
        '<rect width="100%" height="100%" fill="white"/>',
        '<text x="40" y="40" class="title">Pipeline event counts</text>',
        '<text x="40" y="62" class="subtitle">Step-level counts from the packaging logs. Labels include the unit for each metric because removed mRNAs are transcript counts while the other bars are gene counts.</text>',
    ]
    for tick in ticks:
        x = left + (bar_width * tick / axis_max if axis_max else 0.0)
        parts.append(f'<line x1="{x:.2f}" y1="{top - 24}" x2="{x:.2f}" y2="{height - 32}" stroke="#e5e7eb"/>')
        parts.append(f'<text x="{x:.2f}" y="{top - 30}" text-anchor="middle" class="tick">{int(tick) if float(tick).is_integer() else f"{tick:.1f}"}</text>')
    for index, event in enumerate(events):
        y = top + index * bar_gap
        bar_len = 0.0 if axis_max == 0 else bar_width * event.count / axis_max
        color = event.color if event.count > 0 else "#cbd5e1"
        parts.append(f'<text x="40" y="{y + 12}" class="label">{escape(event.label)}</text>')
        parts.append(f'<text x="40" y="{y + 28}" class="unit">{escape(event.unit)}</text>')
        parts.append(f'<rect x="{left}" y="{y}" width="{bar_width}" height="{bar_height}" rx="4" fill="#f8fafc" stroke="#cbd5e1"/>')
        if event.count > 0:
            parts.append(f'<rect x="{left}" y="{y}" width="{bar_len:.2f}" height="{bar_height}" rx="4" fill="{color}"/>')
        parts.append(f'<text x="{left + bar_width + 16}" y="{y + 18}" class="value">{event.count:,}</text>')
    parts.append("</svg>")
    return write_text(output_path, "\n".join(parts) + "\n")


def write_event_counts_pdf(events: list[EventCount], output_path: Path) -> Path:
    width = 1200.0
    height = float(140 + len(events) * 72)
    left = 360.0
    top = 110.0
    bar_height = 26.0
    bar_gap = 72.0
    bar_width = 740.0
    max_count = max((event.count for event in events), default=1)
    axis_max = max(1, int(math.ceil(max_count / 5.0) * 5))
    ticks = [0, axis_max / 4.0, axis_max / 2.0, axis_max * 3.0 / 4.0, axis_max]
    commands = [
        f"1 1 1 rg 0 0 {width:.2f} {height:.2f} re f",
        _pdf_text_command(page_height=height, x=40, y_top=40, text="Pipeline event counts", font="F2", size=CHART_FONT_SIZE_PT, color=TEXT_RGB),
        _pdf_text_command(page_height=height, x=40, y_top=62, text="Step-level counts from packaging logs. Removed mRNAs are transcript counts; the other bars are gene counts.", font="F1", size=CHART_FONT_SIZE_PT, color=MUTED_RGB),
    ]
    for tick in ticks:
        x = left + (bar_width * tick / axis_max if axis_max else 0.0)
        y1 = _pdf_top_to_bottom(height, top - 24)
        y2 = _pdf_top_to_bottom(height, height - 32)
        commands.append(f"{GRID_RGB[0]:.3f} {GRID_RGB[1]:.3f} {GRID_RGB[2]:.3f} RG 1 w {x:.2f} {y1:.2f} m {x:.2f} {y2:.2f} l S")
        tick_text = str(int(tick) if float(tick).is_integer() else round(tick, 1))
        commands.append(_pdf_text_command(page_height=height, x=x - 8, y_top=top - 30, text=tick_text, font="F1", size=CHART_FONT_SIZE_PT, color=MUTED_RGB))
    for index, event in enumerate(events):
        y = top + index * bar_gap
        bar_len = 0.0 if axis_max == 0 else bar_width * event.count / axis_max
        bg_y = _pdf_top_to_bottom(height, y, bar_height)
        commands.append(_pdf_text_command(page_height=height, x=40, y_top=y + 12, text=event.label, font="F2", size=CHART_FONT_SIZE_PT, color=TEXT_RGB))
        commands.append(_pdf_text_command(page_height=height, x=40, y_top=y + 28, text=event.unit, font="F1", size=CHART_FONT_SIZE_PT, color=MUTED_RGB))
        commands.append(f"0.973 0.980 0.988 rg {left:.2f} {bg_y:.2f} {bar_width:.2f} {bar_height:.2f} re f")
        commands.append(f"{GRID_RGB[0]:.3f} {GRID_RGB[1]:.3f} {GRID_RGB[2]:.3f} RG {left:.2f} {bg_y:.2f} {bar_width:.2f} {bar_height:.2f} re S")
        if event.count > 0:
            r, g, b = _hex_to_rgb(event.color)
            commands.append(f"{r:.3f} {g:.3f} {b:.3f} rg {left:.2f} {bg_y:.2f} {bar_len:.2f} {bar_height:.2f} re f")
        commands.append(_pdf_text_command(page_height=height, x=left + bar_width + 16, y_top=y + 18, text=f"{event.count:,}", font="F1", size=CHART_FONT_SIZE_PT, color=MUTED_RGB))
    return write_single_page_pdf(width=width, height=height, commands=commands, output_path=output_path)


def _overlap_rows_for_plot(overlap_rows: tuple[GeneOverlapRow, ...]) -> list[GeneOverlapRow]:
    return list(overlap_rows[:12])


def write_overlap_svg(
    gene_sets: tuple[PipelineGeneSet, ...],
    overlap_rows: tuple[GeneOverlapRow, ...],
    output_path: Path,
) -> Path:
    plot_rows = _overlap_rows_for_plot(overlap_rows)
    width = 1280
    row_height = 38
    matrix_left = 420
    matrix_top = 170
    marker_size = 14
    column_gap = 82
    bar_left = matrix_left + max(1, len(gene_sets) - 1) * column_gap + 90
    bar_width = max(200, width - bar_left - 110)
    height = 220 + max(1, len(plot_rows)) * row_height
    max_count = max((row.count for row in plot_rows), default=1)
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        f'<style>text{{font-family:Helvetica,Arial,sans-serif;fill:#111827}} .title{{font-size:{SVG_FONT_SIZE};font-weight:700}} .subtitle{{font-size:{SVG_FONT_SIZE};fill:#4b5563}} .axis{{font-size:{SVG_FONT_SIZE};fill:#475569}} .count{{font-size:{SVG_FONT_SIZE};fill:#334155}} .label{{font-size:{SVG_FONT_SIZE};fill:#334155}}</style>',
        '<rect width="100%" height="100%" fill="white"/>',
        '<text x="40" y="40" class="title">Changed-gene overlap</text>',
        '<text x="40" y="62" class="subtitle">Exclusive intersections across stage-emitted changed-gene ID sets. Filled squares mark the stages included in each overlap row.</text>',
    ]
    for column, gene_set in enumerate(gene_sets):
        x = matrix_left + column * column_gap + marker_size / 2.0
        parts.append(f'<text x="{x:.2f}" y="128" text-anchor="start" transform="rotate(-35 {x:.2f} 128)" class="label">{escape(gene_set.label)}</text>')
    parts.append(f'<text x="{bar_left:.2f}" y="128" class="label">Exclusive genes</text>')
    for row_index, row in enumerate(plot_rows):
        y = matrix_top + row_index * row_height
        parts.append(f'<line x1="40" y1="{y + marker_size / 2.0:.2f}" x2="{width - 40}" y2="{y + marker_size / 2.0:.2f}" stroke="#f1f5f9"/>')
        included = set(row.member_keys)
        for column, gene_set in enumerate(gene_sets):
            x = matrix_left + column * column_gap
            fill = gene_set.color if gene_set.key in included else "#ffffff"
            stroke = gene_set.color if gene_set.key in included else "#cbd5e1"
            parts.append(f'<rect x="{x:.2f}" y="{y:.2f}" width="{marker_size}" height="{marker_size}" rx="3" fill="{fill}" stroke="{stroke}" stroke-width="1.2"/>')
        bar_len = 0.0 if max_count == 0 else bar_width * row.count / max_count
        parts.append(f'<rect x="{bar_left:.2f}" y="{y - 2:.2f}" width="{bar_width:.2f}" height="{marker_size + 4}" rx="4" fill="#f8fafc" stroke="#cbd5e1"/>')
        parts.append(f'<rect x="{bar_left:.2f}" y="{y - 2:.2f}" width="{bar_len:.2f}" height="{marker_size + 4}" rx="4" fill="#1d4ed8"/>')
        parts.append(f'<text x="{bar_left + bar_width + 16:.2f}" y="{y + 10:.2f}" class="count">{row.count:,}</text>')
    if not plot_rows:
        parts.append('<text x="40" y="170" class="label">No non-empty exclusive overlaps were available for plotting.</text>')
    parts.append("</svg>")
    return write_text(output_path, "\n".join(parts) + "\n")


def write_overlap_pdf(
    gene_sets: tuple[PipelineGeneSet, ...],
    overlap_rows: tuple[GeneOverlapRow, ...],
    output_path: Path,
) -> Path:
    plot_rows = _overlap_rows_for_plot(overlap_rows)
    width = 1280.0
    row_height = 38.0
    matrix_left = 420.0
    matrix_top = 170.0
    marker_size = 14.0
    column_gap = 82.0
    bar_left = matrix_left + max(1, len(gene_sets) - 1) * column_gap + 90.0
    bar_width = max(200.0, width - bar_left - 110.0)
    height = 220.0 + max(1, len(plot_rows)) * row_height
    max_count = max((row.count for row in plot_rows), default=1)
    commands = [
        f"1 1 1 rg 0 0 {width:.2f} {height:.2f} re f",
        _pdf_text_command(page_height=height, x=40, y_top=40, text="Changed-gene overlap", font="F2", size=CHART_FONT_SIZE_PT, color=TEXT_RGB),
        _pdf_text_command(page_height=height, x=40, y_top=62, text="Exclusive intersections across stage-emitted changed-gene ID sets.", font="F1", size=CHART_FONT_SIZE_PT, color=MUTED_RGB),
    ]
    for column, gene_set in enumerate(gene_sets):
        x = matrix_left + column * column_gap
        commands.append(_pdf_text_command(page_height=height, x=x, y_top=128, text=gene_set.label, font="F1", size=CHART_FONT_SIZE_PT, color=MUTED_RGB))
    commands.append(_pdf_text_command(page_height=height, x=bar_left, y_top=128, text="Exclusive genes", font="F1", size=CHART_FONT_SIZE_PT, color=MUTED_RGB))
    for row_index, row in enumerate(plot_rows):
        y = matrix_top + row_index * row_height
        line_y = _pdf_top_to_bottom(height, y + marker_size / 2.0)
        commands.append(f"{GRID_RGB[0]:.3f} {GRID_RGB[1]:.3f} {GRID_RGB[2]:.3f} RG 1 w 40 {line_y:.2f} m {width - 40:.2f} {line_y:.2f} l S")
        included = set(row.member_keys)
        for column, gene_set in enumerate(gene_sets):
            x = matrix_left + column * column_gap
            marker_y = _pdf_top_to_bottom(height, y, marker_size)
            if gene_set.key in included:
                r, g, b = _hex_to_rgb(gene_set.color)
                commands.append(f"{r:.3f} {g:.3f} {b:.3f} rg {x:.2f} {marker_y:.2f} {marker_size:.2f} {marker_size:.2f} re f")
            commands.append(f"{GRID_RGB[0]:.3f} {GRID_RGB[1]:.3f} {GRID_RGB[2]:.3f} RG {x:.2f} {marker_y:.2f} {marker_size:.2f} {marker_size:.2f} re S")
        bar_len = 0.0 if max_count == 0 else bar_width * row.count / max_count
        bg_y = _pdf_top_to_bottom(height, y - 2.0, marker_size + 4.0)
        commands.append(f"0.973 0.980 0.988 rg {bar_left:.2f} {bg_y:.2f} {bar_width:.2f} {marker_size + 4.0:.2f} re f")
        commands.append(f"{GRID_RGB[0]:.3f} {GRID_RGB[1]:.3f} {GRID_RGB[2]:.3f} RG {bar_left:.2f} {bg_y:.2f} {bar_width:.2f} {marker_size + 4.0:.2f} re S")
        commands.append(f"0.114 0.306 0.847 rg {bar_left:.2f} {bg_y:.2f} {bar_len:.2f} {marker_size + 4.0:.2f} re f")
        commands.append(_pdf_text_command(page_height=height, x=bar_left + bar_width + 16.0, y_top=y + 10.0, text=f"{row.count:,}", font="F1", size=CHART_FONT_SIZE_PT, color=MUTED_RGB))
    if not plot_rows:
        commands.append(_pdf_text_command(page_height=height, x=40, y_top=170, text="No non-empty exclusive overlaps were available for plotting.", font="F1", size=CHART_FONT_SIZE_PT, color=MUTED_RGB))
    return write_single_page_pdf(width=width, height=height, commands=commands, output_path=output_path)


def update_plot_manifest(
    manifest_path: Path,
    *,
    artifacts: PipelinePlotArtifacts,
    metrics: PipelinePlotMetrics,
    gene_sets: tuple[PipelineGeneSet, ...],
    overlap_rows: tuple[GeneOverlapRow, ...],
) -> None:
    if manifest_path.exists():
        payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    else:
        payload = {}
    existing_plots = payload.get("plots")
    plots = dict(existing_plots) if isinstance(existing_plots, dict) else {}
    plots["pipeline"] = {
        "enabled": True,
        "stage_wise": True,
        "gene_id_overlap": True,
        "summary_json": str(artifacts.summary_json),
        "summary_tsv": str(artifacts.summary_tsv),
        "gene_flow_tsv": str(artifacts.gene_flow_tsv),
        "gene_flow_svg": str(artifacts.gene_flow_svg),
        "gene_flow_pdf": str(artifacts.gene_flow_pdf),
        "event_counts_tsv": str(artifacts.event_counts_tsv),
        "event_counts_svg": str(artifacts.event_counts_svg),
        "event_counts_pdf": str(artifacts.event_counts_pdf),
        "overlap_tsv": str(artifacts.overlap_tsv),
        "overlap_svg": str(artifacts.overlap_svg),
        "overlap_pdf": str(artifacts.overlap_pdf),
        "metrics": metrics.to_dict(),
        "sources": metrics.sources,
        "gene_sets": {gene_set.key: gene_set.to_dict() for gene_set in gene_sets},
        "overlap": {"row_count": len(overlap_rows)},
    }
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
