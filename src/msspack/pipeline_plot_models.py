from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

SANKEY_COLORS = {
    "start": "#475569",
    "removed": "#dc2626",
    "kept": "#64748b",
    "transcript_changed": "#2563eb",
    "transcript_unchanged": "#94a3b8",
    "merge": "#475569",
    "inframe_updated": "#d97706",
    "inframe_unchanged": "#cbd5e1",
    "padding_updated": "#0891b2",
    "genes_with_stops": "#be123c",
    "padding_unchanged": "#cbd5e1",
    "final_cds": "#059669",
    "final_misc": "#b45309",
}

GENE_SET_SPECS = (
    (
        "duplicate_removed_genes",
        "Duplicate genes removed",
        SANKEY_COLORS["removed"],
        "06.drop-duplicate-coordinate-gene.changed-gene-ids.txt",
    ),
    (
        "transcript_changed_genes",
        "Transcript pruning",
        SANKEY_COLORS["transcript_changed"],
        "07.select-one-mrna.changed-gene-ids.txt",
    ),
    (
        "inframe_updated_genes",
        "Frame correction",
        SANKEY_COLORS["inframe_updated"],
        "09.update-gff-to-inframe.changed-gene-ids.txt",
    ),
    (
        "padding_updated_genes",
        "Padding updated",
        SANKEY_COLORS["padding_updated"],
        "11.gff.updated-genes.txt",
    ),
    (
        "genes_with_stops",
        "Genes with stops",
        SANKEY_COLORS["genes_with_stops"],
        "11.gff.genes-with-stops.txt",
    ),
    (
        "converted_to_misc_genes",
        "Converted to misc_feature",
        SANKEY_COLORS["final_misc"],
        "16.mss-cds-to-misc.changed-gene-ids.txt",
    ),
)


@dataclass(frozen=True)
class ParsedStepRecord:
    path: Path
    step: str
    count_unit: str
    input_total: int | None
    changed_total: int | None
    output_total: int | None
    details: dict[str, object]


@dataclass(frozen=True)
class PipelinePlotMetrics:
    initial_genes: int
    duplicate_removed_genes: int
    genes_after_dedup: int
    transcript_changed_genes: int
    transcript_unchanged_genes: int
    removed_mrnas: int
    genes_after_single_mrna: int
    inframe_updated_genes: int
    inframe_unchanged_genes: int
    genes_after_inframe: int
    padding_updated_genes: int
    genes_with_stops: int
    padding_unchanged_genes: int
    genes_after_padding: int
    converted_to_misc_genes: int
    final_cds_genes: int
    total_cds_input: int
    total_cds_output: int
    misc_feature_output: int
    sources: dict[str, str]

    def to_dict(self) -> dict[str, object]:
        return {
            "initial_genes": self.initial_genes,
            "duplicate_removed_genes": self.duplicate_removed_genes,
            "genes_after_dedup": self.genes_after_dedup,
            "transcript_changed_genes": self.transcript_changed_genes,
            "transcript_unchanged_genes": self.transcript_unchanged_genes,
            "removed_mrnas": self.removed_mrnas,
            "genes_after_single_mrna": self.genes_after_single_mrna,
            "inframe_updated_genes": self.inframe_updated_genes,
            "inframe_unchanged_genes": self.inframe_unchanged_genes,
            "genes_after_inframe": self.genes_after_inframe,
            "padding_updated_genes": self.padding_updated_genes,
            "genes_with_stops": self.genes_with_stops,
            "padding_unchanged_genes": self.padding_unchanged_genes,
            "genes_after_padding": self.genes_after_padding,
            "converted_to_misc_genes": self.converted_to_misc_genes,
            "final_cds_genes": self.final_cds_genes,
            "total_cds_input": self.total_cds_input,
            "total_cds_output": self.total_cds_output,
            "misc_feature_output": self.misc_feature_output,
        }


@dataclass(frozen=True)
class PipelineGeneSet:
    key: str
    label: str
    color: str
    path: Path
    gene_ids: tuple[str, ...]

    @property
    def count(self) -> int:
        return len(self.gene_ids)

    def to_dict(self) -> dict[str, object]:
        return {
            "label": self.label,
            "color": self.color,
            "count": self.count,
            "path": str(self.path),
        }


@dataclass(frozen=True)
class GeneOverlapRow:
    member_keys: tuple[str, ...]
    member_labels: tuple[str, ...]
    count: int
    gene_ids: tuple[str, ...]

    @property
    def key(self) -> str:
        return "+".join(self.member_keys)

    @property
    def label(self) -> str:
        return " + ".join(self.member_labels)

    def to_dict(self) -> dict[str, object]:
        return {
            "member_keys": list(self.member_keys),
            "member_labels": list(self.member_labels),
            "count": self.count,
        }


@dataclass(frozen=True)
class PipelinePlotArtifacts:
    root: Path
    summary_json: Path
    summary_tsv: Path
    gene_flow_tsv: Path
    gene_flow_svg: Path
    gene_flow_pdf: Path
    event_counts_tsv: Path
    event_counts_svg: Path
    event_counts_pdf: Path
    overlap_tsv: Path
    overlap_svg: Path
    overlap_pdf: Path


@dataclass(frozen=True)
class SankeyNode:
    id: str
    label: str
    stage: int
    count: int
    color: str


@dataclass(frozen=True)
class SankeyLink:
    source: str
    target: str
    count: int
    color: str


@dataclass(frozen=True)
class EventCount:
    key: str
    label: str
    count: int
    unit: str
    color: str


@dataclass(frozen=True)
class PipelinePlotDataBundle:
    records: dict[str, ParsedStepRecord]
    metrics: PipelinePlotMetrics
    gene_sets: tuple[PipelineGeneSet, ...]
    overlap_rows: tuple[GeneOverlapRow, ...]
    summary_payload: dict[str, object]
