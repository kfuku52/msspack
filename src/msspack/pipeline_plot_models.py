from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

SANKEY_COLORS = {
    "start": "#334155",
    "removed": "#e11d48",
    "kept": "#6366f1",
    "transcript_changed": "#2563eb",
    "transcript_unchanged": "#0ea5e9",
    "merge": "#475569",
    "inframe_updated": "#f59e0b",
    "inframe_unchanged": "#14b8a6",
    "padding_updated": "#8b5cf6",
    "genes_with_stops": "#e11d48",
    "padding_unchanged": "#22c55e",
    "final_cds": "#059669",
    "final_misc": "#ea580c",
    "annotation_similarity": "#2563eb",
    "annotation_uniref": "#0d9488",
    "annotation_domain": "#7c3aed",
    "annotation_cdd": "#c026d3",
    "annotation_existing": "#ea580c",
    "annotation_missing": "#94a3b8",
    "consistency_consistent": "#1d4ed8",
    "consistency_review": "#d97706",
    "consistency_no_peer": "#64748b",
    "consistency_unannotated": "#cbd5e1",
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
    coordinate_duplicates_tsv: Path
    coordinate_duplicates_svg: Path
    coordinate_duplicates_pdf: Path
    name_consistency_tsv: Path
    name_consistency_svg: Path
    name_consistency_pdf: Path
    source_consistency_tsv: Path
    source_consistency_svg: Path
    source_consistency_pdf: Path


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
class SankeyBuscoSummary:
    label: str
    stage: int
    lineage_dataset: str
    input_sequences: int
    total_buscos: int
    single_copy: int
    duplicated: int
    fragmented: int
    missing: int

    @property
    def complete(self) -> int:
        return self.single_copy + self.duplicated

    @property
    def complete_pct(self) -> float:
        return 100.0 * self.complete / self.total_buscos

    def segment_counts(self) -> tuple[tuple[str, int], ...]:
        return (
            ("single_copy", self.single_copy),
            ("duplicated", self.duplicated),
            ("fragmented", self.fragmented),
            ("missing", self.missing),
        )


@dataclass(frozen=True)
class EventCount:
    key: str
    label: str
    count: int
    unit: str
    color: str


@dataclass(frozen=True)
class FunctionalAnnotationGroup:
    key: str
    source: str
    label: str
    color: str
    locus_tags: tuple[str, ...]

    @property
    def count(self) -> int:
        return len(self.locus_tags)

    def to_dict(self) -> dict[str, object]:
        return {
            "source": self.source,
            "label": self.label,
            "color": self.color,
            "count": self.count,
        }


@dataclass(frozen=True)
class FunctionalAnnotationSummary:
    path: Path
    groups: tuple[FunctionalAnnotationGroup, ...]

    @property
    def total(self) -> int:
        return sum(group.count for group in self.groups)

    def to_dict(self) -> dict[str, object]:
        return {
            "path": str(self.path),
            "total": self.total,
            "groups": {group.key: group.to_dict() for group in self.groups},
        }


@dataclass(frozen=True)
class AnnotationConsistencyGroup:
    key: str
    label: str
    color: str
    locus_tags: tuple[str, ...]

    @property
    def count(self) -> int:
        return len(self.locus_tags)

    def to_dict(self) -> dict[str, object]:
        return {
            "label": self.label,
            "color": self.color,
            "count": self.count,
        }


@dataclass(frozen=True)
class AnnotationConsistencySummary:
    path: Path
    summary_path: Path
    source_pair_path: Path
    groups: tuple[AnnotationConsistencyGroup, ...]
    comparison_tier: str = "family"
    identity_threshold: float = 70.0
    coverage_threshold: float = 80.0

    @property
    def total(self) -> int:
        return sum(group.count for group in self.groups)

    def to_dict(self) -> dict[str, object]:
        return {
            "path": str(self.path),
            "summary_path": str(self.summary_path),
            "source_pair_path": str(self.source_pair_path),
            "comparison_tier": self.comparison_tier,
            "identity_threshold": self.identity_threshold,
            "coverage_threshold": self.coverage_threshold,
            "total": self.total,
            "groups": {group.key: group.to_dict() for group in self.groups},
        }


@dataclass(frozen=True)
class PipelinePlotDataBundle:
    records: dict[str, ParsedStepRecord]
    metrics: PipelinePlotMetrics
    gene_sets: tuple[PipelineGeneSet, ...]
    functional_annotation: FunctionalAnnotationSummary | None
    annotation_consistency: AnnotationConsistencySummary | None
    summary_payload: dict[str, object]
