from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from .utils import default_cache_dir, expand_path


@dataclass
class ProjectConfig:
    name: str
    output_dir: str | None = None


@dataclass
class InputsConfig:
    fasta: str
    gff: str


@dataclass
class SampleConfig:
    locus_tag: str
    locus_tag_digits: int
    scientific_name: str
    strain: str = ""
    isolate: str = ""
    geo_loc_name: str = ""
    collection_date: str = ""
    sex: str = ""
    linkage_evidence: str = "paired-ends"
    genetic_code: str = "1"


@dataclass
class SubmissionConfig:
    datatype: str
    hold_date: str
    bioproject: str
    biosample: str
    sra: list[str] = field(default_factory=list)
    keywords: list[str] = field(default_factory=lambda: ["WGS"])


@dataclass
class SubmitterConfig:
    ab_name: list[str]
    contact: str
    institute: str
    department: str
    country: str
    state: str
    city: str
    street: str
    zip: str
    phone: str
    email: str


@dataclass
class ReferenceConfig:
    title: str
    ab_name: list[str]
    year: int
    status: str = "Unpublished"


@dataclass
class StCommentConfig:
    tagset_id: str = "Genome-Assembly-Data"
    assembly_method: str = ""
    assembly_name: str = ""
    genome_coverage: str = ""
    sequencing_technology: str = ""


@dataclass
class PipelineConfig:
    run_gapjust: bool = True
    gapjust_gap_len: int = 100
    gapjust_min: int = 80
    gapjust_max: int = 120
    gap_assembly: str = "unknown"
    infer_complete: bool = True
    feature_with_gap: str = "misc_feature"
    min_assembly_gap: int = 50
    min_artificial_intron_size: int = 10
    replace_product_patterns: list[str] = field(default_factory=list)
    replace_product_with: str = "hypothetical protein"
    validate_with_parser: bool = True
    validate_with_transchecker: bool = True
    validate_in_parallel: bool = True


@dataclass
class ToolsConfig:
    cache_dir: str | None = None
    java: str = "java"
    gff3sort: str | None = None
    java_heap: str = "16G"


@dataclass
class BuscoConfig:
    command: str = "busco"
    run_cds: bool = True
    run_genome: bool = False
    mode: str = "genome"
    cds_mode: str = "transcriptome"
    lineage_dataset: str = ""
    auto_lineage: bool = True
    auto_lineage_scope: str = "all"
    threads: int = 1
    download_path: str = ""
    offline: bool = False
    force: bool = False
    opt_out_run_stats: bool = True


@dataclass
class MSSPackConfig:
    base_dir: Path
    project: ProjectConfig
    inputs: InputsConfig
    sample: SampleConfig
    submission: SubmissionConfig
    submitter: SubmitterConfig
    reference: ReferenceConfig
    st_comment: StCommentConfig
    pipeline: PipelineConfig
    tools: ToolsConfig
    busco: BuscoConfig = field(default_factory=BuscoConfig)

    @property
    def output_dir(self) -> Path:
        if self.project.output_dir:
            return expand_path(self.project.output_dir, self.base_dir)
        return (self.base_dir / "build" / self.project.name).resolve()

    @property
    def cache_dir(self) -> Path:
        if self.tools.cache_dir:
            return expand_path(self.tools.cache_dir, self.base_dir)
        return default_cache_dir()

    @property
    def fasta_path(self) -> Path:
        return expand_path(self.inputs.fasta, self.base_dir)

    @property
    def gff_path(self) -> Path:
        return expand_path(self.inputs.gff, self.base_dir)
