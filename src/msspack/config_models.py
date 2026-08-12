from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path

from .database_lock import DatabaseLockSettings
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
    tissue_type: str = ""
    isolation_source: str = ""


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
    coordinate_duplicate_policy: str = "longest_valid_cds"
    replace_product_patterns: list[str] = field(default_factory=list)
    replace_product_with: str = "hypothetical protein"
    validate_with_parser: bool = True
    validate_with_transchecker: bool = True
    validate_in_parallel: bool = True


@dataclass
class PlotsConfig:
    coordinate_duplicate_limit: int = 50


@dataclass
class ToolsConfig:
    cache_dir: str | None = None
    java: str = "java"
    gff3sort: str | None = None
    java_heap: str = "16G"


@dataclass
class DatabasesConfig:
    root: str = "msspack_db"
    lock_poll_seconds: float = 5.0
    lock_timeout_seconds: float = 86_400.0
    lock_heartbeat_seconds: float = 60.0
    lock_stale_seconds: float = 900.0


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
class FunctionalAnnotationConsistencyConfig:
    enabled: bool = False
    harmonize_safe_equivalents: bool = False
    auto_resolve_conflicts: bool = True
    near_identical_identity: float = 90.0
    near_identical_coverage: float = 90.0
    family_identity: float = 70.0
    family_coverage: float = 80.0
    broad_identity: float = 40.0
    broad_coverage: float = 60.0
    evalue: float = 1e-10
    name_similarity_threshold: float = 0.45
    source_pair_min_pairs: int = 5


@dataclass
class FunctionalAnnotationTaxonomyConfig:
    enabled: bool = True
    target_taxon_id: int = 0
    resolve_scientific_name: bool = True
    offline: bool = False
    busco_crosscheck: bool = True
    strict: bool = False
    distant_specificity_identity: float = 50.0


@dataclass
class FunctionalAnnotationConfig:
    enabled: bool = False
    diamond_command: str = "diamond"
    hmmscan_command: str = "hmmscan"
    hmmpress_command: str = "hmmpress"
    rpsblast_command: str = "rpsblast"
    rpsbproc_command: str = "rpsbproc"
    threads: int = 8
    sensitivity: str = "sensitive"
    evalue: float = 1e-10
    max_target_seqs: int = 25
    min_bitscore: float = 50.0
    min_identity: float = 35.0
    min_query_coverage: float = 70.0
    min_subject_coverage: float = 70.0
    near_top_bitscore_ratio: float = 0.90
    min_token_score: float = 0.60
    overwrite_existing: bool = False
    swissprot_enabled: bool = True
    swissprot_fasta: str = ""
    swissprot_url: str = (
        "https://ftp.uniprot.org/pub/databases/uniprot/current_release/"
        "knowledgebase/complete/uniprot_sprot.fasta.gz"
    )
    swissprot_weight: float = 1.0
    uniref90_enabled: bool = False
    uniref90_fasta: str = ""
    uniref90_url: str = (
        "https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/"
        "uniref90.fasta.gz"
    )
    uniref90_taxon_id: int = 0
    uniref90_weight: float = 0.5
    reference_proteins: str = ""
    reference_name: str = "reference"
    reference_weight: float = 1.2
    pfam_enabled: bool = True
    pfam_hmm: str = ""
    pfam_url: str = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"
    pfam_max_i_evalue: float = 1e-5
    pfam_min_domain_coverage: float = 0.35
    cdd_enabled: bool = False
    cdd_database: str = ""
    cdd_url: str = "https://ftp.ncbi.nlm.nih.gov/pub/mmdb/cdd/little_endian/Cdd_LE.tar.gz"
    cdd_data_dir: str = ""
    cdd_data_url: str = "https://ftp.ncbi.nlm.nih.gov/pub/mmdb/cdd"
    cdd_evalue: float = 0.01
    taxonomy: FunctionalAnnotationTaxonomyConfig = field(
        default_factory=FunctionalAnnotationTaxonomyConfig
    )
    consistency: FunctionalAnnotationConsistencyConfig = field(
        default_factory=FunctionalAnnotationConsistencyConfig
    )


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
    plots: PlotsConfig = field(default_factory=PlotsConfig)
    databases: DatabasesConfig = field(default_factory=DatabasesConfig)
    busco: BuscoConfig = field(default_factory=BuscoConfig)
    functional_annotation: FunctionalAnnotationConfig = field(
        default_factory=FunctionalAnnotationConfig
    )

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
    def database_dir(self) -> Path:
        override = os.environ.get("MSSPACK_DB_DIR", "").strip()
        if override:
            return expand_path(override, self.base_dir)
        return expand_path(self.databases.root, self.base_dir)

    @property
    def database_lock_settings(self) -> DatabaseLockSettings:
        return DatabaseLockSettings(
            poll_seconds=self.databases.lock_poll_seconds,
            timeout_seconds=self.databases.lock_timeout_seconds,
            heartbeat_seconds=self.databases.lock_heartbeat_seconds,
            stale_seconds=self.databases.lock_stale_seconds,
        )

    @property
    def busco_database_dir(self) -> Path:
        if self.busco.download_path:
            return expand_path(self.busco.download_path, self.base_dir)
        return self.database_dir / "busco"

    @property
    def fasta_path(self) -> Path:
        return expand_path(self.inputs.fasta, self.base_dir)

    @property
    def gff_path(self) -> Path:
        return expand_path(self.inputs.gff, self.base_dir)
