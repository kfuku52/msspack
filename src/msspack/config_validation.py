from __future__ import annotations

import math
import re
from datetime import datetime

from Bio.Data import CodonTable

from .config_errors import ConfigError
from .config_models import (
    BuscoConfig,
    DatabasesConfig,
    FunctionalAnnotationConfig,
    MSSPackConfig,
    PipelineConfig,
    PlotsConfig,
    ProjectConfig,
    ReferenceConfig,
    SampleConfig,
    SubmissionConfig,
    SubmitterConfig,
    ToolsConfig,
)
from .coordinate_duplicates import COORDINATE_DUPLICATE_POLICIES

HOLD_DATE_RE = re.compile(r"^\d{8}$")
JAVA_HEAP_RE = re.compile(r"^\d+[KMGkmg]$")
VALID_GAP_ASSEMBLY = {"known", "unknown"}
VALID_FEATURE_WITH_GAP = {"asis", "misc_feature"}
VALID_BUSCO_MODE = {"genome", "proteins", "transcriptome"}
VALID_BUSCO_AUTO_LINEAGE_SCOPE = {"all", "euk", "prok"}
VALID_DIAMOND_SENSITIVITY = {
    "default",
    "fast",
    "mid-sensitive",
    "sensitive",
    "more-sensitive",
    "very-sensitive",
    "ultra-sensitive",
}


def ensure_nonempty(value: str, key: str) -> None:
    if not value.strip():
        raise ConfigError(f"Config value '{key}' must not be empty")
    if any(character in value for character in ("\x00", "\r", "\n", "\t")):
        raise ConfigError(f"Config value '{key}' must not contain control characters")


def ensure_nonnegative(value: int, key: str) -> None:
    if value < 0:
        raise ConfigError(f"Config value '{key}' must be >= 0")


def ensure_positive(value: int, key: str) -> None:
    if value <= 0:
        raise ConfigError(f"Config value '{key}' must be > 0")


def ensure_date(value: str, key: str, fmt: str) -> None:
    if not value:
        return
    try:
        datetime.strptime(value, fmt)
    except ValueError as exc:
        raise ConfigError(f"Config value '{key}' must match {fmt}") from exc


def ensure_collection_date(value: str, key: str) -> None:
    if not value:
        return
    parts = value.split("/")
    if len(parts) not in {1, 2}:
        raise ConfigError(
            f"Config value '{key}' must match YYYY-MM-DD or YYYY-MM-DD/YYYY-MM-DD"
        )
    try:
        dates = [datetime.strptime(part, "%Y-%m-%d") for part in parts]
    except ValueError as exc:
        raise ConfigError(
            f"Config value '{key}' must match YYYY-MM-DD or YYYY-MM-DD/YYYY-MM-DD"
        ) from exc
    if len(dates) == 2 and dates[0] > dates[1]:
        raise ConfigError(
            f"Config value '{key}' range start must not be after its end"
        )


def ensure_choice(value: str, key: str, choices: set[str]) -> None:
    if value not in choices:
        raise ConfigError(f"Config value '{key}' must be one of: {', '.join(sorted(choices))}")


def ensure_string_list(values: list[str], key: str) -> None:
    if not values:
        raise ConfigError(f"Config value '{key}' must contain at least one value")
    for value in values:
        ensure_nonempty(value, key)


def validate_project_config(project: ProjectConfig) -> None:
    ensure_nonempty(project.name, "project.name")


def validate_sample_config(sample: SampleConfig) -> None:
    ensure_nonempty(sample.locus_tag, "sample.locus_tag")
    ensure_positive(sample.locus_tag_digits, "sample.locus_tag_digits")
    ensure_nonempty(sample.scientific_name, "sample.scientific_name")
    ensure_collection_date(sample.collection_date, "sample.collection_date")
    try:
        genetic_code = int(sample.genetic_code)
    except ValueError as exc:
        raise ConfigError("Config value 'sample.genetic_code' must be an NCBI table number") from exc
    if genetic_code not in CodonTable.unambiguous_dna_by_id:
        raise ConfigError(f"Unknown NCBI genetic code table: {genetic_code}")


def validate_submission_config(submission: SubmissionConfig) -> None:
    if submission.hold_date:
        if not HOLD_DATE_RE.match(submission.hold_date):
            raise ConfigError("Config value 'submission.hold_date' must match YYYYMMDD")
        ensure_date(submission.hold_date, "submission.hold_date", "%Y%m%d")
    ensure_nonempty(submission.bioproject, "submission.bioproject")
    ensure_nonempty(submission.biosample, "submission.biosample")
    ensure_nonempty(submission.datatype, "submission.datatype")


def validate_submitter_config(submitter: SubmitterConfig) -> None:
    ensure_string_list(submitter.ab_name, "submitter.ab_name")
    ensure_nonempty(submitter.contact, "submitter.contact")
    ensure_nonempty(submitter.email, "submitter.email")


def validate_reference_config(reference: ReferenceConfig) -> None:
    ensure_string_list(reference.ab_name, "reference.ab_name")
    ensure_nonempty(reference.title, "reference.title")
    ensure_positive(reference.year, "reference.year")


def validate_pipeline_config(pipeline: PipelineConfig) -> None:
    ensure_nonnegative(pipeline.gapjust_gap_len, "pipeline.gapjust_gap_len")
    ensure_nonnegative(pipeline.gapjust_min, "pipeline.gapjust_min")
    ensure_nonnegative(pipeline.gapjust_max, "pipeline.gapjust_max")
    if pipeline.gapjust_min > pipeline.gapjust_max:
        raise ConfigError("Config value 'pipeline.gapjust_min' must be <= pipeline.gapjust_max")
    ensure_choice(pipeline.gap_assembly, "pipeline.gap_assembly", VALID_GAP_ASSEMBLY)
    ensure_choice(
        pipeline.feature_with_gap,
        "pipeline.feature_with_gap",
        VALID_FEATURE_WITH_GAP,
    )
    ensure_nonnegative(pipeline.min_assembly_gap, "pipeline.min_assembly_gap")
    ensure_nonnegative(
        pipeline.min_artificial_intron_size,
        "pipeline.min_artificial_intron_size",
    )
    ensure_choice(
        pipeline.coordinate_duplicate_policy,
        "pipeline.coordinate_duplicate_policy",
        COORDINATE_DUPLICATE_POLICIES,
    )
    for pattern in pipeline.replace_product_patterns:
        try:
            re.compile(pattern)
        except re.error as exc:
            raise ConfigError(
                f"Invalid regex in 'pipeline.replace_product_patterns': {pattern!r}: {exc}"
            ) from exc


def validate_plots_config(plots: PlotsConfig) -> None:
    ensure_positive(
        plots.coordinate_duplicate_limit,
        "plots.coordinate_duplicate_limit",
    )


def validate_tools_config(tools: ToolsConfig) -> None:
    ensure_nonempty(tools.java, "tools.java")
    if not JAVA_HEAP_RE.match(tools.java_heap):
        raise ConfigError("Config value 'tools.java_heap' must look like 16G or 512M")


def validate_databases_config(databases: DatabasesConfig) -> None:
    ensure_nonempty(databases.root, "databases.root")
    for key, value in (
        ("lock_poll_seconds", databases.lock_poll_seconds),
        ("lock_timeout_seconds", databases.lock_timeout_seconds),
        ("lock_heartbeat_seconds", databases.lock_heartbeat_seconds),
        ("lock_stale_seconds", databases.lock_stale_seconds),
    ):
        if not math.isfinite(value) or value <= 0:
            raise ConfigError(f"Config value 'databases.{key}' must be finite and > 0")
    if databases.lock_stale_seconds <= databases.lock_heartbeat_seconds:
        raise ConfigError(
            "Config value 'databases.lock_stale_seconds' must be greater than "
            "'databases.lock_heartbeat_seconds'"
        )


def validate_busco_config(busco: BuscoConfig) -> None:
    ensure_nonempty(busco.command, "busco.command")
    if not busco.run_cds and not busco.run_genome:
        raise ConfigError("At least one of 'busco.run_cds' or 'busco.run_genome' must be true")
    ensure_choice(busco.mode, "busco.mode", VALID_BUSCO_MODE)
    ensure_choice(busco.cds_mode, "busco.cds_mode", VALID_BUSCO_MODE)
    ensure_choice(
        busco.auto_lineage_scope,
        "busco.auto_lineage_scope",
        VALID_BUSCO_AUTO_LINEAGE_SCOPE,
    )
    ensure_positive(busco.threads, "busco.threads")
    if not busco.auto_lineage and not busco.lineage_dataset.strip():
        raise ConfigError(
            "Config value 'busco.lineage_dataset' is required when busco.auto_lineage = false"
        )


def _ensure_range(value: float, key: str, minimum: float, maximum: float) -> None:
    if not minimum <= value <= maximum:
        raise ConfigError(f"Config value '{key}' must be between {minimum:g} and {maximum:g}")


def validate_functional_annotation_config(annotation: FunctionalAnnotationConfig) -> None:
    ensure_nonempty(annotation.diamond_command, "functional_annotation.diamond_command")
    ensure_nonempty(annotation.hmmscan_command, "functional_annotation.hmmscan_command")
    ensure_nonempty(annotation.hmmpress_command, "functional_annotation.hmmpress_command")
    ensure_nonempty(annotation.rpsblast_command, "functional_annotation.rpsblast_command")
    ensure_nonempty(annotation.rpsbproc_command, "functional_annotation.rpsbproc_command")
    ensure_positive(annotation.threads, "functional_annotation.threads")
    ensure_choice(
        annotation.sensitivity,
        "functional_annotation.sensitivity",
        VALID_DIAMOND_SENSITIVITY,
    )
    if annotation.evalue <= 0:
        raise ConfigError("Config value 'functional_annotation.evalue' must be > 0")
    ensure_positive(annotation.max_target_seqs, "functional_annotation.max_target_seqs")
    if annotation.min_bitscore < 0:
        raise ConfigError("Config value 'functional_annotation.min_bitscore' must be >= 0")
    _ensure_range(annotation.min_identity, "functional_annotation.min_identity", 0, 100)
    _ensure_range(
        annotation.min_query_coverage,
        "functional_annotation.min_query_coverage",
        0,
        100,
    )
    _ensure_range(
        annotation.min_subject_coverage,
        "functional_annotation.min_subject_coverage",
        0,
        100,
    )
    _ensure_range(
        annotation.near_top_bitscore_ratio,
        "functional_annotation.near_top_bitscore_ratio",
        0,
        1,
    )
    _ensure_range(
        annotation.min_token_score,
        "functional_annotation.min_token_score",
        0,
        1,
    )
    if (
        annotation.swissprot_weight <= 0
        or annotation.uniref90_weight <= 0
        or annotation.reference_weight <= 0
    ):
        raise ConfigError("Functional annotation database weights must be > 0")
    if annotation.swissprot_enabled and not (
        annotation.swissprot_fasta.strip() or annotation.swissprot_url.strip()
    ):
        raise ConfigError(
            "Set functional_annotation.swissprot_fasta or functional_annotation.swissprot_url"
        )
    if annotation.reference_proteins.strip():
        ensure_nonempty(annotation.reference_name, "functional_annotation.reference_name")
        if annotation.swissprot_enabled and annotation.reference_name.casefold() == "swissprot":
            raise ConfigError(
                "Config value 'functional_annotation.reference_name' must not be "
                "'swissprot' when Swiss-Prot is enabled"
            )
    if annotation.uniref90_taxon_id < 0:
        raise ConfigError("Config value 'functional_annotation.uniref90_taxon_id' must be >= 0")
    if annotation.uniref90_enabled and not (
        annotation.uniref90_fasta.strip()
        or annotation.uniref90_url.strip()
        or annotation.uniref90_taxon_id > 0
    ):
        raise ConfigError(
            "Set functional_annotation.uniref90_fasta, uniref90_url, or "
            "uniref90_taxon_id"
        )
    if annotation.pfam_max_i_evalue <= 0:
        raise ConfigError("Config value 'functional_annotation.pfam_max_i_evalue' must be > 0")
    _ensure_range(
        annotation.pfam_min_domain_coverage,
        "functional_annotation.pfam_min_domain_coverage",
        0,
        1,
    )
    if annotation.cdd_enabled and not (
        annotation.cdd_database.strip() or annotation.cdd_url.strip()
    ):
        raise ConfigError(
            "Set functional_annotation.cdd_database or functional_annotation.cdd_url"
        )
    if annotation.cdd_enabled and not (
        annotation.cdd_data_dir.strip() or annotation.cdd_data_url.strip()
    ):
        raise ConfigError(
            "Set functional_annotation.cdd_data_dir or functional_annotation.cdd_data_url"
        )
    if annotation.cdd_evalue <= 0:
        raise ConfigError("Config value 'functional_annotation.cdd_evalue' must be > 0")
    if annotation.pfam_enabled and not (annotation.pfam_hmm.strip() or annotation.pfam_url.strip()):
        raise ConfigError("Set functional_annotation.pfam_hmm or functional_annotation.pfam_url")
    if annotation.enabled and not (
        annotation.swissprot_enabled
        or annotation.uniref90_enabled
        or annotation.reference_proteins.strip()
        or annotation.pfam_enabled
        or annotation.cdd_enabled
    ):
        raise ConfigError("Functional annotation is enabled but no annotation source is enabled")

    taxonomy = annotation.taxonomy
    if taxonomy.target_taxon_id < 0:
        raise ConfigError(
            "Config value 'functional_annotation.taxonomy.target_taxon_id' must be >= 0"
        )
    _ensure_range(
        taxonomy.distant_specificity_identity,
        "functional_annotation.taxonomy.distant_specificity_identity",
        0,
        100,
    )

    consistency = annotation.consistency
    for key, value in (
        ("near_identical_identity", consistency.near_identical_identity),
        ("near_identical_coverage", consistency.near_identical_coverage),
        ("family_identity", consistency.family_identity),
        ("family_coverage", consistency.family_coverage),
        ("broad_identity", consistency.broad_identity),
        ("broad_coverage", consistency.broad_coverage),
    ):
        _ensure_range(value, f"functional_annotation.consistency.{key}", 0, 100)
    if consistency.evalue <= 0:
        raise ConfigError(
            "Config value 'functional_annotation.consistency.evalue' must be > 0"
        )
    _ensure_range(
        consistency.name_similarity_threshold,
        "functional_annotation.consistency.name_similarity_threshold",
        0,
        1,
    )
    ensure_positive(
        consistency.source_pair_min_pairs,
        "functional_annotation.consistency.source_pair_min_pairs",
    )
    if not (
        consistency.near_identical_identity >= consistency.family_identity
        >= consistency.broad_identity
    ):
        raise ConfigError(
            "Functional annotation consistency identity thresholds must satisfy "
            "near_identical >= family >= broad"
        )
    if not (
        consistency.near_identical_coverage >= consistency.family_coverage
        >= consistency.broad_coverage
    ):
        raise ConfigError(
            "Functional annotation consistency coverage thresholds must satisfy "
            "near_identical >= family >= broad"
        )
    if consistency.enabled and not annotation.enabled:
        raise ConfigError(
            "functional_annotation.consistency.enabled requires functional_annotation.enabled"
        )


def validate_config(config: MSSPackConfig) -> None:
    validate_project_config(config.project)
    validate_sample_config(config.sample)
    validate_submission_config(config.submission)
    validate_submitter_config(config.submitter)
    validate_reference_config(config.reference)
    validate_pipeline_config(config.pipeline)
    validate_plots_config(config.plots)
    validate_tools_config(config.tools)
    validate_databases_config(config.databases)
    validate_busco_config(config.busco)
    validate_functional_annotation_config(config.functional_annotation)
