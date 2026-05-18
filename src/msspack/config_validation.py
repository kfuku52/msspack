from __future__ import annotations

import re
from datetime import datetime

from .config_errors import ConfigError
from .config_models import (
    BuscoConfig,
    MSSPackConfig,
    PipelineConfig,
    ProjectConfig,
    ReferenceConfig,
    SampleConfig,
    SubmissionConfig,
    SubmitterConfig,
    ToolsConfig,
)

HOLD_DATE_RE = re.compile(r"^\d{8}$")
JAVA_HEAP_RE = re.compile(r"^\d+[KMGkmg]$")
VALID_GAP_ASSEMBLY = {"known", "unknown"}
VALID_FEATURE_WITH_GAP = {"asis", "misc_feature"}
VALID_BUSCO_MODE = {"genome", "proteins", "transcriptome"}
VALID_BUSCO_AUTO_LINEAGE_SCOPE = {"all", "euk", "prok"}


def ensure_nonempty(value: str, key: str) -> None:
    if not value.strip():
        raise ConfigError(f"Config value '{key}' must not be empty")


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
    ensure_date(sample.collection_date, "sample.collection_date", "%Y-%m-%d")


def validate_submission_config(submission: SubmissionConfig) -> None:
    if not HOLD_DATE_RE.match(submission.hold_date):
        raise ConfigError("Config value 'submission.hold_date' must match YYYYMMDD")
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


def validate_tools_config(tools: ToolsConfig) -> None:
    ensure_nonempty(tools.java, "tools.java")
    if not JAVA_HEAP_RE.match(tools.java_heap):
        raise ConfigError("Config value 'tools.java_heap' must look like 16G or 512M")


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


def validate_config(config: MSSPackConfig) -> None:
    validate_project_config(config.project)
    validate_sample_config(config.sample)
    validate_submission_config(config.submission)
    validate_submitter_config(config.submitter)
    validate_reference_config(config.reference)
    validate_pipeline_config(config.pipeline)
    validate_tools_config(config.tools)
    validate_busco_config(config.busco)
