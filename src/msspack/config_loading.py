from __future__ import annotations

from pathlib import Path
from typing import Any

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover
    import tomli as tomllib

from .config_errors import ConfigError
from .config_models import (
    BuscoConfig,
    InputsConfig,
    PipelineConfig,
    ProjectConfig,
    ReferenceConfig,
    SampleConfig,
    StCommentConfig,
    SubmissionConfig,
    SubmitterConfig,
    ToolsConfig,
)


def require(mapping: dict[str, Any], key: str, section: str) -> Any:
    if key not in mapping:
        raise ConfigError(f"Missing required key '{section}.{key}'")
    return mapping[key]


def listify(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, list):
        return [str(item) for item in value]
    return [str(value)]


def section(data: dict[str, Any], name: str) -> dict[str, Any]:
    value = data.get(name, {})
    if isinstance(value, dict):
        return value
    raise ConfigError(f"Section '{name}' must be a table/object")


def load_project_config(data: dict[str, Any]) -> ProjectConfig:
    return ProjectConfig(
        name=str(require(data, "name", "project")),
        output_dir=data.get("output_dir"),
    )


def load_inputs_config(data: dict[str, Any]) -> InputsConfig:
    return InputsConfig(
        fasta=str(require(data, "fasta", "inputs")),
        gff=str(require(data, "gff", "inputs")),
    )


def load_sample_config(data: dict[str, Any]) -> SampleConfig:
    return SampleConfig(
        locus_tag=str(require(data, "locus_tag", "sample")),
        locus_tag_digits=int(data.get("locus_tag_digits", 6)),
        scientific_name=str(require(data, "scientific_name", "sample")),
        strain=str(data.get("strain", "")),
        isolate=str(data.get("isolate", "")),
        geo_loc_name=str(data.get("geo_loc_name", "")),
        collection_date=str(data.get("collection_date", "")),
        sex=str(data.get("sex", "")),
        linkage_evidence=str(data.get("linkage_evidence", "paired-ends")),
        genetic_code=str(data.get("genetic_code", "1")),
    )


def load_submission_config(data: dict[str, Any]) -> SubmissionConfig:
    return SubmissionConfig(
        datatype=str(data.get("datatype", "WGS")),
        hold_date=str(require(data, "hold_date", "submission")),
        bioproject=str(require(data, "bioproject", "submission")),
        biosample=str(require(data, "biosample", "submission")),
        sra=listify(data.get("sra", [])),
        keywords=listify(data.get("keywords", ["WGS"])),
    )


def load_submitter_config(data: dict[str, Any]) -> SubmitterConfig:
    return SubmitterConfig(
        ab_name=listify(require(data, "ab_name", "submitter")),
        contact=str(require(data, "contact", "submitter")),
        institute=str(require(data, "institute", "submitter")),
        department=str(require(data, "department", "submitter")),
        country=str(require(data, "country", "submitter")),
        state=str(require(data, "state", "submitter")),
        city=str(require(data, "city", "submitter")),
        street=str(require(data, "street", "submitter")),
        zip=str(require(data, "zip", "submitter")),
        phone=str(require(data, "phone", "submitter")),
        email=str(require(data, "email", "submitter")),
    )


def load_reference_config(data: dict[str, Any]) -> ReferenceConfig:
    return ReferenceConfig(
        title=str(require(data, "title", "reference")),
        ab_name=listify(require(data, "ab_name", "reference")),
        year=int(require(data, "year", "reference")),
        status=str(data.get("status", "Unpublished")),
    )


def load_st_comment_config(data: dict[str, Any]) -> StCommentConfig:
    return StCommentConfig(
        tagset_id=str(data.get("tagset_id", "Genome-Assembly-Data")),
        assembly_method=str(data.get("assembly_method", "")),
        assembly_name=str(data.get("assembly_name", "")),
        genome_coverage=str(data.get("genome_coverage", "")),
        sequencing_technology=str(data.get("sequencing_technology", "")),
    )


def load_pipeline_config(data: dict[str, Any]) -> PipelineConfig:
    return PipelineConfig(
        run_gapjust=bool(data.get("run_gapjust", True)),
        gapjust_gap_len=int(data.get("gapjust_gap_len", 100)),
        gapjust_min=int(data.get("gapjust_min", 80)),
        gapjust_max=int(data.get("gapjust_max", 120)),
        gap_assembly=str(data.get("gap_assembly", "unknown")),
        infer_complete=bool(data.get("infer_complete", True)),
        feature_with_gap=str(data.get("feature_with_gap", "misc_feature")),
        min_assembly_gap=int(data.get("min_assembly_gap", 50)),
        min_artificial_intron_size=int(data.get("min_artificial_intron_size", 10)),
        replace_product_patterns=listify(data.get("replace_product_patterns", [])),
        replace_product_with=str(data.get("replace_product_with", "hypothetical protein")),
        validate_with_parser=bool(data.get("validate_with_parser", True)),
        validate_with_transchecker=bool(data.get("validate_with_transchecker", True)),
        validate_in_parallel=bool(data.get("validate_in_parallel", True)),
    )


def load_tools_config(data: dict[str, Any]) -> ToolsConfig:
    return ToolsConfig(
        cache_dir=data.get("cache_dir"),
        java=str(data.get("java", "java")),
        gff3sort=data.get("gff3sort"),
        java_heap=str(data.get("java_heap", "16G")),
    )


def load_busco_config(data: dict[str, Any]) -> BuscoConfig:
    return BuscoConfig(
        command=str(data.get("command", "busco")),
        run_cds=bool(data.get("run_cds", True)),
        run_genome=bool(data.get("run_genome", False)),
        mode=str(data.get("mode", "genome")),
        cds_mode=str(data.get("cds_mode", "transcriptome")),
        lineage_dataset=str(data.get("lineage_dataset", "")),
        auto_lineage=bool(data.get("auto_lineage", True)),
        auto_lineage_scope=str(data.get("auto_lineage_scope", "all")),
        threads=int(data.get("threads", 1)),
        download_path=str(data.get("download_path", "")),
        offline=bool(data.get("offline", False)),
        force=bool(data.get("force", False)),
        opt_out_run_stats=bool(data.get("opt_out_run_stats", True)),
    )


def read_config_data(config_path: Path) -> dict[str, Any]:
    return tomllib.loads(config_path.read_text(encoding="utf-8"))


def load_sections(data: dict[str, Any]) -> dict[str, Any]:
    return {
        "project": load_project_config(section(data, "project")),
        "inputs": load_inputs_config(section(data, "inputs")),
        "sample": load_sample_config(section(data, "sample")),
        "submission": load_submission_config(section(data, "submission")),
        "submitter": load_submitter_config(section(data, "submitter")),
        "reference": load_reference_config(section(data, "reference")),
        "st_comment": load_st_comment_config(section(data, "st_comment")),
        "pipeline": load_pipeline_config(section(data, "pipeline")),
        "tools": load_tools_config(section(data, "tools")),
        "busco": load_busco_config(section(data, "busco")),
    }
