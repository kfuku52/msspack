from __future__ import annotations

import tomllib
from pathlib import Path
from typing import Any

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

ExpectedType = type[Any] | tuple[type[Any], ...]


SECTION_TYPES: dict[str, dict[str, ExpectedType]] = {
    "project": {"name": str, "output_dir": str},
    "inputs": {"fasta": str, "gff": str},
    "sample": {
        "locus_tag": str,
        "locus_tag_digits": int,
        "scientific_name": str,
        "strain": str,
        "isolate": str,
        "geo_loc_name": str,
        "collection_date": str,
        "sex": str,
        "linkage_evidence": str,
        "genetic_code": str,
    },
    "submission": {
        "datatype": str,
        "hold_date": str,
        "bioproject": str,
        "biosample": str,
        "sra": (str, list),
        "keywords": (str, list),
    },
    "submitter": {
        "ab_name": (str, list),
        "contact": str,
        "institute": str,
        "department": str,
        "country": str,
        "state": str,
        "city": str,
        "street": str,
        "zip": str,
        "phone": str,
        "email": str,
    },
    "reference": {
        "title": str,
        "ab_name": (str, list),
        "year": int,
        "status": str,
    },
    "st_comment": {
        "tagset_id": str,
        "assembly_method": str,
        "assembly_name": str,
        "genome_coverage": str,
        "sequencing_technology": str,
    },
    "pipeline": {
        "run_gapjust": bool,
        "gapjust_gap_len": int,
        "gapjust_min": int,
        "gapjust_max": int,
        "gap_assembly": str,
        "infer_complete": bool,
        "feature_with_gap": str,
        "min_assembly_gap": int,
        "min_artificial_intron_size": int,
        "replace_product_patterns": (str, list),
        "replace_product_with": str,
        "validate_with_parser": bool,
        "validate_with_transchecker": bool,
        "validate_in_parallel": bool,
    },
    "tools": {
        "cache_dir": str,
        "java": str,
        "gff3sort": str,
        "java_heap": str,
    },
    "busco": {
        "command": str,
        "run_cds": bool,
        "run_genome": bool,
        "mode": str,
        "cds_mode": str,
        "lineage_dataset": str,
        "auto_lineage": bool,
        "auto_lineage_scope": str,
        "threads": int,
        "download_path": str,
        "offline": bool,
        "force": bool,
        "opt_out_run_stats": bool,
    },
}


def _type_label(expected: ExpectedType) -> str:
    if isinstance(expected, tuple):
        return " or ".join(item.__name__ for item in expected)
    return expected.__name__


def _validate_raw_config(data: dict[str, Any]) -> None:
    unknown_sections = sorted(set(data) - set(SECTION_TYPES))
    if unknown_sections:
        raise ConfigError(f"Unknown config section(s): {', '.join(unknown_sections)}")
    for section_name, values in data.items():
        if not isinstance(values, dict):
            raise ConfigError(f"Section '{section_name}' must be a table/object")
        schema = SECTION_TYPES[section_name]
        unknown_keys = sorted(set(values) - set(schema))
        if unknown_keys:
            dotted = ", ".join(f"{section_name}.{key}" for key in unknown_keys)
            raise ConfigError(f"Unknown config key(s): {dotted}")
        for key, value in values.items():
            expected = schema[key]
            if expected is int and isinstance(value, bool):
                valid = False
            else:
                valid = isinstance(value, expected)
            if not valid:
                raise ConfigError(
                    f"Config value '{section_name}.{key}' must be {_type_label(expected)}, "
                    f"got {type(value).__name__}"
                )
            if isinstance(value, list) and any(not isinstance(item, str) for item in value):
                raise ConfigError(
                    f"Config value '{section_name}.{key}' must contain only strings"
                )
            strings = value if isinstance(value, list) else [value]
            if any(
                isinstance(item, str)
                and any(character in item for character in ("\x00", "\r", "\n", "\t"))
                for item in strings
            ):
                raise ConfigError(
                    f"Config value '{section_name}.{key}' must not contain control characters"
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
    try:
        payload: object = tomllib.loads(config_path.read_text(encoding="utf-8"))
    except OSError as exc:
        raise ConfigError(f"Could not read config file {config_path}: {exc}") from exc
    except tomllib.TOMLDecodeError as exc:
        raise ConfigError(f"Invalid TOML in {config_path}: {exc}") from exc
    if not isinstance(payload, dict):
        raise ConfigError(f"TOML root must be a table in {config_path}")
    return payload


def load_sections(data: dict[str, Any]) -> dict[str, Any]:
    _validate_raw_config(data)
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
