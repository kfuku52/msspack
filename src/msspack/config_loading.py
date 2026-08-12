from __future__ import annotations

import tomllib
from pathlib import Path
from typing import Any

from .config_errors import ConfigError
from .config_models import (
    BuscoConfig,
    DatabasesConfig,
    FunctionalAnnotationConfig,
    FunctionalAnnotationConsistencyConfig,
    FunctionalAnnotationTaxonomyConfig,
    InputsConfig,
    PipelineConfig,
    PlotsConfig,
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
        "tissue_type": str,
        "isolation_source": str,
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
        "coordinate_duplicate_policy": str,
        "replace_product_patterns": (str, list),
        "replace_product_with": str,
        "validate_with_parser": bool,
        "validate_with_transchecker": bool,
        "validate_in_parallel": bool,
    },
    "plots": {
        "coordinate_duplicate_limit": int,
    },
    "tools": {
        "cache_dir": str,
        "java": str,
        "gff3sort": str,
        "java_heap": str,
    },
    "databases": {
        "root": str,
        "lock_poll_seconds": (int, float),
        "lock_timeout_seconds": (int, float),
        "lock_heartbeat_seconds": (int, float),
        "lock_stale_seconds": (int, float),
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
    "functional_annotation": {
        "enabled": bool,
        "diamond_command": str,
        "hmmscan_command": str,
        "hmmpress_command": str,
        "rpsblast_command": str,
        "rpsbproc_command": str,
        "threads": int,
        "sensitivity": str,
        "evalue": (int, float),
        "max_target_seqs": int,
        "min_bitscore": (int, float),
        "min_identity": (int, float),
        "min_query_coverage": (int, float),
        "min_subject_coverage": (int, float),
        "near_top_bitscore_ratio": (int, float),
        "min_token_score": (int, float),
        "overwrite_existing": bool,
        "swissprot_enabled": bool,
        "swissprot_fasta": str,
        "swissprot_url": str,
        "swissprot_weight": (int, float),
        "uniref90_enabled": bool,
        "uniref90_fasta": str,
        "uniref90_url": str,
        "uniref90_taxon_id": int,
        "uniref90_weight": (int, float),
        "reference_proteins": str,
        "reference_name": str,
        "reference_weight": (int, float),
        "pfam_enabled": bool,
        "pfam_hmm": str,
        "pfam_url": str,
        "pfam_max_i_evalue": (int, float),
        "pfam_min_domain_coverage": (int, float),
        "cdd_enabled": bool,
        "cdd_database": str,
        "cdd_url": str,
        "cdd_data_dir": str,
        "cdd_data_url": str,
        "cdd_evalue": (int, float),
        "taxonomy": dict,
        "consistency": dict,
    },
}

FUNCTIONAL_ANNOTATION_TAXONOMY_TYPES: dict[str, ExpectedType] = {
    "enabled": bool,
    "target_taxon_id": int,
    "resolve_scientific_name": bool,
    "offline": bool,
    "busco_crosscheck": bool,
    "strict": bool,
    "distant_specificity_identity": (int, float),
}

FUNCTIONAL_ANNOTATION_CONSISTENCY_TYPES: dict[str, ExpectedType] = {
    "enabled": bool,
    "harmonize_safe_equivalents": bool,
    "auto_resolve_conflicts": bool,
    "near_identical_identity": (int, float),
    "near_identical_coverage": (int, float),
    "family_identity": (int, float),
    "family_coverage": (int, float),
    "broad_identity": (int, float),
    "broad_coverage": (int, float),
    "evalue": (int, float),
    "name_similarity_threshold": (int, float),
    "source_pair_min_pairs": int,
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
            if isinstance(value, bool) and (
                expected is int
                or expected is float
                or (isinstance(expected, tuple) and any(item in (int, float) for item in expected))
            ):
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
            if section_name == "functional_annotation" and key in {"taxonomy", "consistency"}:
                assert isinstance(value, dict)
                nested_schema = (
                    FUNCTIONAL_ANNOTATION_TAXONOMY_TYPES
                    if key == "taxonomy"
                    else FUNCTIONAL_ANNOTATION_CONSISTENCY_TYPES
                )
                unknown_consistency_keys = sorted(
                    set(value) - set(nested_schema)
                )
                if unknown_consistency_keys:
                    dotted = ", ".join(
                        f"functional_annotation.{key}.{item}"
                        for item in unknown_consistency_keys
                    )
                    raise ConfigError(f"Unknown config key(s): {dotted}")
                for nested_key, nested_value in value.items():
                    nested_expected = nested_schema[nested_key]
                    nested_numeric_bool = isinstance(nested_value, bool) and (
                        nested_expected is int
                        or nested_expected is float
                        or (
                            isinstance(nested_expected, tuple)
                            and any(item in (int, float) for item in nested_expected)
                        )
                    )
                    if nested_numeric_bool or not isinstance(nested_value, nested_expected):
                        raise ConfigError(
                            f"Config value 'functional_annotation.{key}."
                            f"{nested_key}' must be {_type_label(nested_expected)}, "
                            f"got {type(nested_value).__name__}"
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
        tissue_type=str(data.get("tissue_type", "")),
        isolation_source=str(data.get("isolation_source", "")),
        geo_loc_name=str(data.get("geo_loc_name", "")),
        collection_date=str(data.get("collection_date", "")),
        sex=str(data.get("sex", "")),
        linkage_evidence=str(data.get("linkage_evidence", "paired-ends")),
        genetic_code=str(data.get("genetic_code", "1")),
    )


def load_submission_config(data: dict[str, Any]) -> SubmissionConfig:
    return SubmissionConfig(
        datatype=str(data.get("datatype", "WGS")),
        hold_date=str(data.get("hold_date", "")),
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
        coordinate_duplicate_policy=str(
            data.get("coordinate_duplicate_policy", "longest_valid_cds")
        ),
        replace_product_patterns=listify(data.get("replace_product_patterns", [])),
        replace_product_with=str(data.get("replace_product_with", "hypothetical protein")),
        validate_with_parser=bool(data.get("validate_with_parser", True)),
        validate_with_transchecker=bool(data.get("validate_with_transchecker", True)),
        validate_in_parallel=bool(data.get("validate_in_parallel", True)),
    )


def load_plots_config(data: dict[str, Any]) -> PlotsConfig:
    return PlotsConfig(
        coordinate_duplicate_limit=int(data.get("coordinate_duplicate_limit", 50)),
    )


def load_tools_config(data: dict[str, Any]) -> ToolsConfig:
    return ToolsConfig(
        cache_dir=data.get("cache_dir"),
        java=str(data.get("java", "java")),
        gff3sort=data.get("gff3sort"),
        java_heap=str(data.get("java_heap", "16G")),
    )


def load_databases_config(data: dict[str, Any]) -> DatabasesConfig:
    defaults = DatabasesConfig()
    return DatabasesConfig(
        root=str(data.get("root", defaults.root)),
        lock_poll_seconds=float(
            data.get("lock_poll_seconds", defaults.lock_poll_seconds)
        ),
        lock_timeout_seconds=float(
            data.get("lock_timeout_seconds", defaults.lock_timeout_seconds)
        ),
        lock_heartbeat_seconds=float(
            data.get("lock_heartbeat_seconds", defaults.lock_heartbeat_seconds)
        ),
        lock_stale_seconds=float(
            data.get("lock_stale_seconds", defaults.lock_stale_seconds)
        ),
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


def load_functional_annotation_config(data: dict[str, Any]) -> FunctionalAnnotationConfig:
    defaults = FunctionalAnnotationConfig()
    taxonomy_data = data.get("taxonomy", {})
    if not isinstance(taxonomy_data, dict):
        taxonomy_data = {}
    taxonomy_defaults = FunctionalAnnotationTaxonomyConfig()
    consistency_data = data.get("consistency", {})
    if not isinstance(consistency_data, dict):
        consistency_data = {}
    consistency_defaults = FunctionalAnnotationConsistencyConfig()
    return FunctionalAnnotationConfig(
        enabled=bool(data.get("enabled", defaults.enabled)),
        diamond_command=str(data.get("diamond_command", defaults.diamond_command)),
        hmmscan_command=str(data.get("hmmscan_command", defaults.hmmscan_command)),
        hmmpress_command=str(data.get("hmmpress_command", defaults.hmmpress_command)),
        rpsblast_command=str(data.get("rpsblast_command", defaults.rpsblast_command)),
        rpsbproc_command=str(data.get("rpsbproc_command", defaults.rpsbproc_command)),
        threads=int(data.get("threads", defaults.threads)),
        sensitivity=str(data.get("sensitivity", defaults.sensitivity)),
        evalue=float(data.get("evalue", defaults.evalue)),
        max_target_seqs=int(data.get("max_target_seqs", defaults.max_target_seqs)),
        min_bitscore=float(data.get("min_bitscore", defaults.min_bitscore)),
        min_identity=float(data.get("min_identity", defaults.min_identity)),
        min_query_coverage=float(data.get("min_query_coverage", defaults.min_query_coverage)),
        min_subject_coverage=float(data.get("min_subject_coverage", defaults.min_subject_coverage)),
        near_top_bitscore_ratio=float(
            data.get("near_top_bitscore_ratio", defaults.near_top_bitscore_ratio)
        ),
        min_token_score=float(data.get("min_token_score", defaults.min_token_score)),
        overwrite_existing=bool(data.get("overwrite_existing", defaults.overwrite_existing)),
        swissprot_enabled=bool(data.get("swissprot_enabled", defaults.swissprot_enabled)),
        swissprot_fasta=str(data.get("swissprot_fasta", defaults.swissprot_fasta)),
        swissprot_url=str(data.get("swissprot_url", defaults.swissprot_url)),
        swissprot_weight=float(data.get("swissprot_weight", defaults.swissprot_weight)),
        uniref90_enabled=bool(data.get("uniref90_enabled", defaults.uniref90_enabled)),
        uniref90_fasta=str(data.get("uniref90_fasta", defaults.uniref90_fasta)),
        uniref90_url=str(data.get("uniref90_url", defaults.uniref90_url)),
        uniref90_taxon_id=int(data.get("uniref90_taxon_id", defaults.uniref90_taxon_id)),
        uniref90_weight=float(data.get("uniref90_weight", defaults.uniref90_weight)),
        reference_proteins=str(data.get("reference_proteins", defaults.reference_proteins)),
        reference_name=str(data.get("reference_name", defaults.reference_name)),
        reference_weight=float(data.get("reference_weight", defaults.reference_weight)),
        pfam_enabled=bool(data.get("pfam_enabled", defaults.pfam_enabled)),
        pfam_hmm=str(data.get("pfam_hmm", defaults.pfam_hmm)),
        pfam_url=str(data.get("pfam_url", defaults.pfam_url)),
        pfam_max_i_evalue=float(data.get("pfam_max_i_evalue", defaults.pfam_max_i_evalue)),
        pfam_min_domain_coverage=float(
            data.get("pfam_min_domain_coverage", defaults.pfam_min_domain_coverage)
        ),
        cdd_enabled=bool(data.get("cdd_enabled", defaults.cdd_enabled)),
        cdd_database=str(data.get("cdd_database", defaults.cdd_database)),
        cdd_url=str(data.get("cdd_url", defaults.cdd_url)),
        cdd_data_dir=str(data.get("cdd_data_dir", defaults.cdd_data_dir)),
        cdd_data_url=str(data.get("cdd_data_url", defaults.cdd_data_url)),
        cdd_evalue=float(data.get("cdd_evalue", defaults.cdd_evalue)),
        taxonomy=FunctionalAnnotationTaxonomyConfig(
            enabled=bool(taxonomy_data.get("enabled", taxonomy_defaults.enabled)),
            target_taxon_id=int(
                taxonomy_data.get(
                    "target_taxon_id",
                    taxonomy_defaults.target_taxon_id,
                )
            ),
            resolve_scientific_name=bool(
                taxonomy_data.get(
                    "resolve_scientific_name",
                    taxonomy_defaults.resolve_scientific_name,
                )
            ),
            offline=bool(taxonomy_data.get("offline", taxonomy_defaults.offline)),
            busco_crosscheck=bool(
                taxonomy_data.get(
                    "busco_crosscheck",
                    taxonomy_defaults.busco_crosscheck,
                )
            ),
            strict=bool(taxonomy_data.get("strict", taxonomy_defaults.strict)),
            distant_specificity_identity=float(
                taxonomy_data.get(
                    "distant_specificity_identity",
                    taxonomy_defaults.distant_specificity_identity,
                )
            ),
        ),
        consistency=FunctionalAnnotationConsistencyConfig(
            enabled=bool(consistency_data.get("enabled", consistency_defaults.enabled)),
            harmonize_safe_equivalents=bool(
                consistency_data.get(
                    "harmonize_safe_equivalents",
                    consistency_defaults.harmonize_safe_equivalents,
                )
            ),
            auto_resolve_conflicts=bool(
                consistency_data.get(
                    "auto_resolve_conflicts",
                    consistency_defaults.auto_resolve_conflicts,
                )
            ),
            near_identical_identity=float(
                consistency_data.get(
                    "near_identical_identity",
                    consistency_defaults.near_identical_identity,
                )
            ),
            near_identical_coverage=float(
                consistency_data.get(
                    "near_identical_coverage",
                    consistency_defaults.near_identical_coverage,
                )
            ),
            family_identity=float(
                consistency_data.get("family_identity", consistency_defaults.family_identity)
            ),
            family_coverage=float(
                consistency_data.get("family_coverage", consistency_defaults.family_coverage)
            ),
            broad_identity=float(
                consistency_data.get("broad_identity", consistency_defaults.broad_identity)
            ),
            broad_coverage=float(
                consistency_data.get("broad_coverage", consistency_defaults.broad_coverage)
            ),
            evalue=float(consistency_data.get("evalue", consistency_defaults.evalue)),
            name_similarity_threshold=float(
                consistency_data.get(
                    "name_similarity_threshold",
                    consistency_defaults.name_similarity_threshold,
                )
            ),
            source_pair_min_pairs=int(
                consistency_data.get(
                    "source_pair_min_pairs",
                    consistency_defaults.source_pair_min_pairs,
                )
            ),
        ),
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
        "plots": load_plots_config(section(data, "plots")),
        "tools": load_tools_config(section(data, "tools")),
        "databases": load_databases_config(section(data, "databases")),
        "busco": load_busco_config(section(data, "busco")),
        "functional_annotation": load_functional_annotation_config(
            section(data, "functional_annotation")
        ),
    }
