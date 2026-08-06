from __future__ import annotations

import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Literal, cast

from .config import MSSPackConfig, load_config
from .ddbj_tools import ToolInstallation, require_installed, run_parser, run_transchecker
from .execution import (
    NamedJob,
    append_job_if_needed,
    module_origin,
    path_list,
    run_named_jobs,
)
from .step_logging import count_fasta_records
from .utils import MSSPackError, ensure_dir, write_text

ValidationStatus = Literal["passed", "failed", "not_run"]
VALIDATION_COMPONENTS = ("parser", "transchecker")
VALIDATION_LABELS = {
    "parser": "Parser",
    "transchecker": "transChecker",
}
VALIDATION_SUMMARY_SCHEMA_VERSION = 1


@dataclass(frozen=True)
class ValidationCheckResult:
    component: str
    label: str
    status: ValidationStatus
    version: str | None
    log_path: Path | None
    output_paths: dict[str, Path]
    warning_count: int = 0
    error_count: int = 0
    record_counts: dict[str, int] | None = None
    message: str | None = None

    def to_dict(self) -> dict[str, object]:
        payload: dict[str, object] = {
            "component": self.component,
            "label": self.label,
            "status": self.status,
            "version": self.version,
            "log_path": str(self.log_path) if self.log_path is not None else None,
            "outputs": {key: str(path) for key, path in self.output_paths.items()},
            "warning_count": self.warning_count,
            "error_count": self.error_count,
            "record_counts": dict(self.record_counts or {}),
            "message": self.message,
        }
        return payload

    @classmethod
    def from_dict(cls, payload: object, *, source: Path) -> ValidationCheckResult:
        if not isinstance(payload, dict):
            raise MSSPackError(f"Invalid DDBJ validation check in {source}")
        component = payload.get("component")
        label = payload.get("label")
        status_value = payload.get("status")
        if not isinstance(component, str) or component not in VALIDATION_COMPONENTS:
            raise MSSPackError(f"Invalid DDBJ validation component in {source}")
        if not isinstance(label, str) or not label:
            raise MSSPackError(f"Invalid DDBJ validation label in {source}")
        if status_value not in {"passed", "failed", "not_run"}:
            raise MSSPackError(f"Invalid DDBJ validation status in {source}")
        status = cast(ValidationStatus, status_value)
        version_value = payload.get("version")
        if version_value is not None and not isinstance(version_value, str):
            raise MSSPackError(f"Invalid DDBJ validation version in {source}")
        log_value = payload.get("log_path")
        if log_value is not None and not isinstance(log_value, str):
            raise MSSPackError(f"Invalid DDBJ validation log path in {source}")
        raw_outputs = payload.get("outputs", {})
        if not isinstance(raw_outputs, dict) or not all(
            isinstance(key, str) and isinstance(value, str)
            for key, value in raw_outputs.items()
        ):
            raise MSSPackError(f"Invalid DDBJ validation outputs in {source}")
        raw_records = payload.get("record_counts", {})
        if not isinstance(raw_records, dict) or not all(
            isinstance(key, str)
            and isinstance(value, int)
            and not isinstance(value, bool)
            and value >= 0
            for key, value in raw_records.items()
        ):
            raise MSSPackError(f"Invalid DDBJ validation record counts in {source}")
        warning_count = payload.get("warning_count", 0)
        error_count = payload.get("error_count", 0)
        if (
            not isinstance(warning_count, int)
            or isinstance(warning_count, bool)
            or warning_count < 0
            or not isinstance(error_count, int)
            or isinstance(error_count, bool)
            or error_count < 0
        ):
            raise MSSPackError(f"Invalid DDBJ validation diagnostic counts in {source}")
        message = payload.get("message")
        if message is not None and not isinstance(message, str):
            raise MSSPackError(f"Invalid DDBJ validation message in {source}")
        return cls(
            component=component,
            label=label,
            status=status,
            version=version_value,
            log_path=Path(log_value) if isinstance(log_value, str) else None,
            output_paths={
                cast(str, key): Path(cast(str, value))
                for key, value in raw_outputs.items()
            },
            warning_count=warning_count,
            error_count=error_count,
            record_counts={
                cast(str, key): cast(int, value) for key, value in raw_records.items()
            },
            message=message,
        )


@dataclass(frozen=True)
class ValidationSummary:
    path: Path
    status: ValidationStatus
    checks: tuple[ValidationCheckResult, ...]

    @property
    def attempted(self) -> bool:
        return any(check.status != "not_run" for check in self.checks)

    def to_dict(self) -> dict[str, object]:
        return {
            "schema_version": VALIDATION_SUMMARY_SCHEMA_VERSION,
            "status": self.status,
            "attempted": self.attempted,
            "checks": {check.component: check.to_dict() for check in self.checks},
        }


def load_validation_summary(path: Path) -> ValidationSummary | None:
    if not path.is_file():
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise MSSPackError(f"Invalid DDBJ validation summary JSON: {path}") from exc
    if not isinstance(payload, dict):
        raise MSSPackError(f"Invalid DDBJ validation summary: {path}")
    if payload.get("schema_version") != VALIDATION_SUMMARY_SCHEMA_VERSION:
        raise MSSPackError(f"Unsupported DDBJ validation summary schema: {path}")
    status_value = payload.get("status")
    if status_value not in {"passed", "failed", "not_run"}:
        raise MSSPackError(f"Invalid DDBJ validation summary status: {path}")
    raw_checks = payload.get("checks")
    if not isinstance(raw_checks, dict):
        raise MSSPackError(f"Invalid DDBJ validation checks: {path}")
    if set(raw_checks) != set(VALIDATION_COMPONENTS):
        raise MSSPackError(f"Invalid DDBJ validation check set: {path}")
    checks: list[ValidationCheckResult] = []
    for component in VALIDATION_COMPONENTS:
        if component not in raw_checks:
            raise MSSPackError(f"DDBJ validation summary is missing {component}: {path}")
        check = ValidationCheckResult.from_dict(raw_checks[component], source=path)
        if check.component != component:
            raise MSSPackError(f"Mismatched DDBJ validation component in {path}")
        if check.label != VALIDATION_LABELS[component]:
            raise MSSPackError(f"Mismatched DDBJ validation label in {path}")
        checks.append(check)
    typed_checks = tuple(checks)
    expected_status = _summary_status(typed_checks)
    if status_value != expected_status:
        raise MSSPackError(
            f"DDBJ validation summary status does not match its checks: {path}"
        )
    attempted_value = payload.get("attempted")
    if not isinstance(attempted_value, bool):
        raise MSSPackError(f"Invalid DDBJ validation attempted flag: {path}")
    expected_attempted = any(check.status != "not_run" for check in typed_checks)
    if attempted_value != expected_attempted:
        raise MSSPackError(
            f"DDBJ validation attempted flag does not match its checks: {path}"
        )
    return ValidationSummary(
        path=path,
        status=cast(ValidationStatus, status_value),
        checks=typed_checks,
    )


@dataclass(frozen=True)
class ValidationOptions:
    cache_dir: Path | None
    heap: str
    parallel: bool
    java_cmd: str = "java"
    run_parser: bool = True
    run_transchecker: bool = True

    @classmethod
    def from_config(cls, config: MSSPackConfig | None) -> ValidationOptions:
        if config is None:
            return cls(
                cache_dir=None,
                heap="16G",
                parallel=True,
                java_cmd="java",
            )
        return cls(
            cache_dir=config.cache_dir,
            heap=config.tools.java_heap,
            parallel=config.pipeline.validate_in_parallel,
            java_cmd=config.tools.java,
            run_parser=config.pipeline.validate_with_parser,
            run_transchecker=config.pipeline.validate_with_transchecker,
        )

    @classmethod
    def for_explicit_validate(cls, config: MSSPackConfig | None) -> ValidationOptions:
        options = cls.from_config(config)
        return cls(
            cache_dir=options.cache_dir,
            heap=options.heap,
            parallel=options.parallel,
            java_cmd=options.java_cmd,
            run_parser=True,
            run_transchecker=True,
        )


@dataclass(frozen=True)
class ValidationArtifacts:
    ann_path: Path
    fasta_path: Path
    logs_dir: Path
    output_dir: Path

    @classmethod
    def for_existing_submission(
        cls,
        *,
        ann_path: Path,
        fasta_path: Path,
    ) -> ValidationArtifacts:
        base = ann_path.parent
        return cls(
            ann_path=ann_path,
            fasta_path=fasta_path,
            logs_dir=ensure_dir(base / "logs"),
            output_dir=ensure_dir(base / "validation"),
        )

    @property
    def parser_log(self) -> Path:
        return self.logs_dir / "parser.log"

    @property
    def transchecker_log(self) -> Path:
        return self.logs_dir / "transchecker.log"

    @property
    def aa_fasta(self) -> Path:
        return self.output_dir / "transChecker.aa.fasta"

    @property
    def nuc_fasta(self) -> Path:
        return self.output_dir / "transChecker.nuc.fasta"

    @property
    def parser_result(self) -> Path:
        return self.output_dir / "ddbj-parser-result.json"

    @property
    def transchecker_result(self) -> Path:
        return self.output_dir / "ddbj-transchecker-result.json"

    @property
    def validation_summary(self) -> Path:
        return self.output_dir / "ddbj-validation-summary.json"

    def result_paths(
        self,
        *,
        include_parser: bool,
        include_transchecker: bool,
    ) -> dict[str, Path]:
        paths: dict[str, Path] = {"validation_summary": self.validation_summary}
        if include_parser:
            paths["parser_log"] = self.parser_log
        if include_transchecker:
            paths["transchecker_log"] = self.transchecker_log
            paths["aa_fasta"] = self.aa_fasta
            paths["nuc_fasta"] = self.nuc_fasta
        return paths


def _diagnostic_counts(log_path: Path) -> tuple[int, int]:
    if not log_path.is_file():
        return 0, 0
    warnings = 0
    errors = 0
    for line in log_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if re.search(r"(?:^|:)WAR(?=:)|\[\s*WARN(?:ING)?\s*\]", line, re.IGNORECASE):
            warnings += 1
        if re.search(r"(?:^|:)ERR(?=:)|\[\s*ERROR\s*\]", line, re.IGNORECASE):
            errors += 1
    return warnings, errors


def _count_mss_cds_features(ann_path: Path) -> int:
    count = 0
    with ann_path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            fields = line.rstrip("\r\n").split("\t", 2)
            if len(fields) >= 2 and fields[0] == "" and fields[1] == "CDS":
                count += 1
    return count


def _write_check_result(path: Path, result: ValidationCheckResult) -> Path:
    write_text(path, json.dumps(result.to_dict(), indent=2, sort_keys=True) + "\n")
    return path


def _read_check_result(path: Path) -> ValidationCheckResult | None:
    if not path.is_file():
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise MSSPackError(f"Invalid DDBJ validation result JSON: {path}") from exc
    return ValidationCheckResult.from_dict(payload, source=path)


def _result_for_component(
    *,
    component: str,
    installation: ToolInstallation,
    status: ValidationStatus,
    log_path: Path,
    output_paths: dict[str, Path] | None = None,
    record_counts: dict[str, int] | None = None,
    message: str | None = None,
) -> ValidationCheckResult:
    warning_count, error_count = _diagnostic_counts(log_path)
    return ValidationCheckResult(
        component=component,
        label=VALIDATION_LABELS[component],
        status=status,
        version=installation.version,
        log_path=log_path,
        output_paths=output_paths or {},
        warning_count=warning_count,
        error_count=error_count,
        record_counts=record_counts or {},
        message=message,
    )


def _run_parser_with_result(
    *,
    installation: ToolInstallation,
    artifacts: ValidationArtifacts,
    options: ValidationOptions,
) -> None:
    try:
        run_parser(
            installation,
            ann_path=artifacts.ann_path,
            fasta_path=artifacts.fasta_path,
            heap=options.heap,
            java_cmd=options.java_cmd,
            log_path=artifacts.parser_log,
        )
    except Exception as exc:
        _write_check_result(
            artifacts.parser_result,
            _result_for_component(
                component="parser",
                installation=installation,
                status="failed",
                log_path=artifacts.parser_log,
                message=str(exc) or type(exc).__name__,
            ),
        )
        raise
    _write_check_result(
        artifacts.parser_result,
        _result_for_component(
            component="parser",
            installation=installation,
            status="passed",
            log_path=artifacts.parser_log,
        ),
    )


def _run_transchecker_with_result(
    *,
    installation: ToolInstallation,
    artifacts: ValidationArtifacts,
    options: ValidationOptions,
) -> None:
    outputs = {
        "aa_fasta": artifacts.aa_fasta,
        "nuc_fasta": artifacts.nuc_fasta,
    }
    record_counts: dict[str, int] = {}
    try:
        run_transchecker(
            installation,
            ann_path=artifacts.ann_path,
            fasta_path=artifacts.fasta_path,
            aa_out=artifacts.aa_fasta,
            nuc_out=artifacts.nuc_fasta,
            heap=options.heap,
            java_cmd=options.java_cmd,
            log_path=artifacts.transchecker_log,
        )
        expected_records = _count_mss_cds_features(artifacts.ann_path)
        aa_records = count_fasta_records(artifacts.aa_fasta)
        nuc_records = count_fasta_records(artifacts.nuc_fasta)
        record_counts = {
            "annotation_cds": expected_records,
            "aa_fasta": aa_records,
            "nuc_fasta": nuc_records,
        }
        if aa_records != expected_records or nuc_records != expected_records:
            raise MSSPackError(
                "transChecker output record counts do not match the final annotation: "
                f"annotation CDS={expected_records}, AA={aa_records}, "
                f"nucleotide={nuc_records}"
            )
    except Exception as exc:
        _write_check_result(
            artifacts.transchecker_result,
            _result_for_component(
                component="transchecker",
                installation=installation,
                status="failed",
                log_path=artifacts.transchecker_log,
                output_paths={key: path for key, path in outputs.items() if path.is_file()},
                record_counts=record_counts,
                message=str(exc) or type(exc).__name__,
            ),
        )
        raise
    _write_check_result(
        artifacts.transchecker_result,
        _result_for_component(
            component="transchecker",
            installation=installation,
            status="passed",
            log_path=artifacts.transchecker_log,
            output_paths=outputs,
            record_counts=record_counts,
        ),
    )


def _not_run_result(component: str, *, message: str) -> ValidationCheckResult:
    return ValidationCheckResult(
        component=component,
        label=VALIDATION_LABELS[component],
        status="not_run",
        version=None,
        log_path=None,
        output_paths={},
        message=message,
    )


def _summary_status(checks: tuple[ValidationCheckResult, ...]) -> ValidationStatus:
    if any(check.status == "failed" for check in checks):
        return "failed"
    if any(check.status == "passed" for check in checks):
        return "passed"
    return "not_run"


def _write_validation_summary(
    *,
    options: ValidationOptions,
    artifacts: ValidationArtifacts,
    failure_message: str | None = None,
    ignore_existing_results: bool = False,
) -> ValidationSummary:
    requested = {
        "parser": options.run_parser,
        "transchecker": options.run_transchecker,
    }
    result_paths = {
        "parser": artifacts.parser_result,
        "transchecker": artifacts.transchecker_result,
    }
    checks: list[ValidationCheckResult] = []
    missing_requested: list[str] = []
    for component in VALIDATION_COMPONENTS:
        if not requested[component]:
            checks.append(_not_run_result(component, message="disabled by configuration"))
            continue
        result = (
            None
            if ignore_existing_results
            else _read_check_result(result_paths[component])
        )
        if result is None:
            missing_requested.append(component)
        else:
            checks.append(result)

    has_recorded_failure = any(check.status == "failed" for check in checks)
    for index, component in enumerate(missing_requested):
        if failure_message is not None and (has_recorded_failure or index > 0):
            checks.append(
                _not_run_result(
                    component,
                    message="not executed after an earlier validation failure",
                )
            )
            continue
        checks.append(
            ValidationCheckResult(
                component=component,
                label=VALIDATION_LABELS[component],
                status="failed",
                version=None,
                log_path=None,
                output_paths={},
                error_count=1,
                message=failure_message or "validation did not produce a result",
            )
        )
        has_recorded_failure = True
    checks_by_component = {check.component: check for check in checks}
    checks = [checks_by_component[component] for component in VALIDATION_COMPONENTS]
    typed_checks = tuple(checks)
    summary = ValidationSummary(
        path=artifacts.validation_summary,
        status=_summary_status(typed_checks),
        checks=typed_checks,
    )
    write_text(
        artifacts.validation_summary,
        json.dumps(summary.to_dict(), indent=2, sort_keys=True) + "\n",
    )
    return summary


def write_validation_not_run_summary(
    *,
    artifacts: ValidationArtifacts,
    reason: str = "validation disabled",
) -> ValidationSummary:
    checks = tuple(
        _not_run_result(component, message=reason) for component in VALIDATION_COMPONENTS
    )
    summary = ValidationSummary(
        path=artifacts.validation_summary,
        status="not_run",
        checks=checks,
    )
    write_text(
        artifacts.validation_summary,
        json.dumps(summary.to_dict(), indent=2, sort_keys=True) + "\n",
    )
    return summary


def _build_validation_jobs(
    *,
    options: ValidationOptions,
    artifacts: ValidationArtifacts,
) -> list[NamedJob]:
    requested_components: list[str] = []
    if options.run_parser:
        requested_components.append("parser")
    if options.run_transchecker:
        requested_components.append("transchecker")
    if not requested_components:
        return []

    tools = require_installed(
        requested_components,
        cache_dir=options.cache_dir,
    )
    jobs: list[NamedJob] = []

    if options.run_parser:
        job_count = len(jobs)
        append_job_if_needed(
            jobs,
            name="parser",
            outputs=[artifacts.parser_log, artifacts.parser_result],
            dependencies=path_list(
                artifacts.ann_path,
                artifacts.fasta_path,
                tools["parser"].root,
                module_origin("msspack.validation"),
                module_origin("msspack.ddbj_tools"),
            ),
            cache_key={
                "component": "parser",
                "heap": options.heap,
                "java": options.java_cmd,
            },
            action=lambda: _run_parser_with_result(
                installation=tools["parser"],
                artifacts=artifacts,
                options=options,
            ),
        )
        if len(jobs) > job_count:
            artifacts.parser_result.unlink(missing_ok=True)

    if options.run_transchecker:
        job_count = len(jobs)
        append_job_if_needed(
            jobs,
            name="transchecker",
            outputs=[
                artifacts.transchecker_log,
                artifacts.aa_fasta,
                artifacts.nuc_fasta,
                artifacts.transchecker_result,
            ],
            dependencies=path_list(
                artifacts.ann_path,
                artifacts.fasta_path,
                tools["transchecker"].root,
                module_origin("msspack.validation"),
                module_origin("msspack.ddbj_tools"),
            ),
            cache_key={
                "component": "transchecker",
                "heap": options.heap,
                "java": options.java_cmd,
            },
            action=lambda: _run_transchecker_with_result(
                installation=tools["transchecker"],
                artifacts=artifacts,
                options=options,
            ),
        )
        if len(jobs) > job_count:
            artifacts.transchecker_result.unlink(missing_ok=True)

    return jobs


def run_validation(
    *,
    options: ValidationOptions,
    artifacts: ValidationArtifacts,
) -> dict[str, Path]:
    ensure_dir(artifacts.logs_dir)
    ensure_dir(artifacts.output_dir)
    failure_message: str | None = None
    jobs: list[NamedJob] | None = None
    try:
        jobs = _build_validation_jobs(options=options, artifacts=artifacts)
        run_named_jobs(
            jobs,
            parallel=options.parallel,
        )
    except Exception as exc:
        failure_message = str(exc) or type(exc).__name__
        raise
    finally:
        _write_validation_summary(
            options=options,
            artifacts=artifacts,
            failure_message=failure_message,
            ignore_existing_results=jobs is None,
        )
    return artifacts.result_paths(
        include_parser=options.run_parser,
        include_transchecker=options.run_transchecker,
    )


def validate_existing(
    *,
    config_file: str | Path | None,
    ann_path: str | Path,
    fasta_path: str | Path,
) -> dict[str, Path]:
    config = load_config(config_file) if config_file else None
    return run_validation(
        options=ValidationOptions.for_explicit_validate(config),
        artifacts=ValidationArtifacts.for_existing_submission(
            ann_path=Path(ann_path).expanduser().resolve(),
            fasta_path=Path(fasta_path).expanduser().resolve(),
        ),
    )
