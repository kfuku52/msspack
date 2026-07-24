from __future__ import annotations

import json
import os
import shutil
import tempfile
import time
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path

from .busco import BuscoArtifacts, run_busco_comparison
from .config import load_config
from .databases import DatabaseStatus, collect_database_status, write_database_manifest
from .doctor import doctor_succeeded, render_doctor_report, run_doctor
from .pipeline import PipelineOutputs, run_pipeline
from .pipeline_plots import PipelinePlotArtifacts, run_pipeline_plots
from .report import ReportArtifacts, run_html_report
from .utils import MSSPackError, write_text


@dataclass(frozen=True)
class RunArtifacts:
    pipeline: PipelineOutputs
    busco: BuscoArtifacts | None
    plots: PipelinePlotArtifacts
    report: ReportArtifacts | None
    database_status: DatabaseStatus
    duration_seconds: float


@contextmanager
def database_directory_override(path: str = "") -> Iterator[None]:
    previous = os.environ.get("MSSPACK_DB_DIR")
    if path.strip():
        os.environ["MSSPACK_DB_DIR"] = path
    try:
        yield
    finally:
        if path.strip():
            if previous is None:
                os.environ.pop("MSSPACK_DB_DIR", None)
            else:
                os.environ["MSSPACK_DB_DIR"] = previous


def _invalidate_compute_cache(
    output_root: Path,
    *,
    protected_paths: tuple[Path, ...],
) -> int:
    resolved_root = output_root.expanduser().resolve()
    filesystem_root = Path(resolved_root.anchor)
    unsafe_roots = {
        filesystem_root,
        Path.home().resolve(),
        Path(tempfile.gettempdir()).resolve(),
    }
    if resolved_root in unsafe_roots or len(resolved_root.parts) <= 2:
        raise MSSPackError(
            f"Refusing to recursively clear compute caches from an unsafe output root: "
            f"{resolved_root}"
        )
    for protected_path in protected_paths:
        resolved_protected = protected_path.expanduser().resolve()
        if resolved_root == resolved_protected or resolved_root in resolved_protected.parents:
            raise MSSPackError(
                "Refusing to clear compute caches because the output root contains a "
                f"protected path: output={resolved_root}, protected={resolved_protected}"
            )

    removed = 0
    if not resolved_root.is_dir():
        return removed
    for cache_dir in sorted(resolved_root.rglob(".msspack-cache"), reverse=True):
        if cache_dir.is_symlink():
            cache_dir.unlink()
        elif cache_dir.is_dir():
            shutil.rmtree(cache_dir)
        else:
            continue
        removed += 1
    return removed


def _preflight(
    config_file: str | Path,
    *,
    validate: bool,
    run_busco: bool,
) -> None:
    config = load_config(config_file)
    checks = run_doctor(config)
    if not validate:
        validation_names = {
            "DDBJ validation platform",
            "java",
            "bash",
            "DDBJ parser",
            "DDBJ transchecker",
        }
        for check in checks:
            if check.name in validation_names:
                check.required = False
    busco_names = {"BUSCO (optional)", "BUSCO database root"}
    for check in checks:
        if check.name in busco_names:
            check.required = run_busco
    if not doctor_succeeded(checks):
        raise MSSPackError(
            "Preflight checks failed:\n" + render_doctor_report(config, checks=checks)
        )


def _update_run_manifest(
    path: Path,
    *,
    status: str,
    started_at: datetime,
    duration_seconds: float,
    database_status: DatabaseStatus,
    busco_enabled: bool,
    validation_enabled: bool,
    report_enabled: bool,
    force_compute: bool,
    error: str = "",
) -> None:
    payload: dict[str, object]
    if path.is_file():
        loaded = json.loads(path.read_text(encoding="utf-8"))
        payload = loaded if isinstance(loaded, dict) else {}
    else:
        payload = {}
    payload["run"] = {
        "status": status,
        "started_at": started_at.isoformat(),
        "finished_at": (
            datetime.now(UTC).isoformat() if status in {"completed", "failed"} else None
        ),
        "duration_seconds": duration_seconds,
        "database_root": str(database_status.root),
        "database_mode": database_status.mode,
        "busco_enabled": busco_enabled,
        "validation_enabled": validation_enabled,
        "report_enabled": report_enabled,
        "force_compute": force_compute,
        "error": error or None,
    }
    write_text(path, json.dumps(payload, indent=2, sort_keys=True) + "\n")


def run_all(
    config_file: str | Path,
    *,
    database_dir: str = "",
    force_compute: bool = False,
    run_busco: bool = True,
    validate: bool = True,
    write_report: bool = True,
) -> RunArtifacts:
    started_at = datetime.now(UTC)
    started = time.perf_counter()
    with database_directory_override(database_dir):
        config = load_config(config_file)
        database_status = collect_database_status(config)
        manifest_path = config.output_dir / "build-manifest.json"
        effective_busco = run_busco and (config.busco.run_cds or config.busco.run_genome)
        effective_validation = validate and (
            config.pipeline.validate_with_parser
            or config.pipeline.validate_with_transchecker
        )
        if force_compute:
            _invalidate_compute_cache(
                config.output_dir,
                protected_paths=(
                    config.base_dir,
                    config.database_dir,
                    config.busco_database_dir,
                    config.cache_dir,
                ),
            )
        _update_run_manifest(
            manifest_path,
            status="running",
            started_at=started_at,
            duration_seconds=0.0,
            database_status=database_status,
            busco_enabled=effective_busco,
            validation_enabled=effective_validation,
            report_enabled=write_report,
            force_compute=force_compute,
        )
        try:
            _preflight(
                config_file,
                validate=effective_validation,
                run_busco=effective_busco,
            )

            busco_artifacts: BuscoArtifacts | None = None
            if effective_busco:
                busco_artifacts = run_busco_comparison(
                    config_file,
                    force=force_compute,
                    prepare_only=True,
                )

            pipeline_outputs = run_pipeline(config_file, validate=validate)
            plot_artifacts = run_pipeline_plots(config_file, force=force_compute)
            database_status = collect_database_status(config)
            write_database_manifest(config)
            report_artifacts = (
                run_html_report(config_file, force=False) if write_report else None
            )
            duration_seconds = time.perf_counter() - started
            _update_run_manifest(
                pipeline_outputs.manifest_path,
                status="completed",
                started_at=started_at,
                duration_seconds=duration_seconds,
                database_status=database_status,
                busco_enabled=effective_busco,
                validation_enabled=effective_validation,
                report_enabled=write_report,
                force_compute=force_compute,
            )
            return RunArtifacts(
                pipeline=pipeline_outputs,
                busco=busco_artifacts,
                plots=plot_artifacts,
                report=report_artifacts,
                database_status=database_status,
                duration_seconds=duration_seconds,
            )
        except BaseException as exc:
            try:
                _update_run_manifest(
                    manifest_path,
                    status="failed",
                    started_at=started_at,
                    duration_seconds=time.perf_counter() - started,
                    database_status=database_status,
                    busco_enabled=effective_busco,
                    validation_enabled=effective_validation,
                    report_enabled=write_report,
                    force_compute=force_compute,
                    error=str(exc) or type(exc).__name__,
                )
            except Exception:
                pass
            raise
