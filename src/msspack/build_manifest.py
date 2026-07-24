from __future__ import annotations

import hashlib
import json
import platform
import sys
import time
from collections.abc import Iterable
from dataclasses import dataclass, field
from datetime import UTC, datetime
from pathlib import Path

from . import __version__
from .config import MSSPackConfig
from .utils import ensure_dir, write_text


def _iso_utc_now() -> str:
    return datetime.now(UTC).isoformat()


def _file_info(path: Path) -> dict[str, object]:
    info: dict[str, object] = {"path": str(path)}
    if path.exists():
        stat = path.stat()
        info.update(
            {
                "exists": True,
                "size": stat.st_size,
                "mtime": datetime.fromtimestamp(stat.st_mtime, tz=UTC).isoformat(),
            }
        )
    else:
        info["exists"] = False
    return info


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _read_existing_manifest(path: Path) -> dict[str, object]:
    if not path.exists():
        return {}
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return {}
    return payload if isinstance(payload, dict) else {}


@dataclass
class StageRecord:
    name: str
    outputs: list[str]
    dependencies: list[str]
    ran: bool
    duration_seconds: float


@dataclass
class ManifestRecorder:
    config: MSSPackConfig
    config_path: Path
    output_root: Path
    ann_path: Path
    fasta_path: Path
    started_at: str = field(default_factory=_iso_utc_now)
    started_perf: float = field(default_factory=time.perf_counter, repr=False)
    stages: list[StageRecord] = field(default_factory=list)
    validation: dict[str, object] = field(default_factory=dict)
    status: str = "running"
    error: str | None = None

    def record_stage(
        self,
        *,
        name: str,
        outputs: Iterable[Path],
        dependencies: Iterable[Path],
        ran: bool,
        duration_seconds: float,
    ) -> None:
        self.stages.append(
            StageRecord(
                name=name,
                outputs=[str(path) for path in outputs],
                dependencies=[str(path) for path in dependencies],
                ran=ran,
                duration_seconds=duration_seconds,
            )
        )

    def set_validation(
        self,
        *,
        enabled: bool,
        result_paths: dict[str, Path] | None = None,
        options: dict[str, object] | None = None,
    ) -> None:
        payload: dict[str, object] = {"enabled": enabled}
        if options is not None:
            payload["options"] = options
        if result_paths is not None:
            payload["outputs"] = {key: str(path) for key, path in result_paths.items()}
        self.validation = payload

    def mark_failed(self, exc: Exception) -> None:
        self.status = "failed"
        self.error = str(exc)

    def mark_completed(self) -> None:
        self.status = "completed"
        self.error = None

    def write(self, path: Path) -> Path:
        ensure_dir(path.parent)
        stages_ran = sum(1 for stage in self.stages if stage.ran)
        stages_reused = len(self.stages) - stages_ran
        existing_payload = _read_existing_manifest(path)
        output_payload = {
            "root": str(self.output_root),
            "annotation": _file_info(self.ann_path),
            "fasta": _file_info(self.fasta_path),
        }
        if self.config.functional_annotation.enabled:
            output_payload["functional_annotation_taxonomy"] = _file_info(
                self.output_root / "final" / "functional-annotation-taxonomy.json"
            )
            output_payload["functional_annotation_evidence"] = _file_info(
                self.output_root / "final" / "functional-annotation.tsv"
            )
            output_payload["functional_annotation_name_standardization"] = _file_info(
                self.output_root
                / "final"
                / "functional-annotation-name-standardization.tsv"
            )
            output_payload["functional_domain_search_comparison"] = _file_info(
                self.output_root / "final" / "functional-domain-search-comparison.tsv"
            )
            if self.config.functional_annotation.consistency.enabled:
                output_payload["functional_annotation_consistency"] = _file_info(
                    self.output_root / "final" / "functional-annotation-consistency.tsv"
                )
                output_payload["functional_annotation_families"] = _file_info(
                    self.output_root / "final" / "functional-annotation-families.tsv"
                )
                output_payload["functional_annotation_conflicts"] = _file_info(
                    self.output_root / "final" / "functional-annotation-conflicts.tsv"
                )
                output_payload["functional_annotation_consistency_summary"] = _file_info(
                    self.output_root
                    / "final"
                    / "functional-annotation-consistency-summary.tsv"
                )
        payload = {
            "schema_version": 1,
            "msspack_version": __version__,
            "status": self.status,
            "error": self.error,
            "started_at": self.started_at,
            "finished_at": _iso_utc_now(),
            "duration_seconds": time.perf_counter() - self.started_perf,
            "runtime": {
                "python_version": sys.version.split()[0],
                "python_executable": sys.executable,
                "platform": platform.platform(),
            },
            "project": self.config.project.name,
            "databases": {
                "root": str(self.config.database_dir),
                "mode": (
                    "project"
                    if self.config.database_dir.is_relative_to(self.config.base_dir)
                    else "shared"
                ),
                "busco_root": str(self.config.busco_database_dir),
            },
            "config": {
                "path": str(self.config_path),
                "sha256": _sha256_file(self.config_path),
            },
            "inputs": {
                "fasta": _file_info(self.config.fasta_path),
                "gff": _file_info(self.config.gff_path),
            },
            "outputs": output_payload,
            "stages": [
                {
                    "name": stage.name,
                    "ran": stage.ran,
                    "duration_seconds": stage.duration_seconds,
                    "outputs": stage.outputs,
                    "dependencies": stage.dependencies,
                }
                for stage in self.stages
            ],
            "stage_summary": {
                "count": len(self.stages),
                "ran": stages_ran,
                "reused": stages_reused,
            },
            "validation": self.validation,
        }
        for key, value in existing_payload.items():
            payload.setdefault(key, value)
        write_text(path, json.dumps(payload, indent=2, sort_keys=True) + "\n")
        return path
