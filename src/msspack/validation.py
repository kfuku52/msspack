from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from .config import MSSPackConfig, load_config
from .ddbj_tools import ensure_installed, run_parser, run_transchecker
from .execution import (
    NamedJob,
    append_job_if_needed,
    module_origin,
    path_list,
    run_named_jobs,
)
from .utils import ensure_dir


@dataclass(frozen=True)
class ValidationOptions:
    cache_dir: Optional[Path]
    heap: str
    parallel: bool
    java_cmd: str = "java"
    run_parser: bool = True
    run_transchecker: bool = True

    @classmethod
    def from_config(cls, config: Optional[MSSPackConfig]) -> "ValidationOptions":
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
    def for_explicit_validate(cls, config: Optional[MSSPackConfig]) -> "ValidationOptions":
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
    ) -> "ValidationArtifacts":
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

    def result_paths(
        self,
        *,
        include_parser: bool,
        include_transchecker: bool,
    ) -> dict[str, Path]:
        paths: dict[str, Path] = {}
        if include_parser:
            paths["parser_log"] = self.parser_log
        if include_transchecker:
            paths["transchecker_log"] = self.transchecker_log
            paths["aa_fasta"] = self.aa_fasta
            paths["nuc_fasta"] = self.nuc_fasta
        return paths


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

    tools = ensure_installed(
        requested_components,
        cache_dir=options.cache_dir,
    )
    jobs: list[NamedJob] = []

    if options.run_parser:
        append_job_if_needed(
            jobs,
            name="parser",
            outputs=[artifacts.parser_log],
            dependencies=path_list(
                artifacts.ann_path,
                artifacts.fasta_path,
                tools["parser"].root,
                module_origin("msspack.ddbj_tools"),
            ),
            cache_key={
                "component": "parser",
                "heap": options.heap,
                "java": options.java_cmd,
            },
            action=lambda: run_parser(
                tools["parser"],
                ann_path=artifacts.ann_path,
                fasta_path=artifacts.fasta_path,
                heap=options.heap,
                java_cmd=options.java_cmd,
                log_path=artifacts.parser_log,
            ),
        )

    if options.run_transchecker:
        append_job_if_needed(
            jobs,
            name="transchecker",
            outputs=[
                artifacts.transchecker_log,
                artifacts.aa_fasta,
                artifacts.nuc_fasta,
            ],
            dependencies=path_list(
                artifacts.ann_path,
                artifacts.fasta_path,
                tools["transchecker"].root,
                module_origin("msspack.ddbj_tools"),
            ),
            cache_key={
                "component": "transchecker",
                "heap": options.heap,
                "java": options.java_cmd,
            },
            action=lambda: run_transchecker(
                tools["transchecker"],
                ann_path=artifacts.ann_path,
                fasta_path=artifacts.fasta_path,
                aa_out=artifacts.aa_fasta,
                nuc_out=artifacts.nuc_fasta,
                heap=options.heap,
                java_cmd=options.java_cmd,
                log_path=artifacts.transchecker_log,
            ),
        )

    return jobs


def run_validation(
    *,
    options: ValidationOptions,
    artifacts: ValidationArtifacts,
) -> dict[str, Path]:
    run_named_jobs(
        _build_validation_jobs(options=options, artifacts=artifacts),
        parallel=options.parallel,
    )
    return artifacts.result_paths(
        include_parser=options.run_parser,
        include_transchecker=options.run_transchecker,
    )


def validate_existing(
    *,
    config_file: Optional[str | Path],
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
