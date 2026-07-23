from __future__ import annotations

import gzip
import importlib.util
import platform
from dataclasses import dataclass

from .config import MSSPackConfig
from .ddbj_tools import list_installed
from .fasta import iter_fasta
from .gff import GFFRecord
from .utils import expand_path, which


@dataclass
class Check:
    name: str
    ok: bool
    detail: str
    required: bool = True


def _importable(module_name: str) -> bool:
    return importlib.util.find_spec(module_name) is not None


def _input_checks(config: MSSPackConfig) -> list[Check]:
    checks: list[Check] = []
    fasta_path = config.fasta_path
    gff_path = config.gff_path
    checks.append(Check("input FASTA", fasta_path.is_file(), str(fasta_path)))
    checks.append(Check("input GFF", gff_path.is_file(), str(gff_path)))
    if not (fasta_path.is_file() and gff_path.is_file()):
        return checks
    try:
        fasta_ids: set[str] = set()
        duplicates: set[str] = set()
        empty_records: set[str] = set()
        for fasta_record in iter_fasta(fasta_path):
            if fasta_record.id in fasta_ids:
                duplicates.add(fasta_record.id)
            fasta_ids.add(fasta_record.id)
            if not fasta_record.sequence:
                empty_records.add(fasta_record.id)
        checks.append(
            Check(
                "FASTA records",
                bool(fasta_ids),
                f"{len(fasta_ids)} unique record(s)" if fasta_ids else "no sequence records",
            )
        )
        checks.append(
            Check(
                "unique FASTA IDs",
                not duplicates,
                "unique" if not duplicates else ", ".join(sorted(duplicates)[:10]),
            )
        )
        checks.append(
            Check(
                "non-empty FASTA records",
                not empty_records,
                "all records contain sequence data"
                if not empty_records
                else ", ".join(sorted(empty_records)[:10]),
            )
        )
        gff_seqids: set[str] = set()
        if gff_path.suffix == ".gz":
            gff_handle = gzip.open(gff_path, "rt", encoding="utf-8")
        else:
            gff_handle = gff_path.open("r", encoding="utf-8")
        with gff_handle as handle:
            for line_number, raw_line in enumerate(handle, start=1):
                line = raw_line.rstrip("\n")
                if line == "##FASTA":
                    break
                if not line or line.startswith("#"):
                    continue
                fields = line.split("\t")
                if len(fields) != 9:
                    raise ValueError(f"line {line_number} has {len(fields)} columns")
                gff_record = GFFRecord.from_line(line)
                if gff_record.start < 1 or gff_record.end < gff_record.start:
                    raise ValueError(
                        f"line {line_number} has invalid coordinates "
                        f"{gff_record.start}..{gff_record.end}"
                    )
                gff_seqids.add(gff_record.seqid)
        missing = sorted(gff_seqids - fasta_ids)
        checks.append(
            Check(
                "GFF/FASTA seqids",
                not missing,
                "all GFF seqids are present"
                if not missing
                else "missing from FASTA: " + ", ".join(missing[:10]),
            )
        )
    except (OSError, UnicodeError, ValueError) as exc:
        checks.append(Check("input format", False, str(exc)))
    return checks


def run_doctor(config: MSSPackConfig | None = None) -> list[Check]:
    checks: list[Check] = []

    java_cmd = config.tools.java if config else "java"
    busco_cmd = config.busco.command if config else "busco"
    validation_required = config is None or (
        config.pipeline.validate_with_parser or config.pipeline.validate_with_transchecker
    )

    supported_platform = platform.system() != "Windows"
    checks.append(
        Check(
            "DDBJ validation platform",
            supported_platform,
            "Linux/macOS"
            if supported_platform
            else "native Windows is unsupported; use WSL or DDBJ Windows tools separately",
            required=validation_required,
        )
    )

    checks.append(
        Check("java", which(java_cmd) is not None, java_cmd, required=validation_required)
    )
    checks.append(
        Check("bash", which("bash") is not None, "bash", required=validation_required)
    )
    checks.append(
        Check("BUSCO (optional)", which(busco_cmd) is not None, busco_cmd, required=False)
    )
    annotation = config.functional_annotation if config is not None else None
    annotation_enabled = bool(annotation and annotation.enabled)
    diamond_required = bool(
        annotation_enabled
        and annotation
        and (
            annotation.swissprot_enabled
            or annotation.uniref90_enabled
            or annotation.reference_proteins.strip()
            or annotation.consistency.enabled
        )
    )
    diamond_cmd = annotation.diamond_command if annotation else "diamond"
    checks.append(
        Check(
            "DIAMOND (functional annotation)",
            which(diamond_cmd) is not None,
            diamond_cmd,
            required=diamond_required,
        )
    )
    pfam_required = bool(annotation_enabled and annotation and annotation.pfam_enabled)
    hmmscan_cmd = annotation.hmmscan_command if annotation else "hmmscan"
    hmmpress_cmd = annotation.hmmpress_command if annotation else "hmmpress"
    checks.extend(
        [
            Check(
                "HMMER hmmscan (Pfam fallback)",
                which(hmmscan_cmd) is not None,
                hmmscan_cmd,
                required=pfam_required,
            ),
            Check(
                "HMMER hmmpress (Pfam fallback)",
                which(hmmpress_cmd) is not None,
                hmmpress_cmd,
                required=pfam_required,
            ),
        ]
    )
    cdd_required = bool(annotation_enabled and annotation and annotation.cdd_enabled)
    rpsblast_cmd = annotation.rpsblast_command if annotation else "rpsblast"
    rpsbproc_cmd = annotation.rpsbproc_command if annotation else "rpsbproc"
    checks.extend(
        [
            Check(
                "NCBI RPS-BLAST (CDD fallback)",
                which(rpsblast_cmd) is not None,
                rpsblast_cmd,
                required=cdd_required,
            ),
            Check(
                "NCBI rpsbproc (CDD fallback)",
                which(rpsbproc_cmd) is not None,
                rpsbproc_cmd,
                required=cdd_required,
            ),
        ]
    )
    converter_modules = {
        "Bio": "biopython",
    }
    missing = [
        label for module_name, label in converter_modules.items() if not _importable(module_name)
    ]
    checks.append(
        Check(
            "msspack converter deps",
            not missing,
            ", ".join(missing) if missing else "bundled converter and GFF sorter ready",
        )
    )

    installed = list_installed(config.cache_dir if config else None)
    for component in ("parser", "transchecker"):
        item = installed.get(component)
        required = validation_required and (
            (component == "parser" and (config is None or config.pipeline.validate_with_parser))
            or (
                component == "transchecker"
                and (config is None or config.pipeline.validate_with_transchecker)
            )
        )
        checks.append(
            Check(
                f"DDBJ {component}",
                item is not None,
                str(item.root) if item else "not installed",
                required=required,
            )
        )

    ume = installed.get("ume")
    checks.append(
        Check(
            "DDBJ ume (optional)",
            ume is not None,
            str(ume.root) if ume else "not installed",
            required=False,
        )
    )

    if config is not None:
        checks.extend(_input_checks(config))
        annotation = config.functional_annotation
        if annotation.enabled:
            local_databases = (
                ("Swiss-Prot FASTA", annotation.swissprot_fasta),
                ("UniRef90 FASTA", annotation.uniref90_fasta),
                ("reference protein FASTA", annotation.reference_proteins),
                ("Pfam HMM", annotation.pfam_hmm),
            )
            for label, configured_path in local_databases:
                if not configured_path.strip():
                    continue
                resolved = expand_path(configured_path, config.base_dir)
                checks.append(
                    Check(
                        f"functional annotation {label}",
                        resolved.is_file(),
                        str(resolved),
                    )
                )
            for label, configured_path in (
                ("CDD database", annotation.cdd_database),
                ("CDD data directory", annotation.cdd_data_dir),
            ):
                if not configured_path.strip():
                    continue
                resolved = expand_path(configured_path, config.base_dir)
                exists = (
                    (resolved / "Cdd.aux").is_file()
                    if label == "CDD database" and resolved.is_dir()
                    else resolved.exists()
                )
                checks.append(
                    Check(
                        f"functional annotation {label}",
                        exists,
                        str(resolved),
                    )
                )

    return checks


def doctor_succeeded(checks: list[Check]) -> bool:
    return all(check.ok or not check.required for check in checks)


def render_doctor_report(
    config: MSSPackConfig | None = None,
    *,
    checks: list[Check] | None = None,
) -> str:
    checks = checks if checks is not None else run_doctor(config)
    lines = []
    for check in checks:
        status = "OK" if check.ok else ("MISSING" if check.required else "OPTIONAL")
        lines.append(f"[{status}] {check.name}: {check.detail}")
    return "\n".join(lines)
