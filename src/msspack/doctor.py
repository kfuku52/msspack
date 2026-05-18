from __future__ import annotations

import importlib.util
from dataclasses import dataclass
from typing import List, Optional

from .config import MSSPackConfig
from .ddbj_tools import list_installed
from .utils import which


@dataclass
class Check:
    name: str
    ok: bool
    detail: str


def _importable(module_name: str) -> bool:
    return importlib.util.find_spec(module_name) is not None


def run_doctor(config: Optional[MSSPackConfig] = None) -> List[Check]:
    checks: List[Check] = []

    java_cmd = config.tools.java if config else "java"
    busco_cmd = config.busco.command if config else "busco"

    checks.append(Check("java", which(java_cmd) is not None, java_cmd))
    checks.append(Check("BUSCO (optional)", which(busco_cmd) is not None, busco_cmd))
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
        checks.append(
            Check(
                f"DDBJ {component}",
                item is not None,
                str(item.root) if item else "not installed",
            )
        )

    ume = installed.get("ume")
    checks.append(
        Check(
            "DDBJ ume (optional)",
            True,
            str(ume.root) if ume else "not installed",
        )
    )

    return checks


def render_doctor_report(config: Optional[MSSPackConfig] = None) -> str:
    checks = run_doctor(config)
    lines = []
    for check in checks:
        status = "OK" if check.ok else "MISSING"
        lines.append(f"[{status}] {check.name}: {check.detail}")
    return "\n".join(lines)
