"""Determine whether retained BUSCO artifacts describe the current inputs."""
from __future__ import annotations

import hashlib
import json
from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import Any


def busco_input_provenance(paths: Iterable[Path]) -> dict[str, str]:
    result = {}
    for path in paths:
        with path.open("rb") as handle:
            digest = hashlib.file_digest(handle, "sha256").hexdigest()
        result[str(path.resolve())] = digest
    return result


def available_busco_comparisons(
    output_root: Path, payload: Mapping[str, Any] | None = None,
) -> frozenset[str]:
    legacy_comparisons = frozenset({"cds", "genome"})
    if payload is None:
        manifest = output_root / "build-manifest.json"
        if not manifest.is_file():
            return legacy_comparisons  # Standalone comparisons from older versions.
        try:
            payload = json.loads(manifest.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            return frozenset()
        if not isinstance(payload, dict):
            return frozenset()
    busco = payload.get("busco")
    run = payload.get("run")
    if not isinstance(busco, dict):
        return (
            frozenset() if isinstance(run, dict) and run.get("busco_enabled") is False
            else legacy_comparisons
        )
    comparisons = busco.get("comparisons")
    active = frozenset(comparisons) if isinstance(comparisons, dict) else legacy_comparisons
    if not busco.get("enabled"):
        return frozenset()
    provenance = busco.get("input_provenance")
    if provenance is None:
        return (
            frozenset() if isinstance(run, dict) and run.get("busco_enabled") is False else active
        )
    if not isinstance(provenance, dict) or not provenance:
        return frozenset()
    if any(not isinstance(path, str) or not isinstance(digest, str)
           for path, digest in provenance.items()):
        return frozenset()
    try:
        return active if busco_input_provenance(Path(path) for path in provenance) == provenance else frozenset()
    except OSError:
        return frozenset()


def busco_results_available(
    output_root: Path, payload: Mapping[str, Any] | None = None,
) -> bool:
    return bool(available_busco_comparisons(output_root, payload))
