from __future__ import annotations

import hashlib
import importlib.util
import json
from collections.abc import Callable
from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import cache
from pathlib import Path

from .utils import MSSPackError, write_text

NamedJob = tuple[str, Callable[[], object]]
CACHE_SCHEMA_VERSION = 2
_FILE_FINGERPRINT_CACHE: dict[tuple[str, int, int, int], dict[str, object]] = {}


@cache
def module_origin(module_name: str) -> Path:
    spec = importlib.util.find_spec(module_name)
    if spec is None or spec.origin is None:
        raise MSSPackError(f"Could not resolve module path: {module_name}")
    return Path(spec.origin).resolve()


def path_list(*items: object) -> list[Path]:
    paths: list[Path] = []
    for item in items:
        if item is None:
            continue
        if isinstance(item, list | tuple | set):
            paths.extend(path_list(*item))
        elif isinstance(item, Path):
            paths.append(item)
        else:
            paths.append(Path(str(item)))
    return paths


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _fingerprint(path: Path) -> dict[str, object] | None:
    if not path.exists():
        return None
    resolved = path.resolve()
    if resolved.is_file():
        stat = resolved.stat()
        cache_identity = (
            str(resolved),
            stat.st_size,
            stat.st_mtime_ns,
            stat.st_ctime_ns,
        )
        cached = _FILE_FINGERPRINT_CACHE.get(cache_identity)
        if cached is not None:
            return cached
        result: dict[str, object] = {
            "kind": "file",
            "size": stat.st_size,
            "sha256": _sha256_file(resolved),
        }
        _FILE_FINGERPRINT_CACHE[cache_identity] = result
        return result
    if resolved.is_dir():
        entries: list[dict[str, object]] = []
        for child in sorted(item for item in resolved.rglob("*") if item.is_file()):
            stat = child.stat()
            entries.append(
                {
                    "path": str(child.relative_to(resolved)),
                    "size": stat.st_size,
                    "sha256": _sha256_file(child),
                }
            )
        return {"kind": "directory", "entries": entries}
    return {"kind": "other"}


def _fingerprints(paths: list[Path]) -> dict[str, dict[str, object]] | None:
    fingerprints: dict[str, dict[str, object]] = {}
    for path in paths:
        fingerprint = _fingerprint(path)
        if fingerprint is None:
            return None
        fingerprints[str(path.resolve())] = fingerprint
    return fingerprints


def _normalized_cache_key(cache_key: object | None) -> object | None:
    if cache_key is None:
        return None
    try:
        normalized: object = json.loads(json.dumps(cache_key, sort_keys=True))
        return normalized
    except (TypeError, ValueError) as exc:
        raise MSSPackError(f"Cache key is not JSON serializable: {cache_key!r}") from exc


def _cache_state_path(outputs: list[Path]) -> Path:
    identity = "\0".join(str(path.resolve()) for path in outputs)
    digest = hashlib.sha256(identity.encode("utf-8")).hexdigest()
    return outputs[0].parent / ".msspack-cache" / f"{digest}.json"


def _read_cache_state(path: Path) -> dict[str, object] | None:
    if not path.exists():
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    if not isinstance(payload, dict) or payload.get("schema_version") != CACHE_SCHEMA_VERSION:
        return None
    return payload


def is_up_to_date(
    outputs: list[Path],
    dependencies: list[Path],
    *,
    cache_key: object | None = None,
) -> bool:
    if not outputs:
        return False
    state = _read_cache_state(_cache_state_path(outputs))
    if state is None:
        return False
    output_fingerprints = _fingerprints(outputs)
    dependency_fingerprints = _fingerprints(dependencies)
    if output_fingerprints is None or dependency_fingerprints is None:
        return False
    return (
        state.get("outputs") == output_fingerprints
        and state.get("dependencies") == dependency_fingerprints
        and state.get("cache_key") == _normalized_cache_key(cache_key)
    )


def _write_cache_state(
    outputs: list[Path],
    dependencies: list[Path],
    *,
    cache_key: object | None,
) -> None:
    output_fingerprints = _fingerprints(outputs)
    dependency_fingerprints = _fingerprints(dependencies)
    if output_fingerprints is None:
        missing = [str(path) for path in outputs if not path.exists()]
        raise MSSPackError("Action did not create all declared outputs: " + ", ".join(missing))
    if dependency_fingerprints is None:
        missing = [str(path) for path in dependencies if not path.exists()]
        raise MSSPackError("Action dependencies disappeared: " + ", ".join(missing))
    payload = {
        "schema_version": CACHE_SCHEMA_VERSION,
        "outputs": output_fingerprints,
        "dependencies": dependency_fingerprints,
        "cache_key": _normalized_cache_key(cache_key),
    }
    state_path = _cache_state_path(outputs)
    write_text(state_path, json.dumps(payload, indent=2, sort_keys=True) + "\n")


def run_if_needed(
    *,
    outputs: list[Path],
    dependencies: list[Path],
    action: Callable[[], object],
    cache_key: object | None = None,
) -> bool:
    if not outputs:
        raise MSSPackError("Cached actions must declare at least one output")
    if is_up_to_date(outputs, dependencies, cache_key=cache_key):
        return False
    state_path = _cache_state_path(outputs) if outputs else None
    if state_path is not None:
        state_path.unlink(missing_ok=True)
    dependency_fingerprints_before = _fingerprints(dependencies)
    if dependency_fingerprints_before is None:
        missing = [str(path) for path in dependencies if not path.exists()]
        raise MSSPackError("Action dependencies are missing: " + ", ".join(missing))
    action()
    dependency_fingerprints_after = _fingerprints(dependencies)
    if dependency_fingerprints_before != dependency_fingerprints_after:
        raise MSSPackError("Action dependencies changed while the step was running; rerun the step")
    _write_cache_state(outputs, dependencies, cache_key=cache_key)
    return True


def append_job_if_needed(
    jobs: list[NamedJob],
    *,
    name: str,
    outputs: list[Path],
    dependencies: list[Path],
    action: Callable[[], object],
    cache_key: object | None = None,
) -> None:
    if is_up_to_date(outputs, dependencies, cache_key=cache_key):
        return

    def cached_action() -> object:
        run_if_needed(
            outputs=outputs,
            dependencies=dependencies,
            action=action,
            cache_key=cache_key,
        )
        return None

    jobs.append((name, cached_action))


def run_named_jobs(jobs: list[NamedJob], *, parallel: bool) -> None:
    if not jobs:
        return
    if not parallel or len(jobs) == 1:
        for _, action in jobs:
            action()
        return

    failures: list[tuple[str, Exception]] = []
    with ThreadPoolExecutor(
        max_workers=len(jobs),
        thread_name_prefix="msspack-validate",
    ) as executor:
        future_to_name = {executor.submit(action): name for name, action in jobs}
        for future in as_completed(future_to_name):
            try:
                future.result()
            except Exception as exc:  # pragma: no cover - exercised via callers
                failures.append((future_to_name[future], exc))

    if failures:
        _, first_error = failures[0]
        if isinstance(first_error, MSSPackError):
            raise first_error
        raise MSSPackError(f"Validation job failed: {failures[0][0]}") from first_error
