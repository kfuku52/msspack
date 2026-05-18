from __future__ import annotations

import importlib.util
from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import lru_cache
from pathlib import Path
from typing import Callable

from .utils import MSSPackError

NamedJob = tuple[str, Callable[[], object]]


@lru_cache(maxsize=None)
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
        if isinstance(item, (list, tuple, set)):
            paths.extend(path_list(*item))
        elif isinstance(item, Path):
            paths.append(item)
        else:
            paths.append(Path(str(item)))
    return paths


def is_up_to_date(outputs: list[Path], dependencies: list[Path]) -> bool:
    if not outputs:
        return False
    if any(not path.exists() for path in outputs):
        return False
    if not dependencies:
        return True
    existing_deps = [path for path in dependencies if path.exists()]
    if not existing_deps:
        return True
    newest_input = max(path.stat().st_mtime for path in existing_deps)
    oldest_output = min(path.stat().st_mtime for path in outputs)
    return oldest_output >= newest_input


def run_if_needed(
    *,
    outputs: list[Path],
    dependencies: list[Path],
    action: Callable[[], object],
) -> bool:
    if is_up_to_date(outputs, dependencies):
        return False
    action()
    return True


def append_job_if_needed(
    jobs: list[NamedJob],
    *,
    name: str,
    outputs: list[Path],
    dependencies: list[Path],
    action: Callable[[], object],
) -> None:
    if not is_up_to_date(outputs, dependencies):
        jobs.append((name, action))


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
        future_to_name = {
            executor.submit(action): name
            for name, action in jobs
        }
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
