#!/usr/bin/env python3

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DIRECTORY_PATTERNS = [
    "build",
    "dist",
    ".mypy_cache",
    ".ruff_cache",
    ".pytest_cache",
    ".cache-test",
    "busco_downloads",
    "**/__pycache__",
    "src/*.egg-info",
]
FILE_PATTERNS = [
    "busco_*.log",
]


def _unique_paths(patterns: list[str]) -> list[Path]:
    found: dict[Path, None] = {}
    for pattern in patterns:
        for path in ROOT.glob(pattern):
            resolved = path.resolve()
            try:
                relative = resolved.relative_to(ROOT)
            except ValueError:
                continue
            if relative.parts and relative.parts[0] == ".venv":
                continue
            found[resolved] = None
    return sorted(found.keys())


def _remove_path(path: Path) -> None:
    if path.is_dir() and not path.is_symlink():
        shutil.rmtree(path)
    else:
        path.unlink(missing_ok=True)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Remove repo-local build, cache, and BUSCO artifact files."
    )
    parser.add_argument("--dry-run", action="store_true", help="show what would be removed")
    args = parser.parse_args()

    targets = _unique_paths(DIRECTORY_PATTERNS) + _unique_paths(FILE_PATTERNS)
    if not targets:
        print("nothing to remove")
        return 0

    for path in targets:
        relative = path.relative_to(ROOT)
        if args.dry_run:
            print(relative)
            continue
        _remove_path(path)
        print(relative)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
