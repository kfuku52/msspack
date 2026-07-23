from __future__ import annotations

import importlib.resources
from pathlib import Path

from .utils import MSSPackError, atomic_binary_writer, ensure_dir

DEMO_FILENAMES = (
    "README.md",
    "annotation.gff3",
    "config.functional.toml",
    "config.toml",
    "expected-summary.json",
    "genome.fa",
    "reference.faa",
)


def write_demo_dataset(output: str | Path, *, force: bool = False) -> Path:
    destination = Path(output).expanduser().resolve()
    if destination.exists() and not destination.is_dir():
        raise MSSPackError(f"Demo output exists and is not a directory: {destination}")

    conflicts = tuple(
        destination / filename
        for filename in DEMO_FILENAMES
        if (destination / filename).exists()
    )
    if conflicts and not force:
        names = ", ".join(path.name for path in conflicts)
        raise MSSPackError(
            f"Refusing to overwrite demo files in {destination}: {names}. "
            "Pass --force to replace the bundled files."
        )
    directory_conflicts = tuple(path for path in conflicts if path.is_dir())
    if directory_conflicts:
        names = ", ".join(path.name for path in directory_conflicts)
        raise MSSPackError(f"Demo file targets are directories in {destination}: {names}")

    source_root = importlib.resources.files("msspack").joinpath("demo_data")
    ensure_dir(destination)
    for filename in DEMO_FILENAMES:
        source = source_root.joinpath(filename)
        if not source.is_file():
            raise MSSPackError(f"Bundled demo file is missing: {filename}")
        with atomic_binary_writer(destination / filename) as output_handle:
            output_handle.write(source.read_bytes())
    return destination


__all__ = ["DEMO_FILENAMES", "write_demo_dataset"]
