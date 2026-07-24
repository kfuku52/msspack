#!/usr/bin/env python3

from __future__ import annotations

import argparse
import tarfile
import zipfile
from pathlib import Path


def _wheel_members(path: Path) -> set[str]:
    with zipfile.ZipFile(path) as archive:
        return set(archive.namelist())


def _sdist_members(path: Path) -> set[str]:
    with tarfile.open(path, "r:gz") as archive:
        return set(archive.getnames())


def verify_distribution(path: Path) -> None:
    if path.suffix == ".whl":
        members = _wheel_members(path)
        required_suffixes = (
            "LICENSE",
            "THIRD_PARTY_NOTICES.md",
            "msspack/annotation_consistency.py",
            "msspack/annotation_consistency_plots.py",
            "msspack/annotation_taxonomy.py",
            "msspack/database_lock.py",
            "msspack/databases.py",
            "msspack/functional_annotation.py",
            "msspack/product_names.py",
            "msspack/py.typed",
            "msspack/demo_data/README.md",
            "msspack/demo_data/annotation.gff3",
            "msspack/demo_data/config.functional.toml",
            "msspack/demo_data/config.toml",
            "msspack/demo_data/expected-summary.json",
            "msspack/demo_data/genome.fa",
            "msspack/demo_data/reference.faa",
            "msspack/templates/msspack.example.toml",
            "msspack/workflow.py",
        )
    elif path.name.endswith(".tar.gz"):
        members = _sdist_members(path)
        required_suffixes = (
            "LICENSE",
            "THIRD_PARTY_NOTICES.md",
            "CHANGELOG.md",
            "CONTRIBUTING.md",
            "RELEASE.md",
            "docs/assets/msspack-header.png",
            "docs/assets/sample-busco-cds-comparison.svg",
            "docs/assets/sample-pipeline-gene-flow.sankey.svg",
            "examples/msspack.example.toml",
            "scripts/check_distribution.py",
            "src/msspack/annotation_consistency.py",
            "src/msspack/annotation_consistency_plots.py",
            "src/msspack/annotation_taxonomy.py",
            "src/msspack/database_lock.py",
            "src/msspack/databases.py",
            "src/msspack/functional_annotation.py",
            "src/msspack/product_names.py",
            "src/msspack/workflow.py",
            "src/msspack/demo_data/README.md",
            "src/msspack/demo_data/annotation.gff3",
            "src/msspack/demo_data/config.functional.toml",
            "src/msspack/demo_data/config.toml",
            "src/msspack/demo_data/expected-summary.json",
            "src/msspack/demo_data/genome.fa",
            "src/msspack/demo_data/reference.faa",
            "tests/test_annotation_consistency.py",
            "tests/test_annotation_taxonomy.py",
            "tests/test_database_lock.py",
            "tests/test_databases.py",
            "tests/test_functional_annotation.py",
            "tests/test_product_names.py",
            "tests/test_workflow.py",
            "tests/fixtures/minimal_pack/config.toml",
            "tests/fixtures/minimal_pack/expected.ann.txt",
            "tests/fixtures/minimal_pack/expected.fasta",
            "tests/fixtures/minimal_pack/input.fa",
            "tests/fixtures/minimal_pack/input.gff3",
        )
    else:
        raise ValueError(f"Unsupported distribution file: {path}")
    missing = [
        suffix
        for suffix in required_suffixes
        if not any(member.endswith("/" + suffix) or member == suffix for member in members)
    ]
    if missing:
        raise ValueError(f"{path} is missing required files: {', '.join(missing)}")


def main() -> int:
    parser = argparse.ArgumentParser(description="Verify built msspack distributions")
    parser.add_argument("paths", nargs="*", type=Path, help="wheel/sdist files")
    args = parser.parse_args()
    paths = args.paths or sorted(Path("dist").glob("msspack-*"))
    if not paths:
        raise ValueError("No distribution files were found")
    for path in paths:
        verify_distribution(path)
        print(f"verified {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
