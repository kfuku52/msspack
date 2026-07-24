from __future__ import annotations

import json
import re
from dataclasses import asdict, dataclass
from pathlib import Path

from .config import MSSPackConfig
from .utils import ensure_dir, expand_path, write_text

_CDD_DATA_NAMES = (
    ("cddid.tbl", "cddid.tbl.gz"),
    ("cdtrack.txt", "cdtrack.txt"),
    ("family_superfamily_links", "family_superfamily_links"),
    ("cddannot.dat", "cddannot.dat.gz"),
    ("cddannot_generic.dat", "cddannot_generic.dat.gz"),
    ("bitscore_specific.txt", "bitscore_specific.txt"),
)


@dataclass(frozen=True)
class DatabaseResourceStatus:
    name: str
    path: Path
    enabled: bool
    ready: bool


@dataclass(frozen=True)
class DatabaseStatus:
    root: Path
    mode: str
    resources: tuple[DatabaseResourceStatus, ...]


def _has_files(path: Path) -> bool:
    return path.is_dir() and any(item.is_file() for item in path.rglob("*"))


def _diamond_database_ready(root: Path, name: str) -> bool:
    database_root = root / "diamond" / name
    return (database_root / f"{name}.dmnd").is_file() or any(
        path.is_file()
        for path in database_root.glob(f"*/{name}.dmnd")
    )


def _pfam_database_ready(root: Path) -> bool:
    candidates = [
        root / "pfam" / "Pfam-A.hmm",
        *(root / "pfam" / "objects").glob("*.hmm"),
    ]
    return any(
        path.is_file()
        and all(
            Path(str(path) + suffix).is_file()
            for suffix in (".h3f", ".h3i", ".h3m", ".h3p")
        )
        for path in candidates
    )


def _cdd_prefix_ready(prefix: Path) -> bool:
    if Path(str(prefix) + ".pal").is_file():
        return bool(
            tuple(prefix.parent.glob(prefix.name + ".*.aux"))
            and tuple(prefix.parent.glob(prefix.name + ".*.rps"))
        )
    return all(
        Path(str(prefix) + suffix).is_file()
        for suffix in (".aux", ".freq", ".loo", ".rps")
    )


def _cdd_directory_ready(directory: Path) -> bool:
    return any(
        _cdd_prefix_ready(candidate.with_suffix(""))
        for candidate in (
            *directory.rglob("Cdd.pal"),
            *directory.rglob("Cdd.aux"),
        )
    )


def _cdd_data_directory_ready(directory: Path) -> bool:
    return directory.is_dir() and all(
        (directory / local_name).is_file() or (directory / remote_name).is_file()
        for local_name, remote_name in _CDD_DATA_NAMES
    )


def _cdd_data_ready(config: MSSPackConfig, root: Path) -> bool:
    annotation = config.functional_annotation
    if annotation.cdd_data_dir.strip():
        return _cdd_data_directory_ready(
            expand_path(annotation.cdd_data_dir, config.base_dir)
        )
    cdd_root = root / "cdd"
    if _cdd_data_directory_ready(cdd_root / "data"):
        return True
    versions_root = cdd_root / "data-versions"
    return versions_root.is_dir() and any(
        _cdd_data_directory_ready(candidate)
        and (candidate / "provenance.json").is_file()
        for candidate in versions_root.iterdir()
        if candidate.is_dir()
    )


def _busco_lineage_ready(path: Path, *, explicit_lineage: bool) -> bool:
    if explicit_lineage:
        return (path / "dataset.cfg").is_file()
    return path.is_dir() and any(
        candidate.is_file() for candidate in path.glob("*/dataset.cfg")
    )


def collect_database_status(config: MSSPackConfig) -> DatabaseStatus:
    root = config.database_dir
    annotation = config.functional_annotation
    resources: list[DatabaseResourceStatus] = [
        DatabaseResourceStatus(
            "Swiss-Prot",
            root / "diamond" / "swissprot",
            annotation.enabled and annotation.swissprot_enabled,
            _diamond_database_ready(root, "swissprot"),
        )
    ]
    if annotation.reference_proteins.strip():
        safe_reference_name = re.sub(
            r"[^A-Za-z0-9_.-]+",
            "_",
            annotation.reference_name,
        ).strip("._") or "database"
        reference_path = (
            root
            / "diamond"
            / safe_reference_name
        )
        resources.append(
            DatabaseResourceStatus(
                annotation.reference_name,
                reference_path,
                annotation.enabled,
                _diamond_database_ready(root, safe_reference_name),
            )
        )
    cdd_path = root / "cdd"
    cdd_ready = _cdd_directory_ready(cdd_path)
    if annotation.cdd_database.strip():
        configured_cdd = expand_path(annotation.cdd_database, config.base_dir)
        cdd_prefix = configured_cdd / "Cdd" if configured_cdd.is_dir() else configured_cdd
        cdd_path = configured_cdd
        cdd_ready = (
            _cdd_directory_ready(configured_cdd)
            if configured_cdd.is_dir()
            else _cdd_prefix_ready(cdd_prefix)
        )
    cdd_ready = cdd_ready and _cdd_data_ready(config, root)
    resources.extend(
        [
            DatabaseResourceStatus(
                "UniRef90",
                root / "diamond" / "uniref90",
                annotation.enabled and annotation.uniref90_enabled,
                _diamond_database_ready(root, "uniref90"),
            ),
            DatabaseResourceStatus(
                "Pfam",
                root / "pfam",
                annotation.enabled and annotation.pfam_enabled,
                _pfam_database_ready(root),
            ),
            DatabaseResourceStatus(
                "CDD",
                cdd_path,
                annotation.enabled and annotation.cdd_enabled,
                cdd_ready,
            ),
            DatabaseResourceStatus(
                "taxonomy",
                root / "taxonomy",
                annotation.enabled and annotation.taxonomy.enabled,
                _has_files(root / "taxonomy"),
            ),
        ]
    )
    busco_lineage_path = (
        config.busco_database_dir / "lineages" / config.busco.lineage_dataset
        if config.busco.lineage_dataset
        else config.busco_database_dir / "lineages"
    )
    resources.append(
        DatabaseResourceStatus(
            (
                f"BUSCO {config.busco.lineage_dataset}"
                if config.busco.lineage_dataset
                else "BUSCO lineages"
            ),
            busco_lineage_path,
            config.busco.run_cds or config.busco.run_genome,
            _busco_lineage_ready(
                busco_lineage_path,
                explicit_lineage=bool(config.busco.lineage_dataset),
            ),
        )
    )
    try:
        root.relative_to(config.base_dir)
    except ValueError:
        mode = "shared"
    else:
        mode = "project"
    return DatabaseStatus(root=root, mode=mode, resources=tuple(resources))


def render_database_status(status: DatabaseStatus) -> str:
    lines = [
        f"Database root: {status.root}",
        f"Database mode: {status.mode}",
    ]
    for resource in status.resources:
        if resource.enabled:
            state = "ready" if resource.ready else "missing"
        else:
            state = "disabled (cached)" if resource.ready else "disabled"
        lines.append(f"{resource.name}: {state}\t{resource.path}")
    return "\n".join(lines)


def write_database_manifest(config: MSSPackConfig) -> Path:
    status = collect_database_status(config)
    ensure_dir(status.root)
    manifest_path = status.root / "manifest.json"
    payload = {
        "schema_version": 1,
        "root": str(status.root),
        "mode": status.mode,
        "resources": [
            {
                **asdict(resource),
                "path": str(resource.path),
            }
            for resource in status.resources
        ],
    }
    return write_text(
        manifest_path,
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
    )
