from __future__ import annotations

import hashlib
import json
import os
import re
import shlex
import shutil
import stat
import tarfile
import urllib.request
import zipfile
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from tempfile import TemporaryDirectory, mkstemp
from typing import Dict, Iterable, Optional, Tuple

from .utils import MSSPackError, default_cache_dir, ensure_dir, run_command, write_text

DDBJ_TOOL_INDEX = "https://ddbj.nig.ac.jp/public/ddbj-cib/MSS/"
TRUSTED_ARCHIVE_SHA256 = {
    "UME_unix_V1.66.zip": "62342c07396ee8670486e5d3b7043d839de3ad4158de5a893908eafed2e1c351",
    "Parser_V6.80.tar.gz": "c3ec1cf9f90e5dcf647be42e9cd46f30b83e4a7b58c1dd30b7acafc7b94fdd64",
    "transChecker_V2.26.tar.gz": "2a4a68a44ef4f4a81ba7c72c516e2a7fec2922e2d43a18672d9abcc115657f21",
}


@dataclass(frozen=True)
class ToolInstallation:
    component: str
    version: str
    archive_name: str
    root: Path
    metadata: Optional[dict[str, object]] = None

    @property
    def executable(self) -> Path:
        if self.component == "parser":
            return self.root / "jParser.sh"
        if self.component == "transchecker":
            return self.root / "transChecker.sh"
        if self.component == "ume":
            return self.root / "UME.sh"
        raise MSSPackError(f"Unsupported tool component: {self.component}")

    @property
    def metadata_path(self) -> Path:
        return self.root / ".msspack-install.json"


_PATTERNS = {
    "ume": re.compile(r"(UME_unix_V(?P<version>[\d.]+)\.zip)"),
    "parser": re.compile(r"(Parser_V(?P<version>[\d.]+)\.tar\.gz)"),
    "transchecker": re.compile(r"(transChecker_V(?P<version>[\d.]+)\.tar\.gz)"),
}


def _version_key(version: str) -> Tuple[int, ...]:
    return tuple(int(part) for part in version.split("."))


def _iso_utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def fetch_index_html() -> str:
    with urllib.request.urlopen(DDBJ_TOOL_INDEX, timeout=30) as response:
        return response.read().decode("utf-8")


def resolve_latest_archives(html: str) -> Dict[str, Tuple[str, str]]:
    resolved: Dict[str, Tuple[str, str]] = {}
    for component, pattern in _PATTERNS.items():
        matches = [
            (match.group("version"), match.group(1))
            for match in pattern.finditer(html)
        ]
        if not matches:
            continue
        version, filename = max(matches, key=lambda item: _version_key(item[0]))
        resolved[component] = (version, filename)
    return resolved


def cache_root(path: Optional[str | Path] = None) -> Path:
    if path is None:
        return default_cache_dir()
    return Path(path).expanduser().resolve()


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _installation_tree_sha256(root: Path) -> str:
    digest = hashlib.sha256()
    metadata_path = root / ".msspack-install.json"
    for child in sorted(root.rglob("*"), key=lambda path: str(path.relative_to(root))):
        if child == metadata_path or child.is_dir():
            continue
        if child.is_symlink() or not child.is_file():
            raise MSSPackError(f"Unsupported file in DDBJ tool installation: {child}")
        relative = str(child.relative_to(root)).encode("utf-8")
        digest.update(relative)
        digest.update(b"\0")
        digest.update(f"{stat.S_IMODE(child.stat().st_mode):o}".encode("ascii"))
        digest.update(b"\0")
        with child.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        digest.update(b"\0")
    return digest.hexdigest()


def _write_json(path: Path, payload: dict[str, object]) -> Path:
    write_text(path, json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return path


def _read_json(path: Path) -> Optional[dict[str, object]]:
    if not path.exists():
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    return payload if isinstance(payload, dict) else None


def _download_metadata_path(archive_path: Path) -> Path:
    return archive_path.parent / f"{archive_path.name}.json"


def _download(
    url: str,
    destination: Path,
    *,
    expected_sha256: Optional[str] = None,
) -> Path:
    ensure_dir(destination.parent)
    descriptor, temporary_name = mkstemp(
        prefix=f".{destination.name}.",
        suffix=".part",
        dir=destination.parent,
    )
    temporary = Path(temporary_name)
    try:
        with urllib.request.urlopen(url, timeout=120) as response, os.fdopen(
            descriptor, "wb"
        ) as handle:
            shutil.copyfileobj(response, handle)
        actual_sha256 = _sha256_file(temporary)
        if expected_sha256 is not None and actual_sha256 != expected_sha256:
            raise MSSPackError(
                f"Checksum mismatch for {destination.name}: expected {expected_sha256}, "
                f"got {actual_sha256}"
            )
        os.replace(temporary, destination)
    except BaseException:
        try:
            os.close(descriptor)
        except OSError:
            pass
        temporary.unlink(missing_ok=True)
        raise
    _write_json(
        _download_metadata_path(destination),
        {
            "url": url,
            "path": str(destination),
            "size": destination.stat().st_size,
            "sha256": actual_sha256,
            "expected_sha256": expected_sha256,
            "checksum_verified": expected_sha256 is not None,
            "downloaded_at": _iso_utc_now(),
        },
    )
    return destination


def _safe_member_path(destination: Path, member_name: str) -> Path:
    base = destination.resolve()
    target = (destination / member_name).resolve()
    try:
        target.relative_to(base)
    except ValueError as exc:
        raise MSSPackError(f"Archive member escapes destination: {member_name}") from exc
    return target


def _extract_zip_safely(archive: zipfile.ZipFile, destination: Path) -> None:
    for info in archive.infolist():
        _safe_member_path(destination, info.filename)
        mode = info.external_attr >> 16
        if stat.S_ISLNK(mode):
            raise MSSPackError(f"Unsupported symlink in zip archive: {info.filename}")
        archive.extract(info, destination)


def _extract_tar_safely(archive: tarfile.TarFile, destination: Path) -> None:
    for member in archive.getmembers():
        _safe_member_path(destination, member.name)
        if member.issym() or member.islnk() or member.isdev():
            raise MSSPackError(f"Unsupported tar member: {member.name}")
        archive.extract(member, destination)


def _is_valid_installation(root: Path, component: str) -> bool:
    if not root.exists() or not root.is_dir():
        return False
    executable = ToolInstallation(
        component=component,
        version="",
        archive_name="",
        root=root,
    ).executable
    if not (executable.exists() and executable.is_file()):
        return False
    metadata = read_installation_metadata(root)
    if metadata is None:
        return False
    archive_name = metadata.get("archive_name")
    match = (
        _PATTERNS[component].fullmatch(archive_name)
        if isinstance(archive_name, str)
        else None
    )
    if (
        metadata.get("component") != component
        or metadata.get("version") != root.name
        or match is None
        or match.group("version") != root.name
        or archive_name not in TRUSTED_ARCHIVE_SHA256
    ):
        return False
    expected_executable_sha256 = metadata.get("executable_sha256")
    if not isinstance(expected_executable_sha256, str):
        return False
    if expected_executable_sha256 != _sha256_file(executable):
        return False
    expected_tree_sha256 = metadata.get("installation_tree_sha256")
    if not isinstance(expected_tree_sha256, str):
        return False
    try:
        if expected_tree_sha256 != _installation_tree_sha256(root):
            return False
    except (MSSPackError, OSError):
        return False
    download = metadata.get("download")
    if isinstance(download, dict):
        trusted_sha256 = TRUSTED_ARCHIVE_SHA256[archive_name]
        if (
            download.get("checksum_verified") is not True
            or download.get("expected_sha256") != trusted_sha256
            or download.get("sha256") != trusted_sha256
        ):
            return False
    return True


def read_installation_metadata(root: Path | ToolInstallation) -> Optional[dict[str, object]]:
    metadata_path = root.metadata_path if isinstance(root, ToolInstallation) else root / ".msspack-install.json"
    return _read_json(metadata_path)


def _write_installation_metadata(
    *,
    component: str,
    version: str,
    archive_name: str,
    root: Path,
    archive_path: Optional[Path],
) -> dict[str, object]:
    installation = ToolInstallation(
        component=component,
        version=version,
        archive_name=archive_name,
        root=root,
    )
    archive_metadata = None
    if archive_path is not None:
        archive_metadata = _read_json(_download_metadata_path(archive_path))
    payload: dict[str, object] = {
        "component": component,
        "version": version,
        "archive_name": archive_name,
        "root": str(root),
        "executable": str(installation.executable),
        "executable_sha256": _sha256_file(installation.executable),
        "installation_tree_sha256": _installation_tree_sha256(root),
        "installed_at": _iso_utc_now(),
    }
    if archive_metadata is not None:
        payload["download"] = archive_metadata
    _write_json(installation.metadata_path, payload)
    return payload


def describe_installation(installation: ToolInstallation) -> dict[str, object]:
    payload: dict[str, object] = {
        "component": installation.component,
        "version": installation.version,
        "archive_name": installation.archive_name,
        "root": str(installation.root),
        "executable": str(installation.executable),
    }
    metadata = installation.metadata or read_installation_metadata(installation)
    if metadata is not None:
        payload["metadata"] = metadata
    return payload


def _unpack(archive_path: Path, destination: Path) -> Path:
    ensure_dir(destination)
    if archive_path.suffix == ".zip":
        with zipfile.ZipFile(archive_path) as archive:
            _extract_zip_safely(archive, destination)
    else:
        with tarfile.open(archive_path) as archive:
            _extract_tar_safely(archive, destination)
    entries = [entry for entry in destination.iterdir() if entry.is_dir()]
    if len(entries) != 1:
        raise MSSPackError(f"Unexpected archive layout in {archive_path}")
    return entries[0]


def install_component(
    component: str,
    *,
    cache_dir: Optional[str | Path] = None,
    reinstall: bool = False,
) -> ToolInstallation:
    component = component.lower()
    if component not in _PATTERNS:
        raise MSSPackError(f"Unsupported component: {component}")

    resolved = resolve_latest_archives(fetch_index_html())
    if component not in resolved:
        raise MSSPackError(f"Could not resolve a download for {component}")

    version, archive_name = resolved[component]
    expected_sha256 = TRUSTED_ARCHIVE_SHA256.get(archive_name)
    if expected_sha256 is None:
        raise MSSPackError(
            f"DDBJ archive {archive_name} has no trusted checksum in this msspack release. "
            "Upgrade msspack after the new tool release has been reviewed."
        )
    root = cache_root(cache_dir) / "ddbj-tools" / component / version
    if root.exists() and not reinstall and _is_valid_installation(root, component):
        metadata = read_installation_metadata(root)
        assert metadata is not None
        return ToolInstallation(component, version, archive_name, root, metadata)

    if root.exists():
        shutil.rmtree(root)

    archive_path = cache_root(cache_dir) / "downloads" / archive_name
    _download(
        f"{DDBJ_TOOL_INDEX}{archive_name}",
        archive_path,
        expected_sha256=expected_sha256,
    )
    actual_sha256 = _sha256_file(archive_path)
    if actual_sha256 != expected_sha256:
        raise MSSPackError(
            f"Checksum mismatch for cached archive {archive_name}: expected {expected_sha256}, "
            f"got {actual_sha256}"
        )

    with TemporaryDirectory(prefix=f"msspack-{component}-") as tmp_dir:
        extracted_root = _unpack(archive_path, Path(tmp_dir))
        ensure_dir(root.parent)
        shutil.move(str(extracted_root), str(root))

    metadata = _write_installation_metadata(
        component=component,
        version=version,
        archive_name=archive_name,
        root=root,
        archive_path=archive_path,
    )
    return ToolInstallation(component, version, archive_name, root, metadata)


def list_installed(cache_dir: Optional[str | Path] = None) -> Dict[str, ToolInstallation]:
    installations: Dict[str, ToolInstallation] = {}
    base = cache_root(cache_dir) / "ddbj-tools"
    if not base.exists():
        return installations

    for component in _PATTERNS:
        component_dir = base / component
        if not component_dir.exists():
            continue
        versions = []
        for entry in component_dir.iterdir():
            if not entry.is_dir():
                continue
            try:
                _version_key(entry.name)
            except ValueError:
                continue
            versions.append(entry.name)
        if not versions:
            continue
        version = max(versions, key=_version_key)
        root = component_dir / version
        if not _is_valid_installation(root, component):
            continue
        metadata = read_installation_metadata(root)
        archive_name = ""
        if metadata is not None:
            archive_name = str(metadata.get("archive_name", ""))
        installations[component] = ToolInstallation(
            component=component,
            version=version,
            archive_name=archive_name,
            root=root,
            metadata=metadata,
        )
    return installations


def ensure_installed(
    components: Iterable[str],
    *,
    cache_dir: Optional[str | Path] = None,
) -> Dict[str, ToolInstallation]:
    installed = list_installed(cache_dir)
    result: Dict[str, ToolInstallation] = {}
    for component in components:
        normalized = component.lower()
        if normalized in installed:
            result[normalized] = installed[normalized]
        else:
            result[normalized] = install_component(normalized, cache_dir=cache_dir)
    return result


def _run_with_java(
    *,
    command: list[str],
    cwd: Path,
    log_path: Path,
    java_cmd: str,
) -> None:
    env = dict(os.environ)
    resolved = shutil.which(java_cmd) if "/" not in java_cmd else str(Path(java_cmd).expanduser().resolve())
    if resolved and Path(resolved).name == "java":
        env["PATH"] = str(Path(resolved).parent) + os.pathsep + env.get("PATH", "")
        run_command(command, cwd=cwd, log_path=log_path, env=env)
        return

    with TemporaryDirectory(prefix="msspack-java-") as tmp_dir:
        shim = Path(tmp_dir) / "java"
        shim.write_text(
            "#!/bin/sh\n"
            f"exec {shlex.quote(java_cmd)} \"$@\"\n",
            encoding="utf-8",
        )
        shim.chmod(0o755)
        env["PATH"] = str(Path(tmp_dir)) + os.pathsep + env.get("PATH", "")
        run_command(command, cwd=cwd, log_path=log_path, env=env)


def run_parser(
    installation: ToolInstallation,
    *,
    ann_path: Path,
    fasta_path: Path,
    heap: str,
    java_cmd: str,
    log_path: Path,
) -> None:
    _run_with_java(
        command=[
            "bash",
            str(installation.executable),
            "-x",
            str(ann_path),
            "-s",
            str(fasta_path),
            "-M",
            heap,
        ],
        cwd=installation.root,
        log_path=log_path,
        java_cmd=java_cmd,
    )


def run_transchecker(
    installation: ToolInstallation,
    *,
    ann_path: Path,
    fasta_path: Path,
    aa_out: Path,
    nuc_out: Path,
    heap: str,
    java_cmd: str,
    log_path: Path,
) -> None:
    _run_with_java(
        command=[
            "bash",
            str(installation.executable),
            "-x",
            str(ann_path),
            "-s",
            str(fasta_path),
            "-o",
            str(aa_out),
            "-t",
            str(nuc_out),
            "-M",
            heap,
        ],
        cwd=installation.root,
        log_path=log_path,
        java_cmd=java_cmd,
    )
