from __future__ import annotations

import gzip
import os
import platform
import shutil
import subprocess
import sys
import tempfile
from contextlib import contextmanager
from pathlib import Path
from typing import BinaryIO, Iterator, Mapping, Optional, Sequence, TextIO


class MSSPackError(RuntimeError):
    """Base exception for msspack errors."""


def default_cache_dir() -> Path:
    if platform.system() == "Darwin":
        return Path.home() / "Library" / "Caches" / "msspack"
    return Path.home() / ".cache" / "msspack"


def expand_path(value: str, base_dir: Optional[Path] = None) -> Path:
    path = Path(value).expanduser()
    if path.is_absolute():
        return path
    if base_dir is not None:
        return (base_dir / path).resolve()
    return path.resolve()


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def shell_join(command: Sequence[str]) -> str:
    return " ".join(subprocess.list2cmdline([part]) for part in command)


def run_command(
    command: Sequence[str],
    *,
    cwd: Optional[Path] = None,
    log_path: Optional[Path] = None,
    env: Optional[Mapping[str, str]] = None,
) -> None:
    command = list(command)
    if command:
        resolved = which(command[0])
        if resolved:
            command[0] = resolved

    if log_path is not None:
        ensure_dir(log_path.parent)
        with log_path.open("w", encoding="utf-8") as handle:
            handle.write(f"$ {shell_join(command)}\n\n")
            handle.flush()
            try:
                result = subprocess.run(
                    command,
                    cwd=str(cwd) if cwd else None,
                    env=dict(env) if env else None,
                    stdout=handle,
                    stderr=subprocess.STDOUT,
                    check=False,
                    text=True,
                )
            except FileNotFoundError as exc:
                raise MSSPackError(f"Command not found: {command[0]}") from exc
        if result.returncode != 0:
            raise MSSPackError(
                f"Command failed with exit code {result.returncode}: {shell_join(command)}"
            )
        return

    try:
        result = subprocess.run(
            command,
            cwd=str(cwd) if cwd else None,
            env=dict(env) if env else None,
            check=False,
            text=True,
        )
    except FileNotFoundError as exc:
        raise MSSPackError(f"Command not found: {command[0]}") from exc
    if result.returncode != 0:
        raise MSSPackError(
            f"Command failed with exit code {result.returncode}: {shell_join(command)}"
        )


def copy_or_decompress(source: Path, destination: Path) -> Path:
    with atomic_binary_writer(destination) as dst:
        if source.suffix == ".gz":
            with gzip.open(source, "rb") as src:
                shutil.copyfileobj(src, dst)
        else:
            with source.open("rb") as src:
                shutil.copyfileobj(src, dst)
    return destination


def link_or_copy(source: Path, destination: Path) -> Path:
    ensure_dir(destination.parent)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{destination.name}.",
        suffix=".tmp",
        dir=destination.parent,
    )
    os.close(descriptor)
    temporary = Path(temporary_name)
    temporary.unlink()
    try:
        try:
            os.link(source, temporary)
        except OSError:
            shutil.copy2(source, temporary)
        os.replace(temporary, destination)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise
    return destination


def write_text(path: Path, text: str) -> Path:
    with atomic_text_writer(path) as handle:
        handle.write(text)
    return path


@contextmanager
def atomic_text_writer(path: Path) -> Iterator[TextIO]:
    """Write text to a sibling temporary file and replace the target on success."""
    ensure_dir(path.parent)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.",
        suffix=".tmp",
        dir=path.parent,
        text=True,
    )
    temporary_path = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            yield handle
            handle.flush()
            os.fsync(handle.fileno())
        if path.exists() and not path.is_symlink():
            temporary_path.chmod(path.stat().st_mode & 0o777)
        else:
            temporary_path.chmod(0o644)
        os.replace(temporary_path, path)
    except BaseException:
        temporary_path.unlink(missing_ok=True)
        raise


@contextmanager
def atomic_binary_writer(path: Path) -> Iterator[BinaryIO]:
    """Write bytes to a sibling temporary file and replace the target on success."""
    ensure_dir(path.parent)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.",
        suffix=".tmp",
        dir=path.parent,
    )
    temporary_path = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            yield handle
            handle.flush()
            os.fsync(handle.fileno())
        if path.exists() and not path.is_symlink():
            temporary_path.chmod(path.stat().st_mode & 0o777)
        else:
            temporary_path.chmod(0o644)
        os.replace(temporary_path, path)
    except BaseException:
        temporary_path.unlink(missing_ok=True)
        raise


def which(command: str) -> Optional[str]:
    resolved = shutil.which(command)
    if resolved:
        return resolved

    if "/" not in command:
        for base in (Path(sys.executable).parent, Path(sys.executable).resolve().parent):
            sibling = base / command
            if sibling.exists() and sibling.is_file():
                return str(sibling)

    return None
