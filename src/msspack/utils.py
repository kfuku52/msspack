from __future__ import annotations

import gzip
import os
import platform
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Mapping, Optional, Sequence


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
    ensure_dir(destination.parent)
    if source.suffix == ".gz":
        with gzip.open(source, "rb") as src, destination.open("wb") as dst:
            shutil.copyfileobj(src, dst)
    else:
        shutil.copy2(source, destination)
    return destination


def link_or_copy(source: Path, destination: Path) -> Path:
    ensure_dir(destination.parent)
    if destination.exists() or destination.is_symlink():
        destination.unlink()
    try:
        os.link(source, destination)
    except OSError:
        shutil.copy2(source, destination)
    return destination


def write_text(path: Path, text: str) -> Path:
    ensure_dir(path.parent)
    path.write_text(text, encoding="utf-8")
    return path


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
