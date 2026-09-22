"""Exclusive project writers and publication of complete submission generations."""
from __future__ import annotations

import json
import os
import shutil
import sys
import threading
import uuid
from collections.abc import Callable, Iterator
from contextlib import contextmanager
from functools import wraps
from pathlib import Path
from typing import ParamSpec, TypeVar

from .config import load_config
from .execution import _fingerprint
from .utils import MSSPackError, ensure_dir

P = ParamSpec("P")
T = TypeVar("T")
_guard = threading.Lock()
_locks: dict[Path, threading.RLock] = {}
_local = threading.local()


@contextmanager
def output_directory_lock(root: Path) -> Iterator[None]:
    """Fail promptly on another writer; nested operations in one thread are allowed."""
    root = root.resolve()
    with _guard:
        lock = _locks.setdefault(root, threading.RLock())
    if not lock.acquire(blocking=False):
        raise MSSPackError(f"Another operation is writing output directory: {root}")
    held: set[Path] = getattr(_local, "held", set())
    _local.held = held
    try:
        if root in held:
            yield
            return
        ensure_dir(root)
        with (root / ".msspack-output.lock").open("a+b") as handle:
            try:
                if sys.platform == "win32":
                    import msvcrt

                    handle.seek(0)
                    handle.write(b"0")
                    handle.flush()
                    handle.seek(0)
                    msvcrt.locking(handle.fileno(), msvcrt.LK_NBLCK, 1)
                else:
                    import fcntl

                    fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
            except OSError as exc:
                raise MSSPackError(f"Another operation is writing output directory: {root}") from exc
            held.add(root)
            try:
                yield
            finally:
                held.remove(root)
                if sys.platform == "win32":
                    import msvcrt

                    handle.seek(0)
                    msvcrt.locking(handle.fileno(), msvcrt.LK_UNLCK, 1)
                else:
                    import fcntl

                    fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
    finally:
        lock.release()


def locked_output(function: Callable[P, T]) -> Callable[P, T]:
    @wraps(function)
    def wrapped(*args: P.args, **kwargs: P.kwargs) -> T:
        config_file = args[0] if args else kwargs.get("config_file")
        if not isinstance(config_file, str | Path):
            raise TypeError("Expected a config file path")
        with output_directory_lock(load_config(config_file).output_dir):
            return function(*args, **kwargs)
    return wrapped


def _relocate_json_paths(value: object, source: Path, destination: Path) -> object:
    if isinstance(value, dict):
        return {key: _relocate_json_paths(item, source, destination)
                for key, item in value.items()}
    if isinstance(value, list):
        return [_relocate_json_paths(item, source, destination) for item in value]
    if isinstance(value, str):
        if value == str(source):
            return str(destination)
        prefix = str(source) + os.sep
        if value.startswith(prefix):
            return str(destination / value[len(prefix):])
    return value


def publish_submission(root: Path, files: list[Path]) -> Path:
    """Publish a complete directory, keeping older generations available to readers."""
    root = root.resolve()
    generations = ensure_dir(root / ".msspack-generations")
    final = root / "final"
    backup = root / ".msspack-previous-final"
    # Recover the bounded directory-to-pointer migration (or Windows directory swap).
    if backup.exists():
        if not os.path.lexists(final):
            backup.rename(final)
        else:
            backup.rename(generations / ("previous-" + uuid.uuid4().hex))
    source_signatures = {path.name: _fingerprint(path) for path in files}
    stamp_name = ".msspack-generation.json"
    stamp = final / stamp_name
    if stamp.is_file():
        try:
            previous = json.loads(stamp.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            previous = {}
        if (isinstance(previous, dict) and previous.get("schema_version") == 1
                and previous.get("sources") == source_signatures):
            published_signatures = {path.name: _fingerprint(final / path.name) for path in files}
            if previous.get("published") == published_signatures:
                return final
    generation = generations / uuid.uuid4().hex
    generation.mkdir()
    try:
        for source in files:
            target = generation / source.name
            shutil.copy2(source, target)
            if target.suffix == ".json":
                # Validation results should link to the immutable generation they describe.
                payload = json.loads(target.read_text(encoding="utf-8"))
                relocated = _relocate_json_paths(
                    payload, source.parent, final if sys.platform == "win32" else generation,
                )
                target.write_text(
                    json.dumps(relocated, indent=2) + "\n", encoding="utf-8",
                )
        (generation / stamp_name).write_text(json.dumps({
            "schema_version": 1,
            "sources": source_signatures,
            "published": {path.name: _fingerprint(generation / path.name) for path in files},
        }, sort_keys=True) + "\n", encoding="utf-8")
        if sys.platform == "win32":
            # Directory renames preserve a complete pair on platforms without symlinks.
            if final.exists():
                final.rename(backup)
            try:
                generation.rename(final)
            except BaseException:
                if backup.exists():
                    backup.rename(final)
                raise
        else:
            pointer = root / (".msspack-final-" + uuid.uuid4().hex)
            pointer.symlink_to(generation.relative_to(root), target_is_directory=True)
            try:
                if final.exists() and not final.is_symlink():
                    final.rename(backup)
                try:
                    os.replace(pointer, final)
                except BaseException:
                    if backup.exists() and not os.path.lexists(final):
                        backup.rename(final)
                    raise
            finally:
                pointer.unlink(missing_ok=True)
        if backup.exists():
            backup.rename(generations / ("previous-" + uuid.uuid4().hex))
        return final
    except BaseException:
        # If publication completed, keep the generation even when later housekeeping fails.
        if generation.exists() and (not final.exists() or final.resolve() != generation):
            shutil.rmtree(generation)
        raise
