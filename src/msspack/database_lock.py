from __future__ import annotations

import errno
import json
import math
import os
import socket
import stat
import threading
import time
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path

from .utils import ensure_dir


@dataclass(frozen=True)
class DatabaseLockSettings:
    poll_seconds: float = 5.0
    timeout_seconds: float = 86_400.0
    heartbeat_seconds: float = 60.0
    stale_seconds: float = 900.0

    def validate(self) -> None:
        for name, value in (
            ("poll_seconds", self.poll_seconds),
            ("timeout_seconds", self.timeout_seconds),
            ("heartbeat_seconds", self.heartbeat_seconds),
            ("stale_seconds", self.stale_seconds),
        ):
            if not math.isfinite(value) or value <= 0:
                raise ValueError(f"Database lock {name} must be finite and > 0")
        if self.stale_seconds <= self.heartbeat_seconds:
            raise ValueError(
                "Database lock stale_seconds must be greater than heartbeat_seconds"
            )


def database_lock_path(database_root: Path, resource_name: str) -> Path:
    safe_name = "".join(
        character if character.isalnum() or character in "._-" else "_"
        for character in resource_name.strip()
    ).strip("._")
    return database_root / "locks" / f"{safe_name or 'database'}.lock"


def _assert_regular_file_or_absent(path: Path, *, label: str) -> None:
    if not os.path.lexists(path):
        return
    try:
        result = path.lstat()
    except FileNotFoundError:
        return
    if stat.S_ISLNK(result.st_mode) or not stat.S_ISREG(result.st_mode):
        raise IsADirectoryError(f"{label} path exists but is not a regular file: {path}")


def _local_boot_id() -> str | None:
    path = Path("/proc/sys/kernel/random/boot_id")
    if not path.is_file():
        return None
    try:
        value = path.read_text(encoding="utf-8").strip()
    except OSError:
        return None
    return value or None


def _metadata() -> dict[str, object]:
    return {
        "format": "msspack-database-lock-v1",
        "pid": os.getpid(),
        "hostname": socket.gethostname(),
        "boot_id": _local_boot_id(),
        "created_at": time.time(),
    }


def _read_metadata(path: Path) -> dict[str, object] | None:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    return payload if isinstance(payload, dict) else None


def _describe_owner(metadata: dict[str, object] | None) -> str:
    if metadata is None:
        return "owner=unknown"
    parts = []
    if metadata.get("hostname"):
        parts.append(f"host={metadata['hostname']}")
    if metadata.get("pid"):
        parts.append(f"pid={metadata['pid']}")
    created_at = metadata.get("created_at")
    if isinstance(created_at, int | float):
        parts.append(f"created_at={float(created_at):.0f}")
    return ", ".join(parts) if parts else "owner=unknown"


def _try_create(path: Path) -> bool:
    _assert_regular_file_or_absent(path, label="Database lock")
    try:
        descriptor = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o644)
    except FileExistsError:
        return False
    with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
        json.dump(_metadata(), handle, sort_keys=True)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    return True


def _process_alive(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    except OSError as exc:
        return exc.errno != errno.ESRCH
    return True


def _same_local_boot(metadata: dict[str, object] | None) -> bool:
    if metadata is None or metadata.get("hostname") != socket.gethostname():
        return False
    owner_boot = metadata.get("boot_id")
    local_boot = _local_boot_id()
    if owner_boot and local_boot:
        return owner_boot == local_boot
    return True


def _remove_stale_lock(path: Path, *, settings: DatabaseLockSettings) -> bool:
    if not os.path.lexists(path):
        return False
    _assert_regular_file_or_absent(path, label="Database lock")
    try:
        before = path.stat()
    except FileNotFoundError:
        return False
    metadata = _read_metadata(path)
    reason = ""
    raw_pid = metadata.get("pid") if metadata else None
    if _same_local_boot(metadata) and isinstance(raw_pid, int) and not _process_alive(raw_pid):
        reason = f"same-host owner PID {raw_pid} is not running"
    elif time.time() - before.st_mtime > settings.stale_seconds:
        reason = (
            f"heartbeat expired after {time.time() - before.st_mtime:.0f}s "
            f"({_describe_owner(metadata)})"
        )
    if not reason:
        return False
    try:
        current = path.stat()
    except FileNotFoundError:
        return False
    if (
        before.st_ino,
        before.st_size,
        before.st_mtime_ns,
    ) != (
        current.st_ino,
        current.st_size,
        current.st_mtime_ns,
    ):
        return False
    try:
        path.unlink()
    except FileNotFoundError:
        return False
    print(f"Removed stale database lock: {path} ({reason})", flush=True)
    return True


def _start_heartbeat(path: Path, seconds: float) -> tuple[threading.Event, threading.Thread]:
    stop = threading.Event()

    def heartbeat() -> None:
        while not stop.wait(seconds):
            try:
                os.utime(path, None)
            except OSError:
                return

    thread = threading.Thread(
        target=heartbeat,
        name="msspack-database-lock-heartbeat",
        daemon=True,
    )
    thread.start()
    return stop, thread


@contextmanager
def acquire_database_lock(
    lock_path: Path,
    *,
    label: str,
    settings: DatabaseLockSettings | None = None,
) -> Iterator[None]:
    resolved_settings = settings or DatabaseLockSettings()
    resolved_settings.validate()
    ensure_dir(lock_path.parent)
    _assert_regular_file_or_absent(lock_path, label=label)
    started = time.monotonic()
    reported_wait = False
    while True:
        if _try_create(lock_path):
            try:
                acquired_stat = lock_path.stat()
            except FileNotFoundError:
                continue
            stop, thread = _start_heartbeat(
                lock_path,
                resolved_settings.heartbeat_seconds,
            )
            try:
                yield
            finally:
                stop.set()
                thread.join(timeout=max(1.0, resolved_settings.heartbeat_seconds))
                if os.path.lexists(lock_path):
                    _assert_regular_file_or_absent(lock_path, label=label)
                    try:
                        current_stat = lock_path.stat()
                    except FileNotFoundError:
                        pass
                    else:
                        if (
                            acquired_stat.st_dev,
                            acquired_stat.st_ino,
                        ) == (
                            current_stat.st_dev,
                            current_stat.st_ino,
                        ):
                            try:
                                lock_path.unlink()
                            except FileNotFoundError:
                                pass
            return
        if _remove_stale_lock(lock_path, settings=resolved_settings):
            continue
        if not reported_wait:
            print(
                f"Another process holds {label}; waiting: {lock_path} "
                f"({_describe_owner(_read_metadata(lock_path))})",
                flush=True,
            )
            reported_wait = True
        elapsed = time.monotonic() - started
        if elapsed > resolved_settings.timeout_seconds:
            raise TimeoutError(
                f"Timed out after {resolved_settings.timeout_seconds:.0f}s "
                f"waiting for {label}: {lock_path}"
            )
        time.sleep(resolved_settings.poll_seconds)
