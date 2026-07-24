import os
import tempfile
import time
import unittest
from concurrent.futures import ThreadPoolExecutor
from math import inf, nan
from pathlib import Path
from unittest.mock import patch

from msspack.database_lock import (
    DatabaseLockSettings,
    _assert_regular_file_or_absent,
    acquire_database_lock,
    database_lock_path,
)


class DatabaseLockTests(unittest.TestCase):
    def test_settings_reject_nonfinite_and_unsafe_stale_intervals(self) -> None:
        for value in (nan, inf):
            with self.subTest(value=value), self.assertRaises(ValueError):
                DatabaseLockSettings(poll_seconds=value).validate()
        with self.assertRaises(ValueError):
            DatabaseLockSettings(
                heartbeat_seconds=60,
                stale_seconds=60,
            ).validate()

    def test_disappearing_path_is_treated_as_absent(self) -> None:
        path = Path("/does/not/matter")
        with patch("msspack.database_lock.os.path.lexists", return_value=True), patch.object(
            Path,
            "lstat",
            side_effect=FileNotFoundError,
        ):
            _assert_regular_file_or_absent(path, label="racing lock")

    def test_repeated_thread_handoffs_do_not_raise(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            lock_path = Path(tmp_dir) / "locks" / "resource.lock"
            settings = DatabaseLockSettings(
                poll_seconds=0.001,
                timeout_seconds=10,
                heartbeat_seconds=0.01,
                stale_seconds=1,
            )
            counter = 0

            def worker() -> None:
                nonlocal counter
                for _ in range(25):
                    with acquire_database_lock(
                        lock_path,
                        label="thread stress database",
                        settings=settings,
                    ):
                        counter += 1

            with patch("builtins.print"), ThreadPoolExecutor(max_workers=8) as executor:
                futures = [executor.submit(worker) for _ in range(8)]
                for future in futures:
                    future.result()

            self.assertEqual(counter, 200)
            self.assertFalse(lock_path.exists())

    def test_acquire_creates_heartbeats_and_removes_lock(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            lock_path = database_lock_path(Path(tmp_dir), "Pfam index")
            settings = DatabaseLockSettings(
                poll_seconds=0.01,
                timeout_seconds=1.0,
                heartbeat_seconds=0.01,
                stale_seconds=1.0,
            )

            with acquire_database_lock(
                lock_path,
                label="test database",
                settings=settings,
            ):
                self.assertTrue(lock_path.is_file())
                first_mtime = lock_path.stat().st_mtime_ns
                time.sleep(0.03)
                self.assertGreaterEqual(lock_path.stat().st_mtime_ns, first_mtime)

            self.assertFalse(os.path.lexists(lock_path))

    def test_waiter_times_out_while_owner_is_alive(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            lock_path = Path(tmp_dir) / "locks" / "resource.lock"
            owner_settings = DatabaseLockSettings(
                poll_seconds=0.01,
                timeout_seconds=1.0,
                heartbeat_seconds=0.01,
                stale_seconds=1.0,
            )
            waiter_settings = DatabaseLockSettings(
                poll_seconds=0.01,
                timeout_seconds=0.03,
                heartbeat_seconds=0.01,
                stale_seconds=1.0,
            )

            with acquire_database_lock(
                lock_path,
                label="owner",
                settings=owner_settings,
            ):
                with self.assertRaises(TimeoutError):
                    with acquire_database_lock(
                        lock_path,
                        label="waiter",
                        settings=waiter_settings,
                    ):
                        self.fail("waiter unexpectedly acquired the live lock")

    def test_recovers_same_host_lock_whose_owner_is_dead(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            lock_path = Path(tmp_dir) / "locks" / "resource.lock"
            lock_path.parent.mkdir()
            lock_path.write_text(
                '{"format":"msspack-database-lock-v1","hostname":"test","pid":123}\n',
                encoding="utf-8",
            )
            settings = DatabaseLockSettings(
                poll_seconds=0.01,
                timeout_seconds=1.0,
                heartbeat_seconds=0.01,
                stale_seconds=60.0,
            )

            with patch("msspack.database_lock._same_local_boot", return_value=True), patch(
                "msspack.database_lock._process_alive",
                return_value=False,
            ):
                with acquire_database_lock(
                    lock_path,
                    label="stale database",
                    settings=settings,
                ):
                    self.assertTrue(lock_path.is_file())

            self.assertFalse(lock_path.exists())

    def test_rejects_symlink_lock_path(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            root = Path(tmp_dir)
            target = root / "target"
            target.write_text("do not replace\n", encoding="utf-8")
            lock_path = root / "resource.lock"
            lock_path.symlink_to(target)

            with self.assertRaises(IsADirectoryError):
                with acquire_database_lock(lock_path, label="unsafe lock"):
                    self.fail("symlink lock was accepted")

            self.assertEqual(target.read_text(encoding="utf-8"), "do not replace\n")

    def test_heartbeat_does_not_recreate_removed_lock(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            lock_path = Path(tmp_dir) / "resource.lock"
            settings = DatabaseLockSettings(
                poll_seconds=0.01,
                timeout_seconds=1.0,
                heartbeat_seconds=0.01,
                stale_seconds=1.0,
            )

            with acquire_database_lock(
                lock_path,
                label="removed lock",
                settings=settings,
            ):
                lock_path.unlink()
                time.sleep(0.03)
                self.assertFalse(lock_path.exists())
