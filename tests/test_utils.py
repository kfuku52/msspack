import os
import shlex
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from msspack.utils import default_cache_dir, shell_join


class UtilsTests(unittest.TestCase):
    def test_shell_join_round_trips_posix_arguments(self) -> None:
        command = ["tool", "a b", "x'y", "$HOME", "a&b", ""]

        with patch("msspack.utils.platform.system", return_value="Linux"):
            rendered = shell_join(command)

        self.assertEqual(shlex.split(rendered), command)

    def test_default_cache_dir_respects_xdg_cache_home(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir, patch(
            "msspack.utils.platform.system", return_value="Linux"
        ), patch.dict(os.environ, {"XDG_CACHE_HOME": tmp_dir}):
            cache_dir = default_cache_dir()

        self.assertEqual(cache_dir, Path(tmp_dir) / "msspack")

    def test_default_cache_dir_uses_local_app_data_on_windows(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir, patch(
            "msspack.utils.platform.system", return_value="Windows"
        ), patch.dict(os.environ, {"LOCALAPPDATA": tmp_dir}):
            cache_dir = default_cache_dir()

        self.assertEqual(cache_dir, Path(tmp_dir) / "msspack" / "Cache")


if __name__ == "__main__":
    unittest.main()
