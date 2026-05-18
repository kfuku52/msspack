import unittest
from pathlib import Path

from msspack.cli import _example_config_text


class TemplateSyncTests(unittest.TestCase):
    def test_cli_template_matches_example_file(self) -> None:
        expected = (
            "# Copy values from examples/msspack.example.toml and edit them for your project.\n"
            + Path("/Users/kf/repos/msspack/examples/msspack.example.toml").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(_example_config_text(), expected)
