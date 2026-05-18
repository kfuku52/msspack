import unittest
from pathlib import Path

from msspack.cli import _example_config_text

REPO_ROOT = Path(__file__).resolve().parents[1]


class TemplateSyncTests(unittest.TestCase):
    def test_cli_template_matches_example_file(self) -> None:
        expected = (
            "# Copy values from examples/msspack.example.toml and edit them for your project.\n"
            + (REPO_ROOT / "examples/msspack.example.toml").read_text(encoding="utf-8")
        )
        self.assertEqual(_example_config_text(), expected)
