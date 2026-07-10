import runpy
import tempfile
import unittest
from pathlib import Path


class ScriptSafetyTests(unittest.TestCase):
    def test_benchmark_refuses_to_remove_config_directory(self) -> None:
        namespace = runpy.run_path(str(Path(__file__).parents[1] / "scripts" / "benchmark_pack.py"))
        remove_output_dir = namespace["remove_output_dir"]
        fixture = Path(__file__).parent / "fixtures" / "minimal_pack"
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            config_path = base / "config.toml"
            config_path.write_text(
                (fixture / "config.toml")
                .read_text(encoding="utf-8")
                .replace('output_dir = "build/Fixture"', 'output_dir = "."'),
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ValueError, "protected"):
                remove_output_dir(config_path)

            self.assertTrue(config_path.exists())


if __name__ == "__main__":
    unittest.main()
