import json
import shutil
import tempfile
import unittest
from pathlib import Path

from msspack.pipeline import run_pipeline


class IntegrationPackTests(unittest.TestCase):
    def test_minimal_fixture_matches_expected_outputs(self) -> None:
        fixture_dir = Path(__file__).resolve().parent / "fixtures" / "minimal_pack"
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir) / "minimal_pack"
            shutil.copytree(fixture_dir, base)

            outputs = run_pipeline(base / "config.toml", validate=False)

            self.assertEqual(
                outputs.fasta_path.read_text(encoding="utf-8"),
                (base / "expected.fasta").read_text(encoding="utf-8"),
            )
            self.assertEqual(
                outputs.ann_path.read_text(encoding="utf-8"),
                (base / "expected.ann.txt").read_text(encoding="utf-8"),
            )
            manifest = json.loads((outputs.root / "build-manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(manifest["status"], "completed")
            self.assertIn("duration_seconds", manifest)
            self.assertIn("runtime", manifest)
            self.assertFalse(manifest["validation"]["enabled"])
            self.assertEqual(manifest["stage_summary"]["count"], len(manifest["stages"]))
            self.assertGreaterEqual(manifest["stage_summary"]["ran"], 1)
            self.assertEqual(manifest["stages"][0]["name"], "00.copy-input-fasta")
            self.assertTrue(all("duration_seconds" in stage for stage in manifest["stages"]))
            self.assertTrue(any(stage["name"] == "17.render-final-annotation" for stage in manifest["stages"]))
            copy_log = outputs.logs / "00.copy-input-fasta.log"
            copy_metrics = outputs.logs / "00.copy-input-fasta.metrics.json"
            render_log = outputs.logs / "17.render-final-annotation.log"
            render_metrics = outputs.logs / "17.render-final-annotation.metrics.json"
            self.assertTrue(copy_log.exists())
            self.assertTrue(copy_metrics.exists())
            self.assertTrue(render_log.exists())
            self.assertTrue(render_metrics.exists())
            self.assertIn("Summary table:", copy_log.read_text(encoding="utf-8"))
            self.assertIn("| Metric | Value |", copy_log.read_text(encoding="utf-8"))
            self.assertIn("Input total:", copy_log.read_text(encoding="utf-8"))
            self.assertIn("Changed total:", copy_log.read_text(encoding="utf-8"))
            self.assertIn("Output total:", copy_log.read_text(encoding="utf-8"))
            self.assertIn("Summary table:", render_log.read_text(encoding="utf-8"))
            self.assertIn("| Metric | Value |", render_log.read_text(encoding="utf-8"))
            self.assertIn("Input total:", render_log.read_text(encoding="utf-8"))
            self.assertIn("Changed total:", render_log.read_text(encoding="utf-8"))
            self.assertIn("Output total:", render_log.read_text(encoding="utf-8"))
            self.assertEqual(
                json.loads(copy_metrics.read_text(encoding="utf-8"))["step"],
                "copy-input-fasta",
            )
            self.assertEqual(
                json.loads(render_metrics.read_text(encoding="utf-8"))["step"],
                "render-final-annotation",
            )


if __name__ == "__main__":
    unittest.main()
