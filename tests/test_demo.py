import contextlib
import io
import json
import tempfile
import unittest
from pathlib import Path

from msspack.cli import main
from msspack.config import load_config
from msspack.demo import DEMO_FILENAMES, write_demo_dataset
from msspack.pipeline import run_pipeline
from msspack.utils import MSSPackError


class DemoDatasetTests(unittest.TestCase):
    def test_write_demo_dataset_copies_all_files_and_requires_force(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            destination = Path(tmp_dir) / "demo"

            output = write_demo_dataset(destination)

            self.assertEqual(output, destination.resolve())
            self.assertEqual(
                {path.name for path in destination.iterdir()},
                set(DEMO_FILENAMES),
            )
            marker = destination / "config.toml"
            marker.write_text("do not overwrite\n", encoding="utf-8")
            with self.assertRaisesRegex(MSSPackError, "Refusing to overwrite"):
                write_demo_dataset(destination)
            self.assertEqual(marker.read_text(encoding="utf-8"), "do not overwrite\n")

            write_demo_dataset(destination, force=True)
            self.assertIn("[project]", marker.read_text(encoding="utf-8"))

    def test_demo_metadata_is_deliberately_non_submittable_and_pseudonymized(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            destination = write_demo_dataset(Path(tmp_dir) / "demo")
            text = "\n".join(
                (destination / filename).read_text(encoding="utf-8") for filename in DEMO_FILENAMES
            )

            for forbidden in (
                "Pluau",
                "Plumbago",
                "ptg000",
                "PRJDB42760",
                "SAMD01944293",
                "DRR1076161",
            ):
                self.assertNotIn(forbidden, text)
            self.assertIn("Msspackia fictitia", text)
            self.assertIn("NOT_FOR_SUBMISSION", text)

    def test_demo_core_pipeline_plot_and_report_match_expected_events(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            destination = write_demo_dataset(Path(tmp_dir) / "demo")
            config_path = destination / "config.toml"

            outputs = run_pipeline(config_path, validate=False)
            plots_dir = outputs.root / "plots"
            plots_dir.mkdir()
            legacy_paths = [
                plots_dir / f"pipeline-gene-overlap.{ext}" for ext in ("tsv", "svg", "pdf")
            ]
            for path in legacy_paths:
                path.write_text("obsolete plot\n", encoding="utf-8")
            stdout = io.StringIO()
            with contextlib.redirect_stdout(stdout):
                self.assertEqual(main(["plot", "--config", str(config_path), "--force"]), 0)
                self.assertEqual(main(["report", "--config", str(config_path)]), 0)
            manifest = json.loads(outputs.manifest_path.read_text(encoding="utf-8"))
            plots = manifest["plots"]["pipeline"]
            for key in (
                "summary_json",
                "summary_tsv",
                "gene_flow_tsv",
                "gene_flow_svg",
                "gene_flow_pdf",
                "event_counts_tsv",
                "event_counts_svg",
                "event_counts_pdf",
            ):
                self.assertTrue(Path(plots[key]).is_file(), key)
            for extension in ("tsv", "svg", "pdf"):
                self.assertTrue(Path(plots["coordinate_duplicates"][extension]).is_file())
            self.assertEqual(plots["coordinate_duplicates"]["total_removed_genes"], 1)
            self.assertTrue(all(not path.exists() for path in legacy_paths))
            self.assertIn(plots["gene_flow_svg"], stdout.getvalue())
            self.assertIn(manifest["report"]["index_html"], stdout.getvalue())

            expected = json.loads(
                (destination / "expected-summary.json").read_text(encoding="utf-8")
            )["metrics"]
            actual = json.loads(Path(plots["summary_json"]).read_text(encoding="utf-8"))["metrics"]
            for metric, value in expected.items():
                self.assertEqual(actual[metric], value, metric)
            self.assertTrue(outputs.ann_path.is_file())
            self.assertTrue(outputs.fasta_path.is_file())
            self.assertTrue(Path(manifest["report"]["index_html"]).is_file())
            duplicate_svg = Path(plots["coordinate_duplicates"]["svg"]).read_text(encoding="utf-8")
            self.assertIn("MSSPACK_TEST_GENE_0006", duplicate_svg)
            self.assertIn("MSSPACK_TEST_GENE_0007", duplicate_svg)
            self.assertIn(
                "MSSPACK_TEST_GENE_0010",
                (outputs.logs / "16.mss-cds-to-misc.changed-gene-ids.txt").read_text(
                    encoding="utf-8"
                ),
            )
            self.assertIn(
                "\tmisc_feature\t",
                outputs.ann_path.read_text(encoding="utf-8"),
            )

    def test_functional_config_uses_only_the_local_reference(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            destination = write_demo_dataset(Path(tmp_dir) / "demo")

            config = load_config(destination / "config.functional.toml")

            annotation = config.functional_annotation
            self.assertTrue(annotation.enabled)
            self.assertEqual(annotation.reference_proteins, "reference.faa")
            self.assertTrue((destination / annotation.reference_proteins).is_file())
            self.assertFalse(annotation.swissprot_enabled)
            self.assertFalse(annotation.uniref90_enabled)
            self.assertFalse(annotation.pfam_enabled)
            self.assertFalse(annotation.cdd_enabled)
            self.assertFalse(annotation.taxonomy.enabled)


if __name__ == "__main__":
    unittest.main()
