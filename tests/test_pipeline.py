import json
import os
import tempfile
import threading
import time
import unittest
from pathlib import Path

from msspack.build_manifest import ManifestRecorder
from msspack.config import load_config
from msspack.execution import (
    is_up_to_date,
    run_if_needed,
    run_named_jobs,
)
from msspack.validation import ValidationArtifacts, ValidationOptions


class PipelineCacheTests(unittest.TestCase):
    def test_manifest_write_preserves_out_of_band_sections(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            config_path = base / "config.toml"
            (base / "input.fa").write_text(">chr1\nACGT\n", encoding="utf-8")
            (base / "input.gff").write_text("##gff-version 3\n", encoding="utf-8")
            config_path.write_text(
                "\n".join(
                    [
                        "[project]",
                        'name = "Demo"',
                        'output_dir = "build/Demo"',
                        "",
                        "[inputs]",
                        'fasta = "input.fa"',
                        'gff = "input.gff"',
                        "",
                        "[sample]",
                        'locus_tag = "Demo"',
                        "locus_tag_digits = 6",
                        'scientific_name = "Demo demo"',
                        "",
                        "[submission]",
                        'datatype = "WGS"',
                        'hold_date = "20270401"',
                        'bioproject = "PRJDB000001"',
                        'biosample = "SAMD000001"',
                        "",
                        "[submitter]",
                        'ab_name = ["A. Author"]',
                        'contact = "A. Author"',
                        'institute = "Demo Institute"',
                        'department = "Demo Department"',
                        'country = "Japan"',
                        'state = "Tokyo"',
                        'city = "Tokyo"',
                        'street = "1 Demo Street"',
                        'zip = "100-0001"',
                        'phone = "+81-3-0000-0000"',
                        'email = "demo@example.org"',
                        "",
                        "[reference]",
                        'title = "Demo title"',
                        'ab_name = ["A. Author"]',
                        "year = 2026",
                        "",
                        "[st_comment]",
                        'assembly_method = "demo"',
                        'assembly_name = "demo-v1"',
                        'genome_coverage = "10x"',
                        'sequencing_technology = "ONT"',
                        "",
                        "[pipeline]",
                        "validate_with_parser = false",
                        "validate_with_transchecker = false",
                        "",
                        "[busco]",
                        "run_cds = true",
                        "run_genome = false",
                        'lineage_dataset = "embryophyta_odb12"',
                        "auto_lineage = false",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            config = load_config(config_path)
            output_root = base / "build" / "Demo"
            final_dir = output_root / "final"
            final_dir.mkdir(parents=True)
            ann_path = final_dir / "Demo.ann.txt"
            fasta_path = final_dir / "Demo.fasta"
            ann_path.write_text("", encoding="utf-8")
            fasta_path.write_text("", encoding="utf-8")
            manifest_path = output_root / "build-manifest.json"
            manifest_path.write_text(
                '{"plots":{"pipeline":{"gene_flow_pdf":"/tmp/flow.pdf"}},"busco":{"enabled":true}}\n',
                encoding="utf-8",
            )

            recorder = ManifestRecorder(
                config=config,
                config_path=config_path,
                output_root=output_root,
                ann_path=ann_path,
                fasta_path=fasta_path,
            )
            recorder.mark_completed()
            recorder.write(manifest_path)
            payload = json.loads(manifest_path.read_text(encoding="utf-8"))

        self.assertIn("plots", payload)
        self.assertIn("busco", payload)
        self.assertEqual(payload["plots"]["pipeline"]["gene_flow_pdf"], "/tmp/flow.pdf")

    def test_validation_artifacts_from_existing_creates_consistent_paths(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            ann = base / "final" / "sample.ann.txt"
            fasta = base / "final" / "sample.fasta"
            artifacts = ValidationArtifacts.for_existing_submission(
                ann_path=ann,
                fasta_path=fasta,
            )

            self.assertEqual(artifacts.logs_dir, ann.parent / "logs")
            self.assertEqual(artifacts.output_dir, ann.parent / "validation")
            self.assertTrue(artifacts.logs_dir.exists())
            self.assertTrue(artifacts.output_dir.exists())
            self.assertEqual(
                artifacts.result_paths(
                    include_parser=True,
                    include_transchecker=True,
                ),
                {
                    "parser_log": ann.parent / "logs" / "parser.log",
                    "transchecker_log": ann.parent / "logs" / "transchecker.log",
                    "aa_fasta": ann.parent / "validation" / "transChecker.aa.fasta",
                    "nuc_fasta": ann.parent / "validation" / "transChecker.nuc.fasta",
                },
            )

    def test_validation_options_defaults_when_config_is_missing(self) -> None:
        options = ValidationOptions.from_config(None)

        self.assertIsNone(options.cache_dir)
        self.assertEqual(options.heap, "16G")
        self.assertTrue(options.parallel)
        self.assertTrue(options.run_parser)
        self.assertTrue(options.run_transchecker)

    def test_is_up_to_date_uses_oldest_output_and_newest_dependency(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            dep = base / "dep.txt"
            out_a = base / "out_a.txt"
            out_b = base / "out_b.txt"
            dep.write_text("dep", encoding="utf-8")
            out_a.write_text("a", encoding="utf-8")
            out_b.write_text("b", encoding="utf-8")

            os.utime(dep, (100, 100))
            os.utime(out_a, (200, 200))
            os.utime(out_b, (300, 300))
            self.assertTrue(is_up_to_date([out_a, out_b], [dep]))

            os.utime(dep, (400, 400))
            self.assertFalse(is_up_to_date([out_a, out_b], [dep]))

    def test_run_if_needed_skips_when_outputs_are_fresh(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            dep = base / "dep.txt"
            out = base / "out.txt"
            dep.write_text("dep", encoding="utf-8")
            out.write_text("out", encoding="utf-8")
            os.utime(dep, (100, 100))
            os.utime(out, (200, 200))

            calls = {"count": 0}

            def action() -> None:
                calls["count"] += 1

            self.assertFalse(
                run_if_needed(outputs=[out], dependencies=[dep], action=action)
            )
            self.assertEqual(calls["count"], 0)

            os.utime(dep, (300, 300))
            self.assertTrue(
                run_if_needed(outputs=[out], dependencies=[dep], action=action)
            )
            self.assertEqual(calls["count"], 1)

    def test_run_named_jobs_parallel_overlaps_work(self) -> None:
        current = 0
        max_current = 0
        lock = threading.Lock()

        def job() -> None:
            nonlocal current, max_current
            with lock:
                current += 1
                max_current = max(max_current, current)
            time.sleep(0.05)
            with lock:
                current -= 1

        run_named_jobs(
            [("parser", job), ("transchecker", job)],
            parallel=True,
        )
        self.assertEqual(max_current, 2)

    def test_run_named_jobs_sequential_keeps_single_worker(self) -> None:
        current = 0
        max_current = 0
        lock = threading.Lock()

        def job() -> None:
            nonlocal current, max_current
            with lock:
                current += 1
                max_current = max(max_current, current)
            time.sleep(0.02)
            with lock:
                current -= 1

        run_named_jobs(
            [("parser", job), ("transchecker", job)],
            parallel=False,
        )
        self.assertEqual(max_current, 1)
