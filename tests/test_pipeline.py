import json
import tempfile
import threading
import unittest
from pathlib import Path

from msspack.build_manifest import ManifestRecorder
from msspack.config import load_config
from msspack.execution import (
    is_up_to_date,
    run_if_needed,
    run_named_jobs,
)
from msspack.utils import MSSPackError


class PipelineCacheTests(unittest.TestCase):
    def test_cached_action_requires_a_declared_output(self) -> None:
        with self.assertRaisesRegex(MSSPackError, "at least one output"):
            run_if_needed(outputs=[], dependencies=[], action=lambda: None)

    def test_manifest_write_preserves_out_of_band_sections(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            config_path = base / "config.toml"
            fixture = Path(__file__).parent / "fixtures" / "minimal_pack" / "config.toml"
            config_path.write_text(
                fixture.read_text(encoding="utf-8").replace("Fixture", "Demo"),
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

    def test_run_if_needed_skips_when_outputs_are_fresh(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            dep = base / "dep.txt"
            out = base / "out.txt"
            dep.write_text("dep", encoding="utf-8")

            calls = {"count": 0}

            def action() -> None:
                calls["count"] += 1
                out.write_text(f"out-{calls['count']}", encoding="utf-8")

            self.assertTrue(
                run_if_needed(outputs=[out], dependencies=[dep], action=action)
            )
            self.assertEqual(calls["count"], 1)

            self.assertFalse(
                run_if_needed(outputs=[out], dependencies=[dep], action=action)
            )
            self.assertEqual(calls["count"], 1)

            dep.write_text("changed", encoding="utf-8")
            self.assertTrue(
                run_if_needed(outputs=[out], dependencies=[dep], action=action)
            )
            self.assertEqual(calls["count"], 2)

    def test_failed_action_never_creates_a_reusable_cache_entry(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            dependency = base / "input.txt"
            output = base / "parser.log"
            dependency.write_text("input", encoding="utf-8")
            calls = 0

            def failing_action() -> None:
                nonlocal calls
                calls += 1
                output.write_text("failed validation", encoding="utf-8")
                raise RuntimeError("failed")

            for _ in range(2):
                with self.assertRaisesRegex(RuntimeError, "failed"):
                    run_if_needed(
                        outputs=[output],
                        dependencies=[dependency],
                        action=failing_action,
                    )

            self.assertEqual(calls, 2)
            self.assertFalse(is_up_to_date([output], [dependency]))

    def test_cache_key_invalidates_otherwise_identical_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            dependency = base / "input.txt"
            output = base / "summary.json"
            dependency.write_text("input", encoding="utf-8")

            run_if_needed(
                outputs=[output],
                dependencies=[dependency],
                cache_key={"lineage": "first"},
                action=lambda: output.write_text("first", encoding="utf-8"),
            )

            self.assertTrue(
                is_up_to_date(
                    [output],
                    [dependency],
                    cache_key={"lineage": "first"},
                )
            )
            self.assertFalse(
                is_up_to_date(
                    [output],
                    [dependency],
                    cache_key={"lineage": "second"},
                )
            )

    def test_missing_dependency_never_reuses_output(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            dependency = base / "input.txt"
            output = base / "output.txt"
            dependency.write_text("input", encoding="utf-8")
            run_if_needed(
                outputs=[output],
                dependencies=[dependency],
                action=lambda: output.write_text("output", encoding="utf-8"),
            )
            dependency.unlink()

            self.assertFalse(is_up_to_date([output], [dependency]))
            with self.assertRaisesRegex(MSSPackError, "missing"):
                run_if_needed(
                    outputs=[output],
                    dependencies=[dependency],
                    action=lambda: None,
                )

    def test_run_named_jobs_parallel_overlaps_work(self) -> None:
        rendezvous = threading.Barrier(2, timeout=5)
        run_named_jobs(
            [("parser", rendezvous.wait), ("transchecker", rendezvous.wait)],
            parallel=True,
        )
