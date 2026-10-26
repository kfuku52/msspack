import json
import tempfile
import threading
import unittest
from pathlib import Path
from unittest.mock import patch

import msspack.execution as execution
from msspack.build_manifest import ManifestRecorder
from msspack.config import load_config
from msspack.execution import (
    is_up_to_date,
    run_if_needed,
    run_named_jobs,
)
from msspack.utils import MSSPackError


class PipelineCacheTests(unittest.TestCase):
    def test_source_change_during_action_rejects_reuse(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            helper, output = base / "helper.py", base / "output"
            helper.write_text("value = 1\n")

            def action() -> None:
                output.write_text(helper.read_text())
                helper.write_text("value = 2\n")

            with patch("msspack.execution.__file__", str(base / "execution.py")):
                with self.assertRaisesRegex(MSSPackError, "Implementation changed"):
                    run_if_needed(outputs=[output], dependencies=[], action=action)
                self.assertFalse(is_up_to_date([output], []))

    def test_source_change_after_validation_does_not_retag_old_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            helper, output = base / "helper.py", base / "output"
            helper.write_text("value = 1\n")
            write_state = execution._write_cache_state

            def delayed_write(*args, **kwargs):
                helper.write_text("value = 2\n")
                return write_state(*args, **kwargs)

            with patch("msspack.execution.__file__", str(base / "execution.py")):
                with patch("msspack.execution._write_cache_state", side_effect=delayed_write):
                    run_if_needed(outputs=[output], dependencies=[],
                                  action=lambda: output.write_text(helper.read_text()))
                self.assertEqual(output.read_text(), "value = 1\n")
                self.assertFalse(is_up_to_date([output], []))
                self.assertTrue(run_if_needed(
                    outputs=[output], dependencies=[],
                    action=lambda: output.write_text(helper.read_text()),
                ))
                self.assertTrue(is_up_to_date([output], []))

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
