import contextlib
import io
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from msspack.busco import BuscoArtifacts
from msspack.cli import main
from msspack.databases import DatabaseStatus
from msspack.demo import write_demo_dataset
from msspack.pipeline import PipelineOutputs
from msspack.pipeline_plots import PipelinePlotArtifacts
from msspack.report import ReportArtifacts
from msspack.utils import MSSPackError
from msspack.workflow import _invalidate_compute_cache, run_all


class WorkflowTests(unittest.TestCase):
    def test_run_all_orders_full_workflow_and_preserves_databases_on_force(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            root = Path(tmp_dir)
            output_root = root / "build" / "Demo"
            cache_dir = output_root / "intermediate" / ".msspack-cache"
            cache_dir.mkdir(parents=True)
            (cache_dir / "cached.json").write_text("{}\n", encoding="utf-8")
            database_root = root / "shared-db"
            database_root.mkdir()
            (database_root / "keep.db").write_text("keep\n", encoding="utf-8")
            manifest_path = output_root / "build-manifest.json"
            manifest_path.parent.mkdir(parents=True, exist_ok=True)
            manifest_path.write_text("{}\n", encoding="utf-8")
            pipeline = PipelineOutputs(
                root=output_root,
                intermediate=output_root / "intermediate",
                logs=output_root / "logs",
                final=output_root / "final",
                ann_path=output_root / "final" / "Demo.ann.txt",
                fasta_path=output_root / "final" / "Demo.fasta",
                manifest_path=manifest_path,
            )
            plots = PipelinePlotArtifacts(
                root=output_root / "plots",
                summary_json=output_root / "plots" / "summary.json",
                summary_tsv=output_root / "plots" / "summary.tsv",
                gene_flow_tsv=output_root / "plots" / "flow.tsv",
                gene_flow_svg=output_root / "plots" / "flow.svg",
                gene_flow_pdf=output_root / "plots" / "flow.pdf",
                event_counts_tsv=output_root / "plots" / "events.tsv",
                event_counts_svg=output_root / "plots" / "events.svg",
                event_counts_pdf=output_root / "plots" / "events.pdf",
                coordinate_duplicates_tsv=output_root / "plots" / "duplicates.tsv",
                coordinate_duplicates_svg=output_root / "plots" / "duplicates.svg",
                coordinate_duplicates_pdf=output_root / "plots" / "duplicates.pdf",
                name_consistency_tsv=output_root / "plots" / "names.tsv",
                name_consistency_svg=output_root / "plots" / "names.svg",
                name_consistency_pdf=output_root / "plots" / "names.pdf",
                source_consistency_tsv=output_root / "plots" / "sources.tsv",
                source_consistency_svg=output_root / "plots" / "sources.svg",
                source_consistency_pdf=output_root / "plots" / "sources.pdf",
            )
            busco = BuscoArtifacts(root=output_root / "busco", cds=None, genome=None)
            report = ReportArtifacts(
                root=output_root / "report",
                index_html=output_root / "report" / "index.html",
            )
            status = DatabaseStatus(root=database_root, mode="shared", resources=())
            config = SimpleNamespace(
                output_dir=output_root,
                base_dir=root,
                database_dir=database_root,
                busco_database_dir=database_root / "busco",
                cache_dir=root / "tools-cache",
                busco=SimpleNamespace(run_cds=True, run_genome=False),
                pipeline=SimpleNamespace(
                    validate_with_parser=True,
                    validate_with_transchecker=True,
                ),
            )
            calls: list[str] = []

            with patch("msspack.workflow.load_config", return_value=config), patch(
                "msspack.workflow._preflight",
            ), patch(
                "msspack.workflow.run_busco_comparison",
                side_effect=lambda *_args, **_kwargs: calls.append("busco") or busco,
            ) as busco_mock, patch(
                "msspack.workflow.run_pipeline",
                side_effect=lambda *_args, **_kwargs: calls.append("pipeline") or pipeline,
            ), patch(
                "msspack.workflow.run_pipeline_plots",
                side_effect=lambda *_args, **_kwargs: calls.append("plots") or plots,
            ), patch(
                "msspack.workflow.collect_database_status",
                return_value=status,
            ), patch(
                "msspack.workflow.write_database_manifest",
            ), patch(
                "msspack.workflow.run_html_report",
                side_effect=lambda *_args, **_kwargs: calls.append("report") or report,
            ):
                artifacts = run_all(
                    root / "config.toml",
                    database_dir=str(database_root),
                    force_compute=True,
                )

            self.assertEqual(calls, ["busco", "pipeline", "plots", "report"])
            busco_mock.assert_called_once_with(
                root / "config.toml",
                force=True,
                prepare_only=True,
            )
            self.assertFalse(cache_dir.exists())
            self.assertEqual((database_root / "keep.db").read_text(encoding="utf-8"), "keep\n")
            self.assertEqual(artifacts.report, report)
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            self.assertEqual(manifest["run"]["status"], "completed")
            self.assertEqual(manifest["run"]["database_mode"], "shared")
            self.assertTrue(manifest["run"]["force_compute"])
            self.assertIsNone(manifest["run"]["error"])

    def test_force_compute_refuses_output_root_containing_protected_paths(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            root = Path(tmp_dir)
            cache_dir = root / "project" / "intermediate" / ".msspack-cache"
            cache_dir.mkdir(parents=True)
            (cache_dir / "cached.json").write_text("{}\n", encoding="utf-8")

            with self.assertRaises(MSSPackError):
                _invalidate_compute_cache(
                    root,
                    protected_paths=(root / "project",),
                )

            self.assertTrue(cache_dir.is_dir())

    def test_force_compute_refuses_system_temporary_root(self) -> None:
        with self.assertRaises(MSSPackError):
            _invalidate_compute_cache(
                Path(tempfile.gettempdir()),
                protected_paths=(),
            )

    def test_failed_run_updates_run_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            root = Path(tmp_dir)
            output_root = root / "build" / "Demo"
            manifest_path = output_root / "build-manifest.json"
            pipeline = SimpleNamespace(manifest_path=manifest_path)
            status = DatabaseStatus(root=root / "database", mode="project", resources=())
            config = SimpleNamespace(
                output_dir=output_root,
                base_dir=root,
                database_dir=root / "database",
                busco_database_dir=root / "database/busco",
                cache_dir=root / "tools-cache",
                busco=SimpleNamespace(run_cds=True, run_genome=False),
                pipeline=SimpleNamespace(
                    validate_with_parser=False,
                    validate_with_transchecker=False,
                ),
            )

            with patch("msspack.workflow.load_config", return_value=config), patch(
                "msspack.workflow._preflight",
            ), patch(
                "msspack.workflow.run_pipeline",
                return_value=pipeline,
            ), patch(
                "msspack.workflow.run_pipeline_plots",
                side_effect=RuntimeError("injected plot failure"),
            ), patch(
                "msspack.workflow.collect_database_status",
                return_value=status,
            ):
                with self.assertRaisesRegex(RuntimeError, "injected plot failure"):
                    run_all(
                        root / "config.toml",
                        run_busco=False,
                        validate=False,
                        write_report=False,
                    )

            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            self.assertEqual(manifest["run"]["status"], "failed")
            self.assertEqual(manifest["run"]["error"], "injected plot failure")
            self.assertIsNotNone(manifest["run"]["finished_at"])
            self.assertFalse(manifest["run"]["busco_enabled"])
            self.assertFalse(manifest["run"]["validation_enabled"])

    def test_run_all_skips_busco_when_config_disables_both_modes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            root = write_demo_dataset(Path(tmp_dir) / "demo")
            config_path = root / "config.toml"
            with config_path.open("a", encoding="utf-8") as handle:
                handle.write(
                    '\n[busco]\nrun_cds = false\nrun_genome = false\n'
                    'command = "msspack-test-missing-busco"\n'
                )
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(main(["doctor", "--config", str(config_path)]), 0)
                self.assertEqual(main(["run", "--config", str(config_path)]), 0)
            output_root = root / "build" / "MSSPackDemo"
            manifest = json.loads((output_root / "build-manifest.json").read_text())
            self.assertEqual(manifest["run"]["status"], "completed")
            self.assertFalse(manifest["run"]["busco_enabled"])
            self.assertTrue((output_root / "final" / "MSSPACK_TEST_GENE.ann.txt").is_file())
            self.assertTrue((output_root / "report" / "index.html").is_file())
            self.assertFalse((output_root / "busco").exists())
            # Explicit BUSCO still requires a mode after resolving CLI overrides.
            stderr = io.StringIO()
            with contextlib.redirect_stderr(stderr), patch(
                "msspack.busco._prepare_busco_lineage_database"
            ) as prepare:
                self.assertEqual(main(["busco", "--config", str(config_path)]), 1)
                prepare.assert_not_called()
                self.assertIn("--cds", stderr.getvalue())
                prepare.side_effect = MSSPackError("reached enabled BUSCO preparation")
                self.assertEqual(main(["busco", "--config", str(config_path), "--cds"]), 1)
                prepare.assert_called_once()
                self.assertTrue(prepare.call_args.kwargs["busco"].run_cds)

    def test_doctor_and_run_require_busco_but_no_busco_overrides_preflight(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            root = write_demo_dataset(Path(tmp_dir) / "demo")
            config_path = root / "config.toml"
            with config_path.open("a", encoding="utf-8") as handle:
                handle.write(
                    '\n[busco]\ncommand = "msspack-test-missing-busco"\n'
                    'download_path = "blocked-busco-cache"\n'
                )
            (root / "blocked-busco-cache").write_text("not a directory")
            for command in ("doctor", "run"):
                with self.subTest(command=command):
                    output = io.StringIO()
                    with contextlib.redirect_stdout(output), contextlib.redirect_stderr(output):
                        self.assertEqual(main([command, "--config", str(config_path)]), 1)
                    self.assertIn("[MISSING] BUSCO:", output.getvalue())
                    self.assertIn("[MISSING] BUSCO database root:", output.getvalue())
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(main(["run", "--config", str(config_path), "--no-busco"]), 0)
            manifest_path = root / "build" / "MSSPackDemo" / "build-manifest.json"
            manifest = json.loads(manifest_path.read_text())
            self.assertEqual(manifest["run"]["status"], "completed")
            self.assertFalse(manifest["run"]["busco_enabled"])
