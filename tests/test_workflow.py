import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from msspack.busco import BuscoArtifacts
from msspack.databases import DatabaseStatus
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
            root = Path(tmp_dir)
            config_path = root / "config.toml"
            output_root = root / "build" / "Demo"
            manifest_path = output_root / "build-manifest.json"
            status = DatabaseStatus(root=root / "database", mode="project", resources=())
            config = SimpleNamespace(
                output_dir=output_root,
                base_dir=root,
                database_dir=root / "database",
                busco_database_dir=root / "database/busco",
                cache_dir=root / "tools-cache",
                busco=SimpleNamespace(run_cds=False, run_genome=False),
                pipeline=SimpleNamespace(
                    validate_with_parser=False,
                    validate_with_transchecker=False,
                ),
            )

            with patch("msspack.workflow.load_config", return_value=config), patch(
                "msspack.workflow._preflight",
            ) as preflight, patch(
                "msspack.workflow.run_busco_comparison",
            ) as busco_mock, patch(
                "msspack.workflow.run_pipeline",
                side_effect=RuntimeError("stop after busco decision"),
            ), patch(
                "msspack.workflow.collect_database_status",
                return_value=status,
            ):
                with self.assertRaisesRegex(RuntimeError, "stop after busco decision"):
                    run_all(
                        config_path,
                        validate=False,
                        write_report=False,
                    )

            preflight.assert_called_once_with(
                config_path,
                validate=False,
                run_busco=False,
            )
            busco_mock.assert_not_called()
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            self.assertFalse(manifest["run"]["busco_enabled"])
