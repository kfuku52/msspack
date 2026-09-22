import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from msspack.config import ToolsConfig
from msspack.utils import MSSPackError
from msspack.validation import (
    ValidationArtifacts,
    ValidationOptions,
    load_validation_summary,
    run_validation,
    validate_existing,
)


class ValidationTests(unittest.TestCase):
    def test_run_validation_writes_structured_pass_results(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            artifacts = ValidationArtifacts(
                ann_path=base / "sample.ann.txt",
                fasta_path=base / "sample.fasta",
                logs_dir=base / "logs",
                output_dir=base / "validation",
            )
            options = ValidationOptions(
                cache_dir=None,
                heap="1G",
                parallel=True,
                run_parser=True,
                run_transchecker=True,
            )
            tools = {
                "parser": SimpleNamespace(version="6.80", root=base / "parser"),
                "transchecker": SimpleNamespace(
                    version="2.26", root=base / "transchecker"
                ),
            }
            artifacts.ann_path.write_text(
                "COMMON\n\tCDS\t1..3\tlocus_tag\tDemo_1\n"
                "\tCDS\t4..6\tlocus_tag\tDemo_2\n",
                encoding="utf-8",
            )
            artifacts.fasta_path.write_text(">ctg1\nATGATG\n", encoding="utf-8")

            def fake_parser(*args: object, **kwargs: object) -> None:
                Path(str(kwargs["log_path"])).write_text(
                    "jParser (Ver. 6.80) finished.\n", encoding="utf-8"
                )

            def fake_transchecker(*args: object, **kwargs: object) -> None:
                Path(str(kwargs["log_path"])).write_text(
                    "TransChecker (Ver. 2.26) finished\n", encoding="utf-8"
                )
                Path(str(kwargs["aa_out"])).write_text(">a\nM\n>b\nM\n", encoding="utf-8")
                Path(str(kwargs["nuc_out"])).write_text(
                    ">a\nATG\n>b\nATG\n", encoding="utf-8"
                )

            for installation in tools.values():
                installation.root.mkdir()

            with patch(
                "msspack.validation.require_installed", return_value=tools
            ), patch("msspack.validation.run_parser", side_effect=fake_parser), patch(
                "msspack.validation.run_transchecker", side_effect=fake_transchecker
            ):
                run_validation(options=options, artifacts=artifacts)

            summary = load_validation_summary(artifacts.validation_summary)
            assert summary is not None
            self.assertEqual(summary.status, "passed")
            self.assertEqual([check.status for check in summary.checks], ["passed", "passed"])
            self.assertEqual(summary.checks[0].version, "6.80")
            self.assertEqual(
                summary.checks[1].record_counts,
                {"annotation_cds": 2, "aa_fasta": 2, "nuc_fasta": 2},
            )

    def test_transchecker_fails_when_outputs_do_not_cover_annotation_cds(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            artifacts = ValidationArtifacts(
                ann_path=base / "sample.ann.txt",
                fasta_path=base / "sample.fasta",
                logs_dir=base / "logs",
                output_dir=base / "validation",
            )
            artifacts.ann_path.write_text(
                "COMMON\n\tCDS\t1..3\tlocus_tag\tDemo_1\n",
                encoding="utf-8",
            )
            artifacts.fasta_path.write_text(">ctg1\nATG\n", encoding="utf-8")
            options = ValidationOptions(
                cache_dir=None,
                heap="1G",
                parallel=False,
                run_parser=False,
                run_transchecker=True,
            )
            tools = {
                "transchecker": SimpleNamespace(
                    version="2.26", root=base / "transchecker"
                )
            }

            def empty_transchecker(*args: object, **kwargs: object) -> None:
                Path(str(kwargs["log_path"])).write_text("finished\n", encoding="utf-8")
                Path(str(kwargs["aa_out"])).write_text("", encoding="utf-8")
                Path(str(kwargs["nuc_out"])).write_text("", encoding="utf-8")

            for installation in tools.values():
                installation.root.mkdir()

            with patch(
                "msspack.validation.require_installed", return_value=tools
            ), patch(
                "msspack.validation.run_transchecker", side_effect=empty_transchecker
            ):
                with self.assertRaisesRegex(
                    MSSPackError,
                    "annotation CDS=1, AA=0, nucleotide=0",
                ):
                    run_validation(options=options, artifacts=artifacts)

            summary = load_validation_summary(artifacts.validation_summary)
            assert summary is not None
            self.assertEqual(summary.status, "failed")
            self.assertEqual(summary.checks[0].status, "not_run")
            self.assertEqual(summary.checks[1].status, "failed")
            self.assertEqual(
                summary.checks[1].record_counts,
                {"annotation_cds": 1, "aa_fasta": 0, "nuc_fasta": 0},
            )

    def test_transchecker_allows_zero_outputs_for_zero_cds_annotation(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            artifacts = ValidationArtifacts(
                ann_path=base / "sample.ann.txt",
                fasta_path=base / "sample.fasta",
                logs_dir=base / "logs",
                output_dir=base / "validation",
            )
            artifacts.ann_path.write_text("COMMON\n", encoding="utf-8")
            artifacts.fasta_path.write_text(">ctg1\nATG\n", encoding="utf-8")
            options = ValidationOptions(
                cache_dir=None,
                heap="1G",
                parallel=False,
                run_parser=False,
                run_transchecker=True,
            )
            tools = {
                "transchecker": SimpleNamespace(
                    version="2.26", root=base / "transchecker"
                )
            }

            def empty_transchecker(*args: object, **kwargs: object) -> None:
                Path(str(kwargs["log_path"])).write_text("finished\n", encoding="utf-8")
                Path(str(kwargs["aa_out"])).write_text("", encoding="utf-8")
                Path(str(kwargs["nuc_out"])).write_text("", encoding="utf-8")

            for installation in tools.values():
                installation.root.mkdir()

            with patch(
                "msspack.validation.require_installed", return_value=tools
            ), patch(
                "msspack.validation.run_transchecker", side_effect=empty_transchecker
            ):
                run_validation(options=options, artifacts=artifacts)

            summary = load_validation_summary(artifacts.validation_summary)
            assert summary is not None
            self.assertEqual(summary.status, "passed")
            self.assertEqual(summary.checks[1].status, "passed")
            self.assertEqual(
                summary.checks[1].record_counts,
                {"annotation_cds": 0, "aa_fasta": 0, "nuc_fasta": 0},
            )

    def test_sequential_parser_failure_marks_transchecker_not_run(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            artifacts = ValidationArtifacts(
                ann_path=base / "sample.ann.txt",
                fasta_path=base / "sample.fasta",
                logs_dir=base / "logs",
                output_dir=base / "validation",
            )
            artifacts.ann_path.write_text("COMMON\n", encoding="utf-8")
            artifacts.fasta_path.write_text(">ctg1\nATG\n", encoding="utf-8")
            options = ValidationOptions(
                cache_dir=None,
                heap="1G",
                parallel=False,
                run_parser=True,
                run_transchecker=True,
            )
            tools = {
                "parser": SimpleNamespace(version="6.80", root=base / "parser"),
                "transchecker": SimpleNamespace(
                    version="2.26", root=base / "transchecker"
                ),
            }

            def fail_parser(*args: object, **kwargs: object) -> None:
                Path(str(kwargs["log_path"])).write_text(
                    "JP0038:WAR:STX:SEQ:short sequence\n",
                    encoding="utf-8",
                )
                raise MSSPackError("parser rejected the submission")

            for installation in tools.values():
                installation.root.mkdir()

            with patch(
                "msspack.validation.require_installed", return_value=tools
            ), patch("msspack.validation.run_parser", side_effect=fail_parser), patch(
                "msspack.validation.run_transchecker"
            ) as mocked_transchecker:
                with self.assertRaisesRegex(MSSPackError, "rejected"):
                    run_validation(options=options, artifacts=artifacts)

            mocked_transchecker.assert_not_called()
            summary = load_validation_summary(artifacts.validation_summary)
            assert summary is not None
            self.assertEqual(summary.status, "failed")
            self.assertEqual(summary.checks[0].status, "failed")
            self.assertEqual(summary.checks[0].warning_count, 1)
            self.assertEqual(summary.checks[1].status, "not_run")
            self.assertIn("earlier validation failure", summary.checks[1].message or "")
            self.assertEqual(summary.checks[1].error_count, 0)

    def test_load_validation_summary_rejects_contradictory_aggregate_fields(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            path = Path(tmp_dir) / "summary.json"
            checks = {
                "parser": {
                    "component": "parser",
                    "label": "Parser",
                    "status": "failed",
                    "version": "6.80",
                    "log_path": None,
                    "outputs": {},
                    "warning_count": 0,
                    "error_count": 1,
                    "record_counts": {},
                    "message": "failed",
                },
                "transchecker": {
                    "component": "transchecker",
                    "label": "transChecker",
                    "status": "not_run",
                    "version": None,
                    "log_path": None,
                    "outputs": {},
                    "warning_count": 0,
                    "error_count": 0,
                    "record_counts": {},
                    "message": "skipped",
                },
            }
            path.write_text(
                json.dumps(
                    {
                        "schema_version": 1,
                        "status": "passed",
                        "attempted": False,
                        "checks": checks,
                    }
                ),
                encoding="utf-8",
            )

            with self.assertRaisesRegex(MSSPackError, "status does not match"):
                load_validation_summary(path)

            path.write_text(
                json.dumps(
                    {
                        "schema_version": 1,
                        "status": "failed",
                        "attempted": False,
                        "checks": checks,
                    }
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(MSSPackError, "attempted flag does not match"):
                load_validation_summary(path)

    def test_run_validation_skips_install_when_all_checks_disabled(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            artifacts = ValidationArtifacts(
                ann_path=base / "sample.ann.txt",
                fasta_path=base / "sample.fasta",
                logs_dir=base / "logs",
                output_dir=base / "validation",
            )
            options = ValidationOptions(
                cache_dir=None,
                heap="1G",
                parallel=False,
                java_cmd="java",
                run_parser=False,
                run_transchecker=False,
            )

            with patch("msspack.validation.require_installed") as mocked:
                result = run_validation(options=options, artifacts=artifacts)

            mocked.assert_not_called()
            self.assertEqual(result, {"validation_summary": artifacts.validation_summary})
            payload = artifacts.validation_summary.read_text(encoding="utf-8")
            self.assertIn('"status": "not_run"', payload)

    def test_validate_existing_forces_both_checks_for_explicit_command(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            config_path = base / "msspack.toml"
            fixture = Path(__file__).parent / "fixtures" / "minimal_pack" / "config.toml"
            config_path.write_text(fixture.read_text(encoding="utf-8"), encoding="utf-8")

            with patch("msspack.validation.run_validation", return_value={}) as mocked:
                validate_existing(
                    config_file=config_path,
                    ann_path=base / "sample.ann.txt",
                    fasta_path=base / "sample.fasta",
                )

            options = mocked.call_args.kwargs["options"]
            self.assertTrue(options.run_parser)
            self.assertTrue(options.run_transchecker)

    def test_validation_options_from_config_includes_java_command(self) -> None:
        config = SimpleNamespace(
            cache_dir=Path("/tmp/cache"),
            tools=ToolsConfig(java="/custom/java", java_heap="2G"),
            pipeline=SimpleNamespace(
                validate_in_parallel=False,
                validate_with_parser=True,
                validate_with_transchecker=True,
            ),
        )

        options = ValidationOptions.from_config(config)

        self.assertEqual(options.java_cmd, "/custom/java")
        self.assertEqual(options.heap, "2G")
