import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from msspack.config import ToolsConfig
from msspack.validation import (
    ValidationArtifacts,
    ValidationOptions,
    run_validation,
    validate_existing,
)


class ValidationTests(unittest.TestCase):
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
            self.assertEqual(result, {})

    def test_run_validation_returns_only_enabled_outputs(self) -> None:
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
                run_parser=True,
                run_transchecker=False,
            )

            with patch(
                "msspack.validation.require_installed",
                return_value={
                    "parser": SimpleNamespace(
                        executable=base / "jParser.sh",
                        root=base / "parser",
                    ),
                },
            ), patch(
                "msspack.validation.append_job_if_needed",
                side_effect=lambda jobs, **kwargs: jobs.append((kwargs["name"], lambda: None)),
            ):
                result = run_validation(options=options, artifacts=artifacts)

            self.assertEqual(result, {"parser_log": artifacts.parser_log})

    def test_validate_existing_forces_both_checks_for_explicit_command(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            config_path = base / "msspack.toml"
            config_path.write_text(
                "\n".join(
                    [
                        "[project]",
                        'name = "Demo"',
                        "",
                        "[inputs]",
                        'fasta = "input.fa"',
                        'gff = "input.gff3"',
                        "",
                        "[sample]",
                        'locus_tag = "Demo"',
                        'scientific_name = "Demo species"',
                        "",
                        "[submission]",
                        'hold_date = "20261231"',
                        'bioproject = "PRJDB1"',
                        'biosample = "SAMD1"',
                        "",
                        "[submitter]",
                        'ab_name = ["Doe,J."]',
                        'contact = "Jane Doe"',
                        'institute = "NIG"',
                        'department = "Lab"',
                        'country = "Japan"',
                        'state = "Shizuoka"',
                        'city = "Mishima"',
                        'street = "1111 Yata"',
                        'zip = "411-8540"',
                        'phone = "81-00-0000-0000"',
                        'email = "x@example.org"',
                        "",
                        "[reference]",
                        'title = "Demo sequencing"',
                        'ab_name = ["Doe,J."]',
                        "year = 2026",
                        "",
                        "[pipeline]",
                        "validate_with_parser = false",
                        "validate_with_transchecker = false",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

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
