import gzip
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from msspack.config import load_config
from msspack.doctor import doctor_succeeded, run_doctor


class DoctorTests(unittest.TestCase):
    def test_run_doctor_accepts_gzipped_fasta_and_gff(self) -> None:
        fixture = Path(__file__).parent / "fixtures" / "minimal_pack"
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            with gzip.open(base / "input.fa.gz", "wt", encoding="utf-8") as handle:
                handle.write(">ctg1\nATG\n")
            with gzip.open(base / "input.gff3.gz", "wt", encoding="utf-8") as handle:
                handle.write("ctg1\tsrc\tgene\t1\t3\t.\t+\t.\tID=g1\n")
            config_path = base / "config.toml"
            config_path.write_text(
                (fixture / "config.toml")
                .read_text(encoding="utf-8")
                .replace('fasta = "input.fa"', 'fasta = "input.fa.gz"')
                .replace('gff = "input.gff3"', 'gff = "input.gff3.gz"'),
                encoding="utf-8",
            )
            config = load_config(config_path)
            with patch("msspack.doctor.which", return_value="/usr/bin/tool"), patch(
                "msspack.doctor._importable", return_value=True
            ), patch("msspack.doctor.list_installed", return_value={}):
                checks = run_doctor(config)

        by_name = {check.name: check for check in checks}
        self.assertTrue(by_name["FASTA records"].ok)
        self.assertTrue(by_name["GFF/FASTA seqids"].ok)
        self.assertNotIn("input format", by_name)

    def test_functional_annotation_tools_are_required_only_when_enabled(self) -> None:
        fixture = Path(__file__).parent / "fixtures" / "minimal_pack"
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            config_path = base / "config.toml"
            config_path.write_text(
                (fixture / "config.toml").read_text(encoding="utf-8")
                + "\n[functional_annotation]\nenabled = true\n",
                encoding="utf-8",
            )
            config = load_config(config_path)
            with patch(
                "msspack.doctor.which",
                side_effect=lambda command: None
                if command in {"diamond", "hmmscan", "hmmpress"}
                else "/usr/bin/tool",
            ), patch("msspack.doctor._importable", return_value=True), patch(
                "msspack.doctor.list_installed", return_value={}
            ):
                checks = run_doctor(config)

        by_name = {check.name: check for check in checks}
        self.assertTrue(by_name["DIAMOND (functional annotation)"].required)
        self.assertTrue(by_name["HMMER hmmscan (Pfam fallback)"].required)
        self.assertTrue(by_name["HMMER hmmpress (Pfam fallback)"].required)
        self.assertFalse(doctor_succeeded(checks))

    def test_run_doctor_rejects_native_windows_for_validation(self) -> None:
        with patch("msspack.doctor.platform.system", return_value="Windows"), patch(
            "msspack.doctor.which", return_value="C:/tool.exe"
        ), patch("msspack.doctor._importable", return_value=True), patch(
            "msspack.doctor.list_installed", return_value={}
        ):
            checks = run_doctor()

        by_name = {check.name: check for check in checks}
        self.assertFalse(by_name["DDBJ validation platform"].ok)
        self.assertTrue(by_name["DDBJ validation platform"].required)

    def test_run_doctor_marks_ume_optional(self) -> None:
        with patch(
            "msspack.doctor.which",
            side_effect=lambda command: "/usr/bin/" + command if command in {"java", "bash"} else None,
        ), patch(
            "msspack.doctor._importable",
            return_value=True,
        ), patch("msspack.doctor.list_installed", return_value={}):
            checks = run_doctor()

        by_name = {check.name: check for check in checks}
        self.assertFalse(by_name["BUSCO (optional)"].ok)
        self.assertFalse(by_name["BUSCO (optional)"].required)
        self.assertFalse(by_name["DDBJ parser"].ok)
        self.assertFalse(by_name["DDBJ transchecker"].ok)
        self.assertFalse(by_name["DDBJ ume (optional)"].ok)
        self.assertFalse(by_name["DDBJ ume (optional)"].required)
        self.assertFalse(doctor_succeeded(checks))

    def test_run_doctor_reports_mismatched_input_seqids(self) -> None:
        fixture = Path(__file__).parent / "fixtures" / "minimal_pack"
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            config_path = base / "config.toml"
            config_path.write_text(
                (fixture / "config.toml")
                .read_text(encoding="utf-8")
                .replace('fasta = "input.fa"', f'fasta = "{fixture / "input.fa"}"')
                .replace('gff = "input.gff3"', 'gff = "mismatch.gff3"'),
                encoding="utf-8",
            )
            (base / "mismatch.gff3").write_text(
                "other\tsrc\tgene\t1\t3\t.\t+\t.\tID=g1\n",
                encoding="utf-8",
            )
            config = load_config(config_path)
            with patch("msspack.doctor.which", return_value="/usr/bin/tool"), patch(
                "msspack.doctor._importable", return_value=True
            ), patch("msspack.doctor.list_installed", return_value={}):
                checks = run_doctor(config)

            by_name = {check.name: check for check in checks}
            self.assertFalse(by_name["GFF/FASTA seqids"].ok)

    def test_run_doctor_reports_empty_fasta_records_and_malformed_gff(self) -> None:
        fixture = Path(__file__).parent / "fixtures" / "minimal_pack"
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            (base / "empty.fa").write_text(">empty\n", encoding="utf-8")
            (base / "bad.gff3").write_text(
                "empty\tsrc\tgene\tbad\t3\t.\t+\t.\tID=g1\n",
                encoding="utf-8",
            )
            config_path = base / "config.toml"
            config_path.write_text(
                (fixture / "config.toml")
                .read_text(encoding="utf-8")
                .replace('fasta = "input.fa"', 'fasta = "empty.fa"')
                .replace('gff = "input.gff3"', 'gff = "bad.gff3"'),
                encoding="utf-8",
            )
            config = load_config(config_path)
            with patch("msspack.doctor.which", return_value="/usr/bin/tool"), patch(
                "msspack.doctor._importable", return_value=True
            ), patch("msspack.doctor.list_installed", return_value={}):
                checks = run_doctor(config)

            by_name = {check.name: check for check in checks}
            self.assertFalse(by_name["non-empty FASTA records"].ok)
            self.assertFalse(by_name["input format"].ok)

    def test_run_doctor_rejects_database_root_that_is_a_file(self) -> None:
        fixture = Path(__file__).parent / "fixtures" / "minimal_pack"
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            database_file = base / "not-a-directory"
            database_file.write_text("file\n", encoding="utf-8")
            config_path = base / "config.toml"
            config_path.write_text(
                (fixture / "config.toml").read_text(encoding="utf-8")
                + f'\n[databases]\nroot = "{database_file}"\n',
                encoding="utf-8",
            )
            config = load_config(config_path)
            with patch("msspack.doctor.which", return_value="/usr/bin/tool"), patch(
                "msspack.doctor._importable", return_value=True
            ), patch("msspack.doctor.list_installed", return_value={}):
                checks = run_doctor(config)

        by_name = {check.name: check for check in checks}
        self.assertFalse(by_name["database root"].ok)
        self.assertTrue(by_name["database root"].required)
        self.assertIn("not a directory", by_name["database root"].detail)
