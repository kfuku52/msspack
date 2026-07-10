import hashlib
import io
import json
import tarfile
import tempfile
import unittest
import zipfile
from pathlib import Path
from unittest.mock import patch

from msspack.ddbj_tools import (
    MSSPackError,
    ToolInstallation,
    _download,
    _installation_tree_sha256,
    _sha256_file,
    _unpack,
    describe_installation,
    install_component,
    list_installed,
    read_installation_metadata,
    require_installed,
    resolve_latest_archives,
    run_parser,
)

FAKE_ARCHIVE_SHA256 = hashlib.sha256(b"archive").hexdigest()


class ToolResolutionTests(unittest.TestCase):
    def test_install_component_rejects_native_windows(self) -> None:
        with patch("msspack.ddbj_tools.platform.system", return_value="Windows"):
            with self.assertRaisesRegex(MSSPackError, "Linux and macOS only"):
                install_component("parser")

    def test_resolve_latest_archives(self) -> None:
        html = """
        <a href="Parser_V6.75.tar.gz">Parser_V6.75.tar.gz</a>
        <a href="Parser_V6.80.tar.gz">Parser_V6.80.tar.gz</a>
        <a href="UME_unix_V1.61.zip">UME_unix_V1.61.zip</a>
        <a href="UME_unix_V1.66.zip">UME_unix_V1.66.zip</a>
        <a href="transChecker_V2.20.tar.gz">transChecker_V2.20.tar.gz</a>
        <a href="transChecker_V2.26.tar.gz">transChecker_V2.26.tar.gz</a>
        """
        resolved = resolve_latest_archives(html)
        self.assertEqual(resolved["parser"], ("6.80", "Parser_V6.80.tar.gz"))
        self.assertEqual(resolved["ume"], ("1.66", "UME_unix_V1.66.zip"))
        self.assertEqual(
            resolved["transchecker"],
            ("2.26", "transChecker_V2.26.tar.gz"),
        )

    def test_unpack_rejects_path_traversal_in_zip(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            archive_path = base / "bad.zip"
            with zipfile.ZipFile(archive_path, "w") as archive:
                archive.writestr("../escape.txt", "bad")

            with self.assertRaises(MSSPackError):
                _unpack(archive_path, base / "out")

    def test_unpack_rejects_path_traversal_in_tar(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            archive_path = base / "bad.tar.gz"
            payload = base / "payload.txt"
            payload.write_text("bad", encoding="utf-8")
            with tarfile.open(archive_path, "w:gz") as archive:
                archive.add(payload, arcname="../escape.txt")

            with self.assertRaises(MSSPackError):
                _unpack(archive_path, base / "out")

    def test_list_installed_ignores_invalid_installation_root(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            root = Path(tmp_dir) / "ddbj-tools" / "parser" / "9.99"
            root.mkdir(parents=True)
            self.assertNotIn("parser", list_installed(tmp_dir))

    def test_list_installed_rejects_legacy_install_without_integrity_metadata(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            root = Path(tmp_dir) / "ddbj-tools" / "parser" / "6.80"
            root.mkdir(parents=True)
            (root / "jParser.sh").write_text("#!/bin/sh\n", encoding="utf-8")

            self.assertNotIn("parser", list_installed(tmp_dir))

    def test_require_installed_does_not_download_missing_tools(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir, patch(
            "msspack.ddbj_tools.install_component"
        ) as install:
            with self.assertRaisesRegex(MSSPackError, "msspack tools install parser"):
                require_installed(["parser"], cache_dir=tmp_dir)

        install.assert_not_called()

    def test_list_installed_rejects_tampered_tool_tree(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            root = Path(tmp_dir) / "ddbj-tools" / "parser" / "6.80"
            root.mkdir(parents=True)
            executable = root / "jParser.sh"
            helper = root / "parser.jar"
            executable.write_text("#!/bin/sh\n", encoding="utf-8")
            helper.write_bytes(b"original")
            (root / ".msspack-install.json").write_text(
                json.dumps(
                    {
                        "component": "parser",
                        "version": "6.80",
                        "archive_name": "Parser_V6.80.tar.gz",
                        "executable_sha256": _sha256_file(executable),
                        "installation_tree_sha256": _installation_tree_sha256(root),
                    }
                ),
                encoding="utf-8",
            )
            self.assertIn("parser", list_installed(tmp_dir))

            helper.write_bytes(b"tampered")

            self.assertNotIn("parser", list_installed(tmp_dir))

    def test_install_component_reinstalls_incomplete_root(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            cache_dir = Path(tmp_dir)
            incomplete_root = cache_dir / "ddbj-tools" / "parser" / "9.99"
            incomplete_root.mkdir(parents=True)

            def fake_unpack(_archive_path: Path, destination: Path) -> Path:
                extracted = destination / "Parser"
                extracted.mkdir()
                (extracted / "jParser.sh").write_text("#!/bin/sh\n", encoding="utf-8")
                return extracted

            def fake_download(_url: str, destination: Path, **_kwargs) -> Path:
                destination.parent.mkdir(parents=True, exist_ok=True)
                destination.write_bytes(b"archive")
                return destination

            with patch(
                "msspack.ddbj_tools.fetch_index_html",
                return_value="<html></html>",
            ), patch(
                "msspack.ddbj_tools.resolve_latest_archives",
                return_value={"parser": ("9.99", "Parser_V9.99.tar.gz")},
            ), patch(
                "msspack.ddbj_tools._download",
                side_effect=fake_download,
            ), patch(
                "msspack.ddbj_tools._unpack",
                side_effect=fake_unpack,
            ), patch.dict(
                "msspack.ddbj_tools.TRUSTED_ARCHIVE_SHA256",
                {"Parser_V9.99.tar.gz": FAKE_ARCHIVE_SHA256},
            ):
                installation = install_component("parser", cache_dir=cache_dir)

            self.assertTrue(installation.executable.exists())

    def test_install_component_writes_installation_metadata(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            cache_dir = Path(tmp_dir)

            def fake_unpack(_archive_path: Path, destination: Path) -> Path:
                extracted = destination / "Parser"
                extracted.mkdir()
                (extracted / "jParser.sh").write_text("#!/bin/sh\n", encoding="utf-8")
                return extracted

            def fake_download(_url: str, destination: Path, **_kwargs) -> Path:
                destination.parent.mkdir(parents=True, exist_ok=True)
                destination.write_text("archive", encoding="utf-8")
                return destination

            with patch(
                "msspack.ddbj_tools.fetch_index_html",
                return_value="<html></html>",
            ), patch(
                "msspack.ddbj_tools.resolve_latest_archives",
                return_value={"parser": ("9.99", "Parser_V9.99.tar.gz")},
            ), patch(
                "msspack.ddbj_tools._download",
                side_effect=fake_download,
            ), patch(
                "msspack.ddbj_tools._unpack",
                side_effect=fake_unpack,
            ), patch.dict(
                "msspack.ddbj_tools.TRUSTED_ARCHIVE_SHA256",
                {"Parser_V9.99.tar.gz": FAKE_ARCHIVE_SHA256},
            ):
                installation = install_component("parser", cache_dir=cache_dir)

            metadata = read_installation_metadata(installation)
            self.assertIsNotNone(metadata)
            assert metadata is not None
            self.assertEqual(metadata["component"], "parser")
            self.assertEqual(metadata["version"], "9.99")
            self.assertEqual(metadata["archive_name"], "Parser_V9.99.tar.gz")
            self.assertEqual(
                describe_installation(installation)["metadata"]["version"],
                "9.99",
            )

    def test_download_writes_sidecar_metadata(self) -> None:
        class FakeResponse(io.BytesIO):
            def __enter__(self) -> "FakeResponse":
                return self

            def __exit__(self, exc_type, exc, tb) -> None:
                self.close()

        with tempfile.TemporaryDirectory() as tmp_dir:
            destination = Path(tmp_dir) / "Parser_V9.99.tar.gz"
            with patch(
                "urllib.request.urlopen",
                return_value=FakeResponse(b"archive-bytes"),
            ):
                _download("https://example.test/Parser_V9.99.tar.gz", destination)

            metadata = json.loads(
                Path(f"{destination}.json").read_text(encoding="utf-8")
            )
            self.assertEqual(metadata["url"], "https://example.test/Parser_V9.99.tar.gz")
            self.assertEqual(metadata["size"], len(b"archive-bytes"))
            self.assertIn("sha256", metadata)

    def test_download_rejects_checksum_mismatch_without_publishing_archive(self) -> None:
        class FakeResponse(io.BytesIO):
            def __enter__(self) -> "FakeResponse":
                return self

            def __exit__(self, exc_type, exc, tb) -> None:
                self.close()

        with tempfile.TemporaryDirectory() as tmp_dir:
            destination = Path(tmp_dir) / "Parser_V9.99.tar.gz"
            with patch(
                "urllib.request.urlopen",
                return_value=FakeResponse(b"tampered"),
            ), self.assertRaisesRegex(MSSPackError, "Checksum mismatch"):
                _download(
                    "https://example.test/Parser_V9.99.tar.gz",
                    destination,
                    expected_sha256="0" * 64,
                )

            self.assertFalse(destination.exists())

    def test_install_rejects_unreviewed_archive_version(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir, patch(
            "msspack.ddbj_tools.fetch_index_html",
            return_value="<html></html>",
        ), patch(
            "msspack.ddbj_tools.resolve_latest_archives",
            return_value={"parser": ("9.99", "Parser_V9.99.tar.gz")},
        ):
            with self.assertRaisesRegex(MSSPackError, "trusted checksum"):
                install_component("parser", cache_dir=tmp_dir)

    def test_run_parser_uses_configured_java_binary_directory_when_named_java(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            installation = ToolInstallation(
                component="parser",
                version="1.0",
                archive_name="",
                root=base,
            )
            executable = installation.executable
            executable.write_text("#!/bin/sh\n", encoding="utf-8")
            ann = base / "sample.ann.txt"
            fasta = base / "sample.fasta"
            ann.write_text("", encoding="utf-8")
            fasta.write_text("", encoding="utf-8")

            with patch("msspack.ddbj_tools.run_command") as mocked:
                run_parser(
                    installation,
                    ann_path=ann,
                    fasta_path=fasta,
                    heap="1G",
                    java_cmd="/opt/custom/bin/java",
                    log_path=base / "parser.log",
                )

            env = mocked.call_args.kwargs["env"]
            self.assertTrue(env["PATH"].startswith("/opt/custom/bin"))

    def test_run_parser_creates_java_shim_for_non_java_command_name(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            installation = ToolInstallation(
                component="parser",
                version="1.0",
                archive_name="",
                root=base,
            )
            executable = installation.executable
            executable.write_text("#!/bin/sh\n", encoding="utf-8")
            ann = base / "sample.ann.txt"
            fasta = base / "sample.fasta"
            ann.write_text("", encoding="utf-8")
            fasta.write_text("", encoding="utf-8")

            captured = {}

            def fake_run_command(*_args, **kwargs) -> None:
                env = kwargs["env"]
                shim_dir = env["PATH"].split(":")[0]
                shim_path = Path(shim_dir) / "java"
                captured["shim_path"] = shim_path
                captured["shim_text"] = shim_path.read_text(encoding="utf-8")

            with patch("msspack.ddbj_tools.run_command", side_effect=fake_run_command):
                run_parser(
                    installation,
                    ann_path=ann,
                    fasta_path=fasta,
                    heap="1G",
                    java_cmd="java17",
                    log_path=base / "parser.log",
                )

            self.assertIn('exec java17 "$@"', captured["shim_text"])
