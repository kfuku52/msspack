import json
import os
import shutil
import tempfile
import unittest
from pathlib import Path

from msspack.pipeline import run_pipeline
from msspack.pipeline_plots import run_pipeline_plots


class IntegrationPackTests(unittest.TestCase):
    def test_functional_annotation_updates_ann_with_local_reference(self) -> None:
        fixture_dir = Path(__file__).resolve().parent / "fixtures" / "minimal_pack"
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir) / "minimal_pack"
            shutil.copytree(fixture_dir, base)
            fake_diamond = base / "fake-diamond.py"
            fake_diamond.write_text(
                "#!/usr/bin/env python3\n"
                "import sys\n"
                "from pathlib import Path\n"
                "args = sys.argv[1:]\n"
                "if args[0] == 'makedb':\n"
                "    target = args[args.index('--db') + 1]\n"
                "    Path(target + '.dmnd').write_bytes(b'fake database')\n"
                "elif args[0] == 'blastp':\n"
                "    target = Path(args[args.index('--out') + 1])\n"
                "    if '--max-target-seqs' in args and args[args.index('--max-target-seqs') + 1] == '0':\n"
                "        target.write_text('tx1\\ttx1\\t100\\t100\\t100\\t0\\t100\\t2\\t2\\n')\n"
                "    else:\n"
                "        target.write_text('tx1\\tref1\\t100\\t2\\t2\\t100\\t100\\t1e-30\\t100\\tref1 ATP synthase subunit alpha\\n')\n"
                "else:\n"
                "    raise SystemExit(2)\n",
                encoding="utf-8",
            )
            fake_diamond.chmod(0o755)
            (base / "reference.faa").write_text(
                ">ref1 ATP synthase subunit alpha\nMK\n", encoding="utf-8"
            )
            config_path = base / "config.toml"
            config_text = config_path.read_text(encoding="utf-8")
            config_text = config_text.replace(
                "replace_product_patterns = []",
                'replace_product_patterns = ["test protein"]',
            ).replace(
                'java_heap = "1G"',
                'java_heap = "1G"\ncache_dir = "cache"',
            )
            config_text += (
                "\n[functional_annotation]\n"
                "enabled = true\n"
                f'diamond_command = "{fake_diamond}"\n'
                "swissprot_enabled = false\n"
                'reference_proteins = "reference.faa"\n'
                'reference_name = "close-reference"\n'
                "pfam_enabled = false\n"
                "\n[functional_annotation.taxonomy]\n"
                "enabled = false\n"
                "\n[functional_annotation.consistency]\n"
                "enabled = true\n"
            )
            config_path.write_text(config_text, encoding="utf-8")

            outputs = run_pipeline(config_path, validate=False)

            annotation_text = outputs.ann_path.read_text(encoding="utf-8")
            self.assertIn("\t\t\tproduct\tATP synthase subunit alpha", annotation_text)
            evidence = outputs.final / "functional-annotation.tsv"
            self.assertTrue(evidence.is_file())
            self.assertIn(
                "\tclose-reference\tclose-reference\tref1\thigh\t",
                evidence.read_text(encoding="utf-8"),
            )
            manifest = json.loads(outputs.manifest_path.read_text(encoding="utf-8"))
            self.assertEqual(
                manifest["outputs"]["functional_annotation_evidence"]["path"],
                str(evidence),
            )
            name_standardization = (
                outputs.final / "functional-annotation-name-standardization.tsv"
            )
            self.assertTrue(name_standardization.is_file())
            self.assertEqual(
                manifest["outputs"]["functional_annotation_name_standardization"][
                    "path"
                ],
                str(name_standardization),
            )
            taxonomy = outputs.final / "functional-annotation-taxonomy.json"
            self.assertTrue(taxonomy.is_file())
            self.assertEqual(
                json.loads(taxonomy.read_text(encoding="utf-8"))["status"],
                "disabled",
            )
            self.assertEqual(
                manifest["outputs"]["functional_annotation_taxonomy"]["path"],
                str(taxonomy),
            )
            stage_names = [stage["name"] for stage in manifest["stages"]]
            self.assertIn("14a.functional-annotation-extract-proteins", stage_names)
            self.assertIn("14a2.functional-annotation-resolve-taxonomy", stage_names)
            self.assertIn("14b.functional-annotation-primary-search", stage_names)
            self.assertIn("14c.functional-annotation-uniref90-search", stage_names)
            self.assertIn("14d.functional-annotation-pfam-search", stage_names)
            self.assertIn("14e.functional-annotation-cdd-search", stage_names)
            assignment_stage = (
                "14f.functional-annotation-assign-and-standardize-products"
            )
            self.assertIn(assignment_stage, stage_names)
            self.assertIn("14g.functional-annotation-family-search", stage_names)
            self.assertIn("14h.functional-annotation-consistency-audit", stage_names)
            self.assertLess(
                stage_names.index(assignment_stage),
                stage_names.index("14h.functional-annotation-consistency-audit"),
            )
            consistency = outputs.final / "functional-annotation-consistency.tsv"
            review = outputs.final / "functional-annotation-conflicts.tsv"
            self.assertTrue(consistency.is_file())
            self.assertTrue(review.is_file())
            self.assertIn(
                "no_close_family_peer",
                consistency.read_text(encoding="utf-8"),
            )
            self.assertEqual(
                manifest["outputs"]["functional_annotation_consistency"]["path"],
                str(consistency),
            )

            plot_artifacts = run_pipeline_plots(config_path)
            self.assertTrue(plot_artifacts.name_consistency_svg.is_file())
            self.assertTrue(plot_artifacts.source_consistency_svg.is_file())
            self.assertNotIn(
                "Name consistency",
                plot_artifacts.gene_flow_tsv.read_text(encoding="utf-8"),
            )
            self.assertIn(
                "Name consistency (genes)",
                plot_artifacts.gene_flow_svg.read_text(encoding="utf-8"),
            )
            self.assertIn(
                "each bar totals 100%",
                plot_artifacts.name_consistency_svg.read_text(encoding="utf-8"),
            )

    def test_minimal_fixture_matches_expected_outputs(self) -> None:
        fixture_dir = Path(__file__).resolve().parent / "fixtures" / "minimal_pack"
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir) / "minimal_pack"
            shutil.copytree(fixture_dir, base)

            outputs = run_pipeline(base / "config.toml", validate=False)

            self.assertEqual(
                outputs.fasta_path.read_text(encoding="utf-8"),
                (base / "expected.fasta").read_text(encoding="utf-8"),
            )
            self.assertEqual(
                outputs.ann_path.read_text(encoding="utf-8"),
                (base / "expected.ann.txt").read_text(encoding="utf-8"),
            )
            manifest = json.loads(
                (outputs.root / "build-manifest.json").read_text(encoding="utf-8")
            )
            self.assertEqual(manifest["status"], "completed")
            self.assertIn("duration_seconds", manifest)
            self.assertIn("runtime", manifest)
            self.assertFalse(manifest["validation"]["enabled"])
            self.assertEqual(manifest["stage_summary"]["count"], len(manifest["stages"]))
            self.assertGreaterEqual(manifest["stage_summary"]["ran"], 1)
            self.assertEqual(manifest["stages"][0]["name"], "00.copy-input-fasta")
            self.assertTrue(all("duration_seconds" in stage for stage in manifest["stages"]))
            self.assertTrue(
                any(stage["name"] == "17.render-final-annotation" for stage in manifest["stages"])
            )
            copy_log = outputs.logs / "00.copy-input-fasta.log"
            copy_metrics = outputs.logs / "00.copy-input-fasta.metrics.json"
            render_log = outputs.logs / "17.render-final-annotation.log"
            render_metrics = outputs.logs / "17.render-final-annotation.metrics.json"
            self.assertTrue(copy_log.exists())
            self.assertTrue(copy_metrics.exists())
            self.assertTrue(render_log.exists())
            self.assertTrue(render_metrics.exists())
            self.assertIn("Summary table:", copy_log.read_text(encoding="utf-8"))
            self.assertIn("| Metric | Value |", copy_log.read_text(encoding="utf-8"))
            self.assertIn("Input total:", copy_log.read_text(encoding="utf-8"))
            self.assertIn("Changed total:", copy_log.read_text(encoding="utf-8"))
            self.assertIn("Output total:", copy_log.read_text(encoding="utf-8"))
            self.assertIn("Summary table:", render_log.read_text(encoding="utf-8"))
            self.assertIn("| Metric | Value |", render_log.read_text(encoding="utf-8"))
            self.assertIn("Input total:", render_log.read_text(encoding="utf-8"))
            self.assertIn("Changed total:", render_log.read_text(encoding="utf-8"))
            self.assertIn("Output total:", render_log.read_text(encoding="utf-8"))
            self.assertEqual(
                json.loads(copy_metrics.read_text(encoding="utf-8"))["step"],
                "copy-input-fasta",
            )
            self.assertEqual(
                json.loads(render_metrics.read_text(encoding="utf-8"))["step"],
                "render-final-annotation",
            )

    def test_second_identical_run_reuses_every_pipeline_stage(self) -> None:
        fixture_dir = Path(__file__).resolve().parent / "fixtures" / "minimal_pack"
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir) / "minimal_pack"
            shutil.copytree(fixture_dir, base)

            outputs = run_pipeline(base / "config.toml", validate=False)
            run_pipeline(base / "config.toml", validate=False)

            manifest = json.loads(outputs.manifest_path.read_text(encoding="utf-8"))
            self.assertEqual(manifest["stage_summary"]["ran"], 0)
            self.assertEqual(
                manifest["stage_summary"]["reused"],
                manifest["stage_summary"]["count"],
            )

    def test_same_size_same_mtime_input_edit_invalidates_pipeline_cache(self) -> None:
        fixture_dir = Path(__file__).resolve().parent / "fixtures" / "minimal_pack"
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir) / "minimal_pack"
            shutil.copytree(fixture_dir, base)
            config_path = base / "config.toml"
            input_fasta = base / "input.fa"

            outputs = run_pipeline(config_path, validate=False)
            original_stat = input_fasta.stat()
            original_text = input_fasta.read_text(encoding="utf-8")
            changed_text = original_text.replace("AAA", "AAT", 1)
            self.assertEqual(len(changed_text), len(original_text))
            input_fasta.write_text(changed_text, encoding="utf-8")
            os.utime(
                input_fasta,
                ns=(original_stat.st_atime_ns, original_stat.st_mtime_ns),
            )

            run_pipeline(config_path, validate=False)

            manifest = json.loads(outputs.manifest_path.read_text(encoding="utf-8"))
            stages = {stage["name"]: stage for stage in manifest["stages"]}
            self.assertTrue(stages["00.copy-input-fasta"]["ran"])
            self.assertFalse(stages["00.copy-input-gff"]["ran"])
            self.assertGreater(manifest["stage_summary"]["ran"], 1)


if __name__ == "__main__":
    unittest.main()
