import json
import os
import shutil
import tempfile
import unittest
from pathlib import Path

from msspack.pipeline import run_pipeline
from msspack.pipeline_plots import run_pipeline_plots
from msspack.utils import MSSPackError


class IntegrationPackTests(unittest.TestCase):
    def test_curator_options_apply_to_the_published_pair_and_invalidate_cached_inputs(self) -> None:
        fixture = Path(__file__).resolve().parent / "fixtures" / "minimal_pack"
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir) / "input"
            shutil.copytree(fixture, base)
            config_path = base / "config.toml"
            original_config = config_path.read_text()
            config = (
                original_config.replace('gff = "input.gff3"', 'gff = "input.gff3"\nexclude_contigs = ["foreign"]')
                .replace("locus_tag_digits = 6", 'locus_tag_digits = 5\nlocus_tag_number_prefix = "g"')
                .replace("[pipeline]", "[pipeline]\nomit_redundant_mrna = true")
                .replace('genome_coverage = "10x"', 'genome_coverage = "Unknown"')
            )
            fasta_path = base / "input.fa"
            gff_path = base / "input.gff3"
            fasta = fasta_path.read_text() + ">foreign\nATGCCCTAA\n"
            gff = (
                gff_path.read_text().replace("gene1", "g1").replace("tx1", "g1.t1")
                .replace("test protein", "DNA -binding protein")
                + "##sequence-region foreign 1 9\n"
                + "foreign\tsrc\tgene\t1\t9\t.\t+\t.\tID=g2\n"
                + "##FASTA\n>ctg1\nATGAAATAA\n>foreign\nATGCCCTAA\n"
            )
            fasta_path.write_text(fasta)
            gff_path.write_text(gff)
            original = run_pipeline(config_path, validate=False)
            self.assertIn(">foreign", original.fasta_path.read_text())
            self.assertIn("\tmRNA\t", original.ann_path.read_text())

            config_path.write_text(config)
            updated = run_pipeline(config_path, validate=False)
            self.assertEqual(updated.fasta_path.read_text(), (base / "expected.fasta").read_text())
            ann = updated.ann_path.read_text()
            self.assertIn("\tGenome Coverage\tUnknown\n", ann)
            self.assertIn("\tCDS\t1..9\tlocus_tag\tFix_g00001\n", ann)
            self.assertIn("transcript_id:Fix_g00001.t1", ann)
            self.assertIn("\tproduct\tDNA-binding protein\n", ann)
            self.assertNotIn("\tmRNA\t", ann)
            self.assertNotIn("foreign", ann)
            self.assertEqual(fasta_path.read_text(), fasta)
            self.assertEqual(gff_path.read_text(), gff)
            copied_gff = (updated.intermediate / "00.input.annotations.gff").read_text()
            self.assertNotIn("foreign", copied_gff)
            self.assertIn("##FASTA\n>ctg1\nATGAAATAA\n", copied_gff)
            counts = json.loads((updated.logs / "00.copy-input-fasta.metrics.json").read_text())
            self.assertEqual(counts["details"]["excluded_contigs"], {"foreign": 9})
            counts = json.loads((updated.logs / "00.copy-input-gff.metrics.json").read_text())
            self.assertEqual(counts["details"]["excluded_gff_features"], {"foreign": 1})
            review = (updated.final / "submission-review.tsv").read_text()
            self.assertIn("short_cds\tctg1\tCDS\tFix_g00001", review)
            self.assertNotIn("missing_genome_coverage", review)

            config_path.write_text(config.replace('["foreign"]', '["missing"]'))
            with self.assertRaisesRegex(MSSPackError, "absent from FASTA: missing"):
                run_pipeline(config_path, validate=False)
            self.assertEqual(updated.ann_path.read_text(), ann)

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

            database_files = list(
                (
                    base
                    / "msspack_db"
                    / "diamond"
                    / "close-reference"
                ).glob("*/close-reference.dmnd")
            )
            self.assertEqual(len(database_files), 1)
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
