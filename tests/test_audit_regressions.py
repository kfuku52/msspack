import copy
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path
from unittest.mock import patch

from msspack.config import load_config
from msspack.config_validation import validate_functional_annotation_config
from msspack.doctor import _input_checks
from msspack.execution import run_if_needed
from msspack.fasta import iter_fasta, reverse_complement
from msspack.functional_annotation import (
    _materialize_database_file,
    _prepare_diamond_database,
    write_translated_protein_fasta,
)
from msspack.gff import child_ids, parse_attributes
from msspack.gff_cleanup import drop_duplicate_coordinate_genes
from msspack.output_state import output_directory_lock, publish_submission
from msspack.padding_tools import write_spliced_cds_fasta
from msspack.pipeline import run_pipeline
from msspack.submission_update import prepare_update
from msspack.utils import MSSPackError
from msspack.workflow import run_all

FIXTURE = Path(__file__).parent / "fixtures" / "minimal_pack"


class AuditRegressionTests(unittest.TestCase):
    def setUp(self) -> None:
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.base = Path(temporary.name)
        shutil.copytree(FIXTURE, self.base, dirs_exist_ok=True)
        self.config = self.base / "config.toml"

    def test_shared_descendants_survive_coordinate_deduplication(self) -> None:
        gff = self.base / "input.gff3"
        gff.write_text(
            "ctg1\t.\tgene\t1\t9\t.\t+\t.\tID=g1\n"
            "ctg1\t.\tgene\t1\t9\t.\t+\t.\tID=g2\n"
            "ctg1\t.\tmRNA\t1\t9\t.\t+\t.\tID=t%2C1;Parent=g1\n"
            "ctg1\t.\tmRNA\t1\t9\t.\t+\t.\tID=t2;Parent=g2\n"
            "ctg1\t.\tCDS\t1\t9\t.\t+\t0\tID=shared;Parent=t%2C1,t2\n"
            "ctg1\t.\texon\t1\t9\t.\t+\t.\tParent=t%2C1,t2\n"
            "ctg1\t.\tmisc_feature\t2\t3\t.\t+\t.\tParent=shared\n"
        )
        for policy in ("first", "longest_valid_cds"):
            with self.subTest(policy=policy):
                target = self.base / "dedup.gff"
                drop_duplicate_coordinate_genes(
                    input_path=gff, fasta_path=self.base / "input.fa", output_path=target,
                    log_path=self.base / "dedup.log", selection_policy=policy,
                )
                result = target.read_text()
                self.assertIn("ID=shared;Parent=t%2C1\n", result)
                self.assertIn("\texon\t1\t9\t.\t+\t.\tParent=t%2C1\n", result)
                self.assertIn("Parent=shared", result)
                self.assertNotIn("ID=g2", result)

    def test_parent_lists_keep_escaped_delimiters(self) -> None:
        attrs = parse_attributes("ID=a%3Bb;Parent=tx%2C1,tx%253B2;Alias=x%3Dy,z")
        self.assertEqual(attrs["ID"], "a;b")
        self.assertEqual(child_ids(attrs["Parent"]), ["tx,1", "tx%3B2"])
        self.assertEqual(child_ids(attrs["Alias"]), ["x=y", "z"])
        self.assertEqual(child_ids(copy.deepcopy(attrs)["Parent"]), ["tx,1", "tx%3B2"])

    def test_rootless_and_gene_direct_models_preserve_translation_and_emit_one_cds(self) -> None:
        for parent_type in ("mRNA", "gene"):
            for strand in ("+", "-"):
                for phase in (0, 1, 2):
                    with self.subTest(parent_type=parent_type, strand=strand, phase=phase):
                        coding = "A" * phase + "ATGAAATAA"
                        sequence = coding if strand == "+" else reverse_complement(coding)
                        (self.base / "input.fa").write_text(f">ctg1\n{sequence}\n")
                        # Deliberately put the child before its parent.
                        (self.base / "input.gff3").write_text(
                            f"ctg1\t.\tCDS\t1\t{len(sequence)}\t.\t{strand}\t{phase}\tParent=t\n"
                            f"ctg1\t.\t{parent_type}\t1\t{len(sequence)}\t.\t{strand}\t.\tID=t\n"
                        )
                        outputs = run_pipeline(self.config, validate=False)
                        cds_rows = [row for row in outputs.ann_path.read_text().splitlines()
                                    if "\tCDS\t" in row]
                        self.assertEqual(len(cds_rows), 1, cds_rows)
                        protein = self.base / "protein.fa"
                        write_translated_protein_fasta(
                            fasta_path=outputs.intermediate / "02.gap-normalized.genome.fasta",
                            gff_path=outputs.intermediate / "12.gff.final-sorted.gff",
                            output_path=protein, genetic_code="1", log_path=self.base / "protein.log",
                            metrics_path=self.base / "protein.json",
                        )
                        self.assertEqual([record.sequence for record in iter_fasta(protein)], ["MK"])

    def test_extractors_handle_child_first_and_spliced_phase_on_both_strands(self) -> None:
        for strand in ("+", "-"):
            with self.subTest(strand=strand):
                coding = "AATG" + "CCCC" + "AAATAA"
                sequence = coding if strand == "+" else reverse_complement(coding)
                cds = [(1, 4, 1), (9, 14, 0)] if strand == "+" else [(1, 6, 0), (11, 14, 1)]
                gff = self.base / "input.gff3"
                gff.write_text("".join(
                    f"ctg1\t.\tCDS\t{start}\t{end}\t.\t{strand}\t{phase}\tParent=t%2C1\n"
                    for start, end, phase in cds
                ) + f"ctg1\t.\tmRNA\t1\t14\t.\t{strand}\t.\tID=t%2C1\n")
                fasta = self.base / "input.fa"
                fasta.write_text(f">ctg1\n{sequence}\n")
                write_spliced_cds_fasta(fasta_path=fasta, gff_path=gff, output_path=self.base / "cds.fa",
                                       log_path=self.base / "extract.log")
                records = list(iter_fasta(self.base / "cds.fa"))
                self.assertEqual([(row.id, row.sequence) for row in records], [("t,1", "ATGAAATAA")])
                write_translated_protein_fasta(
                    fasta_path=fasta, gff_path=gff, output_path=self.base / "protein.fa",
                    genetic_code="1", log_path=self.base / "protein.log", metrics_path=self.base / "metrics.json",
                )
                self.assertEqual([(row.id, row.sequence) for row in iter_fasta(self.base / "protein.fa")],
                                 [("t,1", "MK")])

    def test_failed_rebuild_preserves_published_pair_including_validation_failure(self) -> None:
        outputs = run_pipeline(self.config, validate=False)
        original = (outputs.ann_path.read_bytes(), outputs.fasta_path.read_bytes())
        original_fasta = (self.base / "input.fa").read_bytes()
        (self.base / "input.fa").write_text(">different_contig\nATGAAATAA\n")
        with self.assertRaisesRegex(MSSPackError, "missing from the FASTA"):
            run_pipeline(self.config, validate=False)
        self.assertEqual((outputs.ann_path.read_bytes(), outputs.fasta_path.read_bytes()), original)
        (self.base / "input.fa").write_bytes(original_fasta)
        for error in (MSSPackError("validation failed"), KeyboardInterrupt()):
            with self.subTest(error=type(error).__name__):
                with patch("msspack.pipeline.run_validation", side_effect=error):
                    with self.assertRaises(type(error)):
                        run_pipeline(self.config)
                self.assertEqual((outputs.ann_path.read_bytes(), outputs.fasta_path.read_bytes()), original)

    def test_empty_after_trimming_never_publishes(self) -> None:
        (self.base / "input.fa").write_text(">ctg1\nNNN\n\nNNN\n")
        with self.assertRaisesRegex(MSSPackError, "No sequence remains"):
            run_pipeline(self.config, validate=False)
        self.assertFalse((self.base / "build/Fixture/final/Fix.fasta").exists())

    def test_database_reuse_repairs_corrupt_object_and_materialized_file(self) -> None:
        source = self.base / "reference.fa"
        source.write_text(">ref\nMKK\n")
        kwargs = dict(local_value=str(source), url="", base_dir=self.base,
                      cache_dir=self.base / "db", stem="reference", expected_suffix=".fasta",
                      content_addressed=True)
        first, _, digest, _ = _materialize_database_file(**kwargs)
        for target in (first, self.base / "db/reference.fasta"):
            target.write_text(">bad\nXXX\n")
            restored, _, new_digest, _ = _materialize_database_file(**kwargs)
            self.assertEqual(restored.read_bytes(), source.read_bytes())
            self.assertEqual(digest, new_digest)

    def test_diamond_index_checks_content_and_builder_identity(self) -> None:
        source = self.base / "reference.fa"
        source.write_text(">ref\nMKK\n")
        builder = self.base / "diamond"
        builder.write_text("builder version one")
        builder.chmod(0o755)
        kwargs = dict(name="ref", fasta_path=source, weight=1.0, origin=str(source),
                      sha256=hashlib.sha256(source.read_bytes()).hexdigest(), source_provenance={},
                      cache_dir=self.base / "db", database_root=self.base / "db", command=str(builder))
        def build(command: list[str], **kwargs: object) -> None:
            Path(command[-1] + ".dmnd").write_text("valid index")
        with patch("msspack.functional_annotation.run_command", side_effect=build) as invoked:
            first = _prepare_diamond_database(**kwargs)
            _prepare_diamond_database(**kwargs)
            self.assertEqual(invoked.call_count, 1)
            first.database_path.write_text("broken index")
            _prepare_diamond_database(**kwargs)
            self.assertEqual(invoked.call_count, 2)
            self.assertEqual(first.database_path.read_text(), "valid index")
            builder.write_text("builder version two")
            second = _prepare_diamond_database(**kwargs)
            self.assertEqual(invoked.call_count, 3)
            self.assertNotEqual(first.database_path, second.database_path)
            self.assertEqual(first.database_path.read_text(), "valid index")

    def test_update_rejects_common_entry_macros_without_partial_output(self) -> None:
        ann = self.base / "input.ann"
        ann.write_text("COMMON\tmisc_feature\t1..9\tnote\tsequence @@[entry]@@\n"
                       "ctg1\tsource\t1..9\torganism\tTest organism\n")
        mapping = self.base / "mapping.tsv"
        mapping.write_text("entry\taccession\tsubmitter_seqid\nctg1\tAB123456\tctg1\n")
        with self.assertRaisesRegex(MSSPackError, "COMMON entry macros"):
            prepare_update(ann_path=ann, fasta_path=self.base / "input.fa", mapping_path=mapping,
                           output_dir=self.base / "update")
        self.assertFalse((self.base / "update").exists())

    def test_nonfinite_config_rejected_in_toml_and_library(self) -> None:
        original = self.config.read_text()
        annotation = load_config(self.config).functional_annotation
        for field in ("evalue", "min_bitscore", "swissprot_weight", "pfam_max_i_evalue", "cdd_evalue"):
            for number in ("nan", "inf", "-inf"):
                with self.subTest(field=field, number=number):
                    self.config.write_text(original + f"\n[functional_annotation]\n{field} = {number}\n")
                    with self.assertRaisesRegex(MSSPackError, "finite"):
                        load_config(self.config)
                    with self.assertRaisesRegex(MSSPackError, "finite"):
                        validate_functional_annotation_config(replace(annotation, **{field: float(number)}))

    def test_semicolon_repair_precedes_doctor_and_gap_normalization(self) -> None:
        gff = self.base / "input.gff3"
        gff.write_text(gff.read_text().replace("ID=gene1", "ID=gene1;Note=alpha;beta"))
        original = self.config.read_text()
        for enabled in (False, True):
            with self.subTest(gapjust=enabled):
                self.config.write_text(original.replace("run_gapjust = false", f"run_gapjust = {str(enabled).lower()}"))
                self.assertTrue(all(check.ok for check in _input_checks(load_config(self.config))))
                result = run_all(self.config, run_busco=False, validate=False, write_report=False)
                repaired = result.pipeline.intermediate / "04.gff.semicolons-fixed.gff"
                self.assertIn("Note=alpha.beta", repaired.read_text())

    def test_output_lock_excludes_other_process_and_releases_on_failure(self) -> None:
        code = (
            "from pathlib import Path; from msspack.output_state import output_directory_lock; "
            f"guard=output_directory_lock(Path({str(self.base)!r})); guard.__enter__()"
        )
        with output_directory_lock(self.base), output_directory_lock(self.base):
            result = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=10)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("Another operation", result.stderr)
        result = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=10)
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_publication_pointer_failure_retains_legacy_directory(self) -> None:
        final = self.base / "final"
        final.mkdir()
        (final / "old").write_text("old generation")
        source = self.base / "new"
        source.write_text("new generation")
        if sys.platform == "win32":
            self.skipTest("POSIX pointer failure; Windows uses directory rename")
        with patch("msspack.output_state.os.replace", side_effect=OSError("publish failure")):
            with self.assertRaisesRegex(OSError, "publish failure"):
                publish_submission(self.base, [source])
        self.assertEqual((final / "old").read_text(), "old generation")

    def test_publication_reuses_intact_generation_and_repairs_changed_public_file(self) -> None:
        source = self.base / "input.fa"
        final = publish_submission(self.base, [source])
        generation = final.resolve()
        publish_submission(self.base, [source])
        self.assertEqual(final.resolve(), generation)
        (final / source.name).write_text("changed by consumer")
        publish_submission(self.base, [source])
        self.assertEqual((final / source.name).read_bytes(), source.read_bytes())

    def test_publication_relocates_json_paths_without_changing_other_strings(self) -> None:
        for name in ('plain', 'quoted"path', 'back\\slash', '\u65e5\u672c\u8a9e'):
            with self.subTest(name=name):
                root = self.base / name
                work = root / "work"
                work.mkdir(parents=True)
                source = work / "summary.json"
                payload = {
                    "outputs": {"ann": str(work / "input.ann.txt")},
                    "nested": [str(work / "input.fasta"), None, 3],
                    "sibling": str(root / "work-other" / "result"),
                    "message": f"Read {work} for diagnostics",
                }
                source.write_text(json.dumps(payload), encoding="utf-8")
                final = publish_submission(root, [source])
                actual = json.loads((final / source.name).read_text())
                self.assertEqual(actual["outputs"]["ann"], str(final.resolve() / "input.ann.txt"))
                self.assertEqual(actual["nested"], [str(final.resolve() / "input.fasta"), None, 3])
                self.assertEqual(actual["sibling"], payload["sibling"])
                self.assertEqual(actual["message"], payload["message"])
                self.assertEqual(json.loads(source.read_text()), payload)
                generation = final.resolve()
                publish_submission(root, [source])
                self.assertEqual(final.resolve(), generation)

    def test_publication_rebuilds_generations_with_legacy_json_relocation(self) -> None:
        source = self.base / "input.fa"
        final = publish_submission(self.base, [source])
        generation = final.resolve()
        stamp = final / ".msspack-generation.json"
        payload = json.loads(stamp.read_text())
        payload.pop("schema_version")
        stamp.write_text(json.dumps(payload))
        publish_submission(self.base, [source])
        if sys.platform != "win32":
            self.assertNotEqual(final.resolve(), generation)
        self.assertEqual(json.loads((final / stamp.name).read_text())["schema_version"], 1)

    def test_pack_rejects_duplicate_fasta_without_replacing_published_submission(self) -> None:
        outputs = run_pipeline(self.config, validate=False)
        original = outputs.fasta_path.read_bytes()
        fasta = self.base / "input.fa"
        fasta.write_text(fasta.read_text() + fasta.read_text())
        with self.assertRaisesRegex(MSSPackError, "Duplicate.*ctg1"):
            run_pipeline(self.config, validate=False)
        self.assertEqual(outputs.fasta_path.read_bytes(), original)

    def test_transitive_source_change_invalidates_cache_without_version_bump(self) -> None:
        package = self.base / "package"
        package.mkdir()
        helper = package / "helper.py"
        helper.write_text("value = 1\n")
        output = self.base / "result"
        def action() -> None:
            output.write_text(helper.read_text())
        with patch("msspack.execution.__file__", str(package / "execution.py")):
            self.assertTrue(run_if_needed(outputs=[output], dependencies=[], action=action))
            self.assertFalse(run_if_needed(outputs=[output], dependencies=[], action=action))
            stat = helper.stat()
            helper.write_text("value = 2\n")
            os.utime(helper, ns=(stat.st_atime_ns, stat.st_mtime_ns))
            self.assertTrue(run_if_needed(outputs=[output], dependencies=[], action=action))
        self.assertEqual(output.read_text(), "value = 2\n")
