import csv
import json
import shutil
import tarfile
import tempfile
import unittest
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from unittest.mock import patch

from msspack.config_models import FunctionalAnnotationConfig
from msspack.database_lock import DatabaseLockSettings
from msspack.functional_annotation import (
    DiamondHit,
    DiamondMetadata,
    _cdd_assignments,
    _combined_diamond_assignments,
    _diamond_assignments,
    _fasta_taxonomy,
    _materialize_database_file,
    _pfam_verification,
    _prepare_cdd_database,
    _rpsblast_database_prefix,
    _submission_safe_product,
    _taxonomy_adjusted_product,
    _uniprot_verification,
    _uniref90_download_url,
    apply_functional_annotations,
    run_pfam_domain_search,
    write_translated_protein_fasta,
)
from msspack.utils import MSSPackError


class FunctionalAnnotationTests(unittest.TestCase):
    @staticmethod
    def _write_cdd_data(directory: Path, marker: str) -> None:
        directory.mkdir(parents=True)
        for name in (
            "cddid.tbl",
            "cdtrack.txt",
            "family_superfamily_links",
            "cddannot.dat",
            "cddannot_generic.dat",
            "bitscore_specific.txt",
        ):
            (directory / name).write_text(f"{name}\t{marker}\n", encoding="utf-8")

    @staticmethod
    def _write_cdd_archive(path: Path, marker: str) -> None:
        source = path.parent / f"source-{marker}"
        source.mkdir()
        for suffix in (".aux", ".freq", ".loo", ".rps"):
            (source / f"Cdd{suffix}").write_text(
                f"{suffix}\t{marker}\n",
                encoding="utf-8",
            )
        with tarfile.open(path, "w:gz") as archive:
            for database_file in sorted(source.iterdir()):
                archive.add(database_file, arcname=f"little_endian/{database_file.name}")

    def test_content_addressed_database_files_do_not_overwrite_prior_sources(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            source = base / "reference.faa"
            source.write_text(">first\nMA\n", encoding="utf-8")
            first_path, _origin, first_digest, _provenance = _materialize_database_file(
                local_value=str(source),
                url="",
                base_dir=base,
                cache_dir=base / "database",
                database_root=base / "database-root",
                stem="reference",
                expected_suffix=".fasta",
                content_addressed=True,
            )
            source.write_text(">second\nMK\n", encoding="utf-8")
            second_path, _origin, second_digest, _provenance = _materialize_database_file(
                local_value=str(source),
                url="",
                base_dir=base,
                cache_dir=base / "database",
                database_root=base / "database-root",
                stem="reference",
                expected_suffix=".fasta",
                content_addressed=True,
            )

            self.assertNotEqual(first_digest, second_digest)
            self.assertNotEqual(first_path, second_path)
            self.assertEqual(first_path.read_text(encoding="utf-8"), ">first\nMA\n")
            self.assertEqual(second_path.read_text(encoding="utf-8"), ">second\nMK\n")

    def test_cdd_data_sources_get_immutable_version_directories(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            database_prefix = base / "local-cdd" / "Cdd"
            database_prefix.parent.mkdir()
            for suffix in (".aux", ".freq", ".loo", ".rps"):
                Path(str(database_prefix) + suffix).write_text(
                    f"{suffix}\n",
                    encoding="utf-8",
                )
            first_data = base / "data-first"
            second_data = base / "data-second"
            self._write_cdd_data(first_data, "first")
            self._write_cdd_data(second_data, "second")
            cache = base / "cache"

            _prefix, first_version, _provenance = _prepare_cdd_database(
                config=FunctionalAnnotationConfig(
                    cdd_database=str(database_prefix),
                    cdd_data_dir=str(first_data),
                ),
                base_dir=base,
                cache_dir=cache,
            )
            first_contents = (first_version / "cddid.tbl").read_text(encoding="utf-8")
            _prefix, second_version, _provenance = _prepare_cdd_database(
                config=FunctionalAnnotationConfig(
                    cdd_database=str(database_prefix),
                    cdd_data_dir=str(second_data),
                ),
                base_dir=base,
                cache_dir=cache,
            )

            self.assertNotEqual(first_version, second_version)
            self.assertEqual(first_contents, "cddid.tbl\tfirst\n")
            self.assertEqual(
                (first_version / "cddid.tbl").read_text(encoding="utf-8"),
                "cddid.tbl\tfirst\n",
            )
            self.assertEqual(
                (second_version / "cddid.tbl").read_text(encoding="utf-8"),
                "cddid.tbl\tsecond\n",
            )

    def test_cdd_download_sources_do_not_replace_an_active_database(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            data = base / "data"
            self._write_cdd_data(data, "common")
            first_archive = base / "first.tar.gz"
            second_archive = base / "second.tar.gz"
            self._write_cdd_archive(first_archive, "first")
            self._write_cdd_archive(second_archive, "second")
            archives = {
                "https://example.test/first.tar.gz": first_archive,
                "https://example.test/second.tar.gz": second_archive,
            }

            def fake_download(url: str, destination: Path) -> None:
                shutil.copyfile(archives[url], destination)

            cache = base / "cache"
            with patch(
                "msspack.functional_annotation._download",
                side_effect=fake_download,
            ):
                first_prefix, _data_dir, _provenance = _prepare_cdd_database(
                    config=FunctionalAnnotationConfig(
                        cdd_url="https://example.test/first.tar.gz",
                        cdd_data_dir=str(data),
                    ),
                    base_dir=base,
                    cache_dir=cache,
                )
                first_contents = Path(str(first_prefix) + ".aux").read_text(
                    encoding="utf-8"
                )
                second_prefix, _data_dir, _provenance = _prepare_cdd_database(
                    config=FunctionalAnnotationConfig(
                        cdd_url="https://example.test/second.tar.gz",
                        cdd_data_dir=str(data),
                    ),
                    base_dir=base,
                    cache_dir=cache,
                )

            self.assertNotEqual(first_prefix, second_prefix)
            self.assertEqual(first_contents, ".aux\tfirst\n")
            self.assertEqual(
                Path(str(first_prefix) + ".aux").read_text(encoding="utf-8"),
                ".aux\tfirst\n",
            )
            self.assertEqual(
                Path(str(second_prefix) + ".aux").read_text(encoding="utf-8"),
                ".aux\tsecond\n",
            )
            self.assertEqual(len(tuple((cache / "cdd" / "sources").glob("*.json"))), 2)

    def test_concurrent_cdd_data_sources_remain_internally_consistent(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            database_prefix = base / "local-cdd" / "Cdd"
            database_prefix.parent.mkdir()
            for suffix in (".aux", ".freq", ".loo", ".rps"):
                Path(str(database_prefix) + suffix).write_text(
                    f"{suffix}\n",
                    encoding="utf-8",
                )
            first_data = base / "data-first"
            second_data = base / "data-second"
            self._write_cdd_data(first_data, "first")
            self._write_cdd_data(second_data, "second")
            settings = DatabaseLockSettings(
                poll_seconds=0.001,
                timeout_seconds=5,
                heartbeat_seconds=0.01,
                stale_seconds=1,
            )

            def prepare(data_dir: Path) -> Path:
                _prefix, version, _provenance = _prepare_cdd_database(
                    config=FunctionalAnnotationConfig(
                        cdd_database=str(database_prefix),
                        cdd_data_dir=str(data_dir),
                    ),
                    base_dir=base,
                    cache_dir=base / "cache",
                    lock_settings=settings,
                )
                return version

            with patch("builtins.print"), ThreadPoolExecutor(max_workers=2) as executor:
                first_future = executor.submit(prepare, first_data)
                second_future = executor.submit(prepare, second_data)
                first_version = first_future.result()
                second_version = second_future.result()

            self.assertNotEqual(first_version, second_version)
            for version, marker in (
                (first_version, "first"),
                (second_version, "second"),
            ):
                for name in (
                    "cddid.tbl",
                    "cdtrack.txt",
                    "family_superfamily_links",
                    "cddannot.dat",
                    "cddannot_generic.dat",
                    "bitscore_specific.txt",
                ):
                    self.assertEqual(
                        (version / name).read_text(encoding="utf-8"),
                        f"{name}\t{marker}\n",
                    )

    def test_rpsblast_database_prefix_avoids_whitespace(self) -> None:
        with (
            tempfile.TemporaryDirectory(prefix="msspack-cdd-source ") as source_dir,
            tempfile.TemporaryDirectory(prefix="msspack-cdd-alias-") as alias_dir,
        ):
            source_prefix = Path(source_dir) / "Cdd"
            resolved = _rpsblast_database_prefix(source_prefix, Path(alias_dir))

            self.assertNotIn(" ", str(resolved))
            self.assertEqual(resolved.resolve(), source_prefix.resolve())

    def test_parses_uniprot_and_uniref_taxonomy_headers(self) -> None:
        self.assertEqual(
            _fasta_taxonomy(
                "sp|P1|TEST Test protein OS=Arabidopsis thaliana OX=3702 GN=TEST"
            ),
            (3702, "Arabidopsis thaliana"),
        )
        self.assertEqual(
            _fasta_taxonomy(
                "UniRef90_P1 Cluster protein n=3 Tax=Viridiplantae TaxID=33090 RepID=P1"
            ),
            (33090, "Viridiplantae"),
        )

    def test_taxonomy_weight_selects_the_closer_near_top_description(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            hits = base / "hits.tsv"
            metadata = base / "metadata.tsv"
            hits.write_text(
                "database\tqseqid\tsseqid\tpident\tqlen\tslen\tqcovhsp\t"
                "scovhsp\tevalue\tbitscore\n"
                "swissprot\tq1\tplant\t70\t300\t300\t100\t100\t1e-50\t300\n"
                "swissprot\tq1\tanimal\t70\t300\t300\t100\t100\t1e-50\t300\n",
                encoding="utf-8",
            )
            metadata.write_text(
                "database\tsubject_id\tdescription\tweight\tsubject_taxon_id\t"
                "subject_organism\ttaxonomy_relation\ttaxonomy_weight\n"
                "swissprot\tplant\tPlant-specific kinase\t1\t3702\t"
                "Arabidopsis thaliana\tsame_order\t1.2\n"
                "swissprot\tanimal\tAnimal-specific phosphatase\t1\t9606\t"
                "Homo sapiens\tcross_kingdom\t0.72\n",
                encoding="utf-8",
            )

            assignments = _diamond_assignments(
                hit_path=hits,
                metadata_path=metadata,
                config=FunctionalAnnotationConfig(min_token_score=0.5),
            )

            self.assertEqual(assignments["q1"].product, "Plant-specific kinase")
            self.assertEqual(assignments["q1"].subject_taxon_id, 3702)
            self.assertEqual(assignments["q1"].taxonomy_relation, "same_order")

    def test_distant_low_identity_product_loses_lineage_specific_number(self) -> None:
        product, adjustment = _taxonomy_adjusted_product(
            "RNA-binding protein 38",
            hit=DiamondHit(
                query_id="q1",
                subject_id="mouse",
                database="swissprot",
                identity=36.9,
                query_length=300,
                subject_length=300,
                query_coverage=80.8,
                subject_coverage=91.6,
                evalue=1e-30,
                bitscore=123.0,
            ),
            metadata=DiamondMetadata(
                description="RNA-binding protein 38",
                source_weight=1.0,
                subject_taxon_id=10090,
                subject_organism="Mus musculus",
                taxonomy_relation="cross_kingdom",
                taxonomy_weight=0.72,
            ),
            config=FunctionalAnnotationConfig(),
        )

        self.assertEqual(product, "RNA-binding protein")
        self.assertIn("numbering/localization removed", adjustment)

    def test_closer_uniref_assignment_replaces_distant_swissprot_assignment(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            primary_hits = base / "primary.tsv"
            primary_metadata = base / "primary-metadata.tsv"
            uniref_hits = base / "uniref.tsv"
            uniref_metadata = base / "uniref-metadata.tsv"
            header = (
                "database\tqseqid\tsseqid\tpident\tqlen\tslen\tqcovhsp\t"
                "scovhsp\tevalue\tbitscore\n"
            )
            metadata_header = (
                "database\tsubject_id\tdescription\tweight\tsubject_taxon_id\t"
                "subject_organism\ttaxonomy_relation\ttaxonomy_weight\n"
            )
            primary_hits.write_text(
                header
                + "swissprot\tq1\tmouse\t45\t300\t300\t90\t90\t1e-40\t250\n",
                encoding="utf-8",
            )
            primary_metadata.write_text(
                metadata_header
                + "swissprot\tmouse\tRNA-binding protein 38\t1\t10090\t"
                "Mus musculus\tcross_kingdom\t0.72\n",
                encoding="utf-8",
            )
            uniref_hits.write_text(
                header
                + "uniref90\tq1\tplant\t55\t300\t300\t85\t85\t1e-35\t230\n",
                encoding="utf-8",
            )
            uniref_metadata.write_text(
                metadata_header
                + "uniref90\tplant\tLOW QUALITY PROTEIN: Plant RNA-binding protein\t"
                "0.85\t3702\t"
                "Arabidopsis thaliana\tsame_order\t1.2\n",
                encoding="utf-8",
            )

            assignments = _combined_diamond_assignments(
                inputs=[
                    (primary_hits, primary_metadata),
                    (uniref_hits, uniref_metadata),
                ],
                config=FunctionalAnnotationConfig(min_token_score=0.5),
            )

            self.assertEqual(assignments["q1"].source, "uniref90")
            self.assertEqual(assignments["q1"].product, "Plant RNA-binding protein")
            self.assertEqual(assignments["q1"].taxonomy_relation, "same_order")

    def test_generated_products_are_safe_for_ddbj_qualifiers(self) -> None:
        self.assertEqual(
            _submission_safe_product('regulatory subunit B" gamma\\delta\tprotein'),
            "regulatory subunit B gamma/delta protein",
        )

    def test_rejects_download_with_a_checksum_mismatch(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)

            def fake_download(_url: str, destination: Path) -> None:
                destination.write_bytes(b"not the expected database")

            with (
                patch(
                    "msspack.functional_annotation._remote_verification",
                    return_value=("md5", "00000000000000000000000000000000", None, "test"),
                ),
                patch(
                    "msspack.functional_annotation._download",
                    side_effect=fake_download,
                ),
            ):
                with self.assertRaisesRegex(MSSPackError, "checksum mismatch"):
                    _materialize_database_file(
                        local_value="",
                        url="https://example.org/database.fasta.gz",
                        base_dir=base,
                        cache_dir=base / "cache",
                        stem="database",
                        expected_suffix=".fasta",
                    )

            self.assertFalse((base / "cache" / "database.fasta").exists())

    def test_parses_official_database_checksum_metadata(self) -> None:
        metalink = (
            '<metalink xmlns="http://www.metalinker.org/">'
            "<version>2026_02</version><files>"
            '<file name="uniprot_sprot.fasta.gz"><size>123</size><verification>'
            '<hash type="md5">ABCDEF</hash></verification></file>'
            "</files></metalink>"
        )
        with patch("msspack.functional_annotation._fetch_text", return_value=metalink):
            self.assertEqual(
                _uniprot_verification(
                    "https://ftp.uniprot.org/pub/databases/uniprot/current_release/"
                    "knowledgebase/complete/uniprot_sprot.fasta.gz"
                ),
                ("md5", "abcdef", 123, "2026_02"),
            )
        with patch(
            "msspack.functional_annotation._fetch_text",
            return_value="0123456789abcdef  Pfam-A.hmm.gz\n",
        ):
            self.assertEqual(
                _pfam_verification(
                    "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"
                ),
                ("md5", "0123456789abcdef", None, "current_release"),
            )

    def test_builds_taxon_scoped_uniref90_download_url(self) -> None:
        config = FunctionalAnnotationConfig(uniref90_taxon_id=3398)
        url = _uniref90_download_url(config)
        self.assertIn("rest.uniprot.org/uniref/stream", url)
        self.assertIn("identity%3A0.9", url)
        self.assertIn("taxonomy_id%3A3398", url)

    def test_translates_spliced_cds_on_both_strands(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            fasta = base / "genome.fa"
            gff = base / "genes.gff3"
            proteins = base / "proteins.faa"
            fasta.write_text(">ctg1\nATGAAATAATTATTTCAT\n", encoding="utf-8")
            gff.write_text(
                "ctg1\tsrc\tgene\t1\t9\t.\t+\t.\tID=g1\n"
                "ctg1\tsrc\tmRNA\t1\t9\t.\t+\t.\tID=t1;Parent=g1\n"
                "ctg1\tsrc\tCDS\t1\t9\t.\t+\t0\tID=c1;Parent=t1\n"
                "ctg1\tsrc\tgene\t10\t18\t.\t-\t.\tID=g2\n"
                "ctg1\tsrc\tmRNA\t10\t18\t.\t-\t.\tID=t2;Parent=g2\n"
                "ctg1\tsrc\tCDS\t10\t18\t.\t-\t0\tID=c2;Parent=t2\n",
                encoding="utf-8",
            )

            write_translated_protein_fasta(
                fasta_path=fasta,
                gff_path=gff,
                output_path=proteins,
                genetic_code="1",
                log_path=base / "extract.log",
                metrics_path=base / "extract.json",
            )

            self.assertEqual(proteins.read_text(encoding="utf-8"), ">t1\nMK\n>t2\nMK\n")

    def test_ahrd_consensus_preserves_existing_and_uses_pfam_fallback(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            annotation = base / "annotation.tsv"
            diamond_hits = base / "diamond.tsv"
            diamond_metadata = base / "diamond-metadata.tsv"
            pfam_hits = base / "pfam.domtblout"
            pfam_metadata = base / "pfam-metadata.tsv"
            output = base / "annotated.tsv"
            evidence = base / "evidence.tsv"
            name_summary = base / "name-standardization.tsv"
            annotation.write_text(
                "ID\tDescription\tLocus_tag\n"
                "t1\thypothetical protein\tX_1\n"
                "t2\texisting enzyme\tX_2\n"
                "t3\thypothetical protein\tX_3\n"
                "t4\thypothetical protein\tX_4\n",
                encoding="utf-8",
            )
            diamond_hits.write_text(
                "database\tqseqid\tsseqid\tpident\tqlen\tslen\tqcovhsp\tscovhsp\tevalue\tbitscore\n"
                "swissprot\tt1\tsp1\t80\t300\t300\t100\t100\t1e-80\t500\n"
                "reference\tt1\tref1\t78\t300\t300\t99\t99\t1e-75\t490\n"
                "swissprot\tt2\tsp2\t90\t300\t300\t100\t100\t1e-90\t550\n"
                "swissprot\tt4\tsp3\t70\t300\t300\t100\t100\t1e-40\t300\n"
                "swissprot\tt4\tsp4\t70\t300\t300\t100\t100\t1e-40\t300\n",
                encoding="utf-8",
            )
            diamond_metadata.write_text(
                "database\tsubject_id\tdescription\tweight\n"
                "swissprot\tsp1\tATP synthase subunit alpha\t1.0\n"
                "reference\tref1\tATP synthase alpha subunit\t1.2\n"
                "swissprot\tsp2\tmalate dehydrogenase\t1.0\n"
                "swissprot\tsp3\tcellulose synthase\t1.0\n"
                "swissprot\tsp4\tnitrate reductase\t1.0\n",
                encoding="utf-8",
            )
            pfam_metadata.write_text(
                "name\taccession\tdescription\nPkinase\tPF00069\tProtein kinase domain\n",
                encoding="utf-8",
            )
            pfam_hits.write_text(
                "Pkinase PF00069.1 250 t3 - 300 1e-30 100 0 1 1 "
                "1e-25 1e-20 80 0 1 240 10 250 8 255 0.95 Protein kinase domain\n",
                encoding="utf-8",
            )

            apply_functional_annotations(
                annotation_table_path=annotation,
                diamond_hits_path=diamond_hits,
                diamond_metadata_path=diamond_metadata,
                pfam_hits_path=pfam_hits,
                pfam_metadata_path=pfam_metadata,
                output_path=output,
                evidence_path=evidence,
                log_path=base / "assign.log",
                metrics_path=base / "assign.json",
                config=FunctionalAnnotationConfig(enabled=True),
                name_standardization_summary_path=name_summary,
            )

            with output.open("r", encoding="utf-8", newline="") as handle:
                rows = {row["ID"]: row for row in csv.DictReader(handle, delimiter="\t")}
            self.assertIn(
                rows["t1"]["Description"],
                {
                    "ATP synthase subunit alpha",
                    "ATP synthase alpha subunit",
                },
            )
            self.assertEqual(rows["t2"]["Description"], "existing enzyme")
            self.assertEqual(rows["t3"]["Description"], "protein kinase domain-containing protein")
            self.assertEqual(rows["t4"]["Description"], "hypothetical protein")

            with evidence.open("r", encoding="utf-8", newline="") as handle:
                evidence_rows = {row["ID"]: row for row in csv.DictReader(handle, delimiter="\t")}
            self.assertIn(evidence_rows["t1"]["source"], {"swissprot", "reference"})
            self.assertEqual(evidence_rows["t2"]["source"], "existing")
            self.assertEqual(evidence_rows["t3"]["source"], "pfam")
            self.assertEqual(
                evidence_rows["t3"]["proposed_product"],
                "Protein kinase domain-containing protein",
            )
            self.assertEqual(
                evidence_rows["t3"]["name_standardization"],
                "lowercased_initial",
            )
            self.assertEqual(evidence_rows["t4"]["source"], "none")
            self.assertIn(
                "action\tlowercased_initial\t1",
                name_summary.read_text(encoding="utf-8"),
            )

    def test_pfam_fallback_skips_similarity_assignments_and_runs_shards(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            proteins = base / "proteins.fasta"
            hmm = base / "Pfam-A.hmm"
            cached_metadata = base / "Pfam-A.metadata.tsv"
            proteins.write_text(
                "".join(f">t{index}\nMPEPTIDE\n" for index in range(1, 7)),
                encoding="utf-8",
            )
            hmm.write_text("HMMER3/f\n", encoding="utf-8")
            cached_metadata.write_text(
                "name\taccession\tdescription\n",
                encoding="utf-8",
            )
            shard_queries: list[tuple[str, ...]] = []

            def fake_run_command(command: list[str], *, log_path: Path) -> None:
                query_path = Path(command[-1])
                query_ids = tuple(
                    line[1:].strip()
                    for line in query_path.read_text(encoding="utf-8").splitlines()
                    if line.startswith(">")
                )
                shard_queries.append(query_ids)
                Path(command[command.index("-o") + 1]).write_text("", encoding="utf-8")
                Path(command[command.index("--domtblout") + 1]).write_text(
                    "# fake hmmscan output\n",
                    encoding="utf-8",
                )
                log_path.write_text("fake hmmscan\n", encoding="utf-8")

            with (
                patch(
                    "msspack.functional_annotation._diamond_assignments",
                    return_value={"t1": object(), "t2": object()},
                ),
                patch(
                    "msspack.functional_annotation._prepare_pfam_database",
                    return_value=(hmm, cached_metadata, {"release": "test"}),
                ),
                patch(
                    "msspack.functional_annotation.run_command",
                    side_effect=fake_run_command,
                ),
            ):
                run_pfam_domain_search(
                    protein_fasta_path=proteins,
                    output_path=base / "pfam.domtblout",
                    metadata_path=base / "pfam-metadata.tsv",
                    provenance_path=base / "pfam-provenance.json",
                    log_path=base / "pfam.log",
                    metrics_path=base / "pfam.metrics.json",
                    config=FunctionalAnnotationConfig(threads=3),
                    base_dir=base,
                    cache_dir=base / "cache",
                    diamond_hits_path=base / "diamond.tsv",
                    diamond_metadata_path=base / "diamond-metadata.tsv",
                )

            self.assertEqual(len(shard_queries), 3)
            self.assertEqual(
                {identifier for shard in shard_queries for identifier in shard},
                {"t3", "t4", "t5", "t6"},
            )
            self.assertEqual(sorted(map(len, shard_queries)), [1, 1, 2])
            metrics = json.loads((base / "pfam.metrics.json").read_text(encoding="utf-8"))
            self.assertEqual(metrics["details"]["queries_scanned"], 4)
            self.assertEqual(metrics["details"]["similarity_assignments_skipped"], 2)
            self.assertEqual(metrics["details"]["hmmscan_shards"], 3)

    def test_cdd_assignments_prefer_specific_informative_domains(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            metadata = base / "cdd-metadata.tsv"
            results = base / "cdd.tsv"
            metadata.write_text(
                "pssm_id\taccession\tname\tdescription\tmodel_length\n"
                "101\tcd00001\tKinase\tProtein kinase domain\t250\n"
                "102\tcl00002\tLow_complexity\tLow complexity region\t80\n",
                encoding="utf-8",
            )
            results.write_text(
                "DATA\n"
                "SESSION\t1\tRPSBLAST\tCdd\tBLOSUM62\t0.01\n"
                "QUERY\tQuery_1\tPeptide\t300\tt1 translated protein\n"
                "DOMAINS\n"
                "1\tQuery_1\tSuperfamily\t102\t1\t80\t1e-30\t90\tcl00002\tLow_complexity\t-\t-\n"
                "1\tQuery_1\tSpecific\t101\t20\t260\t1e-40\t120\tcd00001\tKinase\t-\t200\n"
                "ENDDOMAINS\nENDQUERY\tQuery_1\nENDSESSION\t1\nENDDATA\n",
                encoding="utf-8",
            )

            assignments = _cdd_assignments(
                result_path=results,
                metadata_path=metadata,
                config=FunctionalAnnotationConfig(cdd_evalue=0.01),
            )

            self.assertEqual(assignments["t1"].source, "cdd")
            self.assertEqual(assignments["t1"].evidence_id, "cd00001")
            self.assertEqual(
                assignments["t1"].product,
                "Protein kinase domain-containing protein",
            )
            self.assertEqual(assignments["t1"].confidence, "moderate")


if __name__ == "__main__":
    unittest.main()
