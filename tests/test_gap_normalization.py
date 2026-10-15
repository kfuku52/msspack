import tempfile
import unittest
from pathlib import Path

from msspack.gap_normalization import normalize_gap_lengths
from msspack.gff import iter_gff_records
from msspack.utils import MSSPackError


class GapNormalizationTests(unittest.TestCase):
    def test_gapjust_maps_unsorted_features_across_expanding_and_shrinking_gaps(self) -> None:
        # Three gaps of lengths 4, 1 and 6 become length 2. Include both gap
        # boundaries, clipped suffixes, adjacent bases and a different contig.
        coordinates = [23, 1, 22, 3, 4, 5, 6, 7, 8, 10, 11, 12, 16, 17, 18, 19]
        expected = [18, 1, 17, 3, 4, 5, 5, 5, 6, 8, 9, 11, 15, 16, 17, 17]
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            (base / "input.fa").write_text(
                ">chr1\nAAANNNNAAANAAAAANNNNNNA\n>chr2\nAAAA\n",
            )
            rows = [
                f"chr1\t.\tregion\t{value}\t{value}\t.\t-\t.\tID=r{index}"
                for index, value in enumerate(coordinates)
            ]
            rows.insert(3, "chr2\t.\tgene\t2\t4\t.\t+\t.\tID=unaffected")
            (base / "input.gff3").write_text("\n".join(rows) + "\n")
            normalize_gap_lengths(
                fasta_path=base / "input.fa", output_fasta_path=base / "out.fa",
                input_gff_path=base / "input.gff3", output_gff_path=base / "out.gff3",
                log_path=base / "log", gap_len=2,
            )
            records = list(iter_gff_records(base / "out.gff3"))
            unchanged = records.pop(3)
            self.assertEqual((unchanged.seqid, unchanged.start, unchanged.end), ("chr2", 2, 4))
            self.assertEqual([(r.start, r.end) for r in records], [(v, v) for v in expected])

    def test_gapjust_maps_boundaries_inside_gaps_and_after_multiple_edits(self) -> None:
        cases = (
            ("ATG" + "N" * 120 + "AAA", 100, [(1, 15), (110, 126)], [(1, 15), (103, 106)]),
            ("AAA" + "N" * 4 + "AAA", 6, [(5, 6), (8, 10)], [(5, 6), (10, 12)]),
            ("AAA" + "N" * 4 + "AAA" + "N" * 4 + "AAA", 2,
             [(6, 14), (15, 17)], [(5, 10), (11, 13)]),
            ("AAA" + "N" * 4 + "AAA", 0, [(1, 6), (5, 10)], [(1, 3), (4, 6)]),
        )
        for sequence, length, before, expected in cases:
            with self.subTest(length=length, sequence=sequence), tempfile.TemporaryDirectory() as tmp:
                base = Path(tmp)
                fasta, gff = base / "input.fa", base / "input.gff3"
                fasta.write_text(f">chr1\n{sequence}\n")
                gff.write_text("".join(
                    f"chr1\t.\tgene\t{start}\t{end}\t.\t+\t.\tID=g{index}\n"
                    for index, (start, end) in enumerate(before)
                ))
                output = base / "out.gff3"
                normalize_gap_lengths(
                    fasta_path=fasta, output_fasta_path=base / "out.fa", log_path=base / "log",
                    gap_len=length, input_gff_path=gff, output_gff_path=output,
                )
                self.assertEqual([(r.start, r.end) for r in iter_gff_records(output)], expected)

    def test_gap_deletion_rejects_features_with_no_remaining_bases(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            fasta, gff = base / "input.fa", base / "input.gff3"
            fasta.write_text(">chr1\nAAANNNNAAA\n")
            gff.write_text("chr1\t.\tgene\t4\t7\t.\t+\t.\tID=g1\n")
            with self.assertRaisesRegex(MSSPackError, "remove every base.*chr1:4..7"):
                normalize_gap_lengths(
                    fasta_path=fasta, output_fasta_path=base / "out.fa", log_path=base / "log",
                    gap_len=0, input_gff_path=gff, output_gff_path=base / "out.gff3",
                )

    def test_gapjust_updates_gff_coordinates(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            fasta = base / "input.fa"
            gff = base / "input.gff3"
            out_fasta = base / "out.fa"
            out_gff = base / "out.gff3"
            log = base / "out.log"
            fasta.write_text(">chr1\nAAANNNNAAA\n", encoding="utf-8")
            gff.write_text(
                "chr1\tsrc\tgene\t8\t10\t.\t+\t.\tID=g1\n",
                encoding="utf-8",
            )

            normalize_gap_lengths(
                fasta_path=fasta,
                output_fasta_path=out_fasta,
                log_path=log,
                gap_len=2,
                gap_just_min=0,
                gap_just_max=10,
                input_gff_path=gff,
                output_gff_path=out_gff,
            )

            self.assertIn("AAANNAAA", out_fasta.read_text(encoding="utf-8"))
            self.assertIn("\t6\t8\t", out_gff.read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
