import unittest
from types import SimpleNamespace

from msspack.mss_converter.core import build_gff_feature_indexes, detect_gap_regions


class MssConverterTests(unittest.TestCase):
    def test_detect_gap_regions_finds_contiguous_n_runs(self) -> None:
        record = SimpleNamespace(seq="AANNNTAANNNNN")
        out, gaps = detect_gap_regions(
            record,
            linkage_evidence="paired-ends",
            min_assembly_gap_size=3,
            gap_estimated_length="unknown",
        )
        self.assertEqual(gaps, [(3, 5), (9, 13)])
        self.assertIn("assembly_gap\t3..5", out)
        self.assertIn("assembly_gap\t9..13", out)

    def test_build_gff_feature_indexes_normalizes_parent_and_phase(self) -> None:
        rows = [
            {
                "seq_id": "chr1",
                "type": "gene",
                "ID": "Gene1",
                "Parent": None,
                "start": 10,
                "end": 50,
                "strand": "+",
                "phase": ".",
                "Name": None,
                "Type": None,
                "anticodon": None,
            },
            {
                "seq_id": "chr1",
                "type": "mRNA",
                "ID": "Gene1-T1",
                "Parent": "Gene1",
                "start": 10,
                "end": 50,
                "strand": "+",
                "phase": ".",
                "Name": "Gene1-T1",
                "Type": None,
                "anticodon": None,
            },
            {
                "seq_id": "chr1",
                "type": "CDS",
                "ID": "Gene1-T1.cds1",
                "Parent": "Gene1-T1",
                "start": 20,
                "end": 40,
                "strand": "+",
                "phase": "2",
                "Name": None,
                "Type": None,
                "anticodon": None,
            },
        ]

        gene_lookup, parent_lookup = build_gff_feature_indexes(rows)

        self.assertEqual(gene_lookup["chr1"][0].id, "Gene1")
        self.assertEqual(parent_lookup["Gene1"][0].id, "Gene1-T1")
        self.assertEqual(parent_lookup["Gene1-T1"][0].phase, 2)
        self.assertEqual(parent_lookup["Gene1-T1"][0].parent, "Gene1-T1")


if __name__ == "__main__":
    unittest.main()
