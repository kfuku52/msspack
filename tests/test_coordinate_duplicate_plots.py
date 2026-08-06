import tempfile
import unittest
from pathlib import Path

from msspack.coordinate_duplicate_plots import (
    build_coordinate_duplicate_plot_data,
    coordinate_duplicate_plot_summary,
    write_coordinate_duplicate_pdf,
    write_coordinate_duplicate_plot_tsv,
    write_coordinate_duplicate_svg,
)
from msspack.coordinate_duplicates import (
    CoordinateDuplicatePair,
    read_coordinate_duplicate_map,
    write_coordinate_duplicate_map,
)


def _write_fixture(
    gff_path: Path,
    map_path: Path,
    *,
    removed_count: int,
) -> None:
    lines = [
        "##gff-version 3",
        "chr1\tsrc\tgene\t1\t100\t.\t+\t.\tID=kept",
        "chr1\tsrc\tmRNA\t1\t100\t.\t+\t.\tID=kept.t1;Parent=kept",
        "chr1\tsrc\texon\t1\t100\t.\t+\t.\tID=kept.ex1;Parent=kept.t1",
        "chr1\tsrc\tCDS\t1\t100\t.\t+\t0\tID=kept.cds1;Parent=kept.t1",
    ]
    pairs: list[CoordinateDuplicatePair] = []
    for index in range(1, removed_count + 1):
        removed_id = f"removed{index:02d}"
        lines.extend(
            [
                f"chr1\tsrc\tgene\t1\t100\t.\t+\t.\tID={removed_id}",
                (
                    f"chr1\tsrc\tmRNA\t1\t100\t.\t+\t.\t"
                    f"ID={removed_id}.t1;Parent={removed_id}"
                ),
                (
                    f"chr1\tsrc\texon\t1\t40\t.\t+\t.\t"
                    f"ID={removed_id}.ex1;Parent={removed_id}.t1"
                ),
                (
                    f"chr1\tsrc\tCDS\t1\t40\t.\t+\t0\t"
                    f"ID={removed_id}.cds1;Parent={removed_id}.t1"
                ),
                (
                    f"chr1\tsrc\texon\t60\t100\t.\t+\t.\t"
                    f"ID={removed_id}.ex2;Parent={removed_id}.t1"
                ),
                (
                    f"chr1\tsrc\tCDS\t60\t100\t.\t+\t0\t"
                    f"ID={removed_id}.cds2;Parent={removed_id}.t1"
                ),
            ]
        )
        pairs.append(
            CoordinateDuplicatePair(
                group_id="coordinate_duplicate_000001",
                seqid="chr1",
                start=1,
                end=100,
                strand="+",
                kept_gene_id="kept",
                removed_gene_id=removed_id,
            )
        )
    gff_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    write_coordinate_duplicate_map(map_path, pairs)


class CoordinateDuplicatePlotTests(unittest.TestCase):
    def test_legacy_duplicate_map_remains_readable(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            map_path = Path(tmp_dir) / "legacy.tsv"
            map_path.write_text(
                "group_id\tseqid\tstart\tend\tstrand\tkept_gene_id\tremoved_gene_id\n"
                "coordinate_duplicate_000001\tchr1\t1\t100\t+\tkept\tremoved\n",
                encoding="utf-8",
            )

            pair = read_coordinate_duplicate_map(map_path)[0]

            self.assertEqual(pair.kept_gene_id, "kept")
            self.assertEqual(pair.selection_policy, "first")
            self.assertEqual(pair.selection_reason, "input_order")
            self.assertEqual(pair.kept_cds_length, 0)
            self.assertFalse(pair.low_confidence)

    def test_default_limit_draws_first_50_but_tsv_keeps_all_rows(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            gff_path = base / "input.gff3"
            map_path = base / "duplicates.tsv"
            _write_fixture(gff_path, map_path, removed_count=55)
            data = build_coordinate_duplicate_plot_data(
                source_gff=gff_path,
                duplicate_map=map_path,
            )
            tsv_path = base / "plot.tsv"
            svg_path = base / "plot.svg"
            pdf_path = base / "plot.pdf"

            write_coordinate_duplicate_plot_tsv(data, tsv_path)
            write_coordinate_duplicate_svg(data, svg_path, limit=50)
            write_coordinate_duplicate_pdf(data, pdf_path, limit=50)

            tsv_text = tsv_path.read_text(encoding="utf-8")
            svg_text = svg_path.read_text(encoding="utf-8")
            pdf_text = pdf_path.read_bytes().decode("latin-1")
            self.assertIn("removed55", tsv_text)
            self.assertIn(
                'data-full-label="REMOVED removed50 / removed50.t1"',
                svg_text,
            )
            self.assertNotIn(
                'data-full-label="REMOVED removed51 / removed51.t1"',
                svg_text,
            )
            self.assertIn("Showing first 50 of 55 removed genes", svg_text)
            self.assertIn('text-anchor="end">100</text>', svg_text)
            self.assertGreater(pdf_text.count("/Type /Page"), 1)
            self.assertEqual(
                coordinate_duplicate_plot_summary(data, 50),
                {
                    "limit": 50,
                    "total_removed_genes": 55,
                    "shown_removed_genes": 50,
                    "truncated": True,
                },
            )

    def test_limit_can_be_increased(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            gff_path = base / "input.gff3"
            map_path = base / "duplicates.tsv"
            _write_fixture(gff_path, map_path, removed_count=55)
            data = build_coordinate_duplicate_plot_data(
                source_gff=gff_path,
                duplicate_map=map_path,
            )
            svg_path = base / "plot.svg"

            write_coordinate_duplicate_svg(data, svg_path, limit=55)

            svg_text = svg_path.read_text(encoding="utf-8")
            self.assertIn(
                'data-full-label="REMOVED removed55 / removed55.t1"',
                svg_text,
            )
            self.assertIn("Showing 55 removed genes", svg_text)


if __name__ == "__main__":
    unittest.main()
