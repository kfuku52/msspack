"""Review emitted MSS features and omit mRNAs that carry no additional annotation."""
from __future__ import annotations

import csv
from collections import defaultdict
from collections.abc import Iterable, Iterator
from dataclasses import dataclass, field
from pathlib import Path

from Bio.SeqFeature import AfterPosition, BeforePosition, ExactPosition, Location

from .trna_annotations import trna_amino_acids
from .utils import MSSPackError, atomic_text_writer

SUBMISSION_REVIEW_NAME = "submission-review.tsv"
_FEATURES_WITH_EXCEPTIONS = {"exception", "pseudo", "pseudogene", "trans_splicing"}


@dataclass(slots=True)
class _Block:
    feature: str
    location: str
    lines: list[str] = field(default_factory=list)
    qualifiers: dict[str, list[str]] = field(default_factory=dict)

    def value(self, qualifier: str) -> str:
        values = self.qualifiers.get(qualifier)
        return values[0] if values else ""


def _entries(lines: Iterable[str]) -> Iterator[tuple[str, list[_Block]]]:
    """Read up to one MSS entry at a time, retaining each original row."""
    entry = ""
    blocks: list[_Block] = []
    for number, line in enumerate(lines, 1):
        if not line.strip():
            if blocks:
                blocks[-1].lines.append(line)
            continue
        fields = line.rstrip("\r\n").split("\t")
        if not 3 <= len(fields) <= 5:
            raise MSSPackError(f"Expected three to five MSS annotation columns at line {number}")
        name, feature, location, qualifier, value = fields + [""] * (5 - len(fields))
        if name and name != entry:
            if blocks:
                yield entry, blocks
            entry, blocks = name, []
        if feature:
            blocks.append(_Block(feature, location))
        if not entry or not blocks:
            raise MSSPackError(f"MSS annotation row without an entry or feature at line {number}")
        blocks[-1].lines.append(line)
        if qualifier:
            blocks[-1].qualifiers.setdefault(qualifier, []).append(value)
    if blocks:
        yield entry, blocks


def _local_location(text: str) -> Location | None:
    try:
        location = Location.fromstring(text)
    except (TypeError, ValueError):
        return None
    if not all(
        not part.ref and part.strand in {1, -1}
        and type(part.start) in {ExactPosition, BeforePosition, AfterPosition}
        and type(part.end) in {ExactPosition, BeforePosition, AfterPosition}
        and int(part.start) < int(part.end)
        for part in location.parts
    ):
        return None
    return location


def _same_spliced_coverage(mrna: _Block, cds: _Block) -> bool:
    # Fuzzy CDS ends are retained even when the RNA has the same numeric boundaries.
    # Fuzzy mRNA ends, in contrast, are information carried only by the mRNA.
    rna_location = _local_location(mrna.location)
    cds_location = _local_location(cds.location)
    if rna_location is None or cds_location is None:
        return False
    if any(getattr(loc, "operator", "join") != "join" for loc in (rna_location, cds_location)):
        return False
    if any(
        type(part.start) is not ExactPosition or type(part.end) is not ExactPosition
        for part in rna_location.parts
    ):
        return False
    return [(int(p.start), int(p.end), p.strand) for p in rna_location.parts] == [
        (int(p.start), int(p.end), p.strand) for p in cds_location.parts
    ]


def omit_redundant_mrna(body: str) -> tuple[str, int]:
    """Omit same-locus mRNAs fully represented by one surviving CDS."""
    kept: list[str] = []
    omitted = 0
    for _entry, blocks in _entries(body.splitlines(keepends=True)):
        by_locus: dict[str, list[_Block]] = defaultdict(list)
        for block in blocks:
            locus = block.value("locus_tag")
            if locus and block.feature in {
                "mRNA", "CDS", "tRNA", "rRNA", "tmRNA", "ncRNA", "misc_RNA", "misc_feature",
            }:
                by_locus[locus].append(block)
        for block in blocks:
            if block.feature == "mRNA":
                group = by_locus.get(block.value("locus_tag"), [])
                if len(group) == 2:
                    other = group[0] if group[1] is block else group[1]
                    if (
                        other.feature == "CDS"
                        and not _FEATURES_WITH_EXCEPTIONS.intersection(other.qualifiers)
                        and all(
                            value in other.qualifiers.get(name, [])
                            for name, values in block.qualifiers.items()
                            for value in values
                        )
                        and _same_spliced_coverage(block, other)
                    ):
                        omitted += 1
                        continue
            kept.extend(block.lines)
    return "".join(kept), omitted


def _trna_mismatch(block: _Block) -> str:
    gene, product = block.value("gene"), block.value("product")
    amino_acids = trna_amino_acids(gene, product)
    if amino_acids is not None and amino_acids[0] != amino_acids[1]:
        return f"{gene} indicates tRNA-{amino_acids[0]}; product={product}"
    return ""


def _location_label(start: int, end: int, strand: int | None) -> str:
    span = str(end) if end - start == 1 else f"{start + 1}..{end}"
    return f"complement({span})" if strand == -1 else span


def _entry_warnings(
    entry: str, blocks: list[_Block],
) -> Iterator[tuple[str, str, str, str, str, str]]:
    exons: dict[tuple[int, int, int | None], set[str]] = defaultdict(set)
    for block in blocks:
        locus = block.value("locus_tag")
        if block.feature == "tRNA" and (detail := _trna_mismatch(block)):
            yield "trna_gene_product_mismatch", entry, "tRNA", locus, block.location, detail
        if block.feature != "CDS":
            continue
        location = _local_location(block.location)
        if location is None:
            yield ("unreviewed_cds_location", entry, "CDS", locus, block.location,
                   "Review the CDS location with the DDBJ validators")
            continue
        if {"pseudo", "pseudogene"}.intersection(block.qualifiers):
            continue
        codon_start = block.value("codon_start") or "1"
        if codon_start not in {"1", "2", "3"}:
            yield ("unreviewed_codon_start", entry, "CDS", locus, block.location,
                   f"Unexpected codon_start={codon_start}")
            continue
        coding_bp = len(location) - int(codon_start) + 1
        partial = codon_start != "1" or any(
            type(part.start) is not ExactPosition or type(part.end) is not ExactPosition
            for part in location.parts
        )
        if coding_bp < 93:
            description = f"{coding_bp} coding bp; review against the 93 bp threshold"
            if partial:
                description += " (partial CDS)"
            yield "short_cds", entry, "CDS", locus, block.location, description
        for part in location.parts:
            start, end = int(part.start), int(part.end)
            if end - start == 1:
                yield ("single_base_cds_segment", entry, "CDS", locus,
                       _location_label(start, end, part.strand), f"CDS {block.location}")
            if (locus and type(part.start) is ExactPosition and type(part.end) is ExactPosition
                    and part.strand in {1, -1}):
                exons[start, end, part.strand].add(locus)
    for (start, end, strand), loci in exons.items():
        if len(loci) > 1:
            yield ("shared_cds_segment", entry, "CDS", ",".join(sorted(loci)),
                   _location_label(start, end, strand), "Exact same-strand segment at distinct loci")


def write_submission_review(annotation_path: Path, output_path: Path) -> int:
    """Write an advisory TSV from an MSS annotation; return the warning count."""
    if annotation_path.resolve() == output_path.resolve():
        raise MSSPackError("The review output must differ from the annotation input")
    if not annotation_path.is_file():
        raise MSSPackError(f"MSS annotation not found: {annotation_path}")
    warnings = 0
    is_wgs = False
    common_coverage = False
    entry_count = 0
    missing_coverage: list[str] = []
    with annotation_path.open(encoding="utf-8") as src, atomic_text_writer(output_path) as dst:
        writer = csv.writer(dst, delimiter="\t", lineterminator="\n")
        writer.writerow(["warning", "entry", "feature", "locus_tag", "location", "detail"])
        for entry, blocks in _entries(src):
            coverage = any(
                b.value("Genome Coverage") for b in blocks if b.feature == "ST_COMMENT"
            )
            if entry == "COMMON":
                is_wgs = any(b.feature == "DATATYPE" and b.value("type") == "WGS" for b in blocks)
                common_coverage = coverage
                continue
            entry_count += 1
            if is_wgs and not (common_coverage or coverage):
                missing_coverage.append(entry)
            for row in _entry_warnings(entry, blocks):
                writer.writerow(row)
                warnings += 1
        if missing_coverage and len(missing_coverage) == entry_count:
            missing_coverage = ["COMMON"]
        for entry in missing_coverage:
            writer.writerow(("missing_genome_coverage", entry, "ST_COMMENT", "", "",
                             "Set Genome Coverage to a depth such as 130x or to Unknown"))
            warnings += 1
    return warnings
