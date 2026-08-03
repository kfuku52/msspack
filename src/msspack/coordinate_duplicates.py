from __future__ import annotations

import csv
from collections import Counter, defaultdict
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path

from Bio.Data import CodonTable
from Bio.Seq import Seq

from .fasta import iter_fasta, reverse_complement
from .gff import GFFRecord, child_ids
from .utils import MSSPackError, atomic_text_writer

COORDINATE_DUPLICATE_POLICIES = {
    "first",
    "keep_all",
    "longest_valid_cds",
}

LEGACY_DUPLICATE_MAP_FIELDS = (
    "group_id",
    "seqid",
    "start",
    "end",
    "strand",
    "kept_gene_id",
    "removed_gene_id",
)

DUPLICATE_MAP_FIELDS = (
    *LEGACY_DUPLICATE_MAP_FIELDS,
    "kept_transcript_id",
    "removed_transcript_id",
    "kept_cds_length",
    "removed_cds_length",
    "kept_intron_count",
    "removed_intron_count",
    "kept_internal_stops",
    "removed_internal_stops",
    "kept_complete",
    "removed_complete",
    "kept_ambiguous_amino_acids",
    "removed_ambiguous_amino_acids",
    "kept_structurally_valid",
    "removed_structurally_valid",
    "kept_splice_motifs",
    "removed_splice_motifs",
    "selection_policy",
    "selection_reason",
    "low_confidence",
)


@dataclass(frozen=True)
class CoordinateGene:
    gene_id: str
    seqid: str
    start: int
    end: int
    strand: str
    input_order: int = 0
    transcript_id: str = ""
    cds_length: int = 0
    intron_count: int = 0
    internal_stop_count: int = 0
    complete_boundary_count: int = 0
    ambiguous_amino_acid_count: int = 0
    structurally_valid: bool = False
    splice_motifs: tuple[str, ...] = ()
    protein_sequence: str = ""

    @property
    def coordinate_key(self) -> tuple[str, int, int, str]:
        return self.seqid, self.start, self.end, self.strand

    @property
    def complete(self) -> bool:
        return self.complete_boundary_count == 2


@dataclass(frozen=True)
class CoordinateDuplicatePair:
    group_id: str
    seqid: str
    start: int
    end: int
    strand: str
    kept_gene_id: str
    removed_gene_id: str
    kept_transcript_id: str = ""
    removed_transcript_id: str = ""
    kept_cds_length: int = 0
    removed_cds_length: int = 0
    kept_intron_count: int = 0
    removed_intron_count: int = 0
    kept_internal_stops: int = 0
    removed_internal_stops: int = 0
    kept_complete: bool = False
    removed_complete: bool = False
    kept_ambiguous_amino_acids: int = 0
    removed_ambiguous_amino_acids: int = 0
    kept_structurally_valid: bool = False
    removed_structurally_valid: bool = False
    kept_splice_motifs: tuple[str, ...] = ()
    removed_splice_motifs: tuple[str, ...] = ()
    selection_policy: str = "first"
    selection_reason: str = "input_order"
    low_confidence: bool = False


@dataclass(frozen=True)
class _TranscriptModel:
    gene: CoordinateGene
    transcript_id: str
    input_order: int
    cds_records: tuple[GFFRecord, ...]


def _candidate_rank(gene: CoordinateGene) -> tuple[int, int, int, int, int, int, int]:
    return (
        int(gene.structurally_valid),
        -gene.internal_stop_count,
        gene.complete_boundary_count,
        -gene.ambiguous_amino_acid_count,
        gene.cds_length,
        -gene.intron_count,
        -gene.input_order,
    )


def _selection_reason(kept: CoordinateGene, removed: CoordinateGene) -> str:
    comparisons = (
        ("valid_structure", int(kept.structurally_valid), int(removed.structurally_valid)),
        ("fewer_internal_stops", -kept.internal_stop_count, -removed.internal_stop_count),
        (
            "more_complete_boundaries",
            kept.complete_boundary_count,
            removed.complete_boundary_count,
        ),
        (
            "fewer_ambiguous_amino_acids",
            -kept.ambiguous_amino_acid_count,
            -removed.ambiguous_amino_acid_count,
        ),
        ("longer_cds", kept.cds_length, removed.cds_length),
        ("fewer_introns", -kept.intron_count, -removed.intron_count),
    )
    for reason, kept_value, removed_value in comparisons:
        if kept_value != removed_value:
            return reason
    return "input_order"


def _pair_from_genes(
    *,
    group_id: str,
    kept: CoordinateGene,
    removed: CoordinateGene,
    policy: str,
) -> CoordinateDuplicatePair:
    reason = _selection_reason(kept, removed) if policy == "longest_valid_cds" else "input_order"
    low_confidence = (
        kept.structurally_valid
        and removed.structurally_valid
        and kept.internal_stop_count == 0
        and removed.internal_stop_count == 0
        and kept.protein_sequence != removed.protein_sequence
        and reason in {"longer_cds", "fewer_introns", "input_order"}
    )
    return CoordinateDuplicatePair(
        group_id=group_id,
        seqid=kept.seqid,
        start=kept.start,
        end=kept.end,
        strand=kept.strand,
        kept_gene_id=kept.gene_id,
        removed_gene_id=removed.gene_id,
        kept_transcript_id=kept.transcript_id,
        removed_transcript_id=removed.transcript_id,
        kept_cds_length=kept.cds_length,
        removed_cds_length=removed.cds_length,
        kept_intron_count=kept.intron_count,
        removed_intron_count=removed.intron_count,
        kept_internal_stops=kept.internal_stop_count,
        removed_internal_stops=removed.internal_stop_count,
        kept_complete=kept.complete,
        removed_complete=removed.complete,
        kept_ambiguous_amino_acids=kept.ambiguous_amino_acid_count,
        removed_ambiguous_amino_acids=removed.ambiguous_amino_acid_count,
        kept_structurally_valid=kept.structurally_valid,
        removed_structurally_valid=removed.structurally_valid,
        kept_splice_motifs=kept.splice_motifs,
        removed_splice_motifs=removed.splice_motifs,
        selection_policy=policy,
        selection_reason=reason,
        low_confidence=low_confidence,
    )


def identify_coordinate_duplicate_pairs(
    genes: Iterable[CoordinateGene],
    *,
    policy: str = "first",
) -> tuple[CoordinateDuplicatePair, ...]:
    if policy not in COORDINATE_DUPLICATE_POLICIES:
        raise MSSPackError(f"Unsupported coordinate duplicate policy: {policy}")
    if policy == "keep_all":
        return ()

    genes_by_coordinate: dict[tuple[str, int, int, str], list[CoordinateGene]] = {}
    for gene in genes:
        genes_by_coordinate.setdefault(gene.coordinate_key, []).append(gene)

    pairs: list[CoordinateDuplicatePair] = []
    group_number = 0
    for coordinate_genes in genes_by_coordinate.values():
        if len(coordinate_genes) < 2:
            continue
        group_number += 1
        group_id = f"coordinate_duplicate_{group_number:06d}"
        kept = (
            coordinate_genes[0]
            if policy == "first"
            else max(coordinate_genes, key=_candidate_rank)
        )
        for removed in coordinate_genes:
            if removed is kept:
                continue
            pairs.append(
                _pair_from_genes(
                    group_id=group_id,
                    kept=kept,
                    removed=removed,
                    policy=policy,
                )
            )
    return tuple(pairs)


def _transcript_sequence(
    model: _TranscriptModel,
    genome_sequence: str,
) -> tuple[str, int, tuple[str, ...], bool]:
    ordered = tuple(sorted(model.cds_records, key=lambda record: (record.start, record.end)))
    valid = (
        bool(ordered)
        and model.gene.strand in {"+", "-"}
        and 1 <= model.gene.start <= model.gene.end <= len(genome_sequence)
    )
    previous_end = 0
    for record in ordered:
        valid = (
            valid
            and record.seqid == model.gene.seqid
            and record.strand == model.gene.strand
            and model.gene.start <= record.start <= record.end <= model.gene.end
            and record.phase in {"0", "1", "2"}
            and (previous_end == 0 or record.start > previous_end)
        )
        previous_end = record.end

    nucleotide_sequence = "".join(
        genome_sequence[record.start - 1 : record.end] for record in ordered
    )
    if model.gene.strand == "-":
        nucleotide_sequence = reverse_complement(nucleotide_sequence)

    splice_motifs: list[str] = []
    for left, right in zip(ordered, ordered[1:], strict=False):
        intron_start = left.end + 1
        intron_end = right.start - 1
        if intron_start > intron_end:
            continue
        intron = genome_sequence[intron_start - 1 : intron_end]
        if model.gene.strand == "-":
            intron = reverse_complement(intron)
        motif = (
            f"{intron[:2].upper()}-{intron[-2:].upper()}"
            if len(intron) >= 4
            else intron.upper()
        )
        splice_motifs.append(motif)
    if model.gene.strand == "-":
        splice_motifs.reverse()

    intron_count = len(splice_motifs)
    valid = valid and bool(nucleotide_sequence) and len(nucleotide_sequence) % 3 == 0
    return nucleotide_sequence.upper(), intron_count, tuple(splice_motifs), valid


def _evaluate_transcript(
    model: _TranscriptModel,
    genome_sequence: str,
    genetic_code: int,
) -> CoordinateGene:
    sequence, intron_count, splice_motifs, structurally_valid = _transcript_sequence(
        model,
        genome_sequence,
    )
    codon_table = CodonTable.unambiguous_dna_by_id[genetic_code]
    start_complete = len(sequence) >= 3 and sequence[:3] in codon_table.start_codons
    stop_complete = len(sequence) >= 3 and sequence[-3:] in codon_table.stop_codons
    protein = ""
    if sequence:
        usable_length = len(sequence) - (len(sequence) % 3)
        protein = str(Seq(sequence[:usable_length]).translate(table=genetic_code))
    internal_protein = protein[:-1] if stop_complete and protein.endswith("*") else protein
    internal_stop_count = internal_protein.count("*")
    protein_sequence = protein[:-1] if protein.endswith("*") else protein
    return CoordinateGene(
        gene_id=model.gene.gene_id,
        seqid=model.gene.seqid,
        start=model.gene.start,
        end=model.gene.end,
        strand=model.gene.strand,
        input_order=model.gene.input_order,
        transcript_id=model.transcript_id,
        cds_length=len(sequence),
        intron_count=intron_count,
        internal_stop_count=internal_stop_count,
        complete_boundary_count=int(start_complete) + int(stop_complete),
        ambiguous_amino_acid_count=protein_sequence.count("X"),
        structurally_valid=structurally_valid,
        splice_motifs=splice_motifs,
        protein_sequence=protein_sequence,
    )


def build_coordinate_gene_candidates(
    *,
    records: Sequence[GFFRecord],
    fasta_path: Path,
    genetic_code: str,
) -> tuple[CoordinateGene, ...]:
    try:
        genetic_code_number = int(genetic_code)
        CodonTable.unambiguous_dna_by_id[genetic_code_number]
    except (KeyError, ValueError) as exc:
        raise MSSPackError(f"Unsupported NCBI genetic code: {genetic_code}") from exc

    genes: list[CoordinateGene] = []
    gene_records: dict[str, GFFRecord] = {}
    gene_order: dict[str, int] = {}
    coordinate_counts: Counter[tuple[str, int, int, str]] = Counter()
    for record in records:
        if record.type != "gene":
            continue
        gene_id = record.attributes.get("ID", "")
        if not gene_id:
            continue
        gene = CoordinateGene(
            gene_id=gene_id,
            seqid=record.seqid,
            start=record.start,
            end=record.end,
            strand=record.strand,
            input_order=len(genes),
        )
        genes.append(gene)
        gene_records[gene_id] = record
        gene_order[gene_id] = gene.input_order
        coordinate_counts[gene.coordinate_key] += 1

    duplicate_gene_ids = {
        gene.gene_id for gene in genes if coordinate_counts[gene.coordinate_key] > 1
    }
    if not duplicate_gene_ids:
        return tuple(genes)

    transcripts_by_gene: dict[str, list[str]] = defaultdict(list)
    transcript_order: dict[str, int] = {}
    cds_by_parent: dict[str, list[GFFRecord]] = defaultdict(list)
    for record_order, record in enumerate(records):
        record_id = record.attributes.get("ID", "")
        if record.type in {"mRNA", "transcript"} and record_id:
            transcript_order.setdefault(record_id, record_order)
            for gene_id in child_ids(record.attributes.get("Parent")):
                if gene_id in duplicate_gene_ids and record_id not in transcripts_by_gene[gene_id]:
                    transcripts_by_gene[gene_id].append(record_id)
        elif record.type == "CDS":
            for parent_id in child_ids(record.attributes.get("Parent")):
                cds_by_parent[parent_id].append(record)

    models_by_seqid: dict[str, list[_TranscriptModel]] = defaultdict(list)
    for gene in genes:
        if gene.gene_id not in duplicate_gene_ids:
            continue
        transcript_ids = list(transcripts_by_gene.get(gene.gene_id, []))
        if cds_by_parent.get(gene.gene_id):
            transcript_ids.append(gene.gene_id)
        for transcript_id in transcript_ids:
            cds_records = tuple(cds_by_parent.get(transcript_id, []))
            if not cds_records:
                continue
            models_by_seqid[gene.seqid].append(
                _TranscriptModel(
                    gene=gene,
                    transcript_id=transcript_id,
                    input_order=transcript_order.get(
                        transcript_id,
                        gene_order[gene.gene_id],
                    ),
                    cds_records=cds_records,
                )
            )

    evaluations_by_gene: dict[str, list[tuple[int, CoordinateGene]]] = defaultdict(list)
    found_seqids: set[str] = set()
    for fasta_record in iter_fasta(fasta_path):
        models = models_by_seqid.get(fasta_record.id, [])
        if not models:
            continue
        found_seqids.add(fasta_record.id)
        for model in models:
            evaluations_by_gene[model.gene.gene_id].append(
                (
                    model.input_order,
                    _evaluate_transcript(
                        model,
                        fasta_record.sequence,
                        genetic_code_number,
                    ),
                )
            )
    missing_seqids = set(models_by_seqid) - found_seqids
    if missing_seqids:
        shown = ", ".join(sorted(missing_seqids)[:10])
        raise MSSPackError(
            f"Coordinate duplicate genes reference seqids missing from FASTA: {shown}"
        )

    summarized: list[CoordinateGene] = []
    for gene in genes:
        evaluations = evaluations_by_gene.get(gene.gene_id)
        if not evaluations:
            summarized.append(gene)
            continue
        evaluations.sort(key=lambda item: item[0])
        summarized.append(max((item[1] for item in evaluations), key=_candidate_rank))
    return tuple(summarized)


def write_coordinate_duplicate_map(
    path: Path,
    pairs: Iterable[CoordinateDuplicatePair],
) -> Path:
    with atomic_text_writer(path) as handle:
        writer = csv.DictWriter(handle, fieldnames=DUPLICATE_MAP_FIELDS, delimiter="\t")
        writer.writeheader()
        for pair in pairs:
            writer.writerow(
                {
                    "group_id": pair.group_id,
                    "seqid": pair.seqid,
                    "start": pair.start,
                    "end": pair.end,
                    "strand": pair.strand,
                    "kept_gene_id": pair.kept_gene_id,
                    "removed_gene_id": pair.removed_gene_id,
                    "kept_transcript_id": pair.kept_transcript_id,
                    "removed_transcript_id": pair.removed_transcript_id,
                    "kept_cds_length": pair.kept_cds_length,
                    "removed_cds_length": pair.removed_cds_length,
                    "kept_intron_count": pair.kept_intron_count,
                    "removed_intron_count": pair.removed_intron_count,
                    "kept_internal_stops": pair.kept_internal_stops,
                    "removed_internal_stops": pair.removed_internal_stops,
                    "kept_complete": str(pair.kept_complete).lower(),
                    "removed_complete": str(pair.removed_complete).lower(),
                    "kept_ambiguous_amino_acids": pair.kept_ambiguous_amino_acids,
                    "removed_ambiguous_amino_acids": pair.removed_ambiguous_amino_acids,
                    "kept_structurally_valid": str(pair.kept_structurally_valid).lower(),
                    "removed_structurally_valid": str(
                        pair.removed_structurally_valid
                    ).lower(),
                    "kept_splice_motifs": ",".join(pair.kept_splice_motifs),
                    "removed_splice_motifs": ",".join(pair.removed_splice_motifs),
                    "selection_policy": pair.selection_policy,
                    "selection_reason": pair.selection_reason,
                    "low_confidence": str(pair.low_confidence).lower(),
                }
            )
    return path


def _parse_int(row: dict[str, str], field: str, line_number: int, path: Path) -> int:
    try:
        return int(row.get(field, "0") or "0")
    except ValueError as exc:
        raise MSSPackError(
            f"Invalid integer {field} at {path}:{line_number}"
        ) from exc


def _parse_bool(row: dict[str, str], field: str, line_number: int, path: Path) -> bool:
    value = (row.get(field, "false") or "false").casefold()
    if value not in {"true", "false"}:
        raise MSSPackError(f"Invalid boolean {field} at {path}:{line_number}")
    return value == "true"


def read_coordinate_duplicate_map(path: Path) -> tuple[CoordinateDuplicatePair, ...]:
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fieldnames = tuple(reader.fieldnames or ())
            if fieldnames not in {LEGACY_DUPLICATE_MAP_FIELDS, DUPLICATE_MAP_FIELDS}:
                raise MSSPackError(f"Invalid coordinate duplicate map header: {path}")
            legacy = fieldnames == LEGACY_DUPLICATE_MAP_FIELDS
            pairs: list[CoordinateDuplicatePair] = []
            for line_number, row in enumerate(reader, start=2):
                pair = CoordinateDuplicatePair(
                    group_id=row.get("group_id", ""),
                    seqid=row.get("seqid", ""),
                    start=_parse_int(row, "start", line_number, path),
                    end=_parse_int(row, "end", line_number, path),
                    strand=row.get("strand", ""),
                    kept_gene_id=row.get("kept_gene_id", ""),
                    removed_gene_id=row.get("removed_gene_id", ""),
                    kept_transcript_id=row.get("kept_transcript_id", ""),
                    removed_transcript_id=row.get("removed_transcript_id", ""),
                    kept_cds_length=_parse_int(row, "kept_cds_length", line_number, path),
                    removed_cds_length=_parse_int(
                        row,
                        "removed_cds_length",
                        line_number,
                        path,
                    ),
                    kept_intron_count=_parse_int(
                        row,
                        "kept_intron_count",
                        line_number,
                        path,
                    ),
                    removed_intron_count=_parse_int(
                        row,
                        "removed_intron_count",
                        line_number,
                        path,
                    ),
                    kept_internal_stops=_parse_int(
                        row,
                        "kept_internal_stops",
                        line_number,
                        path,
                    ),
                    removed_internal_stops=_parse_int(
                        row,
                        "removed_internal_stops",
                        line_number,
                        path,
                    ),
                    kept_complete=_parse_bool(row, "kept_complete", line_number, path),
                    removed_complete=_parse_bool(
                        row,
                        "removed_complete",
                        line_number,
                        path,
                    ),
                    kept_ambiguous_amino_acids=_parse_int(
                        row,
                        "kept_ambiguous_amino_acids",
                        line_number,
                        path,
                    ),
                    removed_ambiguous_amino_acids=_parse_int(
                        row,
                        "removed_ambiguous_amino_acids",
                        line_number,
                        path,
                    ),
                    kept_structurally_valid=_parse_bool(
                        row,
                        "kept_structurally_valid",
                        line_number,
                        path,
                    ),
                    removed_structurally_valid=_parse_bool(
                        row,
                        "removed_structurally_valid",
                        line_number,
                        path,
                    ),
                    kept_splice_motifs=tuple(
                        filter(None, row.get("kept_splice_motifs", "").split(","))
                    ),
                    removed_splice_motifs=tuple(
                        filter(None, row.get("removed_splice_motifs", "").split(","))
                    ),
                    selection_policy=(
                        "first" if legacy else row.get("selection_policy", "")
                    ),
                    selection_reason=(
                        "input_order" if legacy else row.get("selection_reason", "")
                    ),
                    low_confidence=_parse_bool(
                        row,
                        "low_confidence",
                        line_number,
                        path,
                    ),
                )
                if (
                    not pair.group_id
                    or not pair.seqid
                    or not pair.kept_gene_id
                    or not pair.removed_gene_id
                    or pair.start > pair.end
                    or pair.selection_policy not in COORDINATE_DUPLICATE_POLICIES
                ):
                    raise MSSPackError(
                        f"Invalid coordinate duplicate map row at {path}:{line_number}"
                    )
                pairs.append(pair)
    except OSError as exc:
        raise MSSPackError(f"Could not read coordinate duplicate map: {path}") from exc
    return tuple(pairs)
