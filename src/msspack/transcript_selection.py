from __future__ import annotations

from collections import OrderedDict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, TypedDict

from .gff import parse_attributes
from .step_logging import write_id_list


@dataclass
class _TranscriptData:
    gene_id: str | None = None
    mrna_lines: list[str] = field(default_factory=list)
    cds_lines: list[str] = field(default_factory=list)
    five_prime_utr_lines: list[str] = field(default_factory=list)
    three_prime_utr_lines: list[str] = field(default_factory=list)
    exon_lines: list[str] = field(default_factory=list)
    total_cds: int = 0
    total_utr: int = 0
    sort_key: tuple[int, int, str] | None = None


@dataclass
class _GeneData:
    gene_lines: list[str] = field(default_factory=list)
    transcript_ids: list[str] = field(default_factory=list)


class TranscriptSelectionSummary(TypedDict):
    removed_mrnas: int
    changed_genes: int
    changed_gene_ids: list[str]


def _feature_length(start_str: str, end_str: str) -> int:
    return abs(int(end_str) - int(start_str)) + 1


def _ensure_gene(
    genes: OrderedDict[str, _GeneData],
    gene_id: str,
) -> _GeneData:
    if gene_id not in genes:
        genes[gene_id] = _GeneData()
    return genes[gene_id]


def _ensure_transcript(
    transcripts: OrderedDict[str, _TranscriptData],
    transcript_id: str,
) -> _TranscriptData:
    if transcript_id not in transcripts:
        transcripts[transcript_id] = _TranscriptData()
    return transcripts[transcript_id]


def select_one_mrna_per_gene(
    *,
    input_gff_path: str | Path,
    output_gff_path: str | Path,
    changed_gene_ids_path: str | Path | None = None,
) -> TranscriptSelectionSummary:
    genes: OrderedDict[str, _GeneData] = OrderedDict()
    transcripts: OrderedDict[str, _TranscriptData] = OrderedDict()

    with Path(input_gff_path).open("r", encoding="utf-8") as infile:
        for raw_line in infile:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue

            cols = line.split("\t")
            if len(cols) < 9:
                continue

            seqid, source, ftype, start, end, score, strand, phase, attributes_str = cols
            attributes = parse_attributes(attributes_str)

            if ftype == "gene":
                gene_id = attributes.get("ID")
                if gene_id is None:
                    continue
                gene_data = _ensure_gene(genes, gene_id)
                gene_data.gene_lines.append(line)
                continue

            if ftype in ("mRNA", "transcript"):
                gene_id = attributes.get("Parent")
                transcript_id = attributes.get("ID")
                if gene_id is None or transcript_id is None:
                    continue
                gene_data = _ensure_gene(genes, gene_id)
                if transcript_id not in gene_data.transcript_ids:
                    gene_data.transcript_ids.append(transcript_id)
                transcript_data = _ensure_transcript(transcripts, transcript_id)
                if transcript_data.gene_id is None:
                    transcript_data.gene_id = gene_id
                transcript_data.mrna_lines.append(line)
                candidate_key = (int(start), int(end), transcript_id)
                sort_key = transcript_data.sort_key
                if sort_key is None or candidate_key < sort_key:
                    transcript_data.sort_key = candidate_key
                continue

            if ftype not in ("CDS", "five_prime_UTR", "three_prime_UTR", "exon"):
                continue

            parent_str = attributes.get("Parent")
            if parent_str is None:
                continue
            for transcript_id in parent_str.split(","):
                transcript_data = _ensure_transcript(transcripts, transcript_id)
                if transcript_data.sort_key is None:
                    transcript_data.sort_key = (int(start), int(end), transcript_id)

                if ftype == "CDS":
                    transcript_data.cds_lines.append(line)
                    transcript_data.total_cds += _feature_length(start, end)
                elif ftype == "five_prime_UTR":
                    transcript_data.five_prime_utr_lines.append(line)
                    transcript_data.total_utr += _feature_length(start, end)
                elif ftype == "three_prime_UTR":
                    transcript_data.three_prime_utr_lines.append(line)
                    transcript_data.total_utr += _feature_length(start, end)
                elif ftype == "exon":
                    transcript_data.exon_lines.append(line)

    removed_mrna_count = 0
    changed_genes_count = 0
    changed_gene_ids: list[str] = []
    output_lines: list[str] = []

    for gene_id, gene_data in genes.items():
        gene_lines = gene_data.gene_lines
        transcript_ids = gene_data.transcript_ids
        gene_transcripts = OrderedDict(
            (transcript_id, transcripts[transcript_id])
            for transcript_id in transcript_ids
            if transcript_id in transcripts
        )

        if not gene_transcripts:
            output_lines.extend(gene_lines)
            continue

        if len(gene_transcripts) == 1:
            output_lines.extend(gene_lines)
            transcript_data = next(iter(gene_transcripts.values()))
            output_lines.extend(transcript_data.mrna_lines)
            output_lines.extend(transcript_data.exon_lines)
            output_lines.extend(transcript_data.cds_lines)
            output_lines.extend(transcript_data.five_prime_utr_lines)
            output_lines.extend(transcript_data.three_prime_utr_lines)
            continue

        best_tid: Optional[str] = None
        best_cds_len = -1
        best_utr_len = -1
        ranked_transcripts = sorted(
            gene_transcripts.items(),
            key=lambda item: item[1].sort_key
            if item[1].sort_key is not None
            else (10**18, 10**18, item[0]),
        )
        for transcript_id, transcript_data in ranked_transcripts:
            cds_len = transcript_data.total_cds
            utr_len = transcript_data.total_utr
            if cds_len > best_cds_len:
                best_tid = transcript_id
                best_cds_len = cds_len
                best_utr_len = utr_len
            elif cds_len == best_cds_len and utr_len > best_utr_len:
                best_tid = transcript_id
                best_utr_len = utr_len

        removed_here = len(gene_transcripts) - 1
        removed_mrna_count += removed_here
        if removed_here > 0:
            changed_genes_count += 1
            changed_gene_ids.append(gene_id)

        output_lines.extend(gene_lines)
        assert best_tid is not None
        best_transcript = gene_transcripts[best_tid]
        output_lines.extend(best_transcript.mrna_lines)
        output_lines.extend(best_transcript.exon_lines)
        output_lines.extend(best_transcript.cds_lines)
        output_lines.extend(best_transcript.five_prime_utr_lines)
        output_lines.extend(best_transcript.three_prime_utr_lines)

    Path(output_gff_path).write_text("\n".join(output_lines) + "\n", encoding="utf-8")
    if changed_gene_ids_path is not None:
        write_id_list(Path(changed_gene_ids_path), changed_gene_ids)
    return {
        "removed_mrnas": removed_mrna_count,
        "changed_genes": changed_genes_count,
        "changed_gene_ids": changed_gene_ids,
    }
