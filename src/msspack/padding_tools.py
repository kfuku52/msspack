from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional

from Bio.Data import CodonTable

from .fasta import iter_fasta, reverse_complement, write_fasta_record
from .gff import GFFRecord, child_ids, read_gff_document
from .step_logging import write_step_log, write_step_metrics
from .utils import ensure_dir

# Some of the padding logic below is carried forward from earlier internal
# tooling used before msspack unified the MSS packaging workflow.


@dataclass
class TranscriptModel:
    transcript_id: str
    seqid: str
    strand: str
    cds_records: List[GFFRecord]


@dataclass
class _PaddingChoice:
    new_seq: str
    num_stop: int
    headn: int
    tailn: int


def _resolve_codon_table(genetic_code: str | int):
    try:
        code = int(str(genetic_code))
    except ValueError:
        return CodonTable.unambiguous_dna_by_name[str(genetic_code)]
    return CodonTable.unambiguous_dna_by_id[code]


def _count_internal_stop_codons(seq: str, genetic_code: str | int) -> int:
    stop_codons = set(_resolve_codon_table(genetic_code).stop_codons)
    seq_upper = seq.upper()
    internal_stop_limit = len(seq_upper) - 3
    if internal_stop_limit <= 0:
        return 0
    count = 0
    for pos in range(0, internal_stop_limit, 3):
        if seq_upper[pos : pos + 3] in stop_codons:
            count += 1
    return count


def _get_padding_candidates(
    num_stop_input: int,
    num_missing: int,
    seqlen: int,
) -> List[tuple[int, int]]:
    candidates: List[tuple[int, int]] = []
    if num_stop_input:
        if num_missing in (0, 3):
            candidates.extend([(0, 0), (1, 2), (2, 1)])
        elif num_missing == 1:
            candidates.extend([(0, 1), (1, 0), (2, 2)])
        elif num_missing == 2:
            candidates.extend([(0, 2), (2, 0), (1, 1)])
    if (not num_stop_input) and (seqlen % 3):
        candidates.append((0, num_missing))
    return candidates


def _process_padding(
    record_name: str,
    record_seq: str,
    genetic_code: str | int,
    padchar: str = "N",
) -> dict[str, object]:
    clean_seq = record_seq.replace("X", "N")
    seqlen = len(clean_seq)
    adjlen = seqlen if (seqlen % 3 == 0) else ((seqlen // 3) + 1) * 3
    num_missing = adjlen - seqlen
    tailpad_seq = clean_seq.ljust(adjlen, padchar)
    num_stop_input = _count_internal_stop_codons(tailpad_seq, genetic_code)

    if not (num_stop_input or (seqlen % 3)):
        return {
            "new_seq": record_seq,
            "was_padded": False,
            "log": "",
        }

    best: Optional[_PaddingChoice] = None
    for headn, tailn in _get_padding_candidates(num_stop_input, num_missing, seqlen):
        if (headn == 0) and (tailn == num_missing):
            new_seq = tailpad_seq
            num_stop = num_stop_input
        else:
            new_seq = (padchar * headn) + clean_seq + (padchar * tailn)
            num_stop = _count_internal_stop_codons(new_seq, genetic_code)
        if (best is None) or (num_stop < best.num_stop):
            best = _PaddingChoice(
                new_seq=new_seq,
                num_stop=num_stop,
                headn=headn,
                tailn=tailn,
            )

    assert best is not None
    log = (
        f"{record_name}, original_seqlen={seqlen}, "
        f"head_padding={best.headn}, tail_padding={best.tailn}, "
        f"original_num_stop={num_stop_input}, new_num_stop={best.num_stop}\n"
    )
    return {
        "new_seq": best.new_seq,
        "was_padded": bool(best.headn or best.tailn),
        "log": log,
    }


def _build_transcript_models(gff_path: Path) -> List[TranscriptModel]:
    transcript_records: Dict[str, GFFRecord] = {}
    transcript_order: List[str] = []
    cds_by_parent: Dict[str, List[GFFRecord]] = {}

    for record in read_gff_document(gff_path).records:
        record_id = record.attributes.get("ID", "")
        if record.type in ("mRNA", "transcript") and record_id:
            transcript_records[record_id] = record
            transcript_order.append(record_id)
        elif record.type == "CDS":
            for parent_id in child_ids(record.attributes.get("Parent")):
                cds_by_parent.setdefault(parent_id, []).append(record)
                if parent_id not in transcript_records and parent_id not in transcript_order:
                    transcript_order.append(parent_id)

    models: List[TranscriptModel] = []
    for transcript_id in transcript_order:
        cds_records = cds_by_parent.get(transcript_id, [])
        if not cds_records:
            continue
        transcript_record = transcript_records.get(transcript_id)
        seqid = transcript_record.seqid if transcript_record else cds_records[0].seqid
        strand = transcript_record.strand if transcript_record else cds_records[0].strand
        models.append(
            TranscriptModel(
                transcript_id=transcript_id,
                seqid=seqid,
                strand=strand,
                cds_records=sorted(cds_records, key=lambda record: (record.start, record.end)),
            )
        )
    return models


def write_spliced_cds_fasta(
    *,
    fasta_path: Path,
    gff_path: Path,
    output_path: Path,
    log_path: Path,
    metrics_path: Path | None = None,
) -> None:
    command = [
        "msspack",
        "internal",
        "extract-cds",
        "--fasta",
        str(fasta_path),
        "--gff",
        str(gff_path),
        "--out",
        str(output_path),
    ]
    models = _build_transcript_models(gff_path)
    models_by_seqid: Dict[str, List[TranscriptModel]] = defaultdict(list)
    for model in models:
        models_by_seqid[model.seqid].append(model)
    written = 0
    started_at = datetime.now()

    ensure_dir(output_path.parent)
    with output_path.open("w", encoding="utf-8") as handle:
        for record in iter_fasta(fasta_path):
            seq_models = models_by_seqid.get(record.id)
            if not seq_models:
                continue
            contig_seq = record.sequence
            for model in seq_models:
                pieces = [contig_seq[feature.start - 1 : feature.end] for feature in model.cds_records]
                seq = "".join(pieces)
                if model.strand == "-":
                    seq = reverse_complement(seq)
                write_fasta_record(handle, identifier=model.transcript_id, sequence=seq)
                written += 1

    write_step_log(
        log_path=log_path,
        command=command,
        step="extract-cds",
        started_at=started_at,
        count_unit="transcripts",
        input_total=written,
        changed_total=written,
        output_total=written,
        details=[f"Number of transcript CDS records written: {written:,}"],
    )
    if metrics_path is not None:
        write_step_metrics(
            metrics_path=metrics_path,
            step="extract-cds",
            count_unit="transcripts",
            input_total=written,
            changed_total=written,
            output_total=written,
            details={"transcript_cds_records_written": written},
        )


def pad_fasta(
    *,
    input_path: Path,
    output_path: Path,
    log_path: Path,
    genetic_code: str | int = "1",
    padchar: str = "N",
    metrics_path: Path | None = None,
) -> None:
    command = [
        "msspack",
        "internal",
        "pad",
        "-s",
        str(input_path),
        "-o",
        str(output_path),
    ]
    started_at = datetime.now()
    padded = 0
    record_count = 0
    detail_lines: list[str] = []

    ensure_dir(output_path.parent)
    with output_path.open("w", encoding="utf-8") as out_handle:
        for record in iter_fasta(input_path):
            result = _process_padding(record.id, record.sequence, genetic_code, padchar)
            write_fasta_record(
                out_handle,
                identifier=record.id,
                description=record.description,
                sequence=str(result["new_seq"]),
            )
            log_text = str(result["log"])
            if log_text:
                detail_lines.append(log_text.rstrip("\n"))
            if result["was_padded"]:
                padded += 1
            record_count += 1

    detail_lines.append(f"Number of padded sequences: {padded:,} / {record_count:,}")
    write_step_log(
        log_path=log_path,
        command=command,
        step="pad",
        started_at=started_at,
        count_unit="sequences",
        input_total=record_count,
        changed_total=padded,
        output_total=record_count,
        details=detail_lines,
    )
    if metrics_path is not None:
        write_step_metrics(
            metrics_path=metrics_path,
            step="pad",
            count_unit="sequences",
            input_total=record_count,
            changed_total=padded,
            output_total=record_count,
            details={"padded_sequences": padded},
        )


def write_padding_log_for_gff(
    *,
    fasta_path: Path,
    gff_path: Path,
    log_path: Path,
    genetic_code: str | int = "1",
    padchar: str = "N",
    extract_log_path: Optional[Path] = None,
    extract_metrics_path: Optional[Path] = None,
    metrics_path: Optional[Path] = None,
) -> None:
    extract_command = [
        "msspack",
        "internal",
        "extract-cds",
        "--fasta",
        str(fasta_path),
        "--gff",
        str(gff_path),
    ]
    pad_command = [
        "msspack",
        "internal",
        "pad-from-gff",
        "--fasta",
        str(fasta_path),
        "--gff",
        str(gff_path),
    ]
    started_at = datetime.now()
    models = _build_transcript_models(gff_path)
    models_by_seqid: Dict[str, List[TranscriptModel]] = defaultdict(list)
    for model in models:
        models_by_seqid[model.seqid].append(model)

    transcript_count = 0
    padded = 0
    log_payload: List[str] = []

    for record in iter_fasta(fasta_path):
        seq_models = models_by_seqid.get(record.id)
        if not seq_models:
            continue
        for model in seq_models:
            pieces = [record.sequence[feature.start - 1 : feature.end] for feature in model.cds_records]
            seq = "".join(pieces)
            if model.strand == "-":
                seq = reverse_complement(seq)
            result = _process_padding(model.transcript_id, seq, genetic_code, padchar)
            transcript_count += 1
            if result["was_padded"]:
                padded += 1
            if result["log"]:
                log_payload.append(str(result["log"]).rstrip("\n"))

    detail_lines = list(log_payload)
    detail_lines.append(f"Number of padded sequences: {padded:,} / {transcript_count:,}")
    write_step_log(
        log_path=log_path,
        command=pad_command,
        step="pad-from-gff",
        started_at=started_at,
        count_unit="transcripts",
        input_total=transcript_count,
        changed_total=padded,
        output_total=transcript_count,
        details=detail_lines,
    )
    if metrics_path is not None:
        write_step_metrics(
            metrics_path=metrics_path,
            step="pad-from-gff",
            count_unit="transcripts",
            input_total=transcript_count,
            changed_total=padded,
            output_total=transcript_count,
            details={"padded_sequences": padded},
        )

    if extract_log_path is not None:
        write_step_log(
            log_path=extract_log_path,
            command=extract_command,
            step="extract-cds",
            started_at=started_at,
            count_unit="transcripts",
            input_total=transcript_count,
            changed_total=transcript_count,
            output_total=transcript_count,
            details=[f"Number of transcript CDS records written: {transcript_count:,}"],
        )
    if extract_metrics_path is not None:
        write_step_metrics(
            metrics_path=extract_metrics_path,
            step="extract-cds",
            count_unit="transcripts",
            input_total=transcript_count,
            changed_total=transcript_count,
            output_total=transcript_count,
            details={"transcript_cds_records_written": transcript_count},
        )
