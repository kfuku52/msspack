"""Order-independent CDS models shared by extraction and functional annotation."""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

from .fasta import reverse_complement
from .gff import GFFRecord, child_ids, iter_gff_records
from .utils import MSSPackError


@dataclass(frozen=True)
class TranscriptModel:
    transcript_id: str
    seqid: str
    strand: str
    cds_records: tuple[GFFRecord, ...]


def build_transcript_models(gff_path: Path) -> list[TranscriptModel]:
    parents: dict[str, GFFRecord] = {}
    order: dict[str, None] = {}
    cds_by_parent: dict[str, list[GFFRecord]] = defaultdict(list)
    for record in iter_gff_records(gff_path):
        record_id = record.attributes.get("ID", "")
        if record.type in ("gene", "mRNA", "transcript") and record_id:
            parents[record_id] = record
            order.setdefault(record_id, None)
        elif record.type == "CDS":
            for parent_id in child_ids(record.attributes.get("Parent")):
                order.setdefault(parent_id, None)
                cds_by_parent[parent_id].append(record)
    models = []
    for identifier in order:
        cdss = cds_by_parent.get(identifier, [])
        if not cdss:
            continue
        parent = parents.get(identifier, cdss[0])
        if any(cds.seqid != parent.seqid or cds.strand != parent.strand for cds in cdss):
            raise MSSPackError(f"Inconsistent CDS sequence/strand for {identifier}")
        models.append(TranscriptModel(
            identifier, parent.seqid, parent.strand,
            tuple(sorted(cdss, key=lambda row: (row.start, row.end))),
        ))
    return models


def spliced_cds_sequence(sequence: str, model: TranscriptModel) -> str:
    """Apply only the first CDS phase in transcription order, after splicing."""
    for cds in model.cds_records:
        if not 1 <= cds.start <= cds.end <= len(sequence):
            raise MSSPackError(f"CDS outside FASTA bounds: {model.transcript_id}")
        if cds.phase not in {"0", "1", "2", "."}:
            raise MSSPackError(f"Invalid CDS phase for {model.transcript_id}: {cds.phase}")
    spliced = "".join(sequence[cds.start - 1:cds.end] for cds in model.cds_records)
    first = model.cds_records[0]
    if model.strand == "-":
        spliced = reverse_complement(spliced)
        first = model.cds_records[-1]
    return spliced[int(first.phase) if first.phase != "." else 0:]
