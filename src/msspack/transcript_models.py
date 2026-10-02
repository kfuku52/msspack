"""Order-independent CDS models shared by extraction and functional annotation."""
from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import TypeVar

from .fasta import reverse_complement
from .gff import GFFRecord, child_ids, iter_gff_records
from .utils import MSSPackError


@dataclass(frozen=True)
class TranscriptModel:
    transcript_id: str
    seqid: str
    strand: str
    cds_records: tuple[GFFRecord, ...]
    parent_id: str = ""


_CDS = TypeVar("_CDS")


def group_direct_cds(parent_id: str, records: Iterable[tuple[str, _CDS]]) -> dict[str, list[_CDS]]:
    """Repeated CDS IDs join; distinct IDs under a gene are separate coding units.

    Retain the parent identifier for a single coding unit, including anonymous
    CDS segments, so existing extraction and annotation identifiers stay valid.
    """
    groups: dict[str, list[_CDS]] = {}
    for identifier, record in records:
        groups.setdefault(identifier or parent_id, []).append(record)
    if len(groups) == 1:
        return {parent_id: next(iter(groups.values()))}
    return groups


def gene_model_ids(records: Iterable[GFFRecord]) -> tuple[str, ...]:
    """Identify explicit genes and virtual containers for root transcripts.

    Keep one model per Parent gene, or per transcript ID when Parent is absent.
    This is also the unit used by CDS adjustment and pipeline plot metrics.
    """
    identifiers: dict[str, None] = {}
    for record in records:
        identifier = record.attributes.get("ID", "")
        if not identifier:
            continue
        if record.type == "gene":
            identifiers.setdefault(identifier, None)
        elif record.type in {"mRNA", "transcript"}:
            for gene_id in child_ids(record.attributes.get("Parent")) or [identifier]:
                identifiers.setdefault(gene_id, None)
    return tuple(identifiers)


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
        groups = (
            {identifier: cdss} if parent.type in {"mRNA", "transcript"}
            else group_direct_cds(identifier, (
                (cds.attributes.get("ID", ""), cds) for cds in cdss
            ))
        )
        for model_id, segments in groups.items():
            if any(cds.seqid != parent.seqid or cds.strand != parent.strand for cds in segments):
                raise MSSPackError(f"Inconsistent CDS sequence/strand for {model_id}")
            models.append(TranscriptModel(
                model_id, parent.seqid, parent.strand,
                tuple(sorted(segments, key=lambda row: (row.start, row.end))), identifier,
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
