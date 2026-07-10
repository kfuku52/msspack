from __future__ import annotations

import re
from collections import OrderedDict
from dataclasses import dataclass, field
from pathlib import Path
from typing import TypedDict
from urllib.parse import unquote

from .gff import child_ids, parse_attributes
from .step_logging import write_id_list
from .utils import ensure_dir, write_text


@dataclass
class _TranscriptData:
    gene_ids: list[str] = field(default_factory=list)
    total_cds: int = 0
    total_utr: int = 0
    sort_key: tuple[int, int, str] | None = None


@dataclass
class _GeneData:
    transcript_ids: list[str] = field(default_factory=list)


@dataclass(frozen=True)
class _FeatureLine:
    raw: str
    feature_id: str
    parents: tuple[str, ...]


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


def _raw_parent_tokens(attributes_text: str) -> list[tuple[str, str]]:
    for chunk in attributes_text.split(";"):
        if "=" not in chunk:
            continue
        key, value = chunk.split("=", 1)
        if unquote(key.strip()) != "Parent":
            continue
        return [(unquote(token), token) for token in value.split(",") if token]
    return []


def _filter_parent_attribute(line: str, removed_ids: set[str]) -> str | None:
    fields = line.split("\t")
    if len(fields) != 9:
        return line
    tokens = _raw_parent_tokens(fields[8])
    if not tokens:
        return line
    kept_raw = [raw for decoded, raw in tokens if decoded not in removed_ids]
    if not kept_raw:
        return None
    if len(kept_raw) == len(tokens):
        return line
    fields[8] = re.sub(
        r"(?:(?<=;)|^)Parent=[^;]*",
        "Parent=" + ",".join(kept_raw),
        fields[8],
        count=1,
    )
    return "\t".join(fields)


def _selected_transcript_ids(
    genes: OrderedDict[str, _GeneData],
    transcripts: OrderedDict[str, _TranscriptData],
) -> tuple[set[str], set[str], dict[str, set[str]], list[str]]:
    selected: set[str] = set()
    candidates_seen: set[str] = set()
    selected_gene_ids: dict[str, set[str]] = {}
    changed_gene_ids: list[str] = []
    for gene_id, gene_data in genes.items():
        candidates = [
            (transcript_id, transcripts[transcript_id])
            for transcript_id in gene_data.transcript_ids
            if transcript_id in transcripts
        ]
        if not candidates:
            continue
        candidates_seen.update(transcript_id for transcript_id, _ in candidates)
        ranked = sorted(
            candidates,
            key=lambda item: (
                -item[1].total_cds,
                -item[1].total_utr,
                item[1].sort_key
                if item[1].sort_key is not None
                else (10**18, 10**18, item[0]),
            ),
        )
        selected_transcript_id = ranked[0][0]
        selected.add(selected_transcript_id)
        selected_gene_ids.setdefault(selected_transcript_id, set()).add(gene_id)
        if len(ranked) > 1:
            changed_gene_ids.append(gene_id)
    return selected, candidates_seen, selected_gene_ids, changed_gene_ids


def select_one_mrna_per_gene(
    *,
    input_gff_path: str | Path,
    output_gff_path: str | Path,
    changed_gene_ids_path: str | Path | None = None,
) -> TranscriptSelectionSummary:
    genes: OrderedDict[str, _GeneData] = OrderedDict()
    transcripts: OrderedDict[str, _TranscriptData] = OrderedDict()
    declared_transcript_ids: set[str] = set()
    items: list[str | _FeatureLine] = []

    with Path(input_gff_path).open("r", encoding="utf-8") as infile:
        in_fasta = False
        for raw_line in infile:
            line = raw_line.rstrip("\n")
            if in_fasta:
                items.append(line)
                continue
            if line == "##FASTA":
                in_fasta = True
                items.append(line)
                continue
            if not line or line.startswith("#"):
                items.append(line)
                continue
            cols = line.split("\t")
            if len(cols) != 9:
                items.append(line)
                continue
            _, _, feature_type, start, end, _, _, _, attributes_text = cols
            attributes = parse_attributes(attributes_text)
            feature_id = attributes.get("ID", "")
            parents = tuple(child_ids(attributes.get("Parent")))
            items.append(_FeatureLine(line, feature_id, parents))

            if feature_type == "gene" and feature_id:
                _ensure_gene(genes, feature_id)
                continue
            if feature_type in ("mRNA", "transcript") and feature_id:
                declared_transcript_ids.add(feature_id)
                transcript = _ensure_transcript(transcripts, feature_id)
                for gene_id in parents:
                    gene = _ensure_gene(genes, gene_id)
                    if feature_id not in gene.transcript_ids:
                        gene.transcript_ids.append(feature_id)
                    if gene_id not in transcript.gene_ids:
                        transcript.gene_ids.append(gene_id)
                candidate_key = (int(start), int(end), feature_id)
                if transcript.sort_key is None or candidate_key < transcript.sort_key:
                    transcript.sort_key = candidate_key
                continue
            if feature_type in ("CDS", "five_prime_UTR", "three_prime_UTR"):
                for transcript_id in parents:
                    transcript = _ensure_transcript(transcripts, transcript_id)
                    if transcript.sort_key is None:
                        transcript.sort_key = (int(start), int(end), transcript_id)
                    if feature_type == "CDS":
                        transcript.total_cds += _feature_length(start, end)
                    else:
                        transcript.total_utr += _feature_length(start, end)

    selected, candidates_seen, selected_gene_ids, changed_gene_ids = _selected_transcript_ids(
        genes,
        transcripts,
    )
    # Only declared transcript features that actually participated in a gene's
    # candidate set may be removed.  Parent identifiers inferred from CDS/UTR
    # rows can also be gene identifiers and must never be treated as mRNAs.
    removed_transcript_ids = (declared_transcript_ids & candidates_seen) - selected
    removed_ids = set(removed_transcript_ids)

    # Remove descendants only when every parent was removed. Shared features keep
    # their surviving Parent references.
    changed = True
    while changed:
        changed = False
        for item in items:
            if not isinstance(item, _FeatureLine) or not item.feature_id or not item.parents:
                continue
            if item.feature_id in removed_ids:
                continue
            if all(parent in removed_ids for parent in item.parents):
                removed_ids.add(item.feature_id)
                changed = True

    output_lines: list[str] = []
    for item in items:
        if isinstance(item, str):
            output_lines.append(item)
            continue
        if item.feature_id and item.feature_id in removed_ids:
            continue
        if item.feature_id in selected_gene_ids:
            disallowed_gene_ids = set(item.parents) - selected_gene_ids[item.feature_id]
            filtered = _filter_parent_attribute(item.raw, disallowed_gene_ids)
        else:
            filtered = _filter_parent_attribute(item.raw, removed_ids)
        if filtered is not None:
            output_lines.append(filtered)

    output_path = Path(output_gff_path)
    ensure_dir(output_path.parent)
    write_text(output_path, "\n".join(output_lines) + "\n")
    if changed_gene_ids_path is not None:
        write_id_list(Path(changed_gene_ids_path), changed_gene_ids)
    return {
        "removed_mrnas": len(removed_transcript_ids),
        "changed_genes": len(changed_gene_ids),
        "changed_gene_ids": changed_gene_ids,
    }
