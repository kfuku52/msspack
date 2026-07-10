from __future__ import annotations

import re
from collections import defaultdict, deque
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Sequence
from urllib.parse import unquote_to_bytes

from .utils import MSSPackError, atomic_text_writer, ensure_dir

_INVALID_PERCENT_ESCAPE = re.compile(r"%(?![0-9A-Fa-f]{2})")


def _decode_attribute_component(value: str) -> str:
    if _INVALID_PERCENT_ESCAPE.search(value):
        raise ValueError(f"Invalid percent escape in GFF3 attribute: {value!r}")
    try:
        return unquote_to_bytes(value).decode("utf-8")
    except UnicodeDecodeError as exc:
        raise ValueError(f"Invalid UTF-8 escape in GFF3 attribute: {value!r}") from exc


def parse_attributes(text: str) -> Dict[str, str]:
    if text.strip() in ("", "."):
        return {}
    attrs: Dict[str, str] = {}
    for chunk in text.split(";"):
        chunk = chunk.strip()
        if not chunk:
            continue
        if "=" not in chunk:
            raise ValueError(f"Invalid GFF3 attribute without '=': {chunk!r}")
        key, value = chunk.split("=", 1)
        key = _decode_attribute_component(key.strip())
        if not key:
            raise ValueError("Invalid empty GFF3 attribute key")
        if key in attrs:
            raise ValueError(f"Duplicate GFF3 attribute key: {key}")
        decoded_value = _decode_attribute_component(value.strip())
        if any(character in decoded_value for character in ("\x00", "\r", "\n", "\t")):
            raise ValueError(f"GFF3 attribute {key!r} contains unsupported control characters")
        attrs[key] = decoded_value
    return attrs


@dataclass
class GFFRecord:
    seqid: str
    source: str
    type: str
    start: int
    end: int
    score: str
    strand: str
    phase: str
    attributes_text: str
    attributes: Dict[str, str]

    @classmethod
    def from_line(cls, line: str) -> "GFFRecord":
        fields = line.rstrip("\n").split("\t")
        if len(fields) != 9:
            raise ValueError(f"Expected 9 GFF columns, found {len(fields)}: {line!r}")
        return cls(
            seqid=fields[0],
            source=fields[1],
            type=fields[2],
            start=int(fields[3]),
            end=int(fields[4]),
            score=fields[5],
            strand=fields[6],
            phase=fields[7],
            attributes_text=fields[8],
            attributes=parse_attributes(fields[8]),
        )

    def to_line(self) -> str:
        return "\t".join(
            [
                self.seqid,
                self.source,
                self.type,
                str(self.start),
                str(self.end),
                self.score,
                self.strand,
                self.phase,
                self.attributes_text,
            ]
        )


@dataclass
class GFFDocument:
    header_lines: List[str]
    records: List[GFFRecord]
    fasta_lines: List[str] = field(default_factory=list)


@dataclass(frozen=True)
class _SortEntry:
    line: str
    record_id: str
    parents: tuple[str, ...]


def read_gff_document(path: str | Path) -> GFFDocument:
    header_lines: List[str] = []
    records: List[GFFRecord] = []
    fasta_lines: List[str] = []
    in_fasta = False
    with Path(path).open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if in_fasta:
                fasta_lines.append(line)
                continue
            if line == "##FASTA":
                in_fasta = True
                fasta_lines.append(line)
                continue
            if not line:
                header_lines.append(line)
                continue
            if line.startswith("#"):
                header_lines.append(line)
                continue
            records.append(GFFRecord.from_line(line))
    return GFFDocument(
        header_lines=header_lines,
        records=records,
        fasta_lines=fasta_lines,
    )


def write_gff_document(path: str | Path, document: GFFDocument) -> Path:
    out_path = Path(path)
    ensure_dir(out_path.parent)
    with atomic_text_writer(out_path) as handle:
        for line in document.header_lines:
            handle.write(line + "\n")
        for record in document.records:
            handle.write(record.to_line() + "\n")
        for line in document.fasta_lines:
            handle.write(line + "\n")
    return out_path


def iter_gff_records(path: str | Path) -> Iterator[GFFRecord]:
    with Path(path).open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.rstrip("\n") == "##FASTA":
                break
            if (not line.strip()) or line.startswith("#"):
                continue
            yield GFFRecord.from_line(line)


def child_ids(value: str | None) -> List[str]:
    if not value:
        return []
    return [item for item in value.split(",") if item]


def attribute(record: GFFRecord, key: str, default: str = "") -> str:
    return record.attributes.get(key, default)


def iter_records_by_type(
    records: Iterable[GFFRecord],
    feature_type: str,
) -> Iterator[GFFRecord]:
    for record in records:
        if record.type == feature_type:
            yield record


def _sort_same_start_block(entries: Sequence[_SortEntry]) -> list[str]:
    id_to_indices: Dict[str, list[int]] = defaultdict(list)
    for index, entry in enumerate(entries):
        if entry.record_id:
            id_to_indices[entry.record_id].append(index)

    children: Dict[int, list[int]] = defaultdict(list)
    indegree = [0] * len(entries)
    seen_edges: set[tuple[int, int]] = set()
    for child_index, entry in enumerate(entries):
        for parent_id in entry.parents:
            for parent_index in id_to_indices.get(parent_id, []):
                if parent_index == child_index:
                    continue
                edge = (parent_index, child_index)
                if edge in seen_edges:
                    continue
                seen_edges.add(edge)
                children[parent_index].append(child_index)
                indegree[child_index] += 1

    ready = deque(index for index, degree in enumerate(indegree) if degree == 0)
    ordered: list[int] = []
    while ready:
        index = ready.popleft()
        ordered.append(index)
        for child_index in children.get(index, []):
            indegree[child_index] -= 1
            if indegree[child_index] == 0:
                ready.append(child_index)

    if len(ordered) != len(entries):
        emitted = set(ordered)
        ordered.extend(index for index in range(len(entries)) if index not in emitted)

    return [entries[index].line for index in ordered]


def sort_gff_file_precise(
    input_path: str | Path,
    output_path: str | Path,
) -> Path:
    source = Path(input_path)
    destination = Path(output_path)
    ensure_dir(destination.parent)

    pragma_lines: list[str] = []
    fasta_lines: list[str] = []
    chromosome_order: list[str] = []
    blocks: Dict[str, Dict[int, list[_SortEntry]]] = defaultdict(lambda: defaultdict(list))
    with source.open("r", encoding="utf-8") as handle:
        in_fasta = False
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\n")
            if in_fasta:
                fasta_lines.append(line)
                continue
            if line == "##FASTA":
                in_fasta = True
                fasta_lines.append(line)
                continue
            if not line:
                continue
            if line.startswith("#"):
                pragma_lines.append(line)
                continue

            fields = line.split("\t")
            if len(fields) != 9:
                raise MSSPackError(
                    f"Invalid GFF record at {source}:{line_number}: expected 9 columns, "
                    f"found {len(fields)}"
                )
            seqid = fields[0]
            try:
                start = int(fields[3])
                attrs = parse_attributes(fields[8])
            except ValueError as exc:
                raise MSSPackError(f"Invalid GFF record at {source}:{line_number}: {exc}") from exc
            entry = _SortEntry(
                line=line,
                record_id=attrs.get("ID", ""),
                parents=tuple(child_ids(attrs.get("Parent"))),
            )
            blocks[seqid][start].append(entry)
            if seqid not in chromosome_order:
                chromosome_order.append(seqid)

    chromosome_order = sorted(chromosome_order)

    with atomic_text_writer(destination) as out_handle:
        for line in pragma_lines:
            out_handle.write(line + "\n")
        for seqid in chromosome_order:
            for start in sorted(blocks[seqid]):
                entries = blocks[seqid][start]
                if len(entries) == 1:
                    out_handle.write(entries[0].line + "\n")
                    continue
                for line in _sort_same_start_block(entries):
                    out_handle.write(line + "\n")
        for line in fasta_lines:
            out_handle.write(line + "\n")

    return destination
