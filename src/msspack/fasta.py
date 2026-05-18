from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, TextIO

_COMPLEMENT_TABLE = str.maketrans(
    "ACGTRYMKBDHVSWNUacgtrymkbdhvswnu",
    "TGCAYRKMVHDBSWNAtgcayrkmvhdbswna",
)


@dataclass(frozen=True)
class FastaRecord:
    id: str
    description: str
    sequence: str


def _build_record(header: str, sequence_lines: list[str]) -> FastaRecord:
    description = header.strip()
    if not description:
        raise ValueError("Encountered FASTA header without an identifier")
    seq_id = description.split(None, 1)[0]
    return FastaRecord(
        id=seq_id,
        description=description,
        sequence="".join(sequence_lines),
    )


def iter_fasta_handle(handle: TextIO) -> Iterator[FastaRecord]:
    header: str | None = None
    sequence_lines: list[str] = []
    for raw_line in handle:
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                yield _build_record(header, sequence_lines)
            header = line[1:]
            sequence_lines = []
            continue
        if header is None:
            raise ValueError("Encountered FASTA sequence before the first header")
        sequence_lines.append(line)
    if header is not None:
        yield _build_record(header, sequence_lines)


def iter_fasta(path: str | Path) -> Iterator[FastaRecord]:
    with Path(path).open("r", encoding="utf-8") as handle:
        yield from iter_fasta_handle(handle)


def write_wrapped_sequence(
    handle: TextIO,
    sequence: str,
    *,
    width: int = 60,
) -> None:
    if width <= 0:
        handle.write(sequence + "\n")
        return
    if not sequence:
        handle.write("\n")
        return
    for start in range(0, len(sequence), width):
        handle.write(sequence[start : start + width] + "\n")


def write_fasta_record(
    handle: TextIO,
    *,
    identifier: str,
    sequence: str,
    description: str | None = None,
    width: int = 60,
) -> None:
    header = description or identifier
    handle.write(f">{header}\n")
    write_wrapped_sequence(handle, sequence, width=width)


def write_fasta_records(
    path: str | Path,
    records: Iterable[FastaRecord],
    *,
    width: int = 60,
) -> None:
    with Path(path).open("w", encoding="utf-8") as handle:
        for record in records:
            write_fasta_record(
                handle,
                identifier=record.id,
                description=record.description,
                sequence=record.sequence,
                width=width,
            )


def reverse_complement(sequence: str) -> str:
    return sequence.translate(_COMPLEMENT_TABLE)[::-1]
