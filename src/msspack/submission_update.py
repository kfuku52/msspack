"""Prepare annotation-only updates without rerunning sequence/model cleanup."""
from __future__ import annotations

import csv
import hashlib
import json
import re
import tempfile
from collections import Counter
from pathlib import Path

from Bio.SeqFeature import ExactPosition, Location

from .utils import MSSPackError


def _read_mapping(path: Path) -> dict[str, tuple[str, str]]:
    mapping: dict[str, tuple[str, str]] = {}
    accessions: set[str] = set()
    seqids: set[str] = set()
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        if next(reader, None) != ["entry", "accession", "submitter_seqid"]:
            raise MSSPackError("Mapping header must be: entry<TAB>accession<TAB>submitter_seqid")
        for row in reader:
            if len(row) != 3 or any(not x or any(c.isspace() for c in x) for x in row):
                raise MSSPackError(f"Invalid mapping row at line {reader.line_num}")
            entry, accession, seqid = row
            if entry == "COMMON" or not re.fullmatch(r"[A-Z]{1,6}[0-9]{5,12}", accession):
                raise MSSPackError(f"Invalid entry/accession in mapping: {entry}, {accession}")
            if entry in mapping or accession in accessions or seqid in seqids:
                raise MSSPackError("Mapping entries, accessions and submitter_seqids must be unique")
            mapping[entry] = accession, seqid
            accessions.add(accession)
            seqids.add(seqid)
    if not mapping:
        raise MSSPackError("Mapping is empty")
    return mapping


def _qualifier(block: list[list[str]], name: str) -> str:
    return next((r[4] for r in block if r[3] == name), "")


def _covered_by(feature: list[list[str]], mrna: list[list[str]]) -> bool:
    """Only drop locations whose exact, same-strand coverage is represented."""
    try:
        child = Location.fromstring(feature[0][2])
        parent = Location.fromstring(mrna[0][2])
        if child is None or parent is None:
            return False
        if any(type(p.start) is not ExactPosition or type(p.end) is not ExactPosition
               for p in [*child.parts, *parent.parts]):
            return False
        return all(
            not part.ref and any(
                not outer.ref and part.strand == outer.strand
                and int(outer.start) <= int(part.start) < int(part.end) <= int(outer.end)
                for outer in parent.parts
            )
            for part in child.parts
        )
    except (ValueError, TypeError):
        return False


def _write_entry(
    blocks: list[list[list[str]]], mapping: dict[str, tuple[str, str]],
    retain_utr: bool, counts: Counter[str],
) -> str:
    entry = blocks[0][0][0]
    if entry == "COMMON":
        # An entry macro here would change meaning after renaming; require explicit sources.
        if any(r[3] == "submitter_seqid" for b in blocks for r in b):
            raise MSSPackError("Move COMMON submitter_seqid qualifiers to individual source features")
        if any("@@[entry]@@" in r[4] for b in blocks for r in b):
            raise MSSPackError(
                "Expand COMMON entry macros into individual entries before preparing an update"
            )
        return "".join("\t".join(r) + "\n" for b in blocks for r in b)
    if entry not in mapping:
        raise MSSPackError(f"Annotation entry missing from mapping: {entry}")
    accession, seqid = mapping[entry]
    sources = [b for b in blocks if b[0][1] == "source"]
    if len(sources) != 1:
        raise MSSPackError(f"Expected one source feature for {entry}")
    source = sources[0]
    seqid_rows = [r for r in source if r[3] == "submitter_seqid"]
    if len(seqid_rows) > 1:
        raise MSSPackError(f"Duplicate submitter_seqid for {entry}")
    if seqid_rows:
        old = seqid_rows[0][4]
        if old not in {"@@[entry]@@", seqid}:
            raise MSSPackError(f"Mapping conflicts with existing submitter_seqid for {entry}: {old}")
        seqid_rows[0][4] = seqid
    else:
        source.append(["", "", "", "submitter_seqid", seqid])
    mrnas: dict[str, list[list[list[str]]]] = {}
    for b in blocks:
        locus = _qualifier(b, "locus_tag")
        if b[0][1] == "mRNA" and locus:
            mrnas.setdefault(locus, []).append(b)
    result: list[str] = []
    for b in blocks:
        key = b[0][1]
        candidates = mrnas.get(_qualifier(b, "locus_tag"), [])
        removable = key == "exon" or (not retain_utr and key in {"5'UTR", "3'UTR"})
        if removable and any(_covered_by(b, m) for m in candidates):
            counts[f"omitted_{key}"] += 1
            continue
        if removable:
            counts[f"retained_unrepresented_{key}"] += 1
        for row in b:
            if row[0]:
                row[0] = accession
            # Resolve entry macros before renaming so notes/definitions keep their meaning.
            row[4] = row[4].replace("@@[entry]@@", entry)
        counts[f"output_{key}"] += 1
        result.extend("\t".join(r) + "\n" for r in b)
    return "".join(result)


def _rewrite_annotation(
    source: Path, target: Path, mapping: dict[str, tuple[str, str]], retain_utr: bool,
) -> tuple[Counter[str], dict[str, int]]:
    counts: Counter[str] = Counter()
    seen: set[str] = set()
    lengths: dict[str, int] = {}
    blocks: list[list[list[str]]] = []
    current = ""
    with source.open(encoding="utf-8") as src, target.open("w", encoding="utf-8") as dst:
        for number, line in enumerate(src, 1):
            if not line.strip():
                continue
            row = line.rstrip("\r\n").split("\t")
            if not 3 <= len(row) <= 5:
                raise MSSPackError(f"Expected three to five annotation columns at line {number}")
            row += [""] * (5 - len(row))
            if row[0] and row[0] != current:
                if blocks:
                    dst.write(_write_entry(blocks, mapping, retain_utr, counts))
                current = row[0]
                if current in seen:
                    raise MSSPackError(f"Repeated annotation entry: {current}")
                seen.add(current)
                blocks = []
            if row[1] == "source" and current != "COMMON":
                match = re.fullmatch(r"1\.\.([1-9][0-9]*)", row[2])
                if not match or current in lengths:
                    raise MSSPackError(f"Expected a single full-length source for {current}")
                lengths[current] = int(match.group(1))
            if row[1]:
                blocks.append([row])
            elif blocks:
                blocks[-1].append(row)
            else:
                raise MSSPackError(f"Qualifier without a feature at line {number}")
        if blocks:
            dst.write(_write_entry(blocks, mapping, retain_utr, counts))
    seen.discard("COMMON")
    if seen != set(mapping):
        raise MSSPackError("Annotation entries and mapping entries do not match exactly")
    return counts, lengths


def _rewrite_fasta(
    source: Path, target: Path, mapping: dict[str, tuple[str, str]],
) -> dict[str, dict[str, str | int]]:
    sequences: dict[str, dict[str, str | int]] = {}
    current = ""
    digest = hashlib.sha256()
    length = 0
    ended = False
    with source.open(encoding="utf-8") as src, target.open("w", encoding="utf-8") as dst:
        for line in src:
            if line.startswith(">"):
                if current:
                    sequences[current] = {"length": length, "sha256": digest.hexdigest()}
                fields = line[1:].split()
                entry = fields[0] if fields else ""
                if entry not in mapping or entry in sequences:
                    raise MSSPackError(f"Unknown or duplicate FASTA entry: {entry}")
                current = entry
                digest = hashlib.sha256()
                length = 0
                ended = False
                dst.write(f">{mapping[entry][0]}\n")
            else:
                sequence = line.strip()
                if sequence and not current:
                    raise MSSPackError("FASTA sequence before header")
                if sequence == "//":
                    if ended:
                        raise MSSPackError(f"Repeated FASTA terminator: {current}")
                    ended = True
                elif sequence:
                    if ended or not re.fullmatch(r"[ACGTRYSWKMBDHVNacgtryswkmbdhvn]+", sequence):
                        raise MSSPackError(f"Invalid DNA sequence for {current}")
                    digest.update(sequence.encode("ascii"))
                    length += len(sequence)
                dst.write(line)
        if current:
            sequences[current] = {"length": length, "sha256": digest.hexdigest()}
    if set(sequences) != set(mapping) or any(s["length"] == 0 for s in sequences.values()):
        raise MSSPackError("Nonempty FASTA entries must match the mapping exactly")
    return sequences


def prepare_update(
    *, ann_path: Path, fasta_path: Path, mapping_path: Path, output_dir: Path,
    retain_utr_features: bool = False,
) -> Path:
    """Write a new directory atomically; preserve bases, CDS and protein identifiers."""
    if output_dir.exists():
        raise MSSPackError(f"Output directory already exists: {output_dir}")
    mapping = _read_mapping(mapping_path)
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".msspack-update-", dir=output_dir.parent) as tmp:
        stage = Path(tmp) / "result"
        stage.mkdir()
        counts, lengths = _rewrite_annotation(
            ann_path, stage / "update.ann.txt", mapping, retain_utr_features,
        )
        sequences = _rewrite_fasta(fasta_path, stage / "update.fasta", mapping)
        if any(sequences[entry]["length"] != length for entry, length in lengths.items()):
            raise MSSPackError("Annotation source lengths and FASTA lengths differ")
        manifest = {
            "annotation_input": str(ann_path.resolve()),
            "fasta_input": str(fasta_path.resolve()),
            "mapping_input": str(mapping_path.resolve()),
            "feature_counts": dict(counts),
            "sequences": {
                entry: {"accession": mapping[entry][0], "submitter_seqid": mapping[entry][1], **info}
                for entry, info in sequences.items()
            },
            "validation": "Not run. Run msspack validate and verify mapping against public records.",
        }
        (stage / "update-manifest.json").write_text(
            json.dumps(manifest, indent=2) + "\n", encoding="utf-8",
        )
        stage.rename(output_dir)
    return output_dir
