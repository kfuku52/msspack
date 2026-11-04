"""Read annotated organelles without nuclear-model correction or feature selection."""
from __future__ import annotations

import re
from collections import defaultdict
from pathlib import Path
from typing import Any

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from .config_models import OrganelleConfig
from .fasta import iter_fasta
from .gff import GFFRecord, child_ids, iter_gff_records
from .utils import MSSPackError, expand_path

GFF_FEATURES = {
    "gene", "CDS", "mRNA", "tRNA", "rRNA", "ncRNA", "misc_RNA", "misc_feature",
    "exon", "intron", "repeat_region", "regulatory", "pseudogene",
}


def _qualifiers(row: GFFRecord) -> dict[str, list[str]]:
    aliases = {"Note": "note", "Dbxref": "db_xref", "gene_name": "gene"}
    ignored = {"ID", "Parent", "Name", "Ontology_term", "Alias"}
    qualifiers: dict[str, list[str]] = {}
    for key, value in row.attributes.items():
        if key in ignored:
            continue
        if key in {"partial", "start_range", "end_range"}:
            raise MSSPackError(
                f"GFF fuzzy boundaries require explicit GenBank locations: {row.seqid}:{row.start}"
            )
        if not re.fullmatch(r"[a-z][a-z_]*", key) and key not in aliases:
            raise MSSPackError(f"Unsupported organelle GFF attribute {key!r}; use annotated GenBank")
        qualifiers[aliases.get(key, key)] = child_ids(value) if key == "Dbxref" else [value]
    if row.type in {"gene", "pseudogene", "tRNA", "rRNA", "ncRNA"}:
        if "gene" not in qualifiers and row.attributes.get("Name"):
            qualifiers["gene"] = [row.attributes["Name"]]
    return qualifiers


def _attach_gff(records: list[SeqRecord], path: Path) -> None:
    lookup = {str(record.id): record for record in records}
    rows = list(iter_gff_records(path))
    ids: dict[str, GFFRecord] = {}
    for row in rows:
        if row.seqid not in lookup or not 1 <= row.start <= row.end <= len(lookup[row.seqid]):
            raise MSSPackError(f"Organelle GFF coordinate/sequence mismatch: {row.to_line()}")
        if row.type == "source":
            continue
        if row.type not in GFF_FEATURES:
            raise MSSPackError(f"Unsupported organelle GFF feature {row.type!r}; use annotated GenBank")
        identifier = row.attributes.get("ID", "")
        if identifier:
            previous = ids.get(identifier)
            if previous and (previous.type != row.type or previous.seqid != row.seqid):
                raise MSSPackError(f"Ambiguous organelle GFF ID: {identifier}")
            ids[identifier] = row
    groups: dict[tuple[str, str, str], list[GFFRecord]] = defaultdict(list)
    for index, row in enumerate(rows):
        if row.type == "source":
            continue
        parents = child_ids(row.attributes.get("Parent"))
        if len(parents) > 1:
            raise MSSPackError("Multi-parent organelle GFF features require annotated GenBank")
        for parent_key in parents:
            if parent_key not in ids:
                raise MSSPackError(f"Missing organelle GFF parent: {parent_key}")
        key = row.attributes.get("ID", f"row{index}")
        if row.type == "CDS" and parents:
            # GFF3 CDS rows sharing a transcript form one ordered coding feature.
            parent = ids[parents[0]]
            key = parents[0] if parent.type == "mRNA" else row.attributes.get("ID", parents[0])
        groups[(row.seqid, row.type, key)].append(row)
    for (seqid, feature_type, key), parts in groups.items():
        strands = {row.strand for row in parts}
        if len(strands) != 1 or not strands <= {"+", "-"}:
            raise MSSPackError(f"Mixed/unknown GFF strands require annotated GenBank: {key}")
        strand = -1 if parts[0].strand == "-" else 1
        if len(parts) > 1 and any("trans_splicing" in row.attributes for row in parts):
            raise MSSPackError(f"Trans-spliced GFF must specify part order in annotated GenBank: {key}")
        parts.sort(key=lambda row: row.start, reverse=strand == -1)
        locations = [SimpleLocation(row.start - 1, row.end, strand=strand) for row in parts]
        location = locations[0] if len(locations) == 1 else CompoundLocation(locations, "join")
        qualifiers = _qualifiers(parts[0])
        if any(_qualifiers(row) != qualifiers for row in parts[1:]):
            raise MSSPackError(f"Inconsistent qualifiers in multi-part organelle GFF feature: {key}")
        parent_id = child_ids(parts[0].attributes.get("Parent"))
        visited: set[str] = set()
        gene_id = key
        while parent_id:
            parent_key = parent_id[0]
            if parent_key in visited:
                raise MSSPackError(f"Cycle in organelle GFF parents: {key}")
            visited.add(parent_key)
            parent = ids[parent_key]
            gene_id = parent_key
            if "gene" not in qualifiers:
                gene_name = parent.attributes.get("gene") or (
                    parent.attributes.get("Name") if parent.type in {"gene", "pseudogene"} else None
                )
                if gene_name:
                    qualifiers["gene"] = [gene_name]
            parent_id = child_ids(parent.attributes.get("Parent"))
            if len(parent_id) > 1 or (parent_id and parent_id[0] not in ids):
                raise MSSPackError(f"Ambiguous/missing organelle GFF ancestor: {key}")
        if feature_type == "CDS":
            if parts[0].phase not in {"0", "1", "2"}:
                raise MSSPackError(f"Organelle CDS requires an explicit GFF phase: {key}")
            phase = int(parts[0].phase)
            qualifiers["codon_start"] = [str(phase + 1)]
            length = len(locations[0]) - phase
            for row, part in zip(parts[1:], locations[1:], strict=True):
                if row.phase != str((-length) % 3):
                    raise MSSPackError(f"Inconsistent organelle CDS phase: {key}")
                length += len(part)
        if feature_type == "pseudogene":
            feature_type = "gene"
            qualifiers.setdefault("pseudogene", ["unknown"])
        feature = SeqFeature(location, type=feature_type, id=gene_id, qualifiers=qualifiers)
        lookup[seqid].features.append(feature)


def load_organelle_records(item: OrganelleConfig, base: Path) -> list[SeqRecord]:
    if item.genbank:
        with expand_path(item.genbank, base).open(encoding="utf-8") as handle:
            records = list(SeqIO.parse(handle, "genbank"))
    else:
        records = [SeqRecord(Seq(row.sequence), id=row.id, description=row.description)
                   for row in iter_fasta(expand_path(item.fasta, base))]
    if not records:
        raise MSSPackError(f"No sequences for organelle {item.name}")
    seen: set[str] = set()
    for record in records:
        identifier = str(record.id)
        if identifier in seen:
            raise MSSPackError(f"Duplicate organelle sequence ID: {identifier}")
        seen.add(identifier)
        sequence = str(record.seq)
        if not sequence or set(sequence.upper()) - set("ACGTRYSWKMBDHVN"):
            raise MSSPackError(f"Invalid/empty organelle DNA sequence: {identifier}")
        for feature in record.features:
            if feature.type == "source":
                values = feature.qualifiers.get("organelle", [])
                if values and values != [item.organelle]:
                    raise MSSPackError(f"Organelle source conflicts with configuration: {identifier}")
    for name, metadata in (("topology", item.topology), ("completeness", item.completeness)):
        missing = set(metadata) - seen
        if missing:
            raise MSSPackError(f"Unknown sequence ID in organelle {name}: {', '.join(sorted(missing))}")
    if item.gff:
        _attach_gff(records, expand_path(item.gff, base))
    return records


def feature_qualifiers(feature: Any) -> dict[str, list[str]]:
    return {str(key): [str(value) for value in values]
            for key, values in feature.qualifiers.items()}
