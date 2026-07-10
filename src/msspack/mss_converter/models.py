from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Optional

from ..gff import child_ids, iter_gff_records


@dataclass(frozen=True)
class AnnotationEntry:
    product_name: str
    custom_locus_tag: Optional[str]


@dataclass(frozen=True)
class FeatureRecord:
    seq_id: str
    type: str
    id: str
    parent: str
    start: int
    end: int
    strand: str
    phase: int
    name: str
    rna_type: str
    anticodon: str


def _is_missing(value: Any) -> bool:
    return value is None or (isinstance(value, float) and math.isnan(value))


def _normalize_text(value: Any) -> str:
    if _is_missing(value):
        return ""
    return str(value)


def _normalize_phase(value: Any) -> int:
    if _is_missing(value) or value in ("", "."):
        return 0
    return int(value)


def feature_from_row(row: Any) -> FeatureRecord:
    if isinstance(row, FeatureRecord):
        return row
    if isinstance(row, dict):
        getter = row.get
    else:
        def getter(key: str, default: Any = None) -> Any:
            return getattr(row, key, default)
    return FeatureRecord(
        seq_id=_normalize_text(getter("seq_id")),
        type=_normalize_text(getter("type")),
        id=_normalize_text(getter("ID")),
        parent=_normalize_text(getter("Parent")),
        start=int(_normalize_text(getter("start"))),
        end=int(_normalize_text(getter("end"))),
        strand=_normalize_text(getter("strand")),
        phase=_normalize_phase(getter("phase", 0)),
        name=_normalize_text(getter("Name", "")),
        rna_type=_normalize_text(getter("Type", "")),
        anticodon=_normalize_text(getter("anticodon", "")),
    )


def build_gff_indexes(
    rows: Iterable[Any],
) -> tuple[dict[str, list[FeatureRecord]], dict[str, list[FeatureRecord]]]:
    parent_lookup: dict[str, list[FeatureRecord]] = {}
    gene_lookup: dict[str, list[FeatureRecord]] = {}
    for row in rows:
        feature = feature_from_row(row)
        if feature.type == "gene":
            gene_lookup.setdefault(feature.seq_id, []).append(feature)
        if feature.parent:
            for parent_id in child_ids(feature.parent):
                parent_lookup.setdefault(parent_id, []).append(feature)
    for seq_id in gene_lookup:
        gene_lookup[seq_id].sort(key=lambda feature: (feature.start, feature.end, feature.id))
    for parent_id in parent_lookup:
        parent_lookup[parent_id].sort(
            key=lambda feature: (feature.start, feature.end, feature.type, feature.id)
        )
    return gene_lookup, parent_lookup


def load_gff_features(path: str | Path) -> list[FeatureRecord]:
    features: list[FeatureRecord] = []
    for record in iter_gff_records(path):
        features.append(
            FeatureRecord(
                seq_id=record.seqid,
                type=record.type,
                id=record.attributes.get("ID", ""),
                parent=record.attributes.get("Parent", ""),
                start=record.start,
                end=record.end,
                strand=record.strand,
                phase=_normalize_phase(record.phase),
                name=record.attributes.get("Name", ""),
                rna_type=record.attributes.get("Type", ""),
                anticodon=record.attributes.get("anticodon", ""),
            )
        )
    features.sort(key=lambda feature: feature.start)
    return features


def load_annotation_lookup(path: str | Path) -> tuple[dict[str, AnnotationEntry], bool]:
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Annotation table has no header: {path}")
        required = {"ID", "Description"}
        missing = sorted(required.difference(reader.fieldnames))
        if missing:
            raise ValueError(f"Annotation table is missing columns: {', '.join(missing)}")
        has_custom_locus_tag = "Locus_tag" in reader.fieldnames
        lookup: dict[str, AnnotationEntry] = {}
        for row in reader:
            record_id = _normalize_text(row.get("ID"))
            if not record_id:
                continue
            custom_locus_tag = _normalize_text(row.get("Locus_tag")) or None
            lookup[record_id] = AnnotationEntry(
                product_name=_normalize_text(row.get("Description")),
                custom_locus_tag=custom_locus_tag,
            )
    return lookup, has_custom_locus_tag


def load_protein_id_lookup(path: Optional[str | Path]) -> Optional[dict[str, str]]:
    if not path:
        return None
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        if header is None or len(header) < 2:
            raise ValueError(f"Protein ID table must have at least two columns: {path}")
        lookup: dict[str, str] = {}
        for row in reader:
            if len(row) < 2:
                continue
            record_id = row[0].strip()
            protein_id = row[1].strip()
            if record_id and protein_id:
                lookup[record_id] = protein_id
    return lookup
