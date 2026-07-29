from __future__ import annotations

import csv
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from xml.sax.saxutils import escape, quoteattr

from .chart_primitives import (
    GRID_RGB,
    MUTED_RGB,
    SVG_FONT_SIZE,
    TEXT_RGB,
    hex_to_rgb,
    pdf_helvetica_text_width,
    pdf_text_command,
    pdf_top_to_bottom,
    write_multi_page_pdf,
)
from .coordinate_duplicates import (
    CoordinateDuplicatePair,
    CoordinateGene,
    identify_coordinate_duplicate_pairs,
    read_coordinate_duplicate_map,
)
from .gff import GFFRecord, child_ids, iter_gff_records
from .utils import MSSPackError, atomic_text_writer, write_text

PLOT_WIDTH = 518.4
PLOT_LEFT = 166.0
PLOT_RIGHT = 16.0
PLOT_AXIS_WIDTH = PLOT_WIDTH - PLOT_LEFT - PLOT_RIGHT
HEADER_HEIGHT = 72.0
TRACK_HEIGHT = 18.0
GROUP_BASE_HEIGHT = 52.0
GROUP_BOTTOM_GAP = 12.0
PDF_PAGE_HEIGHT = 720.0
PDF_BOTTOM_MARGIN = 24.0

KEPT_COLOR = "#2563eb"
REMOVED_COLOR = "#e11d48"
EXON_COLOR = "#cbd5e1"
UTR_COLOR = "#93c5fd"
GRID_COLOR = "#d1d5db"
TEXT_COLOR = "#111827"
MUTED_COLOR = "#4b5563"

TRANSCRIPT_TYPES = {
    "mrna",
    "transcript",
    "primary_transcript",
    "ncrna",
    "trna",
    "rrna",
    "tmrna",
    "mirna",
    "snrna",
    "snorna",
    "lncrna",
}

PLOT_TSV_FIELDS = (
    "removal_order",
    "group_id",
    "seqid",
    "locus_start",
    "locus_end",
    "strand",
    "kept_gene_id",
    "removed_gene_id",
    "kept_transcript_id",
    "removed_transcript_id",
    "kept_cds_length",
    "removed_cds_length",
    "kept_intron_count",
    "removed_intron_count",
    "kept_internal_stops",
    "removed_internal_stops",
    "kept_complete",
    "removed_complete",
    "kept_splice_motifs",
    "removed_splice_motifs",
    "selection_policy",
    "selection_reason",
    "low_confidence",
    "decision",
    "gene_id",
    "track_id",
    "track_order",
    "feature_type",
    "feature_start",
    "feature_end",
    "phase",
)


@dataclass(frozen=True)
class CoordinateDuplicateFeature:
    removal_order: int
    group_id: str
    seqid: str
    locus_start: int
    locus_end: int
    strand: str
    kept_gene_id: str
    removed_gene_id: str
    kept_transcript_id: str
    removed_transcript_id: str
    kept_cds_length: int
    removed_cds_length: int
    kept_intron_count: int
    removed_intron_count: int
    kept_internal_stops: int
    removed_internal_stops: int
    kept_complete: bool
    removed_complete: bool
    kept_splice_motifs: tuple[str, ...]
    removed_splice_motifs: tuple[str, ...]
    selection_policy: str
    selection_reason: str
    low_confidence: bool
    decision: str
    gene_id: str
    track_id: str
    track_order: int
    feature_type: str
    feature_start: int
    feature_end: int
    phase: str


@dataclass(frozen=True)
class CoordinateDuplicatePlotData:
    pairs: tuple[CoordinateDuplicatePair, ...]
    features: tuple[CoordinateDuplicateFeature, ...]
    source_gff: Path
    duplicate_map: Path | None

    @property
    def total_removed_genes(self) -> int:
        return len(self.pairs)


@dataclass(frozen=True)
class _Track:
    decision: str
    gene_id: str
    track_id: str
    track_order: int
    removal_order: int
    features: tuple[CoordinateDuplicateFeature, ...]


@dataclass(frozen=True)
class _Group:
    group_id: str
    seqid: str
    start: int
    end: int
    strand: str
    kept_gene_id: str
    removed_gene_ids: tuple[str, ...]
    tracks: tuple[_Track, ...]

    @property
    def height(self) -> float:
        return GROUP_BASE_HEIGHT + len(self.tracks) * TRACK_HEIGHT + GROUP_BOTTOM_GAP


def _gene_records(records: list[GFFRecord]) -> tuple[CoordinateGene, ...]:
    genes: list[CoordinateGene] = []
    for record in records:
        if record.type != "gene":
            continue
        gene_id = record.attributes.get("ID")
        if not gene_id:
            continue
        genes.append(
            CoordinateGene(
                gene_id=gene_id,
                seqid=record.seqid,
                start=record.start,
                end=record.end,
                strand=record.strand,
            )
        )
    return tuple(genes)


def _record_id(record: GFFRecord) -> str:
    return record.attributes.get("ID", "")


def _collect_track_records(
    *,
    root: tuple[int, GFFRecord],
    children_by_parent: dict[str, list[tuple[int, GFFRecord]]],
) -> list[tuple[int, GFFRecord]]:
    collected: list[tuple[int, GFFRecord]] = []
    seen_records: set[int] = set()
    queued_ids: list[str] = []

    def add(index: int, record: GFFRecord) -> None:
        if index in seen_records:
            return
        seen_records.add(index)
        collected.append((index, record))
        record_id = _record_id(record)
        if record_id:
            queued_ids.append(record_id)

    add(*root)
    seen_ids: set[str] = set()
    while queued_ids:
        parent_id = queued_ids.pop(0)
        if parent_id in seen_ids:
            continue
        seen_ids.add(parent_id)
        for child in children_by_parent.get(parent_id, []):
            add(*child)
    return sorted(collected, key=lambda item: item[0])


def _tracks_for_gene(
    *,
    gene_id: str,
    gene_record: GFFRecord,
    children_by_parent: dict[str, list[tuple[int, GFFRecord]]],
) -> tuple[tuple[str, tuple[tuple[int, GFFRecord], ...]], ...]:
    direct_children = children_by_parent.get(gene_id, [])
    transcript_roots: list[tuple[int, GFFRecord]] = []
    gene_level_roots: list[tuple[int, GFFRecord]] = []
    for child in direct_children:
        if child[1].type.casefold() in TRANSCRIPT_TYPES and _record_id(child[1]):
            transcript_roots.append(child)
        else:
            gene_level_roots.append(child)

    tracks: list[tuple[str, tuple[tuple[int, GFFRecord], ...]]] = []
    for index, root in transcript_roots:
        track_id = _record_id(root) or f"{gene_id}:transcript:{index}"
        tracks.append(
            (
                track_id,
                tuple(
                    _collect_track_records(
                        root=(index, root),
                        children_by_parent=children_by_parent,
                    )
                ),
            )
        )
    if gene_level_roots:
        gene_level_records: dict[int, GFFRecord] = {}
        for gene_level_root in gene_level_roots:
            for index, record in _collect_track_records(
                root=gene_level_root,
                children_by_parent=children_by_parent,
            ):
                gene_level_records[index] = record
        tracks.append((gene_id, tuple(sorted(gene_level_records.items()))))
    if not tracks:
        tracks.append((gene_id, ((-1, gene_record),)))
    return tuple(tracks)


def _validate_pair_gene(
    pair: CoordinateDuplicatePair,
    gene_id: str,
    gene_record: GFFRecord | None,
    source_gff: Path,
) -> GFFRecord:
    if gene_record is None:
        raise MSSPackError(
            f"Coordinate duplicate gene {gene_id} is missing from {source_gff}"
        )
    observed = (
        gene_record.seqid,
        gene_record.start,
        gene_record.end,
        gene_record.strand,
    )
    expected = (pair.seqid, pair.start, pair.end, pair.strand)
    if observed != expected:
        raise MSSPackError(
            f"Coordinate duplicate map does not match gene {gene_id} in {source_gff}"
        )
    return gene_record


def build_coordinate_duplicate_plot_data(
    *,
    source_gff: Path,
    duplicate_map: Path | None,
) -> CoordinateDuplicatePlotData:
    if duplicate_map is not None and duplicate_map.exists():
        pairs = read_coordinate_duplicate_map(duplicate_map)
        if not pairs and not source_gff.exists():
            return CoordinateDuplicatePlotData(
                pairs=(),
                features=(),
                source_gff=source_gff,
                duplicate_map=duplicate_map,
            )
    else:
        if not source_gff.exists():
            raise MSSPackError(
                "Coordinate duplicate plot requires the pre-deduplication GFF: "
                f"{source_gff}"
            )
        records_for_pairs = list(iter_gff_records(source_gff))
        pairs = identify_coordinate_duplicate_pairs(_gene_records(records_for_pairs))

    if not pairs:
        return CoordinateDuplicatePlotData(
            pairs=(),
            features=(),
            source_gff=source_gff,
            duplicate_map=duplicate_map if duplicate_map and duplicate_map.exists() else None,
        )
    if not source_gff.exists():
        raise MSSPackError(
            f"Coordinate duplicate plot source GFF is missing: {source_gff}"
        )

    records = list(iter_gff_records(source_gff))
    gene_by_id: dict[str, GFFRecord] = {}
    children_by_parent: dict[str, list[tuple[int, GFFRecord]]] = defaultdict(list)
    for index, record in enumerate(records):
        record_id = _record_id(record)
        if record.type == "gene" and record_id:
            gene_by_id[record_id] = record
        for parent_id in child_ids(record.attributes.get("Parent")):
            children_by_parent[parent_id].append((index, record))

    track_cache: dict[
        str,
        tuple[tuple[str, tuple[tuple[int, GFFRecord], ...]], ...],
    ] = {}
    feature_rows: list[CoordinateDuplicateFeature] = []
    for removal_order, pair in enumerate(pairs, start=1):
        for decision, gene_id in (
            ("kept", pair.kept_gene_id),
            ("removed", pair.removed_gene_id),
        ):
            gene_record = _validate_pair_gene(
                pair,
                gene_id,
                gene_by_id.get(gene_id),
                source_gff,
            )
            tracks = track_cache.get(gene_id)
            if tracks is None:
                tracks = _tracks_for_gene(
                    gene_id=gene_id,
                    gene_record=gene_record,
                    children_by_parent=children_by_parent,
                )
                track_cache[gene_id] = tracks
            for track_order, (track_id, track_records) in enumerate(tracks):
                for _record_order, record in track_records:
                    feature_rows.append(
                        CoordinateDuplicateFeature(
                            removal_order=removal_order,
                            group_id=pair.group_id,
                            seqid=pair.seqid,
                            locus_start=pair.start,
                            locus_end=pair.end,
                            strand=pair.strand,
                            kept_gene_id=pair.kept_gene_id,
                            removed_gene_id=pair.removed_gene_id,
                            kept_transcript_id=pair.kept_transcript_id,
                            removed_transcript_id=pair.removed_transcript_id,
                            kept_cds_length=pair.kept_cds_length,
                            removed_cds_length=pair.removed_cds_length,
                            kept_intron_count=pair.kept_intron_count,
                            removed_intron_count=pair.removed_intron_count,
                            kept_internal_stops=pair.kept_internal_stops,
                            removed_internal_stops=pair.removed_internal_stops,
                            kept_complete=pair.kept_complete,
                            removed_complete=pair.removed_complete,
                            kept_splice_motifs=pair.kept_splice_motifs,
                            removed_splice_motifs=pair.removed_splice_motifs,
                            selection_policy=pair.selection_policy,
                            selection_reason=pair.selection_reason,
                            low_confidence=pair.low_confidence,
                            decision=decision,
                            gene_id=gene_id,
                            track_id=track_id,
                            track_order=track_order,
                            feature_type=record.type,
                            feature_start=record.start,
                            feature_end=record.end,
                            phase=record.phase,
                        )
                    )
    return CoordinateDuplicatePlotData(
        pairs=pairs,
        features=tuple(feature_rows),
        source_gff=source_gff,
        duplicate_map=duplicate_map if duplicate_map and duplicate_map.exists() else None,
    )


def write_coordinate_duplicate_plot_tsv(
    data: CoordinateDuplicatePlotData,
    output_path: Path,
) -> Path:
    with atomic_text_writer(output_path) as handle:
        writer = csv.DictWriter(handle, fieldnames=PLOT_TSV_FIELDS, delimiter="\t")
        writer.writeheader()
        for feature in data.features:
            writer.writerow(
                {
                    "removal_order": feature.removal_order,
                    "group_id": feature.group_id,
                    "seqid": feature.seqid,
                    "locus_start": feature.locus_start,
                    "locus_end": feature.locus_end,
                    "strand": feature.strand,
                    "kept_gene_id": feature.kept_gene_id,
                    "removed_gene_id": feature.removed_gene_id,
                    "kept_transcript_id": feature.kept_transcript_id,
                    "removed_transcript_id": feature.removed_transcript_id,
                    "kept_cds_length": feature.kept_cds_length,
                    "removed_cds_length": feature.removed_cds_length,
                    "kept_intron_count": feature.kept_intron_count,
                    "removed_intron_count": feature.removed_intron_count,
                    "kept_internal_stops": feature.kept_internal_stops,
                    "removed_internal_stops": feature.removed_internal_stops,
                    "kept_complete": str(feature.kept_complete).lower(),
                    "removed_complete": str(feature.removed_complete).lower(),
                    "kept_splice_motifs": ",".join(feature.kept_splice_motifs),
                    "removed_splice_motifs": ",".join(feature.removed_splice_motifs),
                    "selection_policy": feature.selection_policy,
                    "selection_reason": feature.selection_reason,
                    "low_confidence": str(feature.low_confidence).lower(),
                    "decision": feature.decision,
                    "gene_id": feature.gene_id,
                    "track_id": feature.track_id,
                    "track_order": feature.track_order,
                    "feature_type": feature.feature_type,
                    "feature_start": feature.feature_start,
                    "feature_end": feature.feature_end,
                    "phase": feature.phase,
                }
            )
    return output_path


def _build_groups(data: CoordinateDuplicatePlotData, limit: int) -> tuple[_Group, ...]:
    if limit <= 0:
        raise MSSPackError("Coordinate duplicate plot limit must be greater than zero")
    selected_pairs = data.pairs[:limit]
    selected_orders = set(range(1, len(selected_pairs) + 1))
    rows_by_group: dict[str, list[CoordinateDuplicateFeature]] = defaultdict(list)
    for feature in data.features:
        if feature.removal_order in selected_orders:
            rows_by_group[feature.group_id].append(feature)

    pairs_by_group: dict[str, list[tuple[int, CoordinateDuplicatePair]]] = defaultdict(list)
    group_order: list[str] = []
    for removal_order, pair in enumerate(selected_pairs, start=1):
        if pair.group_id not in pairs_by_group:
            group_order.append(pair.group_id)
        pairs_by_group[pair.group_id].append((removal_order, pair))

    groups: list[_Group] = []
    for group_id in group_order:
        ordered_pairs = pairs_by_group[group_id]
        first_pair = ordered_pairs[0][1]
        removed_order = {
            pair.removed_gene_id: order for order, pair in ordered_pairs
        }
        track_rows: dict[
            tuple[str, str, str],
            list[CoordinateDuplicateFeature],
        ] = defaultdict(list)
        seen_features: set[
            tuple[str, str, str, str, int, int, str]
        ] = set()
        for feature in rows_by_group[group_id]:
            feature_key = (
                feature.decision,
                feature.gene_id,
                feature.track_id,
                feature.feature_type,
                feature.feature_start,
                feature.feature_end,
                feature.phase,
            )
            if feature_key in seen_features:
                continue
            seen_features.add(feature_key)
            track_rows[(feature.decision, feature.gene_id, feature.track_id)].append(
                feature
            )
        tracks: list[_Track] = []
        for (decision, gene_id, track_id), features in track_rows.items():
            removal_order = 0 if decision == "kept" else removed_order[gene_id]
            tracks.append(
                _Track(
                    decision=decision,
                    gene_id=gene_id,
                    track_id=track_id,
                    track_order=min(feature.track_order for feature in features),
                    removal_order=removal_order,
                    features=tuple(features),
                )
            )
        tracks.sort(
            key=lambda track: (
                0 if track.decision == "kept" else 1,
                track.removal_order,
                track.track_order,
                track.track_id,
            )
        )
        groups.append(
            _Group(
                group_id=group_id,
                seqid=first_pair.seqid,
                start=first_pair.start,
                end=first_pair.end,
                strand=first_pair.strand,
                kept_gene_id=first_pair.kept_gene_id,
                removed_gene_ids=tuple(pair.removed_gene_id for _, pair in ordered_pairs),
                tracks=tuple(tracks),
            )
        )
    return tuple(groups)


def coordinate_duplicate_plot_summary(
    data: CoordinateDuplicatePlotData,
    limit: int,
) -> dict[str, object]:
    shown = min(data.total_removed_genes, limit)
    return {
        "limit": limit,
        "total_removed_genes": data.total_removed_genes,
        "shown_removed_genes": shown,
        "truncated": shown < data.total_removed_genes,
    }


def _subtitle(data: CoordinateDuplicatePlotData, limit: int) -> str:
    shown = min(data.total_removed_genes, limit)
    if shown < data.total_removed_genes:
        gene_label = "gene" if data.total_removed_genes == 1 else "genes"
        return (
            f"Showing first {shown:,} of {data.total_removed_genes:,} removed {gene_label} "
            f"(plot limit {limit:,})."
        )
    gene_label = "gene" if shown == 1 else "genes"
    return f"Showing {shown:,} removed {gene_label} (plot limit {limit:,})."


def _shorten(text: str, length: int = 22) -> str:
    if len(text) <= length:
        return text
    return text[: length - 3] + "..."


def _coordinate_x(position: int, start: int, end: int) -> float:
    span = max(1, end - start)
    fraction = (position - start) / span
    return PLOT_LEFT + max(0.0, min(1.0, fraction)) * PLOT_AXIS_WIDTH


def _ticks(start: int, end: int) -> tuple[int, ...]:
    if start == end:
        return (start,)
    return tuple(round(start + (end - start) * index / 4) for index in range(5))


def _track_label(track: _Track) -> str:
    if track.track_id == track.gene_id:
        return track.gene_id
    return f"{track.gene_id} / {track.track_id}"


def _track_full_label(track: _Track) -> str:
    status = "KEPT" if track.decision == "kept" else "REMOVED"
    return f"{status} {_track_label(track)}"


def _track_span(track: _Track, group: _Group) -> tuple[int, int]:
    transcript_features = [
        feature
        for feature in track.features
        if feature.feature_type.casefold() in TRANSCRIPT_TYPES
    ]
    if transcript_features:
        return (
            min(feature.feature_start for feature in transcript_features),
            max(feature.feature_end for feature in transcript_features),
        )
    visible_features = [
        feature
        for feature in track.features
        if _feature_style(feature.feature_type) is not None
    ]
    if visible_features:
        return (
            min(feature.feature_start for feature in visible_features),
            max(feature.feature_end for feature in visible_features),
        )
    return group.start, group.end


def _feature_style(feature_type: str) -> str | None:
    normalized = feature_type.casefold()
    if normalized == "cds":
        return "cds"
    if normalized == "exon":
        return "exon"
    if "utr" in normalized:
        return "utr"
    if normalized in {"start_codon", "stop_codon"}:
        return "codon"
    return None


def _ordered_visible_features(
    track: _Track,
) -> tuple[CoordinateDuplicateFeature, ...]:
    priority = {"exon": 0, "utr": 1, "cds": 2, "codon": 3}
    return tuple(
        sorted(
            (
                feature
                for feature in track.features
                if _feature_style(feature.feature_type) is not None
            ),
            key=lambda feature: (
                priority[_feature_style(feature.feature_type) or "exon"],
                feature.feature_start,
                feature.feature_end,
            ),
        )
    )


def _svg_feature(
    feature: CoordinateDuplicateFeature,
    *,
    y: float,
    color: str,
    group: _Group,
) -> str:
    style = _feature_style(feature.feature_type)
    start_x = _coordinate_x(feature.feature_start, group.start, group.end)
    end_x = _coordinate_x(feature.feature_end, group.start, group.end)
    width = max(1.2, end_x - start_x)
    if style == "codon":
        return (
            f'<line x1="{start_x:.2f}" y1="{y - 6:.2f}" '
            f'x2="{start_x:.2f}" y2="{y + 6:.2f}" stroke="{color}" stroke-width="1"/>'
        )
    if style == "cds":
        return (
            f'<rect x="{start_x:.2f}" y="{y - 5:.2f}" width="{width:.2f}" height="10" '
            f'fill="{color}" stroke="{color}" stroke-width="0.6"/>'
        )
    if style == "utr":
        return (
            f'<rect x="{start_x:.2f}" y="{y - 3:.2f}" width="{width:.2f}" height="6" '
            f'fill="{UTR_COLOR}" stroke="{color}" stroke-width="0.6"/>'
        )
    return (
        f'<rect x="{start_x:.2f}" y="{y - 3:.2f}" width="{width:.2f}" height="6" '
        f'fill="{EXON_COLOR}" stroke="{color}" stroke-width="0.6"/>'
    )


def _svg_arrow(x: float, y: float, strand: str, color: str) -> str:
    if strand == "-":
        points = f"{x:.2f},{y:.2f} {x + 5:.2f},{y - 3:.2f} {x + 5:.2f},{y + 3:.2f}"
    else:
        points = f"{x:.2f},{y:.2f} {x - 5:.2f},{y - 3:.2f} {x - 5:.2f},{y + 3:.2f}"
    return f'<polygon points="{points}" fill="{color}"/>'


def _svg_group(group: _Group, top: float) -> list[str]:
    lines = [
        (
            f'<text x="12" y="{top + 10:.2f}" class="group-title">'
            f"{escape(group.seqid)}:{group.start:,}-{group.end:,} "
            f"({escape(group.strand)})</text>"
        ),
        (
            f'<text x="{PLOT_WIDTH - PLOT_RIGHT:.2f}" y="{top + 10:.2f}" '
            f'class="muted" text-anchor="end">kept 1; removed '
            f"{len(group.removed_gene_ids)}</text>"
        ),
    ]
    axis_y = top + 29.0
    lines.append(
        f'<line x1="{PLOT_LEFT:.2f}" y1="{axis_y:.2f}" '
        f'x2="{PLOT_WIDTH - PLOT_RIGHT:.2f}" y2="{axis_y:.2f}" '
        f'stroke="{GRID_COLOR}" stroke-width="1"/>'
    )
    ticks = _ticks(group.start, group.end)
    for tick_index, tick in enumerate(ticks):
        x = _coordinate_x(tick, group.start, group.end)
        text_anchor = (
            "start"
            if tick_index == 0
            else "end"
            if tick_index == len(ticks) - 1
            else "middle"
        )
        lines.extend(
            [
                f'<line x1="{x:.2f}" y1="{axis_y - 3:.2f}" x2="{x:.2f}" '
                f'y2="{axis_y + 3:.2f}" stroke="{GRID_COLOR}" stroke-width="0.8"/>',
                f'<text x="{x:.2f}" y="{axis_y - 6:.2f}" class="muted" '
                f'text-anchor="{text_anchor}">{tick:,}</text>',
            ]
        )
    track_top = top + 47.0
    for track_index, track in enumerate(group.tracks):
        y = track_top + track_index * TRACK_HEIGHT
        color = KEPT_COLOR if track.decision == "kept" else REMOVED_COLOR
        label = _shorten(_track_label(track))
        span_start, span_end = _track_span(track, group)
        start_x = _coordinate_x(span_start, group.start, group.end)
        end_x = _coordinate_x(span_end, group.start, group.end)
        lines.extend(
            [
                f'<rect x="12" y="{y - 4:.2f}" width="6" height="6" fill="{color}"/>',
                f'<text x="22" y="{y + 3:.2f}" class="track-label" '
                f"data-full-label={quoteattr(_track_full_label(track))}>{escape(label)}</text>",
                f'<line x1="{start_x:.2f}" y1="{y:.2f}" x2="{end_x:.2f}" '
                f'y2="{y:.2f}" stroke="{color}" stroke-width="1"/>',
                _svg_arrow(
                    start_x if group.strand == "-" else end_x,
                    y,
                    group.strand,
                    color,
                ),
            ]
        )
        for feature in _ordered_visible_features(track):
            lines.append(_svg_feature(feature, y=y, color=color, group=group))
    separator_y = top + group.height - 3.0
    lines.append(
        f'<line x1="12" y1="{separator_y:.2f}" x2="{PLOT_WIDTH - PLOT_RIGHT:.2f}" '
        f'y2="{separator_y:.2f}" stroke="{GRID_COLOR}" stroke-width="0.6"/>'
    )
    return lines


def write_coordinate_duplicate_svg(
    data: CoordinateDuplicatePlotData,
    output_path: Path,
    *,
    limit: int,
) -> Path:
    groups = _build_groups(data, limit)
    height = max(126.0, HEADER_HEIGHT + sum(group.height for group in groups) + 12.0)
    shown = min(data.total_removed_genes, limit)
    svg = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        (
            f'<svg xmlns="http://www.w3.org/2000/svg" width="7.2in" '
            f'height="{height / 72:.3f}in" viewBox="0 0 {PLOT_WIDTH:.2f} {height:.2f}">'
        ),
        "<style>",
        (
            f"text{{font-family:Helvetica,Arial,sans-serif;font-size:{SVG_FONT_SIZE};"
            f"fill:{TEXT_COLOR}}}"
        ),
        ".title{font-weight:bold}",
        f".muted{{fill:{MUTED_COLOR}}}",
        ".group-title{font-weight:bold}",
        "</style>",
        '<text x="12" y="18" class="title">Coordinate duplicate gene models</text>',
        f'<text x="12" y="34" class="muted">{escape(_subtitle(data, limit))}</text>',
        f'<rect x="12" y="46" width="7" height="7" fill="{KEPT_COLOR}"/>',
        '<text x="23" y="53">Kept</text>',
        f'<rect x="68" y="46" width="7" height="7" fill="{REMOVED_COLOR}"/>',
        '<text x="79" y="53">Removed</text>',
        f'<rect x="145" y="47" width="14" height="6" fill="{EXON_COLOR}" '
        f'stroke="{MUTED_COLOR}" stroke-width="0.6"/>',
        '<text x="164" y="53">Exon</text>',
        f'<rect x="208" y="45" width="14" height="10" fill="{KEPT_COLOR}"/>',
        '<text x="227" y="53">CDS</text>',
    ]
    if not groups:
        svg.append(
            '<text x="12" y="94" class="muted">No coordinate duplicate genes detected.</text>'
        )
    else:
        top = HEADER_HEIGHT
        for group in groups:
            svg.extend(_svg_group(group, top))
            top += group.height
    svg.extend(
        [
            (
                f'<metadata>removed_genes_total={data.total_removed_genes};'
                f"removed_genes_shown={shown};limit={limit}</metadata>"
            ),
            "</svg>",
            "",
        ]
    )
    return write_text(output_path, "\n".join(svg))


def _pdf_safe(text: str) -> str:
    return text.encode("latin-1", errors="replace").decode("latin-1")


def _pdf_line(
    *,
    page_height: float,
    x1: float,
    y1: float,
    x2: float,
    y2: float,
    color: tuple[float, float, float],
    width: float = 1.0,
) -> str:
    bottom_y1 = pdf_top_to_bottom(page_height, y1)
    bottom_y2 = pdf_top_to_bottom(page_height, y2)
    return (
        f"{color[0]:.3f} {color[1]:.3f} {color[2]:.3f} RG {width:.2f} w "
        f"{x1:.2f} {bottom_y1:.2f} m {x2:.2f} {bottom_y2:.2f} l S"
    )


def _pdf_rect(
    *,
    page_height: float,
    x: float,
    y: float,
    width: float,
    height: float,
    fill: tuple[float, float, float],
    stroke: tuple[float, float, float] | None = None,
    stroke_width: float = 0.6,
) -> str:
    bottom_y = pdf_top_to_bottom(page_height, y, height)
    command = (
        f"{fill[0]:.3f} {fill[1]:.3f} {fill[2]:.3f} rg "
        f"{x:.2f} {bottom_y:.2f} {width:.2f} {height:.2f} re "
    )
    if stroke is None:
        return command + "f"
    return (
        command
        + f"{stroke[0]:.3f} {stroke[1]:.3f} {stroke[2]:.3f} RG "
        + f"{stroke_width:.2f} w B"
    )


def _pdf_arrow(
    *,
    page_height: float,
    x: float,
    y: float,
    strand: str,
    color: tuple[float, float, float],
) -> str:
    if strand == "-":
        points = ((x, y), (x + 5, y - 3), (x + 5, y + 3))
    else:
        points = ((x, y), (x - 5, y - 3), (x - 5, y + 3))
    converted = [(point_x, pdf_top_to_bottom(page_height, point_y)) for point_x, point_y in points]
    return (
        f"{color[0]:.3f} {color[1]:.3f} {color[2]:.3f} rg "
        f"{converted[0][0]:.2f} {converted[0][1]:.2f} m "
        f"{converted[1][0]:.2f} {converted[1][1]:.2f} l "
        f"{converted[2][0]:.2f} {converted[2][1]:.2f} l h f"
    )


def _pdf_group(group: _Group, top: float, page_height: float) -> list[str]:
    commands = [
        pdf_text_command(
            page_height=page_height,
            x=12,
            y_top=top + 10,
            text=_pdf_safe(
                f"{group.seqid}:{group.start:,}-{group.end:,} ({group.strand})"
            ),
            font="F2",
            size=8,
            color=TEXT_RGB,
        ),
        pdf_text_command(
            page_height=page_height,
            x=PLOT_WIDTH - 103,
            y_top=top + 10,
            text=f"kept 1; removed {len(group.removed_gene_ids)}",
            font="F1",
            size=8,
            color=MUTED_RGB,
        ),
    ]
    axis_y = top + 29.0
    commands.append(
        _pdf_line(
            page_height=page_height,
            x1=PLOT_LEFT,
            y1=axis_y,
            x2=PLOT_WIDTH - PLOT_RIGHT,
            y2=axis_y,
            color=GRID_RGB,
        )
    )
    ticks = _ticks(group.start, group.end)
    for tick_index, tick in enumerate(ticks):
        x = _coordinate_x(tick, group.start, group.end)
        tick_label = f"{tick:,}"
        label_width = pdf_helvetica_text_width(tick_label, size=7)
        label_x = (
            x
            if tick_index == 0
            else x - label_width
            if tick_index == len(ticks) - 1
            else x - label_width / 2
        )
        commands.extend(
            [
                _pdf_line(
                    page_height=page_height,
                    x1=x,
                    y1=axis_y - 3,
                    x2=x,
                    y2=axis_y + 3,
                    color=GRID_RGB,
                    width=0.8,
                ),
                pdf_text_command(
                    page_height=page_height,
                    x=label_x,
                    y_top=axis_y - 6,
                    text=tick_label,
                    font="F1",
                    size=7,
                    color=MUTED_RGB,
                ),
            ]
        )
    track_top = top + 47.0
    exon_rgb = hex_to_rgb(EXON_COLOR)
    utr_rgb = hex_to_rgb(UTR_COLOR)
    for track_index, track in enumerate(group.tracks):
        y = track_top + track_index * TRACK_HEIGHT
        color = hex_to_rgb(KEPT_COLOR if track.decision == "kept" else REMOVED_COLOR)
        commands.extend(
            [
                _pdf_rect(
                    page_height=page_height,
                    x=12,
                    y=y - 4,
                    width=6,
                    height=6,
                    fill=color,
                ),
                pdf_text_command(
                    page_height=page_height,
                    x=22,
                    y_top=y + 3,
                    text=_pdf_safe(_shorten(_track_label(track))),
                    font="F1",
                    size=7,
                    color=TEXT_RGB,
                ),
            ]
        )
        span_start, span_end = _track_span(track, group)
        start_x = _coordinate_x(span_start, group.start, group.end)
        end_x = _coordinate_x(span_end, group.start, group.end)
        commands.extend(
            [
                _pdf_line(
                    page_height=page_height,
                    x1=start_x,
                    y1=y,
                    x2=end_x,
                    y2=y,
                    color=color,
                ),
                _pdf_arrow(
                    page_height=page_height,
                    x=start_x if group.strand == "-" else end_x,
                    y=y,
                    strand=group.strand,
                    color=color,
                ),
            ]
        )
        for feature in _ordered_visible_features(track):
            style = _feature_style(feature.feature_type)
            feature_x = _coordinate_x(feature.feature_start, group.start, group.end)
            feature_end_x = _coordinate_x(feature.feature_end, group.start, group.end)
            feature_width = max(1.2, feature_end_x - feature_x)
            if style == "codon":
                commands.append(
                    _pdf_line(
                        page_height=page_height,
                        x1=feature_x,
                        y1=y - 6,
                        x2=feature_x,
                        y2=y + 6,
                        color=color,
                    )
                )
            elif style == "cds":
                commands.append(
                    _pdf_rect(
                        page_height=page_height,
                        x=feature_x,
                        y=y - 5,
                        width=feature_width,
                        height=10,
                        fill=color,
                    )
                )
            elif style == "utr":
                commands.append(
                    _pdf_rect(
                        page_height=page_height,
                        x=feature_x,
                        y=y - 3,
                        width=feature_width,
                        height=6,
                        fill=utr_rgb,
                        stroke=color,
                    )
                )
            else:
                commands.append(
                    _pdf_rect(
                        page_height=page_height,
                        x=feature_x,
                        y=y - 3,
                        width=feature_width,
                        height=6,
                        fill=exon_rgb,
                        stroke=color,
                    )
                )
    commands.append(
        _pdf_line(
            page_height=page_height,
            x1=12,
            y1=top + group.height - 3,
            x2=PLOT_WIDTH - PLOT_RIGHT,
            y2=top + group.height - 3,
            color=GRID_RGB,
            width=0.6,
        )
    )
    return commands


def _paginate_groups(groups: tuple[_Group, ...]) -> list[list[_Group]]:
    if not groups:
        return [[]]
    max_tracks = max(
        1,
        int(
            (
                PDF_PAGE_HEIGHT
                - PDF_BOTTOM_MARGIN
                - HEADER_HEIGHT
                - GROUP_BASE_HEIGHT
                - GROUP_BOTTOM_GAP
            )
            // TRACK_HEIGHT
        ),
    )
    page_sized_groups: list[_Group] = []
    for group in groups:
        if len(group.tracks) <= max_tracks:
            page_sized_groups.append(group)
            continue
        kept_tracks = tuple(
            track for track in group.tracks if track.decision == "kept"
        )
        removed_tracks = tuple(
            track for track in group.tracks if track.decision == "removed"
        )
        if kept_tracks and len(kept_tracks) < max_tracks:
            removed_chunk_size = max_tracks - len(kept_tracks)
            for offset in range(0, len(removed_tracks), removed_chunk_size):
                removed_chunk = removed_tracks[offset : offset + removed_chunk_size]
                page_sized_groups.append(
                    _Group(
                        group_id=group.group_id,
                        seqid=group.seqid,
                        start=group.start,
                        end=group.end,
                        strand=group.strand,
                        kept_gene_id=group.kept_gene_id,
                        removed_gene_ids=tuple(
                            dict.fromkeys(track.gene_id for track in removed_chunk)
                        ),
                        tracks=(*kept_tracks, *removed_chunk),
                    )
                )
        else:
            for offset in range(0, len(group.tracks), max_tracks):
                track_chunk = group.tracks[offset : offset + max_tracks]
                page_sized_groups.append(
                    _Group(
                        group_id=group.group_id,
                        seqid=group.seqid,
                        start=group.start,
                        end=group.end,
                        strand=group.strand,
                        kept_gene_id=group.kept_gene_id,
                        removed_gene_ids=tuple(
                            dict.fromkeys(
                                track.gene_id
                                for track in track_chunk
                                if track.decision == "removed"
                            )
                        ),
                        tracks=track_chunk,
                    )
                )
    pages: list[list[_Group]] = []
    current: list[_Group] = []
    used = HEADER_HEIGHT
    for group in page_sized_groups:
        if current and used + group.height > PDF_PAGE_HEIGHT - PDF_BOTTOM_MARGIN:
            pages.append(current)
            current = []
            used = HEADER_HEIGHT
        current.append(group)
        used += group.height
    if current:
        pages.append(current)
    return pages


def _pdf_header(
    *,
    data: CoordinateDuplicatePlotData,
    limit: int,
    page_height: float,
    page_number: int,
    page_count: int,
) -> list[str]:
    kept_rgb = hex_to_rgb(KEPT_COLOR)
    removed_rgb = hex_to_rgb(REMOVED_COLOR)
    return [
        pdf_text_command(
            page_height=page_height,
            x=12,
            y_top=18,
            text="Coordinate duplicate gene models",
            font="F2",
            size=10,
            color=TEXT_RGB,
        ),
        pdf_text_command(
            page_height=page_height,
            x=12,
            y_top=34,
            text=_subtitle(data, limit),
            font="F1",
            size=8,
            color=MUTED_RGB,
        ),
        pdf_text_command(
            page_height=page_height,
            x=PLOT_WIDTH - 54,
            y_top=18,
            text=f"{page_number}/{page_count}",
            font="F1",
            size=8,
            color=MUTED_RGB,
        ),
        _pdf_rect(
            page_height=page_height,
            x=12,
            y=46,
            width=7,
            height=7,
            fill=kept_rgb,
        ),
        pdf_text_command(
            page_height=page_height,
            x=23,
            y_top=53,
            text="Kept",
            font="F1",
            size=8,
            color=TEXT_RGB,
        ),
        _pdf_rect(
            page_height=page_height,
            x=68,
            y=46,
            width=7,
            height=7,
            fill=removed_rgb,
        ),
        pdf_text_command(
            page_height=page_height,
            x=79,
            y_top=53,
            text="Removed",
            font="F1",
            size=8,
            color=TEXT_RGB,
        ),
    ]


def write_coordinate_duplicate_pdf(
    data: CoordinateDuplicatePlotData,
    output_path: Path,
    *,
    limit: int,
) -> Path:
    groups = _build_groups(data, limit)
    page_groups = _paginate_groups(groups)
    pages: list[tuple[float, float, list[str]]] = []
    page_count = len(page_groups)
    for page_number, groups_on_page in enumerate(page_groups, start=1):
        content_height = HEADER_HEIGHT + sum(group.height for group in groups_on_page)
        page_height = max(PDF_PAGE_HEIGHT, content_height + PDF_BOTTOM_MARGIN)
        commands = _pdf_header(
            data=data,
            limit=limit,
            page_height=page_height,
            page_number=page_number,
            page_count=page_count,
        )
        if not groups_on_page:
            commands.append(
                pdf_text_command(
                    page_height=page_height,
                    x=12,
                    y_top=94,
                    text="No coordinate duplicate genes detected.",
                    font="F1",
                    size=8,
                    color=MUTED_RGB,
                )
            )
        else:
            top = HEADER_HEIGHT
            for group in groups_on_page:
                commands.extend(_pdf_group(group, top, page_height))
                top += group.height
        pages.append((PLOT_WIDTH, page_height, commands))
    return write_multi_page_pdf(pages=pages, output_path=output_path)
