from __future__ import annotations

import re
from typing import Protocol

from .models import FeatureRecord
from .render import render_gap_feature

GapRegions = list[tuple[int, int]]
GAP_PATTERN = re.compile("N+")


class SequenceLike(Protocol):
    seq: object


def detect_gap_regions(
    sequence: str | SequenceLike,
    *,
    linkage_evidence: str,
    min_assembly_gap_size: int,
    gap_estimated_length: str,
) -> tuple[str, GapRegions]:
    if not isinstance(sequence, str):
        sequence = str(sequence.seq)
    gap_regions: GapRegions = [
        (match.start() + 1, match.end())
        for match in GAP_PATTERN.finditer(sequence)
    ]
    chunks: list[str] = []
    for start, end in gap_regions:
        chunks.append(
            render_gap_feature(
                gap_start=start,
                gap_end=end,
                linkage_evidence=linkage_evidence,
                min_assembly_gap_size=min_assembly_gap_size,
                gap_estimated_length=gap_estimated_length,
            )
        )
    return "".join(chunks), gap_regions


def compare_exon_to_gaps(
    exon_feature: FeatureRecord,
    gap_regions: GapRegions,
    strand: str,
) -> tuple[list[tuple[int, int]], bool, bool, bool]:
    new_starts: list[int] = []
    new_ends: list[int] = []
    gap_flag = False
    gap_start_flag = False
    gap_end_flag = False

    gap_count = 0
    exon_start = exon_feature.start
    exon_end = exon_feature.end

    if any(start <= exon_start <= exon_end <= end for start, end in gap_regions):
        # There is no flanking sequence to split out. Retain this uncertain
        # segment and let the configured gap policy mark the enclosing feature.
        return [(exon_start, exon_end)], True, True, True

    for gap_start, gap_end in gap_regions:
        if gap_end < exon_start or exon_end < gap_start:
            continue
        gap_count += 1
        if exon_start < gap_start:
            new_starts.append(exon_start)
            if exon_start < gap_start and gap_end < exon_end:
                if gap_count > 1:
                    new_ends.pop()
                    new_starts.pop()
                gap_size_mod = (gap_end - gap_start + 1) % 3
                if gap_size_mod == 1:
                    if strand == "-":
                        exon_end_tmp = gap_start - 3
                        exon_start = gap_end + 1
                    else:
                        exon_end_tmp = gap_start - 1
                        exon_start = gap_end + 3
                elif gap_size_mod == 2:
                    if strand == "-":
                        exon_end_tmp = gap_start - 2
                        exon_start = gap_end + 1
                    else:
                        exon_end_tmp = gap_start - 1
                        exon_start = gap_end + 2
                else:
                    exon_end_tmp = gap_start - 1
                    exon_start = gap_end + 1
                new_ends.append(exon_end_tmp)
                new_starts.append(exon_start)
            elif exon_start < gap_start and exon_end <= gap_end:
                if gap_count > 1:
                    new_ends.pop()
                    new_starts.pop()
                exon_end = gap_start - 1
        elif gap_start <= exon_start and gap_end < exon_end:
            exon_start = gap_end + 1
            new_starts.append(exon_start)
            gap_start_flag = True
        new_ends.append(exon_end)
    if gap_count == 0:
        return [(exon_start, exon_end)], gap_flag, gap_start_flag, gap_end_flag
    gap_flag = True
    return (
        list(zip(map(int, new_starts), map(int, new_ends), strict=True)),
        gap_flag,
        gap_start_flag,
        gap_end_flag,
    )


def append_position(
    feature: FeatureRecord,
    count: int,
    position: str,
    gap_regions: GapRegions,
    strand: str,
) -> tuple[str, str, str, bool]:
    """Join surviving segments without relying on the number of input rows."""

    joint = ""
    joint_close = ""
    gap_segments, gap_flag, _, _ = compare_exon_to_gaps(feature, gap_regions, strand)
    for cds_start, cds_end in gap_segments:
        if cds_start > cds_end:
            continue
        segment = str(cds_start) if cds_start == cds_end else f"{cds_start}..{cds_end}"
        position += ("," if position else "") + segment
    if "," in position:
        joint, joint_close = "join(", ")"
    return position, joint, joint_close, gap_flag
