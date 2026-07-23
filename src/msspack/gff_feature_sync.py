from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable
from dataclasses import dataclass

from .gff import child_ids, parse_attributes
from .gff_feature_policy import (
    CDS_TYPES,
    CODON_TYPES,
    EXON_TYPES,
    INTRON_TYPES,
    START_CODON_TYPES,
    STOP_CODON_TYPES,
    UTR_TYPES,
)


@dataclass(frozen=True)
class HierarchyIssue:
    child_id: str
    parent_id: str
    message: str


def _row_id(row: list[str]) -> str:
    return parse_attributes(row[8]).get("ID", "")


def _row_parents(row: list[str]) -> list[str]:
    return child_ids(parse_attributes(row[8]).get("Parent"))


def _transcript_order(rows: Iterable[list[str]], strand: str) -> list[list[str]]:
    return sorted(rows, key=lambda row: (int(row[3]), int(row[4])), reverse=(strand == "-"))


def _terminal_spans(
    cdss: list[list[str]],
    *,
    strand: str,
    five_prime: bool,
    length: int = 3,
) -> list[tuple[int, int]]:
    ordered = _transcript_order(cdss, strand)
    if not five_prime:
        ordered.reverse()
    remaining = length
    spans: list[tuple[int, int]] = []
    for row in ordered:
        start = int(row[3])
        end = int(row[4])
        take = min(remaining, end - start + 1)
        if five_prime:
            span = (start, start + take - 1) if strand == "+" else (end - take + 1, end)
        else:
            span = (end - take + 1, end) if strand == "+" else (start, start + take - 1)
        spans.append(span)
        remaining -= take
        if remaining == 0:
            break
    return spans if remaining == 0 else []


def _sync_terminal_feature(
    rows: list[list[str]],
    *,
    cdss: list[list[str]],
    strand: str,
    five_prime: bool,
    removed_row_ids: set[int],
) -> int:
    if not rows:
        return 0
    spans = _terminal_spans(cdss, strand=strand, five_prime=five_prime)
    if not spans:
        removed_row_ids.update(id(row) for row in rows)
        return len(rows)

    ordered_rows = _transcript_order(rows, strand)
    if len(ordered_rows) < len(spans):
        # A split codon cannot be represented correctly by the available rows.
        # The rows are optional in GFF3 and CDS is authoritative, so removing
        # them is safer than retaining an incomplete or intron-spanning codon.
        removed_row_ids.update(id(row) for row in ordered_rows)
        return len(ordered_rows)
    changed = 0
    for row, (new_start, new_end) in zip(ordered_rows, spans, strict=False):
        if int(row[3]) != new_start or int(row[4]) != new_end:
            changed += 1
        row[3] = str(new_start)
        row[4] = str(new_end)
        row[7] = "0"
    for row in ordered_rows[len(spans) :]:
        removed_row_ids.add(id(row))
        changed += 1
    return changed


def _sync_introns(
    introns: list[list[str]],
    *,
    exons: list[list[str]],
    strand: str,
    removed_row_ids: set[int],
) -> int:
    if not introns:
        return 0
    genomic_exons = sorted(exons, key=lambda row: (int(row[3]), int(row[4])))
    spans = [
        (int(left[4]) + 1, int(right[3]) - 1)
        for left, right in zip(genomic_exons, genomic_exons[1:], strict=False)
        if int(left[4]) + 1 <= int(right[3]) - 1
    ]
    ordered_introns = _transcript_order(introns, strand)
    ordered_spans = list(reversed(spans)) if strand == "-" else spans
    changed = 0
    for row, (new_start, new_end) in zip(ordered_introns, ordered_spans, strict=False):
        if int(row[3]) != new_start or int(row[4]) != new_end:
            changed += 1
        row[3] = str(new_start)
        row[4] = str(new_end)
    for row in ordered_introns[len(ordered_spans) :]:
        removed_row_ids.add(id(row))
        changed += 1
    return changed


def synchronize_transcript_children(
    *,
    transcript_row: list[str],
    child_rows: list[list[str]],
    removed_row_ids: set[int],
) -> int:
    """Synchronize coordinate-dependent child rows after CDS/exon trimming."""

    strand = transcript_row[6]
    active = [row for row in child_rows if id(row) not in removed_row_ids]
    exons = [row for row in active if row[2] in EXON_TYPES]
    cdss = [row for row in active if row[2] in CDS_TYPES]
    if not cdss:
        return 0

    changed = 0
    starts = [int(row[3]) for row in (exons or cdss)]
    ends = [int(row[4]) for row in (exons or cdss)]
    transcript_start = min(starts)
    transcript_end = max(ends)
    if int(transcript_row[3]) != transcript_start or int(transcript_row[4]) != transcript_end:
        changed += 1
    transcript_row[3] = str(transcript_start)
    transcript_row[4] = str(transcript_end)

    for row in active:
        if row[2] in EXON_TYPES | CDS_TYPES | CODON_TYPES | INTRON_TYPES:
            continue
        old_start = int(row[3])
        old_end = int(row[4])
        new_start = max(old_start, transcript_start)
        new_end = min(old_end, transcript_end)
        if new_start > new_end:
            removed_row_ids.add(id(row))
            changed += 1
            continue
        if new_start != old_start or new_end != old_end:
            row[3] = str(new_start)
            row[4] = str(new_end)
            changed += 1

    changed += _sync_terminal_feature(
        [row for row in active if row[2] in START_CODON_TYPES],
        cdss=cdss,
        strand=strand,
        five_prime=True,
        removed_row_ids=removed_row_ids,
    )
    changed += _sync_terminal_feature(
        [row for row in active if row[2] in STOP_CODON_TYPES],
        cdss=cdss,
        strand=strand,
        five_prime=False,
        removed_row_ids=removed_row_ids,
    )
    changed += _sync_introns(
        [row for row in active if row[2] in INTRON_TYPES],
        exons=exons,
        strand=strand,
        removed_row_ids=removed_row_ids,
    )
    return changed


def validate_parent_child_containment(
    rows: Iterable[list[str]],
    *,
    scope_parent_ids: set[str] | None = None,
    removed_row_ids: set[int] | None = None,
) -> list[HierarchyIssue]:
    """Return hierarchy problems without rejecting unrelated orphan features."""

    removed = removed_row_ids or set()
    active = [row for row in rows if id(row) not in removed]
    by_id: dict[str, list[list[str]]] = defaultdict(list)
    for row in active:
        feature_id = _row_id(row)
        if feature_id:
            by_id[feature_id].append(row)

    issues: list[HierarchyIssue] = []
    for child in active:
        child_id = _row_id(child) or f"{child[2]}:{child[0]}:{child[3]}-{child[4]}"
        for parent_id in _row_parents(child):
            if scope_parent_ids is not None and parent_id not in scope_parent_ids:
                continue
            parents = by_id.get(parent_id)
            if not parents:
                continue
            if any(
                child[0] == parent[0]
                and int(parent[3]) <= int(child[3])
                and int(child[4]) <= int(parent[4])
                for parent in parents
            ):
                continue
            issues.append(
                HierarchyIssue(
                    child_id=child_id,
                    parent_id=parent_id,
                    message=(
                        f"{child_id} ({child[0]}:{child[3]}-{child[4]}) is outside "
                        f"parent {parent_id}"
                    ),
                )
            )

    children_by_parent: dict[str, list[list[str]]] = defaultdict(list)
    for child in active:
        for parent_id in _row_parents(child):
            children_by_parent[parent_id].append(child)
    for parent_id, children in children_by_parent.items():
        if scope_parent_ids is not None and parent_id not in scope_parent_ids:
            continue
        cdss = [row for row in children if row[2] in CDS_TYPES]
        if not cdss:
            continue
        for codon_type in START_CODON_TYPES | STOP_CODON_TYPES:
            codons = [row for row in children if row[2] == codon_type]
            if not codons:
                continue
            total_length = sum(int(row[4]) - int(row[3]) + 1 for row in codons)
            if total_length != 3:
                issues.append(
                    HierarchyIssue(
                        child_id=codon_type,
                        parent_id=parent_id,
                        message=(
                            f"{codon_type} rows under {parent_id} span {total_length} bases, "
                            "expected 3"
                        ),
                    )
                )
            for codon in codons:
                if not any(
                    int(cds[3]) <= int(codon[3]) and int(codon[4]) <= int(cds[4])
                    for cds in cdss
                ):
                    issues.append(
                        HierarchyIssue(
                            child_id=_row_id(codon) or codon_type,
                            parent_id=parent_id,
                            message=f"{codon_type} under {parent_id} is outside its CDS rows",
                        )
                    )

        for utr in (row for row in children if row[2] in UTR_TYPES):
            if any(
                max(int(utr[3]), int(cds[3])) <= min(int(utr[4]), int(cds[4]))
                for cds in cdss
            ):
                issues.append(
                    HierarchyIssue(
                        child_id=_row_id(utr) or utr[2],
                        parent_id=parent_id,
                        message=f"{utr[2]} under {parent_id} overlaps a CDS row",
                    )
                )
    return issues
