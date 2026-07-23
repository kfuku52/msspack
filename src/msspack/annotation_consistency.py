from __future__ import annotations

import csv
import hashlib
import html
import json
import math
import re
import subprocess
import tempfile
import unicodedata
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

from .config_models import (
    FunctionalAnnotationConfig,
    FunctionalAnnotationConsistencyConfig,
)
from .fasta import iter_fasta
from .step_logging import write_step_log, write_step_metrics
from .utils import MSSPackError, atomic_text_writer, run_command, write_text

SIMILARITY_FIELDS = (
    "qseqid",
    "sseqid",
    "pident",
    "qcovhsp",
    "scovhsp",
    "evalue",
    "bitscore",
    "qlen",
    "slen",
)

GENE_AUDIT_FIELDS = (
    "ID",
    "Locus_tag",
    "source",
    "specific_product",
    "resolved_product",
    "near_identical_family_id",
    "family_id",
    "broad_family_id",
    "family_size",
    "family_product",
    "name_consistency",
    "resolution_action",
    "conflict_severity",
    "reason_code",
    "conflict_partner_ids",
)

FAMILY_AUDIT_FIELDS = (
    "tier",
    "family_id",
    "threshold",
    "family_size",
    "annotated_members",
    "unannotated_members",
    "source_count",
    "direct_annotated_pairs",
    "exact_pairs",
    "canonical_equivalent_pairs",
    "compatible_granularity_pairs",
    "review_pairs",
    "minimum_name_similarity",
    "status",
    "family_product",
    "member_ids",
)

PAIR_AUDIT_FIELDS = (
    "left_id",
    "right_id",
    "left_locus_tag",
    "right_locus_tag",
    "left_source",
    "right_source",
    "left_product",
    "right_product",
    "identity",
    "query_coverage",
    "subject_coverage",
    "evalue",
    "bitscore",
    "near_identical_edge",
    "family_edge",
    "broad_edge",
    "name_relation",
    "name_similarity",
    "reason_code",
)

REVIEW_FIELDS = (
    "tier",
    "severity",
    *PAIR_AUDIT_FIELDS,
)

SUMMARY_FIELDS = (
    "tier",
    "threshold",
    "total_clusters",
    "multi_gene_clusters",
    "genes_in_multi_gene_clusters",
    "evaluable_clusters",
    "exact_clusters",
    "canonical_equivalent_clusters",
    "compatible_granularity_clusters",
    "review_clusters",
    "not_evaluable_clusters",
    "all_unannotated_clusters",
    "one_annotated_clusters",
    "strict_equivalence_rate",
    "biological_compatibility_rate",
    "annotated_direct_pairs",
    "review_direct_pairs",
)

SOURCE_PAIR_FIELDS = (
    "tier",
    "source_a",
    "source_b",
    "annotated_direct_pairs",
    "exact_pairs",
    "canonical_equivalent_pairs",
    "compatible_granularity_pairs",
    "review_pairs",
    "review_rate",
    "meets_minimum_pairs",
)

_GENERIC_NAME_TOKENS = {
    "a",
    "and",
    "containing",
    "domain",
    "family",
    "homolog",
    "like",
    "of",
    "precursor",
    "probable",
    "protein",
    "putative",
    "related",
    "subfamily",
    "subunit",
    "superfamily",
    "the",
    "type",
}

_LOCALIZATION_TOKENS = {
    "chloroplastic",
    "cytosolic",
    "extracellular",
    "mitochondrial",
    "nuclear",
    "peroxisomal",
    "secreted",
}

_PHRASE_ALIASES = (
    (re.compile(r"\bbtb poz\b"), "btb"),
    (re.compile(r"\bpoz btb\b"), "btb"),
    (re.compile(r"\bcchc type zinc finger\b"), "cchc"),
    (re.compile(r"\bzinc knuckle\b"), "cchc"),
    (re.compile(r"\brna recognition motif\b"), "rrm"),
    (re.compile(r"\brna binding domain rnp 1\b"), "rrm"),
)

_SOURCE_ORDER = {
    "existing": 0,
    "reference": 1,
    "swissprot": 2,
    "uniref90": 3,
    "pfam": 4,
    "cdd": 5,
    "none": 9,
}

_CONFIDENCE_ORDER = {"high": 0, "moderate": 1, "low": 2, "none": 9, "": 9}


@dataclass(frozen=True)
class SimilarityEdge:
    left_id: str
    right_id: str
    identity: float
    query_coverage: float
    subject_coverage: float
    evalue: float
    bitscore: float
    query_length: int
    subject_length: int


@dataclass(frozen=True)
class TierSpec:
    key: str
    prefix: str
    identity: float
    coverage: float

    @property
    def threshold_label(self) -> str:
        return f"identity>={self.identity:g}%; mutual coverage>={self.coverage:g}%"


@dataclass(frozen=True)
class NameDecision:
    relation: str
    similarity: float | None
    reason_code: str


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _diamond_version(command: str) -> str:
    try:
        completed = subprocess.run(
            [command, "version"],
            check=False,
            capture_output=True,
            text=True,
            timeout=15,
        )
    except (OSError, subprocess.SubprocessError):
        return "unknown"
    text = (completed.stdout or completed.stderr).strip()
    return text.splitlines()[0] if text else "unknown"


def _tier_specs(config: FunctionalAnnotationConsistencyConfig) -> tuple[TierSpec, ...]:
    return (
        TierSpec(
            "near_identical",
            "NID",
            config.near_identical_identity,
            config.near_identical_coverage,
        ),
        TierSpec("family", "FAM", config.family_identity, config.family_coverage),
        TierSpec("broad", "BRD", config.broad_identity, config.broad_coverage),
    )


def _edge_passes(edge: SimilarityEdge, tier: TierSpec) -> bool:
    return (
        edge.identity >= tier.identity
        and edge.query_coverage >= tier.coverage
        and edge.subject_coverage >= tier.coverage
    )


def run_family_similarity_search(
    *,
    protein_fasta_path: Path,
    output_path: Path,
    provenance_path: Path,
    log_path: Path,
    metrics_path: Path,
    config: FunctionalAnnotationConfig,
) -> None:
    started_at = datetime.now()
    consistency = config.consistency
    protein_count = sum(1 for _ in iter_fasta(protein_fasta_path))
    if protein_count == 0:
        write_text(output_path, "\t".join(SIMILARITY_FIELDS) + "\n")
        write_text(provenance_path, '{"status": "skipped", "reason": "no proteins"}\n')
        write_step_log(
            log_path=log_path,
            command="msspack functional-annotation family-search --skipped",
            step="functional-annotation-family-search",
            started_at=started_at,
            status="skipped",
            count_unit="directed alignments",
            input_total=0,
            changed_total=0,
            output_total=0,
        )
        write_step_metrics(
            metrics_path=metrics_path,
            step="functional-annotation-family-search",
            status="skipped",
            count_unit="directed alignments",
            input_total=0,
            changed_total=0,
            output_total=0,
        )
        return

    with tempfile.TemporaryDirectory(prefix="msspack-family-search-") as temporary_name:
        temporary_root = Path(temporary_name)
        database_path = temporary_root / "proteins"
        raw_path = temporary_root / "alignments.tsv"
        run_command(
            [
                config.diamond_command,
                "makedb",
                "--in",
                str(protein_fasta_path),
                "--db",
                str(database_path),
            ],
            log_path=temporary_root / "makedb.log",
        )
        command = [
            config.diamond_command,
            "blastp",
            "--query",
            str(protein_fasta_path),
            "--db",
            str(database_path),
            "--out",
            str(raw_path),
            "--outfmt",
            "6",
            *SIMILARITY_FIELDS,
            "--threads",
            str(config.threads),
            "--evalue",
            f"{consistency.evalue:g}",
            "--max-target-seqs",
            "0",
            "--id",
            f"{consistency.broad_identity:g}",
            "--query-cover",
            f"{consistency.broad_coverage:g}",
            "--subject-cover",
            f"{consistency.broad_coverage:g}",
        ]
        if config.sensitivity != "default":
            command.append(f"--{config.sensitivity}")
        run_command(command, log_path=temporary_root / "blastp.log")
        if not raw_path.is_file():
            raise MSSPackError(f"DIAMOND did not create its result file: {raw_path}")
        alignment_count = 0
        query_ids: set[str] = set()
        with (
            raw_path.open("r", encoding="utf-8") as source,
            atomic_text_writer(output_path) as output,
        ):
            output.write("\t".join(SIMILARITY_FIELDS) + "\n")
            for raw_line in source:
                line = raw_line.rstrip("\r\n")
                if not line:
                    continue
                fields = line.split("\t")
                if len(fields) != len(SIMILARITY_FIELDS):
                    raise MSSPackError(
                        f"Unexpected DIAMOND family-search output with {len(fields)} columns"
                    )
                output.write(line + "\n")
                alignment_count += 1
                if fields[0] != fields[1]:
                    query_ids.add(fields[0])

    write_text(
        provenance_path,
        json.dumps(
            {
                "diamond_command": config.diamond_command,
                "diamond_version": _diamond_version(config.diamond_command),
                "protein_fasta": str(protein_fasta_path),
                "protein_fasta_sha256": _sha256(protein_fasta_path),
                "fields": list(SIMILARITY_FIELDS),
                "max_target_seqs": 0,
                "sensitivity": config.sensitivity,
                "evalue": consistency.evalue,
                "broad_identity": consistency.broad_identity,
                "broad_coverage": consistency.broad_coverage,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
    )
    write_step_log(
        log_path=log_path,
        command=command,
        step="functional-annotation-family-search",
        started_at=started_at,
        count_unit="directed alignments",
        input_total=protein_count,
        changed_total=len(query_ids),
        output_total=alignment_count,
        details=[
            "One all-vs-all search was run at the broadest configured threshold.",
            f"Queries with non-self hits: {len(query_ids):,}",
            "Self hits are retained in the raw table and removed during the audit.",
        ],
    )
    write_step_metrics(
        metrics_path=metrics_path,
        step="functional-annotation-family-search",
        count_unit="directed alignments",
        input_total=protein_count,
        changed_total=len(query_ids),
        output_total=alignment_count,
        details={
            "queries_with_non_self_hits": len(query_ids),
            "max_target_seqs": 0,
            "evalue": consistency.evalue,
            "broad_identity": consistency.broad_identity,
            "broad_coverage": consistency.broad_coverage,
        },
    )


def _surface_normalize(value: str) -> str:
    value = html.unescape(unicodedata.normalize("NFKC", value)).casefold()
    value = re.sub(r"[^a-z0-9]+", " ", value)
    return " ".join(value.split())


def _canonical_tokens(value: str) -> frozenset[str]:
    normalized = _surface_normalize(value)
    for pattern, replacement in _PHRASE_ALIASES:
        normalized = pattern.sub(replacement, normalized)
    return frozenset(token for token in normalized.split() if token not in _GENERIC_NAME_TOKENS)


def _weighted_jaccard(
    left: frozenset[str],
    right: frozenset[str],
    idf: dict[str, float],
) -> float:
    union = left | right
    if not union:
        return 0.0
    return sum(idf.get(token, 1.0) for token in left & right) / sum(
        idf.get(token, 1.0) for token in union
    )


def _specific_modifiers(tokens: frozenset[str]) -> frozenset[str]:
    return frozenset(
        token
        for token in tokens
        if token in _LOCALIZATION_TOKENS
        or token.isdigit()
        or re.fullmatch(r"[a-z]+\d+[a-z]*", token)
        or re.fullmatch(r"[a-z]\d+", token)
    )


def classify_product_names(
    left_product: str,
    right_product: str,
    *,
    idf: dict[str, float] | None = None,
    similarity_threshold: float = 0.45,
) -> NameDecision:
    if _surface_normalize(left_product) == _surface_normalize(right_product):
        return NameDecision("exact", 1.0, "normalized_text_match")
    left_tokens = _canonical_tokens(left_product)
    right_tokens = _canonical_tokens(right_product)
    if left_tokens and left_tokens == right_tokens:
        return NameDecision(
            "canonical_equivalent",
            1.0,
            "approved_alias_or_safe_formatting_difference",
        )
    weights = idf or {}
    similarity = _weighted_jaccard(left_tokens, right_tokens, weights)
    left_modifiers = _specific_modifiers(left_tokens)
    right_modifiers = _specific_modifiers(right_tokens)
    if left_modifiers and right_modifiers and left_modifiers != right_modifiers:
        return NameDecision("review", similarity, "specific_modifier_conflict")
    shared_tokens = left_tokens & right_tokens
    if shared_tokens and similarity >= similarity_threshold:
        return NameDecision(
            "compatible_granularity",
            similarity,
            "shared_family_name_with_granularity_difference",
        )
    return NameDecision("review", similarity, "incompatible_family_name_tokens")


def _read_evidence(path: Path) -> dict[str, dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"ID", "Locus_tag", "assigned_product", "source"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise MSSPackError(
                "Functional annotation evidence must contain ID, Locus_tag, "
                f"assigned_product, and source columns: {path}"
            )
        rows = list(reader)
    evidence = {row["ID"]: row for row in rows}
    if len(evidence) != len(rows):
        raise MSSPackError(f"Functional annotation evidence contains duplicate IDs: {path}")
    return evidence


def _read_similarity_edges(
    path: Path,
    identifiers: set[str],
    *,
    evalue: float,
) -> list[SimilarityEdge]:
    deduplicated: dict[tuple[str, str], SimilarityEdge] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != list(SIMILARITY_FIELDS):
            raise MSSPackError(f"Family similarity table has an unexpected header: {path}")
        for row in reader:
            left = row["qseqid"]
            right = row["sseqid"]
            if left == right:
                continue
            if left not in identifiers or right not in identifiers:
                raise MSSPackError(
                    "Family similarity table references an ID absent from functional "
                    f"annotation evidence: {left if left not in identifiers else right}"
                )
            try:
                edge = SimilarityEdge(
                    left_id=min(left, right),
                    right_id=max(left, right),
                    identity=float(row["pident"]),
                    query_coverage=float(row["qcovhsp"]),
                    subject_coverage=float(row["scovhsp"]),
                    evalue=float(row["evalue"]),
                    bitscore=float(row["bitscore"]),
                    query_length=int(row["qlen"]),
                    subject_length=int(row["slen"]),
                )
            except ValueError as exc:
                raise MSSPackError(f"Could not parse family similarity row: {row}") from exc
            if edge.evalue > evalue:
                continue
            key = (edge.left_id, edge.right_id)
            previous = deduplicated.get(key)
            if previous is None or edge.bitscore > previous.bitscore:
                deduplicated[key] = edge
    return sorted(deduplicated.values(), key=lambda edge: (edge.left_id, edge.right_id))


def _family_identifier(prefix: str, members: tuple[str, ...]) -> str:
    digest = hashlib.sha256("\0".join(members).encode("utf-8")).hexdigest()[:12]
    return f"{prefix}_{digest}"


def _components(
    identifiers: set[str],
    edges: list[SimilarityEdge],
    tier: TierSpec,
) -> tuple[dict[str, str], dict[str, tuple[str, ...]], list[SimilarityEdge]]:
    parent = {identifier: identifier for identifier in identifiers}

    def find(identifier: str) -> str:
        while parent[identifier] != identifier:
            parent[identifier] = parent[parent[identifier]]
            identifier = parent[identifier]
        return identifier

    def union(left: str, right: str) -> None:
        left_root = find(left)
        right_root = find(right)
        if left_root != right_root:
            parent[max(left_root, right_root)] = min(left_root, right_root)

    accepted = [edge for edge in edges if _edge_passes(edge, tier)]
    for edge in accepted:
        union(edge.left_id, edge.right_id)
    grouped: dict[str, list[str]] = defaultdict(list)
    for identifier in sorted(identifiers):
        grouped[find(identifier)].append(identifier)
    id_to_family: dict[str, str] = {}
    families: dict[str, tuple[str, ...]] = {}
    for members_list in grouped.values():
        members = tuple(sorted(members_list))
        family_id = _family_identifier(tier.prefix, members)
        families[family_id] = members
        for member in members:
            id_to_family[member] = family_id
    return id_to_family, families, accepted


def _build_idf(evidence: dict[str, dict[str, str]]) -> dict[str, float]:
    document_frequency: Counter[str] = Counter()
    products = [
        row["assigned_product"] for row in evidence.values() if row["source"].casefold() != "none"
    ]
    for product in products:
        document_frequency.update(_canonical_tokens(product))
    total = len(products)
    return {
        token: math.log((total + 1) / (frequency + 1)) + 1.0
        for token, frequency in document_frequency.items()
    }


def _representative_product(
    members: tuple[str, ...],
    evidence: dict[str, dict[str, str]],
) -> str:
    annotated = [
        evidence[member] for member in members if evidence[member]["source"].casefold() != "none"
    ]
    if not annotated:
        return ""

    def numeric(row: dict[str, str], key: str) -> float:
        try:
            return float(row.get(key, ""))
        except ValueError:
            return 0.0

    selected = min(
        annotated,
        key=lambda row: (
            _SOURCE_ORDER.get(row["source"].casefold(), 6),
            _CONFIDENCE_ORDER.get(row.get("confidence", "").casefold(), 8),
            -numeric(row, "consensus_score"),
            -numeric(row, "bitscore"),
            row["assigned_product"].casefold(),
            row["ID"],
        ),
    )
    return selected["assigned_product"]


def _auto_resolved_near_identical_products(
    families: dict[str, tuple[str, ...]],
    evidence: dict[str, dict[str, str]],
    review_partners: dict[str, set[str]],
) -> dict[str, str]:
    """Propagate only a unique better-source product across a direct 90/90 conflict."""
    resolved: dict[str, str] = {}
    for members in families.values():
        annotated = [
            member for member in members if evidence[member]["source"].casefold() != "none"
        ]
        if not any(review_partners.get(member) for member in annotated):
            continue
        best_priority = min(
            _SOURCE_ORDER.get(evidence[member]["source"].casefold(), 6)
            for member in annotated
        )
        best_members = tuple(
            member
            for member in annotated
            if _SOURCE_ORDER.get(evidence[member]["source"].casefold(), 6) == best_priority
        )
        best_products = {
            _surface_normalize(evidence[member]["assigned_product"])
            for member in best_members
        }
        if len(best_products) != 1:
            continue
        representative = _representative_product(best_members, evidence)
        for member in annotated:
            member_priority = _SOURCE_ORDER.get(evidence[member]["source"].casefold(), 6)
            if member_priority <= best_priority:
                continue
            if not set(best_members) & review_partners.get(member, set()):
                continue
            if _surface_normalize(evidence[member]["assigned_product"]) == _surface_normalize(
                representative
            ):
                continue
            resolved[member] = representative
    return resolved


def _format_float(value: float | None) -> str:
    return "" if value is None else f"{value:.6g}"


def _object_int(value: object) -> int:
    if isinstance(value, bool):
        raise TypeError("Boolean is not an integer count")
    if isinstance(value, int | float | str):
        return int(value)
    raise TypeError(f"Expected an integer-compatible value, got {type(value).__name__}")


def _object_float(value: object) -> float:
    if isinstance(value, bool):
        raise TypeError("Boolean is not a floating-point measurement")
    if isinstance(value, int | float | str):
        return float(value)
    raise TypeError(f"Expected a float-compatible value, got {type(value).__name__}")


def _write_tsv(path: Path, fieldnames: tuple[str, ...], rows: list[dict[str, object]]) -> None:
    with atomic_text_writer(path) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def audit_annotation_consistency(
    *,
    annotation_table_path: Path,
    evidence_path: Path,
    similarity_path: Path,
    output_annotation_table_path: Path,
    gene_output_path: Path,
    family_output_path: Path,
    pair_output_path: Path,
    review_output_path: Path,
    summary_output_path: Path,
    source_pair_output_path: Path,
    log_path: Path,
    metrics_path: Path,
    config: FunctionalAnnotationConfig,
) -> None:
    started_at = datetime.now()
    consistency = config.consistency
    evidence = _read_evidence(evidence_path)
    identifiers = set(evidence)
    edges = _read_similarity_edges(
        similarity_path,
        identifiers,
        evalue=consistency.evalue,
    )
    idf = _build_idf(evidence)
    decisions = {
        (edge.left_id, edge.right_id): (
            NameDecision("not_evaluable", None, "unannotated_pair")
            if evidence[edge.left_id]["source"].casefold() == "none"
            or evidence[edge.right_id]["source"].casefold() == "none"
            else classify_product_names(
                evidence[edge.left_id]["assigned_product"],
                evidence[edge.right_id]["assigned_product"],
                idf=idf,
                similarity_threshold=consistency.name_similarity_threshold,
            )
        )
        for edge in edges
    }

    tier_results: dict[
        str,
        tuple[dict[str, str], dict[str, tuple[str, ...]], list[SimilarityEdge]],
    ] = {}
    family_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    family_row_by_id: dict[str, dict[str, object]] = {}

    for tier in _tier_specs(consistency):
        id_to_family, families, tier_edges = _components(identifiers, edges, tier)
        tier_results[tier.key] = (id_to_family, families, tier_edges)
        edges_by_family: dict[str, list[SimilarityEdge]] = defaultdict(list)
        for edge in tier_edges:
            edges_by_family[id_to_family[edge.left_id]].append(edge)
        status_counts: Counter[str] = Counter()
        annotated_direct_pairs = 0
        review_direct_pairs = 0
        genes_in_multi = 0
        multi_gene_clusters = 0
        for family_id, members in sorted(families.items()):
            if len(members) > 1:
                multi_gene_clusters += 1
                genes_in_multi += len(members)
            annotated_members = [
                member for member in members if evidence[member]["source"].casefold() != "none"
            ]
            family_decisions = [
                decisions[(edge.left_id, edge.right_id)]
                for edge in edges_by_family.get(family_id, [])
                if decisions[(edge.left_id, edge.right_id)].relation != "not_evaluable"
            ]
            counts = Counter(decision.relation for decision in family_decisions)
            annotated_direct_pairs += len(family_decisions)
            review_direct_pairs += counts["review"]
            similarities = [
                decision.similarity
                for decision in family_decisions
                if decision.similarity is not None
            ]
            if not annotated_members:
                status = "all_unannotated"
            elif len(annotated_members) == 1:
                status = "one_annotated"
            elif not family_decisions:
                status = "not_evaluable"
            elif counts["review"]:
                status = "review"
            elif counts["compatible_granularity"]:
                status = "compatible_granularity"
            elif counts["canonical_equivalent"]:
                status = "canonical_equivalent"
            else:
                status = "exact"
            status_counts[status] += 1
            family_product = (
                _representative_product(members, evidence)
                if status in {"exact", "canonical_equivalent"}
                else ""
            )
            source_count = len({evidence[member]["source"] for member in annotated_members})
            row: dict[str, object] = {
                "tier": tier.key,
                "family_id": family_id,
                "threshold": tier.threshold_label,
                "family_size": len(members),
                "annotated_members": len(annotated_members),
                "unannotated_members": len(members) - len(annotated_members),
                "source_count": source_count,
                "direct_annotated_pairs": len(family_decisions),
                "exact_pairs": counts["exact"],
                "canonical_equivalent_pairs": counts["canonical_equivalent"],
                "compatible_granularity_pairs": counts["compatible_granularity"],
                "review_pairs": counts["review"],
                "minimum_name_similarity": (
                    _format_float(min(similarities)) if similarities else ""
                ),
                "status": status,
                "family_product": family_product,
                "member_ids": ",".join(members),
            }
            family_rows.append(row)
            family_row_by_id[family_id] = row

        evaluable = sum(
            status_counts[key]
            for key in ("exact", "canonical_equivalent", "compatible_granularity", "review")
        )
        strict = status_counts["exact"] + status_counts["canonical_equivalent"]
        compatible = strict + status_counts["compatible_granularity"]
        summary_rows.append(
            {
                "tier": tier.key,
                "threshold": tier.threshold_label,
                "total_clusters": len(families),
                "multi_gene_clusters": multi_gene_clusters,
                "genes_in_multi_gene_clusters": genes_in_multi,
                "evaluable_clusters": evaluable,
                "exact_clusters": status_counts["exact"],
                "canonical_equivalent_clusters": status_counts["canonical_equivalent"],
                "compatible_granularity_clusters": status_counts["compatible_granularity"],
                "review_clusters": status_counts["review"],
                "not_evaluable_clusters": status_counts["not_evaluable"],
                "all_unannotated_clusters": status_counts["all_unannotated"],
                "one_annotated_clusters": status_counts["one_annotated"],
                "strict_equivalence_rate": _format_float(strict / evaluable if evaluable else 0.0),
                "biological_compatibility_rate": _format_float(
                    compatible / evaluable if evaluable else 0.0
                ),
                "annotated_direct_pairs": annotated_direct_pairs,
                "review_direct_pairs": review_direct_pairs,
            }
        )

    pair_rows: list[dict[str, object]] = []
    review_rows: list[dict[str, object]] = []
    source_pair_counts: Counter[tuple[str, str, str]] = Counter()
    review_partners: dict[str, dict[str, set[str]]] = {
        tier: defaultdict(set) for tier in tier_results
    }
    review_reasons: dict[str, dict[str, set[str]]] = {
        tier: defaultdict(set) for tier in tier_results
    }
    comparable_by_tier: dict[str, set[str]] = {tier: set() for tier in tier_results}

    for edge in edges:
        key = (edge.left_id, edge.right_id)
        decision = decisions[key]
        qualifies = {tier.key: _edge_passes(edge, tier) for tier in _tier_specs(consistency)}
        if not qualifies["broad"]:
            continue
        left = evidence[edge.left_id]
        right = evidence[edge.right_id]
        row = {
            "left_id": edge.left_id,
            "right_id": edge.right_id,
            "left_locus_tag": left["Locus_tag"],
            "right_locus_tag": right["Locus_tag"],
            "left_source": left["source"],
            "right_source": right["source"],
            "left_product": left["assigned_product"],
            "right_product": right["assigned_product"],
            "identity": _format_float(edge.identity),
            "query_coverage": _format_float(edge.query_coverage),
            "subject_coverage": _format_float(edge.subject_coverage),
            "evalue": _format_float(edge.evalue),
            "bitscore": _format_float(edge.bitscore),
            "near_identical_edge": str(qualifies["near_identical"]).lower(),
            "family_edge": str(qualifies["family"]).lower(),
            "broad_edge": "true",
            "name_relation": decision.relation,
            "name_similarity": _format_float(decision.similarity),
            "reason_code": decision.reason_code,
        }
        pair_rows.append(row)
        highest_tier = (
            "near_identical"
            if qualifies["near_identical"]
            else "family"
            if qualifies["family"]
            else "broad"
        )
        if decision.relation == "review":
            severity = {"near_identical": "high", "family": "medium", "broad": "low"}[highest_tier]
            review_rows.append({"tier": highest_tier, "severity": severity, **row})
            for tier_key, accepted in qualifies.items():
                if not accepted:
                    continue
                review_partners[tier_key][edge.left_id].add(edge.right_id)
                review_partners[tier_key][edge.right_id].add(edge.left_id)
                review_reasons[tier_key][edge.left_id].add(decision.reason_code)
                review_reasons[tier_key][edge.right_id].add(decision.reason_code)
        for tier_key, accepted in qualifies.items():
            if accepted and decision.relation != "not_evaluable":
                comparable_by_tier[tier_key].update((edge.left_id, edge.right_id))
        if qualifies["near_identical"] and decision.relation != "not_evaluable":
            source_a, source_b = sorted(
                (left["source"], right["source"]),
                key=lambda source: (_SOURCE_ORDER.get(source.casefold(), 6), source.casefold()),
            )
            source_pair_counts[(source_a, source_b, decision.relation)] += 1

    source_pair_rows: list[dict[str, object]] = []
    source_pairs = sorted({(left, right) for left, right, _relation in source_pair_counts})
    for source_a, source_b in source_pairs:
        counts = Counter(
            {
                relation: source_pair_counts[(source_a, source_b, relation)]
                for relation in (
                    "exact",
                    "canonical_equivalent",
                    "compatible_granularity",
                    "review",
                )
            }
        )
        total = sum(counts.values())
        source_pair_rows.append(
            {
                "tier": "near_identical",
                "source_a": source_a,
                "source_b": source_b,
                "annotated_direct_pairs": total,
                "exact_pairs": counts["exact"],
                "canonical_equivalent_pairs": counts["canonical_equivalent"],
                "compatible_granularity_pairs": counts["compatible_granularity"],
                "review_pairs": counts["review"],
                "review_rate": _format_float(counts["review"] / total if total else 0.0),
                "meets_minimum_pairs": str(total >= consistency.source_pair_min_pairs).lower(),
            }
        )

    near_ids, near_families, _near_edges = tier_results["near_identical"]
    family_ids, family_families, _family_edges = tier_results["family"]
    broad_ids, _broad_families, _broad_edges = tier_results["broad"]
    gene_rows: list[dict[str, object]] = []
    resolved_products: dict[str, str] = {}
    if consistency.auto_resolve_conflicts:
        resolved_products.update(
            _auto_resolved_near_identical_products(
                near_families,
                evidence,
                review_partners["near_identical"],
            )
        )
    auto_resolved_product_changes = len(resolved_products)
    for identifier in sorted(identifiers):
        evidence_row = evidence[identifier]
        source_key = evidence_row["source"].casefold()
        near_family_id = near_ids[identifier]
        family_id = family_ids[identifier]
        family_row = family_row_by_id[family_id]
        if source_key == "none":
            name_consistency = "unannotated"
            resolution_action = "none"
            severity = ""
            reason_code = "unannotated"
        elif identifier in review_partners["family"]:
            review_tier = (
                "near_identical"
                if identifier in review_partners["near_identical"]
                else "family"
            )
            severity = "high" if review_tier == "near_identical" else "medium"
            reason_code = ",".join(sorted(review_reasons[review_tier][identifier]))
            if consistency.auto_resolve_conflicts:
                name_consistency = "resolved"
                if identifier in resolved_products:
                    resolution_action = "higher_priority_near_identical_product_propagated"
                elif review_tier == "near_identical":
                    resolution_action = "independently_supported_specific_variant_retained"
                else:
                    resolution_action = "close_family_specific_variation_retained"
            else:
                name_consistency = "review"
                resolution_action = "manual_review"
        elif identifier in comparable_by_tier["family"]:
            name_consistency = "consistent"
            resolution_action = "none"
            severity = ""
            reason_code = "close_family_names_compatible"
        else:
            name_consistency = "no_close_family_peer"
            resolution_action = "none"
            if identifier in review_partners["broad"]:
                severity = "low"
                reason_code = ",".join(sorted(review_reasons["broad"][identifier]))
            else:
                severity = ""
                reason_code = "no_annotated_close_family_peer"
        partner_tier = (
            "near_identical"
            if review_partners["near_identical"].get(identifier)
            else "family"
            if review_partners["family"].get(identifier)
            else "broad"
        )
        family_product = str(family_row["family_product"])
        near_family_row = family_row_by_id[near_family_id]
        if (
            consistency.harmonize_safe_equivalents
            and source_key != "none"
            and near_family_row["status"] in {"exact", "canonical_equivalent"}
            and _object_int(near_family_row["direct_annotated_pairs"]) > 0
            and near_family_row["family_product"]
        ):
            resolved_products[identifier] = str(near_family_row["family_product"])
        gene_rows.append(
            {
                "ID": identifier,
                "Locus_tag": evidence_row["Locus_tag"],
                "source": evidence_row["source"],
                "specific_product": evidence_row["assigned_product"],
                "resolved_product": resolved_products.get(
                    identifier, evidence_row["assigned_product"]
                ),
                "near_identical_family_id": near_family_id,
                "family_id": family_id,
                "broad_family_id": broad_ids[identifier],
                "family_size": len(family_families[family_id]),
                "family_product": family_product,
                "name_consistency": name_consistency,
                "resolution_action": resolution_action,
                "conflict_severity": severity,
                "reason_code": reason_code,
                "conflict_partner_ids": ",".join(
                    sorted(review_partners[partner_tier].get(identifier, set()))
                ),
            }
        )

    changed_products = 0
    with (
        annotation_table_path.open("r", encoding="utf-8", newline="") as source_handle,
        atomic_text_writer(output_annotation_table_path) as output,
    ):
        reader = csv.DictReader(source_handle, delimiter="\t")
        if reader.fieldnames != ["ID", "Description", "Locus_tag"]:
            raise MSSPackError("Annotation table must have columns: ID, Description, Locus_tag")
        writer = csv.DictWriter(
            output,
            fieldnames=reader.fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        seen_annotation_ids: set[str] = set()
        for table_row in reader:
            identifier = table_row["ID"]
            seen_annotation_ids.add(identifier)
            replacement = resolved_products.get(identifier)
            if replacement and replacement != table_row["Description"]:
                table_row["Description"] = replacement
                changed_products += 1
            writer.writerow(table_row)
    if seen_annotation_ids != identifiers:
        raise MSSPackError(
            "Annotation table and functional annotation evidence contain different IDs"
        )

    _write_tsv(gene_output_path, GENE_AUDIT_FIELDS, gene_rows)
    _write_tsv(family_output_path, FAMILY_AUDIT_FIELDS, family_rows)
    _write_tsv(pair_output_path, PAIR_AUDIT_FIELDS, pair_rows)
    _write_tsv(review_output_path, REVIEW_FIELDS, review_rows)
    _write_tsv(summary_output_path, SUMMARY_FIELDS, summary_rows)
    _write_tsv(source_pair_output_path, SOURCE_PAIR_FIELDS, source_pair_rows)

    near_summary = next(row for row in summary_rows if row["tier"] == "near_identical")
    family_spec = next(tier for tier in _tier_specs(consistency) if tier.key == "family")
    write_step_log(
        log_path=log_path,
        command="msspack functional-annotation consistency-audit",
        step="functional-annotation-consistency-audit",
        started_at=started_at,
        count_unit="genes",
        input_total=len(identifiers),
        changed_total=changed_products,
        output_total=len(gene_rows),
        details=[
            f"Deduplicated non-self similarity pairs: {len(edges):,}",
            "Gene-level name consistency: close family "
            f"({family_spec.threshold_label})",
            "Near-identical evaluable families: "
            f"{_object_int(near_summary['evaluable_clusters']):,}",
            f"Near-identical review families: {_object_int(near_summary['review_clusters']):,}",
            f"Name-conflict diagnostic pairs across all tiers: {len(review_rows):,}",
            "Automatic conflict resolution: unique higher-priority products were propagated "
            "only across direct near-identical conflicts; other specific family variants "
            "were retained."
            if consistency.auto_resolve_conflicts
            else "Automatic conflict resolution was disabled.",
            "Product changes were limited to safe canonical-equivalent near-identical families."
            if consistency.harmonize_safe_equivalents
            else "Audit-only mode: no annotation products were changed.",
        ],
    )
    write_step_metrics(
        metrics_path=metrics_path,
        step="functional-annotation-consistency-audit",
        count_unit="genes",
        input_total=len(identifiers),
        changed_total=changed_products,
        output_total=len(gene_rows),
        details={
            "deduplicated_similarity_pairs": len(edges),
            "conflict_diagnostic_pairs": len(review_rows),
            "near_identical_evaluable_families": _object_int(near_summary["evaluable_clusters"]),
            "near_identical_review_families": _object_int(near_summary["review_clusters"]),
            "near_identical_biological_compatibility_rate": _object_float(
                near_summary["biological_compatibility_rate"]
            ),
            "gene_consistency_tier": "family",
            "gene_consistency_identity_threshold": family_spec.identity,
            "gene_consistency_coverage_threshold": family_spec.coverage,
            "auto_resolve_conflicts": consistency.auto_resolve_conflicts,
            "auto_resolved_product_changes": auto_resolved_product_changes,
            "harmonize_safe_equivalents": consistency.harmonize_safe_equivalents,
        },
    )


__all__ = [
    "audit_annotation_consistency",
    "classify_product_names",
    "run_family_similarity_search",
]
