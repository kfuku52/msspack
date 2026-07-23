from __future__ import annotations

import hashlib
import json
import re
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path

from .config_models import FunctionalAnnotationTaxonomyConfig
from .step_logging import write_step_log, write_step_metrics
from .utils import MSSPackError, ensure_dir, write_text

NCBI_TAXONOMY_API = "https://api.ncbi.nlm.nih.gov/datasets/v2/taxonomy/taxon"
TAXONOMY_CACHE_SCHEMA_VERSION = 1
TAXONOMY_CONTEXT_SCHEMA_VERSION = 1
NCBI_DOMAIN_IDS = frozenset({2, 2157, 2759, 10239})

RELATION_WEIGHTS: dict[str, float] = {
    "same_species": 1.35,
    "same_genus": 1.30,
    "same_family": 1.25,
    "same_order": 1.20,
    "same_class": 1.12,
    "same_phylum": 1.05,
    "same_kingdom": 1.00,
    "same_domain": 0.88,
    "cross_kingdom": 0.72,
    "cross_domain": 0.45,
    "unknown": 0.90,
    "reference": 1.00,
}

_COMPARISON_RANKS = (
    ("GENUS", "same_genus"),
    ("FAMILY", "same_family"),
    ("ORDER", "same_order"),
    ("CLASS", "same_class"),
    ("PHYLUM", "same_phylum"),
    ("KINGDOM", "same_kingdom"),
)

_BUSCO_RANK_PREFERENCE = {
    "DOMAIN": 0,
    "KINGDOM": 1,
    "PHYLUM": 2,
    "CLASS": 3,
    "ORDER": 4,
    "CLADE": 5,
    "FAMILY": 6,
    "GENUS": 7,
    "SPECIES": 8,
}


@dataclass(frozen=True)
class TaxonomyRecord:
    tax_id: int
    scientific_name: str
    rank: str
    lineage_ids: tuple[int, ...]
    blast_name: str = ""

    @classmethod
    def from_dict(cls, data: dict[str, object]) -> TaxonomyRecord:
        raw_tax_id = data.get("tax_id")
        if isinstance(raw_tax_id, bool) or not isinstance(raw_tax_id, int | str):
            raise TypeError("Taxonomy record TaxID must be integer-compatible")
        raw_lineage = data.get("lineage_ids", [])
        lineage = (
            tuple(int(value) for value in raw_lineage)
            if isinstance(raw_lineage, list)
            else ()
        )
        return cls(
            tax_id=int(raw_tax_id),
            scientific_name=str(data.get("scientific_name", "")),
            rank=str(data.get("rank", "")).upper(),
            lineage_ids=lineage,
            blast_name=str(data.get("blast_name", "")),
        )

    def to_dict(self) -> dict[str, object]:
        return {
            "tax_id": self.tax_id,
            "scientific_name": self.scientific_name,
            "rank": self.rank,
            "lineage_ids": list(self.lineage_ids),
            "blast_name": self.blast_name,
        }


def _now_utc() -> str:
    return datetime.now(UTC).isoformat()


def _record_cache_path(cache_dir: Path, tax_id: int) -> Path:
    return cache_dir / "records" / f"{tax_id}.json"


def _query_cache_path(cache_dir: Path, query: str) -> Path:
    digest = hashlib.sha256(query.casefold().strip().encode("utf-8")).hexdigest()
    return cache_dir / "queries" / f"{digest}.json"


def _read_json(path: Path) -> dict[str, object]:
    if not path.is_file():
        return {}
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}
    return value if isinstance(value, dict) else {}


def _request_taxonomy_payload(query: str) -> dict[str, object]:
    encoded = urllib.parse.quote(query, safe=",")
    request = urllib.request.Request(
        f"{NCBI_TAXONOMY_API}/{encoded}",
        headers={
            "Accept": "application/json",
            "User-Agent": "msspack-functional-annotation",
        },
    )
    try:
        with urllib.request.urlopen(request, timeout=60) as response:
            payload: object = json.loads(response.read().decode("utf-8"))
    except (OSError, UnicodeError, ValueError, urllib.error.HTTPError) as exc:
        raise MSSPackError(f"Could not resolve NCBI taxonomy query {query!r}: {exc}") from exc
    if not isinstance(payload, dict):
        raise MSSPackError(f"NCBI taxonomy returned a non-object response for {query!r}")
    return payload


def _records_from_payload(payload: dict[str, object]) -> list[TaxonomyRecord]:
    raw_nodes = payload.get("taxonomy_nodes", [])
    if not isinstance(raw_nodes, list):
        return []
    records: list[TaxonomyRecord] = []
    for raw_node in raw_nodes:
        if not isinstance(raw_node, dict):
            continue
        raw_taxonomy = raw_node.get("taxonomy")
        if not isinstance(raw_taxonomy, dict):
            continue
        raw_tax_id = raw_taxonomy.get("tax_id")
        if not isinstance(raw_tax_id, int) or isinstance(raw_tax_id, bool):
            continue
        raw_lineage = raw_taxonomy.get("lineage", [])
        lineage_ids = (
            tuple(
                value
                for value in raw_lineage
                if isinstance(value, int) and not isinstance(value, bool)
            )
            if isinstance(raw_lineage, list)
            else ()
        )
        records.append(
            TaxonomyRecord(
                tax_id=raw_tax_id,
                scientific_name=str(raw_taxonomy.get("organism_name", "")),
                rank=str(raw_taxonomy.get("rank", "")).upper(),
                lineage_ids=lineage_ids,
                blast_name=str(raw_taxonomy.get("blast_name", "")),
            )
        )
    return records


def _cache_record(cache_dir: Path, record: TaxonomyRecord) -> None:
    path = _record_cache_path(cache_dir, record.tax_id)
    ensure_dir(path.parent)
    write_text(
        path,
        json.dumps(
            {
                "schema_version": TAXONOMY_CACHE_SCHEMA_VERSION,
                "retrieved_at": _now_utc(),
                **record.to_dict(),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
    )


def _cached_record(cache_dir: Path, tax_id: int) -> TaxonomyRecord | None:
    data = _read_json(_record_cache_path(cache_dir, tax_id))
    if data.get("schema_version") != TAXONOMY_CACHE_SCHEMA_VERSION:
        return None
    try:
        return TaxonomyRecord.from_dict(data)
    except (KeyError, TypeError, ValueError):
        return None


def resolve_taxonomy_records(
    tax_ids: set[int],
    *,
    cache_dir: Path,
    offline: bool,
    strict: bool = False,
) -> tuple[dict[int, TaxonomyRecord], list[str]]:
    records: dict[int, TaxonomyRecord] = {}
    warnings: list[str] = []
    missing: list[int] = []
    for tax_id in sorted(value for value in tax_ids if value > 0):
        cached = _cached_record(cache_dir, tax_id)
        if cached is None:
            missing.append(tax_id)
        else:
            records[tax_id] = cached
    if offline:
        if missing:
            warnings.append(
                f"{len(missing):,} taxonomy records were absent from the offline cache"
            )
        return records, warnings

    for start in range(0, len(missing), 100):
        batch = missing[start : start + 100]
        try:
            fetched = _records_from_payload(
                _request_taxonomy_payload(",".join(str(value) for value in batch))
            )
        except MSSPackError as exc:
            if strict:
                raise
            warnings.append(str(exc))
            continue
        for record in fetched:
            records[record.tax_id] = record
            _cache_record(cache_dir, record)
        unresolved = sorted(set(batch) - {record.tax_id for record in fetched})
        if unresolved:
            warnings.append(
                "NCBI taxonomy did not return records for TaxID(s): "
                + ", ".join(map(str, unresolved[:10]))
                + (" ..." if len(unresolved) > 10 else "")
            )
    return records, warnings


def _select_named_record(
    query: str,
    records: list[TaxonomyRecord],
    *,
    busco_lineage: bool,
) -> TaxonomyRecord | None:
    exact = [
        record
        for record in records
        if record.scientific_name.casefold() == query.casefold().strip()
    ]
    candidates = exact or records
    if not candidates:
        return None
    if busco_lineage:
        return min(
            candidates,
            key=lambda record: (
                _BUSCO_RANK_PREFERENCE.get(record.rank, 20),
                record.tax_id,
            ),
        )
    species_like = {"SPECIES", "SUBSPECIES", "STRAIN", "VARIETAS", "FORMA"}
    return min(
        candidates,
        key=lambda record: (
            0 if record.rank in species_like else 1,
            record.tax_id,
        ),
    )


def resolve_named_taxonomy(
    query: str,
    *,
    cache_dir: Path,
    offline: bool,
    busco_lineage: bool = False,
    strict: bool = False,
) -> tuple[TaxonomyRecord | None, list[str]]:
    warnings: list[str] = []
    normalized = query.strip()
    if not normalized:
        return None, warnings
    query_cache = _read_json(_query_cache_path(cache_dir, normalized))
    cached_tax_id = query_cache.get("tax_id")
    if isinstance(cached_tax_id, int) and not isinstance(cached_tax_id, bool):
        cached = _cached_record(cache_dir, cached_tax_id)
        if cached is not None:
            return cached, warnings
    if offline:
        warnings.append(f"Taxonomy query {normalized!r} was absent from the offline cache")
        return None, warnings
    try:
        records = _records_from_payload(_request_taxonomy_payload(normalized))
    except MSSPackError as exc:
        if strict:
            raise
        return None, [str(exc)]
    selected = _select_named_record(
        normalized,
        records,
        busco_lineage=busco_lineage,
    )
    for record in records:
        _cache_record(cache_dir, record)
    if selected is None:
        message = f"NCBI taxonomy did not resolve {normalized!r}"
        if strict:
            raise MSSPackError(message)
        return None, [message]
    query_path = _query_cache_path(cache_dir, normalized)
    ensure_dir(query_path.parent)
    write_text(
        query_path,
        json.dumps(
            {
                "schema_version": TAXONOMY_CACHE_SCHEMA_VERSION,
                "query": normalized,
                "tax_id": selected.tax_id,
                "retrieved_at": _now_utc(),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
    )
    return selected, warnings


def _resolve_target_by_tax_id(
    tax_id: int,
    *,
    cache_dir: Path,
    offline: bool,
    strict: bool,
) -> tuple[TaxonomyRecord | None, list[str]]:
    records, warnings = resolve_taxonomy_records(
        {tax_id},
        cache_dir=cache_dir,
        offline=offline,
        strict=strict,
    )
    record = records.get(tax_id)
    if record is None and strict:
        raise MSSPackError(f"Could not resolve configured target TaxID {tax_id}")
    return record, warnings


def _busco_dataset_query(dataset: str) -> str:
    value = Path(dataset).name
    value = re.sub(r"_odb\d+(?:\.\d+)?(?:.*)?$", "", value, flags=re.IGNORECASE)
    return value.replace("_", " ").strip()


def discover_busco_summary_paths(output_root: Path) -> list[Path]:
    busco_root = output_root / "busco"
    if not busco_root.is_dir():
        return []
    return sorted(
        path
        for path in busco_root.glob("*/*.summary.json")
        if path.is_file()
    )


def _busco_evidence(
    configured_lineage: str,
    summary_paths: list[Path],
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    if configured_lineage.strip():
        rows.append(
            {
                "dataset": configured_lineage.strip(),
                "selection_strategy": "configured-lineage",
                "source": "config",
            }
        )
    for path in summary_paths:
        data = _read_json(path)
        dataset = str(data.get("lineage_dataset", "")).strip()
        if not dataset:
            continue
        rows.append(
            {
                "dataset": dataset,
                "selection_strategy": str(data.get("selection_strategy", "unknown")),
                "source": str(path),
            }
        )
    unique: dict[tuple[str, str], dict[str, str]] = {}
    for row in rows:
        unique.setdefault((row["dataset"], row["selection_strategy"]), row)
    return list(unique.values())


def _ranked_lineage(
    record: TaxonomyRecord,
    *,
    cache_dir: Path,
    offline: bool,
    strict: bool,
) -> tuple[list[dict[str, object]], list[str]]:
    ancestor_ids = set(record.lineage_ids) | {record.tax_id}
    ancestors, warnings = resolve_taxonomy_records(
        ancestor_ids,
        cache_dir=cache_dir,
        offline=offline,
        strict=strict,
    )
    lineage: list[dict[str, object]] = []
    for tax_id in (*record.lineage_ids, record.tax_id):
        ancestor = ancestors.get(tax_id)
        if ancestor is None:
            lineage.append({"tax_id": tax_id, "scientific_name": "", "rank": ""})
        else:
            lineage.append(
                {
                    "tax_id": ancestor.tax_id,
                    "scientific_name": ancestor.scientific_name,
                    "rank": ancestor.rank,
                }
            )
    return lineage, warnings


def _lineage_tax_ids(target_payload: dict[str, object]) -> set[int]:
    raw_lineage = target_payload.get("lineage", [])
    if not isinstance(raw_lineage, list):
        return set()
    result: set[int] = set()
    for item in raw_lineage:
        if not isinstance(item, dict):
            continue
        value = item.get("tax_id")
        if isinstance(value, int) and not isinstance(value, bool):
            result.add(value)
    return result


def _target_payload(
    record: TaxonomyRecord,
    *,
    cache_dir: Path,
    offline: bool,
    strict: bool,
) -> tuple[dict[str, object], list[str]]:
    lineage, warnings = _ranked_lineage(
        record,
        cache_dir=cache_dir,
        offline=offline,
        strict=strict,
    )
    return {
        "tax_id": record.tax_id,
        "scientific_name": record.scientific_name,
        "rank": record.rank,
        "blast_name": record.blast_name,
        "lineage": lineage,
    }, warnings


def _crosscheck_busco_lineages(
    target: dict[str, object] | None,
    evidence: list[dict[str, str]],
    *,
    cache_dir: Path,
    config: FunctionalAnnotationTaxonomyConfig,
) -> tuple[list[dict[str, object]], list[str]]:
    rows: list[dict[str, object]] = []
    warnings: list[str] = []
    target_ids = _lineage_tax_ids(target or {})
    for row in evidence:
        dataset_query = _busco_dataset_query(row["dataset"])
        record, record_warnings = resolve_named_taxonomy(
            dataset_query,
            cache_dir=cache_dir,
            offline=config.offline,
            busco_lineage=True,
            strict=False,
        )
        warnings.extend(record_warnings)
        concordance = "unresolved"
        if target is not None and record is not None:
            concordance = "consistent" if record.tax_id in target_ids else "mismatch"
            if concordance == "mismatch":
                warnings.append(
                    f"BUSCO lineage {row['dataset']} is not an ancestor of target "
                    f"TaxID {target.get('tax_id')}"
                )
        rows.append(
            {
                **row,
                "lineage_query": dataset_query,
                "lineage_tax_id": record.tax_id if record else None,
                "lineage_name": record.scientific_name if record else "",
                "concordance": concordance,
            }
        )
    return rows, warnings


def resolve_annotation_taxonomy(
    *,
    scientific_name: str,
    configured_busco_lineage: str,
    busco_summary_paths: list[Path],
    output_path: Path,
    log_path: Path,
    metrics_path: Path,
    cache_dir: Path,
    config: FunctionalAnnotationTaxonomyConfig,
) -> None:
    started_at = datetime.now()
    warnings: list[str] = []
    source = "disabled"
    target_record: TaxonomyRecord | None = None
    busco_rows = _busco_evidence(configured_busco_lineage, busco_summary_paths)

    if config.enabled:
        if config.target_taxon_id > 0:
            source = "configured_taxon_id"
            target_record, current_warnings = _resolve_target_by_tax_id(
                config.target_taxon_id,
                cache_dir=cache_dir,
                offline=config.offline,
                strict=config.strict,
            )
            warnings.extend(current_warnings)
        elif config.resolve_scientific_name:
            source = "scientific_name"
            target_record, current_warnings = resolve_named_taxonomy(
                scientific_name,
                cache_dir=cache_dir,
                offline=config.offline,
                strict=config.strict,
            )
            warnings.extend(current_warnings)
        if target_record is None and busco_rows:
            source = "busco_lineage_fallback"
            dataset_query = _busco_dataset_query(busco_rows[0]["dataset"])
            target_record, current_warnings = resolve_named_taxonomy(
                dataset_query,
                cache_dir=cache_dir,
                offline=config.offline,
                busco_lineage=True,
                strict=config.strict,
            )
            warnings.extend(current_warnings)

    target: dict[str, object] | None = None
    if target_record is not None:
        target, current_warnings = _target_payload(
            target_record,
            cache_dir=cache_dir,
            offline=config.offline,
            strict=config.strict,
        )
        warnings.extend(current_warnings)
    elif config.enabled and config.strict:
        raise MSSPackError(
            "Taxonomy-aware annotation is strict, but no target taxonomy could be resolved"
        )

    crosschecks: list[dict[str, object]] = []
    if config.enabled and config.busco_crosscheck and busco_rows:
        crosschecks, current_warnings = _crosscheck_busco_lineages(
            target,
            busco_rows,
            cache_dir=cache_dir,
            config=config,
        )
        warnings.extend(current_warnings)
        if config.strict and any(row["concordance"] == "mismatch" for row in crosschecks):
            raise MSSPackError("Resolved sample taxonomy conflicts with the BUSCO lineage")

    status = "disabled" if not config.enabled else "resolved" if target else "unresolved"
    payload = {
        "schema_version": TAXONOMY_CONTEXT_SCHEMA_VERSION,
        "status": status,
        "source": source,
        "scientific_name_query": scientific_name,
        "target": target,
        "busco_crosschecks": crosschecks,
        "resolver": {
            "provider": "NCBI Taxonomy",
            "endpoint": NCBI_TAXONOMY_API,
            "offline": config.offline,
        },
        "relation_weights": RELATION_WEIGHTS,
        "warnings": sorted(set(warnings)),
    }
    write_text(output_path, json.dumps(payload, indent=2, sort_keys=True) + "\n")
    target_tax_id = target.get("tax_id") if target else None
    details = [
        f"Taxonomy status: {status}",
        f"Resolution source: {source}",
        f"Target TaxID: {target_tax_id if target_tax_id is not None else 'unresolved'}",
        f"BUSCO lineage checks: {len(crosschecks):,}",
        f"Warnings: {len(set(warnings)):,}",
    ]
    write_step_log(
        log_path=log_path,
        command="msspack functional-annotation resolve-taxonomy",
        step="functional-annotation-resolve-taxonomy",
        started_at=started_at,
        status="completed" if status != "unresolved" else "warning",
        count_unit="taxonomy contexts",
        input_total=1,
        changed_total=int(target is not None),
        output_total=1,
        details=details,
    )
    write_step_metrics(
        metrics_path=metrics_path,
        step="functional-annotation-resolve-taxonomy",
        status="completed" if status != "unresolved" else "warning",
        count_unit="taxonomy contexts",
        input_total=1,
        changed_total=int(target is not None),
        output_total=1,
        details={
            "status": status,
            "source": source,
            "target_tax_id": target_tax_id,
            "busco_crosschecks": len(crosschecks),
            "warnings": len(set(warnings)),
        },
    )


def load_taxonomy_context(path: Path | None) -> dict[str, object]:
    if path is None:
        return {}
    payload = _read_json(path)
    if payload.get("schema_version") != TAXONOMY_CONTEXT_SCHEMA_VERSION:
        return {}
    return payload


def _target_rank_ids(context: dict[str, object]) -> dict[str, int]:
    target = context.get("target")
    if not isinstance(target, dict):
        return {}
    ranks: dict[str, int] = {}
    raw_lineage = target.get("lineage", [])
    if isinstance(raw_lineage, list):
        for item in raw_lineage:
            if not isinstance(item, dict):
                continue
            rank = str(item.get("rank", "")).upper()
            tax_id = item.get("tax_id")
            if rank and isinstance(tax_id, int) and not isinstance(tax_id, bool):
                ranks[rank] = tax_id
    return ranks


def taxonomy_relation(
    *,
    target_context: dict[str, object],
    subject: TaxonomyRecord | None,
) -> tuple[str, float]:
    if subject is None:
        return "unknown", RELATION_WEIGHTS["unknown"]
    target = target_context.get("target")
    if not isinstance(target, dict):
        return "unknown", RELATION_WEIGHTS["unknown"]
    target_tax_id = target.get("tax_id")
    if isinstance(target_tax_id, int) and target_tax_id == subject.tax_id:
        return "same_species", RELATION_WEIGHTS["same_species"]
    subject_ids = set(subject.lineage_ids) | {subject.tax_id}
    target_ranks = _target_rank_ids(target_context)
    for rank, relation in _COMPARISON_RANKS:
        ancestor = target_ranks.get(rank)
        if ancestor is not None and ancestor in subject_ids:
            return relation, RELATION_WEIGHTS[relation]
    target_domain = target_ranks.get("DOMAIN")
    subject_domain = next((value for value in subject_ids if value in NCBI_DOMAIN_IDS), None)
    if target_domain is not None and subject_domain is not None:
        if target_domain != subject_domain:
            return "cross_domain", RELATION_WEIGHTS["cross_domain"]
        if target_ranks.get("KINGDOM") is not None:
            return "cross_kingdom", RELATION_WEIGHTS["cross_kingdom"]
        return "same_domain", RELATION_WEIGHTS["same_domain"]
    return "unknown", RELATION_WEIGHTS["unknown"]


def subject_taxonomy_annotations(
    tax_ids: set[int],
    *,
    target_context_path: Path | None,
    cache_dir: Path,
    offline: bool,
    strict: bool,
) -> tuple[dict[int, tuple[TaxonomyRecord, str, float]], list[str]]:
    context = load_taxonomy_context(target_context_path)
    if not isinstance(context.get("target"), dict):
        return {}, []
    records, warnings = resolve_taxonomy_records(
        tax_ids,
        cache_dir=cache_dir,
        offline=offline,
        strict=strict,
    )
    annotations = {
        tax_id: (record, *taxonomy_relation(target_context=context, subject=record))
        for tax_id, record in records.items()
    }
    return annotations, warnings


__all__ = [
    "RELATION_WEIGHTS",
    "TaxonomyRecord",
    "discover_busco_summary_paths",
    "load_taxonomy_context",
    "resolve_annotation_taxonomy",
    "resolve_named_taxonomy",
    "resolve_taxonomy_records",
    "subject_taxonomy_annotations",
    "taxonomy_relation",
]
