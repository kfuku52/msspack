from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import re
import shutil
import tarfile
import tempfile
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from collections.abc import Sequence
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

from Bio.Seq import Seq

from .annotation_taxonomy import (
    RELATION_WEIGHTS,
    TaxonomyRecord,
    load_taxonomy_context,
    subject_taxonomy_annotations,
)
from .config_models import FunctionalAnnotationConfig
from .database_lock import (
    DatabaseLockSettings,
    acquire_database_lock,
    database_lock_path,
)
from .fasta import iter_fasta, reverse_complement, write_fasta_record
from .gff import GFFRecord, child_ids, read_gff_document
from .product_names import ProductNameContext, standardize_product_name
from .step_logging import write_step_log, write_step_metrics
from .utils import (
    MSSPackError,
    atomic_binary_writer,
    atomic_text_writer,
    copy_or_decompress,
    ensure_dir,
    expand_path,
    link_or_copy,
    run_command,
    write_text,
)

DIAMOND_FIELDS = (
    "qseqid",
    "sseqid",
    "pident",
    "qlen",
    "slen",
    "qcovhsp",
    "scovhsp",
    "evalue",
    "bitscore",
)

DIAMOND_COMMAND_FIELDS = (*DIAMOND_FIELDS, "stitle")

CDD_METADATA_FIELDS = ("pssm_id", "accession", "name", "description", "model_length")

EVIDENCE_FIELDS = (
    "ID",
    "Locus_tag",
    "original_product",
    "proposed_product",
    "assigned_product",
    "candidate_source",
    "source",
    "evidence_id",
    "confidence",
    "name_standardization",
    "name_warnings",
    "quality_code",
    "reason",
    "identity",
    "query_coverage",
    "subject_coverage",
    "evalue",
    "bitscore",
    "consensus_score",
    "subject_taxon_id",
    "subject_organism",
    "taxonomy_relation",
    "taxonomy_weight",
    "taxonomy_adjustment",
)

_UNINFORMATIVE_DESCRIPTION_PATTERNS = tuple(
    re.compile(pattern, re.IGNORECASE)
    for pattern in (
        r"\bhypothetical protein\b",
        r"\buncharacteri[sz]ed protein\b",
        r"\bunknown protein\b",
        r"\bprotein of unknown function\b",
        r"\bunnamed protein product\b",
        r"\bpredicted protein\b",
        r"\bexpressed protein\b",
        r"\bconserved (?:hypothetical )?protein\b",
        r"\bwhole genome shotgun\b",
        r"\bfragment(?:s)?\b",
        r"^(?:probable|putative)\b",
    )
)

_TOKEN_STOPWORDS = {
    "a",
    "and",
    "chain",
    "containing",
    "domain",
    "family",
    "for",
    "fragment",
    "homolog",
    "isoform",
    "like",
    "of",
    "precursor",
    "probable",
    "protein",
    "putative",
    "subunit",
    "the",
}

_UNINFORMATIVE_PFAM_PATTERNS = tuple(
    re.compile(pattern, re.IGNORECASE)
    for pattern in (
        r"\bcoiled[- ]coil\b",
        r"\blow complexity\b",
        r"\bmotif\b",
        r"\bsignal peptide\b",
        r"\btransmembrane\b",
        r"\bdisordered\b",
        r"\bregion\b",
    )
)


@dataclass(frozen=True)
class TranscriptModel:
    transcript_id: str
    seqid: str
    strand: str
    cds_records: tuple[GFFRecord, ...]


@dataclass(frozen=True)
class DiamondDatabase:
    name: str
    fasta_path: Path
    database_path: Path
    weight: float
    origin: str
    sha256: str
    release: str
    verification_algorithm: str
    verification_digest: str


@dataclass(frozen=True)
class DiamondHit:
    query_id: str
    subject_id: str
    database: str
    identity: float
    query_length: int
    subject_length: int
    query_coverage: float
    subject_coverage: float
    evalue: float
    bitscore: float


@dataclass(frozen=True)
class DiamondMetadata:
    description: str
    source_weight: float
    subject_taxon_id: int | None = None
    subject_organism: str = ""
    taxonomy_relation: str = "unknown"
    taxonomy_weight: float = RELATION_WEIGHTS["unknown"]


@dataclass(frozen=True)
class DiamondCandidate:
    hit: DiamondHit
    description: str
    metadata: DiamondMetadata
    tokens: set[str]
    evidence_weight: float


@dataclass(frozen=True)
class PfamHit:
    query_id: str
    name: str
    accession: str
    description: str
    i_evalue: float
    bitscore: float
    hmm_from: int
    hmm_to: int
    hmm_length: int

    @property
    def domain_coverage(self) -> float:
        if self.hmm_length <= 0:
            return 0.0
        return (self.hmm_to - self.hmm_from + 1) / self.hmm_length


@dataclass(frozen=True)
class CddHit:
    query_id: str
    hit_type: str
    pssm_id: str
    accession: str
    name: str
    description: str
    evalue: float
    bitscore: float
    incomplete: str


@dataclass(frozen=True)
class Assignment:
    product: str
    source: str
    evidence_id: str
    confidence: str
    reason: str
    quality_code: str = "---"
    identity: float | None = None
    query_coverage: float | None = None
    subject_coverage: float | None = None
    evalue: float | None = None
    bitscore: float | None = None
    consensus_score: float | None = None
    subject_taxon_id: int | None = None
    subject_organism: str = ""
    taxonomy_relation: str = ""
    taxonomy_weight: float | None = None
    taxonomy_adjustment: str = ""


def _sha256(path: Path) -> str:
    return _file_digest(path, "sha256")


def _file_digest(path: Path, algorithm: str) -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _fetch_text(url: str) -> str:
    request = urllib.request.Request(url, headers={"User-Agent": "msspack-functional-annotation"})
    try:
        with urllib.request.urlopen(request, timeout=60) as response:
            payload: bytes = response.read()
            return payload.decode("utf-8")
    except (OSError, UnicodeError, ValueError) as exc:
        raise MSSPackError(
            f"Could not retrieve database checksum metadata from {url}: {exc}"
        ) from exc


def _uniprot_verification(url: str) -> tuple[str, str, int | None, str]:
    metadata_url = url.rsplit("/", 1)[0] + "/RELEASE.metalink"
    try:
        root = ET.fromstring(_fetch_text(metadata_url))
    except ET.ParseError as exc:
        raise MSSPackError(f"Could not parse UniProt release metadata: {metadata_url}") from exc
    filename = Path(urllib.parse.urlparse(url).path).name
    release = ""
    for element in root.iter():
        if element.tag.rsplit("}", 1)[-1] == "version" and element.text:
            release = element.text.strip()
            break
    for file_element in root.iter():
        if file_element.tag.rsplit("}", 1)[-1] != "file":
            continue
        if file_element.attrib.get("name") != filename:
            continue
        size: int | None = None
        algorithm = ""
        checksum = ""
        for child in file_element.iter():
            tag = child.tag.rsplit("}", 1)[-1]
            if tag == "size" and child.text:
                size = int(child.text.strip())
            elif tag == "hash" and child.text:
                algorithm = child.attrib.get("type", "").lower().replace("-", "")
                checksum = child.text.strip().lower()
        if algorithm and checksum:
            return algorithm, checksum, size, release
    raise MSSPackError(f"UniProt release metadata has no checksum for {filename}")


def _pfam_verification(url: str) -> tuple[str, str, int | None, str]:
    metadata_url = url.rsplit("/", 1)[0] + "/md5_checksums"
    filename = Path(urllib.parse.urlparse(url).path).name
    for line in _fetch_text(metadata_url).splitlines():
        fields = line.split()
        if len(fields) >= 2 and fields[-1].lstrip("*") == filename:
            return "md5", fields[0].lower(), None, "current_release"
    raise MSSPackError(f"Pfam release metadata has no checksum for {filename}")


def _remote_verification(url: str) -> tuple[str, str, int | None, str] | None:
    parsed = urllib.parse.urlparse(url)
    host = (parsed.hostname or "").lower()
    if host == "ftp.uniprot.org" and (
        "/uniprot/current_release/" in parsed.path or "/uniprot/uniref/" in parsed.path
    ):
        return _uniprot_verification(url)
    if host.endswith("ebi.ac.uk") and "/Pfam/current_release/" in parsed.path:
        return _pfam_verification(url)
    return None


def _read_json(path: Path) -> dict[str, object]:
    if not path.is_file():
        return {}
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}
    return payload if isinstance(payload, dict) else {}


def _download(url: str, destination: Path) -> None:
    request = urllib.request.Request(url, headers={"User-Agent": "msspack-functional-annotation"})
    try:
        with (
            urllib.request.urlopen(request, timeout=300) as response,
            atomic_binary_writer(destination) as output,
        ):
            while chunk := response.read(1024 * 1024):
                output.write(chunk)
    except (OSError, ValueError) as exc:
        raise MSSPackError(f"Could not download annotation database from {url}: {exc}") from exc
    if destination.stat().st_size == 0:
        raise MSSPackError(f"Downloaded annotation database is empty: {url}")


def _copy_binary(source: Path, destination: Path) -> None:
    with source.open("rb") as input_handle, atomic_binary_writer(destination) as output_handle:
        shutil.copyfileobj(input_handle, output_handle, length=1024 * 1024)


def _materialize_database_file_unlocked(
    *,
    local_value: str,
    url: str,
    base_dir: Path,
    cache_dir: Path,
    stem: str,
    expected_suffix: str,
    preserve_compression: bool = False,
) -> tuple[Path, str, str, dict[str, object]]:
    ensure_dir(cache_dir)
    output_path = cache_dir / f"{stem}{expected_suffix}"
    provenance_path = cache_dir / f"{stem}.source.json"
    previous = _read_json(provenance_path)
    verification_algorithm = ""
    verification_digest = ""
    release = ""

    if local_value.strip():
        source_path = expand_path(local_value, base_dir)
        if not source_path.is_file():
            raise MSSPackError(f"Annotation database file not found: {source_path}")
        source_sha256 = _sha256(source_path)
        origin = str(source_path)
        if previous.get("origin_sha256") != source_sha256 or not output_path.is_file():
            if preserve_compression:
                if source_path.suffix.casefold() != ".gz":
                    raise MSSPackError(
                        "A compression-preserving annotation database must be gzip-compressed: "
                        f"{source_path}"
                    )
                _copy_binary(source_path, output_path)
            else:
                copy_or_decompress(source_path, output_path)
    else:
        if not url.strip():
            raise MSSPackError(f"No URL configured for annotation database {stem}")
        origin = url
        if previous.get("origin") != origin or not output_path.is_file():
            verification = _remote_verification(url)
            expected_size: int | None = None
            if verification is not None:
                (
                    verification_algorithm,
                    verification_digest,
                    expected_size,
                    release,
                ) = verification
            archive_name = f"{stem}.download" + (
                ".gz" if preserve_compression or url.lower().endswith(".gz") else ""
            )
            archive_path = cache_dir / archive_name
            _download(url, archive_path)
            observed_size = archive_path.stat().st_size
            if expected_size is not None and observed_size != expected_size:
                archive_path.unlink(missing_ok=True)
                raise MSSPackError(
                    f"Downloaded database size mismatch for {url}: "
                    f"expected {expected_size}, got {observed_size}"
                )
            if verification_algorithm:
                observed_digest = _file_digest(archive_path, verification_algorithm)
                if observed_digest.lower() != verification_digest:
                    archive_path.unlink(missing_ok=True)
                    raise MSSPackError(
                        f"Downloaded database checksum mismatch for {url}: "
                        f"expected {verification_digest}, got {observed_digest}"
                    )
            if preserve_compression:
                os.replace(archive_path, output_path)
            else:
                copy_or_decompress(archive_path, output_path)
                archive_path.unlink(missing_ok=True)
        source_sha256 = str(previous.get("origin_sha256", ""))
        if not verification_algorithm:
            verification_algorithm = str(previous.get("verification_algorithm", ""))
            verification_digest = str(previous.get("verification_digest", ""))
            release = str(previous.get("release", ""))

    if not output_path.is_file() or output_path.stat().st_size == 0:
        raise MSSPackError(f"Prepared annotation database is empty: {output_path}")
    materialized_sha256 = _sha256(output_path)
    if not source_sha256:
        source_sha256 = materialized_sha256
    write_text(
        provenance_path,
        json.dumps(
            {
                "origin": origin,
                "origin_sha256": source_sha256,
                "release": release,
                "verification_algorithm": verification_algorithm,
                "verification_digest": verification_digest,
                "materialized_path": str(output_path),
                "materialized_sha256": materialized_sha256,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
    )
    return output_path, origin, materialized_sha256, _read_json(provenance_path)


def _materialize_database_file(
    *,
    local_value: str,
    url: str,
    base_dir: Path,
    cache_dir: Path,
    database_root: Path | None = None,
    stem: str,
    expected_suffix: str,
    preserve_compression: bool = False,
    content_addressed: bool = False,
    lock_settings: DatabaseLockSettings | None = None,
) -> tuple[Path, str, str, dict[str, object]]:
    resolved_database_root = database_root or cache_dir
    lock_path = database_lock_path(resolved_database_root, f"download-{stem}")
    with acquire_database_lock(
        lock_path,
        label=f"{stem} database preparation",
        settings=lock_settings,
    ):
        path, origin, digest, provenance = _materialize_database_file_unlocked(
            local_value=local_value,
            url=url,
            base_dir=base_dir,
            cache_dir=cache_dir,
            stem=stem,
            expected_suffix=expected_suffix,
            preserve_compression=preserve_compression,
        )
        if not content_addressed:
            return path, origin, digest, provenance
        suffix = "".join(path.suffixes) or expected_suffix
        object_path = cache_dir / "objects" / f"{digest}{suffix}"
        if not object_path.is_file():
            _copy_binary(path, object_path)
        if object_path.stat().st_size == 0:
            raise MSSPackError(f"Content-addressed database file is empty: {object_path}")
        return (
            object_path,
            origin,
            digest,
            {
                **provenance,
                "content_object_path": str(object_path),
            },
        )


def _build_transcript_models(gff_path: Path) -> list[TranscriptModel]:
    transcript_records: dict[str, GFFRecord] = {}
    transcript_order: list[str] = []
    cds_by_parent: dict[str, list[GFFRecord]] = defaultdict(list)
    for record in read_gff_document(gff_path).records:
        record_id = record.attributes.get("ID", "")
        if record.type in ("mRNA", "transcript") and record_id:
            transcript_records[record_id] = record
            transcript_order.append(record_id)
        elif record.type == "CDS":
            for parent_id in child_ids(record.attributes.get("Parent")):
                cds_by_parent[parent_id].append(record)
                if parent_id not in transcript_records and parent_id not in transcript_order:
                    transcript_order.append(parent_id)

    models: list[TranscriptModel] = []
    for transcript_id in transcript_order:
        cds_records = cds_by_parent.get(transcript_id, [])
        if not cds_records:
            continue
        transcript = transcript_records.get(transcript_id)
        models.append(
            TranscriptModel(
                transcript_id=transcript_id,
                seqid=transcript.seqid if transcript else cds_records[0].seqid,
                strand=transcript.strand if transcript else cds_records[0].strand,
                cds_records=tuple(sorted(cds_records, key=lambda item: (item.start, item.end))),
            )
        )
    return models


def write_translated_protein_fasta(
    *,
    fasta_path: Path,
    gff_path: Path,
    output_path: Path,
    genetic_code: str,
    log_path: Path,
    metrics_path: Path,
) -> None:
    started_at = datetime.now()
    models_by_seqid: dict[str, list[TranscriptModel]] = defaultdict(list)
    for model in _build_transcript_models(gff_path):
        models_by_seqid[model.seqid].append(model)

    written = 0
    skipped = 0
    internal_stops = 0
    with atomic_text_writer(output_path) as output:
        for fasta_record in iter_fasta(fasta_path):
            for model in models_by_seqid.get(fasta_record.id, []):
                sequence = "".join(
                    fasta_record.sequence[cds.start - 1 : cds.end] for cds in model.cds_records
                )
                if model.strand == "-":
                    sequence = reverse_complement(sequence)
                usable_length = len(sequence) - (len(sequence) % 3)
                if usable_length < 3:
                    skipped += 1
                    continue
                protein = str(
                    Seq(sequence[:usable_length]).translate(
                        table=int(genetic_code)
                    )
                )
                if protein.endswith("*"):
                    protein = protein[:-1]
                stop_count = protein.count("*")
                internal_stops += stop_count
                protein = protein.replace("*", "X")
                if not protein:
                    skipped += 1
                    continue
                write_fasta_record(output, identifier=model.transcript_id, sequence=protein)
                written += 1

    command = [
        "msspack",
        "functional-annotation",
        "extract-proteins",
        "--fasta",
        str(fasta_path),
        "--gff",
        str(gff_path),
    ]
    write_step_log(
        log_path=log_path,
        command=command,
        step="functional-annotation-extract-proteins",
        started_at=started_at,
        count_unit="proteins",
        input_total=written + skipped,
        changed_total=written,
        output_total=written,
        details=[
            f"Skipped CDS models shorter than one codon: {skipped:,}",
            f"Internal stop codons replaced with X: {internal_stops:,}",
            f"NCBI genetic code: {genetic_code}",
        ],
    )
    write_step_metrics(
        metrics_path=metrics_path,
        step="functional-annotation-extract-proteins",
        count_unit="proteins",
        input_total=written + skipped,
        changed_total=written,
        output_total=written,
        details={"skipped_short_cds": skipped, "internal_stops_replaced": internal_stops},
    )


def _fasta_description(record_description: str) -> str:
    parts = record_description.split(None, 1)
    description = parts[1] if len(parts) == 2 else ""
    taxonomy_marker = re.search(
        r"\s(?:OS|OX|GN|PE|SV|n|Tax|TaxID|RepID)=",
        description,
    )
    if taxonomy_marker:
        description = description[: taxonomy_marker.start()]
    description = re.sub(r"\s+", " ", description).strip().strip(";")
    return description


def _fasta_taxonomy(record_description: str) -> tuple[int | None, str]:
    tax_id_match = re.search(r"\b(?:OX|TaxID)=(\d+)\b", record_description)
    tax_id = int(tax_id_match.group(1)) if tax_id_match else None
    organism = ""
    organism_match = re.search(
        r"\bOS=(.+?)(?=\s(?:OX|GN|PE|SV)=|$)",
        record_description,
    )
    if organism_match is None:
        organism_match = re.search(
            r"\bTax=(.+?)(?=\s(?:TaxID|RepID)=|$)",
            record_description,
        )
    if organism_match is not None:
        organism = re.sub(r"\s+", " ", organism_match.group(1)).strip()
    return tax_id, organism


def _uniref90_download_url(config: FunctionalAnnotationConfig) -> str:
    if config.uniref90_taxon_id <= 0:
        return config.uniref90_url
    query = urllib.parse.quote(
        f"(identity:0.9) AND (taxonomy_id:{config.uniref90_taxon_id})"
    )
    return (
        "https://rest.uniprot.org/uniref/stream?compressed=true&format=fasta&query="
        f"{query}"
    )


def _prepare_diamond_database(
    *,
    name: str,
    fasta_path: Path,
    weight: float,
    origin: str,
    sha256: str,
    source_provenance: dict[str, object],
    cache_dir: Path,
    database_root: Path,
    command: str,
    lock_settings: DatabaseLockSettings | None = None,
) -> DiamondDatabase:
    safe_name = re.sub(r"[^A-Za-z0-9_.-]+", "_", name).strip("._") or "database"
    legacy_dir = cache_dir / safe_name
    legacy_database_path = legacy_dir / f"{safe_name}.dmnd"
    legacy_provenance = _read_json(legacy_dir / "database.json")
    if (
        legacy_database_path.is_file()
        and legacy_provenance.get("fasta_sha256") == sha256
    ):
        return DiamondDatabase(
            name,
            fasta_path,
            legacy_database_path,
            weight,
            origin,
            sha256,
            str(source_provenance.get("release", "")),
            str(source_provenance.get("verification_algorithm", "")),
            str(source_provenance.get("verification_digest", "")),
        )
    version = sha256[:16]
    database_dir = ensure_dir(cache_dir / safe_name / version)
    database_path = database_dir / f"{safe_name}.dmnd"
    provenance_path = database_dir / "database.json"
    with acquire_database_lock(
        database_lock_path(database_root, f"diamond-{safe_name}-{version}"),
        label=f"{name} DIAMOND database build",
        settings=lock_settings,
    ):
        previous = _read_json(provenance_path)
        if previous.get("fasta_sha256") != sha256 or not database_path.is_file():
            temporary_base = database_dir / f".{safe_name}.{os.getpid()}"
            run_command(
                [command, "makedb", "--in", str(fasta_path), "--db", str(temporary_base)],
                log_path=database_dir / "makedb.log",
            )
            built_path = Path(str(temporary_base) + ".dmnd")
            if not built_path.is_file():
                raise MSSPackError(f"DIAMOND did not create its database: {built_path}")
            os.replace(built_path, database_path)
            write_text(
                provenance_path,
                json.dumps(
                    {
                        "database": name,
                        "origin": origin,
                        "fasta": str(fasta_path),
                        "fasta_sha256": sha256,
                        "release": source_provenance.get("release", ""),
                        "verification_algorithm": source_provenance.get(
                            "verification_algorithm", ""
                        ),
                        "verification_digest": source_provenance.get(
                            "verification_digest", ""
                        ),
                        "diamond_command": command,
                    },
                    indent=2,
                    sort_keys=True,
                )
                + "\n",
            )
    return DiamondDatabase(
        name,
        fasta_path,
        database_path,
        weight,
        origin,
        sha256,
        str(source_provenance.get("release", "")),
        str(source_provenance.get("verification_algorithm", "")),
        str(source_provenance.get("verification_digest", "")),
    )


def _diamond_sources(
    *,
    config: FunctionalAnnotationConfig,
    base_dir: Path,
    cache_dir: Path,
    source_group: str,
    lock_settings: DatabaseLockSettings | None = None,
) -> list[DiamondDatabase]:
    downloads_dir = ensure_dir(cache_dir / "downloads")
    database_cache = ensure_dir(cache_dir / "diamond")
    sources: list[DiamondDatabase] = []
    if source_group == "primary" and config.swissprot_enabled:
        fasta, origin, digest, source_provenance = _materialize_database_file(
            local_value=config.swissprot_fasta,
            url=config.swissprot_url,
            base_dir=base_dir,
            cache_dir=downloads_dir / "swissprot",
            database_root=cache_dir,
            stem="uniprot_sprot",
            expected_suffix=".fasta",
            content_addressed=True,
            lock_settings=lock_settings,
        )
        sources.append(
            _prepare_diamond_database(
                name="swissprot",
                fasta_path=fasta,
                weight=config.swissprot_weight,
                origin=origin,
                sha256=digest,
                source_provenance=source_provenance,
                cache_dir=database_cache,
                database_root=cache_dir,
                command=config.diamond_command,
                lock_settings=lock_settings,
            )
        )
    if source_group == "primary" and config.reference_proteins.strip():
        fasta, origin, digest, source_provenance = _materialize_database_file(
            local_value=config.reference_proteins,
            url="",
            base_dir=base_dir,
            cache_dir=downloads_dir / "reference",
            database_root=cache_dir,
            stem="reference_proteins",
            expected_suffix=".fasta",
            content_addressed=True,
            lock_settings=lock_settings,
        )
        sources.append(
            _prepare_diamond_database(
                name=config.reference_name,
                fasta_path=fasta,
                weight=config.reference_weight,
                origin=origin,
                sha256=digest,
                source_provenance=source_provenance,
                cache_dir=database_cache,
                database_root=cache_dir,
                command=config.diamond_command,
                lock_settings=lock_settings,
            )
        )
    if source_group == "uniref90" and config.uniref90_enabled:
        preserve_compression = not config.uniref90_fasta.strip() or Path(
            config.uniref90_fasta
        ).suffix.casefold() == ".gz"
        fasta, origin, digest, source_provenance = _materialize_database_file(
            local_value=config.uniref90_fasta,
            url=_uniref90_download_url(config),
            base_dir=base_dir,
            cache_dir=downloads_dir / "uniref90",
            database_root=cache_dir,
            stem=(
                f"uniref90-taxon-{config.uniref90_taxon_id}"
                if config.uniref90_taxon_id > 0
                else "uniref90"
            ),
            expected_suffix=".fasta.gz" if preserve_compression else ".fasta",
            preserve_compression=preserve_compression,
            content_addressed=True,
            lock_settings=lock_settings,
        )
        sources.append(
            _prepare_diamond_database(
                name="uniref90",
                fasta_path=fasta,
                weight=config.uniref90_weight,
                origin=origin,
                sha256=digest,
                source_provenance=source_provenance,
                cache_dir=database_cache,
                database_root=cache_dir,
                command=config.diamond_command,
                lock_settings=lock_settings,
            )
        )
    if source_group not in {"primary", "uniref90"}:
        raise MSSPackError(f"Unknown DIAMOND annotation source group: {source_group}")
    return sources


def write_empty_diamond_results(
    *,
    output_path: Path,
    metadata_path: Path,
    provenance_path: Path,
    log_path: Path,
    metrics_path: Path,
    reason: str = "No DIAMOND protein database is enabled.",
    step_name: str = "functional-annotation-diamond-search",
) -> None:
    started_at = datetime.now()
    write_text(output_path, "database\t" + "\t".join(DIAMOND_FIELDS) + "\n")
    write_text(
        metadata_path,
        "database\tsubject_id\tdescription\tweight\tsubject_taxon_id\t"
        "subject_organism\ttaxonomy_relation\ttaxonomy_weight\n",
    )
    write_text(provenance_path, '{"sources": [], "status": "skipped"}\n')
    write_step_log(
        log_path=log_path,
        command="msspack functional-annotation diamond-search --skipped",
        step=step_name,
        started_at=started_at,
        status="skipped",
        count_unit="alignments",
        input_total=0,
        changed_total=0,
        output_total=0,
        details=[reason],
    )
    write_step_metrics(
        metrics_path=metrics_path,
        step=step_name,
        status="skipped",
        count_unit="alignments",
        input_total=0,
        changed_total=0,
        output_total=0,
    )


def run_diamond_annotation_search(
    *,
    protein_fasta_path: Path,
    output_path: Path,
    metadata_path: Path,
    provenance_path: Path,
    log_path: Path,
    metrics_path: Path,
    config: FunctionalAnnotationConfig,
    base_dir: Path,
    cache_dir: Path,
    lock_settings: DatabaseLockSettings | None = None,
    taxonomy_context_path: Path | None = None,
    taxonomy_cache_dir: Path | None = None,
    source_group: str = "primary",
    prior_similarity_inputs: Sequence[tuple[Path, Path]] = (),
    step_name: str = "functional-annotation-diamond-search",
) -> None:
    started_at = datetime.now()
    protein_count = sum(1 for _ in iter_fasta(protein_fasta_path))
    if protein_count == 0:
        write_empty_diamond_results(
            output_path=output_path,
            metadata_path=metadata_path,
            provenance_path=provenance_path,
            log_path=log_path,
            metrics_path=metrics_path,
            reason="No translated protein sequences were available for DIAMOND.",
            step_name=step_name,
        )
        return
    sources = _diamond_sources(
        config=config,
        base_dir=base_dir,
        cache_dir=cache_dir,
        source_group=source_group,
        lock_settings=lock_settings,
    )
    prior_assignments = _combined_diamond_assignments(
        inputs=prior_similarity_inputs,
        config=config,
    )
    skipped_ids = {
        query_id
        for query_id, assignment in prior_assignments.items()
        if assignment.confidence == "high"
    }
    alignment_count = 0
    query_ids: set[str] = set()
    command_example: list[str] | None = None
    matched_metadata: dict[
        tuple[str, str],
        tuple[str, float, int | None, str],
    ] = {}

    with tempfile.TemporaryDirectory(prefix="msspack-diamond-") as temporary_dir:
        temporary_root = Path(temporary_dir)
        query_path = protein_fasta_path
        query_count = protein_count
        if skipped_ids:
            query_path = temporary_root / "unassigned.fasta"
            query_count = 0
            with atomic_text_writer(query_path) as query_handle:
                for record in iter_fasta(protein_fasta_path):
                    if record.id in skipped_ids:
                        continue
                    write_fasta_record(
                        query_handle,
                        identifier=record.id,
                        description=record.description,
                        sequence=record.sequence,
                    )
                    query_count += 1
        if query_count == 0:
            write_empty_diamond_results(
                output_path=output_path,
                metadata_path=metadata_path,
                provenance_path=provenance_path,
                log_path=log_path,
                metrics_path=metrics_path,
                reason="All translated proteins already have accepted prior assignments.",
                step_name=step_name,
            )
            return
        with atomic_text_writer(output_path) as combined:
            combined.write("database\t" + "\t".join(DIAMOND_FIELDS) + "\n")
            for index, source in enumerate(sources):
                raw_output = temporary_root / f"{index}.tsv"
                command = [
                    config.diamond_command,
                    "blastp",
                    "--query",
                    str(query_path),
                    "--db",
                    str(source.database_path),
                    "--out",
                    str(raw_output),
                    "--outfmt",
                    "6",
                    *DIAMOND_COMMAND_FIELDS,
                    "--threads",
                    str(config.threads),
                    "--evalue",
                    f"{config.evalue:g}",
                    "--max-target-seqs",
                    str(config.max_target_seqs),
                ]
                if config.sensitivity != "default":
                    command.append(f"--{config.sensitivity}")
                command_example = command
                run_command(command, log_path=temporary_root / f"{index}.log")
                if not raw_output.is_file():
                    raise MSSPackError(f"DIAMOND did not create its result file: {raw_output}")
                with raw_output.open("r", encoding="utf-8") as raw_handle:
                    for raw_line in raw_handle:
                        line = raw_line.rstrip("\n")
                        if not line:
                            continue
                        fields = line.split("\t", maxsplit=len(DIAMOND_FIELDS))
                        if len(fields) != len(DIAMOND_COMMAND_FIELDS):
                            raise MSSPackError(
                                f"Unexpected DIAMOND output with {len(fields)} columns"
                            )
                        core_fields = fields[: len(DIAMOND_FIELDS)]
                        description = _fasta_description(fields[-1])
                        subject_taxon_id, subject_organism = _fasta_taxonomy(fields[-1])
                        combined.write(source.name + "\t" + "\t".join(core_fields) + "\n")
                        matched_metadata[(source.name, fields[1])] = (
                            description,
                            source.weight,
                            subject_taxon_id,
                            subject_organism,
                        )
                        query_ids.add(fields[0])
                        alignment_count += 1

    subject_tax_ids = {
        tax_id
        for _description, _weight, tax_id, _organism in matched_metadata.values()
        if tax_id is not None
    }
    taxonomy_annotations: dict[int, tuple[TaxonomyRecord, str, float]] = {}
    taxonomy_warnings: list[str] = []
    if taxonomy_context_path is not None and taxonomy_cache_dir is not None:
        taxonomy_annotations, taxonomy_warnings = subject_taxonomy_annotations(
            subject_tax_ids,
            target_context_path=taxonomy_context_path,
            cache_dir=taxonomy_cache_dir,
            offline=config.taxonomy.offline,
            strict=config.taxonomy.strict,
        )
    with atomic_text_writer(metadata_path) as metadata:
        metadata.write(
            "database\tsubject_id\tdescription\tweight\tsubject_taxon_id\t"
            "subject_organism\ttaxonomy_relation\ttaxonomy_weight\n"
        )
        for (database, subject_id), (
            description,
            weight,
            subject_taxon_id,
            subject_organism,
        ) in sorted(
            matched_metadata.items()
        ):
            relation = "reference" if database not in {"swissprot", "uniref90"} else "unknown"
            taxonomy_weight = RELATION_WEIGHTS[relation]
            if subject_taxon_id is not None and subject_taxon_id in taxonomy_annotations:
                _record, relation, taxonomy_weight = taxonomy_annotations[subject_taxon_id]
            metadata.write(
                "\t".join(
                    (
                        database,
                        subject_id,
                        description.replace("\t", " "),
                        f"{weight:g}",
                        str(subject_taxon_id or ""),
                        subject_organism.replace("\t", " "),
                        relation,
                        f"{taxonomy_weight:g}",
                    )
                )
                + "\n"
            )
    metadata_count = len(matched_metadata)

    write_text(
        provenance_path,
        json.dumps(
            {
                "diamond_command": config.diamond_command,
                "fields": list(DIAMOND_FIELDS),
                "source_group": source_group,
                "queries_scanned": query_count,
                "prior_assignments_skipped": len(skipped_ids),
                "taxonomy_context": (
                    str(taxonomy_context_path) if taxonomy_context_path is not None else ""
                ),
                "subject_taxon_ids": len(subject_tax_ids),
                "resolved_subject_taxonomies": len(taxonomy_annotations),
                "taxonomy_warnings": taxonomy_warnings,
                "sources": [
                    {
                        "name": source.name,
                        "origin": source.origin,
                        "fasta_sha256": source.sha256,
                        "release": source.release,
                        "verification_algorithm": source.verification_algorithm,
                        "verification_digest": source.verification_digest,
                        "database_path": str(source.database_path),
                        "weight": source.weight,
                    }
                    for source in sources
                ],
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
    )
    write_step_log(
        log_path=log_path,
        command=command_example or "msspack functional-annotation diamond-search --skipped",
        step=step_name,
        started_at=started_at,
        count_unit="alignments",
        input_total=protein_count,
        changed_total=len(query_ids),
        output_total=alignment_count,
        details=[
            "Databases: " + ", ".join(source.name for source in sources),
            f"Queries scanned: {query_count:,}",
            f"Accepted prior assignments skipped: {len(skipped_ids):,}",
            f"Queries with hits: {len(query_ids):,}",
            f"Matched database metadata records: {metadata_count:,}",
            f"Resolved subject taxonomies: {len(taxonomy_annotations):,}/"
            f"{len(subject_tax_ids):,}",
            f"Taxonomy warnings: {len(taxonomy_warnings):,}",
        ],
    )
    write_step_metrics(
        metrics_path=metrics_path,
        step=step_name,
        count_unit="alignments",
        input_total=protein_count,
        changed_total=len(query_ids),
        output_total=alignment_count,
        details={
            "databases": [source.name for source in sources],
            "queries_scanned": query_count,
            "prior_assignments_skipped": len(skipped_ids),
            "queries_with_hits": len(query_ids),
            "metadata_records": metadata_count,
            "subject_taxon_ids": len(subject_tax_ids),
            "resolved_subject_taxonomies": len(taxonomy_annotations),
            "taxonomy_warnings": len(taxonomy_warnings),
        },
    )


def _parse_pfam_hmm_metadata(hmm_path: Path, output_path: Path) -> int:
    count = 0
    current: dict[str, str] = {}
    with (
        hmm_path.open("r", encoding="utf-8", errors="replace") as handle,
        atomic_text_writer(output_path) as output,
    ):
        output.write("name\taccession\tdescription\n")
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line == "//":
                name = current.get("NAME", "")
                accession = current.get("ACC", "").split(".", 1)[0]
                description = current.get("DESC", "")
                if name:
                    output.write(f"{name}\t{accession}\t{description.replace(chr(9), ' ')}\n")
                    count += 1
                current = {}
                continue
            for field in ("NAME", "ACC", "DESC"):
                prefix = field + "  "
                if line.startswith(prefix):
                    current[field] = line[len(prefix) :].strip()
                    break
    return count


def _prepare_pfam_database(
    *,
    config: FunctionalAnnotationConfig,
    base_dir: Path,
    cache_dir: Path,
    lock_settings: DatabaseLockSettings | None = None,
) -> tuple[Path, Path, dict[str, object]]:
    pfam_dir = ensure_dir(cache_dir / "pfam")
    hmm_path, origin, digest, source_provenance = _materialize_database_file(
        local_value=config.pfam_hmm,
        url=config.pfam_url,
        base_dir=base_dir,
        cache_dir=pfam_dir,
        database_root=cache_dir,
        stem="Pfam-A",
        expected_suffix=".hmm",
        content_addressed=True,
        lock_settings=lock_settings,
    )
    metadata_path = hmm_path.with_name(f"{hmm_path.stem}.metadata.tsv")
    prepared_path = hmm_path.with_name(f"{hmm_path.stem}.prepared.json")
    pressed_paths = [Path(str(hmm_path) + suffix) for suffix in (".h3f", ".h3i", ".h3m", ".h3p")]
    with acquire_database_lock(
        database_lock_path(cache_dir, f"pfam-index-{digest[:16]}"),
        label="Pfam hmmpress build",
        settings=lock_settings,
    ):
        previous = _read_json(prepared_path)
        if previous.get("hmm_sha256") != digest or not all(
            path.is_file() for path in pressed_paths
        ):
            for path in pressed_paths:
                path.unlink(missing_ok=True)
            run_command(
                [config.hmmpress_command, "-f", str(hmm_path)],
                log_path=hmm_path.with_name(f"{hmm_path.stem}.hmmpress.log"),
            )
            if not all(path.is_file() for path in pressed_paths):
                raise MSSPackError("hmmpress did not create all four Pfam index files")
            metadata_count = _parse_pfam_hmm_metadata(hmm_path, metadata_path)
            write_text(
                prepared_path,
                json.dumps(
                    {
                        "origin": origin,
                        "hmm_sha256": digest,
                        "metadata_records": metadata_count,
                        "hmmpress_command": config.hmmpress_command,
                    },
                    indent=2,
                    sort_keys=True,
                )
                + "\n",
            )
        elif not metadata_path.is_file():
            metadata_count = _parse_pfam_hmm_metadata(hmm_path, metadata_path)
        else:
            previous_count = previous.get("metadata_records", 0)
            metadata_count = (
                previous_count
                if isinstance(previous_count, int) and not isinstance(previous_count, bool)
                else 0
            )
    return (
        hmm_path,
        metadata_path,
        {
            "origin": origin,
            "hmm_sha256": digest,
            "release": source_provenance.get("release", ""),
            "verification_algorithm": source_provenance.get("verification_algorithm", ""),
            "verification_digest": source_provenance.get("verification_digest", ""),
            "metadata_records": metadata_count,
        },
    )


def write_empty_pfam_results(
    *,
    output_path: Path,
    metadata_path: Path,
    provenance_path: Path,
    log_path: Path,
    metrics_path: Path,
    reason: str = "Pfam fallback is disabled.",
) -> None:
    started_at = datetime.now()
    write_text(output_path, "")
    write_text(metadata_path, "name\taccession\tdescription\n")
    write_text(provenance_path, '{"status": "skipped"}\n')
    write_step_log(
        log_path=log_path,
        command="msspack functional-annotation pfam-search --skipped",
        step="functional-annotation-pfam-search",
        started_at=started_at,
        status="skipped",
        count_unit="domain hits",
        input_total=0,
        changed_total=0,
        output_total=0,
        details=[reason],
    )
    write_step_metrics(
        metrics_path=metrics_path,
        step="functional-annotation-pfam-search",
        status="skipped",
        count_unit="domain hits",
        input_total=0,
        changed_total=0,
        output_total=0,
    )


def run_pfam_domain_search(
    *,
    protein_fasta_path: Path,
    output_path: Path,
    metadata_path: Path,
    provenance_path: Path,
    log_path: Path,
    metrics_path: Path,
    config: FunctionalAnnotationConfig,
    base_dir: Path,
    cache_dir: Path,
    lock_settings: DatabaseLockSettings | None = None,
    diamond_hits_path: Path | None = None,
    diamond_metadata_path: Path | None = None,
    additional_similarity_inputs: Sequence[tuple[Path, Path]] = (),
) -> None:
    started_at = datetime.now()
    protein_count = sum(1 for _ in iter_fasta(protein_fasta_path))
    if protein_count == 0:
        write_empty_pfam_results(
            output_path=output_path,
            metadata_path=metadata_path,
            provenance_path=provenance_path,
            log_path=log_path,
            metrics_path=metrics_path,
            reason="No translated protein sequences were available for Pfam.",
        )
        return
    similarity_assignments: set[str] = set()
    similarity_inputs = list(additional_similarity_inputs)
    if diamond_hits_path is not None and diamond_metadata_path is not None:
        similarity_inputs.insert(0, (diamond_hits_path, diamond_metadata_path))
    if similarity_inputs:
        similarity_assignments = set(
            _combined_diamond_assignments(inputs=similarity_inputs, config=config)
        )
    fallback_query_count = protein_count - len(similarity_assignments)
    if fallback_query_count <= 0:
        write_empty_pfam_results(
            output_path=output_path,
            metadata_path=metadata_path,
            provenance_path=provenance_path,
            log_path=log_path,
            metrics_path=metrics_path,
            reason="All translated proteins already have accepted similarity assignments.",
        )
        return
    hmm_path, cached_metadata_path, provenance = _prepare_pfam_database(
        config=config,
        base_dir=base_dir,
        cache_dir=cache_dir,
        lock_settings=lock_settings,
    )
    copy_or_decompress(cached_metadata_path, metadata_path)
    with tempfile.TemporaryDirectory(prefix="msspack-hmmscan-") as temporary_dir:
        temporary_root = Path(temporary_dir)
        fallback_records = [
            record
            for record in iter_fasta(protein_fasta_path)
            if record.id not in similarity_assignments
        ]
        if len(fallback_records) != fallback_query_count:
            raise MSSPackError(
                "Pfam fallback query count changed while preparing hmmscan shards"
            )
        shard_count = min(config.threads, fallback_query_count)
        base_shard_size, larger_shards = divmod(fallback_query_count, shard_count)
        commands: list[list[str]] = []
        shard_domtbls: list[Path] = []
        record_offset = 0
        for shard_index in range(shard_count):
            shard_size = base_shard_size + int(shard_index < larger_shards)
            records = fallback_records[record_offset : record_offset + shard_size]
            record_offset += shard_size
            shard_fasta = temporary_root / f"pfam-fallback.{shard_index:03d}.fasta"
            with atomic_text_writer(shard_fasta) as query_handle:
                for record in records:
                    write_fasta_record(
                        query_handle,
                        identifier=record.id,
                        description=record.description,
                        sequence=record.sequence,
                    )
            text_output = temporary_root / f"hmmscan.{shard_index:03d}.txt"
            shard_domtbl = temporary_root / f"pfam.{shard_index:03d}.domtblout"
            shard_domtbls.append(shard_domtbl)
            commands.append(
                [
                    config.hmmscan_command,
                    "--cpu",
                    "1",
                    "--cut_ga",
                    "--noali",
                    "-o",
                    str(text_output),
                    "--domtblout",
                    str(shard_domtbl),
                    str(hmm_path),
                    str(shard_fasta),
                ]
            )
        with ThreadPoolExecutor(max_workers=shard_count) as executor:
            futures = [
                executor.submit(
                    run_command,
                    command,
                    log_path=temporary_root / f"hmmscan.{index:03d}.log",
                )
                for index, command in enumerate(commands)
            ]
            for future in futures:
                future.result()

        temporary_domtbl = temporary_root / "pfam.domtblout"
        with atomic_text_writer(temporary_domtbl) as combined:
            for shard_domtbl in shard_domtbls:
                if not shard_domtbl.is_file():
                    raise MSSPackError(
                        f"hmmscan did not create its domtblout: {shard_domtbl}"
                    )
                with shard_domtbl.open("r", encoding="utf-8") as shard_handle:
                    for line in shard_handle:
                        combined.write(line)
        copy_or_decompress(temporary_domtbl, output_path)

    hit_count = 0
    query_ids: set[str] = set()
    with output_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            if not raw_line.startswith("#") and raw_line.strip():
                fields = raw_line.split()
                if len(fields) >= 4:
                    query_ids.add(fields[3])
                    hit_count += 1
    write_text(
        provenance_path,
        json.dumps(
            {
                **provenance,
                "hmmscan_command": config.hmmscan_command,
                "cutoff": "Pfam gathering thresholds (--cut_ga)",
                "queries_scanned": fallback_query_count,
                "similarity_assignments_skipped": len(similarity_assignments),
                "hmmscan_shards": shard_count,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
    )
    write_step_log(
        log_path=log_path,
        command=commands[0],
        step="functional-annotation-pfam-search",
        started_at=started_at,
        count_unit="domain hits",
        input_total=protein_count,
        changed_total=len(query_ids),
        output_total=hit_count,
        details=[
            "Threshold: model-specific Pfam gathering cutoff (--cut_ga)",
            f"Fallback queries scanned: {fallback_query_count:,}",
            f"Accepted similarity assignments skipped: {len(similarity_assignments):,}",
            f"Parallel hmmscan shards: {shard_count:,}",
            f"Queries with domain hits: {len(query_ids):,}",
        ],
    )
    write_step_metrics(
        metrics_path=metrics_path,
        step="functional-annotation-pfam-search",
        count_unit="domain hits",
        input_total=protein_count,
        changed_total=len(query_ids),
        output_total=hit_count,
        details={
            "queries_scanned": fallback_query_count,
            "similarity_assignments_skipped": len(similarity_assignments),
            "hmmscan_shards": shard_count,
            "queries_with_hits": len(query_ids),
            "cutoff": "cut_ga",
        },
    )


_CDD_DATA_FILE_SPECS = (
    ("cddid.tbl.gz", "cddid.tbl", "cddid", ".tbl"),
    ("cdtrack.txt", "cdtrack.txt", "cdtrack", ".txt"),
    (
        "family_superfamily_links",
        "family_superfamily_links",
        "family_superfamily_links",
        "",
    ),
    ("cddannot.dat.gz", "cddannot.dat", "cddannot", ".dat"),
    (
        "cddannot_generic.dat.gz",
        "cddannot_generic.dat",
        "cddannot_generic",
        ".dat",
    ),
    (
        "bitscore_specific.txt",
        "bitscore_specific.txt",
        "bitscore_specific",
        ".txt",
    ),
)


def _safe_extract_tar(archive_path: Path, destination: Path) -> None:
    destination_root = destination.resolve()
    with tarfile.open(archive_path, "r:gz") as archive:
        members = archive.getmembers()
        for member in members:
            target = (destination / member.name).resolve()
            if destination_root not in target.parents and target != destination_root:
                raise MSSPackError(f"CDD archive contains an unsafe path: {member.name}")
            if member.issym() or member.islnk():
                raise MSSPackError(f"CDD archive contains an unsupported link: {member.name}")
        archive.extractall(destination, members=members, filter="data")  # noqa: S202


def _find_cdd_database_prefix(directory: Path) -> Path:
    candidates = sorted(directory.rglob("Cdd.pal"))
    for candidate in candidates:
        prefix = candidate.with_suffix("")
        shard_aux = sorted(prefix.parent.glob(prefix.name + ".*.aux"))
        shard_rps = sorted(prefix.parent.glob(prefix.name + ".*.rps"))
        if shard_aux and shard_rps:
            return prefix
    candidates = sorted(directory.rglob("Cdd.aux"))
    if not candidates:
        raise MSSPackError(
            f"CDD archive did not contain a Cdd.pal or Cdd.aux database: {directory}"
        )
    prefix = candidates[0].with_suffix("")
    required_suffixes = (".aux", ".freq", ".loo", ".rps")
    missing = [suffix for suffix in required_suffixes if not Path(str(prefix) + suffix).is_file()]
    if missing:
        raise MSSPackError(
            f"CDD database prefix {prefix} is incomplete; missing: {', '.join(missing)}"
        )
    return prefix


def _cdd_database_prefix_complete(prefix: Path) -> bool:
    if Path(str(prefix) + ".pal").is_file():
        return bool(
            tuple(prefix.parent.glob(prefix.name + ".*.aux"))
            and tuple(prefix.parent.glob(prefix.name + ".*.rps"))
        )
    return all(
        Path(str(prefix) + suffix).is_file()
        for suffix in (".aux", ".freq", ".loo", ".rps")
    )


def _prepare_cdd_data_files(
    *,
    config: FunctionalAnnotationConfig,
    base_dir: Path,
    cache_dir: Path,
    database_root: Path,
    lock_settings: DatabaseLockSettings | None = None,
) -> tuple[Path, dict[str, object]]:
    staging_dir = ensure_dir(cache_dir / "data")
    versions_dir = ensure_dir(cache_dir / "data-versions")
    local_dir = (
        expand_path(config.cdd_data_dir, base_dir) if config.cdd_data_dir.strip() else None
    )
    provenance: dict[str, object] = {}
    materialized: list[tuple[str, Path, str]] = []
    for remote_name, local_name, stem, suffix in _CDD_DATA_FILE_SPECS:
        local_value = ""
        if local_dir is not None:
            source = local_dir / local_name
            if not source.is_file() and (local_dir / remote_name).is_file():
                source = local_dir / remote_name
            local_value = str(source)
        url = "" if local_dir is not None else f"{config.cdd_data_url.rstrip('/')}/{remote_name}"
        path, origin, digest, _source_provenance = _materialize_database_file(
            local_value=local_value,
            url=url,
            base_dir=base_dir,
            cache_dir=staging_dir,
            database_root=database_root,
            stem=stem,
            expected_suffix=suffix,
            content_addressed=True,
            lock_settings=lock_settings,
        )
        provenance[local_name] = {"origin": origin, "sha256": digest}
        materialized.append((local_name, path, digest))

    version_payload = [
        {"name": local_name, "sha256": digest}
        for local_name, _path, digest in materialized
    ]
    version_digest = hashlib.sha256(
        json.dumps(version_payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    version_dir = versions_dir / version_digest
    required_names = {local_name for _remote, local_name, _stem, _suffix in _CDD_DATA_FILE_SPECS}
    with acquire_database_lock(
        database_lock_path(database_root, f"cdd-data-version-{version_digest[:16]}"),
        label=f"CDD data version {version_digest[:16]}",
        settings=lock_settings,
    ):
        version_ready = (
            version_dir.is_dir()
            and (version_dir / "provenance.json").is_file()
            and all((version_dir / name).is_file() for name in required_names)
        )
        if not version_ready:
            temporary_dir = Path(
                tempfile.mkdtemp(prefix=f".{version_digest}.", dir=versions_dir)
            )
            try:
                for local_name, path, _digest in materialized:
                    link_or_copy(path, temporary_dir / local_name)
                write_text(
                    temporary_dir / "provenance.json",
                    json.dumps(
                        {
                            "version_sha256": version_digest,
                            "files": provenance,
                        },
                        indent=2,
                        sort_keys=True,
                    )
                    + "\n",
                )
                if version_dir.exists():
                    shutil.rmtree(version_dir)
                os.replace(temporary_dir, version_dir)
            finally:
                if temporary_dir.exists():
                    shutil.rmtree(temporary_dir)
    return version_dir, {
        "version_sha256": version_digest,
        "version_path": str(version_dir),
        "files": provenance,
    }


def _prepared_cdd_version(
    manifest: dict[str, object],
    versions_dir: Path,
) -> tuple[Path, str] | None:
    digest = manifest.get("archive_sha256")
    if (
        not isinstance(digest, str)
        or re.fullmatch(r"[0-9a-f]{64}", digest) is None
    ):
        return None
    version_dir = versions_dir / digest
    try:
        prefix = _find_cdd_database_prefix(version_dir)
    except MSSPackError:
        return None
    return prefix, digest


def _rpsblast_database_prefix(database_prefix: Path, temporary_root: Path) -> Path:
    resolved_prefix = database_prefix.resolve()
    if not any(character.isspace() for character in str(resolved_prefix)):
        return resolved_prefix
    alias_directory = temporary_root / "cdd-database"
    source_files = sorted(
        path
        for path in resolved_prefix.parent.glob(f"{resolved_prefix.name}.*")
        if path.is_file()
    )
    if not source_files:
        raise MSSPackError(f"CDD database prefix has no files: {resolved_prefix}")
    try:
        alias_directory.mkdir(parents=True)
        for source_file in source_files:
            link_or_copy(source_file, alias_directory / source_file.name)
    except OSError as exc:
        raise MSSPackError(
            "RPS-BLAST cannot use a CDD database path containing whitespace, and "
            f"msspack could not materialize a temporary whitespace-free alias for "
            f"{resolved_prefix.parent}: {exc}"
        ) from exc
    return alias_directory / resolved_prefix.name


def _prepare_cdd_database(
    *,
    config: FunctionalAnnotationConfig,
    base_dir: Path,
    cache_dir: Path,
    lock_settings: DatabaseLockSettings | None = None,
) -> tuple[Path, Path, dict[str, object]]:
    cdd_dir = ensure_dir(cache_dir / "cdd")
    data_dir, data_provenance = _prepare_cdd_data_files(
        config=config,
        base_dir=base_dir,
        cache_dir=cdd_dir,
        database_root=cache_dir,
        lock_settings=lock_settings,
    )
    if config.cdd_database.strip():
        configured = expand_path(config.cdd_database, base_dir)
        if configured.is_dir():
            prefix = _find_cdd_database_prefix(configured)
        else:
            prefix = configured
            if not _cdd_database_prefix_complete(prefix):
                raise MSSPackError(f"CDD database prefix is incomplete: {prefix}")
        return prefix, data_dir, {
            "origin": str(configured),
            "database_prefix": str(prefix),
            "data_files": data_provenance,
        }

    legacy_database_dir = cdd_dir / "db"
    legacy_manifest_path = cdd_dir / "database.json"
    legacy_manifest = _read_json(legacy_manifest_path)
    if legacy_manifest.get("origin") == config.cdd_url:
        try:
            legacy_prefix = _find_cdd_database_prefix(legacy_database_dir)
        except MSSPackError:
            pass
        else:
            return legacy_prefix, data_dir, {
                **legacy_manifest,
                "database_prefix": str(legacy_prefix),
                "data_files": data_provenance,
                "cache_layout": "legacy-read-only",
            }

    versions_dir = ensure_dir(cdd_dir / "database-versions")
    sources_dir = ensure_dir(cdd_dir / "sources")
    source_key = hashlib.sha256(config.cdd_url.encode("utf-8")).hexdigest()
    source_manifest_path = sources_dir / f"{source_key}.json"
    with acquire_database_lock(
        database_lock_path(cache_dir, f"cdd-source-{source_key[:16]}"),
        label=f"CDD source {source_key[:16]}",
        settings=lock_settings,
    ):
        previous = _read_json(source_manifest_path)
        prepared = (
            _prepared_cdd_version(previous, versions_dir)
            if previous.get("origin") == config.cdd_url
            else None
        )
        if prepared is None:
            archive_path = sources_dir / f".{source_key}.{os.getpid()}.tar.gz"
            _download(config.cdd_url, archive_path)
            archive_sha256 = _sha256(archive_path)
            try:
                version_dir = versions_dir / archive_sha256
                with acquire_database_lock(
                    database_lock_path(
                        cache_dir,
                        f"cdd-database-version-{archive_sha256[:16]}",
                    ),
                    label=f"CDD database version {archive_sha256[:16]}",
                    settings=lock_settings,
                ):
                    try:
                        prefix = _find_cdd_database_prefix(version_dir)
                    except MSSPackError:
                        temporary_dir = Path(
                            tempfile.mkdtemp(
                                prefix=f".{archive_sha256}.",
                                dir=versions_dir,
                            )
                        )
                        try:
                            _safe_extract_tar(archive_path, temporary_dir)
                            temporary_prefix = _find_cdd_database_prefix(temporary_dir)
                            relative_prefix = temporary_prefix.relative_to(temporary_dir)
                            if version_dir.exists():
                                shutil.rmtree(version_dir)
                            os.replace(temporary_dir, version_dir)
                            prefix = version_dir / relative_prefix
                        finally:
                            if temporary_dir.exists():
                                shutil.rmtree(temporary_dir)
            finally:
                archive_path.unlink(missing_ok=True)
            write_text(
                source_manifest_path,
                json.dumps(
                    {
                        "origin": config.cdd_url,
                        "archive_sha256": archive_sha256,
                        "database_prefix": str(prefix),
                        "cache_layout": "content-addressed",
                    },
                    indent=2,
                    sort_keys=True,
                )
                + "\n",
            )
        else:
            prefix, _archive_sha256 = prepared
    return prefix, data_dir, {
        **_read_json(source_manifest_path),
        "database_prefix": str(prefix),
        "data_files": data_provenance,
    }


def _read_cdd_metadata(path: Path) -> dict[str, tuple[str, str, str, int]]:
    metadata: dict[str, tuple[str, str, str, int]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            value = (
                row["accession"],
                row["name"],
                row["description"],
                int(row["model_length"]),
            )
            metadata[row["pssm_id"]] = value
            metadata[row["accession"]] = value
    return metadata


def _cdd_query_identifier(definition: str) -> str:
    identifier = definition.strip().split(None, 1)[0] if definition.strip() else ""
    return identifier.removeprefix("lcl|")


def _read_cdd_hits(result_path: Path, metadata_path: Path) -> dict[str, list[CddHit]]:
    metadata = _read_cdd_metadata(metadata_path)
    hits: dict[str, list[CddHit]] = defaultdict(list)
    query_aliases: dict[str, str] = {}
    in_domains = False
    with result_path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            fields = raw_line.rstrip("\n").split("\t")
            if not fields:
                continue
            token = fields[0]
            if token == "QUERY" and len(fields) >= 5:
                query_aliases[fields[1]] = _cdd_query_identifier(fields[4])
                continue
            if token == "DOMAINS":
                in_domains = True
                continue
            if token == "ENDDOMAINS":
                in_domains = False
                continue
            if not in_domains or len(fields) < 12:
                continue
            query_id = query_aliases.get(fields[1], fields[1])
            pssm_id = fields[3]
            accession, name, description, _model_length = metadata.get(
                pssm_id,
                metadata.get(fields[8], (fields[8], fields[9], fields[9], 0)),
            )
            try:
                hit = CddHit(
                    query_id=query_id,
                    hit_type=fields[2],
                    pssm_id=pssm_id,
                    accession=accession,
                    name=name,
                    description=description,
                    evalue=float(fields[6]),
                    bitscore=float(fields[7]),
                    incomplete=fields[10],
                )
            except ValueError as exc:
                raise MSSPackError("Could not parse numeric CDD rpsbproc fields") from exc
            hits[query_id].append(hit)
    return hits


def _write_cdd_hit_metadata(
    *,
    data_path: Path,
    result_path: Path,
    output_path: Path,
) -> int:
    needed: set[str] = set()
    in_domains = False
    with result_path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            fields = raw_line.rstrip("\n").split("\t")
            if fields and fields[0] == "DOMAINS":
                in_domains = True
                continue
            if fields and fields[0] == "ENDDOMAINS":
                in_domains = False
                continue
            if in_domains and len(fields) >= 9:
                needed.update((fields[3], fields[8]))
    written = 0
    with data_path.open("r", encoding="utf-8", errors="replace") as input_handle, atomic_text_writer(
        output_path
    ) as output_handle:
        output_handle.write("\t".join(CDD_METADATA_FIELDS) + "\n")
        for raw_line in input_handle:
            fields = raw_line.rstrip("\n").split("\t")
            if len(fields) < 5 or (fields[0] not in needed and fields[1] not in needed):
                continue
            output_handle.write("\t".join(fields[:5]) + "\n")
            written += 1
    return written


def write_empty_cdd_results(
    *,
    output_path: Path,
    metadata_path: Path,
    provenance_path: Path,
    log_path: Path,
    metrics_path: Path,
    reason: str = "CDD fallback is disabled.",
) -> None:
    started_at = datetime.now()
    write_text(output_path, "DATA\nENDDATA\n")
    write_text(metadata_path, "\t".join(CDD_METADATA_FIELDS) + "\n")
    write_text(provenance_path, '{"status": "skipped"}\n')
    write_step_log(
        log_path=log_path,
        command="msspack functional-annotation cdd-search --skipped",
        step="functional-annotation-cdd-search",
        started_at=started_at,
        status="skipped",
        count_unit="domain hits",
        input_total=0,
        changed_total=0,
        output_total=0,
        details=[reason],
    )
    write_step_metrics(
        metrics_path=metrics_path,
        step="functional-annotation-cdd-search",
        status="skipped",
        count_unit="domain hits",
        input_total=0,
        changed_total=0,
        output_total=0,
    )


def run_cdd_domain_search(
    *,
    protein_fasta_path: Path,
    output_path: Path,
    metadata_path: Path,
    provenance_path: Path,
    log_path: Path,
    metrics_path: Path,
    config: FunctionalAnnotationConfig,
    base_dir: Path,
    cache_dir: Path,
    lock_settings: DatabaseLockSettings | None = None,
    similarity_inputs: Sequence[tuple[Path, Path]] = (),
) -> None:
    started_at = datetime.now()
    protein_count = sum(1 for _ in iter_fasta(protein_fasta_path))
    if protein_count == 0:
        write_empty_cdd_results(
            output_path=output_path,
            metadata_path=metadata_path,
            provenance_path=provenance_path,
            log_path=log_path,
            metrics_path=metrics_path,
            reason="No translated protein sequences were available for CDD.",
        )
        return
    similarity_assignments = _combined_diamond_assignments(
        inputs=similarity_inputs,
        config=config,
    )
    fallback_records = [
        record
        for record in iter_fasta(protein_fasta_path)
        if record.id not in similarity_assignments
    ]
    if not fallback_records:
        write_empty_cdd_results(
            output_path=output_path,
            metadata_path=metadata_path,
            provenance_path=provenance_path,
            log_path=log_path,
            metrics_path=metrics_path,
            reason="All translated proteins already have accepted similarity assignments.",
        )
        return
    database_prefix, data_dir, provenance = _prepare_cdd_database(
        config=config,
        base_dir=base_dir,
        cache_dir=cache_dir,
        lock_settings=lock_settings,
    )
    with tempfile.TemporaryDirectory(prefix="msspack-rpsblast-") as temporary_dir:
        temporary_root = Path(temporary_dir)
        rpsblast_database_prefix = _rpsblast_database_prefix(
            database_prefix,
            temporary_root,
        )
        query_path = temporary_root / "cdd-fallback.fasta"
        with atomic_text_writer(query_path) as query_handle:
            for record in fallback_records:
                write_fasta_record(
                    query_handle,
                    identifier=record.id,
                    description=record.description,
                    sequence=record.sequence,
                )
        archive_path = temporary_root / "cdd.asn"
        raw_output = temporary_root / "cdd.rpsbproc.tsv"
        rpsblast_command = [
            config.rpsblast_command,
            "-query",
            str(query_path),
            "-db",
            str(rpsblast_database_prefix),
            "-evalue",
            f"{config.cdd_evalue:g}",
            "-outfmt",
            "11",
            "-num_threads",
            str(config.threads),
            "-mt_mode",
            "1",
            "-out",
            str(archive_path),
        ]
        run_command(
            rpsblast_command,
            log_path=log_path.with_name(f"{log_path.stem}.rpsblast.raw.log"),
        )
        if not archive_path.is_file():
            raise MSSPackError(f"RPS-BLAST did not create its archive output: {archive_path}")
        rpsbproc_command = [
            config.rpsbproc_command,
            "-i",
            str(archive_path),
            "-o",
            str(raw_output),
            "-e",
            f"{config.cdd_evalue:g}",
            "-m",
            "rep",
            "-t",
            "doms",
            "-d",
            str(data_dir),
            "-f",
            "-q",
        ]
        run_command(
            rpsbproc_command,
            log_path=log_path.with_name(f"{log_path.stem}.rpsbproc.raw.log"),
        )
        if not raw_output.is_file():
            raise MSSPackError(f"rpsbproc did not create its result file: {raw_output}")
        copy_or_decompress(raw_output, output_path)
        metadata_count = _write_cdd_hit_metadata(
            data_path=data_dir / "cddid.tbl",
            result_path=output_path,
            output_path=metadata_path,
        )

    cdd_hits = _read_cdd_hits(output_path, metadata_path)
    hit_count = sum(len(hits) for hits in cdd_hits.values())
    query_count = len(cdd_hits)
    write_text(
        provenance_path,
        json.dumps(
            {
                **provenance,
                "rpsblast_command": config.rpsblast_command,
                "rpsbproc_command": config.rpsbproc_command,
                "evalue": config.cdd_evalue,
                "redundancy": "representative",
                "queries_scanned": len(fallback_records),
                "similarity_assignments_skipped": len(similarity_assignments),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
    )
    write_step_log(
        log_path=log_path,
        command=rpsblast_command,
        step="functional-annotation-cdd-search",
        started_at=started_at,
        count_unit="domain hits",
        input_total=protein_count,
        changed_total=query_count,
        output_total=hit_count,
        details=[
            f"Fallback queries scanned: {len(fallback_records):,}",
            f"Accepted similarity assignments skipped: {len(similarity_assignments):,}",
            f"Queries with representative CDD hits: {query_count:,}",
            f"Matched CDD metadata records: {metadata_count:,}",
        ],
    )
    write_step_metrics(
        metrics_path=metrics_path,
        step="functional-annotation-cdd-search",
        count_unit="domain hits",
        input_total=protein_count,
        changed_total=query_count,
        output_total=hit_count,
        details={
            "queries_scanned": len(fallback_records),
            "similarity_assignments_skipped": len(similarity_assignments),
            "queries_with_hits": query_count,
            "metadata_records": metadata_count,
            "evalue": config.cdd_evalue,
            "redundancy": "representative",
        },
    )


def _read_diamond_metadata(path: Path) -> dict[tuple[str, str], DiamondMetadata]:
    metadata: dict[tuple[str, str], DiamondMetadata] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            raw_taxon_id = row.get("subject_taxon_id", "").strip()
            raw_taxonomy_weight = row.get("taxonomy_weight", "").strip()
            relation = row.get("taxonomy_relation", "").strip() or "unknown"
            metadata[(row["database"], row["subject_id"])] = DiamondMetadata(
                description=row["description"],
                source_weight=float(row["weight"]),
                subject_taxon_id=int(raw_taxon_id) if raw_taxon_id else None,
                subject_organism=row.get("subject_organism", "").strip(),
                taxonomy_relation=relation,
                taxonomy_weight=(
                    float(raw_taxonomy_weight)
                    if raw_taxonomy_weight
                    else RELATION_WEIGHTS["unknown"]
                ),
            )
    return metadata


def _read_diamond_hits(path: Path) -> dict[str, list[DiamondHit]]:
    hits: dict[str, list[DiamondHit]] = defaultdict(list)
    with path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            try:
                hit = DiamondHit(
                    query_id=row["qseqid"],
                    subject_id=row["sseqid"],
                    database=row["database"],
                    identity=float(row["pident"]),
                    query_length=int(row["qlen"]),
                    subject_length=int(row["slen"]),
                    query_coverage=float(row["qcovhsp"]),
                    subject_coverage=float(row["scovhsp"]),
                    evalue=float(row["evalue"]),
                    bitscore=float(row["bitscore"]),
                )
            except (KeyError, ValueError) as exc:
                raise MSSPackError(f"Could not parse DIAMOND result row: {row}") from exc
            hits[hit.query_id].append(hit)
    return hits


def _clean_description(description: str) -> str | None:
    value = re.sub(r"\s+", " ", description).strip().strip(".;")
    value = re.sub(r"^Cluster:\s*", "", value, flags=re.IGNORECASE)
    value = re.sub(r"^LOW QUALITY PROTEIN:\s*", "", value, flags=re.IGNORECASE)
    value = re.sub(r"^(?:RecName:\s*Full=|SubName:\s*Full=)", "", value)
    value = re.sub(r"\s*\{ECO:[^}]+\}", "", value)
    value = re.sub(r"\s*\((?:Fragment|Fragments)\)$", "", value, flags=re.IGNORECASE)
    value = value.strip().strip(".;")
    if not value or any(pattern.search(value) for pattern in _UNINFORMATIVE_DESCRIPTION_PATTERNS):
        return None
    return value


def _description_tokens(description: str) -> set[str]:
    return {
        token
        for token in re.findall(r"[A-Za-z0-9]+", description.lower())
        if len(token) > 1 and token not in _TOKEN_STOPWORDS
    }


def _taxonomy_adjusted_product(
    product: str,
    *,
    hit: DiamondHit,
    metadata: DiamondMetadata,
    config: FunctionalAnnotationConfig,
) -> tuple[str, str]:
    if (
        not config.taxonomy.enabled
        or metadata.taxonomy_relation not in {"cross_kingdom", "cross_domain"}
        or hit.identity >= config.taxonomy.distant_specificity_identity
    ):
        return product, ""
    adjusted = re.sub(
        r",?\s+(?:chloroplastic|cytosolic|mitochondrial|nuclear|peroxisomal|secreted)"
        r"(?:/[a-z]+)?$",
        "",
        product,
        flags=re.IGNORECASE,
    )
    adjusted = re.sub(
        r"\s+(?:isoform|homolog)\s+[A-Za-z]?\d+(?:[-.]\d+)*$",
        "",
        adjusted,
        flags=re.IGNORECASE,
    )
    adjusted = re.sub(
        r"\s+(?:family\s+)?member\s+[A-Za-z]?\d+(?:[-.]\d+)*$",
        " family protein",
        adjusted,
        flags=re.IGNORECASE,
    )
    adjusted = re.sub(r"\s+\d+(?:[-.]\d+)*$", "", adjusted)
    adjusted = re.sub(r"\s+", " ", adjusted).strip().strip(",;.")
    if not adjusted:
        return product, ""
    return adjusted, (
        "lineage-specific numbering/localization removed from a distant "
        f"{metadata.taxonomy_relation} hit below "
        f"{config.taxonomy.distant_specificity_identity:g}% identity"
        if adjusted != product
        else ""
    )


def _assignment_confidence(database: str, relation: str) -> str:
    if relation == "cross_domain":
        return "low"
    if database == "uniref90" or relation in {"cross_kingdom", "same_domain", "unknown"}:
        return "moderate"
    return "high"


def _diamond_assignments(
    *,
    hit_path: Path,
    metadata_path: Path,
    config: FunctionalAnnotationConfig,
) -> dict[str, Assignment]:
    metadata = _read_diamond_metadata(metadata_path)
    assignments: dict[str, Assignment] = {}
    for query_id, query_hits in _read_diamond_hits(hit_path).items():
        passing = [
            hit
            for hit in query_hits
            if hit.bitscore >= config.min_bitscore
            and hit.identity >= config.min_identity
            and hit.query_coverage >= config.min_query_coverage
            and hit.subject_coverage >= config.min_subject_coverage
            and hit.evalue <= config.evalue
            and metadata.get(
                (hit.database, hit.subject_id),
                DiamondMetadata("", 0.0),
            ).description
        ]
        if not passing:
            continue
        best_bitscore = max(hit.bitscore for hit in passing)
        near_top = [
            hit for hit in passing if hit.bitscore >= best_bitscore * config.near_top_bitscore_ratio
        ]
        candidates: list[DiamondCandidate] = []
        for hit in near_top:
            hit_metadata = metadata[(hit.database, hit.subject_id)]
            description = _clean_description(hit_metadata.description)
            if description is None:
                continue
            tokens = _description_tokens(description)
            if not tokens:
                continue
            coverage_factor = math.sqrt(
                max(0.0, hit.query_coverage / 100.0) * max(0.0, hit.subject_coverage / 100.0)
            )
            evidence_weight = (
                hit_metadata.source_weight
                * hit_metadata.taxonomy_weight
                * (hit.bitscore / best_bitscore)
                * coverage_factor
            )
            candidates.append(
                DiamondCandidate(
                    hit=hit,
                    description=description,
                    metadata=hit_metadata,
                    tokens=tokens,
                    evidence_weight=evidence_weight,
                )
            )
        total_weight = sum(item.evidence_weight for item in candidates)
        if not candidates or total_weight <= 0:
            continue
        token_weights: dict[str, float] = defaultdict(float)
        for candidate in candidates:
            for token in candidate.tokens:
                token_weights[token] += candidate.evidence_weight

        ranked: list[tuple[float, float, float, DiamondCandidate]] = []
        for candidate in candidates:
            consensus_score = sum(
                token_weights[token] / total_weight for token in candidate.tokens
            ) / len(candidate.tokens)
            ranked.append(
                (
                    consensus_score,
                    candidate.evidence_weight,
                    candidate.metadata.source_weight,
                    candidate,
                )
            )
        consensus_score, _evidence_weight, _database_weight, selected_candidate = max(
            ranked,
            key=lambda item: (item[0], item[1], item[2], item[3].hit.bitscore),
        )
        if consensus_score < config.min_token_score:
            continue
        selected = selected_candidate.hit
        selected_metadata = selected_candidate.metadata
        product, taxonomy_adjustment = _taxonomy_adjusted_product(
            selected_candidate.description,
            hit=selected,
            metadata=selected_metadata,
            config=config,
        )
        taxonomy_reason = (
            f"; taxonomy relation={selected_metadata.taxonomy_relation}, "
            f"weight={selected_metadata.taxonomy_weight:g}"
            if selected_metadata.taxonomy_relation
            else ""
        )
        assignments[query_id] = Assignment(
            product=product,
            source=selected.database,
            evidence_id=selected.subject_id,
            confidence=_assignment_confidence(
                selected.database,
                selected_metadata.taxonomy_relation,
            ),
            reason=(
                "AHRD-like weighted UniRef90 description consensus"
                if selected.database == "uniref90"
                else "AHRD-like weighted curated description consensus"
            )
            + taxonomy_reason,
            quality_code=(
                ("*" if selected.bitscore > 50 and selected.evalue < 1e-10 else "-")
                + ("*" if min(selected.query_coverage, selected.subject_coverage) > 60 else "-")
                + ("*" if consensus_score > 0.5 else "-")
            ),
            identity=selected.identity,
            query_coverage=selected.query_coverage,
            subject_coverage=selected.subject_coverage,
            evalue=selected.evalue,
            bitscore=selected.bitscore,
            consensus_score=consensus_score,
            subject_taxon_id=selected_metadata.subject_taxon_id,
            subject_organism=selected_metadata.subject_organism,
            taxonomy_relation=selected_metadata.taxonomy_relation,
            taxonomy_weight=selected_metadata.taxonomy_weight,
            taxonomy_adjustment=taxonomy_adjustment,
        )
    return assignments


def _combined_diamond_assignments(
    *,
    inputs: Sequence[tuple[Path, Path]],
    config: FunctionalAnnotationConfig,
) -> dict[str, Assignment]:
    assignments: dict[str, Assignment] = {}
    for hit_path, metadata_path in inputs:
        for query_id, assignment in _diamond_assignments(
            hit_path=hit_path,
            metadata_path=metadata_path,
            config=config,
        ).items():
            previous = assignments.get(query_id)
            if previous is None or _assignment_rank(assignment) > _assignment_rank(previous):
                assignments[query_id] = assignment
    return assignments


def _assignment_rank(assignment: Assignment) -> tuple[float, ...]:
    confidence_rank = {"high": 3.0, "moderate": 2.0, "low": 1.0}.get(
        assignment.confidence,
        0.0,
    )
    source_rank = {
        "existing": 5.0,
        "reference": 4.0,
        "swissprot": 3.0,
        "uniref90": 2.0,
        "pfam": 1.0,
        "cdd": 0.5,
    }.get(assignment.source.casefold(), 2.5)
    return (
        confidence_rank,
        assignment.taxonomy_weight or 0.0,
        source_rank,
        assignment.consensus_score or 0.0,
        min(assignment.query_coverage or 0.0, assignment.subject_coverage or 0.0),
        assignment.bitscore or 0.0,
    )


def _read_pfam_metadata(path: Path) -> dict[str, tuple[str, str]]:
    metadata: dict[str, tuple[str, str]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            name = row["name"]
            accession = row["accession"].split(".", 1)[0]
            value = (name, row["description"])
            metadata[name] = value
            if accession:
                metadata[accession] = value
    return metadata


def _read_pfam_hits(domtbl_path: Path, metadata_path: Path) -> dict[str, list[PfamHit]]:
    metadata = _read_pfam_metadata(metadata_path)
    hits: dict[str, list[PfamHit]] = defaultdict(list)
    with domtbl_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            if not raw_line.strip() or raw_line.startswith("#"):
                continue
            fields = raw_line.split(maxsplit=22)
            if len(fields) < 22:
                raise MSSPackError("Could not parse Pfam hmmscan domtblout row")
            name = fields[0]
            accession = fields[1].split(".", 1)[0] if fields[1] != "-" else ""
            meta_name, description = metadata.get(
                accession,
                metadata.get(name, (name, fields[22] if len(fields) > 22 else name)),
            )
            try:
                hit = PfamHit(
                    query_id=fields[3],
                    name=meta_name,
                    accession=accession or name,
                    description=description,
                    i_evalue=float(fields[12]),
                    bitscore=float(fields[13]),
                    hmm_from=int(fields[15]),
                    hmm_to=int(fields[16]),
                    hmm_length=int(fields[2]),
                )
            except ValueError as exc:
                raise MSSPackError("Could not parse numeric Pfam domtblout fields") from exc
            hits[hit.query_id].append(hit)
    return hits


def _domain_product(name: str, raw_description: str) -> str | None:
    description = re.sub(r"\s+", " ", raw_description).strip().strip(".;")
    if not description or any(
        pattern.search(description) for pattern in _UNINFORMATIVE_PFAM_PATTERNS
    ):
        return None
    duf = re.search(r"\bDUF\d+\b", f"{name} {description}", re.IGNORECASE)
    if duf:
        return f"{duf.group(0).upper()} domain-containing protein"
    if re.search(r"domain of unknown function", description, re.IGNORECASE):
        return None
    description = re.sub(r"^(?:The\s+)?", "", description).strip()
    if re.search(r"\bdomain$", description, re.IGNORECASE):
        return f"{description}-containing protein"
    if re.search(r"\bfamily$", description, re.IGNORECASE):
        return f"{description} protein"
    if re.search(r"\brepeat$", description, re.IGNORECASE):
        return f"{description}-containing protein"
    return f"{description} domain-containing protein"


def _pfam_product(hit: PfamHit) -> str | None:
    return _domain_product(hit.name, hit.description)


def _submission_safe_product(product: str) -> str:
    """Remove characters forbidden in DDBJ product qualifiers."""
    product = product.replace("\\", "/").replace('"', "")
    return re.sub(r"[\x00-\x20\x7f]+", " ", product).strip()


def _cdd_product(hit: CddHit) -> str | None:
    description = re.split(r"\.\s+", hit.description.strip(), maxsplit=1)[0]
    description = re.sub(r"\s*\[[^]]+]\s*$", "", description).strip().strip(".;")
    description = re.sub(
        r";\s*(?:Reviewed|Provisional)\s*$",
        "",
        description,
        flags=re.IGNORECASE,
    )
    if (
        not description
        or description.casefold() == "n/a"
        or any(pattern.search(description) for pattern in _UNINFORMATIVE_DESCRIPTION_PATTERNS)
        or re.search(r"\bfunction unknown\b", description, re.IGNORECASE)
    ):
        description = ""
    product = _domain_product(hit.name, description) if description else None
    if product is not None and len(product) <= 120:
        return product

    short_name = re.sub(r"[_-]+", " ", hit.name).strip()
    if not short_name or re.fullmatch(
        r"(?:cl|cd|pfam|smart|COG|KOG|PRK)\d+",
        short_name,
        flags=re.IGNORECASE,
    ):
        return None
    product = _domain_product(short_name, short_name)
    return product if product is not None and len(product) <= 120 else None


def _pfam_assignments(
    *,
    domtbl_path: Path,
    metadata_path: Path,
    config: FunctionalAnnotationConfig,
) -> dict[str, Assignment]:
    assignments: dict[str, Assignment] = {}
    for query_id, hits in _read_pfam_hits(domtbl_path, metadata_path).items():
        accepted = [
            (hit, _pfam_product(hit))
            for hit in hits
            if hit.i_evalue <= config.pfam_max_i_evalue
            and hit.domain_coverage >= config.pfam_min_domain_coverage
        ]
        informative = [(hit, product) for hit, product in accepted if product is not None]
        if not informative:
            continue
        hit, product = max(
            informative,
            key=lambda item: (item[0].bitscore, -item[0].i_evalue, item[0].domain_coverage),
        )
        assert product is not None
        assignments[query_id] = Assignment(
            product=product,
            source="pfam",
            evidence_id=hit.accession,
            confidence="moderate",
            reason="informative Pfam domain fallback",
            evalue=hit.i_evalue,
            bitscore=hit.bitscore,
            consensus_score=hit.domain_coverage,
        )
    return assignments


def _cdd_assignments(
    *,
    result_path: Path,
    metadata_path: Path,
    config: FunctionalAnnotationConfig,
) -> dict[str, Assignment]:
    assignments: dict[str, Assignment] = {}
    type_rank = {"Specific": 2, "Superfamily": 1}
    for query_id, hits in _read_cdd_hits(result_path, metadata_path).items():
        informative = [
            (hit, _cdd_product(hit))
            for hit in hits
            if hit.hit_type in type_rank and hit.evalue <= config.cdd_evalue
        ]
        accepted = [(hit, product) for hit, product in informative if product is not None]
        if not accepted:
            continue
        hit, product = max(
            accepted,
            key=lambda item: (
                type_rank[item[0].hit_type],
                item[0].incomplete == "-",
                item[0].bitscore,
                -item[0].evalue,
            ),
        )
        assert product is not None
        assignments[query_id] = Assignment(
            product=product,
            source="cdd",
            evidence_id=hit.accession,
            confidence="moderate" if hit.hit_type == "Specific" else "low",
            reason=f"informative CDD {hit.hit_type.casefold()} domain fallback",
            evalue=hit.evalue,
            bitscore=hit.bitscore,
        )
    return assignments


def _format_number(value: float | None) -> str:
    return "" if value is None else f"{value:.6g}"


def _logged_duration_seconds(path: Path) -> float | None:
    if not path.is_file():
        return None
    matches = re.findall(
        r"^Duration seconds:\s*([0-9.]+)\s*$",
        path.read_text(encoding="utf-8", errors="replace"),
        flags=re.MULTILINE,
    )
    return float(matches[-1]) if matches else None


def write_domain_search_comparison(
    *,
    pfam_hits_path: Path,
    pfam_metadata_path: Path,
    pfam_metrics_path: Path,
    pfam_log_path: Path,
    cdd_hits_path: Path,
    cdd_metadata_path: Path,
    cdd_metrics_path: Path,
    cdd_log_path: Path,
    output_path: Path,
    config: FunctionalAnnotationConfig,
) -> None:
    rows: list[tuple[str, str, Path, Path, int]] = [
        (
            "pfam",
            "HMMER hmmscan",
            pfam_metrics_path,
            pfam_log_path,
            len(
                _pfam_assignments(
                    domtbl_path=pfam_hits_path,
                    metadata_path=pfam_metadata_path,
                    config=config,
                )
            ),
        ),
        (
            "cdd",
            "NCBI RPS-BLAST + rpsbproc",
            cdd_metrics_path,
            cdd_log_path,
            len(
                _cdd_assignments(
                    result_path=cdd_hits_path,
                    metadata_path=cdd_metadata_path,
                    config=config,
                )
            ),
        ),
    ]
    with atomic_text_writer(output_path) as output:
        output.write(
            "backend\tengine\tstatus\tqueries_scanned\tqueries_with_hits\t"
            "domain_hits\tinformative_assignments\tduration_seconds\t"
            "queries_per_second\n"
        )
        for backend, engine, metrics_path, log_path, informative_count in rows:
            metrics = _read_json(metrics_path)
            details = metrics.get("details", {})
            detail_map = details if isinstance(details, dict) else {}
            scanned = int(detail_map.get("queries_scanned", 0))
            queries_with_hits = int(detail_map.get("queries_with_hits", 0))
            raw_domain_hits = metrics.get("output_total", 0)
            domain_hits = (
                raw_domain_hits
                if isinstance(raw_domain_hits, int) and not isinstance(raw_domain_hits, bool)
                else 0
            )
            duration = _logged_duration_seconds(log_path)
            rate = scanned / duration if duration and duration > 0 else None
            output.write(
                "\t".join(
                    (
                        backend,
                        engine,
                        str(metrics.get("status", "unknown")),
                        str(scanned),
                        str(queries_with_hits),
                        str(domain_hits),
                        str(informative_count),
                        _format_number(duration),
                        _format_number(rate),
                    )
                )
                + "\n"
            )


def apply_functional_annotations(
    *,
    annotation_table_path: Path,
    diamond_hits_path: Path,
    diamond_metadata_path: Path,
    uniref90_hits_path: Path | None = None,
    uniref90_metadata_path: Path | None = None,
    pfam_hits_path: Path,
    pfam_metadata_path: Path,
    cdd_hits_path: Path | None = None,
    cdd_metadata_path: Path | None = None,
    output_path: Path,
    evidence_path: Path,
    log_path: Path,
    metrics_path: Path,
    config: FunctionalAnnotationConfig,
    missing_product: str = "hypothetical protein",
    domain_comparison_path: Path | None = None,
    pfam_search_metrics_path: Path | None = None,
    pfam_search_log_path: Path | None = None,
    cdd_search_metrics_path: Path | None = None,
    cdd_search_log_path: Path | None = None,
    taxonomy_context_path: Path | None = None,
    name_standardization_summary_path: Path | None = None,
) -> None:
    started_at = datetime.now()
    diamond_inputs = [(diamond_hits_path, diamond_metadata_path)]
    if uniref90_hits_path is not None and uniref90_metadata_path is not None:
        diamond_inputs.append((uniref90_hits_path, uniref90_metadata_path))
    diamond = _combined_diamond_assignments(
        inputs=diamond_inputs,
        config=config,
    )
    pfam = _pfam_assignments(
        domtbl_path=pfam_hits_path,
        metadata_path=pfam_metadata_path,
        config=config,
    )
    cdd = (
        _cdd_assignments(
            result_path=cdd_hits_path,
            metadata_path=cdd_metadata_path,
            config=config,
        )
        if cdd_hits_path is not None and cdd_metadata_path is not None
        else {}
    )
    row_count = 0
    updated = 0
    preserved = 0
    hypothetical = 0
    functional_assignments = 0
    standardized_rows = 0
    warning_rows = 0
    source_counts: dict[str, int] = defaultdict(int)
    standardization_counts: Counter[str] = Counter()
    warning_counts: Counter[str] = Counter()
    naming_context = ProductNameContext.from_taxonomy_context(
        load_taxonomy_context(taxonomy_context_path)
    )

    with annotation_table_path.open("r", encoding="utf-8", newline="") as input_handle:
        reader = csv.DictReader(input_handle, delimiter="\t")
        if reader.fieldnames != ["ID", "Description", "Locus_tag"]:
            raise MSSPackError("Annotation table must have columns: ID, Description, Locus_tag")
        with (
            atomic_text_writer(output_path) as output_handle,
            atomic_text_writer(evidence_path) as evidence_handle,
        ):
            writer = csv.DictWriter(
                output_handle,
                fieldnames=reader.fieldnames,
                delimiter="\t",
                lineterminator="\n",
            )
            evidence_writer = csv.DictWriter(
                evidence_handle,
                fieldnames=EVIDENCE_FIELDS,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            evidence_writer.writeheader()
            for row in reader:
                row_count += 1
                identifier = row["ID"]
                original = row["Description"]
                is_missing = original.casefold() in {
                    "hypothetical protein",
                    "uncharacterized protein",
                    missing_product.casefold(),
                }
                assignment: Assignment | None = None
                if config.overwrite_existing or is_missing:
                    assignment = (
                        diamond.get(identifier)
                        or pfam.get(identifier)
                        or cdd.get(identifier)
                    )
                if assignment is not None:
                    functional_assignments += 1
                elif not is_missing:
                    assignment = Assignment(
                        product=original,
                        source="existing",
                        evidence_id="",
                        confidence="existing",
                        reason=(
                            "no accepted replacement; existing product preserved"
                            if config.overwrite_existing
                            else "existing non-hypothetical product preserved"
                        ),
                    )
                    preserved += 1
                else:
                    assignment = Assignment(
                        product=original,
                        source="none",
                        evidence_id="",
                        confidence="none",
                        reason="no acceptable similarity or informative domain evidence",
                    )
                candidate_product = assignment.product
                candidate_source = assignment.source
                standardization = standardize_product_name(
                    candidate_product,
                    source=candidate_source,
                    evidence_id=assignment.evidence_id,
                    subject_organism=assignment.subject_organism,
                    context=naming_context,
                )
                row["Description"] = standardization.product
                effective_source = (
                    candidate_source if standardization.informative else "none"
                )
                effective_confidence = (
                    assignment.confidence if standardization.informative else "none"
                )
                effective_reason = assignment.reason
                if not standardization.informative and candidate_source != "none":
                    effective_reason += (
                        "; candidate name rejected by product-name standardization"
                    )
                updated += int(row["Description"] != original)
                standardized_rows += int(bool(standardization.actions))
                warning_rows += int(bool(standardization.warnings))
                hypothetical += int(not standardization.informative)
                source_counts[effective_source] += 1
                standardization_counts.update(standardization.actions)
                warning_counts.update(standardization.warnings)
                writer.writerow(row)
                evidence_writer.writerow(
                    {
                        "ID": identifier,
                        "Locus_tag": row["Locus_tag"],
                        "original_product": original,
                        "proposed_product": candidate_product,
                        "assigned_product": row["Description"],
                        "candidate_source": candidate_source,
                        "source": effective_source,
                        "evidence_id": assignment.evidence_id,
                        "confidence": effective_confidence,
                        "name_standardization": ",".join(standardization.actions),
                        "name_warnings": ",".join(standardization.warnings),
                        "quality_code": assignment.quality_code,
                        "reason": effective_reason,
                        "identity": _format_number(assignment.identity),
                        "query_coverage": _format_number(assignment.query_coverage),
                        "subject_coverage": _format_number(assignment.subject_coverage),
                        "evalue": _format_number(assignment.evalue),
                        "bitscore": _format_number(assignment.bitscore),
                        "consensus_score": _format_number(assignment.consensus_score),
                        "subject_taxon_id": (
                            str(assignment.subject_taxon_id)
                            if assignment.subject_taxon_id is not None
                            else ""
                        ),
                        "subject_organism": assignment.subject_organism,
                        "taxonomy_relation": assignment.taxonomy_relation,
                        "taxonomy_weight": _format_number(assignment.taxonomy_weight),
                        "taxonomy_adjustment": assignment.taxonomy_adjustment,
                    }
                )

    details = [
        f"Rows changed in the final annotation table: {updated:,}",
        f"Rows with functional-evidence candidates: {functional_assignments:,}",
        f"Existing products retained as the candidate source: {preserved:,}",
        f"Rows standardized before family-name consistency: {standardized_rows:,}",
        f"Rows with residual naming warnings: {warning_rows:,}",
        f"Rows remaining without accepted annotation: {hypothetical:,}",
    ]
    details.extend(
        f"Assignments from {source}: {count:,}" for source, count in sorted(source_counts.items())
    )
    write_step_log(
        log_path=log_path,
        command=(
            "msspack functional-annotation assign-and-standardize-products "
            f"--input {annotation_table_path} --output {output_path}"
        ),
        step="functional-annotation-assign-and-standardize-products",
        started_at=started_at,
        count_unit="annotation rows",
        input_total=row_count,
        changed_total=updated,
        output_total=row_count,
        details=details,
    )
    write_step_metrics(
        metrics_path=metrics_path,
        step="functional-annotation-assign-and-standardize-products",
        count_unit="annotation rows",
        input_total=row_count,
        changed_total=updated,
        output_total=row_count,
        details={
            "existing_products_preserved": preserved,
            "functional_evidence_candidates": functional_assignments,
            "standardized_rows": standardized_rows,
            "residual_warning_rows": warning_rows,
            "standardization_actions": dict(sorted(standardization_counts.items())),
            "naming_warnings": dict(sorted(warning_counts.items())),
            "naming_context": {
                "domain": naming_context.domain,
                "kingdom": naming_context.kingdom,
            },
            "unannotated_rows": hypothetical,
            "source_counts": dict(sorted(source_counts.items())),
            "evidence_path": str(evidence_path),
        },
    )
    if name_standardization_summary_path is not None:
        with atomic_text_writer(name_standardization_summary_path) as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=("kind", "code", "count"),
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(
                {"kind": "action", "code": code, "count": count}
                for code, count in sorted(standardization_counts.items())
            )
            writer.writerows(
                {"kind": "warning", "code": code, "count": count}
                for code, count in sorted(warning_counts.items())
            )
    if (
        domain_comparison_path is not None
        and pfam_search_metrics_path is not None
        and pfam_search_log_path is not None
        and cdd_search_metrics_path is not None
        and cdd_search_log_path is not None
        and cdd_hits_path is not None
        and cdd_metadata_path is not None
    ):
        write_domain_search_comparison(
            pfam_hits_path=pfam_hits_path,
            pfam_metadata_path=pfam_metadata_path,
            pfam_metrics_path=pfam_search_metrics_path,
            pfam_log_path=pfam_search_log_path,
            cdd_hits_path=cdd_hits_path,
            cdd_metadata_path=cdd_metadata_path,
            cdd_metrics_path=cdd_search_metrics_path,
            cdd_log_path=cdd_search_log_path,
            output_path=domain_comparison_path,
            config=config,
        )
