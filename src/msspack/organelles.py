"""Prepare a separate MISC submission, preserving organelle sequences and annotations."""
from __future__ import annotations

import hashlib
import html
import json
import re
import tempfile
from collections import Counter
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any, cast

from Bio.SeqFeature import AfterPosition, BeforePosition, ExactPosition, SimpleLocation

from . import __version__
from .config import MSSPackConfig, load_config
from .execution import implementation_fingerprint
from .fasta import write_wrapped_sequence
from .mss_converter.render import render_source_feature
from .organelle_assessment import render_assessment_html
from .organelle_inputs import feature_qualifiers, load_organelle_records
from .organelle_plots import write_organelle_gene_flows
from .output_state import locked_output, publish_submission
from .submission_render import render_header
from .utils import MSSPackError, ensure_dir, expand_path, write_text
from .validation import (
    ValidationArtifacts,
    ValidationOptions,
    ValidationSummary,
    load_validation_summary,
    run_validation,
    write_validation_not_run_summary,
)


@dataclass(frozen=True)
class OrganelleArtifacts:
    final: Path
    ann_path: Path
    fasta_path: Path
    summary_path: Path
    report_path: Path


def _digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _text(value: str) -> str:
    if any(character in value for character in "\x00\t\r\n"):
        raise MSSPackError("Organelle annotation values must not contain control characters")
    return value


def _part_location(part: Any, length: int) -> str:
    if part.ref or part.ref_db:
        raise MSSPackError("Remote organelle locations need curator review; no feature was discarded")
    if not 0 <= int(part.start) < int(part.end) <= length:
        raise MSSPackError(f"Organelle feature location outside sequence: {part}")
    for position in (part.start, part.end):
        if type(position) not in {ExactPosition, BeforePosition, AfterPosition}:
            raise MSSPackError(f"Unsupported fuzzy organelle location: {part}")
    start = ("<" if isinstance(part.start, BeforePosition) else
             ">" if isinstance(part.start, AfterPosition) else "") + str(int(part.start) + 1)
    end = ("<" if isinstance(part.end, BeforePosition) else
           ">" if isinstance(part.end, AfterPosition) else "") + str(int(part.end))
    return f"{start}..{end}"


def _location(feature: Any, length: int) -> str:
    if feature.location is None:
        raise MSSPackError(f"Organelle feature without location: {feature.type}")
    parts = list(feature.location.parts)
    if any(part.strand not in {1, -1} for part in parts):
        raise MSSPackError(f"Organelle feature has unknown strand: {feature.type}")
    operator = getattr(feature.location, "operator", "join")
    if operator not in {"join", "order"}:
        raise MSSPackError(f"Unsupported organelle location operator: {operator}")
    if operator == "join":
        merged: list[Any] = []
        for part in parts:
            if merged and merged[-1].strand == part.strand:
                previous = merged[-1]
                if (part.strand == 1 and type(previous.end) is ExactPosition
                        and type(part.start) is ExactPosition and previous.end == part.start):
                    merged[-1] = SimpleLocation(previous.start, part.end, strand=1)
                    continue
                if (part.strand == -1 and type(previous.start) is ExactPosition
                        and type(part.end) is ExactPosition and previous.start == part.end):
                    merged[-1] = SimpleLocation(part.start, previous.end, strand=-1)
                    continue
            merged.append(part)
        parts = merged
    if all(part.strand == -1 for part in parts):
        segments = [_part_location(part, length) for part in reversed(parts)]
        inner = segments[0] if len(segments) == 1 else f"{operator}({','.join(segments)})"
        return f"complement({inner})"
    segments = [f"complement({_part_location(part, length)})" if part.strand == -1
                else _part_location(part, length) for part in parts]
    return segments[0] if len(segments) == 1 else f"{operator}({','.join(segments)})"


def _gene_for(record: Any, feature: Any, qualifiers: dict[str, list[str]]) -> Any:
    if feature.type == "gene":
        return feature
    name = qualifiers.get("gene", [""])[0]
    tag = qualifiers.get("locus_tag", [""])[0]
    if tag:
        linked = [gene for gene in record.features if gene.type == "gene"
                  and gene.qualifiers.get("locus_tag") == [tag]]
        if len(linked) > 1:
            raise MSSPackError(f"Ambiguous organelle gene locus_tag: {record.id}/{tag}")
        if linked:
            if name and linked[0].qualifiers.get("gene", [name]) != [name]:
                raise MSSPackError(f"Conflicting organelle gene names for locus_tag: {tag}")
            return linked[0]
    genes = [gene for gene in record.features if gene.type == "gene"
             and (not name or gene.qualifiers.get("gene", [""])[0] == name)
             and (not tag or not gene.qualifiers.get("locus_tag")
                  or gene.qualifiers["locus_tag"] == [tag])
             and gene.location is not None and feature.location is not None
             and all(any(g.strand == p.strand and int(g.start) <= int(p.start) and int(p.end) <= int(g.end)
                         for g in gene.location.parts) for p in feature.location.parts)]
    if len(genes) > 1:
        raise MSSPackError(f"Ambiguous organelle gene association: {record.id}/{name}")
    return genes[0] if genes else None


def _gene_identity(record: Any, feature: Any, qualifiers: dict[str, list[str]]) -> str:
    if str(feature.id) != "<unknown id>":
        return str(feature.id)
    gene = _gene_for(record, feature, qualifiers)
    name = (gene.qualifiers if gene else qualifiers).get("gene", [""])[0]
    return f"{name}/{gene.location if gene else feature.location}"


def _render_features(
    record: Any, *, item_name: str, code: str, prefix: str, tags: dict[str, str],
    owners: dict[str, str], reserved_tags: set[str], audit: list[dict[str, object]],
    counts: Counter[str], warnings: list[str],
) -> str:
    chunks: list[str] = []
    for index, feature in enumerate(record.features):
        if feature.type == "source":
            continue
        key = _text(str(feature.type))
        if not re.fullmatch(r"[A-Za-z0-9_'\-]+", key):
            raise MSSPackError(f"Invalid organelle feature key: {key}")
        location = _location(feature, len(record.seq))
        qualifiers = feature_qualifiers(feature)
        original = {qualifier: list(values) for qualifier, values in qualifiers.items()}
        if key == "gene":
            related = [child for child in record.features
                       if child.type in {"CDS", "mRNA", "tRNA", "rRNA", "ncRNA", "misc_RNA"}
                       and _gene_for(record, child, feature_qualifiers(child)) is feature]
            if related:
                audit.append({"organelle": item_name, "source_record": str(record.id),
                              "feature_index": index, "type": key, "location": location,
                              "input_qualifiers": original, "output_qualifiers": {},
                              "representation": "gene qualifiers on related CDS/RNA features",
                              "feature_sequence_sha256": _digest(str(feature.extract(record.seq)).encode())})
                continue
            key = "misc_feature"
            qualifiers.setdefault("note", []).append("gene annotation without a CDS or RNA feature")
        elif key in {"CDS", "mRNA", "tRNA", "rRNA", "ncRNA", "misc_RNA"}:
            gene = _gene_for(record, feature, qualifiers)
            if gene is not None:
                for qualifier in ("gene", "gene_synonym", "locus_tag", "pseudo", "pseudogene", "note"):
                    for value in feature_qualifiers(gene).get(qualifier, []):
                        if value not in qualifiers.setdefault(qualifier, []):
                            qualifiers[qualifier].append(value)
                    if not qualifiers.get(qualifier):
                        qualifiers.pop(qualifier, None)
        # DDBJ requires the supplied protein translation for biological exceptions.
        # Ordinary translations are recomputed by transChecker and retained in the audit.
        if "exception" not in qualifiers:
            qualifiers.pop("translation", None)
        else:
            warnings.append(
                f"{item_name}/{record.id}/{index}: biological exception "
                f"{','.join(qualifiers['exception'])}; review supplied translation and evidence"
            )
        if key == "CDS":
            if "trans_splicing" in qualifiers:
                warnings.append(f"{item_name}/{record.id}/{index}: trans_splicing model retained; review transcript evidence")
            if any(type(part.start) is not ExactPosition or type(part.end) is not ExactPosition
                   for part in feature.location.parts):
                warnings.append(f"{item_name}/{record.id}/{index}: partial CDS boundary retained")
            supplied = qualifiers.get("transl_table", [])
            if supplied and code and supplied != [code]:
                raise MSSPackError(f"Conflicting organelle CDS translation tables: {record.id}")
            table = supplied[0] if supplied else code
            if not table:
                raise MSSPackError(f"Missing genetic code for organelle CDS: {record.id}")
            qualifiers["transl_table"] = [table]
            qualifiers.setdefault("codon_start", ["1"])
            if "exception" not in qualifiers and not ({"pseudo", "pseudogene"} & set(qualifiers)):
                phase = int(qualifiers["codon_start"][0]) - 1
                if phase not in {0, 1, 2}:
                    raise MSSPackError(f"Invalid organelle codon_start: {record.id}")
                sequence = feature.extract(record.seq)[phase:]
                if len(sequence) % 3:
                    warnings.append(f"{item_name}/{record.id}/{index}: CDS length is not divisible by 3")
                protein = str(sequence[:len(sequence) // 3 * 3].translate(table=int(table)))
                if "*" in protein[:-1]:
                    warnings.append(f"{item_name}/{record.id}/{index}: internal stop codon; review annotation")
        if key in {"CDS", "mRNA", "tRNA", "rRNA", "ncRNA", "misc_RNA"} or feature.type == "gene":
            identity = _gene_identity(record, feature, qualifiers)
            tag_key = f"{item_name}/{record.id}/{identity}"
            # Preserve annotated locus tags; assign new ones consistently within each gene.
            if "locus_tag" not in qualifiers:
                if tag_key not in tags:
                    counter = len(tags) + 1
                    while (candidate := f"{prefix}_{item_name}{counter:06d}") in reserved_tags:
                        counter += 1
                    tags[tag_key] = candidate
                    reserved_tags.add(candidate)
                qualifiers["locus_tag"] = [tags[tag_key]]
            for tag in qualifiers["locus_tag"]:
                if tag in owners and owners[tag] != tag_key:
                    raise MSSPackError(f"Conflicting organelle locus_tag: {tag}")
                owners[tag] = tag_key
        qualifier_rows: list[tuple[str, str]] = []
        for qualifier, values in qualifiers.items():
            qualifier = _text(qualifier)
            if not re.fullmatch(r"[a-z][a-z_]*", qualifier):
                raise MSSPackError(f"Invalid organelle qualifier: {qualifier}")
            for value in values or [""]:
                qualifier_rows.append((qualifier, _text(value)))
        if qualifier_rows:
            first_qualifier, first_value = qualifier_rows[0]
            chunks.append(f"\t{key}\t{location}\t{first_qualifier}\t{first_value}\n")
            chunks.extend(f"\t\t\t{qualifier}\t{value}\n" for qualifier, value in qualifier_rows[1:])
        else:
            chunks.append(f"\t{key}\t{location}\t\t\n")
        counts[key] += 1
        audit.append({"organelle": item_name, "source_record": str(record.id), "feature_index": index,
                      "type": key, "location": location, "input_qualifiers": original,
                      "output_qualifiers": qualifiers,
                      "feature_sequence_sha256": _digest(str(feature.extract(record.seq)).encode())})
    return "".join(chunks)


def _write_report(
    work: Path, records: list[dict[str, Any]], warnings: list[str], validation: ValidationSummary | None,
    flows: list[dict[str, object]],
) -> None:
    from reportlab.lib.colors import HexColor  # type: ignore[import-untyped]
    from reportlab.pdfgen.canvas import Canvas  # type: ignore[import-untyped]

    height = 80 + 50 * len(records)
    canvas = Canvas(str(work / "sequences.pdf"), pagesize=(780, height))
    canvas.setFont("Helvetica", 15)
    canvas.drawString(12, height - 24, "Organelle submission sequences (bp)")
    svg_parts = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 780 {height}">',
                 '<text x="12" y="24" font-family="sans-serif" font-size="15">'
                 'Organelle submission sequences (bp)</text>']
    maximum = max(row["length"] for row in records)
    for index, row in enumerate(records):
        y = height - 70 - 50 * index
        width = 390 * row["length"] / maximum
        canvas.setFont("Helvetica", 10)
        canvas.drawString(12, y + 6, row["entry"])
        canvas.setFillColor(HexColor("#167c80"))
        canvas.rect(220, y, width, 20, stroke=0, fill=1)
        canvas.setFillColor(HexColor("#000000"))
        canvas.drawString(620, y + 6, f"{row['length']:,}")
        svg_parts.extend([
            f'<text x="12" y="{height-y-6}" font-family="sans-serif" font-size="10">{html.escape(row["entry"])}</text>',
            f'<rect x="220" y="{height-y-20}" width="{width}" height="20" fill="#167c80"/>',
            f'<text x="620" y="{height-y-6}" font-family="sans-serif" font-size="10">{row["length"]:,}</text>',
        ])
    canvas.save()
    svg_parts.append("</svg>")
    svg = "".join(svg_parts)
    write_text(work / "sequences.svg", svg)
    columns = ("entry", "organelle", "length", "topology", "completeness", "feature_counts")
    rows = "".join("<tr>" + "".join(f"<td>{html.escape(str(row[k]))}</td>" for k in columns)
                   + "</tr>" for row in records)
    notes = "".join(f"<li>{html.escape(note)}</li>" for note in warnings) or "<li>No annotation warnings.</li>"
    checks_html = "".join(
        f"<tr><td>{html.escape(check.label)}</td><td>{check.status}</td>"
        f"<td>{check.error_count}</td><td>{check.warning_count}</td>"
        f"<td>{html.escape(str(check.record_counts or {}))}</td></tr>"
        for check in validation.checks
    ) if validation else ""
    flow_html = "<h2>Pipeline gene flow</h2>" + "".join(
        "<section><h3>" + html.escape(str(flow["title"])) + "</h3>"
        + (work / str(flow["svg"])).read_text(encoding="utf-8")
        + render_assessment_html(cast(dict[str, Any], flow["assessment"]), work) + "</section>" for flow in flows
    )
    write_text(work / "index.html", "<!doctype html><meta charset='utf-8'><title>Organelle MSS review</title>"
               "<style>body{font:16px system-ui;margin:32px;max-width:1200px}"
               "table{border-collapse:collapse}td,th{border:1px solid #bbb;padding:8px}"
               "svg{max-width:100%;height:auto}</style><h1>Organelle MSS candidate review</h1>"
               "<p>MISC submission; nuclear analysis is separate. No external submission performed.</p>"
               + svg + "<table><tr>" + "".join(f"<th>{c}</th>" for c in columns) + "</tr>"
               + rows + "</table>" + flow_html + "<h2>Official DDBJ validation</h2><table>"
               "<tr><th>Tool</th><th>Status</th><th>Errors</th><th>Warnings</th><th>Record counts</th></tr>"
               + checks_html + "</table><h2>Review notes</h2><ul>" + notes + "</ul>")


def _prepare(
    config: MSSPackConfig, *, validate: bool, nuclear_annotation: Path | None,
) -> OrganelleArtifacts:
    if not config.organelles:
        raise MSSPackError("No [[organelles]] inputs configured")
    root = ensure_dir(config.output_dir / "organelles")
    work = Path(tempfile.mkdtemp(prefix="candidate-", dir=ensure_dir(root / ".msspack-work")))
    ann = work / f"{config.sample.locus_tag}.organelles.ann.txt"
    fasta = work / f"{config.sample.locus_tag}.organelles.fasta"
    misc = replace(config, submission=replace(config.submission, datatype="MISC", keywords=[]))
    # Nuclear assembly method/coverage/name do not describe the organelle assemblies.
    header_lines = [line for line in render_header(misc).split("\tST_COMMENT\t", 1)[0].splitlines()
                    if line.split("\t")[1] != "DATATYPE"]
    first = header_lines[0].split("\t")
    first[0] = "COMMON"
    header_lines[0] = "\t".join(first)
    header = "\n".join(header_lines) + "\n"
    records_summary: list[dict[str, Any]] = []
    audit: list[dict[str, object]] = []
    warnings: list[str] = []
    inputs: list[dict[str, object]] = []
    tags: dict[str, str] = {}
    seen_entries: set[str] = set()
    seen_locus_tags: dict[str, str] = {}
    loaded = [(item, load_organelle_records(item, config.base_dir)) for item in config.organelles]
    reserved_tags = {str(tag) for _, records in loaded for record in records
                     for feature in record.features for tag in feature.qualifiers.get("locus_tag", [])}
    with ann.open("w", encoding="utf-8") as annotation, fasta.open("w", encoding="utf-8") as sequence_file:
        annotation.write(header)
        for item, records in loaded:
            for source in (item.fasta, item.gff, item.genbank):
                if source:
                    path = expand_path(source, config.base_dir)
                    inputs.append({"organelle": item.name, "path": str(path), "sha256": _digest(path.read_bytes())})
            for record in records:
                identifier = str(record.id)
                entry = f"{item.name}_{identifier}"
                if not re.fullmatch(r"[A-Za-z0-9_.-]{1,32}", entry) or entry in seen_entries:
                    raise MSSPackError(f"Invalid/duplicate organelle MSS entry name: {entry}")
                seen_entries.add(entry)
                topology = item.topology.get(identifier, str(record.annotations.get("topology", "")))
                if topology not in {"linear", "circular"}:
                    raise MSSPackError(f"Specify topology for organelle sequence {item.name}/{identifier}")
                completeness = item.completeness.get(identifier, "unknown")
                if completeness == "unknown":
                    warnings.append(f"{entry}: assembly completeness is unconfirmed")
                sequence = str(record.seq)
                sequence_file.write(f">{entry}\n")
                write_wrapped_sequence(sequence_file, sequence)
                sequence_file.write("//\n")
                source = render_source_feature(
                    contig_name=entry, length=len(sequence), organism_name=config.sample.scientific_name,
                    strain=config.sample.strain, mol_type="genomic DNA", country=config.sample.geo_loc_name,
                    isolate=config.sample.isolate, collection_date=config.sample.collection_date,
                    sex=config.sample.sex, tissue_type=config.sample.tissue_type,
                    isolation_source=config.sample.isolation_source,
                )
                suffix = ", complete sequence" if completeness == "complete" else ", partial sequence" if completeness == "partial" else ""
                definition = f"@@[organism]@@ {item.organelle} DNA, @@[submitter_seqid]@@{suffix}"
                source = re.sub(r"(?<=ff_definition\t)[^\n]*", definition, source)
                source = source.replace("@@[submitter_seqid]@@", "@@[entry]@@")
                source = "".join(line + "\n" for line in source.splitlines()
                                 if line.split("\t")[3] != "submitter_seqid")
                if topology == "circular":
                    annotation.write(f"{entry}\tTOPOLOGY\t\tcircular\t\n")
                    source = source[len(entry):]
                annotation.write(source)
                annotation.write(f"\t\t\torganelle\t{item.organelle}\n")
                written = {fields[3]: fields[4] for line in source.splitlines()
                           if len(fields := line.split("\t")) == 5}
                written["organelle"] = item.organelle
                for feature in record.features:
                    if feature.type != "source":
                        continue
                    for qualifier, values in feature_qualifiers(feature).items():
                        if qualifier in {"ff_definition", "submitter_seqid"}:
                            continue
                        if qualifier in written:
                            if values != [written[qualifier]]:
                                raise MSSPackError(f"Conflicting organelle source metadata: {entry}/{qualifier}")
                            continue
                        for value in values:
                            annotation.write(f"\t\t\t{_text(qualifier)}\t{_text(value)}\n")
                counts: Counter[str] = Counter()
                annotation.write(_render_features(record, item_name=item.name, code=item.genetic_code,
                                                  prefix=config.sample.locus_tag, tags=tags,
                                                  owners=seen_locus_tags, reserved_tags=reserved_tags, audit=audit,
                                                  counts=counts, warnings=warnings))
                if not counts:
                    warnings.append(f"{entry}: source-only candidate; biological annotation is missing")
                records_summary.append({"entry": entry, "input_id": identifier, "organelle": item.organelle,
                                        "length": len(sequence), "topology": topology,
                                        "completeness": completeness, "sequence_sha256": _digest(sequence.encode()),
                                        "input_feature_counts": dict(Counter(feature.type for feature in record.features
                                                                              if feature.type != "source")),
                                        "feature_counts": dict(counts)})
    nuclear_check: dict[str, str] | None = None
    if nuclear_annotation is not None:
        nuclear_annotation = nuclear_annotation.expanduser().resolve()
        nuclear_tags = {fields[4] for line in nuclear_annotation.read_text(encoding="utf-8").splitlines()
                        if len(fields := line.split("\t")) == 5 and fields[3] == "locus_tag"}
        collisions = nuclear_tags & set(seen_locus_tags)
        if collisions:
            raise MSSPackError("Nuclear/organelle locus_tag collision: " + ", ".join(sorted(collisions)[:10]))
        nuclear_check = {"annotation": str(nuclear_annotation), "sha256": _digest(nuclear_annotation.read_bytes())}
    validation = ValidationArtifacts(ann, fasta, ensure_dir(root / "logs" / work.name), work)
    if validate:
        run_validation(options=replace(ValidationOptions.from_config(config),
                                       review_trans_splicing_overlaps=True), artifacts=validation)
    else:
        write_validation_not_run_summary(artifacts=validation)
    checked = load_validation_summary(validation.validation_summary)
    if checked is None or checked.status != "passed":
        warnings.append("Official DDBJ validation has not passed")
    if checked:
        for check in checked.checks:
            if check.warning_count:
                warnings.append(f"{check.label}: {check.warning_count} official warnings. {check.message or ''}")
    summary = {"schema_version": 1, "datatype": "MISC", "status": "completed",
               "msspack_version": __version__, "implementation_fingerprint": implementation_fingerprint(),
               "inputs": inputs, "records": records_summary, "warnings": warnings,
               "feature_audit_count": len(audit), "sequence_bases_modified": False,
               "review_required": bool(warnings),
               "nuclear_locus_tag_check": nuclear_check,
               "nuclear_analysis_executed": False, "submission_performed": False,
               "validation_status": checked.status if checked else "not_run"}
    write_text(work / "summary.json", json.dumps(summary, indent=2) + "\n")
    write_text(work / "feature-audit.json", json.dumps(audit, indent=2) + "\n")
    flows = write_organelle_gene_flows(config, summary, work, audit=audit)
    _write_report(work, records_summary, warnings, checked, flows)
    final = publish_submission(root, sorted(path for path in work.iterdir() if path.is_file()))
    return OrganelleArtifacts(final, final / ann.name, final / fasta.name,
                             final / "summary.json", final / "index.html")


@locked_output
def prepare_organelles(
    config_file: str | Path, *, validate: bool = True, nuclear_annotation: str | Path | None = None,
) -> OrganelleArtifacts:
    return _prepare(load_config(config_file), validate=validate,
                    nuclear_annotation=Path(nuclear_annotation) if nuclear_annotation is not None else None)
