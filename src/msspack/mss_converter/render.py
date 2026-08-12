from __future__ import annotations

from collections.abc import Iterable


def render_source_feature(
    *,
    contig_name: str,
    length: int,
    organism_name: str,
    strain: str,
    mol_type: str,
    country: str,
    isolate: str,
    collection_date: str,
    sex: str,
    tissue_type: str = "",
    isolation_source: str = "",
) -> str:
    out = f"{contig_name}\tsource\t1..{length}"
    if isolate:
        out += "\tff_definition\t@@[organism]@@ @@[isolate]@@ DNA, @@[submitter_seqid]@@\n"
    else:
        out += "\tff_definition\t@@[organism]@@ DNA, @@[submitter_seqid]@@\n"
    out += f"\t\t\tmol_type\t{mol_type}\n"
    out += f"\t\t\torganism\t{organism_name}\n"
    if strain:
        out += f"\t\t\tstrain\t{strain}\n"
    if isolate:
        out += f"\t\t\tisolate\t{isolate}\n"
    if tissue_type:
        out += f"\t\t\ttissue_type\t{tissue_type}\n"
    if isolation_source:
        out += f"\t\t\tisolation_source\t{isolation_source}\n"
    if sex:
        out += f"\t\t\tsex\t{sex}\n"
    if country:
        out += f"\t\t\tgeo_loc_name\t{country}\n"
    if collection_date:
        out += f"\t\t\tcollection_date\t{collection_date}\n"
    out += "\t\t\tsubmitter_seqid\t@@[entry]@@\n"
    return out


def render_cds_feature(
    *,
    joined_location: str,
    locus_tag_prefix: str,
    locus_tag_counter: int,
    mrna_id: str,
    product_name: str,
    custom_locus_tag: str | None,
    transl_table: str,
    codon_start: int,
    zfill: int = 9,
) -> str:
    if custom_locus_tag is None:
        out = f"\tCDS\t{joined_location}\tlocus_tag\t{locus_tag_prefix}{str(locus_tag_counter).zfill(zfill)}\n"
    else:
        out = f"\tCDS\t{joined_location}\tlocus_tag\t{custom_locus_tag}\n"
    out += f"\t\t\tnote\ttranscript_id:{mrna_id}\n"
    out += f"\t\t\tproduct\t{product_name}\n"
    out += f"\t\t\ttransl_table\t{transl_table}\n"
    out += f"\t\t\tcodon_start\t{codon_start}\n"
    return out


def render_gap_feature(
    *,
    gap_start: int,
    gap_end: int,
    linkage_evidence: str,
    min_assembly_gap_size: int,
    gap_estimated_length: str,
) -> str:
    gap_size = gap_end - gap_start + 1
    if gap_size < min_assembly_gap_size:
        return ""
    if gap_start == gap_end:
        out = f"\tassembly_gap\t{gap_start}"
    else:
        out = f"\tassembly_gap\t{gap_start}..{gap_end}"
    out += f"\testimated_length\t{gap_estimated_length}\n"
    out += "\t\t\tgap_type\twithin scaffold\n"
    out += f"\t\t\tlinkage_evidence\t{linkage_evidence}\n"
    return out


def render_rrna_feature(name: str, position: str) -> str:
    mapping = {
        "18S": "\trRNA\t{position}\tproduct\t18S rRNA\n",
        "5.8S": "\trRNA\t{position}\tproduct\t5.8S rRNA\n",
        "28S": "\trRNA\t{position}\tproduct\t28S rRNA\n",
        "ITS1": "\tmisc_RNA\t{position}\tnote\tinternal transcribed spacer 1\n",
        "ITS2": "\tmisc_RNA\t{position}\tnote\tinternal transcribed spacer 2\n",
    }
    template = mapping.get(name)
    if template is not None:
        return template.format(position=position)
    product = name.strip() or "ribosomal RNA"
    if product.lower() != "ribosomal rna" and not product.lower().endswith("rrna"):
        product += " rRNA"
    return f"\trRNA\t{position}\tproduct\t{product}\n"


def render_trna_feature(
    *,
    position: str,
    locus_tag_prefix: str,
    locus_tag_counter: int,
    product: str,
    anticodon: str,
    note: str,
    zfill: int = 9,
) -> str:
    out = f"\ttRNA\t{position}\tproduct\t{product}\n"
    out += f"\t\t\tlocus_tag\t{locus_tag_prefix}{str(locus_tag_counter).zfill(zfill)}\n"
    if anticodon:
        out += f"\t\t\tanticodon\t{anticodon}\n"
    if note:
        out += f"\t\t\tnote\t{note}\n"
    return out


def render_generic_feature(
    *,
    feature_key: str,
    position: str,
    qualifiers: Iterable[tuple[str, str]],
) -> str:
    qualifier_rows = [(key, value) for key, value in qualifiers if key]
    out = f"\t{feature_key}\t{position}"
    if not qualifier_rows:
        return out + "\n"
    first_key, first_value = qualifier_rows[0]
    out += f"\t{first_key}\t{first_value}\n"
    for key, value in qualifier_rows[1:]:
        out += f"\t\t\t{key}\t{value}\n"
    return out


def format_event_summary(event_counts: dict[str, int], prefix: str) -> str:
    labels = [
        ("start_codon_missing", "start codon missing"),
        ("stop_codon_missing", "stop codon missing"),
        ("small_introns", "small introns"),
        ("gap_artificial_location", "gap -> artificial_location"),
        ("gap_misc_feature", "gap -> misc_feature"),
        ("non_cds_features", "non-CDS features rendered"),
        ("unknown_features", "unknown GFF types -> misc_feature"),
    ]
    parts = [
        f"{label}={event_counts[key]}"
        for key, label in labels
        if event_counts.get(key, 0)
    ]
    if not parts:
        return f"{prefix}: no warnings"
    return f"{prefix}: " + ", ".join(parts)


def append_repeated_qualifier(lines: list[str], qualifier: str, values: Iterable[str]) -> None:
    for value in values:
        lines.append(f"\t\t\t{qualifier}\t{value}")
