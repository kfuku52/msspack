from __future__ import annotations

import re
from collections.abc import Iterable
from pathlib import Path

from .config import MSSPackConfig, StCommentConfig
from .config_validation import normalize_genome_coverage
from .utils import write_text


def _row(
    entry: str = "",
    feature: str = "",
    location: str = "",
    qualifier: str = "",
    value: str = "",
) -> str:
    return "\t".join([entry, feature, location, qualifier, value])


def _rows_for_repeated_qualifiers(
    feature: str,
    qualifier: str,
    values: Iterable[str],
    *,
    entry: str = "",
) -> list[str]:
    values = list(values)
    if not values:
        return []
    rows = [_row(entry, feature, "", qualifier, values[0])]
    for value in values[1:]:
        rows.append(_row("", "", "", qualifier, value))
    return rows


def _append_block(
    lines: list[str],
    *,
    feature: str,
    qualifier: str,
    values: Iterable[str],
) -> None:
    lines.extend(_rows_for_repeated_qualifiers(feature, qualifier, values))


def render_st_comment(config: StCommentConfig, *, entry: str = "") -> str:
    lines = [_row(entry, "ST_COMMENT", "", "tagset_id", config.tagset_id)]
    lines.extend(
        _row("", "", "", qualifier, value)
        for qualifier, value in (
            ("Assembly Method", config.assembly_method),
            ("Assembly Name", config.assembly_name),
            ("Genome Coverage", normalize_genome_coverage(config.genome_coverage)),
            ("Sequencing Technology", config.sequencing_technology),
        )
        if value
    )
    return "\n".join(lines) + "\n"


def render_header(config: MSSPackConfig, *, include_st_comment: bool = True) -> str:
    lines: list[str] = []

    common_entry = "COMMON"
    if config.submission.hold_date:
        lines.append(_row(common_entry, "DATE", "", "hold_date", config.submission.hold_date))
        common_entry = ""
    lines.append(_row(common_entry, "DATATYPE", "", "type", config.submission.datatype))
    lines.append(_row("", "DBLINK", "", "project", config.submission.bioproject))
    lines.append(_row("", "", "", "biosample", config.submission.biosample))
    for sra_id in config.submission.sra:
        lines.append(_row("", "", "", "sequence read archive", sra_id))

    _append_block(lines, feature="KEYWORD", qualifier="keyword", values=config.submission.keywords)
    _append_block(lines, feature="SUBMITTER", qualifier="ab_name", values=config.submitter.ab_name)

    lines.extend(
        [
            _row("", "", "", "contact", config.submitter.contact),
            _row("", "", "", "institute", config.submitter.institute),
            _row("", "", "", "department", config.submitter.department),
            _row("", "", "", "country", config.submitter.country),
            _row("", "", "", "state", config.submitter.state),
            _row("", "", "", "city", config.submitter.city),
            _row("", "", "", "street", config.submitter.street),
            _row("", "", "", "zip", config.submitter.zip),
            _row("", "", "", "phone", config.submitter.phone),
            _row("", "", "", "email", config.submitter.email),
        ]
    )

    lines.append(_row("", "REFERENCE", "", "title", config.reference.title))
    if config.reference.ab_name:
        lines.append(_row("", "", "", "ab_name", config.reference.ab_name[0]))
        for value in config.reference.ab_name[1:]:
            lines.append(_row("", "", "", "ab_name", value))
    lines.extend(
        [
            _row("", "", "", "year", str(config.reference.year)),
            _row("", "", "", "status", config.reference.status),
        ]
    )

    header = "\n".join(lines) + "\n"
    if include_st_comment:
        header += render_st_comment(config.st_comment)
    return header


def normalize_annotation_body(body: str) -> str:
    return re.sub(r"(^\t\t\t)country(\t)", r"\1geo_loc_name\2", body, flags=re.MULTILINE)


def render_final_annotation(header: str, body: str) -> str:
    return header + normalize_annotation_body(body)


def render_final_annotation_file(
    *,
    header: str,
    body_path: str | Path,
    output_path: str | Path,
) -> Path:
    body = Path(body_path).read_text(encoding="utf-8")
    return write_text(Path(output_path), render_final_annotation(header, body))
