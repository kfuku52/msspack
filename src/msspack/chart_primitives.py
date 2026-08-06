from __future__ import annotations

import re
from functools import cache
from io import BytesIO
from pathlib import Path
from typing import Any

import reportlab  # type: ignore[import-untyped]
from reportlab.pdfbase import pdfmetrics  # type: ignore[import-untyped]
from reportlab.pdfbase.ttfonts import TTFont  # type: ignore[import-untyped]
from reportlab.pdfgen.canvas import Canvas  # type: ignore[import-untyped]

from .utils import atomic_binary_writer

TEXT_RGB = (17 / 255.0, 24 / 255.0, 39 / 255.0)
MUTED_RGB = (75 / 255.0, 85 / 255.0, 99 / 255.0)
GRID_RGB = (229 / 255.0, 231 / 255.0, 235 / 255.0)
CHART_FONT_SIZE_PT = 8
SVG_FONT_SIZE = f"{CHART_FONT_SIZE_PT}pt"

_PDF_REGULAR_FONT = "MSSPackVera"
_PDF_BOLD_FONT = "MSSPackVeraBold"
_PDF_FONT_COMMAND_RE = re.compile(r"/(?P<font>F[12]) (?P<size>[0-9.]+) Tf")
_PDF_TEXT_COMMAND_RE = re.compile(r"\((?P<text>(?:\\.|[^\\)])*)\) Tj")


def hex_to_rgb(hex_color: str) -> tuple[float, float, float]:
    return (
        int(hex_color[1:3], 16) / 255.0,
        int(hex_color[3:5], 16) / 255.0,
        int(hex_color[5:7], 16) / 255.0,
    )


def pdf_escape(text: str) -> str:
    return text.replace("\\", "\\\\").replace("(", "\\(").replace(")", "\\)")


@cache
def _register_pdf_fonts() -> tuple[str, str]:
    font_dir = Path(str(reportlab.__file__)).resolve().parent / "fonts"
    registered = set(pdfmetrics.getRegisteredFontNames())
    if _PDF_REGULAR_FONT not in registered:
        pdfmetrics.registerFont(TTFont(_PDF_REGULAR_FONT, font_dir / "Vera.ttf"))
    if _PDF_BOLD_FONT not in registered:
        pdfmetrics.registerFont(TTFont(_PDF_BOLD_FONT, font_dir / "VeraBd.ttf"))
    return _PDF_REGULAR_FONT, _PDF_BOLD_FONT


def pdf_embedded_text_width(text: str, *, size: float, bold: bool = False) -> float:
    regular_font, bold_font = _register_pdf_fonts()
    return float(pdfmetrics.stringWidth(text, bold_font if bold else regular_font, size))


def pdf_helvetica_text_width(text: str, *, size: float, bold: bool = False) -> float:
    """Backward-compatible alias for the embedded PDF font width."""
    return pdf_embedded_text_width(text, size=size, bold=bold)


def pdf_top_to_bottom(page_height: float, y_top: float, box_height: float = 0.0) -> float:
    return page_height - y_top - box_height


def pdf_text_command(
    *,
    page_height: float,
    x: float,
    y_top: float,
    text: str,
    font: str,
    size: int,
    color: tuple[float, float, float],
) -> str:
    y = pdf_top_to_bottom(page_height, y_top)
    return (
        f"BT {color[0]:.3f} {color[1]:.3f} {color[2]:.3f} rg "
        f"/{font} {size} Tf 1 0 0 1 {x:.2f} {y:.2f} Tm ({pdf_escape(text)}) Tj ET"
    )


def write_single_page_pdf(
    *,
    width: float,
    height: float,
    commands: list[str],
    output_path: Path,
) -> Path:
    return write_multi_page_pdf(
        pages=[(width, height, commands)],
        output_path=output_path,
    )


def write_multi_page_pdf(
    *,
    pages: list[tuple[float, float, list[str]]],
    output_path: Path,
) -> Path:
    if not pages:
        raise ValueError("PDF output requires at least one page")
    regular_font_name, bold_font_name = _register_pdf_fonts()
    buffer = BytesIO()
    first_width, first_height, _ = pages[0]
    canvas = Canvas(
        buffer,
        pagesize=(first_width, first_height),
        pageCompression=0,
        invariant=1,
        pdfVersion=(1, 4),
        initialFontName=regular_font_name,
    )
    fonts: dict[str, Any] = {
        "F1": pdfmetrics.getFont(regular_font_name),
        "F2": pdfmetrics.getFont(bold_font_name),
    }
    for width, height, commands in pages:
        canvas.setPageSize((width, height))
        for command in commands:
            canvas.addLiteral(_command_with_embedded_font(command, canvas, fonts))
        canvas.showPage()
    canvas.save()
    with atomic_binary_writer(output_path) as handle:
        handle.write(buffer.getvalue())
    return output_path


def _pdf_unescape(text: str) -> str:
    characters: list[str] = []
    escaped = False
    for character in text:
        if escaped:
            characters.append(character)
            escaped = False
        elif character == "\\":
            escaped = True
        else:
            characters.append(character)
    if escaped:
        characters.append("\\")
    return "".join(characters)


def _command_with_embedded_font(
    command: str,
    canvas: Canvas,
    fonts: dict[str, Any],
) -> str:
    font_match = _PDF_FONT_COMMAND_RE.search(command)
    text_match = _PDF_TEXT_COMMAND_RE.search(command)
    if font_match is None or text_match is None:
        return command
    font_alias = font_match.group("font")
    font_size = font_match.group("size")
    font = fonts[font_alias]
    text = _pdf_unescape(text_match.group("text"))
    chunks = font.splitString(text or " ", canvas._doc)
    rendered_chunks: list[str] = []
    first_internal_name = font.getSubsetInternalName(chunks[0][0], canvas._doc)
    if not text:
        return command.replace(f"/{font_alias}", first_internal_name, 1)
    for index, (subset, encoded) in enumerate(chunks):
        if index > 0:
            internal_name = font.getSubsetInternalName(subset, canvas._doc)
            rendered_chunks.append(f"{internal_name} {font_size} Tf")
        rendered_chunks.append(f"({canvas._escape(encoded)}) Tj")
    command = (
        command[: text_match.start()]
        + " ".join(rendered_chunks)
        + command[text_match.end() :]
    )
    return command.replace(f"/{font_alias}", first_internal_name, 1)
