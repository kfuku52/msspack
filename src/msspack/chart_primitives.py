from __future__ import annotations

from pathlib import Path

from .utils import atomic_binary_writer

TEXT_RGB = (17 / 255.0, 24 / 255.0, 39 / 255.0)
MUTED_RGB = (75 / 255.0, 85 / 255.0, 99 / 255.0)
GRID_RGB = (229 / 255.0, 231 / 255.0, 235 / 255.0)
CHART_FONT_SIZE_PT = 8
SVG_FONT_SIZE = f"{CHART_FONT_SIZE_PT}pt"

_HELVETICA_UPPERCASE_WIDTHS = (
    667,
    667,
    722,
    722,
    667,
    611,
    778,
    722,
    278,
    500,
    667,
    556,
    833,
    722,
    778,
    667,
    778,
    722,
    667,
    611,
    722,
    667,
    944,
    667,
    667,
    611,
)
_HELVETICA_LOWERCASE_WIDTHS = (
    556,
    556,
    500,
    556,
    556,
    278,
    556,
    556,
    222,
    222,
    500,
    222,
    833,
    556,
    556,
    556,
    556,
    333,
    500,
    278,
    556,
    500,
    722,
    500,
    500,
    500,
)
_HELVETICA_BOLD_UPPERCASE_WIDTHS = (
    722,
    722,
    722,
    722,
    667,
    611,
    778,
    722,
    278,
    556,
    722,
    611,
    833,
    722,
    778,
    667,
    778,
    722,
    667,
    611,
    722,
    722,
    944,
    722,
    722,
    611,
)
_HELVETICA_BOLD_LOWERCASE_WIDTHS = (
    556,
    611,
    556,
    611,
    556,
    333,
    611,
    611,
    278,
    278,
    556,
    278,
    889,
    611,
    611,
    611,
    611,
    389,
    556,
    333,
    611,
    556,
    778,
    556,
    556,
    500,
)
_HELVETICA_WIDTHS = {
    **dict(zip("ABCDEFGHIJKLMNOPQRSTUVWXYZ", _HELVETICA_UPPERCASE_WIDTHS, strict=True)),
    **dict(zip("abcdefghijklmnopqrstuvwxyz", _HELVETICA_LOWERCASE_WIDTHS, strict=True)),
    **{digit: 556 for digit in "0123456789"},
    " ": 278,
    ",": 278,
    "-": 333,
    ".": 278,
    "/": 278,
}
_HELVETICA_BOLD_WIDTHS = {
    **dict(zip("ABCDEFGHIJKLMNOPQRSTUVWXYZ", _HELVETICA_BOLD_UPPERCASE_WIDTHS, strict=True)),
    **dict(zip("abcdefghijklmnopqrstuvwxyz", _HELVETICA_BOLD_LOWERCASE_WIDTHS, strict=True)),
    **{digit: 556 for digit in "0123456789"},
    " ": 278,
    ",": 278,
    "-": 333,
    ".": 278,
    "/": 278,
}


def hex_to_rgb(hex_color: str) -> tuple[float, float, float]:
    return (
        int(hex_color[1:3], 16) / 255.0,
        int(hex_color[3:5], 16) / 255.0,
        int(hex_color[5:7], 16) / 255.0,
    )


def pdf_escape(text: str) -> str:
    return text.replace("\\", "\\\\").replace("(", "\\(").replace(")", "\\)")


def pdf_helvetica_text_width(text: str, *, size: float, bold: bool = False) -> float:
    widths = _HELVETICA_BOLD_WIDTHS if bold else _HELVETICA_WIDTHS
    return sum(widths.get(character, 556) for character in text) * size / 1000.0


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
    stream = "\n".join(commands).encode("latin-1")
    objects = [
        b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica >>",
        b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica-Bold >>",
        b"<< /Length %d >>\nstream\n" % len(stream) + stream + b"\nendstream",
        (
            f"<< /Type /Page /Parent 5 0 R /MediaBox [0 0 {width:.2f} {height:.2f}] "
            f"/Resources << /Font << /F1 1 0 R /F2 2 0 R >> >> /Contents 3 0 R >>"
        ).encode("latin-1"),
        b"<< /Type /Pages /Kids [4 0 R] /Count 1 >>",
        b"<< /Type /Catalog /Pages 5 0 R >>",
    ]
    pdf = bytearray(b"%PDF-1.4\n%\xe2\xe3\xcf\xd3\n")
    offsets = [0]
    for index, obj in enumerate(objects, start=1):
        offsets.append(len(pdf))
        pdf.extend(f"{index} 0 obj\n".encode("latin-1"))
        pdf.extend(obj)
        pdf.extend(b"\nendobj\n")
    xref_start = len(pdf)
    pdf.extend(f"xref\n0 {len(objects) + 1}\n".encode("latin-1"))
    pdf.extend(b"0000000000 65535 f \n")
    for offset in offsets[1:]:
        pdf.extend(f"{offset:010d} 00000 n \n".encode("latin-1"))
    pdf.extend(
        (
            f"trailer\n<< /Size {len(objects) + 1} /Root 6 0 R >>\n"
            f"startxref\n{xref_start}\n%%EOF\n"
        ).encode("latin-1")
    )
    with atomic_binary_writer(output_path) as handle:
        handle.write(pdf)
    return output_path
