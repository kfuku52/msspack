from __future__ import annotations

from pathlib import Path

from .utils import ensure_dir

TEXT_RGB = (17 / 255.0, 24 / 255.0, 39 / 255.0)
MUTED_RGB = (75 / 255.0, 85 / 255.0, 99 / 255.0)
GRID_RGB = (229 / 255.0, 231 / 255.0, 235 / 255.0)


def hex_to_rgb(hex_color: str) -> tuple[float, float, float]:
    return (
        int(hex_color[1:3], 16) / 255.0,
        int(hex_color[3:5], 16) / 255.0,
        int(hex_color[5:7], 16) / 255.0,
    )


def pdf_escape(text: str) -> str:
    return text.replace("\\", "\\\\").replace("(", "\\(").replace(")", "\\)")


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
            f"<< /Type /Page /Parent 5 0 R /MediaBox [0 0 {int(width)} {int(height)}] "
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
    ensure_dir(output_path.parent)
    output_path.write_bytes(pdf)
    return output_path
