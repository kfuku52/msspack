from __future__ import annotations

from datetime import datetime
from pathlib import Path

from .gff import iter_gff_records
from .step_logging import write_step_log, write_step_metrics
from .utils import atomic_text_writer, ensure_dir


def build_annotation_table(
    *,
    gff_path: Path,
    output_path: Path,
    locus_tag_prefix: str,
    log_path: Path,
    metrics_path: Path | None = None,
) -> None:
    started_at = datetime.now()
    gene_to_mrnas: dict[str, list[str]] = {}
    transcript_to_product: dict[str, str] = {}
    hypothetical_count = 0
    other_product_count = 0

    for record in iter_gff_records(gff_path):
        if record.type == "gene":
            gene_id = record.attributes.get("ID")
            if gene_id:
                gene_to_mrnas.setdefault(gene_id, [])
        elif record.type in ("mRNA", "transcript"):
            mrna_id = record.attributes.get("ID")
            parent_gene_ids = [item for item in record.attributes.get("Parent", "").split(",") if item]
            if not (parent_gene_ids and mrna_id):
                continue
            for parent_gene_id in parent_gene_ids:
                gene_to_mrnas.setdefault(parent_gene_id, []).append(mrna_id)
            product = record.attributes.get("product")
            if product:
                transcript_to_product[mrna_id] = product
        elif record.type == "CDS":
            product = record.attributes.get("product")
            if product:
                for transcript_id in record.attributes.get("Parent", "").split(","):
                    if transcript_id:
                        transcript_to_product.setdefault(transcript_id, product)

    ensure_dir(output_path.parent)
    with atomic_text_writer(output_path) as out_handle:
        out_handle.write("ID\tDescription\tLocus_tag\n")
        for gene_id, mrna_ids in gene_to_mrnas.items():
            for mrna_id in mrna_ids:
                product = transcript_to_product.get(mrna_id, "hypothetical protein")
                locus_tag = gene_id
                already_prefixed = locus_tag == locus_tag_prefix or locus_tag.startswith(
                    locus_tag_prefix + "_"
                )
                if locus_tag_prefix and not already_prefixed:
                    locus_tag = f"{locus_tag_prefix}_{locus_tag}"
                out_handle.write(f"{mrna_id}\t{product}\t{locus_tag}\n")
                if product == "hypothetical protein":
                    hypothetical_count += 1
                else:
                    other_product_count += 1

    row_count = hypothetical_count + other_product_count
    write_step_log(
        log_path=log_path,
        command=f"msspack internal gff2annotation-table --gff {gff_path} --output {output_path}",
        step="gff2annotation-table",
        started_at=started_at,
        count_unit="annotation rows",
        input_total=row_count,
        changed_total=row_count,
        output_total=row_count,
        details=[
            f"Hypothetical protein count: {hypothetical_count:,}",
            f"Other product count: {other_product_count:,}",
        ],
    )
    if metrics_path is not None:
        write_step_metrics(
            metrics_path=metrics_path,
            step="gff2annotation-table",
            count_unit="annotation rows",
            input_total=row_count,
            changed_total=row_count,
            output_total=row_count,
            details={
                "hypothetical_count": hypothetical_count,
                "other_product_count": other_product_count,
                "locus_tag_prefix": locus_tag_prefix,
            },
        )
