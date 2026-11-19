from __future__ import annotations

import re

from Bio.Data.IUPACData import protein_letters_1to3

_TRNA_GENE = re.compile(r"(?:trn|tRNA)(?:f(?=M))?([A-Z])(?=$|[-_(])", re.IGNORECASE)
_TRNA_PRODUCT = re.compile(r"tRNA[- ]([A-Z]{3})(?=$|\W)", re.IGNORECASE)
_AMINO_ACIDS = frozenset(protein_letters_1to3.values())


def trna_amino_acids(gene: str, product: str) -> tuple[str, str] | None:
    """Read the amino acids named in both qualifiers when their syntax is recognized."""
    gene_match = _TRNA_GENE.match(gene)
    product_match = _TRNA_PRODUCT.match(product)
    if not (gene_match and product_match):
        return None
    gene_amino_acid = protein_letters_1to3.get(gene_match[1].upper())
    product_amino_acid = product_match[1].title()
    if gene_amino_acid and product_amino_acid in _AMINO_ACIDS:
        return gene_amino_acid, product_amino_acid
    return None
