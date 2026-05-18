"""MSS conversion package.

Parts of the conversion logic in this package were originally adapted from
the MIT-licensed GFF2MSS project and later refactored into native msspack
modules. See THIRD_PARTY_NOTICES.md for attribution details.
"""

from .core import ContigSummary, ConversionOptions, ConversionSummary, convert_gff_to_mss
from .render import format_event_summary

__all__ = [
    "ContigSummary",
    "ConversionOptions",
    "ConversionSummary",
    "convert_gff_to_mss",
    "format_event_summary",
]
