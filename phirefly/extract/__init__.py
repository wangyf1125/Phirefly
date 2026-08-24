"""BAM/VCF read-SNP observation extraction."""

from .parity import build_adjacent_parity_edges, write_parity_edges
from .sites import SnpSite

__all__ = ["SnpSite", "build_adjacent_parity_edges", "write_parity_edges"]
