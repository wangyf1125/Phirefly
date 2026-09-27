"""Phirefly package.

This package contains read-SNP phaselet/hyperread compression and QAIA-based
optimization for haplotype phasing.
"""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("phirefly")
except PackageNotFoundError:
    __version__ = "0+uninstalled"

__all__ = [
    "cli",
    "core",
    "extract",
    "hic",
    "hyperread",
    "metrics",
    "msf",
    "pipeline",
    "phaselet",
    "qaia",
    "vcf",
]
