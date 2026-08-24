"""Genomic ID and interval helpers."""

from __future__ import annotations


def parse_region(region: str | None, *, zero_based_start: bool = False) -> tuple[str, int, int] | None:
    if not region:
        return None
    chrom, rest = region.split(":", 1)
    start_s, end_s = rest.replace(",", "").split("-", 1)
    start = int(start_s)
    if zero_based_start:
        start -= 1
    return chrom, start, int(end_s)


def parse_snp_id(snp_id: str) -> tuple[str, int, str, str]:
    chrom, pos_s, alleles = snp_id.split(":", 2)
    ref, alt = alleles.split(">", 1)
    return chrom, int(pos_s), ref.upper(), alt.upper()


def snp_sort_key(snp_id: str) -> tuple[str, int, str]:
    chrom, pos, _ref, _alt = parse_snp_id(snp_id)
    return chrom, pos, snp_id


def parse_chrom_alias(values: list[str]) -> dict[str, str]:
    aliases: dict[str, str] = {}
    for value in values:
        if "=" not in value:
            raise SystemExit(f"--chrom-alias must be phirefly_chrom=bam_chrom, got {value!r}")
        phirefly_chrom, bam_chrom = value.split("=", 1)
        if not phirefly_chrom or not bam_chrom:
            raise SystemExit(f"--chrom-alias must be phirefly_chrom=bam_chrom, got {value!r}")
        aliases[phirefly_chrom] = bam_chrom
    return aliases


def apply_region_alias(
    region: tuple[str, int, int] | None,
    aliases: dict[str, str],
) -> tuple[str, int, int] | None:
    if region is None:
        return None
    chrom, start, end = region
    return aliases.get(chrom, chrom), start, end


__all__ = ["apply_region_alias", "parse_chrom_alias", "parse_region", "parse_snp_id", "snp_sort_key"]
