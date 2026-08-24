"""VCF site loading for read-SNP extraction."""

from __future__ import annotations

from dataclasses import dataclass

import pysam

from ..core.genomics import parse_region


@dataclass(frozen=True)
class SnpSite:
    chrom: str
    pos0: int
    snp_id: str
    ref: str
    alt: str


def is_het_gt(record: pysam.VariantRecord, sample: str) -> bool:
    if sample not in record.samples:
        return False
    gt = record.samples[sample].get("GT")
    if gt is None or len(gt) < 2:
        return False
    if None in gt:
        return False
    return len(set(gt[:2])) == 2


def choose_sample(vcf: pysam.VariantFile, requested: str | None) -> str:
    samples = list(vcf.header.samples)
    if not samples:
        raise SystemExit(f"VCF has no samples: {vcf.filename}")
    if requested:
        if requested not in samples:
            raise SystemExit(f"sample {requested!r} not found in {vcf.filename}. Available: {samples}")
        return requested
    return samples[0]


def load_het_snps(vcf_path: str, region: str, sample: str | None) -> tuple[str, list[SnpSite]]:
    parsed_region = parse_region(region, zero_based_start=True)
    if parsed_region is None:
        raise SystemExit("--region is required for SNP loading")
    chrom, start0, end0 = parsed_region
    vcf = pysam.VariantFile(vcf_path)
    sample = choose_sample(vcf, sample)
    sites: list[SnpSite] = []
    for rec in vcf.fetch(chrom, start0, end0):
        if not is_het_gt(rec, sample):
            continue
        if len(rec.ref) != 1 or len(rec.alts or []) != 1 or len(rec.alts[0]) != 1:
            continue
        sites.append(
            SnpSite(
                chrom=chrom,
                pos0=rec.pos - 1,
                snp_id=f"{chrom}:{rec.pos}:{rec.ref}>{rec.alts[0]}",
                ref=rec.ref.upper(),
                alt=rec.alts[0].upper(),
            )
        )
    return sample, sites


__all__ = ["SnpSite", "choose_sample", "is_het_gt", "load_het_snps"]
