#!/usr/bin/env python3
"""Compute qHap paper-style phasing metrics from VCF outputs.

The output columns match the qHap paper tables: switch error rate, Hamming
error rate, haplotype N50, phasing completeness, SNP completeness, and running
time. This is intended for matched BAM/VCF/truth/region comparisons.
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import pysam

from .core.genomics import parse_region


@dataclass(frozen=True)
class PhaseCall:
    chrom: str
    pos: int
    ref: str
    alt: str
    phase: int
    ps: str

    @property
    def key(self) -> tuple[str, int, str, str]:
        return (self.chrom, self.pos, self.ref, self.alt)


def sample_name(vcf: pysam.VariantFile, requested: str | None) -> str:
    samples = list(vcf.header.samples)
    if not samples:
        raise ValueError("VCF has no samples")
    if requested:
        if requested not in samples:
            raise ValueError(f"Sample {requested!r} not found. Available: {samples}")
        return requested
    return samples[0]


def records(vcf: pysam.VariantFile, region: tuple[str, int, int] | None) -> Iterable[pysam.VariantRecord]:
    if region is None:
        yield from vcf
        return
    chrom, start, end = region
    try:
        yield from vcf.fetch(chrom, start - 1, end)
    except (ValueError, OSError):
        for rec in vcf:
            if rec.chrom == chrom and start <= rec.pos <= end:
                yield rec


def is_biallelic_snp(rec: pysam.VariantRecord) -> bool:
    return (
        rec.alts is not None
        and len(rec.alts) == 1
        and len(rec.ref) == 1
        and len(rec.alts[0]) == 1
    )


def het_phase_for_record(rec: pysam.VariantRecord, sample: str, require_phased: bool) -> tuple[int, str] | None:
    if not is_biallelic_snp(rec):
        return None
    call = rec.samples[sample]
    gt = call.get("GT")
    if gt not in [(0, 1), (1, 0)]:
        return None
    if require_phased and not call.phased:
        return None
    phase = 0 if gt == (0, 1) else 1
    ps = call.get("PS")
    if ps is None:
        ps = rec.chrom if require_phased else "unphased"
    return phase, str(ps)


def count_het_snps(path: Path, sample: str | None, region: tuple[str, int, int] | None) -> int:
    with pysam.VariantFile(str(path)) as vcf:
        sample = sample_name(vcf, sample)
        total = 0
        for rec in records(vcf, region):
            if not is_biallelic_snp(rec):
                continue
            gt = rec.samples[sample].get("GT")
            if gt in [(0, 1), (1, 0)]:
                total += 1
        return total


def load_phase_calls(
    path: Path,
    sample: str | None,
    region: tuple[str, int, int] | None,
    require_phased: bool,
) -> list[PhaseCall]:
    calls: list[PhaseCall] = []
    with pysam.VariantFile(str(path)) as vcf:
        sample = sample_name(vcf, sample)
        for rec in records(vcf, region):
            parsed = het_phase_for_record(rec, sample, require_phased=require_phased)
            if parsed is None:
                continue
            phase, ps = parsed
            calls.append(PhaseCall(rec.chrom, rec.pos, rec.ref, rec.alts[0], phase, ps))
    return calls


def haplotype_n50_kb(calls: list[PhaseCall]) -> float:
    by_block: dict[tuple[str, str], list[int]] = defaultdict(list)
    for call in calls:
        by_block[(call.chrom, call.ps)].append(call.pos)
    spans = []
    for positions in by_block.values():
        if len(positions) < 2:
            continue
        spans.append(max(positions) - min(positions) + 1)
    if not spans:
        return 0.0
    half = sum(spans) / 2
    total = 0
    for span in sorted(spans, reverse=True):
        total += span
        if total >= half:
            return span / 1000
    return 0.0


def phase_error_metrics(truth_calls: list[PhaseCall], pred_calls: list[PhaseCall]) -> dict[str, float | int]:
    truth_by_key = {call.key: call for call in truth_calls}
    blocks: dict[tuple[str, str], list[tuple[PhaseCall, PhaseCall]]] = defaultdict(list)
    for pred in pred_calls:
        truth = truth_by_key.get(pred.key)
        if truth is None:
            continue
        blocks[(pred.chrom, pred.ps)].append((pred, truth))

    switches = 0
    assessed_pairs = 0
    hamming = 0
    hamming_denominator = 0
    covered_variants = 0

    for block_calls in blocks.values():
        block_calls.sort(key=lambda item: item[0].pos)
        if len(block_calls) < 2:
            continue
        diffs = [pred.phase ^ truth.phase for pred, truth in block_calls]
        switches += sum(1 for a, b in zip(diffs, diffs[1:]) if a != b)
        assessed_pairs += len(diffs) - 1
        hamming += min(sum(diffs), len(diffs) - sum(diffs))
        hamming_denominator += len(diffs)
        covered_variants += len(diffs)

    return {
        "covered_variants": covered_variants,
        "assessed_pairs": assessed_pairs,
        "switches": switches,
        "switch_error_rate": switches / assessed_pairs if assessed_pairs else 0.0,
        "hamming": hamming,
        "hamming_denominator": hamming_denominator,
        "hamming_error_rate": hamming / hamming_denominator if hamming_denominator else 0.0,
    }


def parse_elapsed_seconds(path: Path) -> float:
    if not path.exists():
        return 0.0
    total = 0.0
    pattern = re.compile(r"Elapsed \(wall clock\) time .*: (?:(\d+):)?(\d+):(\d+(?:\.\d+)?)$")
    for line in path.read_text(errors="replace").splitlines():
        match = pattern.search(line.strip())
        if not match:
            continue
        hours = int(match.group(1) or 0)
        minutes = int(match.group(2))
        seconds = float(match.group(3))
        total += hours * 3600 + minutes * 60 + seconds
    return total


METRIC_FIELDNAMES = [
    "method",
    "SE_percent",
    "HE_percent",
    "haplotype_N50_kb",
    "phased_block_N50_kb",
    "phasing_completeness_percent",
    "SNP_completeness_percent",
    "running_time_s",
    "switches",
    "assessed_pairs",
    "hamming",
    "hamming_denominator",
    "phased_snps",
    "input_het_snps",
    "truth_het_snps",
]


def compute_metrics_rows(
    input_vcf: Path,
    truth_vcf: Path,
    region_text: str | None,
    input_sample: str | None,
    truth_sample: str | None,
    pred_sample: str | None,
    methods: list[tuple[str, Path, list[Path]]],
) -> list[dict[str, float | int | str]]:
    region = parse_region(region_text)
    input_het_count = count_het_snps(input_vcf, input_sample, region)
    truth_het_count = count_het_snps(truth_vcf, truth_sample, region)
    truth_calls = load_phase_calls(truth_vcf, truth_sample, region, require_phased=True)

    rows = []
    for name, pred_vcf, time_files in methods:
        pred_calls = load_phase_calls(pred_vcf, pred_sample, region, require_phased=True)
        metrics = phase_error_metrics(truth_calls, pred_calls)
        pred_truth_keys = {call.key for call in pred_calls}
        truth_keys = {call.key for call in truth_calls}
        retained_truth_snps = len(pred_truth_keys & truth_keys)
        runtime = sum(parse_elapsed_seconds(path) for path in time_files)
        n50_kb = haplotype_n50_kb(pred_calls)
        rows.append(
            {
                "method": name,
                "SE_percent": 100 * float(metrics["switch_error_rate"]),
                "HE_percent": 100 * float(metrics["hamming_error_rate"]),
                "haplotype_N50_kb": n50_kb,
                "phased_block_N50_kb": n50_kb,
                "phasing_completeness_percent": 100 * len(pred_calls) / input_het_count if input_het_count else 0.0,
                "SNP_completeness_percent": 100 * retained_truth_snps / truth_het_count if truth_het_count else 0.0,
                "running_time_s": runtime,
                "switches": metrics["switches"],
                "assessed_pairs": metrics["assessed_pairs"],
                "hamming": metrics["hamming"],
                "hamming_denominator": metrics["hamming_denominator"],
                "phased_snps": len(pred_calls),
                "input_het_snps": input_het_count,
                "truth_het_snps": truth_het_count,
            }
        )
    return rows


def write_metrics_rows(out_tsv: Path, rows: list[dict[str, float | int | str]]) -> None:
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with out_tsv.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=METRIC_FIELDNAMES, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

def method_arg(value: str) -> tuple[str, Path, list[Path]]:
    parts = value.split(":", 2)
    if len(parts) < 2:
        raise argparse.ArgumentTypeError("--method must be NAME:VCF[:TIME[,TIME...]]")
    name = parts[0]
    vcf = Path(parts[1])
    time_files = [Path(item) for item in parts[2].split(",")] if len(parts) == 3 and parts[2] else []
    return name, vcf, time_files


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-vcf", required=True)
    parser.add_argument("--truth-vcf", required=True)
    parser.add_argument("--region", default=None)
    parser.add_argument("--input-sample", default=None)
    parser.add_argument("--truth-sample", default=None)
    parser.add_argument("--pred-sample", default=None, help="Predicted VCF sample. Default: first sample in each predicted VCF.")
    parser.add_argument("--out-tsv", required=True)
    parser.add_argument("--method", action="append", type=method_arg, required=True)
    args = parser.parse_args()

    out_path = Path(args.out_tsv)
    rows = compute_metrics_rows(
        input_vcf=Path(args.input_vcf),
        truth_vcf=Path(args.truth_vcf),
        region_text=args.region,
        input_sample=args.input_sample,
        truth_sample=args.truth_sample,
        pred_sample=args.pred_sample,
        methods=args.method,
    )
    write_metrics_rows(out_path, rows)
    print(out_path)


if __name__ == "__main__":
    main()


__all__ = [
    "METRIC_FIELDNAMES",
    "PhaseCall",
    "compute_metrics_rows",
    "write_metrics_rows",
    "method_arg",
    "main",
]
