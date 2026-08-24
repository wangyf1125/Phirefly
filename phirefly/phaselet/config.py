"""Argument parsing and input validation for phaselet QAIA runs."""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path


def required_inputs(region_key: str, args: argparse.Namespace) -> dict[str, Path | str]:
    required = [
        ("--observations", args.observations),
        ("--input-vcf-gz", args.input_vcf_gz),
        ("--tau-snp-phases", args.tau_snp_phases),
        ("--region", args.region),
    ]
    missing = [name for name, value in required if value is None]
    if missing:
        raise SystemExit("phaselet QAIA requires explicit input paths: " + ", ".join(missing))
    if (args.truth_vcf is None) != (args.input_vcf_plain is None):
        raise SystemExit("--truth-vcf and --input-vcf-plain must be provided together for benchmark metrics")
    return {
        "label": str(args.region_label or region_key),
        "region": str(args.region),
        "observations": Path(args.observations),
        "input_vcf_gz": Path(args.input_vcf_gz),
        "input_vcf_plain": Path(args.input_vcf_plain) if args.input_vcf_plain else None,
        "truth_bcf": Path(args.truth_vcf) if args.truth_vcf else None,
        "tau_snp_phases": Path(args.tau_snp_phases),
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("--task-id", type=int, default=None, help=argparse.SUPPRESS)
    parser.add_argument("--region-key", default="custom", help="Region label/key used only for output grouping.")
    parser.add_argument("--algorithm", choices=["CFC"], default="CFC")
    parser.add_argument("--phaselet-config-name", help="Phaselet config label for explicit region mode.")
    parser.add_argument("--min-abs-support", type=float, help="Minimum absolute parity support for phaselet edges.")
    parser.add_argument("--min-count", type=int, help="Minimum read count for phaselet edges.")
    parser.add_argument("--min-confidence", type=float, help="Minimum parity confidence for phaselet edges.")
    parser.add_argument("--max-distance-bp", type=int, help="Maximum adjacent SNP edge distance.")
    parser.add_argument("--outroot", type=Path, default=Path("phirefly_output"))
    parser.add_argument("--python", default=sys.executable)
    parser.add_argument("--bgzip", default=os.environ.get("BGZIP", "bgzip"))
    parser.add_argument("--tabix", default=os.environ.get("TABIX", "tabix"))
    parser.add_argument("--observations", type=Path, help="Read-SNP observations in binary NPZ or text TSV(.gz) format.")
    parser.add_argument("--input-vcf-gz", type=Path, help="Unphased input VCF.gz used for VCF export.")
    parser.add_argument("--input-vcf-plain", type=Path, help="Unphased input VCF/VCF.gz used only for benchmark metrics.")
    parser.add_argument("--truth-vcf", type=Path, help="Truth phased VCF/BCF used only for benchmark metrics.")
    parser.add_argument("--sample", default=None, help="Input/predicted VCF sample. Default: first sample in input VCF.")
    parser.add_argument("--truth-sample", default=None, help="Truth VCF sample for benchmark metrics. Default: first sample.")
    parser.add_argument("--tau-snp-phases", type=Path, help="MSF/local backbone SNP phase TSV.")
    parser.add_argument("--parity-edges", type=Path, help="Optional extraction-time adjacent parity edge TSV(.gz).")
    parser.add_argument("--region", help="Genomic interval, e.g. chr6:28510120-33480577.")
    parser.add_argument("--region-label", help="Human-readable region label for output tables.")
    parser.add_argument("--backend", default="gpu-float32")
    parser.add_argument("--batch-size", type=int, default=100)
    parser.add_argument("--n-iter", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=23)
    parser.add_argument("--seed-scale", type=float, default=0.001)
    parser.add_argument("--seed-noise", type=float, default=0.003)
    parser.add_argument("--hyperread-norm", choices=["sum", "mean", "sqrt", "cap", "sqrt_cap", "degree_norm_cap", "adaptive_dom_cap"], default="sqrt_cap")
    parser.add_argument("--hyperread-cap", type=float, default=10.0)
    parser.add_argument("--adaptive-dom-threshold", type=float, default=0.30)
    parser.add_argument("--bridge-mode", choices=["confidence", "risk"], default="risk")
    parser.add_argument("--risk-threshold", type=float, default=0.45)
    parser.add_argument("--risk-alpha-conf", type=float, default=0.40)
    parser.add_argument("--risk-beta-count", type=float, default=0.15)
    parser.add_argument("--risk-gamma-conflict", type=float, default=0.15)
    parser.add_argument("--risk-delta-distance", type=float, default=0.10)
    parser.add_argument("--risk-epsilon-domination", type=float, default=0.20)
    parser.add_argument("--risk-distance-norm-bp", type=int, default=1_000_000)
    parser.add_argument("--soft-bridge-mode", choices=["none", "phaselet"], default="phaselet")
    parser.add_argument("--soft-min-abs-support", type=float, default=1.0)
    parser.add_argument("--soft-min-count", type=int, default=1)
    parser.add_argument("--soft-min-confidence", type=float, default=0.55)
    parser.add_argument("--soft-risk-threshold", type=float, default=0.75)
    parser.add_argument("--soft-edge-scale", type=float, default=0.25)
    parser.add_argument("--soft-edge-cap", type=float, default=10.0)
    parser.add_argument("--validation-fraction", type=float, default=0.0)
    parser.add_argument("--validation-weight", type=float, default=1.0)
    parser.add_argument("--validation-salt", default="phirefly_phaselet_v3")
    parser.add_argument("--max-qaia-nodes", type=int, default=46_000)
    parser.add_argument("--weighted", action="store_true")
    parser.add_argument("--debug-output", action="store_true", help="Write detailed phaselet/hyperread debug TSVs.")
    parser.add_argument("--force", action="store_true")
    return parser


def parse_args() -> argparse.Namespace:
    args = build_parser().parse_args()
    if args.task_id is not None:
        raise SystemExit("The built-in benchmark task grid was removed from public Phirefly; pass explicit input paths instead.")
    missing = [
        name
        for name, value in [
            ("--phaselet-config-name", args.phaselet_config_name),
            ("--min-abs-support", args.min_abs_support),
            ("--min-count", args.min_count),
            ("--min-confidence", args.min_confidence),
            ("--max-distance-bp", args.max_distance_bp),
        ]
        if value is None
    ]
    if missing:
        raise SystemExit("phaselet QAIA requires " + ", ".join(missing))
    return args


__all__ = ["build_parser", "parse_args", "required_inputs"]
