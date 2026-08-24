"""One-command Phirefly pipeline.

Common public workflow:

    BAM + unphased VCF -> observations -> MSF -> phaselet QAIA -> VCF

Truth is optional and used only to append benchmark metrics.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from types import SimpleNamespace

from .core.matrix import load_observations
from .extract.cli import extract_observations
from .extract.sites import load_het_snps
from .msf import solve_msf_backbone, write_solution as write_msf_solution
from .phaselet.runner import run_phaselet_qaia


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run the Phirefly phaselet-hyperread QAIA pipeline.")
    parser.add_argument("--bam", required=True)
    parser.add_argument("--vcf", required=True, help="Input unphased heterozygous-site VCF.gz.")
    parser.add_argument("--sample", default=None, help="Input sample. Default: first sample in VCF.")
    parser.add_argument("--region", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--truth-vcf", default=None, help="Optional truth VCF/BCF for benchmark metrics only.")
    parser.add_argument("--truth-sample", default=None)
    parser.add_argument("--python", default=sys.executable)
    parser.add_argument("--bgzip", default="bgzip")
    parser.add_argument("--tabix", default="tabix")
    parser.add_argument("--backend", default="cpu-float32")
    parser.add_argument("--algorithm", choices=["CFC"], default="CFC")
    parser.add_argument("--batch-size", type=int, default=100)
    parser.add_argument("--n-iter", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=23)
    parser.add_argument("--seed-scale", type=float, default=0.001)
    parser.add_argument("--seed-noise", type=float, default=0.003)
    parser.add_argument("--weighted", action="store_true")
    parser.add_argument("--min-baseq", type=int, default=13)
    parser.add_argument("--min-mapq", type=int, default=20)
    parser.add_argument("--epsilon-floor", type=float, default=1e-4)
    parser.add_argument("--epsilon-ceiling", type=float, default=0.49)
    parser.add_argument("--mapq-cap", type=int, default=60)
    parser.add_argument("--default-baseq", type=int, default=None)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--chunk-size-bp", type=int, default=0)
    parser.add_argument("--phaselet-config-name", default="risk0p45_soft0p75_s0p25")
    parser.add_argument("--min-abs-support", type=float, default=2.0)
    parser.add_argument("--min-count", type=int, default=2)
    parser.add_argument("--min-confidence", type=float, default=0.65)
    parser.add_argument("--max-distance-bp", type=int, default=1_000_000)
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
    parser.add_argument("--hyperread-norm", choices=["sum", "mean", "sqrt", "cap", "sqrt_cap", "degree_norm_cap", "adaptive_dom_cap"], default="sqrt_cap")
    parser.add_argument("--hyperread-cap", type=float, default=10.0)
    parser.add_argument("--adaptive-dom-threshold", type=float, default=0.30)
    parser.add_argument("--validation-fraction", type=float, default=0.0)
    parser.add_argument("--validation-weight", type=float, default=1.0)
    parser.add_argument("--validation-salt", default="phirefly_phaselet_v3")
    parser.add_argument("--max-qaia-nodes", type=int, default=46_000)
    parser.add_argument("--debug-output", action="store_true")
    parser.add_argument("--force", action="store_true")
    return parser


def phaselet_args(args: argparse.Namespace, observations_npz: Path, parity_edges: Path, tau_phases: Path, phaselet_out: Path, sample: str | None) -> SimpleNamespace:
    return SimpleNamespace(
        task_id=None,
        region_key="phirefly",
        region_label=None,
        algorithm=args.algorithm,
        phaselet_config_name=args.phaselet_config_name,
        min_abs_support=args.min_abs_support,
        min_count=args.min_count,
        min_confidence=args.min_confidence,
        max_distance_bp=args.max_distance_bp,
        outroot=phaselet_out,
        python=args.python,
        bgzip=args.bgzip,
        tabix=args.tabix,
        observations=observations_npz,
        input_vcf_gz=Path(args.vcf),
        input_vcf_plain=Path(args.vcf) if args.truth_vcf else None,
        truth_vcf=Path(args.truth_vcf) if args.truth_vcf else None,
        sample=sample,
        truth_sample=args.truth_sample,
        tau_snp_phases=tau_phases,
        parity_edges=parity_edges,
        region=args.region,
        backend=args.backend,
        batch_size=args.batch_size,
        n_iter=args.n_iter,
        seed=args.seed,
        seed_scale=args.seed_scale,
        seed_noise=args.seed_noise,
        hyperread_norm=args.hyperread_norm,
        hyperread_cap=args.hyperread_cap,
        adaptive_dom_threshold=args.adaptive_dom_threshold,
        bridge_mode=args.bridge_mode,
        risk_threshold=args.risk_threshold,
        risk_alpha_conf=args.risk_alpha_conf,
        risk_beta_count=args.risk_beta_count,
        risk_gamma_conflict=args.risk_gamma_conflict,
        risk_delta_distance=args.risk_delta_distance,
        risk_epsilon_domination=args.risk_epsilon_domination,
        risk_distance_norm_bp=args.risk_distance_norm_bp,
        soft_bridge_mode=args.soft_bridge_mode,
        soft_min_abs_support=args.soft_min_abs_support,
        soft_min_count=args.soft_min_count,
        soft_min_confidence=args.soft_min_confidence,
        soft_risk_threshold=args.soft_risk_threshold,
        soft_edge_scale=args.soft_edge_scale,
        soft_edge_cap=args.soft_edge_cap,
        validation_fraction=args.validation_fraction,
        validation_weight=args.validation_weight,
        validation_salt=args.validation_salt,
        max_qaia_nodes=args.max_qaia_nodes,
        weighted=args.weighted,
        debug_output=args.debug_output,
        force=args.force,
    )


def run_pipeline(args: argparse.Namespace) -> None:
    out_dir = Path(args.out_dir)
    extract_prefix = out_dir / "extract" / "phirefly_reads"
    msf_prefix = out_dir / "msf" / "msf_backbone"
    phaselet_out = out_dir / "phaselet_qaia"
    extract_prefix.parent.mkdir(parents=True, exist_ok=True)
    msf_prefix.parent.mkdir(parents=True, exist_ok=True)
    phaselet_out.mkdir(parents=True, exist_ok=True)

    observations_npz = extract_prefix.with_suffix(".observations.npz")
    parity_edges = extract_prefix.with_suffix(".adjacent_parity_edges.tsv.gz")
    tau_phases = msf_prefix.with_suffix(".snp_phases.tsv")

    sample = args.sample
    if args.force or not observations_npz.exists() or not parity_edges.exists():
        sample, sites = load_het_snps(args.vcf, args.region, args.sample)
        if not sites:
            raise SystemExit(f"No heterozygous biallelic SNVs found for {sample} in {args.region}")
        extract_observations(
            bam_path=args.bam,
            sample=sample,
            sites=sites,
            out_prefix=extract_prefix,
            min_baseq=args.min_baseq,
            min_mapq=args.min_mapq,
            epsilon_floor=args.epsilon_floor,
            epsilon_ceiling=args.epsilon_ceiling,
            mapq_cap=args.mapq_cap,
            default_baseq=args.default_baseq,
            scan_mode="read",
            output_format="npz",
            threads=max(1, int(args.threads)),
            chunk_size_bp=max(0, int(args.chunk_size_bp)),
            write_parity_edges_flag=True,
            parity_max_distance_bp=max(0, int(args.max_distance_bp)),
            parity_weighted=args.weighted,
        )

    if args.force or not tau_phases.exists():
        obs = load_observations(str(observations_npz), args.weighted, "weight")
        read_ids, snp_ids, _matrix, sigma, delta, score, total_weight = solve_msf_backbone(obs, weighted=args.weighted)
        write_msf_solution(str(msf_prefix), read_ids, snp_ids, sigma, delta, score, total_weight)

    run_phaselet_qaia(phaselet_args(args, observations_npz, parity_edges, tau_phases, phaselet_out, sample))


def main() -> None:
    run_pipeline(build_parser().parse_args())


if __name__ == "__main__":
    main()


__all__ = ["build_parser", "phaselet_args", "run_pipeline"]
