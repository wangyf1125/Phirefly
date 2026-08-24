"""Command-line runner for phaselet/hyperread Phirefly-QAIA.

This module is deliberately input-path driven. The public package does not
include cluster-specific built-in benchmark paths; Slurm examples pass all
BAM/VCF-derived inputs explicitly.
"""

from __future__ import annotations

import time
from pathlib import Path

import numpy as np

from ..core.matrix import build_matrix, load_npz_graph, load_observations
from ..core.spins import sign_keep_zero
from ..hyperread.build import build_hyperread_matrix, build_read_entries, split_read_entries
from ..qaia import algorithm_dt, solve_qaia
from .build import build_phaselets, load_parity_edges
from .config import parse_args, required_inputs
from .output import (
    copy_solution_outputs,
    export_component_vcf,
    load_phases,
    read_key_value,
    read_method_metrics,
    run_benchmark_metrics,
    write_config_metrics,
    write_debug_tables,
    write_solution,
)
from .soft_bridge import build_soft_phaselet_coupling

PACKAGE_DIR = Path(__file__).resolve().parents[1]


def run_phaselet_qaia(args) -> None:
    t0 = time.monotonic()
    region_key = str(args.region_key or args.region_label or "custom")
    algorithm = args.algorithm
    cfg_name = str(args.phaselet_config_name)
    min_abs_support = float(args.min_abs_support)
    min_count = int(args.min_count)
    min_confidence = float(args.min_confidence)
    max_distance_bp = int(args.max_distance_bp)
    inputs = required_inputs(region_key, args)
    dt = float(algorithm_dt(algorithm))
    outdir = args.outroot / region_key / f"alg_{algorithm}" / f"phaselet_{cfg_name}"
    phirefly_dir = outdir / "phirefly"
    logs_dir = outdir / "logs"
    paper_dir = outdir / "paper_metrics"
    for directory in [phirefly_dir, logs_dir, paper_dir, args.outroot / "tables"]:
        directory.mkdir(parents=True, exist_ok=True)

    for key in ["observations", "input_vcf_gz", "tau_snp_phases"]:
        path = Path(inputs[key])
        if not path.exists():
            raise SystemExit(f"missing {key}: {path}")
    for key in ["input_vcf_plain", "truth_bcf"]:
        path = inputs[key]
        if path is not None and not Path(path).exists():
            raise SystemExit(f"missing {key}: {path}")

    prefix = phirefly_dir / "phaselet_hyperread_qaia"
    solve_time = logs_dir / "solve.time.txt"
    component_vcf = outdir / "phased.component_ps.vcf"
    component_vcf_gz = outdir / "phased.component_ps.vcf.gz"
    metrics_tsv = paper_dir / "phaselet_hyperread.metrics.tsv"

    if args.force or not prefix.with_suffix(".summary.tsv").exists():
        time_start = time.monotonic()
        observations_path = str(inputs["observations"])
        if observations_path.endswith(".npz"):
            read_ids, snp_ids, read_matrix, read_entries, observation_count = load_npz_graph(
                observations_path,
                weighted=args.weighted,
            )
        else:
            obs = load_observations(observations_path, weighted=args.weighted, weight_column="weight")
            read_ids, snp_ids, read_matrix = build_matrix(obs)
            read_entries = build_read_entries(obs, read_ids, snp_ids)
            observation_count = len(obs)
        train_entries, train_mask, val_mask = split_read_entries(
            read_ids,
            read_entries,
            validation_fraction=float(args.validation_fraction),
            validation_salt=str(args.validation_salt),
        )
        train_read_matrix = read_matrix[train_mask]
        validation_read_matrix = read_matrix[val_mask] if np.any(val_mask) else None
        if args.parity_edges is not None and np.any(val_mask):
            raise SystemExit("--parity-edges cannot be combined with --validation-fraction > 0 unless generated from the same train split")
        tau = load_phases(Path(inputs["tau_snp_phases"]))
        missing_tau = [snp_id for snp_id in snp_ids if snp_id not in tau]
        if missing_tau:
            raise SystemExit(f"missing tau phases for {len(missing_tau)} observed SNPs")
        tau_vec = np.fromiter((tau[snp_id] for snp_id in snp_ids), dtype=np.int8, count=len(snp_ids))

        precomputed_parity_edges = load_parity_edges(args.parity_edges) if args.parity_edges is not None else None
        phaselet_ids, phaselets, phaselet_stats, phaselet_edge_rows, soft_snp_edges = build_phaselets(
            train_entries,
            snp_ids=snp_ids,
            tau_vec=tau_vec,
            min_abs_support=min_abs_support,
            min_count=min_count,
            min_confidence=min_confidence,
            max_distance_bp=max_distance_bp,
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
            parity_edges=precomputed_parity_edges,
        )
        hyper_matrix, hyperread_rows, hyperread_edge_rows, hyper_stats = build_hyperread_matrix(
            train_entries,
            tau_vec=tau_vec,
            phaselet_ids=phaselet_ids,
            n_phaselets=len(phaselets),
            hyperread_norm=args.hyperread_norm,
            hyperread_cap=args.hyperread_cap,
            adaptive_dom_threshold=args.adaptive_dom_threshold,
        )
        soft_phaselet_coupling, soft_phaselet_edge_rows, soft_phaselet_stats = build_soft_phaselet_coupling(
            soft_snp_edges,
            phaselet_ids=phaselet_ids,
            tau_vec=tau_vec,
            n_phaselets=len(phaselets),
            soft_edge_scale=args.soft_edge_scale,
            soft_edge_cap=args.soft_edge_cap,
        )
        compressed_nodes = hyper_matrix.shape[0] + hyper_matrix.shape[1]
        if compressed_nodes > args.max_qaia_nodes:
            raise SystemExit(
                f"phaselet-hyperread graph has {compressed_nodes} nodes, above --max-qaia-nodes {args.max_qaia_nodes}; "
                "use stricter phaselet merging or component-wise execution"
            )

        qaia_result = solve_qaia(
            hyper_matrix,
            phaselet_coupling=soft_phaselet_coupling,
            algorithm=algorithm,
            backend=args.backend,
            batch_size=args.batch_size,
            n_iter=args.n_iter,
            dt=dt,
            seed=args.seed,
            seed_scale=args.seed_scale,
            seed_noise=args.seed_noise,
        )
        z_samples = qaia_result.phaselet_spins
        delta_samples = tau_vec[:, None] * z_samples[phaselet_ids, :]
        full_read_fields = np.asarray(read_matrix @ delta_samples, dtype=np.float64)
        profiled_scores = np.abs(full_read_fields).sum(axis=0)
        train_scores = np.abs(np.asarray(train_read_matrix @ delta_samples, dtype=np.float64)).sum(axis=0)
        if validation_read_matrix is not None:
            validation_scores = np.abs(np.asarray(validation_read_matrix @ delta_samples, dtype=np.float64)).sum(axis=0)
            selection_scores = train_scores + float(args.validation_weight) * validation_scores
            selection_objective = "train_plus_heldout_profiled_F"
        else:
            validation_scores = np.zeros_like(profiled_scores)
            selection_scores = profiled_scores
            selection_objective = "exact_original_read_profiled_F"
        best_idx = int(np.argmax(selection_scores))
        best_delta = delta_samples[:, best_idx].astype(np.int8, copy=False)
        full_sigma = sign_keep_zero(np.asarray(read_matrix @ best_delta, dtype=np.float64).ravel())
        total_weight = float(abs(read_matrix).sum())
        solve_elapsed = time.monotonic() - time_start

        with solve_time.open("w") as out:
            out.write(f"Elapsed (wall clock) time (h:mm:ss or m:ss): 0:{solve_elapsed:05.2f}\n")

        if args.debug_output:
            write_debug_tables(outdir, phaselets, phaselet_edge_rows, soft_phaselet_edge_rows, hyperread_rows, hyperread_edge_rows)

        phaselet_sizes = np.asarray([row["n_snps"] for row in phaselets], dtype=np.int64)
        summary = {
            "method": "phirefly_phaselet_hyperread_qaia",
            "algorithm": algorithm,
            "backend": args.backend,
            "batch_size": args.batch_size,
            "n_iter": args.n_iter,
            "dt": dt,
            "seed": args.seed,
            "seed_scale": args.seed_scale,
            "seed_noise": args.seed_noise,
            "hyperread_norm": args.hyperread_norm,
            "hyperread_cap": args.hyperread_cap,
            "adaptive_dom_threshold": args.adaptive_dom_threshold,
            "weighted": args.weighted,
            "debug_output": args.debug_output,
            "parity_edge_source": str(args.parity_edges) if args.parity_edges is not None else "from_observations",
            "phaselet_config": cfg_name,
            "min_abs_support": min_abs_support,
            "min_count": min_count,
            "min_confidence": min_confidence,
            "max_distance_bp": max_distance_bp,
            "bridge_mode": args.bridge_mode,
            "risk_threshold": args.risk_threshold,
            "risk_alpha_conf": args.risk_alpha_conf,
            "risk_beta_count": args.risk_beta_count,
            "risk_gamma_conflict": args.risk_gamma_conflict,
            "risk_delta_distance": args.risk_delta_distance,
            "risk_epsilon_domination": args.risk_epsilon_domination,
            "risk_distance_norm_bp": args.risk_distance_norm_bp,
            "soft_bridge_mode": args.soft_bridge_mode,
            "soft_min_abs_support": args.soft_min_abs_support,
            "soft_min_count": args.soft_min_count,
            "soft_min_confidence": args.soft_min_confidence,
            "soft_risk_threshold": args.soft_risk_threshold,
            "soft_edge_scale": args.soft_edge_scale,
            "soft_edge_cap": args.soft_edge_cap,
            "validation_fraction": args.validation_fraction,
            "validation_weight": args.validation_weight,
            "validation_salt": args.validation_salt,
            "train_reads": int(train_mask.sum()),
            "validation_reads": int(val_mask.sum()),
            "reads": len(read_ids),
            "snps": len(snp_ids),
            "observations": observation_count,
            "original_spins": len(read_ids) + len(snp_ids),
            "phaselets": len(phaselets),
            "largest_phaselet_snps": int(phaselet_sizes.max()) if phaselet_sizes.size else 0,
            "median_phaselet_snps": float(np.median(phaselet_sizes)) if phaselet_sizes.size else 0.0,
            "hyperreads": hyper_stats["hyperreads"],
            "single_phaselet_reads": hyper_stats["single_phaselet_reads"],
            "multi_phaselet_reads": hyper_stats["multi_phaselet_reads"],
            "compressed_nodes": compressed_nodes,
            "hyper_matrix_nnz": hyper_matrix.nnz,
            "coupling_nnz": qaia_result.coupling_nnz,
            "qaia_runtime_s": f"{qaia_result.runtime_s:.2f}",
            **phaselet_stats,
            **soft_phaselet_stats,
            "single_phaselet_constant": f"{hyper_stats['single_phaselet_constant']:.8f}",
            "max_phaselet_hyperread_domination": hyper_stats["max_phaselet_hyperread_domination"],
            "mean_phaselet_hyperread_domination": hyper_stats["mean_phaselet_hyperread_domination"],
            "compressed_best_score": f"{float(qaia_result.scores[best_idx]):.8f}",
            "selection_score_best": f"{float(selection_scores[best_idx]):.8f}",
            "profiled_F_train_best": f"{float(train_scores[best_idx]):.8f}",
            "profiled_F_validation_best": f"{float(validation_scores[best_idx]):.8f}",
            "profiled_F_best": f"{float(profiled_scores[best_idx]):.8f}",
            "profiled_disagreement": f"{(total_weight - float(profiled_scores[best_idx])) / 2.0:.8f}",
            "total_abs_observation_weight": f"{total_weight:.8f}",
            "best_candidate_index": best_idx,
            "candidate_count": z_samples.shape[1],
            "solver_runtime_s": f"{solve_elapsed:.2f}",
            "postprocess": "none",
            "selection_objective": selection_objective,
        }
        write_solution(prefix, read_ids, snp_ids, full_sigma, best_delta, summary)

    copy_solution_outputs(prefix, outdir)

    if args.force or not component_vcf_gz.exists():
        export_component_vcf(args, inputs, PACKAGE_DIR, logs_dir, outdir, component_vcf, component_vcf_gz)

    have_truth = inputs["truth_bcf"] is not None and inputs["input_vcf_plain"] is not None
    if have_truth and (args.force or not metrics_tsv.exists()):
        run_benchmark_metrics(args, inputs, PACKAGE_DIR, logs_dir, metrics_tsv, component_vcf_gz, solve_time)

    summary = read_key_value(outdir / "qaia_summary.tsv")
    component_metrics = read_method_metrics(metrics_tsv, "phirefly_phaselet_hyperread") if have_truth else {}
    write_config_metrics(
        outdir=outdir,
        inputs=inputs,
        args=args,
        region_key=region_key,
        algorithm=algorithm,
        cfg_name=cfg_name,
        dt=dt,
        have_truth=have_truth,
        summary=summary,
        component_metrics=component_metrics,
        component_vcf_gz=component_vcf_gz,
        metrics_tsv=metrics_tsv,
        solve_time=solve_time,
        t0=t0,
    )


def main() -> None:
    run_phaselet_qaia(parse_args())


if __name__ == "__main__":
    main()
