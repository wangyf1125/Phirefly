"""Solution, diagnostic tables and run summaries for phaselet phasing."""

from __future__ import annotations

import csv
import shutil
import sys
from pathlib import Path

import numpy as np

from ..core.genomics import snp_sort_key
from ..qaia import algorithm_dt


RUN_SETTINGS = (
    "algorithm", "backend", "batch_size", "n_iter", "seed", "seed_scale",
    "seed_noise", "hyperread_norm", "hyperread_cap", "adaptive_dom_threshold",
    "weighted", "debug_output", "min_abs_support", "min_count", "min_confidence",
    "max_distance_bp", "bridge_mode", "risk_threshold", "risk_alpha_conf",
    "risk_beta_count", "risk_gamma_conflict", "risk_delta_distance",
    "risk_epsilon_domination", "risk_distance_norm_bp", "soft_bridge_mode",
    "soft_min_abs_support", "soft_min_count", "soft_min_confidence",
    "soft_risk_threshold", "soft_edge_scale", "soft_edge_cap",
    "validation_fraction", "validation_weight", "validation_salt", "trace_every",
)
CONFIG_FIELDS = (
    "min_abs_support", "min_count", "min_confidence", "max_distance_bp",
    "hyperread_norm", "hyperread_cap", "adaptive_dom_threshold", "debug_output",
    "parity_edge_source", "bridge_mode", "risk_threshold", "soft_bridge_mode",
    "soft_min_abs_support", "soft_min_count", "soft_min_confidence",
    "soft_risk_threshold", "soft_edge_scale", "soft_edge_cap",
    "validation_fraction", "validation_weight", "trace_every", "trace_output",
    "train_reads", "validation_reads",
)
GRAPH_FIELDS = (
    "reads", "snps", "observations", "original_spins", "phaselets",
    "largest_phaselet_snps", "median_phaselet_snps", "hyperreads",
    "single_phaselet_reads", "multi_phaselet_reads", "compressed_nodes",
    "hyper_matrix_nnz", "coupling_nnz", "parity_edges", "retained_phaselet_edges",
    "soft_bridge_edges", "weak_phaselet_edges", "conflicting_phaselet_edges",
    "mean_retained_bridge_risk", "max_retained_bridge_risk",
    "mean_soft_bridge_risk", "max_soft_bridge_risk",
    "soft_phaselet_coupling_edges", "soft_phaselet_coupling_nnz",
    "soft_phaselet_coupling_abs_sum", "max_phaselet_hyperread_domination",
    "mean_phaselet_hyperread_domination", "selection_score_best",
    "profiled_F_train_best", "profiled_F_validation_best", "profiled_F_best",
    "profiled_disagreement",
)
BENCHMARK_FIELDS = (
    "SE_percent", "HE_percent", "haplotype_N50_kb", "phased_block_N50_kb",
    "phasing_completeness_percent", "SNP_completeness_percent", "switches",
    "hamming", "hamming_denominator",
)
DEBUG_TABLES = {
    "phaselets": ("phaselet_id", "chrom", "start", "end", "n_snps"),
    "phaselet_edges": ("snp_i", "snp_j", "status"),
    "soft_phaselet_edges": ("left_phaselet", "right_phaselet", "soft_coupling"),
    "hyperreads": ("hyperread_id", "path", "signs", "support_reads"),
    "hyperread_edges": ("hyperread_id", "phaselet_id", "normalized_value"),
}


def summarize_solution(args, data, solution, phaselets, metadata) -> dict:
    """Score the selected solution and collect statistics without writing files."""
    raw_score = float(np.abs(data.matrix @ solution.delta).sum())
    train_score = float(np.abs(data.train_matrix @ solution.delta).sum())
    validation_score = (float(np.abs(data.validation_matrix @ solution.delta).sum())
                        if data.validation_matrix is not None else 0.0)
    total_weight = float(abs(data.matrix).sum())
    compressed_score = float(np.abs(solution.hyper @ solution.z).sum()
                             + .5 * solution.z @ (solution.soft @ solution.z))
    sizes = np.asarray([p["n_snps"] for p in phaselets], dtype=np.int64)
    hyper = solution.hyper_stats
    return {
        **metadata,
        **{key: getattr(args, key) for key in RUN_SETTINGS},
        "dt": float(algorithm_dt(args.algorithm)),
        "phaselet_config": str(args.phaselet_config_name),
        "parity_edge_source": str(args.parity_edges) if args.parity_edges is not None else "from_observations",
        "method": "phirefly_phaselet_hyperread_qaia",
        "pipeline_revision": "single_cfc_v1",
        "ps_policy": "retained_nonzero_couplings",
        "retained_graph_components": len(np.unique(solution.ps)),
        "qaia_calls": 1,
        "qaia_relaxed_runtime_s": solution.qaia.runtime_s,
        "qaia_runtime_s": f"{solution.qaia.runtime_s:.6f}",
        "train_reads": int(data.train_mask.sum()),
        "validation_reads": int(data.validation_mask.sum()),
        "reads": len(data.read_ids), "snps": len(data.snp_ids),
        "observations": data.observation_count,
        "original_spins": len(data.read_ids) + len(data.snp_ids),
        "phaselets": len(phaselets),
        "largest_phaselet_snps": int(sizes.max()) if sizes.size else 0,
        "median_phaselet_snps": float(np.median(sizes)) if sizes.size else 0.0,
        **{key: hyper[key] for key in ("hyperreads", "single_phaselet_reads", "multi_phaselet_reads",
                                      "max_phaselet_hyperread_domination", "mean_phaselet_hyperread_domination")},
        "compressed_nodes": sum(solution.hyper.shape),
        "hyper_matrix_nnz": solution.hyper.nnz,
        "coupling_nnz": solution.qaia.coupling_nnz,
        "retained_phaselet_edges": len(data.snp_ids) - len(phaselets),
        **solution.soft_stats,
        "single_phaselet_constant": f"{hyper['single_phaselet_constant']:.8f}",
        "compressed_best_score": f"{compressed_score:.8f}",
        "selection_score_best": f"{train_score + float(args.validation_weight) * validation_score:.8f}",
        "profiled_F_train_best": f"{train_score:.8f}",
        "profiled_F_validation_best": f"{validation_score:.8f}",
        "profiled_F_best": f"{raw_score:.8f}",
        "profiled_disagreement": f"{(total_weight - raw_score) / 2.0:.8f}",
        "total_abs_observation_weight": f"{total_weight:.8f}",
        "best_candidate_index": "per_component",
        "candidate_count": solution.qaia.phaselet_spins.shape[1],
        "solver_trace_tsv": solution.qaia.trace_path or "",
        "postprocess": "none", "selection_objective": solution.selection_mode,
    }


def write_solution(outdir: Path, read_ids, snp_ids, sigma, delta, summary, phase_sets) -> None:
    prefix = outdir / "phirefly" / "phaselet_hyperread_qaia"
    with prefix.with_suffix(".read_haplotags.tsv").open("w") as out:
        out.write("read_id\tsigma\n")
        for read_id, value in zip(read_ids, sigma):
            out.write(f"{read_id}\t{int(value)}\n")
    with prefix.with_suffix(".snp_phases.tsv").open("w") as out:
        out.write("snp_id\tdelta\tps\n")
        for i in sorted(range(len(snp_ids)), key=lambda i: snp_sort_key(snp_ids[i])):
            out.write(f"{snp_ids[i]}\t{int(delta[i])}\t{int(phase_sets[i])}\n")
    with prefix.with_suffix(".summary.tsv").open("w") as out:
        out.write("metric\tvalue\n")
        for key, value in summary.items():
            out.write(f"{key}\t{value}\n")
    # Preserve published output paths used by downstream Hi-C and collectors.
    for suffix, name in ((".summary.tsv", "qaia_summary.tsv"),
                         (".snp_phases.tsv", "snp_phases.tsv"),
                         (".read_haplotags.tsv", "read_haplotags.tsv")):
        shutil.copyfile(prefix.with_suffix(suffix), outdir / name)


def write_rows(path: Path, rows: list[dict], fallback_fields=()) -> None:
    with path.open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=list(rows[0]) if rows else fallback_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_debug_output(outdir, data, solution, ids, phaselets, edges, margins) -> None:
    from scipy import sparse

    sparse.save_npz(outdir / "hyperread_matrix.npz", solution.hyper)
    sparse.save_npz(outdir / "soft_phaselet_coupling.npz", solution.soft)
    np.savez_compressed(outdir / "phaselet_mapping.npz", phaselet_ids=ids,
                        tau=data.tau, snp_ids=np.asarray(data.snp_ids), ps=solution.ps)
    np.savez_compressed(outdir / "read_consistency.npz", snp_ids=np.asarray(data.snp_ids), margin=margins)
    rows = (phaselets, edges, solution.soft_rows, solution.hyper_rows, solution.hyper_edges)
    for (name, fields), values in zip(DEBUG_TABLES.items(), rows):
        write_rows(outdir / f"{name}.tsv", values, fields)


def write_config_metrics(outdir, inputs, args, summary, component_metrics, elapsed) -> None:
    have_truth = inputs["truth_bcf"] is not None
    config = {
        "region": inputs["label"],
        "region_key": str(args.region_key or args.region_label or "custom"),
        "region_interval": inputs["region"],
        "sample": args.sample or "", "truth_sample": args.truth_sample or "",
        "benchmark_metrics": int(have_truth),
        **{key: summary[key] for key in ("algorithm", "backend", "batch_size", "n_iter", "dt", "phaselet_config")},
        **{key: summary.get(key, "") for key in CONFIG_FIELDS},
        "postprocess": "none", "solver_runtime_s": summary["solver_runtime_s"],
        "config_wall_s": f"{elapsed:.2f}",
        **{key: summary.get(key, "") for key in GRAPH_FIELDS},
        **{key: component_metrics.get(key, "") for key in BENCHMARK_FIELDS},
        "output_dir": outdir,
        "phased_vcf": outdir / "phased.component_ps.vcf.gz",
        "snp_phases": outdir / "snp_phases.tsv",
        "read_haplotags": outdir / "read_haplotags.tsv",
        "solver_trace_tsv": summary.get("solver_trace_tsv", ""),
        **{f"{name}_tsv": outdir / f"{name}.tsv" if args.debug_output else "" for name in DEBUG_TABLES},
        "metrics_tsv": outdir / "paper_metrics/phaselet_hyperread.metrics.tsv" if have_truth else "",
        "solve_time_log": outdir / "logs/solve.time.txt",
    }
    path = outdir / "config_metrics.tsv"
    write_rows(path, [config])
    print(path)
    sys.stdout.write(path.read_text())
