"""Load evidence, build the final partition, solve once and export the result."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import time

import numpy as np
from scipy import sparse

from ..core.io import read_key_value, require_fresh_output
from ..core.runtime import check_runtime
from ..core.matrix import build_matrix, load_npz_graph, load_observations
from ..core.spins import sign_keep_zero
from ..hyperread.build import build_read_entries, split_read_entries
from ..vcf import export_component_vcf, load_phases, resolve_sample
from .build import build_phaselets, load_parity_edges
from .config import parse_args, required_inputs
from .consistency import relax_phaselets, solve_stage
from .output import summarize_solution, write_config_metrics, write_debug_output, write_solution


@dataclass
class ReadData:
    read_ids: list[str]
    snp_ids: list[str]
    matrix: sparse.csr_matrix
    train_entries: list
    train_mask: np.ndarray
    validation_mask: np.ndarray
    train_matrix: sparse.csr_matrix
    validation_matrix: sparse.csr_matrix | None
    tau: np.ndarray
    observation_count: int


def load_read_data(args, inputs) -> ReadData:
    path = str(inputs["observations"])
    if path.endswith(".npz"):
        read_ids, snp_ids, matrix, entries, count = load_npz_graph(path, args.weighted)
    else:
        obs = load_observations(path, args.weighted, "weight")
        read_ids, snp_ids, matrix = build_matrix(obs)
        entries, count = build_read_entries(obs, read_ids, snp_ids), len(obs)
    entries, train_mask, val_mask = split_read_entries(
        read_ids, entries, validation_fraction=float(args.validation_fraction),
        validation_salt=str(args.validation_salt))
    if args.parity_edges is not None and np.any(val_mask):
        raise SystemExit("--parity-edges cannot be combined with --validation-fraction > 0 unless generated from the same train split")
    phases = load_phases(inputs["tau_snp_phases"])
    missing = sum(snp_id not in phases for snp_id in snp_ids)
    if missing:
        raise SystemExit(f"missing tau phases for {missing} observed SNPs")
    tau = np.fromiter((phases[s] for s in snp_ids), dtype=np.int8, count=len(snp_ids))
    return ReadData(read_ids, snp_ids, matrix, entries, train_mask, val_mask,
                    matrix[train_mask], matrix[val_mask] if np.any(val_mask) else None, tau, count)


def build_initial_phaselets(data, args):
    parity = load_parity_edges(args.parity_edges, snp_ids=data.snp_ids) if args.parity_edges is not None else None
    return build_phaselets(
        data.train_entries, snp_ids=data.snp_ids, tau_vec=data.tau,
        min_abs_support=float(args.min_abs_support), min_count=int(args.min_count),
        min_confidence=float(args.min_confidence), max_distance_bp=int(args.max_distance_bp),
        bridge_mode=args.bridge_mode, risk_threshold=args.risk_threshold,
        risk_alpha_conf=args.risk_alpha_conf, risk_beta_count=args.risk_beta_count,
        risk_gamma_conflict=args.risk_gamma_conflict, risk_delta_distance=args.risk_delta_distance,
        risk_epsilon_domination=args.risk_epsilon_domination, risk_distance_norm_bp=args.risk_distance_norm_bp,
        soft_bridge_mode=args.soft_bridge_mode, soft_min_abs_support=args.soft_min_abs_support,
        soft_min_count=args.soft_min_count, soft_min_confidence=args.soft_min_confidence,
        soft_risk_threshold=args.soft_risk_threshold, parity_edges=parity)


def prepare_output(args, region_key) -> Path:
    outdir = Path(args.outroot) / region_key / f"alg_{args.algorithm}" / f"phaselet_{args.phaselet_config_name}"
    cached = outdir / "phirefly/phaselet_hyperread_qaia.summary.tsv"
    if not args.force and cached.exists() and read_key_value(cached).get("pipeline_revision") != "single_cfc_v1":
        raise ValueError("Existing output is not single-round CFC; use a fresh directory or --force")
    require_fresh_output(outdir, args.force)
    for name in ("phirefly", "logs", "paper_metrics"):
        (outdir / name).mkdir(parents=True, exist_ok=True)
    return outdir


def evaluate_run(args, inputs, outdir) -> dict:
    """Optional benchmark; metrics never participate in inference or selection."""
    if inputs["truth_bcf"] is None:
        return {}
    from ..metrics import compute_metrics_rows, write_metrics_rows

    path = outdir / "paper_metrics/phaselet_hyperread.metrics.tsv"
    rows = compute_metrics_rows(
        input_vcf=inputs["input_vcf_plain"], truth_vcf=inputs["truth_bcf"],
        region_text=inputs["region"], input_sample=args.sample or None,
        truth_sample=args.truth_sample or None, pred_sample=args.sample or None,
        methods=[("phirefly_phaselet_hyperread", outdir / "phased.component_ps.vcf.gz",
                  [outdir / "logs/solve.time.txt"])])
    write_metrics_rows(path, rows)
    (outdir / "logs/metrics.stdout.txt").write_text(str(path) + "\n")
    (outdir / "logs/metrics.stderr.txt").write_text("")
    return rows[0]


def run_phaselet_qaia(args) -> None:
    started = time.monotonic()
    region_key = str(args.region_key or args.region_label or "custom")
    inputs = required_inputs(region_key, args)
    tools = check_runtime(args.bgzip, args.tabix)
    args.bgzip, args.tabix = tools["bgzip"], tools["tabix"]
    args.sample = resolve_sample(inputs["input_vcf_gz"], args.sample)
    if inputs["truth_bcf"] is not None:
        args.truth_sample = resolve_sample(inputs["truth_bcf"], args.truth_sample)
    outdir = prepare_output(args, region_key)

    solve_started = time.monotonic()
    data = load_read_data(args, inputs)
    metadata = {"t_load_inputs": time.monotonic() - solve_started}
    stage_started = time.monotonic()
    _, initial_phaselets, initial_stats, edges, soft_edges = build_initial_phaselets(data, args)
    metadata.update(t_build_phaselets=time.monotonic() - stage_started,
                    initial_phaselets=len(initial_phaselets))
    metadata.update({f"initial_{key}": value for key, value in initial_stats.items()})

    stage_started = time.monotonic()
    ids, phaselets, margins = relax_phaselets(data.train_matrix, data.tau, data.snp_ids, edges)
    metadata["t_read_consistency"] = time.monotonic() - stage_started
    trace = (args.trace_output or outdir / "logs/solver_trace.tsv") if args.trace_every > 0 else None
    row_weights = (np.where(data.train_mask, 1.0, float(args.validation_weight))
                   if data.validation_matrix is not None else None)
    solution = solve_stage(
        data.train_entries, data.matrix, data.snp_ids, data.tau, ids, phaselets,
        soft_edges, args, row_weights=row_weights, trace_output=trace, stage="relaxed")
    metadata.update(solution.timings)

    released = 0
    for edge in edges:
        if edge["status"] in {"hard_retained", "hard_cycle"} and (
                margins[int(edge["snp_i"])] < 0 or margins[int(edge["snp_j"])] < 0):
            edge["status"] = "read_consistency_released"
            released += 1
    metadata.update(suspect_snps=int(np.count_nonzero(margins < 0)),
                    read_consistency_released_edges=released, trace_output=str(trace) if trace else "")
    summary = summarize_solution(args, data, solution, phaselets, metadata)
    sigma = sign_keep_zero(np.asarray(data.matrix @ solution.delta).ravel())
    elapsed = time.monotonic() - solve_started
    summary["solver_runtime_s"] = f"{elapsed:.2f}"
    (outdir / "logs/solve.time.txt").write_text(
        f"Elapsed (wall clock) time (h:mm:ss or m:ss): 0:{elapsed:05.2f}\n")
    if args.debug_output:
        write_debug_output(outdir, data, solution, ids, phaselets, edges, margins)
    write_solution(outdir, data.read_ids, data.snp_ids, sigma, solution.delta, summary, solution.ps)
    export_component_vcf(args, inputs, outdir)
    metrics = evaluate_run(args, inputs, outdir)
    write_config_metrics(outdir, inputs, args, summary, metrics, time.monotonic() - started)


def main() -> None:
    try:
        run_phaselet_qaia(parse_args())
    except ValueError as error:
        raise SystemExit(f"phirefly: {error}") from None


if __name__ == "__main__":
    main()
