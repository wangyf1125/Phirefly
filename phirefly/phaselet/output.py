"""Output, VCF export, and summary helpers for phaselet QAIA runs."""

from __future__ import annotations

import argparse
import csv
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

from ..metrics import compute_metrics_rows, write_metrics_rows
from ..vcf import export_phased_vcf


def snp_sort_key(snp_id: str) -> tuple[str, int, str]:
    chrom, pos_s, rest = snp_id.split(":", 2)
    return chrom, int(pos_s), rest


def load_phases(path: Path) -> dict[str, int]:
    out: dict[str, int] = {}
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            out[row["snp_id"]] = int(row["delta"])
    return out


def run_command(cmd: list[str], stdout_path: Path | None = None, stderr_path: Path | None = None) -> None:
    stdout = stdout_path.open("w") if stdout_path else None
    stderr = stderr_path.open("w") if stderr_path else None
    try:
        subprocess.run(cmd, check=True, stdout=stdout, stderr=stderr)
    finally:
        if stdout:
            stdout.close()
        if stderr:
            stderr.close()


def bgzip_and_index(bgzip: str, tabix: str, plain_vcf: Path, gz_vcf: Path) -> None:
    with gz_vcf.open("wb") as gz_out:
        subprocess.run([bgzip, "-f", "-c", str(plain_vcf)], check=True, stdout=gz_out)
    run_command([tabix, "-f", "-p", "vcf", str(gz_vcf)])


def parse_time_seconds(path: Path) -> float:
    if not path.exists():
        return 0.0
    pattern = re.compile(r"Elapsed \(wall clock\) time .*: (?:(\d+):)?(\d+):(\d+(?:\.\d+)?)$")
    for line in path.read_text(errors="replace").splitlines():
        match = pattern.search(line.strip())
        if match:
            return int(match.group(1) or 0) * 3600 + int(match.group(2)) * 60 + float(match.group(3))
    return 0.0


def read_key_value(path: Path) -> dict[str, str]:
    if not path.exists():
        return {}
    out: dict[str, str] = {}
    with path.open() as fh:
        reader = csv.reader(fh, delimiter="\t")
        next(reader, None)
        for row in reader:
            if len(row) >= 2:
                out[row[0]] = row[1]
    return out


def read_method_metrics(path: Path, method: str) -> dict[str, str]:
    if not path.exists():
        return {}
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row.get("method") == method:
                return row
    return {}


def write_solution(prefix: Path, read_ids: list[str], snp_ids: list[str], sigma, delta, summary) -> None:
    with prefix.with_suffix(".read_haplotags.tsv").open("w") as out:
        out.write("read_id\tsigma\n")
        for read_id, value in zip(read_ids, sigma):
            out.write(f"{read_id}\t{int(value)}\n")
    with prefix.with_suffix(".snp_phases.tsv").open("w") as out:
        out.write("snp_id\tdelta\n")
        for snp_id, value in sorted(zip(snp_ids, delta), key=lambda item: snp_sort_key(item[0])):
            out.write(f"{snp_id}\t{int(value)}\n")
    with prefix.with_suffix(".summary.tsv").open("w") as out:
        out.write("metric\tvalue\n")
        for key, value in summary.items():
            out.write(f"{key}\t{value}\n")


def copy_solution_outputs(prefix: Path, outdir: Path) -> None:
    shutil.copyfile(prefix.with_suffix(".summary.tsv"), outdir / "qaia_summary.tsv")
    shutil.copyfile(prefix.with_suffix(".snp_phases.tsv"), outdir / "snp_phases.tsv")
    shutil.copyfile(prefix.with_suffix(".read_haplotags.tsv"), outdir / "read_haplotags.tsv")


def export_component_vcf(
    args,
    inputs: dict[str, Path | str],
    package_dir: Path,
    logs_dir: Path,
    outdir: Path,
    component_vcf: Path,
    component_vcf_gz: Path,
) -> None:
    written, skipped, output_vcf = export_phased_vcf(
        input_vcf=Path(inputs["input_vcf_gz"]),
        pred_snp_phases=outdir / "snp_phases.tsv",
        region=str(inputs["region"]),
        output_vcf=component_vcf,
        sample=str(args.sample) if args.sample else None,
        ps_start=1,
        observations=Path(inputs["observations"]),
    )
    with (logs_dir / "vcf.stdout.txt").open("w") as out:
        out.write(f"written\t{written}\n")
        out.write(f"skipped\t{skipped}\n")
        out.write(f"output_vcf\t{output_vcf}\n")
    (logs_dir / "vcf.stderr.txt").write_text("")
    bgzip_and_index(args.bgzip, args.tabix, component_vcf, component_vcf_gz)


def run_benchmark_metrics(
    args,
    inputs: dict[str, Path | str],
    package_dir: Path,
    logs_dir: Path,
    metrics_tsv: Path,
    component_vcf_gz: Path,
    solve_time: Path,
) -> None:
    rows = compute_metrics_rows(
        input_vcf=Path(inputs["input_vcf_plain"]),
        truth_vcf=Path(inputs["truth_bcf"]),
        region_text=str(inputs["region"]),
        input_sample=str(args.sample) if args.sample else None,
        truth_sample=str(args.truth_sample) if args.truth_sample else None,
        pred_sample=str(args.sample) if args.sample else None,
        methods=[("phirefly_phaselet_hyperread", component_vcf_gz, [solve_time])],
    )
    write_metrics_rows(metrics_tsv, rows)
    with (logs_dir / "metrics.stdout.txt").open("w") as out:
        out.write(str(metrics_tsv) + "\n")
    (logs_dir / "metrics.stderr.txt").write_text("")

def write_rows(path: Path, rows: list[dict[str, object]], fallback_fields: list[str]) -> None:
    fieldnames = list(rows[0]) if rows else fallback_fields
    with path.open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_debug_tables(
    outdir: Path,
    phaselets: list[dict[str, object]],
    phaselet_edge_rows: list[dict[str, object]],
    soft_phaselet_edge_rows: list[dict[str, object]],
    hyperread_rows: list[dict[str, object]],
    hyperread_edge_rows: list[dict[str, object]],
) -> None:
    write_rows(outdir / "phaselets.tsv", phaselets, ["phaselet_id", "chrom", "start", "end", "n_snps"])
    write_rows(outdir / "phaselet_edges.tsv", phaselet_edge_rows, ["snp_i", "snp_j", "status"])
    write_rows(outdir / "soft_phaselet_edges.tsv", soft_phaselet_edge_rows, ["left_phaselet", "right_phaselet", "soft_coupling"])
    write_rows(outdir / "hyperreads.tsv", hyperread_rows, ["hyperread_id", "path", "signs", "support_reads"])
    write_rows(outdir / "hyperread_edges.tsv", hyperread_edge_rows, ["hyperread_id", "phaselet_id", "normalized_value"])


def write_config_metrics(
    outdir: Path,
    inputs: dict[str, Path | str],
    args: argparse.Namespace,
    region_key: str,
    algorithm: str,
    cfg_name: str,
    dt: float,
    have_truth: bool,
    summary: dict[str, str],
    component_metrics: dict[str, str],
    component_vcf_gz: Path,
    metrics_tsv: Path,
    solve_time: Path,
    t0: float,
) -> None:
    config_metrics = {
        "region": inputs["label"],
        "region_key": region_key,
        "region_interval": inputs["region"],
        "sample": args.sample or "",
        "truth_sample": args.truth_sample or "",
        "benchmark_metrics": int(have_truth),
        "algorithm": algorithm,
        "backend": args.backend,
        "batch_size": args.batch_size,
        "n_iter": args.n_iter,
        "dt": dt,
        "phaselet_config": cfg_name,
        "min_abs_support": summary.get("min_abs_support", ""),
        "min_count": summary.get("min_count", ""),
        "min_confidence": summary.get("min_confidence", ""),
        "max_distance_bp": summary.get("max_distance_bp", ""),
        "hyperread_norm": summary.get("hyperread_norm", ""),
        "hyperread_cap": summary.get("hyperread_cap", ""),
        "adaptive_dom_threshold": summary.get("adaptive_dom_threshold", ""),
        "debug_output": summary.get("debug_output", ""),
        "parity_edge_source": summary.get("parity_edge_source", ""),
        "bridge_mode": summary.get("bridge_mode", ""),
        "risk_threshold": summary.get("risk_threshold", ""),
        "soft_bridge_mode": summary.get("soft_bridge_mode", ""),
        "soft_min_abs_support": summary.get("soft_min_abs_support", ""),
        "soft_min_count": summary.get("soft_min_count", ""),
        "soft_min_confidence": summary.get("soft_min_confidence", ""),
        "soft_risk_threshold": summary.get("soft_risk_threshold", ""),
        "soft_edge_scale": summary.get("soft_edge_scale", ""),
        "soft_edge_cap": summary.get("soft_edge_cap", ""),
        "validation_fraction": summary.get("validation_fraction", ""),
        "validation_weight": summary.get("validation_weight", ""),
        "train_reads": summary.get("train_reads", ""),
        "validation_reads": summary.get("validation_reads", ""),
        "postprocess": "none",
        "solver_runtime_s": f"{parse_time_seconds(solve_time):.2f}",
        "config_wall_s": f"{time.monotonic() - t0:.2f}",
        "reads": summary.get("reads", ""),
        "snps": summary.get("snps", ""),
        "observations": summary.get("observations", ""),
        "original_spins": summary.get("original_spins", ""),
        "phaselets": summary.get("phaselets", ""),
        "largest_phaselet_snps": summary.get("largest_phaselet_snps", ""),
        "median_phaselet_snps": summary.get("median_phaselet_snps", ""),
        "hyperreads": summary.get("hyperreads", ""),
        "single_phaselet_reads": summary.get("single_phaselet_reads", ""),
        "multi_phaselet_reads": summary.get("multi_phaselet_reads", ""),
        "compressed_nodes": summary.get("compressed_nodes", ""),
        "hyper_matrix_nnz": summary.get("hyper_matrix_nnz", ""),
        "coupling_nnz": summary.get("coupling_nnz", ""),
        "parity_edges": summary.get("parity_edges", ""),
        "retained_phaselet_edges": summary.get("retained_phaselet_edges", ""),
        "soft_bridge_edges": summary.get("soft_bridge_edges", ""),
        "weak_phaselet_edges": summary.get("weak_phaselet_edges", ""),
        "conflicting_phaselet_edges": summary.get("conflicting_phaselet_edges", ""),
        "mean_retained_bridge_risk": summary.get("mean_retained_bridge_risk", ""),
        "max_retained_bridge_risk": summary.get("max_retained_bridge_risk", ""),
        "mean_soft_bridge_risk": summary.get("mean_soft_bridge_risk", ""),
        "max_soft_bridge_risk": summary.get("max_soft_bridge_risk", ""),
        "soft_phaselet_coupling_edges": summary.get("soft_phaselet_coupling_edges", ""),
        "soft_phaselet_coupling_nnz": summary.get("soft_phaselet_coupling_nnz", ""),
        "soft_phaselet_coupling_abs_sum": summary.get("soft_phaselet_coupling_abs_sum", ""),
        "max_phaselet_hyperread_domination": summary.get("max_phaselet_hyperread_domination", ""),
        "mean_phaselet_hyperread_domination": summary.get("mean_phaselet_hyperread_domination", ""),
        "selection_score_best": summary.get("selection_score_best", ""),
        "profiled_F_train_best": summary.get("profiled_F_train_best", ""),
        "profiled_F_validation_best": summary.get("profiled_F_validation_best", ""),
        "profiled_F_best": summary.get("profiled_F_best", ""),
        "profiled_disagreement": summary.get("profiled_disagreement", ""),
        "SE_percent": component_metrics.get("SE_percent", ""),
        "HE_percent": component_metrics.get("HE_percent", ""),
        "haplotype_N50_kb": component_metrics.get("haplotype_N50_kb", ""),
        "phased_block_N50_kb": component_metrics.get("phased_block_N50_kb", ""),
        "phasing_completeness_percent": component_metrics.get("phasing_completeness_percent", ""),
        "SNP_completeness_percent": component_metrics.get("SNP_completeness_percent", ""),
        "switches": component_metrics.get("switches", ""),
        "hamming": component_metrics.get("hamming", ""),
        "hamming_denominator": component_metrics.get("hamming_denominator", ""),
        "output_dir": outdir,
        "phased_vcf": component_vcf_gz,
        "snp_phases": outdir / "snp_phases.tsv",
        "read_haplotags": outdir / "read_haplotags.tsv",
        "phaselets_tsv": outdir / "phaselets.tsv" if args.debug_output else "",
        "phaselet_edges_tsv": outdir / "phaselet_edges.tsv" if args.debug_output else "",
        "soft_phaselet_edges_tsv": outdir / "soft_phaselet_edges.tsv" if args.debug_output else "",
        "hyperreads_tsv": outdir / "hyperreads.tsv" if args.debug_output else "",
        "hyperread_edges_tsv": outdir / "hyperread_edges.tsv" if args.debug_output else "",
        "metrics_tsv": metrics_tsv if have_truth else "",
        "solve_time_log": solve_time,
    }
    with (outdir / "config_metrics.tsv").open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=list(config_metrics), delimiter="\t")
        writer.writeheader()
        writer.writerow(config_metrics)
    print(outdir / "config_metrics.tsv")
    with (outdir / "config_metrics.tsv").open() as handle:
        sys.stdout.write(handle.read())


__all__ = [
    "bgzip_and_index",
    "copy_solution_outputs",
    "export_component_vcf",
    "load_phases",
    "parse_time_seconds",
    "read_key_value",
    "read_method_metrics",
    "run_benchmark_metrics",
    "run_command",
    "snp_sort_key",
    "write_config_metrics",
    "write_debug_tables",
    "write_solution",
]
