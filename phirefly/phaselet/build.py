"""Phaselet construction from read-induced SNP-SNP parity evidence."""

from __future__ import annotations

import csv
import gzip
from collections import defaultdict
from pathlib import Path

import numpy as np

from .risk import PhaseEdge, classify_parity_edges, edge_row, soft_edge_rows


class DSU:
    def __init__(self, n: int):
        self.parent = list(range(n))
        self.size = [1] * n

    def find(self, x: int) -> int:
        root = x
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[x] != x:
            nxt = self.parent[x]
            self.parent[x] = root
            x = nxt
        return root

    def union(self, a: int, b: int) -> bool:
        ra = self.find(a)
        rb = self.find(b)
        if ra == rb:
            return False
        if self.size[ra] < self.size[rb]:
            ra, rb = rb, ra
        self.parent[rb] = ra
        self.size[ra] += self.size[rb]
        return True


def build_phaselet_components(n_snps: int, hard_edges: list[PhaseEdge]) -> tuple[DSU, int, list[float], list[tuple[PhaseEdge, str]]]:
    """Merge low-risk parity edges into conservative phaselet components."""

    dsu = DSU(n_snps)
    retained = 0
    retained_risks: list[float] = []
    hard_status: list[tuple[PhaseEdge, str]] = []
    for edge in sorted(hard_edges, key=lambda item: abs(item.score), reverse=True):
        if dsu.union(edge.snp_i, edge.snp_j):
            retained += 1
            retained_risks.append(edge.risk)
            status = "hard_retained"
        else:
            status = "hard_cycle"
        hard_status.append((edge, status))
    return dsu, retained, retained_risks, hard_status


def build_phaselet_table(dsu: DSU, snp_ids: list[str], positions: np.ndarray) -> tuple[np.ndarray, list[dict[str, object]]]:
    root_to_members: dict[int, list[int]] = defaultdict(list)
    for idx in range(len(snp_ids)):
        root_to_members[dsu.find(idx)].append(idx)
    ordered_roots = sorted(root_to_members, key=lambda root: root_to_members[root][0])
    block_ids = np.empty(len(snp_ids), dtype=np.int64)
    phaselets: list[dict[str, object]] = []
    for block_id, root in enumerate(ordered_roots):
        members = root_to_members[root]
        for idx in members:
            block_ids[idx] = block_id
        phaselets.append(
            {
                "phaselet_id": block_id,
                "chrom": snp_ids[members[0]].split(":", 1)[0],
                "start": int(positions[members[0]]),
                "end": int(positions[members[-1]]),
                "n_snps": len(members),
            }
        )
    return block_ids, phaselets


def build_phaselets(
    read_entries: list[list[tuple[int, int, float]]],
    snp_ids: list[str],
    tau_vec: np.ndarray,
    min_abs_support: float,
    min_count: int,
    min_confidence: float,
    max_distance_bp: int,
    bridge_mode: str = "confidence",
    risk_threshold: float = 0.55,
    risk_alpha_conf: float = 0.40,
    risk_beta_count: float = 0.15,
    risk_gamma_conflict: float = 0.15,
    risk_delta_distance: float = 0.10,
    risk_epsilon_domination: float = 0.20,
    risk_distance_norm_bp: int = 1_000_000,
    soft_bridge_mode: str = "none",
    soft_min_abs_support: float = 1.0,
    soft_min_count: int = 1,
    soft_min_confidence: float = 0.55,
    soft_risk_threshold: float = 0.75,
    parity_edges: dict[tuple[int, int], tuple[float, float, int, float]] | None = None,
):
    positions = np.fromiter((int(snp_id.split(":", 2)[1]) for snp_id in snp_ids), dtype=np.int64, count=len(snp_ids))
    if parity_edges is None:
        parity_edges = accumulate_adjacent_parity_edges(read_entries, positions, max_distance_bp=max_distance_bp)

    hard_edges, soft_edges, edge_rows, class_stats = classify_parity_edges(
        parity_edges,
        snp_ids=snp_ids,
        tau_vec=tau_vec,
        positions=positions,
        max_distance_bp=max_distance_bp,
        min_abs_support=min_abs_support,
        min_count=min_count,
        min_confidence=min_confidence,
        bridge_mode=bridge_mode,
        risk_threshold=risk_threshold,
        risk_alpha_conf=risk_alpha_conf,
        risk_beta_count=risk_beta_count,
        risk_gamma_conflict=risk_gamma_conflict,
        risk_delta_distance=risk_delta_distance,
        risk_epsilon_domination=risk_epsilon_domination,
        risk_distance_norm_bp=risk_distance_norm_bp,
        soft_bridge_mode=soft_bridge_mode,
        soft_min_abs_support=soft_min_abs_support,
        soft_min_count=soft_min_count,
        soft_min_confidence=soft_min_confidence,
        soft_risk_threshold=soft_risk_threshold,
    )
    dsu, retained, retained_risks, hard_status = build_phaselet_components(len(snp_ids), hard_edges)
    edge_rows.extend(edge_row(edge, snp_ids, positions, status) for edge, status in hard_status)
    block_ids, phaselets = build_phaselet_table(dsu, snp_ids, positions)
    soft_risks = [edge.risk for edge in soft_edges]
    stats = {
        "bridge_mode": bridge_mode,
        "risk_threshold": risk_threshold,
        "soft_bridge_mode": soft_bridge_mode,
        "soft_min_abs_support": soft_min_abs_support,
        "soft_min_count": soft_min_count,
        "soft_min_confidence": soft_min_confidence,
        "soft_risk_threshold": soft_risk_threshold,
        "parity_edges": len(parity_edges),
        "retained_phaselet_edges": retained,
        "soft_bridge_edges": len(soft_edges),
        **class_stats,
        "mean_retained_bridge_risk": f"{float(np.mean(retained_risks)):.8f}" if retained_risks else "0.00000000",
        "max_retained_bridge_risk": f"{float(np.max(retained_risks)):.8f}" if retained_risks else "0.00000000",
        "mean_soft_bridge_risk": f"{float(np.mean(soft_risks)):.8f}" if soft_risks else "0.00000000",
        "max_soft_bridge_risk": f"{float(np.max(soft_risks)):.8f}" if soft_risks else "0.00000000",
    }
    return block_ids, phaselets, stats, edge_rows, soft_edge_rows(soft_edges)


def accumulate_adjacent_parity_edges(
    read_entries: list[list[tuple[int, int, float]]],
    positions: np.ndarray,
    max_distance_bp: int,
) -> dict[tuple[int, int], tuple[float, float, int, float]]:
    edge_score: dict[tuple[int, int], float] = defaultdict(float)
    edge_abs: dict[tuple[int, int], float] = defaultdict(float)
    edge_count: dict[tuple[int, int], int] = defaultdict(int)
    edge_max_abs: dict[tuple[int, int], float] = defaultdict(float)

    for entries in read_entries:
        if len(entries) < 2:
            continue
        for left, right in zip(entries, entries[1:]):
            i, b_i, w_i = left
            j, b_j, w_j = right
            if i == j:
                continue
            if abs(int(positions[j]) - int(positions[i])) > max_distance_bp:
                continue
            if i > j:
                i, j = j, i
                b_i, b_j = b_j, b_i
                w_i, w_j = w_j, w_i
            gamma = min(float(w_i), float(w_j))
            key = (i, j)
            edge_score[key] += gamma * int(b_i) * int(b_j)
            edge_abs[key] += abs(gamma)
            edge_count[key] += 1
            edge_max_abs[key] = max(edge_max_abs[key], abs(gamma))

    return {
        key: (float(edge_score[key]), float(edge_abs[key]), int(edge_count[key]), float(edge_max_abs[key]))
        for key in edge_score
    }


def open_maybe_gzip(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return path.open()


def load_parity_edges(path: Path) -> dict[tuple[int, int], tuple[float, float, int, float]]:
    edges: dict[tuple[int, int], tuple[float, float, int, float]] = {}
    with open_maybe_gzip(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            i = int(row["snp_i"])
            j = int(row["snp_j"])
            if i == j:
                continue
            if i > j:
                i, j = j, i
            score = float(row["score"])
            support = float(row["support"])
            count = int(row["count"])
            max_abs_raw = row.get("max_abs_contribution") or row.get("max_abs")
            max_abs = float(max_abs_raw) if max_abs_raw not in {None, ""} else support / max(1, count)
            edges[(i, j)] = (score, support, count, max_abs)
    return edges


__all__ = [
    "DSU",
    "accumulate_adjacent_parity_edges",
    "build_phaselet_components",
    "build_phaselet_table",
    "build_phaselets",
    "load_parity_edges",
]
