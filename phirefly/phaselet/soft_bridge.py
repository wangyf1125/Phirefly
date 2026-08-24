"""Soft phaselet-edge coupling from medium-risk SNP parity links."""

from __future__ import annotations

from collections import defaultdict

import numpy as np
from scipy import sparse


def build_soft_phaselet_coupling(
    soft_snp_edges: list[dict[str, object]],
    phaselet_ids: np.ndarray,
    tau_vec: np.ndarray,
    n_phaselets: int,
    soft_edge_scale: float = 1.0,
    soft_edge_cap: float = 10.0,
):
    coupling: dict[tuple[int, int], float] = defaultdict(float)
    support: dict[tuple[int, int], int] = defaultdict(int)
    raw_abs: dict[tuple[int, int], float] = defaultdict(float)
    for edge in soft_snp_edges:
        i = int(edge["snp_i"])
        j = int(edge["snp_j"])
        left = int(phaselet_ids[i])
        right = int(phaselet_ids[j])
        if left == right:
            continue
        if left > right:
            left, right = right, left
        coeff = float(edge["score"]) * int(tau_vec[i]) * int(tau_vec[j]) * float(soft_edge_scale)
        if soft_edge_cap > 0:
            coeff = float(np.clip(coeff, -float(soft_edge_cap), float(soft_edge_cap)))
        key = (left, right)
        coupling[key] += coeff
        support[key] += int(edge["count"])
        raw_abs[key] += abs(float(edge["score"]))

    rows = []
    cols = []
    data = []
    edge_rows = []
    for (left, right), coeff in sorted(coupling.items()):
        if coeff == 0.0:
            continue
        rows.extend([left, right])
        cols.extend([right, left])
        data.extend([coeff, coeff])
        edge_rows.append(
            {
                "left_phaselet": left,
                "right_phaselet": right,
                "soft_coupling": f"{coeff:.8f}",
                "soft_support_count": support[(left, right)],
                "soft_raw_abs_support": f"{raw_abs[(left, right)]:.8f}",
            }
        )
    matrix = sparse.coo_matrix(
        (
            np.asarray(data, dtype=np.float64),
            (np.asarray(rows, dtype=np.int64), np.asarray(cols, dtype=np.int64)),
        ),
        shape=(n_phaselets, n_phaselets),
    ).tocsr()
    matrix.sum_duplicates()
    stats = {
        "soft_phaselet_coupling_edges": len(edge_rows),
        "soft_phaselet_coupling_nnz": matrix.nnz,
        "soft_phaselet_coupling_abs_sum": f"{float(abs(matrix).sum()) / 2.0:.8f}",
    }
    return matrix, edge_rows, stats


__all__ = ["build_soft_phaselet_coupling"]
