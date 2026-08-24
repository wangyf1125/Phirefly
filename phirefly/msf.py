#!/usr/bin/env python3
"""Maximum-spanning-forest SNP backbone for Phirefly."""

from __future__ import annotations

import argparse
from collections import defaultdict

import numpy as np
from scipy import sparse

from .core.matrix import build_matrix, load_observations
from .core.spins import sign_keep_zero

Observation = tuple[str, str, int, float]
ParityEdges = dict[tuple[int, int], float]


def group_observations_by_read(obs: list[Observation], snp_ids: list[str]) -> dict[str, list[tuple[int, int, float]]]:
    """Group read-SNP observations as sorted SNP-index entries per read."""

    snp_index = {snp_id: i for i, snp_id in enumerate(snp_ids)}
    read_entries: dict[str, list[tuple[int, int, float]]] = defaultdict(list)
    for read_id, snp_id, allele, weight in obs:
        read_entries[read_id].append((snp_index[snp_id], int(allele), float(weight)))
    for entries in read_entries.values():
        entries.sort(key=lambda item: item[0])
    return read_entries


def build_snp_parity_edges(read_entries: dict[str, list[tuple[int, int, float]]]) -> ParityEdges:
    """Build adjacent SNP-SNP parity edges induced by reads."""

    edge_scores: ParityEdges = defaultdict(float)
    for entries in read_entries.values():
        if len(entries) < 2:
            continue
        for (i, b_i, w_i), (j, b_j, w_j) in zip(entries, entries[1:]):
            if i == j:
                continue
            if i > j:
                i, j = j, i
                b_i, b_j = b_j, b_i
                w_i, w_j = w_j, w_i
            edge_scores[(i, j)] += min(float(w_i), float(w_j)) * int(b_i) * int(b_j)
    return dict(edge_scores)


def build_msf_backbone(edge_scores: ParityEdges, n_snps: int) -> np.ndarray | None:
    """Solve the signed maximum spanning forest and return SNP signs tau."""

    if not edge_scores:
        return None

    parent = list(range(n_snps))
    rank = [0] * n_snps
    parity = [0] * n_snps

    def find_no_compress(x: int) -> tuple[int, int]:
        acc = 0
        while parent[x] != x:
            acc ^= parity[x]
            x = parent[x]
        return x, acc

    def union(a: int, b: int, relation: int) -> None:
        ra, pa = find_no_compress(a)
        rb, pb = find_no_compress(b)
        if ra == rb:
            return
        if rank[ra] < rank[rb]:
            ra, rb = rb, ra
            pa, pb = pb, pa
        parent[rb] = ra
        parity[rb] = pa ^ pb ^ relation
        if rank[ra] == rank[rb]:
            rank[ra] += 1

    for (i, j), score in sorted(edge_scores.items(), key=lambda item: abs(item[1]), reverse=True):
        union(i, j, 0 if score >= 0 else 1)

    delta = np.ones(n_snps, dtype=np.int8)
    for i in range(n_snps):
        _root, p = find_no_compress(i)
        delta[i] = 1 if p == 0 else -1
    return delta


def mst_delta_from_observations(obs: list[Observation], snp_ids: list[str]) -> np.ndarray | None:
    """Compatibility wrapper: observations -> parity graph -> MSF tau."""

    read_entries = group_observations_by_read(obs, snp_ids)
    return build_msf_backbone(build_snp_parity_edges(read_entries), len(snp_ids))


def infer_read_labels(matrix: sparse.spmatrix, delta: np.ndarray) -> np.ndarray:
    """Given SNP phases, infer read haplotype signs by best response."""

    return sign_keep_zero(np.asarray(matrix @ delta, dtype=np.float64).ravel())


def score_read_snp_objective(matrix: sparse.spmatrix, sigma: np.ndarray, delta: np.ndarray) -> float:
    """Evaluate the original read-SNP Phirefly objective."""

    return float(np.asarray((matrix @ delta) * sigma, dtype=np.float64).sum())


def solve_msf_backbone(obs: list[Observation], weighted: bool = False):
    """Build the read-SNP matrix and MSF backbone from loaded observations."""

    read_ids, snp_ids, matrix = build_matrix(obs)
    delta = mst_delta_from_observations(obs, snp_ids)
    if delta is None:
        delta = np.ones(len(snp_ids), dtype=np.int8)
    sigma = infer_read_labels(matrix, delta)
    score = score_read_snp_objective(matrix, sigma, delta)
    total_weight = float(abs(matrix).sum())
    return read_ids, snp_ids, matrix, sigma, delta, score, total_weight


def write_solution(
    prefix: str,
    read_ids: list[str],
    snp_ids: list[str],
    sigma: np.ndarray,
    delta: np.ndarray,
    score: float,
    total_weight: float,
) -> None:
    with open(prefix + ".read_haplotags.tsv", "w") as out:
        out.write("read_id\tsigma\n")
        for read_id, value in zip(read_ids, sigma):
            out.write(f"{read_id}\t{int(value)}\n")
    with open(prefix + ".snp_phases.tsv", "w") as out:
        out.write("snp_id\tdelta\n")
        for snp_id, value in zip(snp_ids, delta):
            out.write(f"{snp_id}\t{int(value)}\n")
    with open(prefix + ".summary.tsv", "w") as out:
        out.write("metric\tvalue\n")
        out.write("method\tmsf_backbone\n")
        out.write(f"reads\t{len(read_ids)}\n")
        out.write(f"snps\t{len(snp_ids)}\n")
        out.write(f"best_score\t{score:.8f}\n")
        out.write(f"max_possible_score\t{total_weight:.8f}\n")
        out.write(f"weighted_disagreement\t{(total_weight - score) / 2.0:.8f}\n")


def main() -> None:
    parser = argparse.ArgumentParser(description="Build Phirefly's MSF SNP backbone.")
    parser.add_argument("--observations", required=True)
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument("--weighted", action="store_true")
    parser.add_argument("--weight-column", default="weight")
    args = parser.parse_args()

    obs = load_observations(args.observations, args.weighted, args.weight_column)
    read_ids, snp_ids, _matrix, sigma, delta, score, total_weight = solve_msf_backbone(obs, weighted=args.weighted)
    write_solution(args.out_prefix, read_ids, snp_ids, sigma, delta, score, total_weight)
    print(f"reads\t{len(read_ids)}")
    print(f"snps\t{len(snp_ids)}")
    print(f"best_score\t{score:.8f}")
    print(f"summary\t{args.out_prefix}.summary.tsv")


if __name__ == "__main__":
    main()


__all__ = [
    "build_msf_backbone",
    "build_snp_parity_edges",
    "group_observations_by_read",
    "infer_read_labels",
    "main",
    "mst_delta_from_observations",
    "score_read_snp_objective",
    "solve_msf_backbone",
    "write_solution",
]
