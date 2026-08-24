"""Read-SNP observation loading and sparse matrix construction."""

from __future__ import annotations

import gzip

import numpy as np
from scipy import sparse


def open_maybe_gzip(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def load_observations(path: str, weighted: bool, weight_column: str):
    """Load text TSV(.gz) or binary NPZ observations."""

    if path.endswith(".npz"):
        return load_npz_observations(path, weighted=weighted)
    obs = []
    with open_maybe_gzip(path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        if weighted and weight_column not in idx:
            raise SystemExit(f"--weighted requested but column {weight_column!r} is absent from {path}")
        for line in handle:
            if not line.strip():
                continue
            row = line.rstrip("\n").split("\t")
            weight = float(row[idx[weight_column]]) if weighted else 1.0
            obs.append((row[idx["read_id"]], row[idx["snp_id"]], int(row[idx["b_ki"]]), weight))
    return obs


def load_npz_observations(path: str, weighted: bool):
    arr = np.load(path, allow_pickle=False)
    read_ids = arr["read_ids"]
    snp_ids = arr["snp_ids"]
    read_idx = arr["read_idx"]
    snp_idx = arr["snp_idx"]
    b_ki = arr["b_ki"]
    weights = arr["weight"] if weighted and "weight" in arr.files else np.ones_like(b_ki, dtype=np.float64)
    return [
        (str(read_ids[int(read_i)]), str(snp_ids[int(snp_i)]), int(b), float(weight))
        for read_i, snp_i, b, weight in zip(read_idx, snp_idx, b_ki, weights)
    ]


def build_matrix(obs):
    read_ids = sorted({read_id for read_id, _, _, _ in obs})
    snp_ids = sorted({snp_id for _, snp_id, _, _ in obs}, key=lambda x: (x.split(":")[0], int(x.split(":")[1])))
    read_index = {read_id: i for i, read_id in enumerate(read_ids)}
    snp_index = {snp_id: i for i, snp_id in enumerate(snp_ids)}
    rows = np.fromiter((read_index[read_id] for read_id, _, _, _ in obs), dtype=np.int64, count=len(obs))
    cols = np.fromiter((snp_index[snp_id] for _, snp_id, _, _ in obs), dtype=np.int64, count=len(obs))
    data = np.fromiter((weight * allele for _, _, allele, weight in obs), dtype=np.float64, count=len(obs))
    matrix = sparse.coo_matrix((data, (rows, cols)), shape=(len(read_ids), len(snp_ids))).tocsr()
    matrix.sum_duplicates()
    return read_ids, snp_ids, matrix


def load_npz_graph(path: str, weighted: bool):
    """Load NPZ observations directly into a sparse read-SNP matrix."""

    arr = np.load(path, allow_pickle=False)
    read_ids = [str(value) for value in arr["read_ids"]]
    snp_ids = [str(value) for value in arr["snp_ids"]]
    rows = np.asarray(arr["read_idx"], dtype=np.int64)
    cols = np.asarray(arr["snp_idx"], dtype=np.int64)
    alleles = np.asarray(arr["b_ki"], dtype=np.int8)
    weights = np.asarray(arr["weight"], dtype=np.float64) if weighted and "weight" in arr.files else np.ones(len(rows))
    data = weights * alleles
    matrix = sparse.coo_matrix((data, (rows, cols)), shape=(len(read_ids), len(snp_ids))).tocsr()
    matrix.sum_duplicates()
    read_entries: list[list[tuple[int, int, float]]] = [[] for _ in read_ids]
    for read_idx, snp_idx, allele, weight in zip(rows, cols, alleles, weights):
        read_entries[int(read_idx)].append((int(snp_idx), int(allele), float(weight)))
    for entries in read_entries:
        entries.sort(key=lambda item: item[0])
    return read_ids, snp_ids, matrix, read_entries, int(len(rows))


__all__ = [
    "build_matrix",
    "load_npz_graph",
    "load_npz_observations",
    "load_observations",
    "open_maybe_gzip",
]
