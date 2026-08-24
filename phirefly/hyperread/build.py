"""Hyperread construction for compressed Phirefly graphs.

Single-phaselet reads are profiled out as constants. Reads spanning two or
more phaselets are collapsed by phaselet path and sign pattern.
"""

from __future__ import annotations

import hashlib
import math
from collections import defaultdict
from dataclasses import dataclass

import numpy as np
from scipy import sparse

ReadEntry = tuple[int, int, float]
PhaseletVector = list[tuple[int, float]]
HyperreadKey = tuple[tuple[int, ...], tuple[int, ...]]


@dataclass(frozen=True)
class ProfiledReads:
    single_count: int
    single_constant: float
    multi_vectors: list[PhaseletVector]


def build_read_entries(obs, read_ids: list[str], snp_ids: list[str]):
    read_index = {read_id: idx for idx, read_id in enumerate(read_ids)}
    snp_index = {snp_id: idx for idx, snp_id in enumerate(snp_ids)}
    entries: list[list[ReadEntry]] = [[] for _ in read_ids]
    for read_id, snp_id, b_ki, weight in obs:
        entries[read_index[read_id]].append((snp_index[snp_id], int(b_ki), float(weight)))
    for row in entries:
        row.sort(key=lambda item: item[0])
    return entries


def stable_fraction(label: str, salt: str) -> float:
    digest = hashlib.blake2b(f"{salt}:{label}".encode(), digest_size=8).digest()
    return int.from_bytes(digest, byteorder="little", signed=False) / float(1 << 64)


def split_read_entries(
    read_ids: list[str],
    read_entries: list[list[ReadEntry]],
    validation_fraction: float,
    validation_salt: str,
) -> tuple[list[list[ReadEntry]], np.ndarray, np.ndarray]:
    if validation_fraction <= 0.0:
        all_mask = np.ones(len(read_ids), dtype=bool)
        none_mask = np.zeros(len(read_ids), dtype=bool)
        return read_entries, all_mask, none_mask
    if validation_fraction >= 0.5:
        raise ValueError("--validation-fraction must be < 0.5 so the training graph remains well supported")
    val_mask = np.fromiter(
        (stable_fraction(read_id, validation_salt) < validation_fraction for read_id in read_ids),
        dtype=bool,
        count=len(read_ids),
    )
    train_mask = ~val_mask
    if not np.any(train_mask):
        raise ValueError("validation split left no training reads")
    train_entries = [entries for entries, keep in zip(read_entries, train_mask) if keep]
    return train_entries, train_mask, val_mask


def project_read_to_phaselets(
    entries: list[ReadEntry],
    tau_vec: np.ndarray,
    phaselet_ids: np.ndarray,
) -> PhaseletVector:
    by_phaselet: dict[int, float] = defaultdict(float)
    for snp_idx, b_ki, weight in entries:
        phaselet = int(phaselet_ids[snp_idx])
        by_phaselet[phaselet] += float(weight) * int(b_ki) * int(tau_vec[snp_idx])
    return [(phaselet, value) for phaselet, value in sorted(by_phaselet.items()) if value != 0.0]


def project_reads_to_phaselets(
    read_entries: list[list[ReadEntry]],
    tau_vec: np.ndarray,
    phaselet_ids: np.ndarray,
) -> list[PhaseletVector]:
    """Project read-SNP observations into read-by-phaselet signed support vectors."""

    return [project_read_to_phaselets(entries, tau_vec, phaselet_ids) for entries in read_entries]


def profile_single_phaselet_reads(vectors: list[PhaseletVector]) -> ProfiledReads:
    """Remove single-phaselet reads exactly and keep only bridge-read vectors."""

    single_count = 0
    single_constant = 0.0
    multi_vectors: list[PhaseletVector] = []
    for values in vectors:
        if len(values) <= 1:
            if values:
                single_count += 1
                single_constant += abs(values[0][1])
            continue
        multi_vectors.append(values)
    return ProfiledReads(single_count, single_constant, multi_vectors)


def canonicalize_hyperread_pattern(values: PhaseletVector) -> tuple[HyperreadKey, list[tuple[int, float]]]:
    """Canonicalize global read sign so equivalent patterns collapse together."""

    path = tuple(phaselet for phaselet, _ in values)
    signs = [1 if value > 0 else -1 for _, value in values]
    global_sign = signs[0]
    canonical_signs = tuple(sign * global_sign for sign in signs)
    canonical_values = [(phaselet, global_sign * value) for phaselet, value in values]
    return (path, canonical_signs), canonical_values


def collapse_hyperreads(vectors: list[PhaseletVector]):
    """Collapse multi-phaselet reads by phaselet path and canonical sign pattern."""

    hyper: dict[HyperreadKey, dict[int, float]] = {}
    hyper_support: dict[HyperreadKey, int] = defaultdict(int)
    for values in vectors:
        key, canonical_values = canonicalize_hyperread_pattern(values)
        if key not in hyper:
            hyper[key] = defaultdict(float)
        for phaselet, value in canonical_values:
            hyper[key][phaselet] += value
        hyper_support[key] += 1
    return hyper, hyper_support


def normalize_hyperread_value(
    value: float,
    support_reads: int,
    hyperread_degree: int,
    phaselet_degree: int,
    domination_fraction: float,
    hyperread_norm: str,
    hyperread_cap: float,
    adaptive_dom_threshold: float,
) -> tuple[float, int]:
    capped_by_dom = 0
    if hyperread_norm == "mean":
        value = value / max(1, support_reads)
    elif hyperread_norm == "sqrt":
        value = value / math.sqrt(max(1, support_reads))
    elif hyperread_norm == "cap":
        value = float(np.clip(value, -float(hyperread_cap), float(hyperread_cap)))
    elif hyperread_norm == "sqrt_cap":
        value = value / math.sqrt(max(1, support_reads))
        value = float(np.clip(value, -float(hyperread_cap), float(hyperread_cap)))
    elif hyperread_norm == "degree_norm_cap":
        denom = math.sqrt(max(1, hyperread_degree) * max(1, phaselet_degree))
        value = value / denom
        value = float(np.clip(value, -float(hyperread_cap), float(hyperread_cap)))
    elif hyperread_norm == "adaptive_dom_cap":
        if domination_fraction >= float(adaptive_dom_threshold):
            value = float(np.clip(value, -float(hyperread_cap), float(hyperread_cap)))
            capped_by_dom = 1
    elif hyperread_norm != "sum":
        raise ValueError(f"unknown hyperread_norm: {hyperread_norm}")
    return value, capped_by_dom


def build_hyperread_matrix(
    read_entries: list[list[ReadEntry]],
    tau_vec: np.ndarray,
    phaselet_ids: np.ndarray,
    n_phaselets: int,
    hyperread_norm: str = "sum",
    hyperread_cap: float = 10.0,
    adaptive_dom_threshold: float = 0.30,
):
    vectors = project_reads_to_phaselets(read_entries, tau_vec, phaselet_ids)
    profiled = profile_single_phaselet_reads(vectors)
    hyper, hyper_support = collapse_hyperreads(profiled.multi_vectors)
    sorted_keys = sorted(hyper, key=lambda item: (item[0], item[1]))

    phaselet_degree: dict[int, int] = defaultdict(int)
    phaselet_total_abs: dict[int, float] = defaultdict(float)
    for key in sorted_keys:
        for phaselet, value in hyper[key].items():
            if value == 0.0:
                continue
            phaselet_degree[phaselet] += 1
            phaselet_total_abs[phaselet] += abs(value)

    rows = []
    cols = []
    data = []
    hyperread_rows = []
    hyperread_edge_rows = []
    entry_domination_scores = []
    for row_idx, key in enumerate(sorted_keys):
        path, signs = key
        support_reads = hyper_support[key]
        raw_abs_weight = float(sum(abs(value) for value in hyper[key].values()))
        hyperread_degree = sum(1 for value in hyper[key].values() if value != 0.0)
        normalized_abs_weight = 0.0
        hyperread_rows.append(
            {
                "hyperread_id": row_idx,
                "path": ",".join(map(str, path)),
                "signs": ",".join(map(str, signs)),
                "support_reads": support_reads,
                "degree": hyperread_degree,
                "raw_abs_weight": f"{raw_abs_weight:.8f}",
            }
        )
        for phaselet, raw_value in sorted(hyper[key].items()):
            if raw_value == 0.0:
                continue
            total_abs = phaselet_total_abs[phaselet]
            dom = abs(raw_value) / total_abs if total_abs > 0 else 1.0
            entry_domination_scores.append(dom)
            value, capped_by_dom = normalize_hyperread_value(
                raw_value,
                support_reads=support_reads,
                hyperread_degree=hyperread_degree,
                phaselet_degree=phaselet_degree[phaselet],
                domination_fraction=dom,
                hyperread_norm=hyperread_norm,
                hyperread_cap=hyperread_cap,
                adaptive_dom_threshold=adaptive_dom_threshold,
            )
            normalized_abs_weight += abs(value)
            rows.append(row_idx)
            cols.append(phaselet)
            data.append(value)
            hyperread_edge_rows.append(
                {
                    "hyperread_id": row_idx,
                    "phaselet_id": phaselet,
                    "raw_value": f"{raw_value:.8f}",
                    "normalized_value": f"{value:.8f}",
                    "abs_raw_value": f"{abs(raw_value):.8f}",
                    "abs_normalized_value": f"{abs(value):.8f}",
                    "support_reads": support_reads,
                    "hyperread_degree": hyperread_degree,
                    "phaselet_degree": phaselet_degree[phaselet],
                    "phaselet_total_abs_raw": f"{total_abs:.8f}",
                    "phaselet_domination_fraction": f"{dom:.8f}",
                    "capped_by_domination": capped_by_dom,
                }
            )
        hyperread_rows[-1]["normalized_abs_weight"] = f"{normalized_abs_weight:.8f}"

    matrix = sparse.coo_matrix(
        (
            np.asarray(data, dtype=np.float64),
            (np.asarray(rows, dtype=np.int64), np.asarray(cols, dtype=np.int64)),
        ),
        shape=(len(hyperread_rows), n_phaselets),
    ).tocsr()
    matrix.sum_duplicates()
    stats = {
        "single_phaselet_reads": profiled.single_count,
        "single_phaselet_constant": profiled.single_constant,
        "multi_phaselet_reads": len(profiled.multi_vectors),
        "hyperreads": len(hyperread_rows),
        "max_phaselet_hyperread_domination": (
            f"{float(max(entry_domination_scores)):.8f}" if entry_domination_scores else "0.00000000"
        ),
        "mean_phaselet_hyperread_domination": (
            f"{float(np.mean(entry_domination_scores)):.8f}" if entry_domination_scores else "0.00000000"
        ),
    }
    return matrix, hyperread_rows, hyperread_edge_rows, stats


__all__ = [
    "ProfiledReads",
    "build_hyperread_matrix",
    "build_read_entries",
    "canonicalize_hyperread_pattern",
    "collapse_hyperreads",
    "normalize_hyperread_value",
    "profile_single_phaselet_reads",
    "project_read_to_phaselets",
    "project_reads_to_phaselets",
    "split_read_entries",
    "stable_fraction",
]
