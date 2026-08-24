"""Adjacent SNP-SNP parity edge writers generated during extraction."""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

from .sites import SnpSite
from .store import open_text

ObservationRow = tuple[int, int, int, int, int, float, float]


def build_adjacent_parity_edges(
    rows: list[ObservationRow],
    sites: list[SnpSite],
    max_distance_bp: int = 0,
    weighted: bool = False,
) -> list[dict[str, object]]:
    by_read: dict[int, list[tuple[int, int, float]]] = defaultdict(list)
    for read_idx, snp_idx, allele, _baseq, _mapq, _eps, weight in rows:
        by_read[int(read_idx)].append((int(snp_idx), int(allele), float(weight)))

    edge_score: dict[tuple[int, int], float] = defaultdict(float)
    edge_abs: dict[tuple[int, int], float] = defaultdict(float)
    edge_count: dict[tuple[int, int], int] = defaultdict(int)
    edge_max_abs: dict[tuple[int, int], float] = defaultdict(float)

    for entries in by_read.values():
        if len(entries) < 2:
            continue
        entries.sort(key=lambda item: item[0])
        for left, right in zip(entries, entries[1:]):
            i, b_i, w_i = left
            j, b_j, w_j = right
            if i == j:
                continue
            if max_distance_bp > 0 and abs(sites[j].pos0 - sites[i].pos0) > max_distance_bp:
                continue
            if i > j:
                i, j = j, i
                b_i, b_j = b_j, b_i
                w_i, w_j = w_j, w_i
            gamma = min(float(w_i), float(w_j)) if weighted else 1.0
            key = (i, j)
            edge_score[key] += gamma * int(b_i) * int(b_j)
            edge_abs[key] += abs(gamma)
            edge_count[key] += 1
            edge_max_abs[key] = max(edge_max_abs[key], abs(gamma))

    edge_rows = []
    for (i, j), score in sorted(edge_score.items()):
        support = edge_abs[(i, j)]
        confidence = abs(score) / support if support > 0 else 0.0
        edge_rows.append(
            {
                "snp_i": i,
                "snp_j": j,
                "snp_id_i": sites[i].snp_id,
                "snp_id_j": sites[j].snp_id,
                "chrom": sites[i].chrom,
                "pos_i": sites[i].pos0 + 1,
                "pos_j": sites[j].pos0 + 1,
                "distance_bp": abs(sites[j].pos0 - sites[i].pos0),
                "score": f"{score:.8f}",
                "support": f"{support:.8f}",
                "count": edge_count[(i, j)],
                "confidence": f"{confidence:.8f}",
                "max_abs_contribution": f"{edge_max_abs[(i, j)]:.8f}",
            }
        )
    return edge_rows


def write_parity_edges(path: Path, edge_rows: list[dict[str, object]]) -> None:
    fieldnames = [
        "snp_i",
        "snp_j",
        "snp_id_i",
        "snp_id_j",
        "chrom",
        "pos_i",
        "pos_j",
        "distance_bp",
        "score",
        "support",
        "count",
        "confidence",
        "max_abs_contribution",
    ]
    with open_text(path) as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(edge_rows)


__all__ = ["ObservationRow", "build_adjacent_parity_edges", "write_parity_edges"]
