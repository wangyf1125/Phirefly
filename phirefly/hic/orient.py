#!/usr/bin/env python3
"""Orient Phirefly components from signed Hi-C/Pore-C edges using MSF only."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from ..core.union_find import UnionFind


def read_edges(path: str | Path) -> list[dict[str, str]]:
    with Path(path).open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def retained_edges(rows: list[dict[str, str]], min_contacts: int, min_confidence: float, min_abs_score: float):
    out = []
    for row in rows:
        contacts = int(float(row["contacts"]))
        confidence = float(row["confidence"])
        score = float(row["score"])
        if contacts < min_contacts or confidence < min_confidence or abs(score) < min_abs_score:
            continue
        out.append((row["component_i"], row["component_j"], 1 if score >= 0 else -1, abs(score), row))
    return out


def _component_sort_key(component_id: str) -> tuple[int, str]:
    try:
        return int(component_id), component_id
    except ValueError:
        return 10**18, component_id


def orient_components_msf(
    edges: list[tuple[str, str, int, float, dict[str, str]]],
) -> tuple[dict[str, int], dict[str, str], list[dict[str, str]]]:
    """Maximum-spanning forest over signed component edges."""

    uf = UnionFind()
    selected = []
    for comp_i, comp_j, sign, weight, row in sorted(edges, key=lambda item: item[3], reverse=True):
        if uf.find(comp_i) == uf.find(comp_j):
            continue
        uf.union(comp_i, comp_j)
        selected.append((comp_i, comp_j, sign, row))

    adjacency: dict[str, list[tuple[str, int]]] = {}
    for comp_i, comp_j, sign, _row in selected:
        adjacency.setdefault(comp_i, []).append((comp_j, sign))
        adjacency.setdefault(comp_j, []).append((comp_i, sign))

    orientation: dict[str, int] = {}
    phase_sets: dict[str, str] = {}
    for node in sorted(adjacency):
        if node in orientation:
            continue
        orientation[node] = 1
        stack = [node]
        members = []
        while stack:
            cur = stack.pop()
            members.append(cur)
            for nxt, sign in adjacency.get(cur, []):
                value = orientation[cur] * sign
                if nxt not in orientation:
                    orientation[nxt] = value
                    stack.append(nxt)
        phase_set = min(members, key=_component_sort_key)
        for member in members:
            phase_sets[member] = phase_set

    selected_rows = []
    for comp_i, comp_j, sign, row in selected:
        out = dict(row)
        out["selected"] = 1
        out["edge_sign"] = sign
        selected_rows.append(out)
    return orientation, phase_sets, selected_rows


def write_orientations(path: Path, orientation: dict[str, int], phase_sets: dict[str, str] | None = None) -> None:
    phase_sets = phase_sets or {}
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=["component_id", "orientation", "hic_ps"])
        writer.writeheader()
        for component_id in sorted(orientation, key=lambda x: (len(x), x)):
            writer.writerow(
                {
                    "component_id": component_id,
                    "orientation": int(orientation[component_id]),
                    "hic_ps": phase_sets.get(component_id, component_id),
                }
            )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--edges", required=True)
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument("--solver", choices=["msf"], default="msf")
    parser.add_argument("--min-contacts", type=int, default=20)
    parser.add_argument("--min-confidence", type=float, default=0.8)
    parser.add_argument("--min-abs-score", type=float, default=0.0)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    rows = read_edges(args.edges)
    kept = retained_edges(rows, args.min_contacts, args.min_confidence, args.min_abs_score)
    orientation, phase_sets, selected = orient_components_msf(kept)
    prefix = Path(args.out_prefix)
    out_path = Path(f"{prefix}.component_orientations.tsv")
    write_orientations(out_path, orientation, phase_sets)
    with Path(f"{prefix}.summary.tsv").open("w") as out:
        out.write("metric\tvalue\n")
        out.write("solver\tmsf\n")
        out.write(f"input_edges\t{len(rows)}\n")
        out.write(f"retained_edges\t{len(kept)}\n")
        out.write(f"selected_msf_edges\t{len(selected)}\n")
        out.write(f"oriented_components\t{len(orientation)}\n")
    if selected:
        fieldnames = list(selected[0])
        with Path(f"{prefix}.selected_edges.tsv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(selected)
    print(out_path)


if __name__ == "__main__":
    main()


__all__ = ["orient_components_msf", "read_edges", "retained_edges", "write_orientations"]
