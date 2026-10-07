#!/usr/bin/env python3
"""Build component-component orientation edges from Hi-C/Pore-C alignments."""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import pysam

from .contact import alignment_anchor


def remap_site_chroms(site_map, site_positions, aliases: dict[str, str]):
    if not aliases:
        return site_map, site_positions
    remapped = {}
    for chrom, by_pos in site_map.items():
        remapped[aliases.get(chrom, chrom)] = by_pos
    return remapped, {chrom: sorted(by_pos) for chrom, by_pos in remapped.items()}


def collect_anchors(contact_bam: str | Path, region: tuple[str, int, int], site_map, site_positions, args):
    anchors_by_name: dict[str, list[dict]] = defaultdict(list)
    chrom, start0, end0 = region
    with pysam.AlignmentFile(str(contact_bam), "rb", threads=getattr(args, 'threads', 1)) as bam:
        if not bam.has_index():
            raise ValueError(f'Hi-C BAM must be indexed: {contact_bam}')
        for aln in bam.fetch(chrom, start0, end0):
            anchor = alignment_anchor(aln, site_map, site_positions, args)
            if anchor is not None:
                anchors_by_name[anchor["query_name"]].append(anchor)
    return anchors_by_name


def component_edges_from_anchors(anchors_by_name: dict[str, list[dict]], alpha: float) -> list[dict[str, object]]:
    accum: dict[tuple[str, str], dict[str, float | int]] = {}
    for anchors in anchors_by_name.values():
        by_component = {}
        for anchor in anchors:
            component = str(anchor["component_id"])
            previous = by_component.get(component)
            if previous is None or float(anchor["score"]) > float(previous["score"]):
                by_component[component] = anchor
        items = sorted(by_component.values(), key=lambda row: str(row["component_id"]))
        for i, left in enumerate(items):
            for right in items[i + 1 :]:
                comp_i = str(left["component_id"])
                comp_j = str(right["component_id"])
                key = (comp_i, comp_j) if comp_i < comp_j else (comp_j, comp_i)
                sign = int(left["sign"]) * int(right["sign"])
                weight = max(float(alpha), min(float(left["score"]), float(right["score"])))
                if key not in accum:
                    accum[key] = {
                        "score": 0.0,
                        "contacts": 0,
                        "same_support": 0.0,
                        "opposite_support": 0.0,
                    }
                row = accum[key]
                row["score"] = float(row["score"]) + sign * weight
                row["contacts"] = int(row["contacts"]) + 1
                if sign > 0:
                    row["same_support"] = float(row["same_support"]) + weight
                else:
                    row["opposite_support"] = float(row["opposite_support"]) + weight

    rows = []
    for (comp_i, comp_j), values in sorted(accum.items()):
        same = float(values["same_support"])
        opposite = float(values["opposite_support"])
        total = same + opposite
        score = float(values["score"])
        rows.append(
            {
                "component_i": comp_i,
                "component_j": comp_j,
                "score": f"{score:.8f}",
                "abs_score": f"{abs(score):.8f}",
                "contacts": int(values["contacts"]),
                "same_support": f"{same:.8f}",
                "opposite_support": f"{opposite:.8f}",
                "confidence": f"{abs(score) / total if total else 0.0:.8f}",
            }
        )
    return rows


def write_rows(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


__all__ = ["component_edges_from_anchors", "collect_anchors", "remap_site_chroms"]
