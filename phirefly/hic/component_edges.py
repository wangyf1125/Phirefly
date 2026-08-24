#!/usr/bin/env python3
"""Build component-component orientation edges from Hi-C/Pore-C alignments."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import pysam

from ..core.genomics import apply_region_alias, parse_chrom_alias, parse_region
from .components import load_component_site_map
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
    with pysam.AlignmentFile(str(contact_bam), "rb") as bam:
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


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--observations", required=True)
    parser.add_argument("--haplotags", default=None, help="Accepted for run provenance; not used by this edge builder.")
    parser.add_argument("--snp-phases", required=True)
    parser.add_argument("--contact-bam", required=True)
    parser.add_argument("--region", required=True)
    parser.add_argument("--contact-region", default=None)
    parser.add_argument("--chrom-alias", action="append", default=[])
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument("--min-mapq", type=int, default=20)
    parser.add_argument("--mapq-cap", type=int, default=60)
    parser.add_argument("--min-sites", type=int, default=1)
    parser.add_argument("--min-local-score", type=float, default=1.0)
    parser.add_argument("--alpha", type=float, default=0.1)
    parser.add_argument("--threads", type=int, default=1, help="Reserved for future chunk parallelism.")
    parser.add_argument("--chunk-size-bp", type=int, default=0, help="Reserved for future chunk parallelism.")
    parser.add_argument("--keep-secondary", action="store_true")
    parser.add_argument("--keep-supplementary", action="store_true")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    aliases = parse_chrom_alias(args.chrom_alias)
    site_map, site_positions = load_component_site_map(args.observations, args.snp_phases)
    site_map, site_positions = remap_site_chroms(site_map, site_positions, aliases)
    base_region = parse_region(args.region, zero_based_start=True)
    if base_region is None:
        raise SystemExit("--region is required")
    contact_region = (
        parse_region(args.contact_region, zero_based_start=True)
        if args.contact_region
        else apply_region_alias(base_region, aliases)
    )
    if contact_region is None:
        raise SystemExit("--contact-region could not be resolved")

    anchors = collect_anchors(args.contact_bam, contact_region, site_map, site_positions, args)
    edge_rows = component_edges_from_anchors(anchors, alpha=float(args.alpha))
    prefix = Path(args.out_prefix)
    edge_path = prefix.with_suffix(".block_edges.tsv")
    write_rows(
        edge_path,
        edge_rows,
        ["component_i", "component_j", "score", "abs_score", "contacts", "same_support", "opposite_support", "confidence"],
    )
    manifest = prefix.with_suffix(".manifest.tsv")
    with manifest.open("w") as out:
        out.write("key\tvalue\n")
        out.write(f"component_edges\t{len(edge_rows)}\n")
        out.write(f"contact_molecules\t{len(anchors)}\n")
        out.write(f"region\t{args.region}\n")
        out.write(f"contact_region\t{contact_region[0]}:{contact_region[1] + 1}-{contact_region[2]}\n")
        out.write(f"edge_tsv\t{edge_path}\n")
    print(manifest)


if __name__ == "__main__":
    main()


__all__ = ["component_edges_from_anchors", "collect_anchors", "main", "remap_site_chroms"]
