#!/usr/bin/env python3
"""Truth-independent QC for Hi-C/Pore-C component edges."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_edges(path: str | Path) -> list[dict[str, str]]:
    with Path(path).open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def edge_summary(rows: list[dict[str, str]], min_contacts: int, min_confidence: float, min_abs_score: float) -> dict[str, object]:
    retained = []
    for row in rows:
        contacts = int(float(row["contacts"]))
        confidence = float(row["confidence"])
        abs_score = float(row["abs_score"])
        if contacts >= min_contacts and confidence >= min_confidence and abs_score >= min_abs_score:
            retained.append(row)
    confidences = [float(row["confidence"]) for row in retained]
    contacts = [int(float(row["contacts"])) for row in retained]
    return {
        "input_edges": len(rows),
        "retained_edges": len(retained),
        "retained_fraction": len(retained) / len(rows) if rows else 0.0,
        "mean_confidence": sum(confidences) / len(confidences) if confidences else 0.0,
        "mean_contacts": sum(contacts) / len(contacts) if contacts else 0.0,
        "max_contacts": max(contacts) if contacts else 0,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--edges", required=True)
    parser.add_argument("--min-contacts", type=int, default=20)
    parser.add_argument("--min-confidence", type=float, default=0.8)
    parser.add_argument("--min-abs-score", type=float, default=0.0)
    parser.add_argument("--out-tsv", default=None)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    summary = edge_summary(read_edges(args.edges), args.min_contacts, args.min_confidence, args.min_abs_score)
    if args.out_tsv:
        with Path(args.out_tsv).open("w") as out:
            out.write("metric\tvalue\n")
            for key, value in summary.items():
                out.write(f"{key}\t{value}\n")
        print(args.out_tsv)
    else:
        for key, value in summary.items():
            print(f"{key}\t{value}")


if __name__ == "__main__":
    main()


__all__ = ["edge_summary", "main", "read_edges"]
