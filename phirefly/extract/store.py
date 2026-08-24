"""Observation store readers and writers."""

from __future__ import annotations

import csv
import gzip
from collections.abc import Iterator
from pathlib import Path

import numpy as np

from .sites import SnpSite


def open_text(path: Path | str, mode: str = "wt"):
    path = Path(path)
    return gzip.open(path, mode) if path.suffix == ".gz" else path.open(mode)


def ordered_read_names(read_ids: dict[str, int]) -> list[str]:
    return [read_id for read_id, _ in sorted(read_ids.items(), key=lambda item: item[1])]


def write_sidecar_tables(out_prefix: Path, sites: list[SnpSite], read_ids: dict[str, int]) -> tuple[Path, Path]:
    snp_path = out_prefix.with_suffix(".snps.tsv")
    read_path = out_prefix.with_suffix(".reads.tsv")
    with snp_path.open("w") as out:
        out.write("snp_id\tchrom\tpos1\tref\talt\n")
        for site in sites:
            out.write(f"{site.snp_id}\t{site.chrom}\t{site.pos0 + 1}\t{site.ref}\t{site.alt}\n")
    with read_path.open("w") as out:
        out.write("read_id\tread_idx\n")
        for read_id, read_idx in sorted(read_ids.items(), key=lambda x: x[1]):
            out.write(f"{read_id}\t{read_idx}\n")
    return snp_path, read_path


def write_npz_observations(
    path: Path,
    sites: list[SnpSite],
    read_ids: dict[str, int],
    read_idx: list[int],
    snp_idx: list[int],
    b_ki: list[int],
    baseq: list[int],
    mapq: list[int],
    epsilon: list[float],
    weight: list[float],
) -> None:
    np.savez_compressed(
        path,
        format_version=np.asarray(["phirefly_observations_v1"], dtype=str),
        read_idx=np.asarray(read_idx, dtype=np.int32),
        snp_idx=np.asarray(snp_idx, dtype=np.int32),
        b_ki=np.asarray(b_ki, dtype=np.int8),
        baseq=np.asarray(baseq, dtype=np.uint8),
        mapq=np.asarray(mapq, dtype=np.uint16),
        epsilon=np.asarray(epsilon, dtype=np.float32),
        weight=np.asarray(weight, dtype=np.float32),
        read_ids=np.asarray(ordered_read_names(read_ids), dtype=str),
        snp_ids=np.asarray([site.snp_id for site in sites], dtype=str),
        chrom=np.asarray([site.chrom for site in sites], dtype=str),
        pos1=np.asarray([site.pos0 + 1 for site in sites], dtype=np.int64),
        ref=np.asarray([site.ref for site in sites], dtype=str),
        alt=np.asarray([site.alt for site in sites], dtype=str),
    )


def iter_observation_pairs(path: str | Path) -> Iterator[tuple[str, str]]:
    """Yield ``(read_id, snp_id)`` from TSV(.gz) or NPZ observations."""

    path = Path(path)
    if path.suffix == ".npz":
        arr = np.load(path, allow_pickle=False)
        read_ids = arr["read_ids"]
        snp_ids = arr["snp_ids"]
        for read_idx, snp_idx in zip(arr["read_idx"], arr["snp_idx"]):
            yield str(read_ids[int(read_idx)]), str(snp_ids[int(snp_idx)])
        return
    with open_text(path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            yield row["read_id"], row["snp_id"]


def count_observations(path: str | Path) -> int:
    path = Path(path)
    if path.suffix == ".npz":
        arr = np.load(path, allow_pickle=False)
        return int(len(arr["read_idx"]))
    with open_text(path, "rt") as handle:
        next(handle, None)
        return sum(1 for line in handle if line.strip())


__all__ = [
    "count_observations",
    "iter_observation_pairs",
    "open_text",
    "ordered_read_names",
    "write_npz_observations",
    "write_sidecar_tables",
]
