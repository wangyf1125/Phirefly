#!/usr/bin/env python3
"""Command-line read-SNP observation extraction."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from .parity import build_adjacent_parity_edges, write_parity_edges
from .bam import extract_observations_chunked
from .sites import SnpSite, load_het_snps
from ..core.io import open_text
from .store import ordered_read_names, write_npz_observations, write_sidecar_tables


def extract_observations(
    bam_path: str,
    sample: str,
    sites: list[SnpSite],
    out_prefix: Path,
    min_baseq: int,
    min_mapq: int,
    epsilon_floor: float,
    epsilon_ceiling: float,
    mapq_cap: int,
    default_baseq: int | None,
    scan_mode: str,
    output_format: str,
    threads: int,
    chunk_size_bp: int,
    write_parity_edges_flag: bool,
    parity_max_distance_bp: int,
    parity_weighted: bool,
) -> None:
    obs_path = out_prefix.with_suffix(".observations.tsv.gz")
    npz_path = out_prefix.with_suffix(".observations.npz")
    parity_path = out_prefix.with_suffix(".adjacent_parity_edges.tsv.gz")

    read_ids, rows, n_ref, n_alt, n_other, n_chunks = extract_observations_chunked(
        bam_path,
        sites=sites,
        min_baseq=min_baseq,
        min_mapq=min_mapq,
        epsilon_floor=epsilon_floor,
        epsilon_ceiling=epsilon_ceiling,
        mapq_cap=mapq_cap,
        default_baseq=default_baseq,
        scan_mode=scan_mode,
        threads=threads,
        chunk_size_bp=chunk_size_bp,
    )

    if output_format in {"tsv", "both"}:
        ordered_reads = ordered_read_names(read_ids)
        with open_text(obs_path, "wt") as out:
            out.write("read_id\tsnp_id\tb_ki\tbaseq\tmapq\tepsilon\tweight\n")
            for read_idx, snp_idx, allele, baseq, mapq, eps, weight in rows:
                out.write(
                    f"{ordered_reads[read_idx]}\t{sites[snp_idx].snp_id}\t{allele}\t{baseq}\t{mapq}"
                    f"\t{eps:.8g}\t{weight:.8g}\n"
                )

    if output_format in {"npz", "both"}:
        write_npz_observations(
            npz_path,
            sites=sites,
            read_ids=read_ids,
            read_idx=[row[0] for row in rows],
            snp_idx=[row[1] for row in rows],
            b_ki=[row[2] for row in rows],
            baseq=[row[3] for row in rows],
            mapq=[row[4] for row in rows],
            epsilon=[row[5] for row in rows],
            weight=[row[6] for row in rows],
        )

    parity_edges = []
    if write_parity_edges_flag:
        parity_edges = build_adjacent_parity_edges(
            rows,
            sites,
            max_distance_bp=parity_max_distance_bp,
            weighted=parity_weighted,
        )
        write_parity_edges(parity_path, parity_edges)

    snp_path, read_path = write_sidecar_tables(out_prefix, sites, read_ids)
    summary = dict(bam=bam_path, sample=sample, region=sites[0].chrom if sites else 'NA',
                   scan_mode=scan_mode, output_format=output_format, threads=threads,
                   chunks=n_chunks, chunk_size_bp=chunk_size_bp, snps=len(sites),
                   reads=len(read_ids), observations=len(rows), ref_observations=n_ref,
                   alt_observations=n_alt, other_base_skipped=n_other,
                   observations_tsv=obs_path if output_format in {'tsv', 'both'} else '',
                   observations_npz=npz_path if output_format in {'npz', 'both'} else '',
                   parity_edges=len(parity_edges), parity_weighted=int(parity_weighted),
                   parity_edges_tsv=parity_path if write_parity_edges_flag else '',
                   snps_tsv=snp_path, reads_tsv=read_path)
    summary_path = out_prefix.with_suffix('.extract_summary.tsv')
    with summary_path.open('w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(('metric', 'value'))
        writer.writerows(summary.items())
    print(f"Extracted {len(rows):,} observations from {len(read_ids):,} reads; details: {summary_path}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", required=True)
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--sample", default=None, help="Sample to extract. Default: first sample in VCF.")
    parser.add_argument("--region", default="chr11:1-10000000")
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument("--min-baseq", type=int, default=13)
    parser.add_argument("--min-mapq", type=int, default=20)
    parser.add_argument("--epsilon-floor", type=float, default=1e-4)
    parser.add_argument("--epsilon-ceiling", type=float, default=0.49)
    parser.add_argument("--mapq-cap", type=int, default=60)
    parser.add_argument("--scan-mode", choices=["read", "pileup"], default="read")
    parser.add_argument("--output-format", choices=["npz", "tsv", "both"], default="both")
    parser.add_argument("--default-baseq", type=int, default=None)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--chunk-size-bp", type=int, default=0)
    parser.add_argument("--write-parity-edges", action="store_true")
    parser.add_argument("--parity-max-distance-bp", type=int, default=0)
    parser.add_argument("--parity-weighted", action="store_true")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    sample, sites = load_het_snps(args.vcf, args.region, args.sample)
    if not sites:
        raise SystemExit(f"No heterozygous biallelic SNVs found for {sample} in {args.region}")
    extract_observations(
        bam_path=args.bam,
        sample=sample,
        sites=sites,
        out_prefix=Path(args.out_prefix),
        min_baseq=args.min_baseq,
        min_mapq=args.min_mapq,
        epsilon_floor=args.epsilon_floor,
        epsilon_ceiling=args.epsilon_ceiling,
        mapq_cap=args.mapq_cap,
        default_baseq=args.default_baseq,
        scan_mode=args.scan_mode,
        output_format=args.output_format,
        threads=max(1, int(args.threads)),
        chunk_size_bp=max(0, int(args.chunk_size_bp)),
        write_parity_edges_flag=args.write_parity_edges,
        parity_max_distance_bp=max(0, int(args.parity_max_distance_bp)),
        parity_weighted=args.parity_weighted,
    )


if __name__ == "__main__":
    main()
