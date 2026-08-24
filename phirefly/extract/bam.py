"""BAM to read-SNP observation extraction."""

from __future__ import annotations

import bisect
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass

import pysam

from ..core.alignment import aligned_snp_query_positions
from .parity import ObservationRow
from .sites import SnpSite
from .store import ordered_read_names
from .weights import observation_weight


@dataclass
class ExtractionResult:
    site_indices: list[int]
    read_names: list[str]
    rows: list[ObservationRow]
    n_ref: int
    n_alt: int
    n_other: int


def extract_observations_pileup(
    bam: pysam.AlignmentFile,
    sites: list[SnpSite],
    min_baseq: int,
    min_mapq: int,
    epsilon_floor: float,
    epsilon_ceiling: float,
    mapq_cap: int,
    default_baseq: int | None,
):
    site_by_pos = {site.pos0: (idx, site) for idx, site in enumerate(sites)}
    chrom = sites[0].chrom
    start0 = min(site.pos0 for site in sites)
    end0 = max(site.pos0 for site in sites) + 1
    cache: dict[tuple[int, int], tuple[float, float]] = {}
    read_ids: dict[str, int] = {}
    rows: list[ObservationRow] = []
    n_ref = n_alt = n_other = 0

    for pileup_col in bam.pileup(
        chrom,
        start0,
        end0,
        truncate=True,
        stepper="samtools",
        min_base_quality=min_baseq,
        ignore_overlaps=False,
        ignore_orphans=False,
    ):
        site_record = site_by_pos.get(pileup_col.reference_pos)
        if site_record is None:
            continue
        site_idx, site = site_record
        for pileup_read in pileup_col.pileups:
            aln = pileup_read.alignment
            if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
                continue
            if aln.mapping_quality < min_mapq:
                continue
            if pileup_read.is_del or pileup_read.is_refskip:
                continue
            qpos = pileup_read.query_position
            if qpos is None or aln.query_sequence is None:
                continue
            baseq = int(aln.query_qualities[qpos]) if aln.query_qualities is not None else int(default_baseq or min_baseq)
            if baseq < min_baseq:
                continue
            base = aln.query_sequence[qpos].upper()
            if base == site.ref:
                allele = 1
                n_ref += 1
            elif base == site.alt:
                allele = -1
                n_alt += 1
            else:
                n_other += 1
                continue
            read_idx = read_ids.setdefault(aln.query_name, len(read_ids))
            eps, weight = observation_weight(baseq, aln.mapping_quality, epsilon_floor, epsilon_ceiling, mapq_cap, cache)
            rows.append((read_idx, site_idx, allele, baseq, int(aln.mapping_quality), eps, weight))
    return read_ids, rows, n_ref, n_alt, n_other


def extract_observations_read_scan(
    bam: pysam.AlignmentFile,
    sites: list[SnpSite],
    min_baseq: int,
    min_mapq: int,
    epsilon_floor: float,
    epsilon_ceiling: float,
    mapq_cap: int,
    default_baseq: int | None,
):
    chrom = sites[0].chrom
    positions = [site.pos0 for site in sites]
    site_by_pos = {site.pos0: (idx, site) for idx, site in enumerate(sites)}
    start0 = min(positions)
    end0 = max(positions) + 1
    cache: dict[tuple[int, int], tuple[float, float]] = {}
    read_ids: dict[str, int] = {}
    rows: list[ObservationRow] = []
    n_ref = n_alt = n_other = 0

    for aln in bam.fetch(chrom, start0, end0):
        if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
            continue
        if aln.mapping_quality < min_mapq:
            continue
        if aln.reference_start is None or aln.reference_end is None or aln.query_sequence is None:
            continue
        left = bisect.bisect_left(positions, max(start0, int(aln.reference_start)))
        right = bisect.bisect_left(positions, min(end0, int(aln.reference_end)))
        if left >= right:
            continue
        candidate_positions = positions[left:right]
        read_idx_value = None
        seen_positions: set[int] = set()
        for ref_pos, qpos in aligned_snp_query_positions(aln, candidate_positions):
            if ref_pos in seen_positions:
                continue
            seen_positions.add(ref_pos)
            site_idx, site = site_by_pos[ref_pos]
            baseq = int(aln.query_qualities[qpos]) if aln.query_qualities is not None else int(default_baseq or min_baseq)
            if baseq < min_baseq:
                continue
            base = aln.query_sequence[qpos].upper()
            if base == site.ref:
                allele = 1
                n_ref += 1
            elif base == site.alt:
                allele = -1
                n_alt += 1
            else:
                n_other += 1
                continue
            if read_idx_value is None:
                read_idx_value = read_ids.setdefault(aln.query_name, len(read_ids))
            eps, weight = observation_weight(baseq, aln.mapping_quality, epsilon_floor, epsilon_ceiling, mapq_cap, cache)
            rows.append((read_idx_value, site_idx, allele, baseq, int(aln.mapping_quality), eps, weight))
    return read_ids, rows, n_ref, n_alt, n_other


def chunk_site_indices(sites: list[SnpSite], chunk_size_bp: int) -> list[list[int]]:
    if chunk_size_bp <= 0 or not sites:
        return [list(range(len(sites)))]
    chunks: list[list[int]] = []
    current: list[int] = []
    chunk_end = sites[0].pos0 + int(chunk_size_bp)
    for idx, site in enumerate(sites):
        while current and site.pos0 >= chunk_end:
            chunks.append(current)
            current = []
            chunk_end += int(chunk_size_bp)
        while site.pos0 >= chunk_end:
            chunk_end += int(chunk_size_bp)
        current.append(idx)
    if current:
        chunks.append(current)
    return chunks


def _extract_chunk(payload) -> ExtractionResult:
    bam_path, all_sites, site_indices, min_baseq, min_mapq, epsilon_floor, epsilon_ceiling, mapq_cap, default_baseq, scan_mode = payload
    chunk_sites = [all_sites[idx] for idx in site_indices]
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        extractor = extract_observations_read_scan if scan_mode == "read" else extract_observations_pileup
        read_ids, rows, n_ref, n_alt, n_other = extractor(
            bam,
            sites=chunk_sites,
            min_baseq=min_baseq,
            min_mapq=min_mapq,
            epsilon_floor=epsilon_floor,
            epsilon_ceiling=epsilon_ceiling,
            mapq_cap=mapq_cap,
            default_baseq=default_baseq,
        )
    return ExtractionResult(list(site_indices), ordered_read_names(read_ids), rows, n_ref, n_alt, n_other)


def merge_chunk_results(results: list[ExtractionResult]) -> tuple[dict[str, int], list[ObservationRow], int, int, int]:
    read_ids: dict[str, int] = {}
    best_row_by_key: dict[tuple[int, int], ObservationRow] = {}
    n_other = 0
    for result in results:
        n_other += result.n_other
        for row in result.rows:
            local_read_idx, local_snp_idx, allele, baseq, mapq, eps, weight = row
            read_name = result.read_names[int(local_read_idx)]
            read_idx = read_ids.setdefault(read_name, len(read_ids))
            snp_idx = int(result.site_indices[int(local_snp_idx)])
            merged = (read_idx, snp_idx, allele, baseq, mapq, eps, weight)
            key = (read_idx, snp_idx)
            previous = best_row_by_key.get(key)
            if previous is None or abs(float(weight)) > abs(float(previous[6])):
                best_row_by_key[key] = merged
    rows = sorted(best_row_by_key.values(), key=lambda item: (item[0], item[1]))
    return read_ids, rows, sum(1 for row in rows if row[2] == 1), sum(1 for row in rows if row[2] == -1), n_other


def extract_observations_chunked(
    bam_path: str,
    sites: list[SnpSite],
    min_baseq: int,
    min_mapq: int,
    epsilon_floor: float,
    epsilon_ceiling: float,
    mapq_cap: int,
    default_baseq: int | None,
    scan_mode: str,
    threads: int,
    chunk_size_bp: int,
) -> tuple[dict[str, int], list[ObservationRow], int, int, int, int]:
    chunks = chunk_site_indices(sites, chunk_size_bp)
    if threads <= 1 or len(chunks) == 1:
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            extractor = extract_observations_read_scan if scan_mode == "read" else extract_observations_pileup
            read_ids, rows, n_ref, n_alt, n_other = extractor(
                bam,
                sites=sites,
                min_baseq=min_baseq,
                min_mapq=min_mapq,
                epsilon_floor=epsilon_floor,
                epsilon_ceiling=epsilon_ceiling,
                mapq_cap=mapq_cap,
                default_baseq=default_baseq,
            )
        return read_ids, rows, n_ref, n_alt, n_other, 1

    payloads = [
        (bam_path, sites, chunk, min_baseq, min_mapq, epsilon_floor, epsilon_ceiling, mapq_cap, default_baseq, scan_mode)
        for chunk in chunks
    ]
    with ProcessPoolExecutor(max_workers=int(threads)) as pool:
        results = list(pool.map(_extract_chunk, payloads))
    read_ids, rows, n_ref, n_alt, n_other = merge_chunk_results(results)
    return read_ids, rows, n_ref, n_alt, n_other, len(chunks)


__all__ = [
    "ExtractionResult",
    "extract_observations_chunked",
    "extract_observations_pileup",
    "extract_observations_read_scan",
]
