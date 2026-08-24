"""Contact-BAM anchor extraction for Hi-C/Pore-C component edges."""

from __future__ import annotations

import bisect
from collections import Counter

from ..core.alignment import aligned_snp_query_positions


def normalize_query_name(query_name: str) -> str:
    if query_name.endswith("/1") or query_name.endswith("/2"):
        return query_name[:-2]
    return query_name


def alignment_anchor(aln, site_map, site_positions, args):
    """Return the strongest component/sign anchor supported by one alignment."""

    if aln.is_unmapped:
        return None
    if aln.is_secondary and not args.keep_secondary:
        return None
    if aln.is_supplementary and not args.keep_supplementary:
        return None
    if aln.mapping_quality < args.min_mapq or aln.query_sequence is None:
        return None
    if aln.reference_start is None or aln.reference_end is None:
        return None
    chrom = aln.reference_name
    if chrom not in site_map:
        return None

    positions = site_positions.get(chrom, [])
    left = bisect.bisect_left(positions, aln.reference_start)
    right = bisect.bisect_left(positions, aln.reference_end)
    if left >= right:
        return None

    seq = aln.query_sequence.upper()
    component_score: Counter[str] = Counter()
    n_sites = 0
    for ref_pos, qpos in aligned_snp_query_positions(aln, positions[left:right]):
        base = seq[qpos]
        for site in site_map[chrom][ref_pos]:
            if base == site.ref:
                allele = 1
            elif base == site.alt:
                allele = -1
            else:
                continue
            component_score[site.component] += allele * site.delta
            n_sites += 1

    if n_sites < args.min_sites or not component_score:
        return None
    component, score_sum = max(component_score.items(), key=lambda item: abs(item[1]))
    if abs(score_sum) < args.min_local_score:
        return None
    mapq_scale = min(float(aln.mapping_quality), float(args.mapq_cap)) / float(args.mapq_cap)
    return {
        "query_name": normalize_query_name(aln.query_name),
        "component_id": component,
        "sign": 1 if score_sum > 0 else -1,
        "score": abs(float(score_sum)) * max(mapq_scale, 0.05),
        "n_sites": int(n_sites),
        "chrom": chrom,
        "start": int(aln.reference_start) + 1,
        "end": int(aln.reference_end),
        "mapq": int(aln.mapping_quality),
    }


__all__ = ["alignment_anchor", "normalize_query_name"]
