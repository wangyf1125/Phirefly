"""Low-level alignment helpers."""

from __future__ import annotations

from collections.abc import Iterable


def aligned_snp_query_positions(aln, candidate_positions: list[int]) -> Iterable[tuple[int, int]]:
    """Yield (ref_pos0, query_pos0) for candidate SNP positions covered by an alignment."""
    if not candidate_positions or aln.cigartuples is None:
        return

    qpos = 0
    rpos = aln.reference_start
    cand_i = 0
    n_candidates = len(candidate_positions)

    for op, length in aln.cigartuples:
        if op in (0, 7, 8):  # M, =, X consume query and reference.
            block_start = rpos
            block_end = rpos + length
            while cand_i < n_candidates and candidate_positions[cand_i] < block_start:
                cand_i += 1
            j = cand_i
            while j < n_candidates and candidate_positions[j] < block_end:
                ref_pos = candidate_positions[j]
                yield ref_pos, qpos + (ref_pos - block_start)
                j += 1
            cand_i = j
            qpos += length
            rpos += length
        elif op in (1, 4):  # I, S consume query only.
            qpos += length
        elif op in (2, 3):  # D, N consume reference only.
            rpos += length
        elif op in (5, 6):  # H, P consume neither.
            continue


__all__ = ["aligned_snp_query_positions"]
