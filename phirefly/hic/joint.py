"""Add quality-filtered contact evidence to the phaselet orientation problem."""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import connected_components

from ..core.genomics import apply_region_alias, parse_chrom_alias, parse_region, parse_snp_id
from .components import ComponentSite
from .component_edges import collect_anchors, component_edges_from_anchors, remap_site_chroms


@dataclass(frozen=True)
class HiCConfig:
    min_mapq: int = 20
    min_baseq: int = 13
    min_contacts: int = 1
    min_molecules: int = 20
    scale: float = 0.1
    cap: float = 10.0
    threads: int = 1
    keep_secondary: bool = False
    keep_supplementary: bool = False
    min_sites: int = 1
    min_local_score: float = 1.0
    mapq_cap: int = 60

    @classmethod
    def from_args(cls, args):
        return cls(min_mapq=args.hic_min_mapq, min_baseq=args.hic_min_baseq,
            min_contacts=args.hic_min_contacts, min_molecules=args.hic_min_molecules,
            scale=args.hic_scale, cap=args.hic_cap, threads=getattr(args, 'threads', 1))


def add_hic_arguments(parser):
    parser.add_argument('--hic-bam', type=Path, help='Optional same-sample indexed Hi-C BAM; solved jointly in the single CFC batch.')
    parser.add_argument('--hic-chrom-alias', action='append', default=[], help='Explicit chromosome alias, e.g. chr1=NC_000001.11.')
    config = HiCConfig()
    for name, kind, help_text in [
        ('min-mapq', int, 'Minimum Hi-C mapping quality.'),
        ('min-baseq', int, 'Minimum Hi-C SNP base quality; missing qualities are excluded.'),
        ('min-contacts', int, 'Minimum distinct read names supporting a phaselet pair.'),
        ('min-molecules', int, 'Minimum external contact read names per original long-read component.'),
        ('scale', float, 'Scale signed contact support relative to original read evidence.'),
        ('cap', float, 'Absolute contact-coupling cap.')]:
        parser.add_argument('--hic-'+name, type=kind,
            default=getattr(config, name.replace('-', '_')), help=help_text)


@dataclass(frozen=True)
class HiCEvidence:
    coupling: sparse.csr_matrix
    edges: list[dict]
    stats: dict


def build_hic_coupling(bam, region, matrix, snp_ids, tau, ids, config, aliases=()):
    """Keep local phaselets mutable while requiring supported long-range joins."""
    if (min(config.min_contacts, config.min_molecules, config.threads) < 1 or
            min(config.min_mapq, config.min_baseq) < 0 or
            not np.isfinite([config.scale, config.cap]).all() or min(config.scale, config.cap) <= 0):
        raise ValueError('Hi-C counts, scale and cap must be positive; qualities must be nonnegative')
    n_phaselets = int(ids.max())+1
    read_graph = sparse.bmat([[None, abs(matrix)], [abs(matrix.T), None]], format='csr')
    _, labels = connected_components(read_graph, directed=False)
    labels = labels[matrix.shape[0]:]
    groups = np.full(n_phaselets, -1, dtype=int)
    sites = {}
    for i, name in enumerate(snp_ids):
        k = int(ids[i])
        if groups[k] not in (-1, labels[i]):
            raise ValueError('Phaselet crosses original read components')
        groups[k] = labels[i]
        chrom, pos, ref, alt = parse_snp_id(name)
        sites.setdefault(chrom, {}).setdefault(pos-1, []).append(ComponentSite(ref, alt, int(tau[i]), str(k), name))
    positions = {chrom: sorted(values) for chrom, values in sites.items()}
    alias_map = parse_chrom_alias(list(aliases))
    sites, positions = remap_site_chroms(sites, positions, alias_map)
    interval = apply_region_alias(parse_region(region, zero_based_start=True), alias_map)
    anchors = collect_anchors(bam, interval, sites, positions, config)
    external = Counter()
    bridging = 0
    for rows in anchors.values():
        members = {int(groups[int(row['component_id'])]) for row in rows}
        if len(members) > 1:
            external.update(members)
            bridging += 1
    all_edges = component_edges_from_anchors(anchors, 0.1)
    retained, row, col, values = [], [], [], []
    for edge in all_edges:
        i, j = int(edge['component_i']), int(edge['component_j'])
        score = float(edge['score'])
        if int(edge['contacts']) < config.min_contacts or score == 0:
            continue
        if groups[i] != groups[j] and min(external[groups[i]], external[groups[j]]) < config.min_molecules:
            continue
        weight = float(np.clip(config.scale*score, -config.cap, config.cap))
        row.extend([i, j])
        col.extend([j, i])
        values.extend([weight, weight])
        retained.append(dict(edge, normalized_coupling=weight))
    coupling = sparse.coo_matrix((values, (row, col)), shape=(n_phaselets, n_phaselets)).tocsr()
    stats = dict(hic_bam=str(bam), hic_anchor_read_names=len(anchors),
        hic_external_read_names=bridging, hic_candidate_edges=len(all_edges),
        hic_retained_couplings=len(retained), hic_scale=config.scale, hic_cap=config.cap,
        hic_min_contacts=config.min_contacts, hic_min_molecules=config.min_molecules,
        hic_min_mapq=config.min_mapq, hic_min_baseq=config.min_baseq)
    return HiCEvidence(coupling, retained, stats)
