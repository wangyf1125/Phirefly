"""Read-consistency phaselets and separable candidate selection."""
from dataclasses import dataclass
import time

import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import connected_components

from ..hyperread.build import build_hyperread_matrix
from ..qaia import build_ising_coupling, solve_qaia, QAIAResult
from .build import DSU, build_phaselet_table
from .soft_bridge import build_soft_phaselet_coupling


def project_matrix(matrix, phaselet_ids, tau):
    transform = sparse.csr_matrix((tau.astype(float),
        (np.arange(len(tau)), phaselet_ids)), shape=(len(tau), int(phaselet_ids.max())+1))
    projected = (matrix @ transform).tocsr()
    projected.sum_duplicates()
    projected.eliminate_zeros()
    projected.sort_indices()
    return projected


def read_consistency_margin(matrix, tau):
    """Each observation is predicted from the other sites on its read; ties abstain."""
    coo = matrix.tocoo()
    fields = np.asarray(matrix @ tau).ravel()
    labels = np.sign(fields[coo.row] - coo.data * tau[coo.col])
    return tau * np.bincount(coo.col, weights=coo.data * labels, minlength=len(tau))


def relax_phaselets(matrix, tau, snp_ids, edge_rows):
    margin = read_consistency_margin(matrix, tau)
    suspect = margin < 0
    forest = DSU(len(tau))
    for edge in edge_rows:
        i, j = int(edge['snp_i']), int(edge['snp_j'])
        if edge['status'] in {'hard_retained', 'hard_cycle'} and not (suspect[i] or suspect[j]):
            forest.union(i, j)
    positions = np.array([int(s.split(':')[1]) for s in snp_ids])
    ids, phaselets = build_phaselet_table(forest, snp_ids, positions)
    return ids, phaselets, margin


def select_component_candidates(projected, samples, coupling, n_hyper, row_weights=None):
    """Select existing candidates, never change a spin within a candidate component.

    A held-out read can couple training components. In that case use global
    selection rather than incorrectly assuming the validation objective separates.
    """
    count, labels = connected_components(coupling, directed=False)
    phase_labels = labels[n_hyper:]
    fields = np.abs(projected @ samples)
    if row_weights is not None:
        fields *= np.asarray(row_weights)[:, None]
    row_labels = np.full(projected.shape[0], -1, dtype=np.int64)
    for r in range(projected.shape[0]):
        ids = projected.indices[projected.indptr[r]:projected.indptr[r+1]]
        if len(ids):
            components = phase_labels[ids]
            if np.any(components != components[0]):
                best = int(fields.sum(axis=0).argmax())
                return samples[:, best], phase_labels, 'global_nonseparable_validation'
            row_labels[r] = components[0]
    scores = np.zeros((count, samples.shape[1]), dtype=float)
    active = row_labels >= 0
    np.add.at(scores, row_labels[active], fields[active])
    choices = scores.argmax(axis=1)
    selected = samples[np.arange(len(phase_labels)), choices[phase_labels]]
    observed = np.abs(projected @ selected)
    if row_weights is not None:
        observed *= row_weights
    if not np.isclose(observed.sum(), scores.max(axis=1).sum()):
        raise AssertionError('Component objective does not separate')
    return selected, phase_labels, 'component_profiled_F'


@dataclass
class PhaseletSolution:
    hyper: sparse.csr_matrix
    hyper_rows: list
    hyper_edges: list
    hyper_stats: dict
    soft: sparse.csr_matrix
    soft_rows: list
    soft_stats: dict
    qaia: QAIAResult
    delta: np.ndarray
    z: np.ndarray
    ps: np.ndarray
    selection_mode: str
    timings: dict


def solve_stage(entries, matrix, snp_ids, tau, ids, phaselets, soft_edges, args,
                row_weights=None, trace_output=None, stage='relaxed'):
    started = time.monotonic()
    hyper, rows, edges, hs = build_hyperread_matrix(entries, tau, ids, len(phaselets),
        hyperread_norm=args.hyperread_norm, hyperread_cap=args.hyperread_cap,
        adaptive_dom_threshold=args.adaptive_dom_threshold)
    soft, soft_rows, soft_stats = build_soft_phaselet_coupling(soft_edges, ids, tau,
        len(phaselets), args.soft_edge_scale, args.soft_edge_cap)
    coupling = build_ising_coupling(hyper, soft)
    if coupling.shape[0] > args.max_qaia_nodes:
        raise ValueError(f'Compressed graph has {coupling.shape[0]} nodes; limit is {args.max_qaia_nodes}')
    build_seconds = time.monotonic()-started
    result = solve_qaia(hyper, soft, args.algorithm, args.backend, args.batch_size,
        args.n_iter, .1, args.seed, args.seed_scale, args.seed_noise,
        trace_every=args.trace_every, trace_output=trace_output,
        trace_metadata={'region': str(args.region_key), 'solver': args.algorithm,
                        'seed': args.seed, 'batch_size': args.batch_size,
                        'nodes': coupling.shape[0]})
    started = time.monotonic()
    projected = project_matrix(matrix, ids, tau)
    samples = result.phaselet_spins
    z, labels, mode = select_component_candidates(projected, samples, coupling, hyper.shape[0], row_weights)
    positions = np.array([int(s.split(':')[1]) for s in snp_ids])
    minima = np.full(labels.max()+1, np.iinfo(np.int64).max, dtype=np.int64)
    np.minimum.at(minima, labels[ids], positions)
    return PhaseletSolution(hyper, rows, edges, hs, soft, soft_rows, soft_stats, result,
        tau*z[ids], z, minima[labels[ids]], mode,
        {f't_{stage}_build_graph': build_seconds,
         f't_{stage}_candidate_scoring': time.monotonic()-started})
