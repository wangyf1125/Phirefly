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


def select_component_candidates(projected, samples, coupling, n_hyper, row_weights=None, pair_coupling=None):
    """Select existing candidates, never change a spin within a candidate component.

    A held-out read can couple training components. In that case use global
    selection rather than incorrectly assuming the validation objective separates.
    """
    graph = coupling
    if pair_coupling is not None and pair_coupling.nnz:
        graph = abs(coupling).tolil()
        graph[n_hyper:, n_hyper:] = graph[n_hyper:, n_hyper:]+abs(pair_coupling)
        graph = graph.tocsr()
    count, labels = connected_components(graph, directed=False)
    phase_labels = labels[n_hyper:]
    fields = np.abs(projected @ samples)
    if row_weights is not None:
        fields *= np.asarray(row_weights)[:, None]
    pair_fields = (.5*samples*(pair_coupling @ samples)) if pair_coupling is not None else None
    row_labels = np.full(projected.shape[0], -1, dtype=np.int64)
    for r in range(projected.shape[0]):
        ids = projected.indices[projected.indptr[r]:projected.indptr[r+1]]
        if len(ids):
            components = phase_labels[ids]
            if np.any(components != components[0]):
                scores = fields.sum(axis=0)
                if pair_fields is not None:
                    scores += pair_fields.sum(axis=0)
                best = int(scores.argmax())
                return samples[:, best], phase_labels, 'global_nonseparable_validation'
            row_labels[r] = components[0]
    scores = np.zeros((count, samples.shape[1]), dtype=float)
    active = row_labels >= 0
    np.add.at(scores, row_labels[active], fields[active])
    if pair_fields is not None:
        np.add.at(scores, phase_labels, pair_fields)
    choices = scores.argmax(axis=1)
    selected = samples[np.arange(len(phase_labels)), choices[phase_labels]]
    observed = np.abs(projected @ selected)
    if row_weights is not None:
        observed *= row_weights
    objective = observed.sum()
    if pair_coupling is not None:
        objective += .5*selected @ (pair_coupling @ selected)
    if not np.isclose(objective, scores.max(axis=1).sum()):
        raise AssertionError('Component objective does not separate')
    return selected, phase_labels, 'component_joint_evidence' if pair_coupling is not None else 'component_profiled_F'


def select_energy_candidates(coupling, eta, z, pair_coupling=None):
    """Choose measured joint spins by energy; ties retain the first candidate."""
    spins = np.vstack([eta, z])
    if spins.ndim != 2 or not spins.shape[1] or not np.isin(spins, [-1, 1]).all():
        raise ValueError('Expected a nonempty batch of binary spins')
    if coupling.shape != (len(spins), len(spins)):
        raise ValueError('Coupling and sample dimensions differ')
    if (coupling-coupling.T).nnz or np.any(coupling.diagonal()):
        raise ValueError('Expected symmetric couplings without self-interactions')
    graph = coupling.copy()
    graph.eliminate_zeros()
    if pair_coupling is not None and pair_coupling.nnz:
        # Keep the existing Hi-C phase-set boundaries when terms cancel.
        graph = abs(graph).tolil()
        graph[len(eta):, len(eta):] += abs(pair_coupling)
        graph = graph.tocsr()
    count, labels = connected_components(graph, directed=False)
    energies = np.zeros((count, spins.shape[1]))
    np.add.at(energies, labels, -.5*spins*(coupling @ spins))
    choices = energies.argmin(axis=1)
    selected = spins[np.arange(len(spins)), choices[labels]]
    energy = float(energies[np.arange(count), choices].sum())
    return selected[len(eta):], labels[len(eta):], energy


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
    contact: sparse.csr_matrix | None = None
    selected_energy: float | None = None


def solve_stage(entries, matrix, snp_ids, tau, ids, phaselets, soft_edges, args,
                row_weights=None, trace_output=None, stage='relaxed', contact_coupling=None):
    started = time.monotonic()
    hyper, rows, edges, hs = build_hyperread_matrix(entries, tau, ids, len(phaselets),
        hyperread_norm=args.hyperread_norm, hyperread_cap=args.hyperread_cap,
        adaptive_dom_threshold=args.adaptive_dom_threshold)
    soft, soft_rows, soft_stats = build_soft_phaselet_coupling(soft_edges, ids, tau,
        len(phaselets), args.soft_edge_scale, args.soft_edge_cap)
    extra = soft+contact_coupling if contact_coupling is not None else soft
    coupling = build_ising_coupling(hyper, extra)
    if coupling.shape[0] > args.max_qaia_nodes:
        raise ValueError(f'Compressed graph has {coupling.shape[0]} nodes; limit is {args.max_qaia_nodes}')
    build_seconds = time.monotonic()-started
    result = solve_qaia(hyper, extra, args.algorithm, args.backend, args.batch_size,
        args.n_iter, .1, args.seed, args.seed_scale, args.seed_noise,
        trace_every=args.trace_every, trace_output=trace_output,
        trace_metadata={'region': str(args.region_key), 'solver': args.algorithm,
                        'seed': args.seed, 'batch_size': args.batch_size,
                        'nodes': coupling.shape[0]})
    started = time.monotonic()
    z, labels, energy = select_energy_candidates(coupling, result.hyperread_spins,
        result.phaselet_spins, pair_coupling=contact_coupling)
    positions = np.array([int(s.split(':')[1]) for s in snp_ids])
    minima = np.full(labels.max()+1, np.iinfo(np.int64).max, dtype=np.int64)
    np.minimum.at(minima, labels[ids], positions)
    return PhaseletSolution(hyper, rows, edges, hs, soft, soft_rows, soft_stats, result,
        tau*z[ids], z, minima[labels[ids]], 'component_minimum_energy',
        {f't_{stage}_build_graph': build_seconds,
         f't_{stage}_candidate_scoring': time.monotonic()-started}, contact_coupling, energy)
