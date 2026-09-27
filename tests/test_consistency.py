"""Regression tests for read-only phaselet relaxation and candidate selection."""
import numpy as np
from scipy import sparse

from phirefly.phaselet.consistency import (
    project_matrix, read_consistency_margin, relax_phaselets,
    select_component_candidates,
)


def test_independent_component_selection_combines_existing_candidates():
    # Candidate 0 wins the first component; candidate 1 wins the second.
    projected = sparse.csr_matrix([[3., 2., 0., 0.], [0., 0., 1., 3.]])
    samples = np.array([[1, 1], [1, -1], [1, 1], [-1, 1]], dtype=np.int8)
    coupling = sparse.csr_matrix([[0, 1, 0, 0], [1, 0, 0, 0],
                                  [0, 0, 0, 1], [0, 0, 1, 0]])
    selected, labels, mode = select_component_candidates(projected, samples, coupling, 0)
    assert mode == 'component_profiled_F'
    assert np.array_equal(selected, [1, 1, 1, 1])
    assert np.abs(projected @ selected).sum() == 9
    assert np.abs(projected @ samples).sum(axis=0).max() == 7
    assert labels[0] != labels[2]


def test_cross_component_validation_does_not_assume_separability():
    projected = sparse.csr_matrix([[1., 1.]])
    samples = np.array([[1, 1], [-1, 1]], dtype=np.int8)
    selected, _, mode = select_component_candidates(
        projected, samples, sparse.csr_matrix((2, 2)), 0)
    assert mode == 'global_nonseparable_validation'
    assert np.array_equal(selected, [1, 1])


def test_leave_one_site_out_does_not_vote_for_itself():
    matrix = sparse.csr_matrix([[1., 1., -1.], [0., 0., -100.]])
    tau = np.ones(3, dtype=np.int8)
    # Single-site read abstains even though its weight is very large.
    assert np.array_equal(read_consistency_margin(matrix, tau), [0, 0, -1])


def test_relaxation_preserves_good_merge_and_releases_suspect_site():
    matrix = sparse.csr_matrix([[1., 1., -1.]])
    tau = np.ones(3, dtype=np.int8)
    snps = ['chr1:10:A:C', 'chr1:20:A:C', 'chr1:30:A:C']
    edges = [dict(snp_i=0, snp_j=1, status='hard_retained'),
             dict(snp_i=1, snp_j=2, status='hard_retained')]
    ids, phaselets, margins = relax_phaselets(matrix, tau, snps, edges)
    assert ids[0] == ids[1] != ids[2]
    assert len(phaselets) == 2
    assert np.array_equal(tau, [1, 1, 1])  # Never edit the phase directly.
    assert margins[2] < 0
    assert np.array_equal(project_matrix(matrix, ids, tau).toarray(), [[2., -1.]])


def test_component_selection_matches_exhaustive_component_choices():
    rng = np.random.default_rng(3)
    matrix = sparse.block_diag([np.array([[1., 3.], [2., -1.]]),
                                np.array([[4., 2.]])], format='csr')
    coupling = sparse.block_diag([np.ones((2, 2)), np.ones((2, 2))], format='csr')
    samples = rng.choice([-1, 1], size=(4, 9))
    z, _, _ = select_component_candidates(matrix, samples, coupling, 0)
    best = max(np.abs(matrix @ np.r_[samples[:2, a], samples[2:, b]]).sum()
               for a in range(9) for b in range(9))
    assert np.abs(matrix @ z).sum() == best
