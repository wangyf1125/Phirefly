from phirefly.metrics import PhaseCall, haplotype_n50_kb, phased_block_n50_kb, phase_error_metrics


def calls(positions, ps, phase=0, chrom='chr1'):
    return [PhaseCall(chrom, p, 'A', 'G', phase, ps) for p in positions]


def test_interleaved_phase_sets_do_not_inflate_pb_n50():
    data = calls([1, 101, 901, 1001], 'A') + calls([401, 501], 'B') + calls([601], 'singleton')
    assert haplotype_n50_kb(data) == 1.001
    assert phased_block_n50_kb(data) == 0.100


def test_chromosomes_and_singletons_remain_separate():
    data = calls([1, 101], 'A') + calls([20, 30], 'A', chrom='chr2')
    assert phased_block_n50_kb(data) == .100
    assert phased_block_n50_kb(calls([100], 'A')) == 0.
    assert phased_block_n50_kb([]) == 0.


def test_hamming_gauge_uses_joint_truth_and_prediction_blocks():
    truth = calls([10, 20], 'T1') + calls([30, 40], 'T2', phase=1)
    prediction = calls([10, 20, 30, 40], 'joined')
    metrics = phase_error_metrics(truth, prediction)
    assert metrics['hamming'] == metrics['switches'] == 0
    assert metrics['covered_variants'] == 4
    assert metrics['assessed_pairs'] == 2
