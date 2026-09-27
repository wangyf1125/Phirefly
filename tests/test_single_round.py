"""The release runner must solve only the read-consistency partition."""
import numpy as np
import pytest

from phirefly.phaselet import runner, consistency
from phirefly.core.io import read_key_value
from phirefly.pipeline import build_parser, phaselet_args
from phirefly.qaia import QAIAResult, build_ising_coupling


def test_single_call_after_release_and_old_cache_rejected(tmp_path, monkeypatch):
    snps = ['chr1:10:A>C', 'chr1:20:A>C', 'chr1:30:A>C']
    obs = tmp_path / 'obs.npz'
    np.savez(obs, read_ids=['r1'], snp_ids=snps, read_idx=[0, 0, 0],
             snp_idx=[0, 1, 2], b_ki=[1, 1, -1])
    tau = tmp_path / 'tau.tsv'
    tau.write_text('snp_id\tdelta\n' + ''.join(f'{s}\t1\n' for s in snps))
    vcf = tmp_path / 'input.vcf'
    vcf.touch()
    base = build_parser().parse_args(['--bam', 'unused', '--vcf', str(vcf),
        '--region', 'chr1:1-100', '--out-dir', str(tmp_path), '--batch-size', '4'])
    args = phaselet_args(base, obs, None, tau, tmp_path / 'solve', 'SAMPLE')
    events = []
    relax = consistency.relax_phaselets

    def release(*values):
        events.append('release')
        return relax(*values)

    def qaia(hyper, soft, algorithm, backend, batch, *values, **kw):
        events.append('CFC')
        return QAIAResult(np.ones((hyper.shape[1], batch), dtype=np.int8),
            np.ones((hyper.shape[0], batch), dtype=np.int8), np.zeros(batch),
            build_ising_coupling(hyper, soft).nnz, 0.125)

    monkeypatch.setattr(runner, 'relax_phaselets', release)
    monkeypatch.setattr(consistency, 'solve_qaia', qaia)
    monkeypatch.setattr(runner, 'export_component_vcf', lambda *a: None)
    runner.run_phaselet_qaia(args)
    assert events == ['release', 'CFC']
    out = args.outroot / 'phirefly/alg_CFC/phaselet_risk0p45_soft0p75_s0p25'
    summary = read_key_value(out / 'qaia_summary.tsv')
    assert summary['pipeline_revision'] == 'single_cfc_v1'
    assert summary['qaia_calls'] == '1'
    assert summary['candidate_count'] == '4'
    assert float(summary['qaia_runtime_s']) == 0.125
    assert 'qaia_initial_runtime_s' not in summary
    cache = out / 'phirefly/phaselet_hyperread_qaia.summary.tsv'
    cache.write_text(cache.read_text().replace('single_cfc_v1', 'read_consistency_v1'))
    with pytest.raises(ValueError, match='not single-round'):
        runner.run_phaselet_qaia(args)
