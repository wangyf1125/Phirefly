"""Contact filters and the single-batch joint long-read/Hi-C workflow."""
import csv
import json
import shutil

import numpy as np
import pysam
import pytest
from scipy import sparse

from phirefly.hic.components import ComponentSite
from phirefly.hic.contact import alignment_anchor
from phirefly.hic.joint import HiCConfig, build_hic_coupling
from phirefly.pipeline import build_parser, run_pipeline


SITES = (10, 20, 40, 50, 80)
HEADER = {'HD': {'VN': '1.6', 'SO': 'coordinate'},
          'SQ': [{'SN': 'chr1', 'LN': 100}]}


def alignment(name, start, base, flag=0, quality=30):
    aln = pysam.AlignedSegment(pysam.AlignmentHeader.from_dict(HEADER))
    aln.query_name = name
    aln.query_sequence = base*11
    aln.reference_id = 0
    aln.reference_start = start
    aln.flag = flag
    aln.mapping_quality = 60
    aln.cigarstring = '11M'
    aln.query_qualities = [quality]*11 if quality is not None else None
    return aln


def contacts(path, n=20, quality=30):
    with pysam.AlignmentFile(path, 'wb', header=HEADER) as out:
        for start, base, mate in ((9, 'A', 1), (39, 'C', 2)):
            for i in range(n):
                out.write(alignment(f'pair{i}/{mate}', start, base, quality=quality))
    pysam.index(str(path))
    return path


def test_anchor_requires_quality_and_primary_nonduplicate_alignment():
    config = HiCConfig()
    sites = {'chr1': {9: [ComponentSite('A', 'C', 1, '0', 'chr1:10:A>C')]}}
    positions = {'chr1': [9]}
    assert alignment_anchor(alignment('x', 9, 'A'), sites, positions, config)['sign'] == 1
    for flag in (256, 512, 1024, 2048):
        assert alignment_anchor(alignment('x', 9, 'A', flag), sites, positions, config) is None
    for quality in (None, 5):
        assert alignment_anchor(alignment('x', 9, 'A', quality=quality), sites, positions, config) is None


def test_joint_contacts_require_external_component_support(tmp_path):
    bam = contacts(tmp_path/'contacts.bam')
    matrix = sparse.csr_matrix([[1, 1, 0, 0, 0], [0, 0, 1, 1, 0], [0, 0, 0, 0, 1]])
    snps = [f'chr1:{p}:A>C' for p in SITES]
    ids = np.array([0, 0, 1, 1, 2])
    evidence = build_hic_coupling(bam, 'chr1:1-100', matrix, snps, np.ones(5), ids, HiCConfig())
    assert evidence.coupling.shape == (3, 3)
    assert evidence.stats['hic_external_read_names'] == 20
    assert evidence.stats['hic_retained_couplings'] == 1
    assert evidence.coupling[0, 1] == -4
    assert evidence.coupling[2].nnz == 0
    rejected = build_hic_coupling(bam, 'chr1:1-100', matrix, snps, np.ones(5), ids,
                                  HiCConfig(min_molecules=21))
    assert rejected.coupling.nnz == 0
    # Internal contact links do not require evidence for a new long-range join.
    linked = sparse.csr_matrix([[1, 1, 1, 1, 0], [0, 0, 0, 0, 1]])
    local = build_hic_coupling(bam, 'chr1:1-100', linked, snps, np.ones(5), ids,
                              HiCConfig(min_molecules=100))
    assert local.coupling[0, 1] == -4


@pytest.mark.qaia
def test_joint_bam_to_vcf_solves_once_and_preserves_isolated_phase_set(tmp_path):
    pytest.importorskip('mindquantum')
    if not shutil.which('bgzip') or not shutil.which('tabix'):
        pytest.skip('htslib commands are required')
    bam = tmp_path/'reads.bam'
    with pysam.AlignmentFile(bam, 'wb', header=HEADER) as out:
        for start, name in ((9, 'left'), (39, 'right')):
            for i in range(2):
                out.write(alignment(name+str(i), start, 'A'))
        out.write(alignment('isolate', 79, 'A'))
    pysam.index(str(bam))
    vcf = tmp_path/'input.vcf'
    vcf.write_text('##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100>\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n'+
        ''.join(f'chr1\t{p}\t.\tA\tC\t.\tPASS\t.\tGT\t0/1\n' for p in SITES))
    pysam.tabix_compress(str(vcf), str(vcf)+'.gz')
    pysam.tabix_index(str(vcf)+'.gz', preset='vcf')
    hic = contacts(tmp_path/'contacts.bam')
    out = tmp_path/'result'
    args = build_parser().parse_args(['--bam', str(bam), '--vcf', str(vcf)+'.gz',
        '--hic-bam', str(hic), '--region', 'chr1:1-100', '--out-dir', str(out),
        '--backend', 'cpu-float32', '--batch-size', '20', '--n-iter', '200', '--debug-output'])
    run_pipeline(args)
    with pysam.VariantFile(next(out.rglob('phased.component_ps.vcf.gz'))) as result:
        calls = {r.pos: r.samples['SAMPLE'] for r in result.fetch('chr1', 0, 100)}
    assert [calls[p]['PS'] for p in SITES] == [10, 10, 10, 10, 80]
    assert calls[10]['GT'] == calls[20]['GT']
    assert calls[40]['GT'] == calls[50]['GT']
    assert calls[10]['GT'] != calls[40]['GT']
    with next(out.rglob('config_metrics.tsv')).open() as handle:
        summary = next(csv.DictReader(handle, delimiter='\t'))
    assert summary['hic_retained_couplings'] == '1'
    assert summary['hic_min_molecules'] == '20'
    assert 'qaia_calls\t1\n' in next(out.rglob('qaia_summary.tsv')).read_text()
    manifest = json.loads((out/'run_manifest.json').read_text())
    assert manifest['pipeline_revision'] == 'single_cfc_hic_joint_v1'
    assert manifest['inputs']['hic_bam']['path'] == str(hic)
    assert len(list(out.rglob('hic_coupling.npz'))) == 1
