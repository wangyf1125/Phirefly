"""Exercise contact BAM -> retained edges -> orientation -> VCF via installed CLI."""
import csv
import subprocess
import sys

import pysam


def test_hic_cli_preserves_unconnected_phase_sets(tmp_path):
    sites = (10, 20, 40, 50, 80)
    phases = tmp_path / 'phases.tsv'
    phases.write_text('snp_id\tdelta\tps\n' + ''.join(
        f'chr1:{p}:A>C\t1\t{10 if p < 30 else 40 if p < 60 else 80}\n' for p in sites))
    vcf = tmp_path / 'input.vcf'
    vcf.write_text('##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100>\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n' +
        ''.join(f'chr1\t{p}\t.\tA\tC\t.\tPASS\t.\tGT\t0/1\n' for p in sites))
    bam = tmp_path / 'contacts.bam'
    with pysam.AlignmentFile(bam, 'wb', header={
            'HD': {'VN': '1.6', 'SO': 'coordinate'}, 'SQ': [{'SN': 'chr1', 'LN': 100}]}) as out:
        for start, base, mate in ((9, 'A', 1), (39, 'C', 2)):
            for i in range(20):
                aln = pysam.AlignedSegment()
                aln.query_name = f'pair{i}/{mate}'
                aln.query_sequence = base * 11
                aln.reference_id = 0
                aln.reference_start = start
                aln.flag = 0
                aln.mapping_quality = 60
                aln.cigarstring = '11M'
                aln.query_qualities = pysam.qualitystring_to_array('I' * 11)
                out.write(aln)
    pysam.index(str(bam))
    prefix = tmp_path / 'contact'

    def cli(*args):
        subprocess.run([sys.executable, '-m', 'phirefly.cli', *map(str, args)],
                       check=True, capture_output=True, text=True)

    # Explicit PS records avoid reconstructing components from observations.
    obs = tmp_path / 'unused-observations.npz'
    cli('hic-edges', '--observations', obs, '--snp-phases', phases, '--contact-bam', bam,
        '--region', 'chr1:1-100', '--out-prefix', prefix)
    edges = prefix.with_suffix('.block_edges.tsv')
    with edges.open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    assert len(rows) == 1
    assert rows[0]['contacts'] == '20'
    assert float(rows[0]['score']) < 0
    cli('hic-orient', '--edges', edges, '--out-prefix', tmp_path / 'oriented')
    output = tmp_path / 'hic.vcf.gz'
    cli('hic-apply', '--input-vcf', vcf, '--sample', 'SAMPLE', '--observations', obs,
        '--snp-phases', phases, '--component-orientations',
        tmp_path / 'oriented.component_orientations.tsv', '--region', 'chr1:1-100',
        '--output-vcf', output)
    pysam.tabix_index(str(output), preset='vcf')
    with pysam.VariantFile(output) as result:
        calls = {r.pos: r.samples['SAMPLE'] for r in result.fetch('chr1', 0, 100)}
    assert [calls[p]['PS'] for p in sites] == [10, 10, 10, 10, 80]
    assert calls[10]['GT'] == (0, 1)
    assert calls[40]['GT'] == (1, 0)
    assert calls[80]['GT'] == (0, 1)
