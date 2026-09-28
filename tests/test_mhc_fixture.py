"""The distributed truth is indexed, phased and overlaps the smoke input."""
from pathlib import Path

import pysam


def test_bundled_mhc_truth():
    data = Path(__file__).resolve().parents[1] / 'examples/mhc_smoke/data'
    with pysam.VariantFile(str(data / 'truth.region.bcf')) as truth:
        assert list(truth.header.samples) == ['HG002']
        calls = [r for r in truth.fetch('chr6', 28510119, 33480577)
                 if len(r.ref) == 1 and r.alts and len(r.alts) == 1
                 and len(r.alts[0]) == 1 and r.samples['HG002']['GT'] in ((0, 1), (1, 0))]
    assert len(calls) == 13131
    assert all(r.samples['HG002'].phased for r in calls)
    keys = {(r.chrom, r.pos, r.ref, r.alts) for r in calls}
    with pysam.VariantFile(str(data / 'shared.snps.vcf.gz')) as source:
        shared = {(r.chrom, r.pos, r.ref, r.alts) for r in source.fetch('chr6', 28510119, 33480577)}
    assert shared and shared <= keys
