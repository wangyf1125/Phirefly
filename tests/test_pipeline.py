import json
import csv
import shutil

import pysam
import pytest

from phirefly.core.io import require_fresh_output
from phirefly.pipeline import build_parser, run_pipeline


def test_existing_outputs_require_explicit_recompute(tmp_path):
    require_fresh_output(tmp_path)
    marker = tmp_path / "old-result"
    marker.write_text("preserve")
    with pytest.raises(ValueError, match="not empty"):
        require_fresh_output(tmp_path)
    require_fresh_output(tmp_path, force=True)
    assert marker.read_text() == "preserve"


def toy_inputs(tmp_path):
    bam = tmp_path / "reads.bam"
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": "chr1", "LN": 100}]}
    reads = [("left1", 9, 11, "A"), ("left2", 9, 11, "A"),
             ("bridge", 19, 21, "A"),
             ("right1", 39, 11, "C"), ("right2", 39, 11, "C")]
    with pysam.AlignmentFile(bam, "wb", header=header) as out:
        for name, start, length, base in reads:
            aln = pysam.AlignedSegment()
            aln.query_name = name
            sequence = base * length
            if name == "bridge":
                sequence = sequence[:-1] + "C"
            aln.query_sequence = sequence
            aln.flag = 0
            aln.reference_id = 0
            aln.reference_start = start
            aln.mapping_quality = 60
            aln.cigarstring = f"{length}M"
            aln.query_qualities = pysam.qualitystring_to_array("I" * length)
            out.write(aln)
    pysam.index(str(bam))
    vcf = tmp_path / "sites.vcf"
    vcf.write_text('##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100>\n'
                   '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                   '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n' +
                   ''.join(f'chr1\t{pos}\t.\tA\tC\t.\tPASS\t.\tGT\t0/1\n'
                           for pos in (10, 20, 40, 50)))
    gz = str(vcf) + ".gz"
    pysam.tabix_compress(str(vcf), gz)
    pysam.tabix_index(gz, preset="vcf")
    return bam, gz


@pytest.mark.qaia
@pytest.mark.parametrize("threads", [1, 2])
def test_bam_to_indexed_vcf_without_truth(tmp_path, threads):
    pytest.importorskip("mindquantum")
    if not shutil.which("bgzip") or not shutil.which("tabix"):
        pytest.skip("htslib commands are required")
    bam, vcf = toy_inputs(tmp_path)
    out = tmp_path / "out"
    args = build_parser().parse_args([
        "--bam", str(bam), "--vcf", vcf, "--sample", "SAMPLE",
        "--region", "chr1:1-100", "--out-dir", str(out),
        "--threads", str(threads), "--chunk-size-bp", "20",
        "--batch-size", "4", "--n-iter", "20", "--backend", "cpu-float32"])
    run_pipeline(args)
    outputs = list(out.rglob("phased.component_ps.vcf.gz"))
    assert len(outputs) == 1
    with pysam.VariantFile(outputs[0]) as result:
        calls = list(result.fetch("chr1", 0, 100))
    assert len(calls) == 4
    assert all(r.samples["SAMPLE"].phased for r in calls)
    assert len({r.samples["SAMPLE"]["PS"] for r in calls}) == 1
    summary = next(out.rglob("qaia_summary.tsv")).read_text()
    assert "qaia_calls\t1\n" in summary
    assert "phaselets\t2\n" in summary
    assert "hyperreads\t1\n" in summary
    manifest = json.loads((out / "run_manifest.json").read_text())
    assert manifest["pipeline_revision"] == "single_cfc_v1"
    assert manifest["arguments"]["truth_vcf"] is None
    assert manifest["source_sha256"]["qaia.py"]
    assert not list(out.rglob("hyperread_matrix.npz"))
    assert not list(out.rglob("phaselet_edges.tsv"))
    assert not list(out.rglob("phaselet_hyperread.metrics.tsv"))
    with pytest.raises(ValueError, match="not empty"):
        run_pipeline(args)


@pytest.mark.qaia
def test_optional_benchmark_and_debug_outputs(tmp_path):
    pytest.importorskip("mindquantum")
    if not shutil.which("bgzip") or not shutil.which("tabix"):
        pytest.skip("htslib commands are required")
    bam, vcf = toy_inputs(tmp_path)
    truth = tmp_path / "truth.vcf"
    truth.write_text((tmp_path / "sites.vcf").read_text().replace("0/1", "0|1"))
    out = tmp_path / "evaluated"
    args = build_parser().parse_args([
        "--bam", str(bam), "--vcf", vcf, "--sample", "SAMPLE",
        "--region", "chr1:1-100", "--out-dir", str(out),
        "--truth-vcf", str(truth), "--truth-sample", "SAMPLE", "--debug-output",
        "--batch-size", "4", "--n-iter", "20", "--backend", "cpu-float32"])
    run_pipeline(args)
    assert len(list(out.rglob("hyperread_matrix.npz"))) == 1
    assert len(list(out.rglob("phaselet_edges.tsv"))) == 1
    with next(out.rglob("phaselet_hyperread.metrics.tsv")).open() as handle:
        metrics = next(csv.DictReader(handle, delimiter="\t"))
    with next(out.rglob("config_metrics.tsv")).open() as handle:
        config = next(csv.DictReader(handle, delimiter="\t"))
    assert config["benchmark_metrics"] == "1"
    assert float(config["HE_percent"]) == float(metrics["HE_percent"])
    assert int(metrics["phased_snps"]) == 4
