# MHC smoke

This staged HG002 MHC example contains 12,259 reads, 10,200 observed SNPs and
498,889 observations. It exercises single-round phaselet/hyperread solving and
indexed VCF export, not BAM extraction.

~~~sh
python -I -m pip install '.[qaia]'
bash examples/mhc_smoke/run_mhc_smoke.sh
~~~

The default is CFC on cpu-float32, 20 candidates, 200 iterations and seed 23.
Truth is not used by the solver. The included truth data are for separate evaluation only.
The script runs the installed package in Python isolated mode, ignoring inherited
PYTHONPATH/PYTHONHOME and checkout imports. Activate the intended environment
first, or set PHIREFLY_PYTHON to its Python executable. bgzip and tabix must be
on PATH. Missing prerequisites fail before computation or output creation.
The old reference table predates the current inference revision and was removed;
do not compare its phasing errors as a regression target.

Use a new output directory for another run:

~~~sh
bash examples/mhc_smoke/run_mhc_smoke.sh out/mhc_retry
~~~

The script never deletes or silently reuses an earlier output. Its short terminal
summary omits unevaluated SE/HE; full settings and counts are in config_metrics.tsv.
The workflow timing starts from staged observations, not from BAM extraction.

## Hamming and switch errors

~~~sh
python -I -m pip install '.[qaia,benchmark]'
PHIREFLY_BENCHMARK=1 bash examples/mhc_smoke/run_mhc_smoke.sh out/mhc_evaluated
~~~

This first completes the same phasing run, then evaluates it against the bundled
truth. It prints HE/SE and their error counts/denominators, and writes
benchmark_metrics.tsv beside the phased VCF. Evaluation is not included in
phasing time. The phasing-only config_metrics.tsv remains free of truth metrics.
The y1fei Conda package already includes WhatsHap 2.8 and these example files.
It does not need a separate pip installation or repository checkout:

~~~sh
PHIREFLY_BENCHMARK=1 bash \
    "$CONDA_PREFIX/share/phirefly/examples/mhc_smoke/run_mhc_smoke.sh" out/mhc_evaluated
~~~

The truth files are data/truth.region.bcf and its .csi index: HG002 (NA24385),
GRCh38, chr6:28510120-33480577 (1-based inclusive), extracted from
GIAB/NIST HG002_GRCh38_v5.0q_smvar.vcf.gz. This assembly-based draft reference
is distinct from the older GIAB v4.2.1 benchmark. The slice has 29,000 records,
including 13,131 phased biallelic heterozygous SNPs. Evaluation uses matching
CHROM/POS/REF/ALT SNPs and joint prediction/reference phase blocks; it is a smoke
test, not a claim of full MHC coverage or an independent benchmark.

Source: [GIAB v5.0q release](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/v5.0q/).
SHA256 of the bundled BCF:
e5c94d0a08ab16bcd222e7ce13f80579a7b4fdc050d89ffbfd17f1e5a7ec740b.

## GPU

~~~sh
PHIREFLY_BACKEND=gpu-float32 PHIREFLY_BATCH_SIZE=100 PHIREFLY_N_ITER=1000 \
    bash examples/mhc_smoke/run_mhc_smoke.sh out/mhc_gpu
~~~

GPU execution requires the gpu extra and a configured GPU environment.
Inspect qaia_summary.tsv for pipeline_revision=single_cfc_v1 and qaia_calls=1.
Exact candidate scores can vary with the numerical backend. The qaia-marked
tests also cover fresh synthetic BAM-to-VCF extraction and export.
