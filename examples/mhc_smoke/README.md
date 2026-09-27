# MHC smoke

This staged HG002 MHC example contains 12,259 reads, 10,200 observed SNPs and
498,889 observations. It exercises single-round phaselet/hyperread solving and
indexed VCF export, not BAM extraction.

~~~sh
python -m pip install '.[qaia]'
bash examples/mhc_smoke/run_mhc_smoke.sh
~~~

The default is CFC on cpu-float32, 20 candidates, 200 iterations and seed 23.
Truth is not used. The included truth data are for separate evaluation only.
The old reference table predates the current inference revision and was removed;
do not compare its phasing errors as a regression target.

Use a new output directory for another run:

~~~sh
PHIREFLY_BACKEND=gpu-float32 PHIREFLY_BATCH_SIZE=100 PHIREFLY_N_ITER=1000 \
    bash examples/mhc_smoke/run_mhc_smoke.sh out/mhc_gpu
~~~

GPU execution requires the gpu extra and a configured GPU environment.
Inspect qaia_summary.tsv for pipeline_revision=single_cfc_v1 and qaia_calls=1.
Exact candidate scores can vary with the numerical backend. The qaia-marked
tests also cover fresh synthetic BAM-to-VCF extraction and export.
