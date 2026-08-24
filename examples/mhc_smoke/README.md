# MHC Smoke Test

This example runs the phaselet-hyperread QAIA path on a small HG002 MHC
read-SNP observation graph.

The bundled inputs are a small MHC smoke dataset:

```text
phirefly_reads.observations.npz
shared.snps.vcf / shared.snps.vcf.gz
truth.region.bcf
phirefly_reads.msf.snp_phases.tsv
```

Run from the repository root:

```bash
python -m pip install -e '.[qaia]'
bash examples/mhc_smoke/run_mhc_smoke.sh
```

By default the local smoke test uses `cpu-float32`, `batch_size=20`, and
`n_iter=200`.  It still requires MindQuantum QAIA to be installed.  A plain
system Python with only the standard library is not enough; the minimum Python
dependencies are `numpy`, `scipy`, `pysam`, and `mindquantum`.

For a GPU QAIA smoke run, install the GPU extra first because MindQuantum's
`gpu-float32` backend requires PyTorch:

```bash
python -m pip install -e '.[gpu]'
```

Then use:

```bash
PHIREFLY_BACKEND=gpu-float32 PHIREFLY_BATCH_SIZE=100 PHIREFLY_N_ITER=1000 \
  bash examples/mhc_smoke/run_mhc_smoke.sh
```

The smoke run on H100 with the packaged MHC smoke data produced:

```text
SE %       0.255503
HE %       0.196098
N50 kb     308.971
SNP %      77.678775
reads      12259
snps       10200
phaselets  65
hyperreads 75
nodes      140
```

The full expected key/value record is in:

```text
examples/mhc_smoke/expected/mhc_smoke_reference.tsv
```

Exact QAIA scores can vary across backend/version settings. The reference file
tracks the H100 `gpu-float32`, `CFC`, `batch_size=100`, `n_iter=1000` smoke run
for this staged example data.
