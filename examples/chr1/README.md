# HG002 chromosome 1

This is a full-chromosome staged-input example, not a subsampled MHC example.
It contains 565,291 reads, 168,950 SNPs and 8,293,236 read-SNP observations.
The final graph has 1,919 phaselets, 1,965 hyperreads and 7,995 unique couplings.
Input files are checked against SHA256SUMS before solving.

From a cloned repository, with the Phirefly 0.2.1 Conda environment active:

~~~sh
PHIREFLY_PYTHON="$CONDA_PREFIX/bin/python" PHIREFLY_BENCHMARK=1 \
  /usr/bin/time -v -o chr1.resource.log \
  bash examples/chr1/run_chr1.sh out/chr1-example
~~~

Use a fresh output directory. The example defaults to CPU CFC, batch 100,
1000 iterations and seed 23. Set PHIREFLY_BACKEND=gpu-float32 only with a
configured GPU installation. Evaluation is optional and never feeds the solver.
The result and benchmark_metrics.tsv are under
out/chr1-example/chr1/alg_CFC/phaselet_risk0p45_soft0p75_s0p25/.
This example is in the Git repository and the 0.2.2 source archive. It is not
bundled in Conda 0.2.1 build 0. The 0.2.2 recipe installs it under
share/phirefly/examples/chr1; package publication is a separate release step.

## Scope and resources

The input NPZ and MSF signs were produced by the frozen single-round workflow
on 24 September 2026. The script rebuilds phaselets/hyperreads, runs one CFC
batch, selects candidates and exports VCF. It does NOT measure BAM extraction
or recomputation of the MSF. `/usr/bin/time -v` includes optional evaluation;
config_metrics.tsv reports the staged phasing time separately.

The final graph has 3,884 variables but the raw problem has 734,241 variables.
Graph-only resource estimates do not describe the observation pipeline's peak
memory. See resources.tsv and resources.csv for a dated CPU example measurement;
hardware/backend and timing scope matter, so do not extrapolate from MHC alone.

On an Intel Xeon Gold 6526Y CPU, the staged example plus evaluation took
86.97 s with a measured peak RSS of 2.76 GiB. The timed CFC call took 30.383 s;
the staged phaselet workflow took 67.63 s. Four CPUs and 16 GiB were allocated;
these are one run's measurements, not minimum resource requirements or a GPU
speed comparison. The sample SE/HE counts are in benchmark_reference.tsv and
benchmark_reference.csv. They are a reproducibility record, not strict
cross-platform pass/fail thresholds.

To measure complete BAM-to-VCF resources, provide the original indexed BAM:

~~~sh
/usr/bin/time -v -o chr1.e2e.resource.log \
  phirefly run --bam /path/to/HG002.Revio.bam \
  --vcf examples/chr1/data/shared.snps.vcf.gz --sample HG002 \
  --region chr1:1-248956422 --out-dir out/chr1-e2e --threads 8
~~~

## Data and graph identity

HG002 / NA24385, GRCh38 chromosome 1, interval 1-248956422. Long reads are
[PacBio HiFi Revio, nominal 48x, 2023-10-31](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/PacBio_HiFi-Revio_20231031/HG002_PacBio-HiFi-Revio_20231031_48x_GRCh38-GIABv3.bam).
Truth is the [GIAB/NIST v5.0q assembly-based draft](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/v5.0q/),
not v4.2.1. Benchmark input genotypes/sites were derived from this reference
and unphased, so this is not an independent variant-calling benchmark.

- phirefly_reads.observations.npz: integer-indexed read-SNP observations.
- msf_backbone.snp_phases.tsv: inference-only local signs; not truth phases.
- shared.snps.vcf.gz and .tbi: unphased HG002 input and index.
- truth.region.bcf and .csi: reference for separate HE/SE evaluation.
- hyperread_matrix.npz: saved final hyperread-by-phaselet coupling matrix A.
- soft_phaselet_coupling.npz: symmetric final phaselet coupling matrix S.
- phaselet_mapping.npz: mapping back to SNPs and local signs.

The last three files are for graph/resource inspection; the example rebuilds
its inference graph from observations rather than loading cached solutions.
For Ising inspection the symmetric matrix is J = [[0, A], [A.T, S]], ordered
hyperreads first, with energy -0.5 * spins.T @ J @ spins. Disconnected
components can be submitted separately; a width audit alone does not validate
routing, circuit depth, noise, shot budget or QPU feasibility.
