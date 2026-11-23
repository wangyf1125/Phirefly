# Phirefly

Phirefly phases heterozygous SNPs from long reads and extends phase sets with
Hi-C contacts. It compresses read evidence into phaselets and hyperreads, then
samples phaselet orientations and hyperread labels together. For each connected
component it keeps the sampled assignment with the lowest Ising energy. Hi-C
contacts, when supplied, are included in that energy. Phasing does not use a
truth VCF.

## Installation

To run the publication version (v0.2.3) with the CPU CFC sampler:

~~~sh
unset PYTHONPATH PYTHONHOME
export PYTHONNOUSERSITE=1
git clone --branch v0.2.3 https://github.com/wangyf1125/Phirefly.git
cd Phirefly
conda env create -f environment.yml
conda activate phirefly
python -m pip install '.[qaia,benchmark]'
phirefly --version
PHIREFLY_PYTHON="$CONDA_PREFIX/bin/python" PHIREFLY_BENCHMARK=1 bash \
    examples/mhc_smoke/run_mhc_smoke.sh out/mhc
~~~

The older 0.2.2 installation uses a different selection rule. In v0.2.3,
`qaia_summary.tsv` reports `selection_objective=component_minimum_energy`
(as a key/value row). The bundled MHC example includes HG002 truth for HE/SE checks.
See [installation notes](docs/install.md) for source installation and GPU support,
and [packaging/conda](packaging/conda/README.md) for recipes.

Examples: [MHC installation smoke](examples/mhc_smoke) and
[full HG002 chromosome 1](examples/chr1), with truth, checksums and resource
measurement commands. The chromosome example starts from observations; the
separate BAM-to-VCF command measures the complete workflow.

## Usage

`phirefly run --help` lists common options; `--help-all` includes advanced
settings and defaults. See [parameter notes](docs/parameters.md).

Input: an indexed, coordinate-sorted BAM and an indexed biallelic heterozygous
SNP VCF on the same reference. Coordinates are 1-based, inclusive.

~~~sh
phirefly run --bam reads.bam --vcf het_sites.vcf.gz --sample SAMPLE \
    --region chr20:1-64444167 --out-dir out/chr20 \
    --backend cpu-float32 --threads 8
~~~

The indexed result is:

~~~text
out/chr20/phaselet_qaia/phirefly/alg_CFC/phaselet_risk0p45_soft0p75_s0p25/phased.component_ps.vcf.gz
~~~

Use a new output directory for each configuration. Existing results are never
reused implicitly; --force recomputes the workflow. The run_manifest.json records
arguments, dependency versions and source checksums. Use --debug-output for graph
tables. The terminal prints a short summary; config_metrics.tsv keeps the full
record. GPU execution requires --backend gpu-float32 and a configured GPU.

## Hi-C and Evaluation

Supply `--hic-bam hic.bam` to solve long-read and contact evidence jointly in
one CFC batch. Contacts can change phaselet orientations, not their internal
SNP phases. Retained contact links can join
previously disconnected long-read components. See the [Hi-C example](docs/hic.md).
Use a fresh output directory when changing versions. Use `phirefly benchmark
--help` for truth-based evaluation.
See [metric definitions](docs/reproducibility.md) and [publication data and offline quantum replay](publication/README.md).

## Development

~~~sh
python -m pip install -e '.[qaia,test]'
pytest -q
~~~

Start with phirefly/pipeline.py; it follows extraction, MSF, phaselet/hyperread
construction, one CFC solve and VCF export. The extract/, phaselet/,
hyperread/ and hic/ directories contain the corresponding algorithms.
See [release checks](docs/release.md). MIT license.
