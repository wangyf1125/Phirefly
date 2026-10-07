# Phirefly

Phirefly phases heterozygous SNPs from long reads and extends phase sets with
Hi-C contacts. It compresses read evidence into phaselets and hyperreads, then
uses one CFC batch to sample the remaining orientations. Selection uses the
original read objective, with contact support when Hi-C is supplied. Phasing does
not require a truth VCF.

## Installation

Linux x86-64, Python 3.11, CPU CFC:

~~~sh
unset PYTHONPATH PYTHONHOME
export PYTHONNOUSERSITE=1
conda create -n phirefly --override-channels -c y1fei -c conda-forge \
    -c bioconda python=3.11 y1fei::phirefly=0.2.2
conda activate phirefly
phirefly --version
PHIREFLY_PYTHON="$CONDA_PREFIX/bin/python" PHIREFLY_BENCHMARK=1 bash \
    "$CONDA_PREFIX/share/phirefly/examples/mhc_smoke/run_mhc_smoke.sh" out/mhc
~~~

The [y1fei personal channel](https://anaconda.org/y1fei/phirefly) provides
Phirefly and its required MindQuantum/Rich builds. Conda installs all dependencies,
including htslib and WhatsHap; no separate pip step is needed. This is not a
Bioconda release. The bundled MHC example includes HG002 truth for HE/SE checks.
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
one CFC batch. Local phaselets remain mutable; retained contact links can join
previously disconnected long-read components. See the [Hi-C example](docs/hic.md).
The joint interface is included in version 0.2.2; older versions use a different
Hi-C model. Use fresh output directories when changing versions.
Use phirefly benchmark --help for truth-based evaluation. The Conda package
includes its dependencies; pip installations need the benchmark extra.
See [metric definitions](docs/reproducibility.md).

## Development

~~~sh
python -m pip install -e '.[qaia,test]'
pytest -q
~~~

Start with phirefly/pipeline.py; it follows extraction, MSF, phaselet/hyperread
construction, one CFC solve and VCF export. The extract/, phaselet/,
hyperread/ and hic/ directories contain the corresponding algorithms.
See [release checks](docs/release.md). MIT license.
