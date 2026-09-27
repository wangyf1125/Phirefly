# Phirefly

Phirefly phases heterozygous SNPs from long reads and extends phase sets with
Hi-C contacts. It compresses read evidence into phaselets and hyperreads, then
uses one CFC batch to sample the remaining orientations. Selection uses the
original read objective. Phasing does not require a truth VCF.

## Installation

Python 3.10/3.11 and htslib are required. On Linux:

~~~sh
git clone --branch v0.2.0 https://github.com/wangyf1125/Phirefly.git
cd Phirefly
conda env create -f environment.yml
conda activate phirefly
python -m pip install '.[qaia]'
phirefly --version
bash examples/mhc_smoke/run_mhc_smoke.sh
~~~

This uses Conda for the environment and pip for Phirefly/MindQuantum.
Phirefly is not yet on Bioconda; **conda install phirefly is not available**.
See [installation notes](docs/install.md) for GPU and Bioconda details.

## Usage

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
tables. GPU execution requires --backend gpu-float32 and a configured GPU.

## Hi-C and Evaluation

Hi-C is part of Phirefly: hic-edges, hic-orient and hic-apply extend the
long-read phase sets using contact evidence. See the [Hi-C example](docs/hic.md).
Use phirefly benchmark --help for truth-based evaluation; install the
benchmark extra first. See [metric definitions](docs/reproducibility.md).

## Development

~~~sh
python -m pip install -e '.[qaia,test]'
pytest -q
~~~

Start with phirefly/pipeline.py; it follows extraction, MSF, phaselet/hyperread
construction, one CFC solve and VCF export. The extract/, phaselet/,
hyperread/ and hic/ directories contain the corresponding algorithms.
See [release checks](docs/release.md). MIT license.
