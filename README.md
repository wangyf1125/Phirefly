# Phirefly

Phirefly is a haplotype phasing tool for long reads.  It builds conservative
local phaselets, compresses spanning reads into hyperreads, and solves the
remaining block-orientation problem with a quantum-inspired QAIA solver.

Truth data are not needed for phasing.  Truth VCFs are used only for
benchmarking.

## Install

Phirefly supports Python 3.10 and 3.11.

```bash
git clone https://github.com/wangyf1125/Phirefly.git
cd Phirefly
python -m pip install -e .
```

Install the QAIA dependency when running the solver:

```bash
python -m pip install -e '.[qaia]'
```

For the full research environment:

```bash
python -m pip install -e '.[all]'
```

## Smoke Test

```bash
bash examples/mhc_smoke/run_mhc_smoke.sh
```

This runs a small packaged MHC example and checks that extraction-free
phaselet/hyperread phasing works.

## Basic Run

Run Phirefly from BAM and heterozygous SNP VCF:

```bash
phirefly run \
  --bam reads.bam \
  --vcf het_sites.vcf.gz \
  --sample SAMPLE \
  --region chr20:1-64444167 \
  --out-dir out/chr20 \
  --backend cpu-float32 \
  --threads 8
```

The main output is:

```text
out/chr20/phased.component_ps.vcf.gz
```

## Benchmark

Benchmarking is separate from phasing:

```bash
phirefly benchmark \
  --input-vcf het_sites.vcf.gz \
  --truth-vcf truth.phased.vcf.gz \
  --region chr20:1-64444167 \
  --pred-sample SAMPLE \
  --out-tsv out/benchmark.metrics.tsv \
  --method phirefly:out/chr20/phased.component_ps.vcf.gz
```

Reported metrics include switch error, Hamming error, phase-block N50, SNP
completeness, and runtime.

## Hi-C / Pore-C

Hi-C and Pore-C are optional final-method extensions for orienting Phirefly
components into longer phase blocks.  The public release path uses conservative
component-edge filtering and MSF orientation.

```bash
phirefly hic-edges \
  --observations out/chr20/extract/phirefly_reads.observations.npz \
  --snp-phases out/chr20/phaselet_qaia/phirefly/alg_CFC/phaselet_risk0p45_soft0p75_s0p25/snp_phases.tsv \
  --contact-bam hic_or_porec.bam \
  --region chr20:1-64444167 \
  --out-prefix out/chr20/hic/phirefly_hic

phirefly hic-orient \
  --edges out/chr20/hic/phirefly_hic.block_edges.tsv \
  --out-prefix out/chr20/hic/phirefly_hic.msf \
  --min-contacts 20 \
  --min-confidence 0.8

phirefly hic-apply \
  --input-vcf het_sites.vcf.gz \
  --sample SAMPLE \
  --snp-phases out/chr20/phaselet_qaia/phirefly/alg_CFC/phaselet_risk0p45_soft0p75_s0p25/snp_phases.tsv \
  --observations out/chr20/extract/phirefly_reads.observations.npz \
  --component-orientations out/chr20/hic/phirefly_hic.msf.component_orientations.tsv \
  --region chr20:1-64444167 \
  --output-vcf out/chr20/phirefly.hic.vcf
```

## Commands

```text
phirefly run              one-command phasing workflow
phirefly extract          BAM/VCF to read-SNP observations
phirefly msf              MSF backbone
phirefly phaselet-qaia    phaselet/hyperread QAIA solver
phirefly vcf              VCF export
phirefly benchmark        metrics against truth
phirefly hic-edges        Hi-C/Pore-C component edges
phirefly hic-orient       component orientation
phirefly hic-apply        apply component orientations
phirefly hic-qc           Hi-C edge QC
```

The same tools are also installed with `phirefly-*` entry points, for example
`phirefly-run` and `phirefly-phaselet-qaia`.

## Code Map

Start here:

```text
phirefly/pipeline.py          top-level workflow
phirefly/extract/             observation extraction
phirefly/msf.py               MSF backbone
phirefly/phaselet/runner.py   phaselet/hyperread solver path
phirefly/hyperread/build.py   hyperread compression
phirefly/qaia.py              QAIA wrapper
phirefly/hic/                 Hi-C/Pore-C extension
phirefly/vcf.py               VCF export
phirefly/metrics.py           benchmark metrics
```

## Development

```bash
python -m pip install -e '.[test]'
pytest -q
```

## License

MIT.
