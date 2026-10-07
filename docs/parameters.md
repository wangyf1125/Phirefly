# Parameters

Start with `phirefly run --bam reads.bam --vcf sites.vcf.gz --sample SAMPLE
--region chr1:1-248956422 --out-dir out/chr1`. No tuning is needed for the
release configuration. `phirefly run --help` shows common options;
`phirefly run --help-all` describes every accepted option and its default.
The short-help interface is a source update after 0.2.1 build 0.

| Common option | Meaning |
| --- | --- |
| `--bam` | Coordinate-sorted, indexed long-read alignments |
| `--hic-bam` | Optional indexed same-sample Hi-C alignments for a joint CFC solve |
| `--vcf` | Indexed heterozygous biallelic SNP input, on the same reference |
| `--sample` | VCF sample to phase; defaults to the first sample |
| `--region` | Chromosome interval, 1-based and inclusive |
| `--out-dir` | New output directory |
| `--threads` | BAM extraction workers, not a global solver thread limit |
| `--backend` | `cpu-float32`, or `gpu-float32` with a configured GPU |
| `--seed` | Random initialization seed, default 23 |
| `--debug-output` | Optional detailed graphs and tables |
| `--force` | Recompute in an existing directory, not resume |

## Advanced settings

- Candidate budget: `--batch-size 100 --n-iter 1000`. Memory grows with batch
  size and graph size; more updates cost more time but do not guarantee better
  phasing. `--seed-scale` and `--seed-noise` control continuous initialization.
- Extraction: `--min-baseq 13 --min-mapq 20`. `--weighted` enables log-odds
  weights; the default release method uses unit observation weights.
- Hard links: `--bridge-mode risk --risk-threshold 0.45`. A lower threshold
  is more conservative. Risk is a heuristic score, not an error probability.
  `--min-confidence` applies to confidence mode, not the default risk filter.
- Soft links: `--soft-bridge-mode phaselet --soft-risk-threshold 0.75
  --soft-edge-scale 0.25`. These retain eligible uncertain links without hard
  merging; the scale changes their influence in the compressed objective.
- Hyperreads: `--hyperread-norm sqrt_cap --hyperread-cap 10` normalizes support
  by the square root of the collapsed read count and caps coupling magnitude.
- `--phaselet-config-name` is an output label, not a preset loader.
- `--trace-every` enables instrumentation and adds overhead. Do not report
  traced elapsed time as an uninstrumented solver speed measurement.
- Truth is only for evaluation. Prefer a separate `phirefly benchmark` command
  rather than `run --truth-vcf` when measuring inference time.
- Hi-C: `--hic-min-mapq 20 --hic-min-baseq 13`; primary nonduplicate,
  QC-passing alignments only. `--hic-min-contacts 1` applies to each phaselet
  pair. New links between original long-read components also require
  `--hic-min-molecules 20` external contact read names per component.
  `--hic-scale 0.1 --hic-cap 10` controls their normalized influence. These
  are support filters, not calibrated correctness probabilities.

Hidden options remain accepted for reproducibility. Changing visibility does
not change defaults or inference. Direct `phaselet-qaia` is an expert command
for staged observations and a precomputed MSF; ordinary users should use `run`.
