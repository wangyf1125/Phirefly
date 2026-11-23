# Hi-C

Provide indexed, coordinate-sorted long-read and Hi-C BAMs from the same sample
and reference. Use a fresh output directory for each run.

```sh
phirefly run --bam reads.bam --hic-bam hic.bam \
    --vcf het_sites.vcf.gz --sample HG002 --region chr20:1-64444167 \
    --out-dir out/chr20-hic --backend cpu-float32 --threads 8
```

Hi-C anchors require primary, nonduplicate, QC-passing alignments, mapping
quality >=20 and SNP base quality >=13. Missing base qualities are excluded.
Each alignment contributes its strongest phaselet/sign anchor; distinct read
names are counted once per phaselet pair. Paired names ending in /1 or /2 are
normalized to the same name.

Anchor strength is the absolute local allele sum times capped MQ/60. A pair
uses the smaller endpoint strength (floor 0.1) and their sign product.
Signed aggregate support is scaled by 0.1 and capped at absolute weight 10.
Links between original long-read components require at least 20 external
contact read names at each component. This is not a PCR-independence guarantee.
The corresponding options appear in `phirefly run --help-all`.

One CFC batch samples the joint graph. Each component keeps the sampled
assignment with the lowest energy, including its Hi-C couplings. A phaselet
can change orientation, but the relative phases of its SNPs stay fixed.
Phase sets follow proposal and contact connectivity even when weights cancel.

Use `--debug-output` for `hic_contacts.tsv` and `hic_coupling.npz`.
`qaia_summary.tsv`, `config_metrics.tsv` and `run_manifest.json` record the
contact configuration and source identity. Do not confuse a long PS span with
overlap-split phase-block N50 or independently validated long-range accuracy.
