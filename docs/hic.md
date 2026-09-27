# Hi-C

Run after long-read phasing. The contact BAM must use the same reference and
sample. Components follow the supplied phase-set identifiers, not a separately
reconstructed read graph. Only groups joined by retained contact evidence share
a new PS; disconnected groups are not merged to increase contiguity.

~~~sh
P=out/chr20/phaselet_qaia/phirefly/alg_CFC/phaselet_risk0p45_soft0p75_s0p25
O=out/chr20/extract/phirefly_reads.observations.npz
H=out/chr20/hic/phirefly
phirefly hic-edges --observations "$O" --snp-phases "$P/snp_phases.tsv" \
    --contact-bam hic.bam --region chr20:1-64444167 --out-prefix "$H"
phirefly hic-orient --edges "$H.block_edges.tsv" --out-prefix "$H.msf" \
    --min-contacts 20 --min-confidence 0.8
phirefly hic-apply --input-vcf het_sites.vcf.gz --sample SAMPLE \
    --observations "$O" --snp-phases "$P/snp_phases.tsv" \
    --component-orientations "$H.msf.component_orientations.tsv" \
    --region chr20:1-64444167 --output-vcf out/chr20/phirefly.hic.vcf
bgzip -c out/chr20/phirefly.hic.vcf > out/chr20/phirefly.hic.vcf.gz
tabix -p vcf out/chr20/phirefly.hic.vcf.gz
~~~

Orientation uses a maximum-spanning forest. This optimizes the selected forest,
not the complete cyclic contact graph. A long PS span is not independently
validated long-range accuracy. Report source libraries, filtering, join evidence
and evaluation scope with any Hi-C result.
