# Reproducibility

The inference revision is single_cfc_v1: reads and original MSF signs first
release suspect internal constraints, then one CFC batch samples the resulting
graph. Candidate selection uses the original read objective, independently only
where it separates. The public defaults remain 100 candidates, 1000 iterations
and seed 23. Solver-control experiments do not change these defaults.

Use a fresh output directory after changing inputs or parameters. The --force
option recomputes all stages; it is not a resume operation. The output manifest
records source SHA256 values, package/dependency versions, parameters and input
path/size/modification time. Input stat records are not content hashes: archival
reproduction additionally requires input checksums and the matching source tag.

## Metrics

- span_N50_kb: N50 of inclusive, raw phase-set spans.
- haplotype_N50_kb: legacy column name for that same span metric.
- phased_block_N50_kb: N50 of WhatsHap's nonoverlapping split-block spans.
- Phase-block NG50: uses the chromosome/reference length as denominator;
  it is not another name for phase-block N50.
- Input phased (%): phased input heterozygous SNPs divided by input sites.
- Reference SNP coverage (%): phased predictions matching reference
  heterozygous SNPs divided by reference sites. Neither fraction is read depth.

Package benchmarking uses WhatsHap 2.8 for block statistics. The manuscript
additionally uses a pinned evaluator with explicit comparison denominators.
Report evaluator version, region, prediction/reference samples and PS policy.
Do not substitute a span N50 when block statistics are missing.

## Timing and Validation

Paper E2E begins with staged BAM/VCF and ends with indexed VCF, excluding queueing,
alignment, variant calling and evaluation. Phirefly Core measures the solver;
comparator Core scopes differ. Pipeline timings include optional evaluation if
--truth-vcf is supplied, so do not use that mode for the paper's phasing E2E.

The single-round HG002 baseline and fresh Hi-C workflows completed 22 autosomes.
A separate CFC/LQA follow-up completed 220 fresh executions at 3000 iterations
and five seeds per solver. It does not replace the default 1000-iteration result.
Validation on independent samples with this revision remains outstanding.

IBM hardware experiments sample the fixed final graph using QAOA and retain
classical preprocessing and selection. They are research experiments outside
the release CLI, not evidence of quantum advantage or a dependency of Phirefly.

## Output records

The terminal is a short human-readable summary, not a TSV interface.
config_metrics.tsv retains the detailed machine-readable record. Columns named
initial_* describe parity classification before the read-consistency split;
phaselets, hyperreads and couplings describe the final graph. For compatibility,
parity_edges, soft_bridge_edges, weak_phaselet_edges, conflicting_phaselet_edges
and the bridge-risk summaries alias their initial_* values. They are not final
soft-coupling counts. retained_phaselet_edges describes the final forest size.

unique_couplings counts each undirected term once; coupling_nnz counts both
entries of the symmetric matrix. qaia_runtime_s is CFC-only; solver_runtime_s is
the phaselet computation stage including input loading/build/selection, and
config_wall_s includes phaselet output/export. Neither is BAM-to-VCF E2E; the
one-command workflow records t_total_e2e in timings.json. Missing truth metrics
are blank, not zero. The resolved input sample is always recorded.
