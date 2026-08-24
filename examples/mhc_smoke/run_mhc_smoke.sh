#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
DATA="$ROOT/examples/mhc_smoke/data"
OUT="${1:-$ROOT/examples/mhc_smoke/output}"

mkdir -p "$OUT"

export PYTHONPATH="$ROOT${PYTHONPATH:+:$PYTHONPATH}"

python -m phirefly.phaselet \
  --region-key mhc \
  --region-label MHC \
  --region chr6:28510120-33480577 \
  --observations "$DATA/phirefly_reads.observations.npz" \
  --input-vcf-gz "$DATA/shared.snps.vcf.gz" \
  --input-vcf-plain "$DATA/shared.snps.vcf" \
  --truth-vcf "$DATA/truth.region.bcf" \
  --tau-snp-phases "$DATA/phirefly_reads.msf.snp_phases.tsv" \
  --phaselet-config-name smoke_risk_soft \
  --min-abs-support 2.0 \
  --min-count 2 \
  --min-confidence 0.65 \
  --max-distance-bp 1000000 \
  --algorithm CFC \
  --backend "${PHIREFLY_BACKEND:-cpu-float32}" \
  --batch-size "${PHIREFLY_BATCH_SIZE:-20}" \
  --n-iter "${PHIREFLY_N_ITER:-200}" \
  --seed 23 \
  --seed-scale 0.001 \
  --seed-noise 0.003 \
  --bridge-mode risk \
  --risk-threshold 0.45 \
  --soft-bridge-mode phaselet \
  --soft-min-abs-support 1.0 \
  --soft-min-count 1 \
  --soft-min-confidence 0.55 \
  --soft-risk-threshold 0.75 \
  --soft-edge-scale 0.25 \
  --hyperread-norm sqrt_cap \
  --hyperread-cap 10 \
  --outroot "$OUT" \
  --force
