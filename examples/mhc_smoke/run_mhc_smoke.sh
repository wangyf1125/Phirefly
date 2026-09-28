#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
DATA="$ROOT/examples/mhc_smoke/data"
OUT="${1:-$ROOT/examples/mhc_smoke/output}"

PYTHON="${PHIREFLY_PYTHON:-python}"
BENCHMARK="${PHIREFLY_BENCHMARK:-0}"
# Test the installed distribution, never inherited PYTHONPATH or checkout imports.
"$PYTHON" -I -c 'import phirefly' || {
  printf '%s\n' 'Activate the Phirefly environment and install .[qaia] before running smoke.' >&2
  exit 1
}
case "$BENCHMARK" in
  0) ;;
  1) "$PYTHON" -I -c 'import whatshap' || {
       printf '%s\n' 'Benchmark smoke also requires WhatsHap 2.8: install .[qaia,benchmark].' >&2
       exit 1
     } ;;
  *) printf '%s\n' 'PHIREFLY_BENCHMARK must be 0 or 1.' >&2; exit 1 ;;
esac
RESULT="$OUT/mhc/alg_CFC/phaselet_smoke_risk_soft"
if [[ -d "$RESULT" ]] && [[ -n "$(ls -A "$RESULT")" ]]; then
  printf '%s\n' "Smoke output already exists: $RESULT" \
    'Keep existing results; pass a new output directory as the first argument.' >&2
  exit 1
fi

"$PYTHON" -I -m phirefly.phaselet \
  --region-key mhc \
  --region-label MHC \
  --region chr6:28510120-33480577 \
  --observations "$DATA/phirefly_reads.observations.npz" \
  --input-vcf-gz "$DATA/shared.snps.vcf.gz" \
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
  --outroot "$OUT"

if [[ "$BENCHMARK" == 1 ]]; then
  "$PYTHON" -I -m phirefly.cli benchmark \
    --input-vcf "$DATA/shared.snps.vcf.gz" \
    --truth-vcf "$DATA/truth.region.bcf" \
    --input-sample HG002 --truth-sample HG002 --pred-sample HG002 \
    --region chr6:28510120-33480577 \
    --method "Phirefly:$RESULT/phased.component_ps.vcf.gz" \
    --out-tsv "$RESULT/benchmark_metrics.tsv"
  "$PYTHON" -I - "$RESULT/benchmark_metrics.tsv" <<'PY'
import csv
import sys
with open(sys.argv[1]) as handle:
    row = next(csv.DictReader(handle, delimiter='\t'))
if int(row['assessed_pairs']) == 0 or int(row['hamming_denominator']) == 0:
    raise SystemExit('No evaluable truth overlap; benchmark smoke failed.')
print(f"SE {float(row['SE_percent']):.3f}% ({row['switches']}/{row['assessed_pairs']})  "
      f"HE {float(row['HE_percent']):.3f}% ({row['hamming']}/{row['hamming_denominator']})")
PY
fi
