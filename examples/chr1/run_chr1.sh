#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA="$HERE/data"
OUT="${1:?Usage: bash run_chr1.sh NEW_OUTPUT_DIRECTORY}"
PYTHON="${PHIREFLY_PYTHON:-${VIRTUAL_ENV:-${CONDA_PREFIX:-}}/bin/python}"
if [[ "$PYTHON" == /bin/python ]]; then PYTHON=python; fi
BENCHMARK="${PHIREFLY_BENCHMARK:-0}"
case "$BENCHMARK" in 0|1) ;; *) printf '%s\n' 'PHIREFLY_BENCHMARK must be 0 or 1.' >&2; exit 1 ;; esac
"$PYTHON" -I -c 'import phirefly' || {
  printf '%s\n' "Cannot import Phirefly using: $PYTHON; set PHIREFLY_PYTHON." >&2
  exit 1
}
if [[ "$BENCHMARK" == 1 ]]; then "$PYTHON" -I -c 'import whatshap'; fi
(cd "$DATA" && sha256sum --check ../SHA256SUMS)

"$PYTHON" -I -m phirefly.phaselet \
  --region-key chr1 --region-label chr1 --region chr1:1-248956422 \
  --observations "$DATA/phirefly_reads.observations.npz" \
  --input-vcf-gz "$DATA/shared.snps.vcf.gz" --sample HG002 \
  --tau-snp-phases "$DATA/msf_backbone.snp_phases.tsv" \
  --phaselet-config-name risk0p45_soft0p75_s0p25 \
  --min-abs-support 2 --min-count 2 --min-confidence 0.65 \
  --max-distance-bp 1000000 --algorithm CFC \
  --backend "${PHIREFLY_BACKEND:-cpu-float32}" \
  --batch-size 100 --n-iter 1000 --seed 23 \
  --bridge-mode risk --risk-threshold 0.45 \
  --soft-bridge-mode phaselet --soft-risk-threshold 0.75 --soft-edge-scale 0.25 \
  --hyperread-norm sqrt_cap --hyperread-cap 10 --outroot "$OUT"

RESULT="$OUT/chr1/alg_CFC/phaselet_risk0p45_soft0p75_s0p25"
if [[ "$BENCHMARK" == 1 ]]; then
  "$PYTHON" -I -m phirefly.cli benchmark \
    --input-vcf "$DATA/shared.snps.vcf.gz" --input-sample HG002 \
    --truth-vcf "$DATA/truth.region.bcf" --truth-sample HG002 --pred-sample HG002 \
    --region chr1:1-248956422 \
    --method "Phirefly:$RESULT/phased.component_ps.vcf.gz" \
    --out-tsv "$RESULT/benchmark_metrics.tsv"
fi
