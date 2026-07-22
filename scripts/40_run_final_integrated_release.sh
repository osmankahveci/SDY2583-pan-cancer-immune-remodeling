#!/usr/bin/env bash
set -euo pipefail

: "${SDY2583_ALL10_MATRIX:?Set SDY2583_ALL10_MATRIX}"
: "${SDY2583_CP24_FEATURES:?Set SDY2583_CP24_FEATURES}"
: "${SDY2583_CP8_FEATURES:?Set SDY2583_CP8_FEATURES}"
: "${SDY2583_CP10_FEATURES:?Set SDY2583_CP10_FEATURES}"
: "${SDY2583_CP16_FEATURES:?Set SDY2583_CP16_FEATURES}"
: "${SDY2583_CP22_FEATURES:?Set SDY2583_CP22_FEATURES}"
: "${SDY2583_CP23_FEATURES:?Set SDY2583_CP23_FEATURES}"
: "${SDY2583_CP26_FEATURES:?Set SDY2583_CP26_FEATURES}"
: "${SDY2583_CP28_FEATURES:?Set SDY2583_CP28_FEATURES}"

RESULTS_DIR="${SDY2583_FINAL_RESULTS_DIR:-results/local-final}"
BOOTSTRAPS="${SDY2583_BOOTSTRAP_RESAMPLES:-2000}"
BOOTSTRAP_SEED="${SDY2583_BOOTSTRAP_SEED:-2583}"
mkdir -p "$RESULTS_DIR"

python scripts/41_run_integrated_systems_analysis.py \
  --matrix "$SDY2583_ALL10_MATRIX" \
  --output-dir "$RESULTS_DIR"

python scripts/42_run_frequency_matching_table1.py \
  --matrix "$SDY2583_ALL10_MATRIX" \
  --cp24 "$SDY2583_CP24_FEATURES" \
  --cp8 "$SDY2583_CP8_FEATURES" \
  --cp10 "$SDY2583_CP10_FEATURES" \
  --cp16 "$SDY2583_CP16_FEATURES" \
  --cp22 "$SDY2583_CP22_FEATURES" \
  --cp23 "$SDY2583_CP23_FEATURES" \
  --cp26 "$SDY2583_CP26_FEATURES" \
  --cp28 "$SDY2583_CP28_FEATURES" \
  --output-dir "$RESULTS_DIR"

python scripts/43_run_bootstrap_component_stability.py \
  --matrix "$SDY2583_ALL10_MATRIX" \
  --output-dir "$RESULTS_DIR" \
  --resamples "$BOOTSTRAPS" \
  --seed "$BOOTSTRAP_SEED"

if command -v Rscript >/dev/null 2>&1; then
  Rscript R/figures/figure7_integrated_systems_v6.R
else
  printf '%s\n' 'Rscript not found; skipped Figure 7 regeneration.' >&2
fi

printf '%s\n' "Final integrated analysis completed. Aggregate outputs: $RESULTS_DIR"
