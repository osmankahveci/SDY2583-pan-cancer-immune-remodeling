# Cross-panel convergence analysis

This analysis tests whether independently generated SDY2583 immune-remodeling scores converge at the subject level across cytometry panels.

## Input

Provide the participant-level ALL10 score matrix through:

```bash
export SDY2583_CROSS_PANEL_MATRIX_FILE="/path/to/SDY2583_integrated_clinical_immune_score_matrix_ALL10_with_CP23.csv"
```

Alternatively, place a compatible integrated score matrix under `SDY2583_INTEGRATED_DIR`; the script prioritizes filenames containing `ALL10`, `with_CP23`, or `ALL_PANELS`.

## Run

```bash
Rscript scripts/40_run_cross_panel_convergence.R
```

## Statistical design

The script calculates all pairwise Spearman correlations among analyzable composite scores and distinguishes cross-panel from within-panel pairs.

- Full-cohort scores are residualized for age, sex, and disease group.
- Cancer-only and healthy-only scores are residualized for age and sex.
- Benjamini–Hochberg FDR is calculated globally and within the cross-panel/within-panel analysis scope.
- A minimum of 50 pairwise complete observations is required by default; this can be changed with `SDY2583_CROSS_PANEL_MIN_N`.
- Prespecified convergence summaries cover CP7–CP24, CP7–CP28, CP24–CP28, CP8–CP25, CP10–CP16, CP10–CP23, CP16–CP23, and CP26–CP28.

## Main outputs

The default output folder is `09_cross_panel_convergence` under `SDY2583_INTEGRATED_DIR`.

- full pairwise correlation table;
- cross-panel-only correlation table;
- panel-pair convergence summary;
- prespecified biological-axis summary;
- robustness table comparing full, cancer-only, and healthy-only directions;
- one automatically selected integrated/core score per panel;
- panel-level heatmap in PNG and PDF;
- principal-score heatmap in PNG and PDF;
- principal-score network in PNG and PDF when qualifying edges are present;
- manifest, text summary, and R session information.

Participant-level matrices and generated outputs remain local and are excluded from Git.
