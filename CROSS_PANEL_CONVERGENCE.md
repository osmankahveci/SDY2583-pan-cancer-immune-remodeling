# Cross-panel convergence analysis

This project contains two complementary cross-panel analyses.

## 1. Full-cohort convergence analysis

The original integrated analysis tests whether independently generated
SDY2583 immune-remodeling scores converge at the subject level across cytometry
panels.

Provide the participant-level ALL10 score matrix through:

```bash
export SDY2583_CROSS_PANEL_MATRIX_FILE="/path/to/SDY2583_integrated_clinical_immune_score_matrix_ALL10_with_CP23.csv"
Rscript scripts/40_run_cross_panel_convergence.R
```

The full-cohort analysis residualizes principal scores for age, sex and disease,
calculates pairwise Spearman correlations, applies BH FDR, and evaluates
cancer-only and healthy-only directional robustness. In the locked v1.0.0
results, 38/45 principal-score relationships were FDR-significant.

## 2. v1.1.0 locked held-out architecture

The locked internal validation asks a different question: whether the
cross-panel correlation architecture learned under a training-frozen workflow
is retained in the original held-out VALIDATION cohort.

Run:

```bash
Rscript scripts/44_run_locked_internal_validation.R
```

For each principal score, age/sex/disease residualization coefficients are fit
only in TRAINING and projected into VALIDATION. Pairwise Spearman correlations
are then calculated separately in TRAINING and VALIDATION.

Aggregate reference results:

- 45 principal-score panel pairs tested;
- 40/45 retained the same correlation direction;
- 35/45 were BH-FDR significant in VALIDATION;
- TRAINING-versus-VALIDATION rho concordance: Pearson r = 0.700, p = 8.74e-8;
- Spearman concordance: rho = 0.676, p = 3.54e-7.

This is evidence of internal architectural generalization, not an external
cohort validation.

## Public outputs

Aggregate cross-panel tables are under `results/final/`. Participant-level
residuals, participant identifiers, score matrices and individual PCA
coordinates remain local and are excluded from Git.
