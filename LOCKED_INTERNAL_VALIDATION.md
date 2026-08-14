# Locked internal held-out validation

## Purpose

Version 1.1.0 adds a leakage-reduced **locked internal held-out validation** of the 66-score SDY2583 immune-remodeling framework. This layer is distinct from the repository's earlier **downstream archive validation**, which tests reproducibility of archived panel outputs.

The original SDY2583 cohort labels are retained:

- TRAINING: 503 participants
- VALIDATION: 347 participants
- total: 850 participants

No new random split was generated.

## Locked design

Composite-score membership is fixed from the existing 66-score project registry. For each raw or nested component:

1. the cancer-versus-healthy coefficient is estimated only in TRAINING with an age- and sex-adjusted linear model using HC3 covariance;
2. the component direction is set from the TRAINING cancer coefficient;
3. the mean and standard deviation are estimated only in TRAINING;
4. the frozen direction, mean and standard deviation are applied unchanged to VALIDATION;
5. hierarchical integrated scores are recursively rebuilt from training-frozen component scores;
6. validation data are never used to select direction, centering, scaling or PCA loadings.

The primary analysis is the ten prespecified principal panel scores. The 66-score catalogue, frozen PCA and cross-panel correlation architecture are secondary validation layers.

## Primary result

Nine of ten principal scores retained the cancer-associated direction and remained significant after Benjamini-Hochberg correction across the ten validation tests.

| Panel | TRAINING beta | VALIDATION beta | VALIDATION 95% CI | VALIDATION BH-q | Status |
|---|---:|---:|---:|---:|---|
| CP7 | 0.438 | 0.244 | 0.122 to 0.366 | 1.30e-4 | replicated |
| CP8 | 0.586 | 0.440 | 0.324 to 0.556 | 3.80e-13 | replicated |
| CP10 | 0.620 | 0.294 | 0.160 to 0.427 | 2.76e-5 | replicated |
| CP16 | 0.395 | 0.587 | 0.439 to 0.735 | 8.39e-14 | replicated |
| CP22 | 0.474 | 0.225 | 0.164 to 0.286 | 8.77e-13 | replicated |
| CP23 | 0.287 | 0.189 | 0.043 to 0.335 | 0.0125 | replicated |
| CP24 | 0.244 | -0.022 | -0.224 to 0.180 | 0.830 | **not replicated** |
| CP25 | 0.413 | 0.654 | 0.485 to 0.823 | 1.47e-13 | replicated |
| CP26 | 0.358 | 0.290 | 0.168 to 0.412 | 6.63e-6 | replicated |
| CP28 | 0.356 | 0.184 | 0.084 to 0.283 | 3.75e-4 | replicated |

CP24 is retained transparently as a non-replicated principal axis; it is not re-tuned on the validation cohort.

Across all 66 locked scores, 60/66 retained direction and 44/66 were validation-significant after 66-test BH correction.

## Frozen PCA

The ten locked principal scores were standardized with TRAINING means and standard deviations. PCA was fitted only in TRAINING, the PC1 sign was fixed in TRAINING, and VALIDATION participants were projected without refitting.

- TRAINING complete cases: 485
- VALIDATION complete cases: 344
- PC1 variance explained in TRAINING: 39.43%
- TRAINING cancer beta: 2.148 (95% CI 1.815 to 2.481)
- VALIDATION cancer beta: 1.569 (95% CI 1.232 to 1.905), p = 6.01e-20

This supports generalization of a global multidimensional remodeling axis without validation-cohort PCA refitting.

## Cross-panel architecture

Age, sex and disease effects were fitted in TRAINING and projected into VALIDATION before pairwise Spearman analysis of the ten principal scores.

- 40/45 panel-pair correlations retained direction;
- 35/45 were validation-significant after BH correction;
- TRAINING-versus-VALIDATION correlation-effect concordance: Pearson r = 0.700, p = 8.74e-8;
- Spearman concordance: rho = 0.676, p = 3.54e-7.

## Interpretation boundary

This analysis should be described as **locked internal held-out validation**, **training-frozen internal validation**, or a **prespecified held-out cohort analysis**.

It is **not external validation**. Score membership itself was established during the earlier full-cohort project. The v1.1.0 procedure removes validation-cohort leakage from component orientation, centering/scaling, higher-level score rebuilding and PCA fitting, but it does not retroactively convert original component membership into de novo discovery-only feature selection.

## Reproducibility and privacy

Run:

```bash
Rscript scripts/44_run_locked_internal_validation.R
```

Provide either a common local input root:

```bash
export SDY2583_LOCKED_INPUT_ROOT="/path/to/local/SDY2583/analysis/archive"
```

or explicit file variables such as `SDY2583_LOCKED_MATRIX_FILE` and `SDY2583_LOCKED_CP7_FILE` through `SDY2583_LOCKED_CP28_FILE`.

The public implementation writes aggregate outputs only. Participant identifiers, participant-level locked scores, individual PCA coordinates, raw FCS files and serialized participant objects are not committed.

## Public reference outputs

Aggregate reference results are under `results/final/` with the `locked_validation_` prefix. Exact score membership is under `config/locked_validation/`.
