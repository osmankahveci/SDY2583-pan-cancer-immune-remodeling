# Final integrated analysis release

## Scope

This release completes the manuscript-level analyses added after the ten-panel
downstream release. It is designed to be run on local participant-level files
that are excluded from Git. Public outputs contain aggregate statistics only.

## Executable workflows

### `scripts/41_run_integrated_systems_analysis.py`

- ten prespecified principal integrated scores;
- age/sex/disease residualization;
- pairwise-complete Spearman and Pearson convergence analyses;
- cancer-only and healthy-only direction checks;
- robust immunotype models in the full, training, validation, cancer-only, and
  treatment-adjusted cancer-only cohorts;
- adjusted immunotype marginal means;
- disease-adjusted PCA, immunotype centroids, and training-validation loading
  and centroid replication;
- k-means, Ward, and Gaussian-mixture clustering diagnostics;
- Hopkins, silhouette, Calinski-Harabasz, Davies-Bouldin, ARI, NMI, AIC, and BIC
  summaries.

### `scripts/42_run_frequency_matching_table1.py`

- eight prespecified raw percentage outcomes;
- raw-percentage OLS with HC3 covariance;
- empirical-logit OLS with HC3 covariance;
- beta regression with logit link and HC3 covariance;
- same-sex nearest-age matching without replacement using 5- and 10-year
  calipers;
- pre/post standardized mean differences and age-gap distributions;
- matched-pair-aware standard errors and FDR correction;
- demographic/clinical Table 1 generation;
- panel-level FCS availability summary.

### `scripts/43_run_bootstrap_component_stability.py`

- 2,000 disease-stratified resamples by default;
- complete repetition of covariate residualization, scaling, and PCA;
- component matching and sign alignment;
- loading cosine congruence and two-component canonical correlations;
- explained-variance and immunotype-centroid stability;
- age/sex-adjusted cancer effects for all 66 scores with HC3 covariance;
- direction-stability, bootstrap interval, and FDR-retention probabilities.

### `R/figures/figure7_integrated_systems_v6.R`

Locked v6 code for the integrated manuscript figure:

- Panel A: adjusted immunotype profiles;
- Panel B: selected biologically coherent cross-panel convergence network;
- Panel C: disease-adjusted PCA immunotype centroids.

The script generates 600-dpi PNG and vector PDF outputs locally. Generated
figures are not committed automatically.

## Locked manuscript findings

- 38/45 principal-score Spearman relationships were FDR-significant.
- 39/45 were FDR-significant in Pearson sensitivity analysis.
- All 38 significant Spearman relationships retained direction in the healthy
  and cancer strata.
- All ten principal scores showed significant immunotype effects in the full,
  training, validation, cancer-only, and treatment-adjusted cancer-only models.
- Disease-adjusted PCA: PC1 = 30.5%; PC2 = 16.6%.
- Integrated-score space was structured but not sufficiently separated to
  support new discrete subtype claims.
- All 8/8 bounded outcomes retained direction and FDR significance across three
  model scales and both matched analyses.
- Matching: 265 pairs under the 5-year caliper and 273 pairs under the 10-year
  caliper.
- Bootstrap: 62/66 scores retained direction in at least 95% of resamples; all
  ten principal intervals excluded zero.

## Interpretation boundary

The release supports internal reproducibility, convergent validity, and
biological resolution of the original SDY2583 immunotypes. It does not provide
an independent external cohort, functional validation, tumor-tissue
correspondence, treatment-response prediction, or prognostic validation.

Cancer-subgroup findings, including the breast-cancer CP10 pattern, remain
hypothesis-generating. Treatment annotations do not justify claims of treatment
independence because stage, burden, timing, prior exposure, and treatment line
were incompletely observed.

## Data boundary

Do not commit:

- raw FCS files;
- local ImmPort extracts;
- participant-level feature or score matrices;
- participant identifiers or matched-pair IDs;
- participant-level PCA coordinates;
- local path files or credentials.
