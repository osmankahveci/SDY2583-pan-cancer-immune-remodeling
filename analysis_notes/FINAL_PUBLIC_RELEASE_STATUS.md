# Final public release status

## Completed scope

The `cross-panel-convergence` release branch now contains the complete final
analysis layer used in the manuscript:

- adjusted cross-panel Spearman and Pearson convergence;
- healthy-only and cancer-only direction checks;
- full, training, validation, cancer-only, and treatment-adjusted immunotype
  models;
- disease-adjusted PCA and training-validation replication;
- k-means, Ward, and Gaussian-mixture clustering diagnostics;
- bounded-percentage linear-HC3, empirical-logit, and beta-regression models;
- same-sex age matching, balance diagnostics, and matched-pair inference;
- cohort Table 1 and panel availability summaries;
- consolidated panel threshold/event-count robustness;
- 2,000-resample bootstrap stability for all 66 scores and the ten principal
  scores;
- locked Figure 7 v6 source;
- curated aggregate-only result tables.

## Numerical checks

The modular final Python release reproduced the locked manuscript results on the
local analysis files:

- 38/45 FDR-significant adjusted Spearman relationships;
- 39/45 FDR-significant Pearson relationships;
- all 38 significant Spearman relationships directionally preserved in both
  disease strata;
- PC1 = 30.5% and PC2 = 16.6%;
- training-validation loading similarities = 0.953 and 0.884;
- projected centroid correlations = 0.984 and 0.977;
- maximum non-GMM silhouette = 0.207;
- 8/8 bounded outcomes concordant and FDR-significant across model scales;
- 265 and 273 matched pairs with 8/8 paired outcomes FDR-significant;
- 2,000 bootstrap resamples: 62/66 scores directionally stable in at least 95%
  of resamples and 10/10 principal intervals excluding zero.

## Public-data boundary

No raw FCS files, participant-level metadata, subject identifiers,
participant-level scores, matched-pair identifiers, or individual PCA
coordinates are included. Source data remain available from ImmPort under
SDY2583.

## Validation boundary

All ten panel downstream analyses are archive-validated. Fresh raw-FCS
end-to-end re-execution was not performed. The final integrated analyses are
internally reproduced but do not constitute an independent external cohort or
functional validation.
