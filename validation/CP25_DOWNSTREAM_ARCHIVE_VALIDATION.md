# CP25 downstream archive validation

Date: 2026-07-20

## Scope

Archived CP25 participant data, 77-feature primary and sensitivity models, nine subject-level composite scores, age-stratified analyses, same-sex age matching, disease-by-age interactions, and the exact 73-outcome threshold-sensitivity family were independently recalculated. This validation exposed incomplete age and threshold target families in the initial reconstructed source; those sections and the panel validator were corrected.

This validates downstream reproducibility from archived participant-level tables. It does not replace local execution against the 850 raw CP25 FCS files.

## Results

| Layer | Archived rows | Result | Maximum absolute numerical difference |
|---|---:|---|---:|
| Primary age/sex-adjusted feature models | 77 | Pass | 1.95 × 10^-13 |
| Binary-sex sensitivity feature models | 77 | Pass | 1.96 × 10^-13 |
| Primary-gate event-QC feature models | 77 | Pass | 1.25 × 10^-13 |
| Subject-level composite scores | 850 subjects × 9 scores | Pass | 1.42 × 10^-14 |
| Primary composite-score models | 9 | Pass | 7.66 × 10^-15 |
| Age-stratified adjusted models | 258 | Pass | 1.38 × 10^-13 |
| Same-sex nearest-age subject matching | 265 five-year pairs; 273 ten-year pairs | Exact pair match | 0 |
| Age/sex-adjusted matched models | 172 | Pass | 5.71 × 10^-13 |
| Disease-by-age interaction models | 86 | Pass | 2.75 × 10^-14 |
| Threshold-set subject-level scores | 2,550 rows × 9 scores | Pass | Machine precision |
| Targeted threshold-sensitivity models | 219 | Pass | 2.36 × 10^-13 |
| Threshold robustness summary | 73 outcomes | Exact classification match | 0 |

## Archive-specific rules confirmed

1. CP25 contains 77 modeled CD4 regulatory/checkpoint features and nine composite scores.
2. Primary models use `outcome ~ disease_group + age_for_model + sex` in 832 model-ready subjects.
3. Binary-sex sensitivity excludes the 15 `Not Specified` sex records; event-QC sensitivity requires `n_cd3_cd4_primary >= 300`.
4. The nine composite scores are constructed from oriented raw-feature z scores; the integrated score is the mean of standardized domain scores. Subject-level values reproduce the archive exactly.
5. All nine composite scores are higher in cancer and FDR-significant. The archived integrated adjusted beta is 0.4392101623.
6. Age sensitivity uses the complete 86-outcome family: 77 raw features plus nine scores. Each of the three age strata therefore contains 86 models, for 258 rows total.
7. Age-stratum FDR is calculated separately within Young <40, Middle 40-59, and Older 60+.
8. Same-sex nearest-age matching processes cancer subjects in ascending sex/age/subject order. For equal age distance, the lower-age control is selected first, followed by subject ID. This exactly reproduces 265 five-year and 273 ten-year pairs.
9. Matched inference remains age/sex adjusted; it is not a matched-pair fixed-effect model.
10. Disease-by-age interaction testing uses all 86 outcomes. Age is standardized separately within each outcome's complete-case dataset before fitting `outcome ~ disease_group * age_z + sex`.
11. The threshold-sensitivity family is an exact archived list of 64 features plus nine scores. It cannot be reconstructed by merely excluding broad modules and taking the first 64 features.
12. Positive-marker thresholds shift by -0.2/+0.2 for permissive/stringent extraction, whereas DUMP-low and IL7RA-low thresholds shift in the opposite direction.
13. Threshold robustness classification is exactly reproduced: 48 outcomes preserve direction and FDR, 23 preserve direction without FDR in every set, and two do not preserve direction.

## Corrections applied to reconstructed CP25 source

- Added the exact archived 64-feature threshold target list to the panel manifest.
- Expanded age sensitivity from a selected 20-outcome subset to the full 86-outcome archive family.
- Replaced subject-ID-first matching ties with lower-age-first tie resolution, reproducing the exact 265/273 archived pairs.
- Changed age-stratum FDR from a pooled adjustment to separate within-stratum adjustments.
- Changed interaction age standardization from a single panel-wide z score to outcome-specific complete-case standardization.
- Rewrote threshold outputs with the archived 73-outcome family and archive-compatible robustness labels.
- Expanded the CP25 validator to cover primary, sensitivity, score, age-stratified, exact-pair, matched, interaction, and threshold layers.

## Remaining validation

- Execute Steps 1–2 on the local 850 CP25 FCS files and compare inventory, compensation/transformation, gating, and feature tables with the archive.
- Run the corrected complete CP25 pipeline with `SDY2583_CP25_REFERENCE_DIR` configured and retain its validation report and session information.
- Validate optional clinical-annotation outputs if the required clinical source tables are available.

CP25 is therefore **downstream archive-validated after age/threshold reconstruction corrections, but raw-FCS and optional clinical-extension validation remain pending**.
