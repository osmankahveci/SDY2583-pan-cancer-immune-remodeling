# CP28 downstream archive validation

Date: 2026-07-20

## Scope

Archived CP28 participant data, the 65-row T/NK-interface feature-model family, nine subject-level composite scores, module- and score-specific event-QC sensitivities, age-stratified analyses, same-sex age matching, disease-by-age interactions, and the exact 80-outcome threshold-sensitivity family were independently recalculated. This validation exposed substantive structural gaps in the initial reconstructed source; the manifest, feature-model, event-QC, score, age-sensitivity, threshold, and validation code were corrected accordingly.

This validates downstream numerical reproducibility from archived participant-level tables. It does not replace local execution against the 850 raw CP28 FCS files.

## Results

| Layer | Archived rows | Result | Maximum absolute numerical difference |
|---|---:|---|---:|
| Primary age/sex-adjusted module-feature models | 65 rows / 62 unique features | Pass | 1.54 × 10^-13 |
| Binary-sex sensitivity module-feature models | 65 | Pass | 9.68 × 10^-14 |
| Module-specific event-QC feature models | 65 | Pass | 1.28 × 10^-13 |
| Subject-level composite scores | 850 subjects × 9 scores | Pass | 1.02 × 10^-14 |
| Primary composite-score models | 9 | Pass | 1.02 × 10^-14 |
| Score-specific event-QC models | 9 | Pass | 1.13 × 10^-14 |
| Age-stratified adjusted models | 156 | Pass | 9.14 × 10^-14 |
| Same-sex nearest-age subject matching | 265 five-year pairs; 273 ten-year pairs | Exact pair match | 0 |
| Age/sex-adjusted matched models | 104 | Pass | 5.39 × 10^-13 |
| Disease-by-age interaction models | 52 | Pass | 2.00 × 10^-14 |
| Threshold-set subject-level scores | 2,550 rows × 9 scores | Pass | Machine precision |
| Targeted threshold-sensitivity models | 240 | Pass | 1.84 × 10^-13 |
| Threshold robustness summary | 80 outcomes | Exact classification match | 0 |

All numerical differences are at floating-point/machine-precision scale.

## Archive-specific rules confirmed

1. CP28 contains 850 subjects and 832 valid-age/model-ready subjects.
2. The primary analysis contains 65 module-feature rows but 62 unique features. Three CD3-negative/CD56-positive NK-like composition endpoints belong to both the global T/NK composition module and the NK-like module.
3. The archived module sizes are 13 composition, eight CD8 differentiation, 12 CD8 effector/NK-like, 26 CD3 innate-like T-cell, and six CD3-negative/CD56-positive NK-like rows.
4. Modeling duplicated module membership by joining a unique-feature result table back to a nonunique module map creates 71 rows and is incorrect. Each module-feature row must be modeled and FDR-adjusted in its own module context.
5. Primary and binary-sex feature inference use `outcome ~ disease_group + age_for_model + sex`.
6. Event-QC sensitivity uses a 300-event threshold, but the gate denominator is module-specific: no additional filter for composition, `n_cd3_cd8_pos` for both CD8 modules, `n_cd3_pos_dump_low` for the CD3 innate-like module, and `n_cd3neg_cd56pos_nk_like` for the NK-like module.
7. The corresponding archived event-QC model sizes are 832, 823, 823, 829, and 810, respectively.
8. CP28 contains eight domain scores plus the integrated T/NK-interface remodeling score. Scores are constructed across all 850 post-extraction subjects and modeled in the 832 valid-age subjects.
9. The archived integrated adjusted beta is 0.256754715641472. The CD8-enrichment-within-T score has an adjusted beta of approximately 0.4793; CP28 remains a T/NK-interface panel rather than a pure NK panel.
10. Score event-QC is also domain-specific. The integrated score requires all three relevant gates (`n_cd3_pos_dump_low`, `n_cd3_cd8_pos`, and `n_cd3neg_cd56pos_nk_like`) to have at least 300 events, yielding 808 model records.
11. Age sensitivity uses an exact 52-outcome family: nine scores and 43 selected raw features. Each of three age strata contains 52 models, for 156 rows total.
12. Age-stratum FDR is calculated separately within Young <40, Middle 40-59, and Older 60+.
13. Same-sex nearest-age matching processes cancer subjects in ascending sex/age/subject order. For equal age distance, the lower-age healthy control is selected first, followed by subject ID. This exactly reproduces 265 five-year and 273 ten-year pairs.
14. Matched inference remains age/sex adjusted and uses the complete 52-outcome family, producing 104 rows.
15. Disease-by-age interaction testing uses the same 52 outcomes. Age is standardized separately within each outcome's complete-case dataset before fitting `outcome ~ disease_group * age_z + sex`.
16. Threshold sensitivity uses an exact 80-outcome family: 71 raw features and nine scores. Threshold-set scores are standardized independently within each all-850 extraction.
17. Positive-marker thresholds shift by -0.2/+0.2 for permissive/stringent extraction, while the dump-low threshold shifts in the opposite direction.
18. Threshold robustness classification is exactly reproduced: 35 outcomes preserve direction and FDR, 34 preserve direction without FDR in every set, and 11 do not preserve direction.

## Corrections applied to reconstructed CP28 source

- Added exact 52-outcome age and 80-outcome threshold target lists to the panel manifest.
- Added module-specific feature event-QC and score-specific event-QC maps with the archived 300-event threshold.
- Replaced the generic many-to-many module join with row-wise modeling of the exact 65 module-feature rows.
- Replaced the single `n_cd3_cd8_pos >= 1000` sensitivity with the archived module-specific 300-event rules.
- Added the integrated score's three-gate event-QC intersection.
- Expanded age sensitivity from a small selected subset to the full 52-outcome family.
- Replaced subject-ID-first matching ties with lower-age-first tie resolution, reproducing the exact 265/273 archived pairs.
- Changed age-stratum FDR to separate within-stratum adjustments and interaction age standardization to outcome-specific complete-case standardization.
- Expanded threshold sensitivity from the initial reconstructed family to the exact 80 archived outcomes.
- Replaced the permissive validator with fixed feature, score, event-QC, age, exact-pair, matched, interaction, threshold, and optional archive-table checks.

## Repository validation gate

`scripts/24_validate_cp28_reconstruction.R` now checks:

- the 850-file feature-extraction and event-count benchmarks;
- the exact 65-row/62-unique-feature model family;
- all module-specific event-QC model sizes;
- all nine subject-level and primary score outputs, including the integrated beta;
- score-specific event-QC and the integrated 808-record intersection;
- the full 156 age-stratified, 104 matched, and 52 interaction rows;
- exact 265/273 subject matching;
- all 2,550 threshold score rows, 240 models, 80 outcomes, and 35/34/11 robustness classes;
- optional archive-table comparisons when `SDY2583_CP28_REFERENCE_DIR` is configured.

## Remaining validation

- Execute inventory, compensation, fixed-logicle transformation, gating, and feature extraction on the local 850 CP28 FCS files.
- Compare regenerated feature and QC tables with the archived reference tree.
- Run the corrected complete CP28 pipeline with `SDY2583_CP28_REFERENCE_DIR` configured and retain its validation report and session information.
- Validate optional clinical-annotation outputs if the required clinical source tables are available.

CP28 is therefore **downstream archive-validated after structural model/event-QC/age/threshold corrections, but raw-FCS and optional clinical-extension validation remain pending**.
