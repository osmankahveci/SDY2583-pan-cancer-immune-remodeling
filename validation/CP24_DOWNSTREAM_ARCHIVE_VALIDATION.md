# CP24 downstream archive validation

Date: 2026-07-20

## Scope

The official CP24 850-subject clean-age feature table, 55-row outcome-statistics table, three-set event-count sensitivity tables, participant-level score dataset, seven primary composite-score results, and three-set composite sensitivity results were independently recalculated. This validation exposed substantive gaps in the initial reconstructed CP24 source; the manifest, extraction, modeling, scoring, figure, and validation code were corrected accordingly.

This establishes downstream numerical reproducibility from the archived participant-level tables. It does not replace local execution against the 850 raw CP24 FCS files.

## Results

| Layer | Archived rows | Result | Maximum absolute numerical difference |
|---|---:|---|---:|
| Primary clean-age outcome models | 55 rows / 54 unique outcomes | Pass | 9.42 × 10^-14 across model quantities |
| Module and global FDR values | 55 | Pass | 1.69 × 10^-14 |
| Event-count sensitivity models | 165 | Pass | 3.34 × 10^-13 across model quantities |
| Event-count sensitivity FDR values | 165 | Pass | 2.09 × 10^-13 |
| Subject-level hierarchical composite scores | 850 subjects × 7 scores | Pass | 1.69 × 10^-14 |
| Primary composite-score models | 7 | Pass | 4.33 × 10^-15 |
| Composite event-count sensitivity models | 21 | Pass | 5.55 × 10^-14 across model quantities |
| Composite sensitivity FDR values | 21 | Pass | 9.44 × 10^-15 |

All numerical differences are at floating-point/machine-precision scale.

## Archive-specific rules confirmed

1. CP24 contains 850 subjects: 408 healthy controls and 442 cancer patients.
2. Clean age is defined by retaining ages from 18 through 100 years; 832 subjects meet this age criterion.
3. Fifteen subjects have missing sex in the archived feature table. Four of these have a valid clean age.
4. The official table reports `n_model = 832`, defined by nonmissing outcome and clean age. The actual `lm()` coefficient is estimated after complete-case removal of missing sex and therefore uses 828 records. Treating missing sex as a third category does not reproduce the archive.
5. The official model family contains 55 module-outcome rows and 54 unique outcomes because `median_PD1_in_CD3CD8` belongs to both `marker_medians` and `pd1_module`.
6. The official module sizes are 4 composition, 9 marker-median, 8 CD45RA/CD62L quadrant, 8 CD27-refined TEMRA, 7 extended-differentiation, and 19 PD-1 rows.
7. Sixteen of 55 rows are significant after global FDR and 17 after within-module FDR.
8. Official event-count sensitivity uses all subjects, `n_cd3_cd8_pos >= 500`, and `n_cd3_cd8_pos >= 1000`. The extraction-set sizes are 850, 849, and 844; the corresponding reported clean-age model sizes are 832, 831, and 826.
9. All 55 outcome directions are preserved across the three event-count sets.
10. CP24 has seven hierarchical scores. Five domain scores are means of oriented raw-feature z scores. Core differentiation is the arithmetic mean of the marker-loss, quadrant-shift, CD27-refined terminal, and extended-effector scores. Integrated remodeling is the arithmetic mean of core differentiation and the PD-1-associated terminal score.
11. The archived adjusted betas are 0.195445486256139 for marker loss, 0.17236847025103 for quadrant shift, 0.178728624529361 for CD27-refined terminal differentiation, 0.125120779018976 for extended effector remodeling, 0.0597846945686508 for the PD-1-associated terminal score, 0.167915840013875 for core differentiation, and 0.113850267291264 for integrated remodeling.
12. Six of seven scores are FDR-significant. The PD-1-associated terminal score is not significant and must not be presented as a robust primary finding.
13. The same model-reporting convention and event-count sets are used for the 21 composite-score sensitivity models.
14. Same-sex matching, disease-by-age interaction, and ±0.2 threshold-shift analyses in the repository are explicitly reconstructed extensions; no direct official CP24 archive tables were available for those extensions in this validation pass.

## Corrections applied to reconstructed CP24 source

- Expanded the module map from 46 to the official 55 rows.
- Added the missing PD-1 state-specific, PD-1-high, and PD-1/co-differentiation endpoints.
- Corrected `pct_pd1_pos_within_temra` to the archived `pct_pd1_pos_within_temra_like` endpoint.
- Replaced the generic module-map join with row-wise panel-specific modeling so the duplicated median PD-1 endpoint is represented exactly twice rather than being accidentally multiplied by a many-to-many join.
- Reproduced the official distinction between descriptive counts, reported clean-age `n_model`, and complete-case `lm()` records.
- Replaced the initial generic composite construction with the exact archived hierarchical score formulas.
- Updated age and threshold extensions to use all seven official scores.
- Updated quantitative figures to use the official score family and FDR column names.
- Replaced the permissive rounded validator with fixed subject, feature, model, score, event-sensitivity, significance, and optional archive-table checks.

## Repository validation gate

`scripts/26_validate_cp24_reconstruction.R` now checks:

- 850 successful FCS-derived rows;
- archived event-count and clean-age benchmarks;
- the exact 55-row/54-unique-outcome family;
- global and module significance counts;
- all 165 event-sensitivity rows and reported model sizes;
- all seven hierarchical score columns, primary score results, and 21 score-sensitivity rows;
- exact core and integrated score betas;
- optional participant-, outcome-, and score-table comparisons when `SDY2583_CP24_REFERENCE_DIR` is configured.

## Remaining validation

- Execute inventory, compensation, fixed-logicle transformation, gating, and feature extraction on the local 850 CP24 FCS files.
- Compare the regenerated feature table and QC summaries with the official archive.
- Run the complete corrected CP24 pipeline with `SDY2583_CP24_REFERENCE_DIR` configured and retain its validation report and session information.
- Evaluate the reconstructed age-matching, interaction, and ±0.2 threshold extensions separately without describing them as recovered official analyses.

CP24 is therefore **downstream archive-validated after substantive reconstruction corrections, but raw-FCS validation and reconstructed-extension validation remain pending**.
