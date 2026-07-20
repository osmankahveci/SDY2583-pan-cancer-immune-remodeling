# CP22 downstream archive validation

Date: 2026-07-20

## Scope

Archived CP22 post-extraction participant data, feature-model outputs, composite definitions and scores, same-sex age matching, disease-by-age interactions, and threshold-sensitivity tables were independently recalculated. This validates the downstream numerical core of recovered `SAFE` Steps 3B–5. It does not replace local execution of Steps 1–2 against the 850 raw CP22 FCS files.

## Results

| Layer | Archived rows | Result | Maximum absolute numerical difference |
|---|---:|---|---:|
| Primary age/sex-adjusted feature models | 48 | Pass | 2.74 × 10^-13 |
| Binary-sex sensitivity feature models | 48 | Pass | 1.65 × 10^-13 |
| B-cell event-QC feature models | 48 | Pass | 1.74 × 10^-13 |
| Standard-dump sensitivity feature models | 48 | Pass | 1.74 × 10^-13 |
| Subject-level composite scores | 850 subjects × 7 scores | Pass | 3.11 × 10^-15 |
| Primary and sensitivity composite models | 28 | Pass | 6.77 × 10^-15 |
| Same-sex nearest-age matching | 265 five-year pairs; 273 ten-year pairs | Archived counts and age-difference summaries reproduced | machine precision |
| Age/sex-adjusted matched composite models | 14 | Pass | 5.11 × 10^-15 |
| Disease-by-age interaction models | 7 | Pass | 7.55 × 10^-15 |
| Threshold-set subject-level scores | 2,550 rows × 7 scores | Pass | 5.33 × 10^-15 |
| Threshold-sensitivity models | 165 | Pass | 2.74 × 10^-13 |
| Threshold robustness summary | 55 variables | Exact classification match | 0 |

All numerical differences are at floating-point/machine-precision scale.

## Archive-specific rules confirmed

1. CP22 contains 850 FCS files/subjects; all 850 archived files were successfully processed and 832 subjects were model-ready after age/sex QC.
2. The primary analysis tests 48 B-cell/humoral features under `feature ~ disease_group + age_for_model + sex`.
3. Prespecified feature sensitivities are binary sex, `n_cd19_b >= 300`, and restriction to the standard `Viability_CD3_CD7_CD13` dump-marker configuration.
4. Forty-four of 48 feature directions are preserved across all four model sets; global FDR remains below 0.05 for 31 features in every set.
5. Seven phenotype-oriented composite scores are standardized across all 850 post-extraction subjects.
6. The integrated humoral B-cell remodeling score has an archived adjusted beta of 0.4072877744339053 and is higher in cancer.
7. Same-sex age matching processes cancer subjects in ascending-age order and selects the nearest unused same-sex healthy control, with lower-age tie resolution. Archived 265/273 matched-pair counts and age-difference summaries were reproduced.
8. Matched inference uses `score ~ disease_group + age_for_model + sex`; it is not a matched-pair fixed-effect model.
9. Disease-by-age interaction testing is restricted to the seven composite scores and uses `score ~ disease_group * age_z + sex`.
10. Threshold sensitivity tests 55 variables: 48 feature-dictionary variables and seven composites. BH FDR is calculated separately within each threshold set, both globally and within biological modules.
11. Threshold robustness classification is exactly reproduced: 39 variables preserve direction and global FDR, six preserve direction without FDR in all sets, and ten do not preserve direction.
12. `CP22_Ig_isotype_repatterning_score` is a directed B-cell immunoglobulin-isotype architecture/repatterning composite. It is not a total IgG- or total IgA-abundance score, and it must not be described as a global increase or decrease in IgG+ or IgA+ B cells.

## Repository changes

- Added `scripts/29_validate_cp22_recovered.R`, a CP22-specific fixed-benchmark and optional archive-table validation gate.
- Updated `scripts/17_run_cp22_recovered.R` to invoke the panel-specific validator.
- Recovered Steps 3B–5 were retained unchanged because their downstream algorithms reproduce the archive exactly.

## Remaining validation

- Execute the recovered Steps 1–2 on the local 850 CP22 FCS files.
- Compare regenerated inventory, marker/channel variants, compensation/transformation, gating, and feature tables with the archived reference tree.
- Run the complete CP22 pipeline with `SDY2583_CP22_REFERENCE_DIR` configured and retain the validation report and session information.
- Validate optional clinical-annotation extensions where the necessary clinical source tables are available.

CP22 is therefore **downstream archive-validated but raw-FCS and optional clinical-extension validation pending**.
