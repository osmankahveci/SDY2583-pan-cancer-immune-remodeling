# CP16 downstream archive validation

Date: 2026-07-20

## Scope

Archived CP16 post-extraction participant data, feature-model outputs, composite definitions and scores, same-sex age matching, disease-by-age interactions, and threshold-sensitivity tables were independently recalculated. This validates the downstream numerical core of recovered `SAFE` Steps 3B–5. It does not replace local execution of Steps 1–2 against the 847 raw CP16 FCS files.

## Results

| Layer | Archived rows | Result | Maximum absolute numerical difference |
|---|---:|---|---:|
| Primary age/sex-adjusted feature models | 54 | Pass | 6.82 × 10^-13 |
| Binary-sex sensitivity feature models | 54 | Pass | 3.37 × 10^-13 |
| CD45 dump-low event-QC feature models | 54 | Pass | 8.81 × 10^-13 |
| HLA-DR APC-core event-QC feature models | 54 | Pass | 4.12 × 10^-13 |
| Technical mismatch-exclusion feature models | 54 | Pass | 4.83 × 10^-13 |
| Subject-level composite scores | 847 subjects × 6 scores | Pass | 2.22 × 10^-15 |
| Primary and sensitivity composite models | 30 | Pass | 9.99 × 10^-15 |
| Same-sex nearest-age matching | 264 five-year pairs; 272 ten-year pairs | Exact pair match | 0 |
| Age/sex-adjusted matched composite models | 12 | Pass | 3.89 × 10^-15 |
| Disease-by-age interaction models | 6 | Pass | 3.66 × 10^-15 |
| Threshold-set subject-level scores | 2,541 rows × 6 scores | Pass | 2.31 × 10^-14 |
| Threshold-sensitivity models | 180 | Pass | 6.82 × 10^-13 |
| Threshold robustness summary | 60 variables | Exact classification match | 0 |

All numerical differences are at floating-point/machine-precision scale.

## Archive-specific rules confirmed

1. The CP16 cohort contains 847 FCS files/subjects, not 850. All 847 archived FCS files were successfully processed; 829 subjects were model-ready after age/sex QC.
2. The primary analysis tests 54 features under `feature ~ disease_group + age_for_model + sex`.
3. Prespecified feature sensitivities are binary sex, `n_cd45_dump_low >= 300`, `n_hladr_apc_core >= 300`, and exclusion of channel/marker mismatch files.
4. Fifty-one of 54 feature directions are preserved across all five model sets; global FDR remains below 0.05 for 36 features in every set.
5. Six phenotype-oriented composite scores are standardized over all 847 post-extraction subjects. A score requires at least half of its defined components to be available.
6. The integrated APC/DC-like myeloid remodeling score has an archived adjusted beta of 0.447860154036832 and is higher in cancer.
7. Same-sex age matching processes sex/age/subject-ID ordered cancer records and selects the nearest unused control, with lower-age and subject-ID tie resolution. The exact 264/272 archived pairs were reproduced.
8. Matched inference uses `score ~ disease_group + age_for_model + sex`; it is not a matched-pair fixed-effect model.
9. Disease-by-age interaction testing is restricted to the six composites and uses `score ~ disease_group * age_z + sex`.
10. Threshold sensitivity tests 60 variables: 54 feature-dictionary variables and six composites. BH FDR is calculated separately by threshold set and separately for feature versus composite families.
11. Threshold robustness classification is exactly reproduced: 33 variables preserve direction and global FDR, 14 preserve direction without FDR in all sets, and 13 do not preserve direction.
12. CP16 terminology remains phenotype-based: APC-like myeloid, monocyte-like, cDC1-like, cDC2-like, and pDC-like. The analysis does not establish definitive cell identities.

## Repository changes

- Added `scripts/28_validate_cp16_recovered.R`, a CP16-specific fixed-benchmark and optional archive-table validation gate.
- Updated `scripts/16_run_cp16_recovered.R` to invoke the panel-specific validator.
- Recovered Steps 3B–5 were retained unchanged because their downstream algorithms reproduce the archive exactly.

## Remaining validation

- Execute the recovered Steps 1–2 on the local 847 CP16 FCS files.
- Compare regenerated inventory, mismatch, compensation/transformation, gating, and feature tables with the archived reference tree.
- Run the complete CP16 pipeline with `SDY2583_CP16_REFERENCE_DIR` configured and retain the validation report and session information.
- Validate optional clinical-annotation extensions where the necessary clinical source tables are available.

CP16 is therefore **downstream archive-validated but raw-FCS and optional clinical-extension validation pending**.
