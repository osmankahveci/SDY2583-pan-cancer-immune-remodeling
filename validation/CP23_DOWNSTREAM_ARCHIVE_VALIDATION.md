# CP23 downstream archive validation

Date: 2026-07-20

## Scope

Archived CP23 post-extraction participant data, feature-model outputs, composite definitions and scores, same-sex age matching, disease-by-age interactions, and targeted threshold-sensitivity tables were independently recalculated. This validates the downstream numerical core of recovered `SAFE` Steps 3B–5. It does not replace local execution of Steps 1–2 against the 850 raw CP23 FCS files.

## Results

| Layer | Archived rows | Result | Maximum absolute numerical difference |
|---|---:|---|---:|
| Primary age/sex-adjusted feature models | 69 | Pass | 1.79 × 10^-13 |
| Binary-sex sensitivity feature models | 69 | Pass | 8.17 × 10^-14 |
| CD45 event-QC feature models | 69 | Pass | 1.31 × 10^-13 |
| CD33+HLA-DR+ myeloid-like event-QC feature models | 69 | Pass | 2.20 × 10^-13 |
| Technical mismatch-exclusion feature models | 69 | Pass | 1.10 × 10^-13 |
| Subject-level composite scores | 850 subjects × 5 scores | Pass | 4.44 × 10^-15 |
| Primary and sensitivity composite models | 25 | Pass | 7.66 × 10^-15 |
| Same-sex nearest-age matching | 265 five-year pairs; 273 ten-year pairs | Exact subject and pair match | 0 |
| Age/sex-adjusted matched composite models | 10 | Pass | 3.16 × 10^-15 |
| Disease-by-age interaction models | 5 | Pass | 1.33 × 10^-15 |
| Threshold-set subject-level scores | 2,550 rows × 5 scores | Pass | 4.44 × 10^-15 |
| Targeted threshold-sensitivity models | 60 | Pass | 1.79 × 10^-13 |
| Threshold robustness summary | 20 variables | Exact classification match | 0 |

All numerical differences are at floating-point/machine-precision scale.

## Archive-specific rules confirmed

1. CP23 contains 850 FCS files/subjects; all 850 archived files were successfully processed and 832 subjects were model-ready after age/sex QC.
2. The primary analysis tests 69 threshold-defined monocyte/macrophage-like myeloid features under `feature ~ disease_group + age_for_model + sex`.
3. Prespecified feature sensitivities are binary sex, `n_cd45_dump_low >= 1000`, `n_cd33_hladr_myeloid_like >= 300`, and exclusion of channel/marker-order mismatch files.
4. Fifteen of 69 features preserve direction and global FDR across all five model sets; 48 preserve direction without FDR in every set, and six do not preserve direction across every sensitivity.
5. Five phenotype-oriented composite scores are standardized across all 850 post-extraction subjects.
6. The integrated monocyte/macrophage-like myeloid remodeling score has an archived adjusted beta of 0.2476159779193443 and is higher in cancer.
7. Same-sex age matching uses Female/Male-coded subjects, processes cancer subjects in ascending age order, and selects the nearest unused same-sex healthy control with lower-age tie resolution. The exact archived 265/273 subject pairs were reproduced.
8. Matched inference uses `score ~ disease_group + age_for_model + sex_match`; it is not a matched-pair fixed-effect model.
9. Disease-by-age interaction testing is restricted to the five composite scores and uses `score ~ disease_group * age_z + sex`.
10. Targeted threshold sensitivity tests 20 variables: 15 unique composite-component features and five composite scores. BH FDR is calculated separately within each main, permissive, and stringent threshold set.
11. Threshold robustness classification is exactly reproduced: 16 variables preserve direction and global FDR across all sets, four preserve direction without FDR in every set, and none fail direction preservation.
12. CP23 labels remain phenotype-based. Threshold-defined monocyte/macrophage-like, APC-like, and myeloid-like populations are not treated as definitive tissue macrophages, MDSCs, dendritic cells, or polarization states.

## Repository changes

- Added `scripts/31_validate_cp23_recovered.R`, a CP23-specific fixed-benchmark and optional archive-table validation gate.
- Updated `scripts/18_run_cp23_recovered.R` to invoke the panel-specific validator.
- The CP23 validator excludes only environment-dependent absolute path columns during archive-table comparison; all shared identifiers and scientific columns remain compared.
- Recovered Steps 3B–5 were retained unchanged because their downstream algorithms reproduce the archive exactly.

## Remaining validation

- Execute recovered Steps 1–2 on the local 850 CP23 FCS files.
- Compare regenerated inventory, mismatch diagnostics, compensation/transformation, gating, and feature tables with the archived reference tree.
- Run the complete CP23 pipeline with `SDY2583_CP23_REFERENCE_DIR` configured and retain the validation report and session information.
- Validate optional clinical-annotation and time-from-start repair extensions where the necessary clinical source tables are available.

CP23 is therefore **downstream archive-validated but raw-FCS and optional clinical-extension validation pending**.
