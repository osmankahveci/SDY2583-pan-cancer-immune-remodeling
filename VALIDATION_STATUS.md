# Validation status

## Accepted release-validation boundary — 2026-07-20

All ten manuscript panels have completed downstream archive validation from archived participant-level feature, metadata, score, matching, interaction, and threshold-sensitivity tables. Numerical agreement is at floating-point/machine-precision scale, with exact subject-pair reproduction where archived matching tables were available.

A fresh raw-FCS end-to-end rerun was not performed. FCS inventory, compensation, fixed-logicle transformation, threshold gating, and regeneration of participant-level feature tables are explicitly outside the accepted validation scope for this release. This omission must remain visible in repository and manuscript-facing reproducibility wording and must not be described as completed raw-data validation.

Accepted description:

> Provenance-tracked all-panel reproducibility release with all ten panel downstream analyses independently archive-validated. Raw-FCS end-to-end re-execution was not performed and is outside the accepted validation scope of this release.

The detailed panel records below preserve the original evidence and the previously identified optional raw-FCS follow-up tasks. Those tasks are no longer release blockers under the accepted scope decision.

## CP7 downstream archive validation — 2026-07-20

### Scope

The archived CP7 participant-level feature, metadata, score, age-sensitivity, and threshold-sensitivity CSV tables were independently recalculated outside the reconstructed R pipeline. This validation starts from archived post-extraction tables; it does **not** replace a raw-FCS execution of Steps 1–2.

### Results

| Layer | Archived rows | Result | Maximum absolute numerical difference |
|---|---:|---|---:|
| Primary age/sex-adjusted feature models | 55 | Pass | 1.52 × 10^-13 |
| Subject-level component and integrated scores | 832 subjects × 6 scores | Pass | 7.99 × 10^-15 |
| Composite-score adjusted models | 6 | Pass | 9.95 × 10^-14 |
| Age-stratified adjusted models | 96 | Pass | 5.80 × 10^-14 |
| Same-sex nearest-age subject matching | 302 five-year pairs; 332 ten-year pairs | Exact pair match | 0 |
| Matched pair-fixed-effect models | 64 | Pass | 7.99 × 10^-14 |
| Disease-by-age interaction models | 32 | Pass | 2.35 × 10^-14 |
| Targeted threshold-sensitivity models | 33 | Pass | 7.55 × 10^-14 |

All differences are at floating-point/machine-precision scale.

### Archive-specific rules confirmed

1. Primary CP7 component and integrated scores are standardized in the **832-subject model-ready set**, not in all 850 feature-extracted subjects.
2. Same-sex nearest-age matching assigns row IDs after binary-sex filtering, processes cancer subjects from **oldest to youngest**, and greedily selects the nearest unused same-sex healthy control.
3. Matched inference uses `outcome ~ disease_group + matched_pair_id`, equivalently a paired-difference model.
4. Age sensitivity uses an exact 32-outcome family.
5. Targeted threshold sensitivity uses an exact 11-outcome family; BH FDR is calculated separately within permissive, main, and stringent threshold sets.
6. Threshold-set composite scores are standardized within each all-850 threshold extraction, which explains the small expected difference between the primary integrated beta (0.351119363746345) and the main-threshold integrated beta (approximately 0.351140).

### Corrections applied to reconstructed CP7 source

- primary composite standardization changed from 850 subjects to the archived 832-subject model-ready set;
- age-sensitivity outcome family changed to the archived 32 outcomes;
- matching order changed from ascending to descending cancer age;
- matched model changed from age/sex adjustment to matched-pair fixed effects;
- official targeted-threshold FDR family changed from 28 candidate outcomes to the archived 11 outcomes;
- CP7 validation expanded to compare subject-level scores, matching, age strata, interactions, and threshold tables.

### Optional future raw-FCS follow-up

- execute FCS inventory, compensation, fixed-logicle transformation, gating, and feature extraction on the local 850 CP7 FCS files;
- compare the generated Step 1–2 outputs with the archived feature table and QC summaries;
- execute the corrected full CP7 R pipeline and retain its validation report plus session information.

CP7 is **downstream archive-validated**; fresh raw-FCS re-execution was not performed and is outside the accepted release scope.

## CP8 downstream archive validation — 2026-07-20

### Scope

The archived CP8 post-extraction participant, score, binary-sex, age-sensitivity, matching, interaction, and threshold-sensitivity CSV tables were independently recalculated outside the reconstructed R pipeline. This establishes downstream reproducibility from archived participant-level tables but does **not** represent a fresh raw-FCS execution of inventory, compensation, transformation, gating, and feature extraction.

### Results

| Layer | Archived rows | Result | Maximum absolute numerical difference |
|---|---:|---|---:|
| Primary age/sex-adjusted feature models | 44 | Pass | 2.69 × 10^-13 |
| Subject-level component and integrated scores | 850 subjects × 6 scores | Pass | 1.16 × 10^-14 |
| Composite-score adjusted models | 6 | Pass | 6.93 × 10^-14 |
| Binary-sex sensitivity feature models | 44 | Pass | 1.05 × 10^-13 |
| Age-stratified adjusted models | 93 | Pass | 7.46 × 10^-14 |
| Same-sex nearest-age subject matching | 265 five-year pairs; 273 ten-year pairs | Exact pair match | 0 |
| Age/sex-adjusted matched models | 62 | Pass | 4.91 × 10^-13 |
| Disease-by-age interaction models | 31 | Pass | 2.49 × 10^-14 |
| Threshold-set subject-level scores | 2,550 rows × 6 scores | Pass | 1.16 × 10^-14 |
| Targeted threshold-sensitivity models | 84 | Pass | 2.08 × 10^-13 |
| Threshold robustness summary | 28 outcomes | Pass: direction and FDR preserved for all 28 | Exact classification match |

All numerical differences are at floating-point/machine-precision scale.

### Archive-specific rules confirmed

1. CP8 component and integrated scores are standardized across **all 850 post-extraction subjects**. This differs from the CP7 primary-score rule.
2. Age-stratified, matched, and interaction analyses use an exact **31-outcome** family.
3. Age-stratum BH FDR is calculated separately within Young <40, Middle 40–59, and Older 60+.
4. Same-sex nearest-age matching sorts cancer and healthy subjects by sex, ascending age, and subject ID. When two controls have the same age distance, the first control in ascending-age order is selected; this reproduces the archived lower-age-first tie rule.
5. Matched CP8 inference uses `outcome ~ disease_group + age_for_model + sex`; it is not a matched-pair fixed-effect model.
6. Disease-by-age interaction models use age standardized over the complete valid-age CP8 set and fit `outcome ~ disease_group * age_z + sex`.
7. Threshold sensitivity uses an exact **28-outcome** family. Positive-marker thresholds use −0.2/+0.2 permissive/stringent shifts, whereas the IL7RA-low gate uses the reversed +0.2/−0.2 rule.
8. Threshold-set component and integrated scores are standardized independently within each all-850 permissive, main, and stringent extraction. BH FDR is calculated separately within each 28-outcome threshold set.

### Corrections applied to reconstructed CP8 source

- added explicit archived 31-outcome age-sensitivity and 28-outcome threshold-sensitivity families to the panel manifest;
- corrected same-sex nearest-age tie handling from subject-ID-first to lower-age-first selection, reproducing the exact 265/273 archived pairs;
- changed age-stratum FDR from a pooled 93-test adjustment to separate 31-test adjustments within each age group;
- retained the archive-correct age/sex-adjusted matched model rather than applying CP7’s pair-fixed-effect rule;
- corrected Step 5 to load the Step 1 FCS inventory before using `fcs_files`;
- expanded threshold outputs to include the full data-with-metadata, scored-data, 28-outcome statistics, and robustness tables;
- expanded the CP8 validator to cover subject-level scores, binary-sex models, age strata, exact pairs, matched models, interactions, and threshold outputs;
- corrected the fixed main-feature benchmark from 45 to the archived 44 modeled features.

### Optional future raw-FCS follow-up

- execute Steps 1–2 on the local 850 CP8 FCS files and compare the generated feature table and QC summaries with the archive;
- execute the corrected full CP8 R pipeline with `SDY2583_CP8_REFERENCE_DIR` configured and retain its validation report plus session information;
- evaluate the reconstructed event-count sensitivity extension during a local run; no separate archived CP8 event-count result table was available for independent table-level comparison in this validation pass.

CP8 is **downstream archive-validated**; fresh raw-FCS re-execution was not performed and is outside the accepted release scope.

## CP10 downstream archive validation — 2026-07-20

### Scope

The archived CP10 post-extraction participant table, feature-model outputs, composite definitions and scores, age-matched datasets, disease-by-age interaction results, and threshold-sensitivity tables were independently recalculated. Recovered `SAFE` Steps 2–5 were also inspected against the archive-specific algorithms. This validates the numerical downstream core but does **not** represent local execution of the reconstructed Step 1/1B and recovered Step 2 against the 850 raw FCS files.

### Results

| Layer | Archived rows | Result | Maximum absolute numerical difference |
|---|---:|---|---:|
| Primary age/sex-adjusted feature models | 37 | Pass | 2.65 × 10^-13 |
| Binary-sex sensitivity feature models | 37 | Pass | 2.10 × 10^-13 |
| CD45 event-QC sensitivity feature models | 37 | Pass | 1.39 × 10^-13 |
| Technical mismatch-exclusion feature models | 37 | Pass | 1.28 × 10^-13 |
| Subject-level composite scores | 850 subjects × 6 scores | Pass | 1.78 × 10^-15 |
| Primary and sensitivity composite models | 24 | Pass | 8.22 × 10^-15 |
| Same-sex nearest-age subject matching | 265 five-year pairs; 273 ten-year pairs | Exact pair match | 0 |
| Age/sex-adjusted matched composite models | 12 | Pass | 5.00 × 10^-15 |
| Disease-by-age interaction models | 6 | Pass | 1.64 × 10^-14 |
| Threshold-set subject-level scores | 2,550 rows × 6 scores | Pass | 2.66 × 10^-15 |
| Targeted threshold-sensitivity models | 99 | Pass | 2.65 × 10^-13 |
| Threshold robustness summary | 33 variables | Exact classification match | 0 |

All numerical differences are at floating-point/machine-precision scale.

### Archive-specific rules confirmed

1. CP10 tests an exact **37-feature** family under the primary model and three prespecified sensitivities: binary sex, `n_cd45_viable >= 1000`, and exclusion of the two channel/marker-mismatch files.
2. All 37 feature directions are preserved across the four model sets; global FDR remains below 0.05 for 28 features in all sets.
3. Six phenotype-oriented composite scores are standardized across all **850 post-extraction subjects**. Five are significantly higher in cancer; the CCR3/eosinophil-like composite is not significant in the primary model.
4. The integrated myeloid/granulocytic remodeling score has an archived adjusted beta of **0.4893007219212893**.
5. Same-sex age matching sorts healthy and cancer subjects by sex, ascending age, and subject ID, then greedily selects the nearest unused control with lower-age-first and subject-ID tie resolution. The exact archived 265/273 pairs were reproduced.
6. Matched inference remains age/sex adjusted with `score ~ disease_group + age_for_model + sex`; it is not a pair-fixed-effect model.
7. Disease-by-age interaction models are restricted to the six composite scores and fit `score ~ disease_group * age_z + sex` in the 832-subject model-ready set.
8. Threshold sensitivity uses **33 variables**: 27 constituent/targeted features and six composites. BH FDR is calculated separately within each main, permissive, and stringent 33-variable set.
9. Threshold robustness classification is exactly reproduced: 29 variables preserve direction and global FDR, three preserve direction without FDR across all sets, and one does not preserve direction.
10. CP10 labels remain phenotype-based. The code does not convert granulocyte-like, eosinophil-like, APC-like, or monocyte-like flow phenotypes into definitive cell-identity claims.

### Repository changes made during validation

- added `scripts/27_validate_cp10_recovered.R`, a CP10-specific fixed-benchmark and optional archive-table validation gate;
- changed the CP10 runner to invoke the panel-specific validator rather than the generic recovered-panel structural validator;
- configured the validator to compare feature, model, subject-level composite, exact matching, interaction, and threshold tables when `SDY2583_CP10_REFERENCE_DIR` is provided;
- retained recovered Steps 2–5 unchanged because their downstream algorithms already reproduce the archive exactly.

### Optional future raw-FCS and clinical follow-up

- execute the reconstructed Step 1/1B and recovered Step 2 on the local 850 CP10 FCS files;
- compare the regenerated inventory, marker/channel mismatch, compensation/transformation, gating, and feature tables with the archived outputs;
- run the complete CP10 pipeline with `SDY2583_CP10_REFERENCE_DIR` configured and retain the panel-specific validation report and session information;
- separately inspect optional clinical-annotation and cancer-subtype outputs where the necessary clinical source tables are available.

CP10 is **downstream archive-validated**; fresh raw-FCS re-execution and optional clinical-extension validation were not performed and are outside the accepted release scope.
