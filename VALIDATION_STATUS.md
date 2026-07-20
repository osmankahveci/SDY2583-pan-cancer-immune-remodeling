# Validation status

This document records completed validation evidence without promoting provisional reconstructed source beyond the evidence actually available.

## CP7 downstream archive validation — 2026-07-20

### Scope

The archived CP7 participant-level feature, metadata, score, age-sensitivity, and threshold-sensitivity CSV tables were independently recalculated outside the reconstructed R pipeline. This validation starts from archived post-extraction tables; it does **not** replace the pending raw-FCS execution of Steps 1–2.

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

### Remaining CP7 validation

- execute FCS inventory, compensation, fixed-logicle transformation, gating, and feature extraction on the local 850 CP7 FCS files;
- compare the generated Step 1–2 outputs with the archived feature table and QC summaries;
- execute the corrected full CP7 R pipeline and retain its validation report plus session information.

Until those FCS-level checks pass, CP7 is best described as **downstream archive-validated but raw-FCS validation pending**.

## CP8 downstream archive validation — 2026-07-20

### Scope

The archived CP8 post-extraction participant, score, binary-sex, age-sensitivity, matching, interaction, and threshold-sensitivity CSV tables were independently recalculated outside the reconstructed R pipeline. As for CP7, this establishes downstream reproducibility from archived participant-level tables but does **not** replace pending raw-FCS execution of inventory, compensation, transformation, gating, and feature extraction.

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

### Remaining CP8 validation

- execute Steps 1–2 on the local 850 CP8 FCS files and compare the generated feature table and QC summaries with the archive;
- execute the corrected full CP8 R pipeline with `SDY2583_CP8_REFERENCE_DIR` configured and retain its validation report plus session information;
- evaluate the reconstructed event-count sensitivity extension during the local run; no separate archived CP8 event-count result table was available for independent table-level comparison in this validation pass.

CP8 is therefore **downstream archive-validated but raw-FCS validation pending**.
