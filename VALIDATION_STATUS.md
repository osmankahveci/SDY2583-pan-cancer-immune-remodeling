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

## CP8 downstream archive validation — in progress, 2026-07-20

### Completed layers

The archived CP8 post-extraction participant table was used to independently refit the primary models with `outcome ~ disease_group + age_for_model + sex`, retaining all archived sex categories and complete-case model sets.

| Layer | Archived rows | Result | Maximum absolute numerical difference |
|---|---:|---|---:|
| Primary age/sex-adjusted feature models | 44 | Pass | 2.69 × 10^-13 |
| Composite-score adjusted models | 6 | Pass | 6.93 × 10^-14 for coefficient/t statistics; p-value differences below 1.53 × 10^-18 |

The primary and composite model outputs therefore reproduce at floating-point/machine-precision scale.

### Remaining CP8 validation

- independently reconstruct and compare subject-level composite-score definitions;
- validate age-stratified models, same-sex nearest-age matching, matched models, and disease-by-age interactions;
- validate binary-sex, event-count, direction, and targeted threshold-sensitivity families;
- inspect the reconstructed CP8 source for any archive-specific implementation differences and apply corrections;
- execute raw-FCS extraction locally when the 850 FCS files are available.

CP8 is not yet marked complete.