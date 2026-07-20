# CP26 downstream archive validation

Date: 2026-07-20

## Scope

Archived CP26 participant data, 15-feature primary and sensitivity models, five subject-level NK composite scores, same-sex age matching, disease-by-age interactions, and the 20-variable threshold-sensitivity analysis were independently recalculated. This validation identified an incorrect main-score standardization universe and an incorrect integrated-score construction in the initial reconstructed source; the manifest, score, threshold, matching-output, and validation code were corrected.

This validates downstream reproducibility from archived participant-level tables. It does not replace local execution against the 850 raw CP26 FCS files.

## Results

| Layer | Archived rows | Result | Maximum absolute numerical difference |
|---|---:|---|---:|
| Primary age/sex-adjusted feature models | 15 | Pass | 2.72 × 10^-13 |
| Subject-level composite scores | 832 subjects × 5 scores | Pass | 9.77 × 10^-15 |
| Primary composite-score models | 5 | Pass | 7.66 × 10^-15 |
| Same-sex nearest-age matching | 263 five-year pairs; 272 ten-year pairs | Exact pair match | 0 |
| Five-year matched models | 20 | Pass | 9.45 × 10^-13 |
| Ten-year matched models | 20 | Pass | 1.36 × 10^-12 |
| Disease-by-age interaction models | 20 | Pass | 2.09 × 10^-14 |
| Threshold-set subject-level scores | 2,550 rows × 5 scores | Pass | 9.33 × 10^-15 |
| Targeted threshold-sensitivity models | 60 | Pass | 2.65 × 10^-13 |
| Threshold robustness summary | 20 variables | Exact classification match | 0 |

## Archive-specific rules confirmed

1. CP26 contains 850 FCS files/subjects and 832 model-ready subjects after valid-age QC.
2. Seven annotation-variant files require the dump/viability channel fallback to `BV510-A`; repaired extraction contains no remaining failed files.
3. CD45 is measured on `BB700-A`. The threshold-defined NK-like gate is CD45-positive, dump-low, and CD56-positive or CD16-positive.
4. The primary model family contains exactly 15 prespecified NK composition, CD56/CD16-balance, receptor-phenotype, and activation/maturation outcomes.
5. Primary feature inference uses `outcome ~ disease_group + age_for_model + sex`.
6. Main composite scores are standardized only within the 832 model-ready subjects, not across all 850 post-extraction rows.
7. Four domain scores are means of oriented raw-feature z scores.
8. The integrated NK-remodeling score is the mean of ten oriented raw-feature z scores: one positively oriented NKG2C endpoint and nine negatively oriented composition/maturation endpoints. It is not a mean of re-standardized domain scores.
9. The archived integrated adjusted beta is 0.364755491769086.
10. Same-sex age matching uses the standard ascending sex/age/subject ordering and subject-ID tie resolution, reproducing exactly 263 five-year and 272 ten-year pairs.
11. Matched inference remains age/sex adjusted and tests the exact 20-variable family: 15 raw features plus five scores.
12. Disease-by-age interaction testing uses the same 20 variables and fits `outcome ~ disease_group * age_z + sex`.
13. Threshold-set scores are standardized independently within each all-850 main, permissive, and stringent extraction.
14. Threshold sensitivity tests the same 20 variables. Positive-marker thresholds shift by -0.2/+0.2, while the dump-low threshold shifts in the opposite direction.
15. Threshold robustness classification is exactly reproduced: ten variables preserve direction and FDR, seven preserve direction without FDR in every set, and three do not preserve direction.

## Corrections applied to reconstructed CP26 source

- Added the exact ten-feature integrated-score definition to the panel manifest.
- Changed main score standardization from all 850 rows to the archived 832-subject model-ready set.
- Replaced the generic integrated domain-score calculation with the archived oriented raw-feature calculation.
- Retained all-850, within-threshold-set standardization for threshold sensitivity.
- Added archive-compatible threshold scored-data and robustness outputs.
- Rebuilt matching outputs with exact `pair_n` identifiers, age differences, separate five- and ten-year model tables, and archive-compatible interaction columns.
- Expanded the CP26 validator to cover feature extraction, score maps, subject-level scores, exact matching, matched models, interactions, threshold scored data, and robustness classes.

## Remaining validation

- Execute inventory, compensation, fixed-logicle transformation, repaired dump-channel mapping, gating, and feature extraction on the local 850 CP26 FCS files.
- Confirm that exactly the annotation-variant files use the `BV510-A` fallback and that no required-marker failures remain.
- Run the corrected complete CP26 pipeline with `SDY2583_CP26_REFERENCE_DIR` configured and retain its validation report and session information.
- Validate optional clinical-annotation outputs where source clinical tables are available.

CP26 is therefore **downstream archive-validated after score-construction corrections, but raw-FCS and optional clinical-extension validation remain pending**.
