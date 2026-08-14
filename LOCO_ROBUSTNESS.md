# Locked leave-one-cancer-type-out robustness

This analysis asks whether the pan-cancer peripheral immune-remodeling signal is driven by one large tumor subgroup.

## Design

The analysis starts from the **training-frozen principal scores** produced by the locked internal held-out validation workflow. Score membership, component direction, centering, scaling, and weighting are not recalculated after subgroup exclusion.

Eight cancer groups are omitted one at a time: sarcoma, other cancers, breast, colorectal, pancreatic, lung, skin, and prostate. After each omission, the remaining cancer participants are compared with all healthy controls using:

```text
locked principal score ~ cancer status + age + sex
```

Inference uses HC3 robust covariance. Benjamini-Hochberg correction is applied across the ten principal scores within each omission. A second, more stringent correction is also applied across all 80 omission-by-score models.

## Aggregate result

- 80/80 cancer coefficients retained the positive cancer-associated direction.
- 79/80 models remained significant after the global 80-test BH correction.
- The only inferential exception was CP24 after sarcoma omission: beta approximately 0.067, 95% CI approximately -0.068 to 0.201, within-omission q approximately 0.329.
- All other nine principal axes remained FDR-significant after sarcoma omission.
- Removing any of the other seven cancer groups left all ten principal axes FDR-significant.

The interpretation is that the convergent pan-cancer architecture is not attributable to any single cancer subgroup, including the largest subgroup. CP24 is less stable: its direction remains positive after sarcoma removal, but its effect is materially attenuated and loses inferential support. This does **not** establish CP24 as sarcoma-specific.

## Run

Provide the same local inputs used for the locked internal validation, then run:

```bash
Rscript scripts/45_run_locked_loco_robustness.R
```

The runner first reconstructs the locked scores using TRAINING-only orientation/scaling and then performs LOCO in the same session. Participant-level locked scores are not written by this workflow. Aggregate outputs are written under `results/local/locked_loco_robustness/` by default and are ignored by Git.

## Public reference outputs

Aggregate-only reference summaries are distributed under `results/final/`:

- `locked_loco_baseline_results.csv`
- `locked_loco_score_summary.csv`
- `locked_loco_omission_summary.csv`
- `locked_loco_cancer_group_counts.csv`
- `locked_loco_beta_matrix.csv`

The runner additionally regenerates the full 80-model table locally as `locked_loco_principal_results.csv`; that detailed local output is not required as a public reference table.

These results are a cohort-composition robustness analysis, not an independent external validation.
