# Original five-immunotype association with integrated remodeling scores

Status: **analysis completed; retain for manuscript integration only after all planned additional analyses are complete. No figure should be generated yet.**

## Dataset and design

Input: `SDY2583_integrated_clinical_immune_score_matrix_ALL10_with_CP23.csv`

- 850 participants
- 5 original immunotypes: G1_Naive (n=178), G2_Primed (n=294), G3_Progressive (n=205), G4_Chronic (n=95), G5_Suppressive (n=78)
- Training cohort: n=503; validation cohort: n=347
- Primary outcomes: one principal integrated remodeling score per each of the 10 cytometry panels
- Full-cohort model: score ~ immunotype + age + sex + disease group
- Training and validation models: same covariate structure, analyzed separately
- Cancer-only model: score ~ immunotype + age + sex + cancer subgroup
- Treatment sensitivity: cancer-only model additionally adjusted for treatment-status category
- Heteroskedasticity-consistent HC3 covariance used for univariate models
- Omnibus immunotype tests controlled across the 10 principal scores with Benjamini-Hochberg FDR
- Pairwise post-hoc contrasts controlled within score

## Confounding structure requiring adjustment

Immunotype composition was strongly non-random with respect to age and disease group:

- G1_Naive: 44 cancer / 134 healthy; mean age 39.8 years
- G2_Primed: 134 cancer / 160 healthy; mean age 56.3 years
- G3_Progressive: 128 cancer / 77 healthy; mean age 56.2 years
- G4_Chronic: 64 cancer / 31 healthy; mean age 61.9 years
- G5_Suppressive: 72 cancer / 6 healthy; mean age 62.7 years

Therefore, unadjusted immunotype comparisons should not be used for primary inference.

## Primary results

The adjusted omnibus immunotype effect was significant for all 10 principal integrated scores in the full cohort, training cohort, validation cohort, cancer-only cohort, and cancer-only treatment-adjusted cohort.

### Full-cohort omnibus FDR values

- CP28 integrated T/NK-interface remodeling: q = 5.15e-59
- CP8 integrated CD4 helper/regulatory remodeling: q = 2.24e-45
- CP24 integrated CD8 remodeling: q = 8.76e-25
- CP10 integrated myeloid/granulocytic remodeling: q = 8.23e-23
- CP25 integrated CD4 regulatory-checkpoint remodeling: q = 4.30e-20
- CP26 integrated NK remodeling: q = 2.69e-18
- CP7 integrated checkpoint remodeling: q = 1.94e-17
- CP23 integrated monocyte/macrophage-like myeloid remodeling: q = 3.23e-09
- CP22 integrated humoral B-cell remodeling: q = 3.68e-08
- CP16 integrated APC/DC-like myeloid remodeling: q = 7.09e-08

A multivariate model across the 10 principal scores also showed a strong global immunotype effect (Pillai trace = 0.716; approximate F = 17.67; P < 0.001).

## Incremental variance associated with immunotype

After accounting for age, sex, and disease group, addition of immunotype increased model R-squared by:

- CP28: 0.236
- CP26: 0.161
- CP10: 0.158
- CP8: 0.151
- CP24: 0.112
- CP7: 0.086
- CP25: 0.084
- CP23: 0.051
- CP16: 0.050
- CP22: 0.038

These values are descriptive incremental R-squared estimates, not causal attributable fractions.

## Adjusted immunotype profiles

### G1_Naive

G1 occupied the lowest adjusted position across nearly all integrated remodeling axes. It is the low-remodeling reference state, particularly for CD4 helper/regulatory, CD4 regulatory-checkpoint, CD8 differentiation, myeloid, APC/DC-like, humoral, and checkpoint-related scores.

### G2_Primed

G2 generally remained close to G1 but showed an early increase in CD4 helper/regulatory remodeling. The adjusted G2-G1 difference for CP8 was +0.343 score units and replicated in training and validation cohorts. Most NK, CD8/T-NK, and myeloid axes remained comparatively low.

### G3_Progressive

G3 was characterized most consistently by intermediate-to-high CD4 remodeling. Compared with G1, CP8 increased by +0.635 and CP25 by +0.542 score units. Humoral B-cell remodeling and myeloid/granulocytic remodeling were also elevated, while terminal CD8/T-NK remodeling remained below G4.

### G4_Chronic

G4 showed the strongest CD8/T-NK-oriented phenotype. It had the highest adjusted CP24 integrated CD8 remodeling and CP28 integrated T/NK-interface remodeling scores. Compared with G1, the adjusted differences were +0.866 for CP24 and +0.729 for CP28. G4 was also elevated for CP7 checkpoint remodeling and CP25 regulatory-checkpoint remodeling, but did not show the marked CP10/CP23/CP26 pattern observed in G5.

### G5_Suppressive

G5 showed the broadest innate/myeloid and NK remodeling phenotype. It had the highest adjusted CP10 myeloid/granulocytic, CP16 APC/DC-like, CP22 humoral B-cell, CP23 monocyte/macrophage-like, CP26 NK, CP7 checkpoint, and CP8 CD4 helper/regulatory scores. The largest G5-G1 differences were CP10 +1.308, CP26 +0.913, CP8 +0.829, CP16 +0.602, CP7 +0.574, and CP23 +0.560 score units. G5 remained below G4 on CP24 and CP28, distinguishing suppressive/myeloid-NK remodeling from the chronic terminal CD8/T-NK phenotype.

## Internal training-validation replication

All 10 principal score omnibus effects were FDR-significant independently in both the training and validation subsets.

Among the 40 prespecified G1-referenced contrasts, 21 were FDR-significant with the same direction in the full, training, and validation analyses. Strong internally replicated contrasts included:

- G4 vs G1: CP28 and CP24
- G3, G4, and G5 vs G1: CP8
- G3, G4, and G5 vs G1: CP25, with the G5 training estimate weaker but still replicated under the predefined FDR criterion
- G5 vs G1: CP10, CP16, CP23, CP26, CP7, CP22, and CP28
- G3 vs G1: CP10, CP22, CP24, and CP25

This is an internal split-sample replication of score-immunotype mapping, not an external validation cohort.

## Cancer-only and treatment sensitivity

All 10 principal integrated scores retained FDR-significant immunotype effects in cancer-only models adjusted for age, sex, and cancer subgroup. Adding treatment-status category produced negligible changes in pairwise contrast estimates: median absolute change = 0.002 score units; maximum absolute change = 0.015 score units. Thus, the immunotype-remodeling profiles were not materially explained by recorded treatment status.

## Full 66-score supplementary screen

- Full cohort: 63/66 score-level omnibus tests FDR-significant
- Training: 52/66
- Validation: 57/66
- Cancer-only: 59/66

The 66-score screen is supplementary because many component scores overlap and are not statistically independent.

## Interpretive boundary

The original immunotypes and the present remodeling scores were derived from the same underlying SDY2583 flow-cytometry resource. Therefore, these analyses should be framed as biological resolution, decomposition, and internal replication of the original immunotypes—not as independent validation of their existence. Strong claims about external reproducibility, prognosis, treatment prediction, or causal immune states are not supported by this analysis alone.

## Working Results text

To biologically resolve the five previously defined peripheral immunotypes, one principal integrated remodeling score from each cytometry panel was modeled as a function of immunotype after adjustment for age, sex, and disease group. Immunotype was associated with all ten remodeling axes after false-discovery-rate correction, with the largest incremental contributions observed for T/NK-interface, NK, myeloid/granulocytic, and CD4 helper/regulatory remodeling. The adjusted profiles formed biologically distinct rather than uniformly ordered states. G1_Naive occupied the lowest remodeling position across most axes, whereas G2_Primed showed an early shift dominated by CD4 helper/regulatory remodeling. G3_Progressive exhibited intermediate CD4 regulatory, humoral, and myeloid remodeling. G4_Chronic was distinguished by the highest CD8 differentiation and T/NK-interface remodeling, including marked elevations in CP24 and CP28 scores. In contrast, G5_Suppressive showed the strongest myeloid/granulocytic, monocyte/macrophage-like, APC/DC-like, and NK remodeling, together with broad CD4 and checkpoint-associated changes. All ten omnibus associations were reproduced separately in the training and validation subsets and remained significant in cancer-only models adjusted for cancer subgroup. Additional adjustment for treatment status had negligible effects on the estimated immunotype contrasts, indicating that the observed profiles were not materially explained by recorded treatment exposure.

## Working Discussion text

Mapping the composite remodeling framework onto the original five immunotypes provided an interpretable biological decomposition of the cluster labels. The profiles did not represent a simple monotonic continuum from G1 to G5. Instead, two advanced but distinct remodeling configurations emerged. G4_Chronic was dominated by CD8 differentiation, terminal-effector/checkpoint, and T/NK-interface remodeling, whereas G5_Suppressive was characterized by myeloid/granulocytic and monocyte/macrophage-like remodeling together with marked NK and broader multi-compartment changes. G2_Primed and G3_Progressive occupied intermediate positions, with progressively stronger CD4 helper/regulatory and humoral remodeling. These relationships remained after demographic and disease-group adjustment and were reproduced across the original training and validation partitions, supporting the internal stability of the biological mapping. Nevertheless, because both the immunotype labels and the remodeling scores originate from the same flow-cytometry resource, this analysis should be interpreted as mechanistic annotation and dimensional resolution of the original immunotypes rather than independent external validation.

## Figure policy

Do not generate an immunotype-specific figure now. Retain all results for the final integrated figure decision after PCA, clustering/stability, and transcriptomic analyses are complete. The likely final network may encode immunotype-specific loading or enrichment around biological-domain nodes, but only if this improves clarity without exceeding the one-figure, maximum-two-figure manuscript strategy.
