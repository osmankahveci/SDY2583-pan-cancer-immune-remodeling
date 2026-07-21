# Principal integrated-score PCA: global remodeling structure

Status: **analysis completed; retain for later manuscript integration after all planned additional analyses are complete. Do not generate a figure yet.**

## Dataset and analysis population

Input: `SDY2583_integrated_clinical_immune_score_matrix_ALL10_with_CP23.csv`

- 850 total participants
- 829 complete cases across the 10 principal integrated scores, age, sex, disease group, immunotype, and cohort
- Complete cases: training n = 485; validation n = 344
- Outcomes: one principal integrated remodeling score from each of CP7, CP8, CP10, CP16, CP22, CP23, CP24, CP25, CP26, and CP28

Two PCA spaces were evaluated:

1. **Age/sex-adjusted PCA:** each score was residualized for age and sex, standardized, and entered into PCA. This space retains the observed cancer-versus-healthy remodeling signal.
2. **Disease-adjusted PCA:** each score was residualized for age, sex, and disease group before standardization and PCA. This sensitivity analysis tests whether the global structure and immunotype mapping persist beyond the cancer-versus-healthy mean shift.

PCA signs were oriented for interpretability:

- positive PC1 = increasing multi-compartment remodeling;
- positive PC2 = myeloid/monocyte/NK predominance relative to CD8/T-NK and regulatory-checkpoint predominance.

## Variance structure

### Age/sex-adjusted PCA

- PC1: 35.6%
- PC2: 15.8%
- PC1–PC2 cumulative: 51.4%
- PC1–PC4 cumulative: 68.9%

### Disease-adjusted PCA

- PC1: 30.5%
- PC2: 16.6%
- PC1–PC2 cumulative: 47.1%
- PC1–PC4 cumulative: 66.3%

Thus, removal of the cancer-versus-healthy mean shift reduced PC1 variance only modestly and preserved the two-axis architecture.

## Biological interpretation of the principal axes

### PC1: broad multi-compartment remodeling gradient

All ten disease-adjusted PC1 loadings were positive (range 0.205–0.400). The largest loadings were:

- CP8 CD4 helper/regulatory remodeling: 0.400
- CP28 T/NK-interface remodeling: 0.393
- CP7 checkpoint remodeling: 0.348
- CP25 CD4 regulatory-checkpoint remodeling: 0.338
- CP26 NK remodeling: 0.317

The uniform loading direction supports interpretation of PC1 as a general peripheral immune-remodeling burden rather than a single-cell-lineage axis.

### PC2: myeloid–NK versus CD8/T-NK remodeling configuration

Disease-adjusted PC2 had positive loadings for:

- CP23 monocyte/macrophage-like remodeling: +0.465
- CP10 myeloid/granulocytic remodeling: +0.446
- CP26 NK remodeling: +0.330

and negative loadings for:

- CP24 CD8 remodeling: −0.436
- CP28 T/NK-interface remodeling: −0.328
- CP25 CD4 regulatory-checkpoint remodeling: −0.310

PC2 therefore distinguishes two biologically different advanced remodeling configurations rather than measuring overall remodeling magnitude.

### PC3: APC/NK versus CD4 regulatory-humoral configuration

Disease-adjusted PC3 was driven positively by CP16 APC/DC-like remodeling (+0.533) and CP26 NK remodeling (+0.321), and negatively by CP8 CD4 helper/regulatory remodeling (−0.413), CP25 CD4 regulatory-checkpoint remodeling (−0.380), and CP22 humoral B-cell remodeling (−0.285). PC3 accounted for 10.2% of variance and should remain secondary because training-validation stability was lower than for PC1–PC2.

## Disease and immunotype contributions to global score space

Permutation-based multivariate variance partitioning used 4,999 permutations.

- Disease group accounted for 6.8% of the age/sex-adjusted 10-score variance (permutation P = 0.0002).
- Immunotype accounted for an additional 11.4% of score-space variance after controlling for disease group (permutation P = 0.0002).
- In the independently disease-adjusted score space, immunotype accounted for 11.7% of multivariate variance (permutation P = 0.0002).

Therefore, the immunotype organization is not explained solely by the cancer-versus-healthy contrast.

## Immunotype centroids in disease-adjusted PCA space

| Immunotype | n | PC1 mean | PC2 mean | PC3 mean |
|---|---:|---:|---:|---:|
| G1_Naive | 170 | −0.766 | +0.082 | +0.300 |
| G2_Primed | 288 | −0.499 | +0.124 | −0.065 |
| G3_Progressive | 202 | +0.125 | −0.073 | −0.396 |
| G4_Chronic | 93 | +1.088 | −1.252 | +0.378 |
| G5_Suppressive | 76 | +1.941 | +1.073 | +0.168 |

PC1 ordered the adjusted immunotype centroids from G1 through G5, consistent with increasing global remodeling. PC2 then separated the two highest-remodeling immunotypes:

- G4_Chronic: CD8/T-NK-oriented pole
- G5_Suppressive: myeloid/monocyte/NK-oriented pole

Key HC3 contrasts in the disease-adjusted PCA space:

- PC1, G5 versus G1: difference +2.707, P < 1e-20
- PC1, G4 versus G1: difference +1.854, P < 1e-18
- PC1, G5 versus G4: difference +0.853, P = 0.0042
- PC2, G4 versus G1: difference −1.333, P < 1e-18
- PC2, G5 versus G1: difference +0.991, P < 1e-6
- PC2, G5 versus G4: difference +2.325, P < 1e-20

## PC-level immunotype associations

After disease adjustment in the age/sex-adjusted PCA:

- PC1 immunotype incremental R² = 0.169; partial R² = 0.207
- PC2 immunotype incremental R² = 0.178; partial R² = 0.179
- PC3 immunotype incremental R² = 0.069; partial R² = 0.070

In the disease-adjusted PCA space:

- immunotype R² for PC1 = 0.226
- immunotype R² for PC2 = 0.174
- immunotype R² for PC3 = 0.075

All three associations were FDR-significant. PC1 and PC2 contain the principal interpretable immunotype structure.

## Continuous versus categorical organization

Despite significant centroid differences, immunotype silhouette values were negative:

- age/sex-adjusted 10-score space: −0.031
- disease-adjusted 10-score space: −0.031
- disease-adjusted PC1–PC2 space: −0.065

This means that the compressed 10-score representation does not form five non-overlapping islands. Instead, the immunotypes occupy distinct mean positions within broadly overlapping continuous remodeling space. This does **not** invalidate the original immunotypes, which were derived from a much larger cellular-state feature set. It indicates that the present integrated scores biologically annotate the immunotypes as overlapping configurations rather than perfectly separable classes.

## Training-validation stability

Separate disease-adjusted PCA models were fitted in the training and validation subsets.

- PC1 variance explained: training 31.3%; validation 31.4%
- PC2 variance explained: training 16.4%; validation 16.0%
- PC1 loading cosine similarity: 0.953
- PC2 loading cosine similarity after sign alignment: 0.884
- First-two-PC subspace principal angles: 11.4° and 24.3° (cosines 0.980 and 0.911)

A PCA fitted only in the training cohort was then used to project validation participants. Across the five immunotype centroids:

- PC1 centroid Pearson correlation, training versus validation: r = 0.984
- PC2 centroid Pearson correlation, training versus validation: r = 0.977
- PC1 immunotype rank ordering was identical in training and validation.

The first two remodeling axes are therefore internally stable. Higher-order components were less stable and should not drive the primary narrative.

## Primary interpretation

The integrated score space is best described by two complementary dimensions:

1. **PC1 — global multi-compartment remodeling burden**, spanning adaptive, humoral, NK, and myeloid domains.
2. **PC2 — remodeling configuration**, separating a CD8/T-NK chronic-effector pattern from a myeloid/monocyte/NK suppressive pattern.

G1_Naive and G2_Primed occupy the low-PC1 region, G3_Progressive is intermediate, and G4_Chronic and G5_Suppressive both occupy high-PC1 positions but separate strongly along PC2. This provides a compact systems-level explanation for why G4 and G5 are both advanced immunotypes yet biologically different.

## Interpretive boundary

- PCA is descriptive and does not establish temporal progression, causality, or prognosis.
- The G1-to-G5 PC1 ordering should be described as an adjusted centroid gradient, not a demonstrated longitudinal trajectory.
- Negative silhouette values argue against presenting the five immunotypes as sharply discrete in the 10-score compressed space.
- Training-validation analyses are internal split-sample replication, not external validation.
- No PCA figure should be generated now. PCA elements should be considered only when constructing the final one-figure, maximum-two-figure integrated manuscript strategy.

## Working Results text

Principal component analysis was used to determine whether the ten integrated panel scores formed a lower-dimensional systemic remodeling structure. After age and sex adjustment, the first principal component accounted for 35.6% of total score variance and showed uniformly positive contributions from all adaptive, humoral, innate-cytotoxic, and myeloid axes, supporting its interpretation as a broad multi-compartment remodeling gradient. The second component accounted for 15.8% and contrasted myeloid/granulocytic, monocyte/macrophage-like, and NK remodeling with CD8 differentiation, T/NK-interface, and regulatory-checkpoint remodeling. This two-axis structure remained after additional removal of the cancer-versus-healthy mean shift, with PC1 and PC2 explaining 30.5% and 16.6% of disease-adjusted variance, respectively. Immunotype accounted for 11.7% of the disease-adjusted multivariate score variance in permutation analysis. The adjusted immunotype centroids were ordered along PC1 from G1_Naive and G2_Primed through G3_Progressive to G4_Chronic and G5_Suppressive. However, G4 and G5 occupied opposing PC2 positions: G4 was shifted toward the CD8/T-NK pole, whereas G5 was shifted toward the myeloid/monocyte/NK pole. Separate training and validation analyses reproduced the variance contribution and loading structure of the first two components, and validation projection preserved the immunotype-centroid organization. Negative silhouette values indicated substantial individual-level overlap, suggesting that the immunotypes represent distinct but continuous remodeling configurations rather than non-overlapping classes in the compressed integrated-score space.

## Working Discussion text

The PCA results provide a systems-level bridge between the cross-panel convergence analysis and the original immunotype labels. A single broad remodeling component captured coordinated variation across all ten panels, indicating that adaptive, humoral, innate-cytotoxic, and myeloid changes share a common systemic dimension. This global axis was not sufficient to explain the complete architecture, however. A second, internally reproducible component separated CD8/T-NK-oriented remodeling from myeloid/monocyte/NK-oriented remodeling and thereby distinguished G4_Chronic from G5_Suppressive despite their similarly elevated global remodeling positions. The immunotype structure is therefore better represented as a shared remodeling gradient combined with alternative biological configurations than as a simple unidimensional sequence. At the individual level, the immunotype distributions remained substantially overlapping, emphasizing that the categories summarize dominant immune organizations rather than sharply isolated biological states. Because the PCA and immunotypes were derived from the same underlying SDY2583 resource, these results constitute dimensional resolution and internal stability assessment rather than external validation or evidence of temporal progression.
