# Final additional-analysis integration plan

Status: **decision locked. Integrate into the manuscript only after the main text is opened for revision. Do not generate figures before the final integrated design is approved.**

Scope decision:
- Include four completed analysis packages: cross-panel convergence, original-immunotype biological resolution, global PCA structure, and unsupervised clustering/stability.
- Exclude transcriptomic integration and predictive/random-forest analyses.
- Add one new integrative main figure; use a second only if the first becomes unreadable.

## 1. Main Results priorities

### A. Cross-panel convergence — primary new systems-level result

Role: establishes that the ten panel-level remodeling axes are coordinated rather than isolated.

Main-text content:
- 10 principal integrated scores; 45 pairwise relationships.
- 39/45 FDR-significant after age, sex, and disease-group adjustment.
- Emphasize strongest biologically coherent links: CP24–CP28, CP8–CP25, CP7–CP28, CP10–CP23, CP16–CP23, and CP26–CP28.
- State that prespecified associations retained positive direction in cancer-only and healthy-only analyses.

Placement: first additional-analysis subsection after the panel-specific Results and before immunotype/PCA synthesis.

Approximate manuscript allocation: one focused Results paragraph plus figure callout.

### B. Original five-immunotype biological resolution — primary translational interpretation result

Role: converts the published G1–G5 labels into interpretable multi-panel biological architectures.

Main-text content:
- Immunotype effect significant for all 10 principal integrated scores after age, sex, and disease-group adjustment.
- Training and validation subsets independently reproduce all 10 omnibus effects.
- Cancer-only and treatment-adjusted sensitivity remains significant.
- Biological profiles:
  - G1_Naive: low-remodeling reference.
  - G2_Primed: early CD4 helper/regulatory shift.
  - G3_Progressive: intermediate CD4 regulatory, humoral, and myeloid remodeling.
  - G4_Chronic: strongest CD8 differentiation and T/NK-interface configuration.
  - G5_Suppressive: strongest myeloid/granulocytic, monocyte/macrophage-like, APC/DC-like, and NK configuration.
- Emphasize the distinction between G4 and G5 as two advanced but biologically different remodeling states.

Placement: immediately after cross-panel convergence.

Approximate manuscript allocation: one substantial Results paragraph plus one short robustness sentence.

Interpretive boundary: biological resolution/internal replication, not independent validation of the original immunotypes.

### C. PCA/global remodeling structure — primary dimensional-synthesis result

Role: provides the compact systems-level model connecting convergence and immunotypes.

Main-text content:
- PC1 represents broad multi-compartment remodeling burden.
- PC2 separates myeloid/monocyte/NK remodeling from CD8/T-NK/regulatory-checkpoint remodeling.
- Disease-adjusted PCA preserves the architecture.
- G1/G2 occupy low PC1, G3 is intermediate, and G4/G5 occupy high PC1 but opposing PC2 poles.
- Training-validation loading and centroid stability support internal reproducibility.
- Immunotypes are distinct in mean position but substantially overlapping at the individual level.

Placement: directly after immunotype resolution, functioning as the synthesis subsection.

Approximate manuscript allocation: one Results paragraph.

Interpretive boundary: adjusted centroid gradient, not longitudinal progression or prognosis.

## 2. Supplementary Results and tables

### Cross-panel convergence supplementary material
- Full 66 × 66 correlation matrix or clustered heatmap.
- Full pairwise table with n, rho, P, and qFDR.
- Panel-pair median-correlation summary.
- Full/cancer-only/healthy-only robustness table.
- Optional bootstrap confidence intervals for the prespecified panel pairs.

### Immunotype supplementary material
- Full 66-score omnibus screen.
- Adjusted marginal means for all five immunotypes and ten principal scores.
- All pairwise immunotype contrasts.
- Training-validation replication table.
- Cancer-only and treatment-adjusted sensitivity tables.
- Demographic/disease composition of immunotypes to document confounding structure.

### PCA supplementary material
- Complete loading table for all components.
- Variance-explained table.
- Training/validation loading similarity and projected-centroid stability.
- Multivariate permutation variance partitioning.
- Silhouette statistics documenting overlap.

## 3. Methodological-support result

### Unsupervised clustering and stability

Role: protects the manuscript from overclaiming and justifies representing the data as continuous dimensions rather than inventing new subtypes.

Main-text handling:
- Include no dedicated main Results subsection unless required by reviewers.
- Add one concise sentence at the end of the PCA subsection:
  `Unsupervised clustering did not support additional sharply discrete and reproducible remodeling subtypes; the only highly stable solution primarily dichotomized the continuous global-remodeling component.`

Supplementary handling:
- Methods: k-means, Ward clustering, GMM information criteria, Hopkins statistic, subsampling stability, and training-validation reproducibility.
- Results: internal-validity metrics, ARI/Jaccard results, and the finding that k=2 is nearly fully recoverable from PC1.
- State explicitly that k>=3 did not meet reproducibility requirements for subtype claims.

Figure policy: no clustering-specific figure.

## 4. Final main-figure strategy

Preferred final design: **one integrated systems figure**, network-centered but not network-only.

Working title:
`Coordinated cross-panel remodeling resolves the five peripheral immunotypes along shared and divergent systemic immune axes`

Recommended structure:

### Central layer — biological-domain network
- Ten principal panel scores or a reduced set of biological-domain nodes.
- Nodes organized into CD8/T-cell, CD4/regulatory, B-cell/humoral, NK/T-NK, and myeloid/APC/monocyte domains.
- Edges restricted to strong, FDR-supported, biologically interpretable relationships rather than displaying all 39 significant pairs.
- Edge width proportional to adjusted rho.
- The network should communicate coordinated multi-compartment architecture.

### Immunotype layer
- Encode the adjusted G1–G5 profile around or beneath the network using compact standardized profile glyphs, mini-bars, or a five-row domain heat strip.
- Make the G4 CD8/T-NK pole and G5 myeloid/monocyte/NK pole visually unmistakable.
- Avoid ten separate radar plots or dense marker labels.

### Dimensional layer
- Add a compact PC1–PC2 conceptual axis or centroid inset:
  - PC1: low to high global remodeling burden.
  - PC2: CD8/T-NK pole versus myeloid/monocyte/NK pole.
- Plot only immunotype centroids with confidence regions or directional labels, not all 829 participants, unless readability remains excellent.

### Robustness notation
- Use a small symbol or annotation to indicate training-validation stability and cancer/healthy directional concordance.
- Do not add a separate forest plot unless the network-plus-centroid design cannot display robustness clearly.

Maximum two-figure contingency:
- Figure 1: network plus immunotype profiles.
- Figure 2: PC1–PC2 centroid map with training-validation replication.
- This split should be used only if a single figure becomes visually overloaded.

## 5. Main Discussion integration

Discussion should synthesize three claims:

1. Independent cytometry panels converge on a coordinated systemic remodeling architecture.
2. The original immunotypes can be biologically resolved into distinct multi-compartment configurations, especially the divergence between G4_Chronic and G5_Suppressive.
3. The architecture is best represented as a global remodeling gradient plus alternative biological configurations, not as newly discovered sharply discrete subtypes.

Required cautions:
- Same-resource dimensional annotation, not external validation.
- No causal, prognostic, predictive, or longitudinal claims.
- Composite-score overlap and potential shared technical covariance should be acknowledged.
- G1-to-G5 ordering on PC1 is a centroid gradient, not temporal progression.

## 6. Final hierarchy of evidentiary weight

1. Cross-panel convergence — primary mechanistic/systemic evidence.
2. Immunotype biological resolution — primary translational interpretation.
3. PCA — primary dimensional synthesis and internal stability evidence.
4. Clustering — negative model-selection and anti-overclaiming support.

This hierarchy should control word count, figure emphasis, abstract inclusion, and Discussion ordering.
