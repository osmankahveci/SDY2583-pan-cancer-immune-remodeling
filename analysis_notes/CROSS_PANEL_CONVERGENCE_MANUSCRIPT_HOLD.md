# Cross-panel convergence: manuscript hold note

Status: **retain for later integration after all planned additional analyses are complete.**

Figure strategy: **do not generate or finalize figures until every additional analysis is complete. The final manuscript should receive one integrative figure, or at most two figures, covering all new analyses. A network-centered synthesis is currently preferred, potentially combined with compact effect-size or robustness summaries. The earlier four-panel figure concept is only a component library, not the final figure plan.**

Input matrix: `SDY2583_integrated_clinical_immune_score_matrix_ALL10_with_CP23.csv`

Validated matrix dimensions:
- 850 participants
- 127 total columns
- 66 composite/integrated immune-remodeling scores across 10 panels

Primary analysis:
- one integrated score per panel;
- full cohort residualized for age, sex, and disease group;
- cancer-only and healthy-only analyses residualized for age and sex;
- pairwise Spearman correlations;
- Benjamini–Hochberg FDR correction across the 45 principal score pairs.

Key results:
- 39/45 principal integrated-score pairs were FDR-significant.
- CP24 integrated remodeling vs CP28 integrated T/NK-interface remodeling: rho = 0.563, q = 2.67e-69; cancer rho = 0.612; healthy rho = 0.480.
- CP8 integrated CD4 helper/regulatory remodeling vs CP25 integrated CD4 regulatory-checkpoint remodeling: rho = 0.487, q = 1.32e-49; cancer rho = 0.539; healthy rho = 0.425.
- CP7 integrated checkpoint remodeling vs CP28 integrated T/NK-interface remodeling: rho = 0.466, q = 5.80e-45; cancer rho = 0.492; healthy rho = 0.435.
- CP10 integrated myeloid/granulocytic remodeling vs CP23 integrated monocyte/macrophage-like remodeling: rho = 0.370, q = 2.21e-27; cancer rho = 0.427; healthy rho = 0.274.
- CP7 integrated checkpoint remodeling vs CP24 integrated remodeling: rho = 0.261, q = 6.04e-14; cancer rho = 0.289; healthy rho = 0.208.
- CP16 integrated APC/DC-like myeloid remodeling vs CP23 integrated monocyte/macrophage-like remodeling: rho = 0.239, q = 7.53e-12; cancer rho = 0.197; healthy rho = 0.270.
- CP26 integrated NK remodeling vs CP28 integrated T/NK-interface remodeling: rho = 0.256, q = 1.69e-13; cancer rho = 0.250; healthy rho = 0.286.
- CP10 integrated myeloid/granulocytic remodeling vs CP16 integrated APC/DC-like myeloid remodeling was weaker but retained FDR significance: rho = 0.098, q = 0.0060.

Interpretive boundary:
- Primary inference should rely on the 10 principal integrated scores because the full 66-score matrix contains overlapping and non-independent composites.
- Full 66-score results may be presented in supplementary material as evidence of broad connectivity, not as 1,951 independent biological tests.
- Correlations indicate convergent organization, not causality, common cellular identity, or a single unidimensional immune state.

## Results paragraph — working manuscript text

To determine whether the panel-specific findings represented isolated associations or a coordinated systemic architecture, cross-panel convergence was evaluated at the subject level. One principal integrated remodeling score was selected from each of the ten cytometry panels, and pairwise Spearman correlations were calculated after residualizing the scores for age, sex, and disease group in the full cohort. Thirty-nine of the 45 pairwise comparisons remained significant after Benjamini–Hochberg correction. The strongest convergence was observed between CP24 CD8 differentiation remodeling and CP28 T/NK-interface remodeling (rho = 0.563, qFDR = 2.67 × 10^-69), between CP8 CD4 helper/regulatory remodeling and CP25 CD4 regulatory-checkpoint remodeling (rho = 0.487, qFDR = 1.32 × 10^-49), and between CP7 checkpoint remodeling and CP28 T/NK-interface remodeling (rho = 0.466, qFDR = 5.80 × 10^-45). Myeloid convergence was also evident between CP10 myeloid/granulocytic remodeling and CP23 monocyte/macrophage-like remodeling (rho = 0.370, qFDR = 2.21 × 10^-27), whereas CP16 APC/DC-like remodeling correlated with CP23 monocyte/macrophage-like remodeling (rho = 0.239, qFDR = 7.53 × 10^-12). Additional within-domain associations included CP7–CP24 CD8 remodeling (rho = 0.261, qFDR = 6.04 × 10^-14) and CP26–CP28 NK/T–NK-interface remodeling (rho = 0.256, qFDR = 1.69 × 10^-13). All prespecified cross-panel associations retained the same positive direction in the cancer-only and healthy-only analyses after age and sex adjustment, indicating that the observed convergence was not attributable solely to the cancer-versus-healthy group separation.

## Figure strategy — deferred until all analyses are complete

Do not create the final figure now. After all additional analyses are complete, select only the most informative and nonredundant results for one integrated summary figure, or at most two figures.

Preferred current concept:
- a central immune-remodeling network linking the ten panel-level integrated scores or the final reduced biological axes;
- node organization by CD8/T-cell, CD4/regulatory, B-cell/humoral, NK/T–NK, and myeloid/APC/monocyte domains;
- edge width proportional to the final selected association magnitude;
- restrained annotation of the strongest validated relationships;
- compact peripheral elements, only if needed, showing robustness across full, cancer-only, and healthy-only analyses or showing the outputs of later clustering/transcriptomic/immunotype analyses.

Potential components to retain for later selection, not as mandatory separate panels:
- analysis schematic: 850 participants, 66 scores, 10 panels, covariate residualization, and selection of principal integrated axes;
- 10 × 10 principal-score heatmap;
- prespecified convergence effect-size plot with bootstrap 95% confidence intervals;
- integrated-score network;
- final patient-level phenotype, immunotype, transcriptomic, or clustering overlays arising from subsequent analyses.

The final design decision must be made only after reviewing all new analyses together. The objective is to avoid adding multiple fragmented figures and instead produce a single clear visual synthesis of the translational contribution.

Suggested supplementary material, to be reconsidered after all analyses:
- full 66 × 66 clustered score heatmap;
- panel-pair median-correlation matrix;
- full pairwise correlation table with n, rho, P, and qFDR;
- robustness table comparing full, cancer-only, and healthy-only estimates;
- bootstrap confidence intervals and alternative partial-correlation sensitivity analysis.

## Discussion interpretation — working manuscript text

The cross-panel analysis extended the panel-specific results by demonstrating that the major remodeling scores were not statistically isolated. The strongest associations occurred between independently designed panels interrogating overlapping biological domains, including CD8 differentiation and the T/NK interface across CP7, CP24, and CP28; CD4 helper/regulatory and regulatory-checkpoint remodeling across CP8 and CP25; and myeloid, APC/DC-like, and monocyte/macrophage-like remodeling across CP10, CP16, and CP23. Importantly, these relationships persisted after adjustment for disease group in the full cohort and remained directionally concordant within both the cancer and healthy strata. The resulting pattern therefore cannot be explained solely by parallel mean shifts between patients with cancer and healthy controls. Instead, it supports a coordinated systemic architecture in which adaptive, innate-cytotoxic, humoral, and myeloid remodeling axes vary together across individuals. This convergence should not be interpreted as evidence of a single unidimensional immune state or a common causal mechanism, because shared biological regulation, technical covariance, and partial overlap among composite components may all contribute. Nevertheless, the reproducible cross-panel organization suggests that multi-compartment immune-remodeling summaries may provide a more stable translational framework than isolated marker-level abnormalities.
