# Cross-panel convergence: manuscript hold note

Status: **retain for later integration after all planned additional analyses are complete.**

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

## Figure plan — working hold

Proposed main figure title: **Cross-panel convergence defines a coordinated peripheral immune-remodeling architecture**

- **Panel A — Analysis schematic:** 850 participants, 66 scores, 10 panels, covariate residualization, and selection of one principal integrated score per panel.
- **Panel B — Principal-score heatmap:** 10 × 10 adjusted Spearman correlation matrix for the ten integrated scores, with rho values printed and FDR-significant cells marked.
- **Panel C — Prespecified convergence forest/dot plot:** CP7–CP24, CP7–CP28, CP24–CP28, CP8–CP25, CP10–CP16, CP10–CP23, CP16–CP23, and CP26–CP28 displayed for full, cancer-only, and healthy-only analyses. Bootstrap 95% confidence intervals should be added before finalization.
- **Panel D — Integrated-score network:** ten panel nodes; edges restricted to qFDR < 0.05, with edge width proportional to |rho| and edge style indicating direction. Nodes should be organized by CD8/T-cell, CD4/regulatory, B-cell/humoral, NK/T–NK, and myeloid/APC domains.

Suggested supplementary material:
- Supplementary full 66 × 66 clustered score heatmap;
- panel-pair median-correlation matrix;
- full pairwise correlation table with n, rho, P, and qFDR;
- robustness table comparing full, cancer-only, and healthy-only estimates;
- bootstrap confidence intervals and alternative partial-correlation sensitivity analysis.

## Discussion interpretation — working manuscript text

The cross-panel analysis extended the panel-specific results by demonstrating that the major remodeling scores were not statistically isolated. The strongest associations occurred between independently designed panels interrogating overlapping biological domains, including CD8 differentiation and the T/NK interface across CP7, CP24, and CP28; CD4 helper/regulatory and regulatory-checkpoint remodeling across CP8 and CP25; and myeloid, APC/DC-like, and monocyte/macrophage-like remodeling across CP10, CP16, and CP23. Importantly, these relationships persisted after adjustment for disease group in the full cohort and remained directionally concordant within both the cancer and healthy strata. The resulting pattern therefore cannot be explained solely by parallel mean shifts between patients with cancer and healthy controls. Instead, it supports a coordinated systemic architecture in which adaptive, innate-cytotoxic, humoral, and myeloid remodeling axes vary together across individuals. This convergence should not be interpreted as evidence of a single unidimensional immune state or a common causal mechanism, because shared biological regulation, technical covariance, and partial overlap among composite components may all contribute. Nevertheless, the reproducible cross-panel organization suggests that multi-compartment immune-remodeling summaries may provide a more stable translational framework than isolated marker-level abnormalities.
