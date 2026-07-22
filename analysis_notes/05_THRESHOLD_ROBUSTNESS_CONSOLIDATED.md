# Consolidated threshold and analysis-set robustness

Status: **completed; aggregate panel-level summary only.**

## Interpretation

Threshold sensitivity was consolidated to distinguish stability of the principal biological axes from instability of individual threshold-defined features. The principal integrated scores for CP7, CP8, CP10, CP16, CP22, CP23, CP25, CP26, and CP28 retained the expected cancer-versus-healthy direction and FDR significance across the main, permissive, and stringent threshold configurations. CP24 used event-count analysis-set sensitivity rather than the same marker-threshold-shift framework in its final package and is therefore reported separately.

| Panel | Analysis type | Targeted variables | Direction + FDR preserved | Direction preserved; FDR not all sets | Direction not preserved | Composite robustness | Principal integrated score |
|---|---|---:|---:|---:|---:|---|---|
| CP7 | Main/permissive/stringent thresholds | 11 | 11 | 0 | 0 | 4/4 | Robust |
| CP8 | Main/permissive/stringent thresholds | 28 | 28 | 0 | 0 | 6/6 | Robust |
| CP10 | Main/permissive/stringent thresholds | 33 | 29 | 3 | 1 | 5/6 | Robust |
| CP16 | Main/permissive/stringent thresholds | 60 | 33 | 14 | 13 | 4/6 | Robust |
| CP22 | Main/permissive/stringent thresholds | 55 | 39 | 6 | 10 | 6/7 | Robust |
| CP23 | Main/permissive/stringent thresholds | 20 | 16 | 4 | 0 | 5/5 | Robust |
| CP24 | Event-count analysis sets: all, >=500, >=1000 | 16 | 8 | 8 | 0 | Not reported in this summary | Not directly tested in this summary |
| CP25 | Main/permissive/stringent thresholds | 73 | 48 | 23 | 2 | 9/9 | Robust |
| CP26 | Main/permissive/stringent thresholds | 20 | 10 | 7 | 3 | Not separately classified | Robust |
| CP28 | Main/permissive/stringent thresholds | 80 | 35 | 34 | 11 | Not separately classified | Robust |

## Manuscript wording

Threshold sensitivity was summarized at the panel level to distinguish the stability of the principal biological axes from the behavior of individual threshold-defined features. The principal integrated scores for CP7, CP8, CP10, CP16, CP22, CP23, CP25, CP26, and CP28 retained the expected cancer-versus-healthy direction and false-discovery-rate significance across the main, permissive, and stringent threshold configurations. Robustness was strongest for CP7 and CP8, in which all targeted outcomes preserved both direction and FDR significance. Across the remaining panels, individual low-frequency or narrowly defined subpopulations were less stable than the integrated composites, whereas the principal integrated scores remained robust. CP24 used event-count analysis-set sensitivity rather than the same marker-threshold-shift framework in its final package; all 16 selected outcomes retained direction across the full, >=500-event, and >=1000-event analysis sets, and eight retained module-level FDR significance throughout. These findings support using the eight principal biological axes in the main text while retaining the full feature-level robustness catalogue in the Supplementary Results.

## Required guardrails

- Feature-level instability must not be hidden; it is expected to be most prominent for rare or narrowly gated subpopulations.
- The main inference rests on prespecified integrated axes, not on every individual threshold-defined feature.
- CP24 is not directly comparable to the marker-threshold analyses because its final package varied event-count eligibility.
- CP22 remains a humoral/B-cell repatterning construct and must not be described as total IgG or total IgA abundance.
