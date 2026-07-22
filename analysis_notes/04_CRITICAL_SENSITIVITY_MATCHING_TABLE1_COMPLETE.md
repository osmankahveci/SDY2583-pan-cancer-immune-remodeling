# Critical sensitivity analyses, age-matching diagnostics, Table 1, and pan-cancer heterogeneity

Status: **completed; aggregate results only. Participant-level raw feature tables and matched-pair identifiers are not stored in GitHub.**

## Completion summary

- Eight representative raw percentage outcomes were prespecified, one for each principal biological axis.
- Each was analyzed using:
  1. raw-percentage ordinary least squares with HC3 standard errors;
  2. empirical-logit ordinary least squares with HC3 standard errors;
  3. beta regression with a logit link and HC3 covariance.
- All 8/8 outcomes retained the same cancer-versus-healthy direction and remained FDR-significant under all three model families.
- Same-sex nearest-age matching without replacement was reconstructed with 5-year and 10-year calipers.
- Five-year matching produced 265 pairs; absolute age SMD decreased from 0.919 before matching to 0.066, with a mean absolute within-pair difference of 0.99 years.
- Ten-year matching produced 273 pairs; absolute age SMD was 0.064, with a mean absolute within-pair difference of 1.29 years.
- Exact sex matching reduced the sex SMD from 0.072 to 0.
- All eight representative outcomes retained direction and FDR significance in paired analyses under both calipers.
- The complete cohort/treatment Table 1 was generated, including cancer subgroups, treatment-status availability, treatment modalities, treatment-line availability, and panel-level FCS counts.
- CP16 contained 847 FCS files from 847 unique subjects; the other nine panels contained 850 files/subjects.

## Methods — proportion-scale sensitivity analyses

To evaluate whether the principal frequency-based findings were sensitive to the bounded and potentially heteroskedastic distribution of percentage outcomes, one biologically representative raw frequency measure was prespecified for each of the eight principal remodeling axes. Each outcome was analyzed using three complementary models: (i) ordinary least-squares regression on the original percentage scale with heteroskedasticity-consistent HC3 standard errors; (ii) ordinary least-squares regression after empirical-logit transformation, also using HC3 standard errors; and (iii) beta regression with a logit link and HC3 covariance. Values at 0 or 100% were mapped to the open unit interval using the Smithson–Verkuilen adjustment before empirical-logit and beta-regression analyses. All models included cancer status, age, and sex. Benjamini–Hochberg correction was applied across the eight prespecified outcomes within each model family.

## Methods — age-matching diagnostics and paired sensitivity analyses

Age sensitivity was evaluated using greedy nearest-age matching without replacement, restricted to participants of the same recorded binary sex. Cancer participants were processed in ascending age order and matched to the closest available healthy control within prespecified 5-year and 10-year calipers. Covariate balance was quantified using standardized mean differences before and after matching. Absolute within-pair age differences and their distribution were summarized. Disease-associated effects in the matched samples were estimated from within-pair cancer-minus-healthy differences; standard errors and confidence intervals were calculated from the between-pair distribution of these differences, thereby retaining the matched-pair structure. False-discovery-rate correction was applied across the eight representative outcomes for each caliper analysis.

## Results — bounded-outcome robustness

All eight prespecified representative frequency outcomes retained the same cancer-versus-healthy direction across the raw-percentage HC3 model, empirical-logit HC3 model, and beta-regression model, and all remained significant after false-discovery-rate correction within each model family. On the original percentage scale, cancer was associated with higher TEMRA-like CD8 representation in CP24 (+4.22 percentage points), increased CP8 CD25+IL7RA-low regulatory-like representation (+0.68 points), increased switched-memory-like representation within B cells in CP22 (+5.86 points), reduced CP26 NK-like representation (−4.56 points), increased CP28 CD8 representation within CD3+ cells (+5.93 points), increased CP10 CD13+CD66b+ granulocyte-like representation (+15.61 points), reduced CP16 HLA-DR+ APC-core representation (−4.01 points), and increased CP23 CD45+ dump-low representation (+2.75 points). The complete concordance across modeling scales indicates that the main directional conclusions were not artifacts of applying linear regression directly to bounded percentage outcomes.

## Results — matching balance and paired outcomes

Before matching, the standardized mean difference for age was 0.919. Five-year same-sex matching produced 265 pairs and reduced the absolute age standardized mean difference to 0.066; the mean absolute within-pair age difference was 0.99 years, the median was 0.0 years, and the maximum was 5 years. Ten-year matching produced 273 pairs, with an absolute age standardized mean difference of 0.064 and a mean absolute difference of 1.29 years. Exact matching on sex reduced the sex standardized mean difference from 0.072 to 0. All eight representative outcomes retained their primary direction and remained FDR-significant in paired analyses under both calipers.

## Discussion — pan-cancer heterogeneity and interpretive limits

The pooled pan-cancer coefficients should be interpreted as estimates of a shared disease-associated remodeling axis rather than evidence that every malignancy exhibits an identical effect size or cellular configuration. Quantitative heterogeneity across cancer types is biologically expected. In particular, the comparatively lower CP10 integrated myeloid/granulocytic, granulocyte-like, monocyte-like HLA-DR phenotype, and myeloid-to-lymphoid balance scores observed in breast cancer relative to the other malignancies should be considered hypothesis-generating rather than definitive evidence of a breast-cancer-specific immune state. Adjustment for cancer subgroup and recorded treatment status reduced, but could not eliminate, residual confounding related to disease stage, tumor burden, treatment timing, incompletely observed treatment line, and other clinical factors that were not uniformly available. Accordingly, the present study does not establish treatment response, resistance, prognosis, survival associations, or temporal progression of the identified immune-remodeling profiles.

## Table 1 reporting boundaries

- Total cohort: 850 participants; 408 healthy controls and 442 cancer patients.
- Clean age available for 832 participants; 18 invalid or missing ages were excluded from age-adjusted models.
- Treatment metadata were available for 335/442 cancer participants.
- Treatment modalities were non-mutually exclusive and use the treatment-known denominator.
- Treatment line was available for 118/442 cancer participants and unavailable for 324.
- The pooled pan-cancer coefficient is not interpreted as a uniform cancer-type effect.
- The CP22 representative switched-memory-like frequency is one component of humoral repatterning and must not be interpreted as total IgG or total IgA abundance.
