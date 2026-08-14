# SDY2583 pan-cancer immune remodeling

Reproducible code and aggregate validation outputs for the secondary multi-panel
reanalysis of ImmPort **SDY2583 — Blood Immunotypes**. The project covers ten
flow-cytometry panels (CP7, CP8, CP10, CP16, CP22, CP23, CP24, CP25, CP26, and
CP28) and the integrated analyses used in the associated manuscript,
*Convergent Peripheral Immune Remodeling Across Lymphoid and Myeloid
Compartments in Cancer*.

## Release status

**Version 1.1.0 — locked internal validation release**

Version 1.1.0 retains the complete v1.0.0 reproducibility package and adds a
training-frozen held-out validation layer using the original 503-participant
TRAINING and 347-participant VALIDATION assignment.

This repository contains:

- all ten panel-specific downstream analysis pipelines;
- the final cross-panel convergence workflow;
- mapping of the ten principal scores to the five original SDY2583 immunotypes;
- disease-adjusted PCA and clustering-stability analyses;
- bounded-percentage model sensitivity analyses;
- same-sex nearest-age matching and matched-pair inference;
- Table 1 and panel-availability generation;
- 2,000-resample bootstrap component-stability analysis;
- the locked v6 code for the integrated Figure 7;
- **locked internal held-out validation of the 66-score framework**;
- curated aggregate results supporting the manuscript.

Raw FCS files, participant-level metadata, participant-level feature/score
matrices, matched subject identifiers, serialized workspaces, individual PCA
coordinates, and generated participant-level outputs are not distributed.

### Validation layers

The repository distinguishes two different uses of the word validation.

**Downstream archive validation** tests whether archived panel-level analytical
outputs can be independently reconstructed from archived post-extraction tables.
All ten manuscript panels completed this layer. A fresh raw-FCS-to-results rerun
was not performed; compensation, transformation, primary threshold regeneration
and regeneration of participant-level feature tables remain outside that scope.

**Locked internal held-out validation** is the v1.1.0 analysis. The original
503/347 TRAINING/VALIDATION assignment is retained. Component orientation,
centering/scaling and PCA fitting are learned only in TRAINING and then frozen
before application to VALIDATION. Score membership remains fixed from the
existing 66-score registry.

See [`LOCKED_INTERNAL_VALIDATION.md`](LOCKED_INTERNAL_VALIDATION.md),
[`VALIDATION_STATUS.md`](VALIDATION_STATUS.md),
[`PROVENANCE.md`](PROVENANCE.md), and
[`validation/ALL_PANEL_DOWNSTREAM_VALIDATION_SUMMARY.md`](validation/ALL_PANEL_DOWNSTREAM_VALIDATION_SUMMARY.md).

## Study overview

The local master matrix contains 850 participants:

- 408 healthy controls;
- 442 patients with cancer;
- 503 participants in the original training cohort;
- 347 participants in the original validation cohort.

The full integrated catalogue contains 66 composite scores. The main cross-panel
narrative is based on ten prespecified principal scores representing eight
biological axes:

1. CD8 differentiation/checkpoint remodeling — CP7 and CP24
2. CD4 helper/regulatory-like remodeling — CP8 and CP25
3. B-cell/humoral repatterning — CP22
4. NK-cell remodeling/attenuation — CP26
5. T/NK-interface remodeling — CP28
6. Myeloid/granulocytic remodeling — CP10
7. APC/DC-like remodeling — CP16
8. Monocyte/macrophage-like remodeling — CP23

The CP22 immunoglobulin-isotype score represents denominator-specific B-cell
isotype architecture/repatterning. It is not a measure of total IgG, total IgA,
or circulating immunoglobulin concentration.

## Locked internal held-out results

With direction, centering/scaling and hierarchical score rebuilding frozen from
TRAINING, **9/10 principal scores replicated** in VALIDATION with the same
direction and BH-FDR < 0.05. CP24 did not replicate (validation beta -0.022,
95% CI -0.224 to 0.180; BH-q 0.830) and is retained transparently as a
non-replicated axis.

Across all 66 locked scores:

- 60/66 retained the TRAINING direction in VALIDATION;
- 44/66 were significant in VALIDATION after 66-test BH correction.

A PCA fitted only in TRAINING explained 39.43% of variance on PC1. Projection
into VALIDATION without refitting retained a strong cancer-associated PC1 effect
(beta 1.569, 95% CI 1.232 to 1.905; p = 6.01e-20).

For the 45 principal-score panel pairs, 40/45 correlations retained direction
and 35/45 were validation-significant after BH correction. TRAINING-versus-
VALIDATION correlation-effect concordance was Pearson r = 0.700 (p = 8.74e-8).

These findings are **internal held-out validation, not external validation**.
Score membership predates the locked split analysis and remains an explicit
limitation.

## Original final integrated results

The v1.0.0 integrated evidence remains available under [`results/final/`](results/final/):

- 38/45 age-, sex-, and disease-adjusted Spearman relationships among the ten
  principal scores were FDR-significant;
- 39/45 were FDR-significant in Pearson sensitivity analysis;
- all 38 Spearman-significant relationships retained direction in both the
  cancer-only and healthy-only analyses;
- all ten principal scores differed across the five original immunotypes in the
  full, training, validation, cancer-only, and treatment-adjusted cancer-only analyses;
- disease-adjusted PC1 explained 30.5% and PC2 explained 16.6% of variance;
- maximum non-GMM silhouette was approximately 0.207;
- all eight representative bounded percentage outcomes retained direction and
  FDR significance across linear-HC3, empirical-logit, beta-regression, and
  age-matched paired analyses;
- 5-year matching produced 265 pairs and 10-year matching produced 273 pairs;
- across 2,000 disease-stratified bootstrap resamples, 62/66 scores retained
  direction in at least 95% of resamples and all ten principal score intervals
  excluded zero.

## Data access

Obtain the source data directly from ImmPort:

- Study: **SDY2583 — Blood Immunotypes**
- Dataset DOI: <https://doi.org/10.21430/M3A0B9RD5T>
- Primary publication: Dyikanov et al., *Cancer Cell* (2024),
  <https://doi.org/10.1016/j.ccell.2024.04.008>

Use remains subject to ImmPort terms. See [`DATA_ACCESS.md`](DATA_ACCESS.md).

## Repository structure

```text
R/
  integration/                     integrated and locked-validation workflows
  figures/                         final integrated Figure 7 code
  shared/                          portable paths and workflow engines
  panels/CP*/                      panel-specific SAFE/RECONSTRUCTED source
config/
  paths.example.R
  final_analysis_paths.example.env
  locked_validation/               public 66-score membership registries
scripts/
  01_run_all_panels.R
  10-19 panel runners
  20-28 panel validators
  30_validate_all_panels.R
  40_run_final_integrated_release.sh
  41_run_integrated_systems_analysis.py
  42_run_frequency_matching_table1.py
  43_run_bootstrap_component_stability.py
  44_run_locked_internal_validation.R
results/final/                     public aggregate-only final results
validation/                        downstream archive-validation evidence
analysis_notes/                    locked scientific interpretation records
```

## Software setup

### R

```bash
Rscript scripts/00_install_dependencies.R
```

Copy `config/paths.example.R` to `config/paths.R` and configure the local
ImmPort/FCS/output paths. The local file is ignored by Git.

### Python

Python 3.10 or later is recommended.

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements-final.txt
```

Copy `config/final_analysis_paths.example.env` to a local file, edit the paths,
and source it. Do not commit participant-level paths or data.

## Running the panel workflows

```bash
Rscript scripts/01_run_all_panels.R
```

For selected panels:

```bash
SDY2583_PANELS=CP7,CP8,CP25 Rscript scripts/01_run_all_panels.R
```

## Running the final integrated analyses

```bash
bash scripts/40_run_final_integrated_release.sh
```

## Running locked internal validation

Set a common local archive root or explicit per-file environment variables:

```bash
export SDY2583_LOCKED_INPUT_ROOT="/path/to/local/SDY2583/analysis/archive"
Rscript scripts/44_run_locked_internal_validation.R
```

The locked-validation script writes aggregate outputs only. See
[`LOCKED_INTERNAL_VALIDATION.md`](LOCKED_INTERNAL_VALIDATION.md) for the exact
methodological and interpretation boundary.

## Public-output boundary

The `results/final/` directory contains aggregate tables only. Public outputs do
not include participant identifiers, participant-level scores or predictions,
matched-pair IDs, individual PCA coordinates, raw event-level data, or FCS files.

## Static checks

```bash
Rscript scripts/02_static_source_audit.R
```

The same data-free audit runs through GitHub Actions.

## Citation and license

Cite this repository together with the ImmPort dataset and the primary SDY2583
publication. Machine-readable metadata are in [`CITATION.cff`](CITATION.cff).
Code is released under the [MIT License](LICENSE); ImmPort data are not covered
by that license.
