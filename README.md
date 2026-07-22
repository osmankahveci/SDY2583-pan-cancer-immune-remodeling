# SDY2583 pan-cancer immune remodeling

Reproducible code and aggregate validation outputs for the secondary multi-panel
reanalysis of ImmPort **SDY2583 — Blood Immunotypes**. The project covers ten
flow-cytometry panels (CP7, CP8, CP10, CP16, CP22, CP23, CP24, CP25, CP26, and
CP28) and the final integrated analyses used in the associated manuscript,
*Convergent Peripheral Immune Remodeling Across Lymphoid and Myeloid
Compartments in Cancer*.

## Release status

**Version 1.0.0 — final public reproducibility release**

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
- curated aggregate results supporting the manuscript.

Raw FCS files, participant-level metadata, participant-level feature/score
matrices, matched subject identifiers, serialized workspaces, and generated
participant-level outputs are not distributed.

### Accepted validation boundary

All ten manuscript panels completed independent **downstream archive
validation** from archived participant-level feature, metadata, score, matching,
interaction, and threshold-sensitivity tables. A fresh raw-FCS-to-results rerun
was not performed. Compensation, logicle transformation, primary threshold
regeneration, and regeneration of participant-level feature tables therefore
remain outside the accepted validation scope.

Accurate release description:

> Provenance-tracked all-panel reproducibility release with downstream analyses independently archive-validated and final integrated analyses numerically executed on the local ALL10 matrix. Fresh raw-FCS end-to-end re-execution was not performed.

Source provenance remains explicit:

- **`SAFE`** — recovered archived source with portability/non-interactive edits only.
- **`RECONSTRUCTED`** — source rebuilt from archived thresholds, schemas,
  composite definitions, Methods/Results records, and numerical output tables.

See [`VALIDATION_STATUS.md`](VALIDATION_STATUS.md),
[`PROVENANCE.md`](PROVENANCE.md), and
[`validation/ALL_PANEL_DOWNSTREAM_VALIDATION_SUMMARY.md`](validation/ALL_PANEL_DOWNSTREAM_VALIDATION_SUMMARY.md).

## Study overview

The local master matrix contains 850 participants:

- 408 healthy controls;
- 442 patients with cancer;
- 503 participants in the original training cohort;
- 347 participants in the original validation cohort.

The full integrated catalogue contains 66 directionally oriented composite
scores. The main cross-panel narrative is based on ten prespecified principal
scores representing eight biological axes:

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

## Final integrated results

Curated aggregate tables are under [`results/final/`](results/final/).
Manuscript-level checks include:

- 38/45 age-, sex-, and disease-adjusted Spearman relationships among the ten
  principal scores were FDR-significant;
- 39/45 were FDR-significant in Pearson sensitivity analysis;
- all 38 Spearman-significant relationships retained direction in both the
  cancer-only and healthy-only analyses;
- all ten principal scores differed across the five original immunotypes in the
  full, training, validation, cancer-only, and treatment-adjusted cancer-only
  analyses;
- disease-adjusted PC1 explained 30.5% and PC2 explained 16.6% of variance;
- the integrated space showed non-random structure but weak discrete-cluster
  separation (maximum non-GMM silhouette approximately 0.207);
- all eight representative bounded percentage outcomes retained direction and
  FDR significance across linear-HC3, empirical-logit, beta-regression, and
  age-matched paired analyses;
- 5-year matching produced 265 pairs and 10-year matching produced 273 pairs;
- across 2,000 disease-stratified bootstrap resamples, 62/66 scores retained
  direction in at least 95% of resamples and all ten principal score intervals
  excluded zero.

These are internal reproducibility and convergent-validity results, not an
independent external validation cohort.

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
  integration/                     metadata and cross-panel builders
  figures/                         final integrated Figure 7 code
  shared/                          portable paths and workflow engines
  panels/CP*/                      panel-specific SAFE/RECONSTRUCTED source
config/
  paths.example.R
  final_analysis_paths.example.env
scripts/
  01_run_all_panels.R
  10-19 panel runners
  20-28 panel validators
  30_validate_all_panels.R
  40_run_final_integrated_release.sh
  41_run_integrated_systems_analysis.py
  42_run_frequency_matching_table1.py
  43_run_bootstrap_component_stability.py
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

Run all ten panel pipelines:

```bash
Rscript scripts/01_run_all_panels.R
```

Run selected panels:

```bash
SDY2583_PANELS=CP7,CP8,CP25 Rscript scripts/01_run_all_panels.R
```

The all-panel release remains subject to the validation boundary above; it must
not be described as a newly completed raw-FCS end-to-end execution.

## Running the final integrated analyses

After setting the local input paths:

```bash
bash scripts/40_run_final_integrated_release.sh
```

The runner executes:

1. integrated convergence, immunotype, PCA, and clustering analyses;
2. bounded-outcome models, matching diagnostics, and Table 1 generation;
3. the 2,000-resample bootstrap stability analysis;
4. the v6 integrated Figure 7 workflow when `Rscript` is available.

Individual commands are documented in the scripts and in
[`docs/FINAL_INTEGRATED_ANALYSIS_RELEASE.md`](docs/FINAL_INTEGRATED_ANALYSIS_RELEASE.md).

## Public-output boundary

The `results/final/` directory contains aggregate tables only. Public outputs do
not include:

- participant identifiers;
- participant-level scores or predictions;
- matched-pair subject IDs;
- individual PCA coordinates;
- raw event-level data;
- FCS files.

## Static checks

A data-free audit parses R source, checks runner references and panel coverage,
rejects author-specific executable paths and interactive file selection, and
screens for accidentally committed data:

```bash
Rscript scripts/02_static_source_audit.R
```

The same audit runs through GitHub Actions.

## Citation and license

Cite this repository together with the ImmPort dataset and the primary SDY2583
publication. Machine-readable metadata are in [`CITATION.cff`](CITATION.cff).
Code is released under the [MIT License](LICENSE); ImmPort data are not covered
by that license.
