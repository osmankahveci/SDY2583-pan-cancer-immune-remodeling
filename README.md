# SDY2583 pan-cancer immune remodeling

Reproducible analysis code for the secondary multi-panel reanalysis of ImmPort
SDY2583 *Blood Immunotypes*. The repository covers CP7, CP8, CP10, CP16, CP22,
CP23, CP24, CP25, CP26, and CP28, including demographic adjustment, composite
phenotypes, robustness analyses, clinical annotation, and figure generation.

## Release status

This is a provenance-tracked all-panel reproducibility release. Raw FCS files,
participant-level metadata and matrices, serialized workspaces, and generated
outputs are excluded from Git.

All ten manuscript panels have completed independent downstream archive
validation from archived participant-level feature, metadata, score, matching,
interaction, and threshold-sensitivity tables. Numerical agreement is at
floating-point/machine-precision scale, with exact subject-pair reproduction
where archived matching tables were available.

A fresh raw-FCS end-to-end rerun was not performed. FCS inventory,
compensation, fixed-logicle transformation, threshold gating, and regeneration
of participant-level feature tables are explicitly outside the accepted
validation scope of this release.

Accurate release description:

> Provenance-tracked all-panel reproducibility release with all ten panel downstream analyses independently archive-validated. Raw-FCS end-to-end re-execution was not performed and is outside the accepted validation scope of this release.

Two source classes are kept separate:

- **`SAFE`**: original analysis source recovered from the project archive and
  modified only for portable paths and non-interactive execution.
- **`RECONSTRUCTED`**: source rebuilt where the original script was not visible,
  using archived thresholds, marker maps, feature schemas, composite
  definitions, Methods/Results records, and numerical output tables.

Do not remove the provenance suffixes. See [VALIDATION_STATUS.md](VALIDATION_STATUS.md)
and [validation/ALL_PANEL_DOWNSTREAM_VALIDATION_SUMMARY.md](validation/ALL_PANEL_DOWNSTREAM_VALIDATION_SUMMARY.md)
for the accepted validation boundary and panel-level evidence.

## Panel coverage

| Panel | Included source | Status |
|---|---|---|
| CP7 | Reconstructed Steps 1–7B | Downstream archive-validated after score, matching, model, and threshold-family corrections. |
| CP8 | Reconstructed Steps 1–7B | Downstream archive-validated after matching, age-family, FDR, and threshold-family corrections. |
| CP10 | Reconstructed Steps 1/1B + recovered Steps 2–8 | Downstream archive-validated; recovered downstream algorithms retained. |
| CP16 | Recovered Steps 1–7B | Downstream archive-validated; recovered sequence retained. |
| CP22 | Recovered Steps 1–7B | Downstream archive-validated; recovered sequence retained. |
| CP23 | Recovered Steps 1–7B, including Step 6B | Downstream archive-validated; recovered sequence retained. |
| CP24 | Recovered pilot codebook + reconstructed full-cohort Steps 1–7B | Downstream archive-validated after outcome-family, PD-1, model-count, and hierarchical-score corrections. |
| CP25 | Reconstructed Steps 1–7B | Downstream archive-validated after age-family, matching, interaction, and exact threshold-family corrections. |
| CP26 | Reconstructed Steps 1–7B | Downstream archive-validated after score-universe and integrated-score corrections; includes BV510-A dump fallback. |
| CP28 | Reconstructed Steps 1–7B | Downstream archive-validated after duplicate-module, event-QC, age-family, matching, and threshold-family corrections. |

See [PROVENANCE.md](PROVENANCE.md) for panel-level source evidence and
limitations.

## Data access

Source data are not redistributed. Obtain **SDY2583 — Blood Immunotypes** from
ImmPort. Dataset DOI: <https://doi.org/10.21430/M3A0B9RD5T>. Use remains
subject to ImmPort terms. See [DATA_ACCESS.md](DATA_ACCESS.md).

## Repository structure

```text
R/
  integration/        # metadata and cross-panel score builders
  shared/             # portable paths and common workflow engines
  panels/CP*/         # panel-specific recovered/reconstructed source
config/paths.example.R
scripts/
  01_run_all_panels.R
  02_static_source_audit.R
  10-19 panel runners
  20-28 panel validators
  30_validate_all_panels.R
  99_session_info.R
validation/            # panel-level archive-validation evidence
```

## Setup

Install dependencies:

```bash
Rscript scripts/00_install_dependencies.R
```

Copy `config/paths.example.R` to `config/paths.R` and configure local ImmPort,
FCS, output, metadata, and optional archived-reference paths. The local config
is ignored by Git.

## Metadata reconstruction

`R/integration/STEP0_build_subject_metadata_matrix_RECONSTRUCTED.R` rebuilds the
common subject/FCS metadata matrix from the unpacked ImmPort tabular package.
It links subject, arm, and flow-result tables and standardizes:

- DBG-style subject identifiers and panel/file mappings;
- cancer versus healthy-control group;
- age with the valid modeling range of 18–100 years and age strata;
- sex and binary-sex sensitivity fields.

A local cancer-subgroup/therapy matrix may be supplied with
`SDY2583_CLINICAL_ANNOTATION_FILE`. Participant-level merged outputs remain
local and excluded from Git.

## Run all panels

From the repository root:

```bash
Rscript scripts/01_run_all_panels.R
```

The master workflow builds/loads metadata, runs each panel in dependency order,
performs panel validation, constructs the cross-panel immune-score matrix,
records session information, and applies the all-panel release gate.

Run a subset with:

```bash
SDY2583_PANELS=CP7,CP8,CP25 Rscript scripts/01_run_all_panels.R
```

Individual runners:

```text
CP7   scripts/10_run_cp7_reconstructed.R
CP8   scripts/11_run_cp8_reconstructed.R
CP25  scripts/12_run_cp25_reconstructed.R
CP26  scripts/13_run_cp26_reconstructed.R
CP28  scripts/14_run_cp28_reconstructed.R
CP10  scripts/15_run_cp10_pipeline.R
CP16  scripts/16_run_cp16_recovered.R
CP22  scripts/17_run_cp22_recovered.R
CP23  scripts/18_run_cp23_recovered.R
CP24  scripts/19_run_cp24_full_reconstructed.R
```

The CP24 runner executes the reconstructed 850-subject full-cohort sequence;
the recovered pilot/QC codebook remains available separately.

## Validation

Panel-specific validators check fixed non-sensitive archive benchmarks,
including file counts where available, participant-table dimensions,
feature/score counts, selected adjusted coefficients, matching counts,
robustness classifications, and required outputs.

When `SDY2583_<PANEL>_REFERENCE_DIR` is configured, generated CSVs can also be
compared with archived tables by stable keys and numerical tolerance. Reports
are written under `outputs/<PANEL>/validation/`.

The completed release validation reported here is downstream archive
validation from archived participant-level tables. It does not claim a fresh
raw-FCS-to-results execution.

## Static checks

A data-free audit parses every R file, validates runner references and all ten
panel directories, rejects author-specific paths and interactive file selection,
and screens for accidentally committed data:

```bash
Rscript scripts/02_static_source_audit.R
```

The same audit runs through GitHub Actions in
`.github/workflows/static-r-audit.yml`.

## Safeguards

- Raw data and participant-level outputs are not committed.
- `.RData`, `.rds`, FCS, local configuration, and outputs are ignored.
- Recovered and reconstructed source remain visibly distinguished.
- Session information and validation reports can be retained locally as
  additional execution evidence.
- CP22 immunoglobulin-isotype scores represent B-cell isotype
  architecture/repatterning, not total IgG or total IgA abundance.

## Citation and license

Cite this repository together with the ImmPort dataset and the primary SDY2583
publication. Machine-readable metadata are in [CITATION.cff](CITATION.cff).
Code is released under the [MIT License](LICENSE); the ImmPort data are not
covered by that license.
