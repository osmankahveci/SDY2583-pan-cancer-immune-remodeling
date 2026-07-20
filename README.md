# SDY2583 pan-cancer immune remodeling

Reproducible analysis code for the secondary multi-panel reanalysis of ImmPort
SDY2583 *Blood Immunotypes*. The repository covers CP7, CP8, CP10, CP16, CP22,
CP23, CP24, CP25, CP26, and CP28, including demographic adjustment, composite
phenotypes, robustness analyses, clinical annotation, and figure generation.

## Release status

This is a provenance-tracked working reproducibility release. Raw FCS files,
participant-level metadata and matrices, serialized workspaces, and generated
outputs are excluded from Git.

Two source classes are kept separate:

- **`SAFE`**: original analysis source recovered from the project archive and
  modified only for portable paths and non-interactive execution.
- **`RECONSTRUCTED`**: source rebuilt where the original script was not visible,
  using archived thresholds, marker maps, feature schemas, composite
  definitions, Methods/Results records, and numerical output tables.

Reconstructed code remains provisional until it executes locally and passes
archived-output validation. Do not remove the provenance suffixes.

## Panel coverage

| Panel | Included source | Status |
|---|---|---|
| CP7 | Reconstructed Steps 1–7B | Full provisional pipeline, runner, validator. |
| CP8 | Reconstructed Steps 1–7B | Full provisional pipeline, runner, validator. |
| CP10 | Reconstructed Steps 1/1B + recovered Steps 2–8 | End-to-end sequence assembled. |
| CP16 | Recovered Steps 1–7B | Complete recovered sequence and runner. |
| CP22 | Recovered Steps 1–7B | Complete recovered sequence and runner. |
| CP23 | Recovered Steps 1–7B, including Step 6B | Complete recovered sequence and runner. |
| CP24 | Recovered pilot codebook + reconstructed full-cohort Steps 1–7B | Pilot and 850-subject analysis kept separate. |
| CP25 | Reconstructed Steps 1–7B | Full provisional pipeline, runner, validator. |
| CP26 | Reconstructed Steps 1–7B | Includes BV510-A dump fallback for annotation-variant files. |
| CP28 | Reconstructed Steps 1–7B | Full provisional T/NK-interface pipeline. |

See [PROVENANCE.md](PROVENANCE.md) for panel-level evidence and limitations.

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
  20-26 panel validators
  30_validate_all_panels.R
  99_session_info.R
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

Reconstructed-panel validators check fixed non-sensitive archive benchmarks,
including file counts, successful extractions, event-count summaries, 832
valid-age records where applicable, feature/score counts, selected adjusted
coefficients, matching counts, robustness classifications, and required
figures.

When `SDY2583_<PANEL>_REFERENCE_DIR` is configured, generated CSVs are also
compared with archived tables by stable keys and numerical tolerance. Reports
are written under `outputs/<PANEL>/validation/`.

Recovered-panel runners use structural execution validation and optional local
archive inventory checks. A panel should be called verified only after all
required scripts execute and its validation report passes.

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
- Session information and validation reports are retained locally as release
  evidence.
- CP22 immunoglobulin-isotype scores represent B-cell isotype
  architecture/repatterning, not total IgG or total IgA abundance.

## Citation and license

Cite this repository together with the ImmPort dataset and the primary SDY2583
publication. Machine-readable metadata are in [CITATION.cff](CITATION.cff).
Code is released under the [MIT License](LICENSE); the ImmPort data are not
covered by that license.
