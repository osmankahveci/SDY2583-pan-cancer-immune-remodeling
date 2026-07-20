# SDY2583 pan-cancer immune remodeling

Reproducible analysis code for a secondary, multi-panel reanalysis of the
publicly available ImmPort SDY2583 *Blood Immunotypes* flow-cytometry study.
The project evaluates cancer-associated peripheral immune remodeling across
complementary lymphoid and myeloid panels using demographic adjustment,
composite phenotypes, and prespecified robustness analyses.

## Repository status

This is a working reproducibility release. Raw FCS files, participant-level
clinical data, derived participant-level matrices, serialized workspaces, and
analysis outputs are intentionally excluded from Git.

Two source-code provenance classes are kept separate:

- **Recovered scripts** are original analysis scripts found in the project
  archive. Machine-specific paths and interactive file selection were replaced
  with portable configuration helpers, but analytical definitions were not
  intentionally changed. These files generally contain `SAFE` in their names.
- **Reconstructed scripts** were rebuilt from archived output schemas,
  thresholds, Methods/Results records, and recovered cross-panel conventions
  when the original source file was not visible in the archive. These files
  contain `RECONSTRUCTED` and remain provisional until their generated outputs
  pass the archived benchmark comparisons.

Do not remove the provenance suffixes. They prevent recovered source from being
confused with code reconstructed after the analysis.

## Current panel coverage

| Panel | Code currently included | Status and limitation |
|---|---|---|
| CP10 | Recovered Steps 2–8 | Step 1/1B source remains unlocated. |
| CP16 | Recovered Steps 1–7B | Complete archived panel sequence located. |
| CP22 | Recovered Steps 1–7B | Complete archived panel sequence located. |
| CP23 | Recovered Steps 1–7B, including Step 6B repair | Complete archived panel sequence located. |
| CP24 | Recovered pilot/QC codebook | Not the final full-cohort production sequence. |
| CP7 | Reconstructed Steps 1, 2, 3A, 3B, and 4 | Requires local execution and archive-output validation before promotion to verified status. |
| CP8, CP25, CP26, CP28 | Not yet reconstructed | Archived outputs exist, but source scripts were not visible in the inspected Drive trees. |

Panel-specific clinical-annotation scripts are included for CP10, CP16, CP22,
and CP23. A standalone all-panel clinical-integration source script has not yet
been recovered or reconstructed. See [PROVENANCE.md](PROVENANCE.md) for the
full coverage statement.

## Data access

The source data are not redistributed here. Obtain SDY2583 directly from
ImmPort:

- Study: **SDY2583 — Blood Immunotypes**
- Dataset DOI: <https://doi.org/10.21430/M3A0B9RD5T>
- ImmPort: <https://www.immport.org/>

Use of the data remains subject to ImmPort terms and the documentation supplied
with the study. See [DATA_ACCESS.md](DATA_ACCESS.md).

## Repository layout

```text
R/
  shared/bootstrap.R
  panels/
    CP7/                    # reconstructed pipeline, currently Steps 1–4
    CP10/                   # recovered pipeline, Steps 2–8
    CP16/                   # recovered complete sequence
    CP22/                   # recovered complete sequence
    CP23/                   # recovered complete sequence
    CP24/                   # recovered pilot/QC codebook
config/
  paths.example.R
scripts/
  00_install_dependencies.R
  10_run_cp7_reconstructed.R
  20_validate_cp7_reconstruction.R
  99_session_info.R
```

## Software setup

The exact R and package versions from the original analysis environment were
not preserved in the Drive archive. Install the dependencies used by the
included scripts from the repository root:

```r
source("scripts/00_install_dependencies.R")
```

After running analyses, capture the actual environment:

```r
source("scripts/99_session_info.R")
```

The generated `session-info.txt` is written under `outputs/`, ignored by Git,
and should be retained with the reproducibility record.

## Path configuration

Scripts do not depend on author-specific Desktop paths. By default, raw data
are expected under `data/raw/<PANEL>` and outputs are written under
`outputs/<PANEL>`.

For data stored elsewhere, copy `config/paths.example.R` to `config/paths.R`,
edit the local paths, and source it before running scripts. The same values can
be supplied through environment variables.

For reconstructed CP7, configure at minimum:

```r
Sys.setenv(
  SDY2583_CP7_FCS_DIR = "/local/path/to/CP7/FCS",
  SDY2583_METADATA_MATRIX_FILE = "/local/path/to/subject_metadata_matrix.csv"
)
```

The metadata matrix must contain a DBG-style `subject_id` and compatible
disease, age, and sex columns. Reconstructing the standalone raw-ImmPort
metadata-integration step remains a separate required task for full end-to-end
reproduction from the downloaded tabular package alone.

## Running reconstructed CP7

From the repository root:

```bash
Rscript scripts/10_run_cp7_reconstructed.R
```

The runner executes:

1. FCS inventory and marker/channel QC;
2. compensation, fixed-logicle transformation, and feature extraction;
3. metadata merge and valid-age QC;
4. age- and sex-adjusted feature models, global/module FDR, binary-sex
   sensitivity, and event-count sensitivity;
5. CP7 composite-score construction and adjusted models;
6. automated benchmark validation.

The runner currently stops after Step 4. CP7 age-stratified/matching,
disease-by-age interaction, threshold-sensitivity, clinical annotation, and
figure-generation scripts still need to be reconstructed and validated.

## CP7 validation

The validation script always checks non-sensitive archived benchmarks,
including:

- 850 CP7 files and 850 unique subjects;
- 850 successful Step 2 extractions;
- median total and CD3+CD8+ event counts;
- 832 valid-age model records;
- 55 main CP7 features;
- selected adjusted coefficients and the integrated composite coefficient.

For full table-level validation, set the local archived-output root:

```r
Sys.setenv(
  SDY2583_CP7_REFERENCE_DIR = "/local/path/to/SDY2583_CP7_FULL_850_ANALYSIS"
)
source("scripts/20_validate_cp7_reconstruction.R")
```

Generated and archived tables are matched by stable keys and compared within a
configurable numeric tolerance. Validation reports are written to
`outputs/CP7/validation/` and are not committed.

## Reproducibility and privacy safeguards

- No raw FCS files or participant-level tables are committed.
- `.RData`, `.rds`, FCS, local path configuration, and output directories are
  excluded by `.gitignore`.
- Recovered scripts and reconstructed scripts are explicitly distinguished.
- Reconstructed code is not described as verified until archived-output
  comparisons pass.
- Missing steps remain documented rather than silently inferred.

## Citation

Cite this software repository together with the ImmPort dataset and the primary
SDY2583 publication. Machine-readable citation metadata are provided in
[CITATION.cff](CITATION.cff).

Key sources:

- Bhattacharya, S., Dunn, P., Thomas, C. G., et al. (2018). ImmPort, toward
  repurposing of open access immunological assay data for translational and
  clinical research. *Scientific Data, 5*, 180015.
  <https://doi.org/10.1038/sdata.2018.15>
- Dyikanov, D., et al. (2024). Comprehensive peripheral blood
  immunoprofiling reveals five immunotypes with immunotherapy response
  characteristics in patients with cancer. *Cancer Cell, 42*(5), 759–779.e12.
  <https://doi.org/10.1016/j.ccell.2024.04.008>

## License

Code is released under the [MIT License](LICENSE). The source dataset is not
covered by this software license.
