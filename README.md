# SDY2583 pan-cancer immune remodeling

Reproducible code archive for a secondary, multi-panel reanalysis of the
publicly available ImmPort SDY2583 *Blood Immunotypes* flow-cytometry study.
The analysis quantifies cancer-associated peripheral immune remodeling across
complementary lymphoid and myeloid panels, with demographic adjustment and
prespecified robustness analyses.

## Repository status

This repository is a provenance-checked working release. It contains only
source scripts recovered from the analysis archive. Raw FCS files,
participant-level clinical data, derived participant-level matrices, and
analysis outputs are intentionally excluded.

The archived source-code coverage is not identical across panels:

| Panel | Archived code included | Important limitation |
|---|---|---|
| CP10 | Steps 2–7B | Step 1/1B source scripts were not located in the archive. |
| CP16 | Steps 1–7B | Complete archived panel sequence was located. |
| CP22 | Steps 1B–7B | Step 1 source script was not located in the archive. |
| CP23 | Steps 1–7B, including the time-from-start repair | Complete archived panel sequence was located. |
| CP24 | Pilot/QC codebook | This is not represented as the final full-cohort production pipeline. |
| CP7, CP8, CP25, CP26, CP28 | Not included | Output artifacts were present, but the corresponding `.R` source scripts were not located. |

Panel-specific clinical-annotation scripts are included for CP10, CP16, CP22,
and CP23. A standalone all-panel clinical-integration source script was not
located. See [PROVENANCE.md](PROVENANCE.md) for the detailed code-coverage
statement.

## Data access

The source data are not redistributed here. Obtain SDY2583 directly from
ImmPort:

- Study: **SDY2583 — Blood Immunotypes**
- Dataset DOI: <https://doi.org/10.21430/M3A0B9RD5T>
- ImmPort: <https://www.immport.org/>

Use of the data remains subject to the ImmPort terms and the documentation
distributed with the study. Additional details are provided in
[DATA_ACCESS.md](DATA_ACCESS.md).

## Repository layout

```text
R/
  shared/bootstrap.R       # portable path helpers
  panels/CP10/             # verified archived CP10 scripts
  panels/CP16/             # verified archived CP16 scripts
  panels/CP22/             # verified archived CP22 scripts
  panels/CP23/             # verified archived CP23 scripts
  panels/CP24/             # verified pilot/QC codebook only
config/paths.example.R     # local path configuration template
scripts/00_install_dependencies.R
scripts/99_session_info.R
```

## Software setup

The exact R and package versions used in the original analysis environment
were not preserved in the Drive archive and are therefore not inferred here.
To install the packages used by the recovered scripts, run from the repository
root:

```r
source("scripts/00_install_dependencies.R")
```

After reproducing an analysis, capture the actual environment with:

```r
source("scripts/99_session_info.R")
```

The resulting `session-info.txt` is written under `outputs/` and is ignored by
Git. It should be retained with the analysis record and reported in the final
software-availability metadata.

## Path configuration

Scripts no longer depend on author-specific Desktop paths. Run them from the
repository root. By default, raw data are expected under `data/raw/<PANEL>` and
outputs are written under `outputs/<PANEL>`.

For data stored elsewhere, copy `config/paths.example.R` to
`config/paths.R`, edit the local paths, and source it before running a panel
script:

```r
source("config/paths.R")
source("R/panels/CP16/STEP1_fcs_inventory_marker_QC_SAFE.R")
```

`config/paths.R` is ignored by Git. The same configuration can be supplied
directly through the environment variables listed in the example file.

## Execution order

Run scripts within each panel in filename step order. CP16 and CP23 contain the
complete archived sequences and can begin with Step 1 after the corresponding
FCS directory is configured. CP10 begins at Step 2 and CP22 at Step 1B; these
scripts require the expected Step 1 `.RData` object to have been produced by
the original, currently unlocated source step.

Clinical-annotation and figure scripts also require the derived inputs named
inside each script. These inputs are deliberately not versioned because they
contain participant-level data or are reproducible outputs from earlier steps.

## Reproducibility and privacy safeguards

- No raw FCS files or participant-level tables are committed.
- `.RData`, `.rds`, FCS, local path configuration, and output directories are
  excluded by `.gitignore`.
- Recovered scripts were modified only to replace machine-specific paths and
  interactive file selection with portable configuration helpers. Gating
  thresholds, feature definitions, score directions, statistical models, and
  figure logic were not intentionally altered.
- Missing archived source steps are documented rather than reconstructed from
  outputs.

## Citation

If this code is used, cite the software repository together with the ImmPort
dataset and the primary SDY2583 publication. Machine-readable software
citation metadata are provided in [CITATION.cff](CITATION.cff).

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
