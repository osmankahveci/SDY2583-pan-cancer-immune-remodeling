# Supplementary flow-cytometry figure workflow

This directory contains the reconstructed R workflow used to generate the ten
representative event-level supplementary flow atlases for the SDY2583
pan-cancer manuscript.

## Scope

| Supplementary figure | Panel | Main visualization scope |
|---|---|---|
| S1 | CP7 | CD8 differentiation and checkpoint co-expression |
| S2 | CP8 | CD4 helper/regulatory-like phenotypes |
| S3 | CP10 | Myeloid/granulocytic phenotypes |
| S4 | CP16 | APC/DC-like phenotypes |
| S5 | CP22 | B-cell and humoral architecture |
| S6 | CP23 | Monocyte/myeloid phenotype maps |
| S7 | CP24 | CD8 differentiation and effector-like phenotypes |
| S8 | CP25 | CD4 regulatory-checkpoint phenotypes |
| S9 | CP26 | NK-like receptor and maturation phenotypes |
| S10 | CP28 | T/NK-interface phenotypes |

The shared workflow produces white-background pseudocolor density plots with
analytical threshold lines, parent-population event counts, and percentages
calculated from all events in the stated parent. Only the displayed point layer
is downsampled and rasterized; text, axes, threshold lines, and labels remain
vector elements in the PDF.

## Public-data boundary

Raw FCS files, participant identifiers, participant-specific filenames,
participant-specific validation values, and generated event-level outputs are
not included in this repository. Obtain the source files directly from ImmPort
under SDY2583.

The workflow receives each representative FCS path through a local environment
variable. Do not replace the portable input variables with author-specific
absolute paths before committing changes.

## Configuration

Copy the example configuration and edit only the local copy:

```bash
cp config/supplementary_flow_paths.example.env \
   config/supplementary_flow_paths.local.env
```

The local file is excluded through the repository's `*.local.env` rule. Source
the configuration before running the figures:

```bash
set -a
source config/supplementary_flow_paths.local.env
set +a
```

Set `SDY2583_FLOW_OUTPUT_DIR` to a local output directory. When it is not set,
the workflow writes to `outputs/supplementary_flow/`, which is excluded from
Git.

## Run

Run one panel:

```bash
SDY2583_FLOW_PANEL=CP7 \
Rscript R/figures/supplementary_flow/RECONSTRUCTED/run_supplementary_flow_atlas.R
```

Run all configured panels:

```bash
Rscript scripts/50_run_supplementary_flow_figures.R
```

Run selected panels:

```bash
SDY2583_FLOW_PANELS=CP7,CP22,CP24 \
Rscript scripts/50_run_supplementary_flow_figures.R
```

A selected panel stops with an informative error when its required local FCS
environment variable is unset, its file is missing, its marker map cannot be
resolved, its spillover matrix cannot be applied where required, or a required
parent population has insufficient events.

## Preprocessing and interpretation

Most panels are read without automatic transformation, compensated using the
file-level spillover matrix, and transformed using a fixed logicle transform
with `W = 0.5`, `T = 262144`, `M = 4.5`, and `A = 0`. CP22 follows its completed
analysis workflow and does not apply an additional compensation step when the
embedded spillover information is not directly usable as a matrix.

The plotted gates are the fixed analytical thresholds used by the automated
feature-extraction workflows. They document the analytical architecture and are
not exact reproductions of manually drawn FlowJo polygons.

Panel-specific interpretation limits remain explicit in the code and captions.
Examples include the absence of definitive DC identity in CP16, the
phenotype-like rather than functional interpretation of CP22 and CP23
populations, the absence of FOXP3 in CP25, and the absence of CD16 for canonical
NK-subset classification in CP28.

## Provenance and validation status

The workflow is marked `RECONSTRUCTED`. It was rebuilt from archived channel
maps, threshold registries, completed Methods/Results records, and local
representative-figure workflows. The public source has been screened to exclude
raw data, participant identifiers, participant-specific metrics, and
author-specific paths.

Repository CI can perform data-free source checks, but it cannot execute this
workflow without local raw FCS files. A submission-ready figure should be
accepted only after local execution, channel-audit review, parent-event review,
percentage-table review, and visual inspection of the exported PDF.
