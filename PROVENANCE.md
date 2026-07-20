# Source-code provenance and coverage

## Provenance classes

The repository deliberately separates two classes of analysis source.

### Recovered source

Recovered files were found as `.R`, `.Rmd`, or `.qmd` source in the project
Google Drive archive. Large raw-data archives, serialized workspaces,
participant-level matrices, and output folders were excluded.

Recovered scripts were preserved in archived step order. Machine-specific
`~/Desktop` and `~/Downloads` paths and interactive file selection were
replaced with shared helpers in `R/shared/bootstrap.R`. Gating thresholds,
feature definitions, score directions, statistical models, and figure logic
were not intentionally rewritten. Recovered portable files generally contain
`SAFE` in their names.

### Reconstructed source

Reconstructed files were created only where the corresponding original source
was not visible in the inspected archive. Reconstruction is grounded in the
available analysis record, including:

- archived output-table schemas and numerical benchmarks;
- marker/channel maps and threshold tables;
- Methods and Results records;
- composite-score feature definitions;
- recovered conventions from analytically parallel panels.

Reconstructed files contain `RECONSTRUCTED` in their names. They are not
represented as recovered original code and are not promoted to verified status
until their generated outputs pass archived benchmark comparisons.

## Panel coverage

### CP10

Recovered: feature extraction (Step 2), metadata/age QC (Step 3A), adjusted
statistics (Step 3B), composite scores (Step 4), age sensitivity (Step 4B),
threshold sensitivity (Step 5), clinical annotation (Step 6), figure
generation (Steps 7A–7B), and manuscript export (Step 8).

Not located: Step 1 FCS inventory/QC and Step 1B mismatch-diagnostic source.

### CP16

Recovered: Steps 1, 1B, 2, 3A, 3B, 4, 4B, 5, 6, 7A, and 7B.

### CP22

Recovered: Steps 1, 1B, 2, 3A, 3B, 4, 4B, 5, 6, 7A, and 7B. Step 1 was
stored directly in the panel's `FULL_850_ANALYSIS` root rather than its script
subfolder and was recovered during the recursive archive audit.

### CP23

Recovered: Steps 1, 1B, 2, 3A, 3B, 4, 4B, 5, 6, 6B, 7A, and 7B. A superseded
threshold-sensitivity file and a recovery utility were excluded in favor of
the archived final `v2` step and the explicit Step 6B repair.

### CP24

Recovered: a commented pilot/QC codebook. The Drive archive contained a final
analysis package and output folders, but a corresponding final full-cohort R
source sequence was not located. The codebook must not be cited as the complete
CP24 production pipeline.

### CP7

Original `.R`, `.Rmd`, or `.qmd` source was not visible in the recursively
inspected Drive analysis tree. The following provisional reconstructed steps
are now included:

- Step 1: FCS inventory and marker/channel QC;
- Step 2: compensation, fixed-logicle transformation, and CD8/checkpoint
  feature extraction;
- Step 3A: metadata merge and valid-age quality control;
- Step 3B: age- and sex-adjusted feature models, global/module FDR,
  binary-sex sensitivity, and CD3+CD8+ event-count sensitivity;
- Step 4: directed CP7 composite scores and adjusted score models.

The reconstruction is tied to archived CP7 thresholds, the 850-subject feature
schema, the 55-feature module map, composite feature sets, and archived
statistical benchmarks. `scripts/20_validate_cp7_reconstruction.R` performs
fixed benchmark checks and optional full table-level comparison against a local
copy of the archived outputs.

CP7 remains **reconstructed, not yet verified** because R execution and archive
comparison have not been completed in this environment. Remaining CP7 steps to
reconstruct include age stratification/matching/interactions, targeted
threshold sensitivity, clinical annotation, and figure generation.

### CP8, CP25, CP26, and CP28

No original `.R`, `.Rmd`, or `.qmd` source scripts were visible in the
recursively inspected `FULL_850_ANALYSIS` folder trees. Output artifacts and/or
serialized objects were present. Reconstruction will proceed panel by panel
using the same explicit provenance and validation rules applied to CP7.

## Cross-panel clinical integration

Panel-specific clinical-annotation scripts were recovered for CP10, CP16,
CP22, and CP23. The archive included integrated matrices and panel `.RData`
inputs, but a standalone all-panel integration source script was not located.
Participant-level integration files are not included in this repository.

The reconstructed CP7 metadata step can consume a local compatible integrated
metadata matrix. Full reproduction directly from the downloaded ImmPort
tabular package still requires a standalone subject/FCS/clinical integration
builder to be recovered or reconstructed and validated.

## Reproducibility claim

At the current stage, the repository should be described as a
**provenance-tracked working reproducibility release containing recovered and
provisional reconstructed code**. It must not yet be described as a fully
validated end-to-end reproduction of every manuscript panel.

A panel can be promoted to verified reconstructed status only after:

1. all required steps execute from a clean local configuration;
2. file counts, feature schemas, and QC benchmarks match the archive;
3. adjusted coefficients, confidence intervals, p values, and FDR values match
   within the declared numerical tolerance;
4. composite scores and robustness analyses match archived outputs;
5. the R and package environment is captured with `scripts/99_session_info.R`.
