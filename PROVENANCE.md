# Source-code provenance and coverage

## Archive source

The files in `R/panels/` were recovered from the project Google Drive archive
for the SDY2583 reanalysis. Only files identifiable as R source code were
included. Large raw-data archives, serialized workspaces, participant-level
matrices, and output folders were excluded.

The recovered scripts were preserved in their archived step order. For
portability, machine-specific `~/Desktop`/`~/Downloads` paths and interactive
file selection were replaced with the shared helpers in
`R/shared/bootstrap.R`. The analytical definitions and numerical procedures
were not intentionally rewritten.

## Verified coverage

### CP10

Recovered: feature extraction (Step 2), metadata/age QC (Step 3A), adjusted
statistics (Step 3B), composite scores (Step 4), age sensitivity (Step 4B),
threshold sensitivity (Step 5), clinical annotation (Step 6), and figure
generation (Steps 7A–7B).

Not located: Step 1 FCS inventory/QC and Step 1B mismatch-diagnostic source.

### CP16

Recovered: Steps 1, 1B, 2, 3A, 3B, 4, 4B, 5, 6, 7A, and 7B.

### CP22

Recovered: Steps 1B, 2, 3A, 3B, 4, 4B, 5, 6, 7A, and 7B.

Not located: Step 1 FCS inventory/QC source.

### CP23

Recovered: Steps 1, 1B, 2, 3A, 3B, 4, 4B, 5, 6, 6B, 7A, and 7B. A
superseded threshold-sensitivity file and a recovery utility were excluded in
favor of the archived final `v2` step and the explicit Step 6B repair.

### CP24

Recovered: a commented pilot/QC codebook. The Drive archive contained a final
analysis package and output folders, but a corresponding final full-cohort R
source sequence was not located. The codebook is therefore labeled as such and
must not be cited as the complete CP24 production pipeline.

### CP7, CP8, CP25, CP26, and CP28

No `.R` source scripts were located in the inspected Drive analysis folders.
Output artifacts and/or serialized objects were present, but they are not a
substitute for source code and are not used here to reconstruct an unverified
pipeline.

### Cross-panel clinical integration

Panel-specific clinical-annotation scripts were recovered for CP10, CP16,
CP22, and CP23. The archive included integrated matrices and panel `.RData`
inputs, but a standalone all-panel integration source script was not located.
Participant-level integration files are not included in this repository.

## Implication for reproducibility claims

The repository supports transparent inspection and reuse of the recovered
analytical code. Full end-to-end reproduction of every manuscript panel will
require the unlocated source steps to be recovered from the original local
analysis environment. Until that occurs, documentation and manuscript code-
availability statements should describe this as a provenance-verified partial
code archive rather than a complete end-to-end workflow.
