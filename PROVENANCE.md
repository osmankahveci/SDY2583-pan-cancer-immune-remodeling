# Source-code provenance and coverage

## Provenance classes

The repository separates recovered source from reconstructed source.

### Recovered source (`SAFE`)

Recovered scripts were located as R source in the project Drive archive. Their
archived step order and analytical definitions were retained. Author-specific
paths and interactive file selection were replaced with portable helpers; raw
and participant-level data were excluded.

### Reconstructed source (`RECONSTRUCTED`)

Reconstructed scripts were created only where original source was not visible.
They are grounded in archived marker/channel maps, thresholds, feature schemas,
Methods/Results records, composite definitions, figures, and numerical output
tables. They are never represented as recovered originals.

### Verified status

Neither the `SAFE` nor `RECONSTRUCTED` suffix alone means that a clean rerun has
passed. Verified status requires successful local execution, archived-output
comparison, and captured session information.

## Panel coverage

### CP7

Original source was not visible. Reconstructed Steps 1–7B now cover FCS QC,
feature extraction, metadata/age QC, adjusted statistics, composite scores, age
sensitivity, threshold sensitivity, clinical annotation, quantitative figures,
and representative real-FCS gating. A runner and benchmark/table validator are
included.

### CP8

Original source was not visible. Reconstructed Steps 1–7B reproduce the CD4
helper/regulatory workflow, including the archived 46-feature schema,
IL7RA-low/CD25 and CCR4/CCR6/CXCR5 axes, five component scores, integrated
score, age sensitivity, threshold sensitivity, annotation, and figures.

### CP10

Recovered Steps 2–8 include feature extraction through manuscript export.
Original Step 1/1B source was not located; reconstructed inventory/QC and
mismatch-diagnostic steps now provide the RData interface expected by recovered
Step 2. The combined runner therefore begins from CP10 FCS files while keeping
the provenance boundary explicit.

### CP16

Recovered Steps 1, 1B, 2, 3A, 3B, 4, 4B, 5, 6, 7A, and 7B. An ordered runner
and structural/archive validator are included.

### CP22

Recovered Steps 1, 1B, 2, 3A, 3B, 4, 4B, 5, 6, 7A, and 7B. Step 1 was found in
the panel root during the recursive audit. The isotype composite must be
interpreted as B-cell immunoglobulin-isotype architecture/repatterning, not
total IgG or total IgA abundance.

### CP23

Recovered Steps 1, 1B, 2, 3A, 3B, 4, 4B, 5, 6, 6B, 7A, and 7B. The final v2
threshold step and explicit Step 6B time-from-start repair were retained.

### CP24

The recovered source is a pilot/QC codebook and remains labeled as such. A
separate 850-subject final analysis output tree was present, but its final R
sequence was not. Reconstructed full-cohort Steps 1–7B now cover the archived
CD8 differentiation workflow, clean-age models, official event-count
sensitivity, composite scores, additional explicitly supplemental threshold
sensitivity, annotation, and figures. The pilot and final pipelines are never
conflated.

### CP25

Original source was not visible. Reconstructed Steps 1–7B cover the CD4
regulatory-checkpoint workflow, dump-low primary gate, 79-feature Step 2 table,
main modules, nine composite scores, event-QC sensitivity at 300 primary-gate
events, age sensitivity, the archived 73-outcome threshold family, annotation,
and figures.

### CP26

Original source was not visible. Reconstructed Steps 1–7B cover the NK panel.
The archived BV510-A fallback for the dump channel is preserved because seven
files lacked the expected marker annotation. The main 15-feature model family,
five scores, NK-like event-QC threshold, matching, threshold robustness,
annotation, and figures are included.

### CP28

Original source was not visible. Reconstructed Steps 1–7B cover the T/NK
interface, including T-cell composition, CD8 differentiation, innate-like
T-cell axes, CD3-negative CD56-positive NK-like phenotypes, nine scores, age
and threshold sensitivity, annotation, and figures.

## Integration source

A reconstructed metadata builder now links the unpacked ImmPort subject, arm,
and flow-result tables into a local subject/FCS/panel matrix. A second builder
integrates locally generated panel scores into a cross-panel matrix. Both write
participant-level outputs only to ignored local directories.

Clinical subgroup and therapy fields can be supplied through a separate local
annotation matrix. The repository does not redistribute these participant-level
outputs.

## Execution and validation architecture

- Every panel has an ordered runner.
- Reconstructed panels have fixed benchmark and optional table-level validators.
- Recovered panels have ordered execution and structural/archive validation.
- `scripts/01_run_all_panels.R` orchestrates metadata, panel pipelines,
  cross-panel integration, session capture, and the all-panel release gate.
- `scripts/02_static_source_audit.R` parses all R files and audits paths,
  references, panel coverage, and accidental data inclusion.
- GitHub Actions runs the data-free static audit on pushes and pull requests.

## Current reproducibility claim

The repository should currently be described as a **provenance-tracked working
reproducibility release containing recovered and provisional reconstructed
code for all manuscript panels**. It must not yet be described as a fully
validated end-to-end release because the newly reconstructed sequences have not
been executed against the local FCS/archive environment in this session.

Promotion to verified status requires:

1. clean execution of all required steps;
2. matching file counts, feature schemas, and QC benchmarks;
3. matching coefficients, confidence intervals, p values, and FDR values within
   declared tolerance;
4. matching composite scores and robustness analyses;
5. successful figure generation;
6. an all-panel passing validation summary;
7. captured R and package session information.
