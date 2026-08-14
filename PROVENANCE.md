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

Neither the `SAFE` nor `RECONSTRUCTED` suffix alone implies raw-FCS end-to-end
verification. The accepted public release validation is downstream archive
validation from archived post-extraction participant tables. Fresh raw-FCS
re-execution remains outside the accepted release boundary.

## Panel coverage

### CP7
Original source was not visible. Reconstructed Steps 1–7B cover FCS QC, feature
extraction, metadata/age QC, adjusted statistics, composite scores, age and
threshold sensitivity, clinical annotation, figures, and representative gating.

### CP8
Original source was not visible. Reconstructed Steps 1–7B cover the CD4
helper/regulatory workflow, including IL7RA-low/CD25, CCR4/CCR6/CXCR5 axes,
composite scores, age/threshold sensitivity, annotation, and figures.

### CP10
Recovered Steps 2–8 include feature extraction through manuscript export.
Original Step 1/1B source was not located; reconstructed inventory/QC and
mismatch-diagnostic steps provide the interface expected by recovered Step 2.

### CP16
Recovered Steps 1, 1B, 2, 3A, 3B, 4, 4B, 5, 6, 7A, and 7B. An ordered runner
and structural/archive validator are included.

### CP22
Recovered Steps 1, 1B, 2, 3A, 3B, 4, 4B, 5, 6, 7A, and 7B. The isotype
composite represents B-cell immunoglobulin-isotype architecture/repatterning,
not total IgG or total IgA abundance.

### CP23
Recovered Steps 1, 1B, 2, 3A, 3B, 4, 4B, 5, 6, 6B, 7A, and 7B. The final v2
threshold step and explicit Step 6B time-from-start repair were retained.

### CP24
The recovered source is a pilot/QC codebook and remains labeled as such. The
850-subject final R sequence was reconstructed from the archived final output
tree. Pilot and final pipelines are not conflated.

### CP25
Original source was not visible. Reconstructed Steps 1–7B cover the CD4
regulatory-checkpoint workflow, nine composite scores, event-QC sensitivity,
age sensitivity, threshold analysis, annotation, and figures.

### CP26
Original source was not visible. Reconstructed Steps 1–7B cover the NK panel,
including the archived dump-channel fallback, five scores, event-QC, matching,
threshold robustness, annotation, and figures.

### CP28
Original source was not visible. Reconstructed Steps 1–7B cover the T/NK
interface, including T-cell composition, CD8 differentiation, innate-like axes,
nine scores, age/threshold sensitivity, annotation, and figures.

## Integration source

Reconstructed builders link ImmPort subject, arm, flow-result and local
panel-score tables into a participant-level integrated matrix. Participant-level
outputs are written only to ignored local directories and are not redistributed.

## v1.1.0 locked internal validation provenance

`R/integration/locked_internal_validation.R` is a new public analytical layer.
It is **not recovered historical source**. It implements the locked held-out
analysis defined after the original full-cohort score framework had been
completed.

Its provenance boundary is explicit:

- exact score membership comes from the public 66-score registry;
- original SDY2583 TRAINING/VALIDATION labels are used unchanged;
- component direction and mean/SD are learned only in TRAINING;
- higher-level scores are rebuilt recursively from training-frozen components;
- PCA is fitted only in TRAINING and VALIDATION is projected without refitting;
- cross-panel residualization coefficients are estimated only in TRAINING;
- only aggregate validation outputs are public.

Because score membership predates this analysis, v1.1.0 is correctly described
as **locked internal held-out validation**, not external validation and not a
fully de novo discovery/validation experiment.

## Execution and validation architecture

- Every panel has an ordered runner.
- Reconstructed panels have fixed benchmark and optional table-level validators.
- Recovered panels have ordered execution and structural/archive validation.
- `scripts/01_run_all_panels.R` orchestrates panel workflows and integration.
- `scripts/44_run_locked_internal_validation.R` runs the v1.1.0 held-out layer.
- `scripts/02_static_source_audit.R` performs data-free source/privacy auditing.
- GitHub Actions runs the static audit on pushes and pull requests.

## Current reproducibility claim

The public repository is a provenance-tracked, aggregate-only reproducibility
release with all ten panel downstream analyses archive-validated and final
integrated analyses numerically executed on the local ALL10 matrix. Version
1.1.0 additionally publishes the design, aggregate reference results and
portable R implementation for the locked 503/347 internal held-out validation.

Fresh raw-FCS end-to-end re-execution remains outside the accepted release
scope. None of the internal analyses should be described as independent
external validation or functional validation.
