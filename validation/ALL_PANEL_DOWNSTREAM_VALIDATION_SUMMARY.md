# All-panel downstream archive-validation summary

Date: 2026-07-20

## Status

All ten manuscript panels have completed participant-table/downstream archive validation:

| Panel | Source provenance | Downstream archive status | Raw-FCS status |
|---|---|---|---|
| CP7 | Reconstructed | Pass after score, matching, model, and threshold-family corrections | Not performed; accepted as outside the release-validation scope |
| CP8 | Reconstructed | Pass after matching, age-family, FDR, and threshold-family corrections | Not performed; accepted as outside the release-validation scope |
| CP10 | Reconstructed Step 1/1B + recovered Steps 2 onward | Pass; recovered downstream algorithms retained | Not performed; accepted as outside the release-validation scope |
| CP16 | Recovered | Pass; recovered downstream algorithms retained | Not performed; accepted as outside the release-validation scope |
| CP22 | Recovered | Pass; recovered downstream algorithms retained | Not performed; accepted as outside the release-validation scope |
| CP23 | Recovered | Pass; recovered downstream algorithms retained | Not performed; accepted as outside the release-validation scope |
| CP24 | Reconstructed | Pass after substantive outcome-family, PD-1, model-count, and hierarchical-score corrections | Not performed; accepted as outside the release-validation scope |
| CP25 | Reconstructed | Pass after age-family, matching, interaction, and exact threshold-family corrections | Not performed; accepted as outside the release-validation scope |
| CP26 | Reconstructed | Pass after score-universe and integrated-score corrections | Not performed; accepted as outside the release-validation scope |
| CP28 | Reconstructed | Pass after duplicate-module, event-QC, age-family, matching, and threshold-family corrections | Not performed; accepted as outside the release-validation scope |

## Meaning of “downstream archive-validated”

For each panel, archived participant-level feature/metadata/score tables and archived model outputs were used to independently reproduce the applicable downstream layers: adjusted models, sensitivity models, composite scores, matching, age interactions, and threshold robustness. Numeric agreement is at floating-point/machine-precision scale, and exact subject-pair matching was required where archived pair tables were available.

The release-validation boundary intentionally begins at the archived participant-level tables. A fresh end-to-end rerun from the local raw FCS files was not performed because those files were not available through the connected environment at executable scale. This omitted layer includes FCS inventory, compensation, transformation, gating, and regeneration of the participant-level feature tables.

This omission is explicitly documented rather than represented as completed validation. It does not change the verified downstream result: all ten archived analysis pipelines and their panel-specific numerical outputs were reproduced after the corrections recorded in the panel validation files.

## Accepted release wording

> Provenance-tracked all-panel reproducibility release with all ten panel downstream analyses independently archive-validated. Raw-FCS end-to-end re-execution was not performed and is outside the accepted validation scope of this release.

The pull request may proceed to review on this basis. It must not be described as having completed a fresh raw-FCS-to-results end-to-end execution.

## Panel-level records

Detailed scope, numerical comparisons, archive-specific rules, corrections, and remaining limitations are recorded in the panel-specific files under `validation/`.
