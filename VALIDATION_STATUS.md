# Validation status

## Accepted release-validation boundary — 2026-07-20

All ten manuscript panels have completed downstream archive validation from archived participant-level feature, metadata, score, matching, interaction, and threshold-sensitivity tables. Numerical agreement is at floating-point/machine-precision scale, with exact subject-pair reproduction where archived matching tables were available.

A fresh raw-FCS end-to-end rerun was not performed. FCS inventory, compensation, fixed-logicle transformation, threshold gating, and regeneration of participant-level feature tables are therefore explicitly outside the accepted validation scope for this release. This omission must remain visible in repository and manuscript-facing reproducibility wording and must not be described as completed raw-data validation.

Accepted description:

> Provenance-tracked all-panel reproducibility release with all ten panel downstream analyses independently archive-validated. Raw-FCS end-to-end re-execution was not performed and is outside the accepted validation scope of this release.

Detailed panel-level records are retained under `validation/`, including the exact archive-specific rules, numerical comparisons, and corrections applied to reconstructed source.
