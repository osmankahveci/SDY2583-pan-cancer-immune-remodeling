# Raw-FCS validation scope decision

Date: 2026-07-20

The accepted release-validation boundary begins at archived participant-level feature, metadata, score, matching, interaction, and threshold-sensitivity tables.

All ten manuscript panels passed independent downstream archive validation with floating-point/machine-precision numerical agreement and exact subject-pair reproduction where archived pair tables were available.

A fresh raw-FCS end-to-end rerun was not performed. The following layers are therefore outside the accepted validation scope of this release:

- local FCS inventory;
- compensation-matrix reapplication;
- fixed-logicle transformation;
- threshold-gate regeneration;
- regeneration of participant-level feature tables from raw events.

This decision does not permit the repository to be described as freshly validated from raw FCS files. The accurate release statement is:

> Provenance-tracked all-panel reproducibility release with all ten panel downstream analyses independently archive-validated. Raw-FCS end-to-end re-execution was not performed and is outside the accepted validation scope of this release.
