# All-panel downstream archive-validation summary

Date: 2026-07-20

## Status

All ten manuscript panels have completed participant-table/downstream archive validation:

| Panel | Source provenance | Downstream archive status | Raw-FCS status |
|---|---|---|---|
| CP7 | Reconstructed | Pass after score, matching, model, and threshold-family corrections | Pending local execution |
| CP8 | Reconstructed | Pass after matching, age-family, FDR, and threshold-family corrections | Pending local execution |
| CP10 | Reconstructed Step 1/1B + recovered Steps 2 onward | Pass; recovered downstream algorithms retained | Pending local execution |
| CP16 | Recovered | Pass; recovered downstream algorithms retained | Pending local execution |
| CP22 | Recovered | Pass; recovered downstream algorithms retained | Pending local execution |
| CP23 | Recovered | Pass; recovered downstream algorithms retained | Pending local execution |
| CP24 | Reconstructed | Pass after substantive outcome-family, PD-1, model-count, and hierarchical-score corrections | Pending local execution |
| CP25 | Reconstructed | Pass after age-family, matching, interaction, and exact threshold-family corrections | Pending local execution |
| CP26 | Reconstructed | Pass after score-universe and integrated-score corrections | Pending local execution |
| CP28 | Reconstructed | Pass after duplicate-module, event-QC, age-family, matching, and threshold-family corrections | Pending local execution |

## Meaning of “downstream archive-validated”

For each panel, archived participant-level feature/metadata/score tables and archived model outputs were used to independently reproduce the applicable downstream layers: adjusted models, sensitivity models, composite scores, matching, age interactions, and threshold robustness. Numeric agreement is at floating-point/machine-precision scale, and exact subject-pair matching was required where archived pair tables were available.

This status does **not** mean that the entire repository has been validated end to end from raw FCS files. The remaining validation boundary is the local raw-data execution layer:

1. inventory all local panel FCS files;
2. verify marker/channel annotations and compensation matrices;
3. apply the fixed-logicle transformation;
4. reproduce threshold-defined gates and participant-level feature tables;
5. compare regenerated Step 1–2 outputs with archived reference trees;
6. run each corrected panel pipeline with its reference-directory environment variable configured;
7. retain all panel validation reports and R session information;
8. run the master all-panel integration and final release gate.

## Release wording

Accurate current description:

> Provenance-tracked all-panel reproducibility release with all ten panel downstream analyses archive-validated; raw-FCS end-to-end validation remains pending.

The repository must not yet be described as a fully end-to-end validated release. The pull request therefore remains a draft.

## Panel-level records

Detailed scope, numerical comparisons, archive-specific rules, corrections, and remaining checks are recorded in the panel-specific files under `validation/`.
