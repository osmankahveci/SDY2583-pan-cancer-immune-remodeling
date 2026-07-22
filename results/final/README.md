# Curated final aggregate results

These files are aggregate-only outputs from the final integrated analysis
release. They contain no participant identifiers, participant-level scores,
matched-pair IDs, or individual PCA coordinates.

## Cross-panel, immunotype, PCA, and clustering evidence

- `manuscript_locked_key_results.csv`
- `integrated_systems_analysis_summary.json`
- `cross_panel_convergence_summary.csv`
- `immunotype_principal_score_omnibus.csv`
- `immunotype_adjusted_means_full_cohort.csv`
- `pca_principal_score_loadings.csv`
- `pca_immunotype_centroids.csv`
- `pca_training_validation_replication.csv`
- `clustering_stability_metrics.csv`

The final code additionally regenerates all-subset adjusted marginal means and
full intermediate correlation tables locally. Participant-level residuals and
PCA coordinates are intentionally excluded.

The locked scientific results and interpretation boundaries are also recorded
in:

- `analysis_notes/CROSS_PANEL_CONVERGENCE_MANUSCRIPT_HOLD.md`
- `analysis_notes/IMMUNOTYPE_REMODELING_ASSOCIATION_HOLD.md`
- `analysis_notes/PCA_GLOBAL_REMODELING_HOLD.md`

## Bounded-outcome, matching, and cohort summaries

- `selected_frequency_model_sensitivity.csv`
- `age_matching_diagnostics.csv`
- `age_matched_paired_outcome_results.csv`
- `table1_cohort_characteristics.csv`
- `panel_fcs_availability.csv`
- `consolidated_threshold_robustness.csv`
- `frequency_matching_table1_summary.json`

## Bootstrap stability

- `bootstrap_global_summary.csv`
- `bootstrap_pca_loading_stability.csv`
- `bootstrap_principal_score_effect_stability.csv`
- `bootstrap_all66_component_effect_stability.csv`
- `bootstrap_panel_component_summary.csv`
- `bootstrap_component_stability_summary.json`

Values are tied to the documented final scripts, input schemas, software
implementation, random seed, and accepted validation boundary. Internal
stability is not external validation.
