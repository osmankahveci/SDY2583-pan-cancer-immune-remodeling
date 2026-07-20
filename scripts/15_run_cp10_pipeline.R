# Run CP10 from reconstructed Step 1/1B through recovered Steps 2–8.
source(file.path(Sys.getenv("SDY2583_REPO_ROOT",unset=getwd()),"R","shared","pipeline_runner.R"))
scripts<-c(
"R/panels/CP10/STEP1_fcs_inventory_marker_QC_RECONSTRUCTED.R",
"R/panels/CP10/STEP1B_mismatch_diagnostic_RECONSTRUCTED.R",
"R/panels/CP10/STEP2_feature_extraction_myeloid_granulocyte_SAFE.R",
"R/panels/CP10/STEP3A_metadata_merge_age_QC_SAFE.R",
"R/panels/CP10/STEP3B_age_sex_adjusted_statistics_SAFE.R",
"R/panels/CP10/STEP4_composite_scores_SAFE.R",
"R/panels/CP10/STEP4B_age_sensitivity_caliper_interaction_SAFE.R",
"R/panels/CP10/STEP5_threshold_sensitivity_TARGETED_SAFE.R",
"R/panels/CP10/STEP6_clinical_annotation_SAFE.R",
"R/panels/CP10/STEP7A_main_figures_SAFE_v2.R",
"R/panels/CP10/STEP7B_flowjo_style_representative_gating_REAL_FCS_SAFE_v2.R",
"R/panels/CP10/STEP8_manuscript_export_SAFE.R")
run_sdy2583_script_sequence("CP10",scripts,"scripts/25_validate_recovered_panels.R")
