# Run the recovered complete CP23 myeloid/macrophage-like pipeline.
source(file.path(Sys.getenv("SDY2583_REPO_ROOT",unset=getwd()),"R","shared","pipeline_runner.R"))
scripts<-c(
"R/panels/CP23/SDY2583_CP23_STEP1_fcs_inventory_marker_QC_SAFE.R",
"R/panels/CP23/STEP1B_mismatch_diagnostic_SAFE.R",
"R/panels/CP23/STEP2_feature_extraction_myeloid_macrophage_like_SAFE.R",
"R/panels/CP23/STEP3A_metadata_merge_age_QC_SAFE.R",
"R/panels/CP23/STEP3B_age_sex_adjusted_statistics_SAFE.R",
"R/panels/CP23/STEP4_composite_scores_SAFE_v2.R",
"R/panels/CP23/STEP4B_age_sensitivity_caliper_interaction_SAFE.R",
"R/panels/CP23/STEP5_threshold_sensitivity_TARGETED_SAFE_v2.R",
"R/panels/CP23/STEP6_clinical_annotation_SAFE.R",
"R/panels/CP23/STEP6B_time_from_start_repair_SAFE.R",
"R/panels/CP23/STEP7A_PATCH_age_threshold_figures_SAFE.R",
"R/panels/CP23/STEP7B_flowjo_style_representative_gating_SAFE.R")
run_sdy2583_script_sequence("CP23",scripts,"scripts/25_validate_recovered_panels.R")
