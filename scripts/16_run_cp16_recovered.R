# Run the recovered complete CP16 pipeline.
source(file.path(Sys.getenv("SDY2583_REPO_ROOT",unset=getwd()),"R","shared","pipeline_runner.R"))
scripts<-c(
"R/panels/CP16/STEP1_fcs_inventory_marker_QC_SAFE.R",
"R/panels/CP16/STEP1B_mismatch_diagnostic_SAFE.R",
"R/panels/CP16/STEP2_feature_extraction_APC_DC_myeloid_SAFE.R",
"R/panels/CP16/STEP3A_metadata_merge_age_QC_SAFE.R",
"R/panels/CP16/STEP3B_age_sex_adjusted_statistics_SAFE.R",
"R/panels/CP16/STEP4_composite_scores_SAFE.R",
"R/panels/CP16/STEP4B_age_sensitivity_caliper_interaction_SAFE.R",
"R/panels/CP16/STEP5_threshold_sensitivity_TARGETED_SAFE.R",
"R/panels/CP16/STEP6_clinical_annotation_SAFE.R",
"R/panels/CP16/STEP7A_quantitative_figures_SAFE.R",
"R/panels/CP16/STEP7B_flowjo_style_representative_gating_SAFE.R")
run_sdy2583_script_sequence("CP16",scripts,"scripts/25_validate_recovered_panels.R")
