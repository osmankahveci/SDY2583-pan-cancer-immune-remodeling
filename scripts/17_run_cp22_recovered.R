# Run the recovered complete CP22 B-cell/isotype pipeline.
source(file.path(Sys.getenv("SDY2583_REPO_ROOT",unset=getwd()),"R","shared","pipeline_runner.R"))
scripts<-c(
"R/panels/CP22/STEP1_FCS_inventory_marker_QC_SAFE_v2.R",
"R/panels/CP22/STEP1B_mismatch_diagnostic.R",
"R/panels/CP22/STEP2_feature_extraction_Bcell_humoral_SAFE.R",
"R/panels/CP22/STEP3A_metadata_merge_age_QC_SAFE_v2.R",
"R/panels/CP22/STEP3B_age_sex_adjusted_statistics_SAFE.R",
"R/panels/CP22/STEP4_composite_scores_SAFE.R",
"R/panels/CP22/STEP4B_age_sensitivity_caliper_interaction_SAFE.R",
"R/panels/CP22/STEP5_threshold_sensitivity_TARGETED_SAFE.R",
"R/panels/CP22/STEP6_clinical_annotation_SAFE_v3.R",
"R/panels/CP22/STEP7A_quantitative_figures_SAFE_v2.R",
"R/panels/CP22/STEP7B_FlowJo_style_Bcell_figures_SAFE_v2.R")
run_sdy2583_script_sequence("CP22",scripts,"scripts/25_validate_recovered_panels.R")
