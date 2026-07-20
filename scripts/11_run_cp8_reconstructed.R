# Run the reconstructed CP8 pipeline from the repository root.
root<-normalizePath(path.expand(Sys.getenv("SDY2583_REPO_ROOT",unset=getwd())),mustWork=FALSE); Sys.setenv(SDY2583_REPO_ROOT=root)
if(file.exists(file.path(root,"config","paths.R"))) source(file.path(root,"config","paths.R"))
steps<-c("STEP1_inventory_marker_QC_RECONSTRUCTED.R","STEP2_feature_extraction_CD4_helper_regulatory_RECONSTRUCTED.R","STEP3A_metadata_merge_age_QC_RECONSTRUCTED.R","STEP3B_age_sex_adjusted_statistics_RECONSTRUCTED.R","STEP4_composite_scores_RECONSTRUCTED.R","STEP4B_age_sensitivity_caliper_interaction_RECONSTRUCTED.R","STEP5_threshold_sensitivity_TARGETED_RECONSTRUCTED.R","STEP6_clinical_annotation_RECONSTRUCTED.R","STEP7A_quantitative_figures_RECONSTRUCTED.R","STEP7B_flowjo_style_representative_gating_RECONSTRUCTED.R")
for(s in steps){cat("\n=== CP8",s,"===\n"); source(file.path(root,"R","panels","CP8",s),local=new.env(parent=globalenv()))}
source(file.path(root,"scripts","21_validate_cp8_reconstruction.R"),local=new.env(parent=globalenv()))
