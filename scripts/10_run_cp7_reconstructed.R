# Run the complete reconstructed CP7 pipeline in dependency order.
repo_root<-normalizePath(path.expand(Sys.getenv("SDY2583_REPO_ROOT",unset=getwd())),mustWork=FALSE); Sys.setenv(SDY2583_REPO_ROOT=repo_root)
local_config<-file.path(repo_root,"config","paths.R"); if(file.exists(local_config))source(local_config)
steps<-c(
"R/panels/CP7/STEP1_fcs_inventory_marker_QC_RECONSTRUCTED.R",
"R/panels/CP7/STEP2_feature_extraction_CD8_checkpoint_RECONSTRUCTED.R",
"R/panels/CP7/STEP3A_metadata_merge_age_QC_RECONSTRUCTED.R",
"R/panels/CP7/STEP3B_age_sex_adjusted_statistics_RECONSTRUCTED.R",
"R/panels/CP7/STEP4_composite_scores_RECONSTRUCTED.R",
"R/panels/CP7/STEP4B_age_sensitivity_caliper_interaction_RECONSTRUCTED.R",
"R/panels/CP7/STEP5_threshold_sensitivity_TARGETED_RECONSTRUCTED.R",
"R/panels/CP7/STEP6_clinical_annotation_RECONSTRUCTED.R",
"R/panels/CP7/STEP7A_quantitative_figures_RECONSTRUCTED.R",
"R/panels/CP7/STEP7B_flowjo_style_representative_gating_RECONSTRUCTED.R")
for(relative_path in steps){script_path<-file.path(repo_root,relative_path);if(!file.exists(script_path))stop("Pipeline step is missing: ",script_path);cat("\n============================================================\nRunning:",relative_path,"\n============================================================\n");source(script_path,local=new.env(parent=globalenv()))}
source(file.path(repo_root,"scripts","20_validate_cp7_reconstruction.R"),local=new.env(parent=globalenv())); cat("\nCP7 reconstructed pipeline completed.\n")
