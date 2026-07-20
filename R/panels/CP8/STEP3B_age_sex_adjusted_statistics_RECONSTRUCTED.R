# CP8 Step 3B reconstructed adjusted statistics and QC sensitivities.
rm(list=ls())
source(file.path(Sys.getenv("SDY2583_REPO_ROOT", unset="."), "R", "shared", "reconstructed_panel_framework.R"))
rp_install_and_load(c("dplyr","readr","tibble","purrr"))
source(file.path(sd_repo_root(),"R","panels","CP8","MANIFEST_RECONSTRUCTED.R"))
analysis_dir<-sd_analysis_dir("CP8"); rdata_dir<-file.path(analysis_dir,"11_RData"); out_dir<-file.path(analysis_dir,"06_statistics"); sens_dir<-file.path(analysis_dir,"07_sensitivity")
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE); dir.create(sens_dir,recursive=TRUE,showWarnings=FALSE)
load(file.path(rdata_dir,"SDY2583_CP8_STEP3A_metadata_merge_age_QC_RECONSTRUCTED.RData"))
module_map<-CP8_MANIFEST$module_map
main_statistics<-rp_fit_feature_set(analysis_data,module_map)
binary_sex_statistics<-rp_fit_feature_set(analysis_data,module_map,"binary_sex_only",TRUE)
event_statistics<-dplyr::bind_rows(lapply(c(0,500,1000,2000),function(n) rp_fit_feature_set(analysis_data,module_map,paste0("min_CD3CD4_events_",n),FALSE,CP8_MANIFEST$primary_event_col,n)))
readr::write_csv(module_map,file.path(out_dir,"SDY2583_CP8_feature_module_map_RECONSTRUCTED.csv")); readr::write_csv(main_statistics,file.path(out_dir,"SDY2583_CP8_FULL850_age_sex_adjusted_statistics_MAIN_RECONSTRUCTED.csv")); readr::write_csv(main_statistics |> dplyr::filter(fdr_all<.05),file.path(out_dir,"SDY2583_CP8_FULL850_significant_global_FDR_0p05_MAIN_RECONSTRUCTED.csv")); readr::write_csv(main_statistics |> dplyr::filter(fdr_within_module<.05),file.path(out_dir,"SDY2583_CP8_FULL850_significant_module_FDR_0p05_MAIN_RECONSTRUCTED.csv")); readr::write_csv(binary_sex_statistics,file.path(sens_dir,"SDY2583_CP8_binary_sex_sensitivity_statistics_RECONSTRUCTED.csv")); readr::write_csv(event_statistics,file.path(sens_dir,"SDY2583_CP8_event_QC_sensitivity_statistics_RECONSTRUCTED.csv"))
save(main_statistics,binary_sex_statistics,event_statistics,module_map,file=file.path(rdata_dir,"SDY2583_CP8_STEP3B_statistics_RECONSTRUCTED.RData")); cat("CP8 features modeled:",nrow(main_statistics),"\n")
