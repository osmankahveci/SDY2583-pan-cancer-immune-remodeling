# CP8 Step 3A reconstructed metadata merge and age QC.
rm(list=ls())
source(file.path(Sys.getenv("SDY2583_REPO_ROOT", unset="."), "R", "shared", "reconstructed_panel_framework.R"))
rp_install_and_load(c("dplyr","readr","stringr","tibble","purrr"))
source(file.path(sd_repo_root(), "R", "panels", "CP8", "MANIFEST_RECONSTRUCTED.R"))
analysis_dir <- sd_analysis_dir("CP8"); rdata_dir <- file.path(analysis_dir,"11_RData"); out_dir <- file.path(analysis_dir,"04_metadata_merge"); age_dir <- file.path(analysis_dir,"05_age_QC")
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE); dir.create(age_dir,recursive=TRUE,showWarnings=FALSE)
load(file.path(rdata_dir,"SDY2583_CP8_STEP2_feature_extraction_RECONSTRUCTED.RData"))
merged <- rp_merge_metadata("CP8", feature_table, CP8_MANIFEST$module_map$feature)
analysis_data <- merged$data; metadata <- merged$metadata
summary <- tibble::tibble(n_rows=nrow(analysis_data), n_disease=sum(!is.na(analysis_data$disease_group)), n_valid_age=sum(!is.na(analysis_data$age_for_model)), n_invalid_age=sum(is.na(analysis_data$age_for_model)), n_healthy=sum(analysis_data$disease_group=="Healthy control",na.rm=TRUE), n_cancer=sum(analysis_data$disease_group=="Cancer patient",na.rm=TRUE), metadata_source=merged$source_file)
age_summary <- analysis_data |> dplyr::group_by(disease_group) |> dplyr::summarise(n_total=dplyr::n(), n_valid_age=sum(!is.na(age_for_model)), mean_age=mean(age_for_model,na.rm=TRUE), sd_age=sd(age_for_model,na.rm=TRUE), median_age=median(age_for_model,na.rm=TRUE), .groups="drop")
readr::write_csv(analysis_data,file.path(out_dir,"SDY2583_CP8_FULL_850_analysis_data_with_metadata_RECONSTRUCTED.csv")); readr::write_csv(summary,file.path(out_dir,"SDY2583_CP8_metadata_merge_summary_RECONSTRUCTED.csv")); readr::write_csv(age_summary,file.path(age_dir,"SDY2583_CP8_age_QC_summary_RECONSTRUCTED.csv"))
save(analysis_data,metadata,summary,age_summary,file=file.path(rdata_dir,"SDY2583_CP8_STEP3A_metadata_merge_age_QC_RECONSTRUCTED.RData")); print(summary)
