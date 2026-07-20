# CP24 full-cohort Step 3A reconstructed metadata merge and clean-age QC.

rm(list = ls())
source(file.path(
  Sys.getenv("SDY2583_REPO_ROOT", unset = "."),
  "R", "shared", "reconstructed_panel_framework.R"
))
rp_install_and_load(c("dplyr", "readr", "stringr", "tibble", "purrr"))
source(file.path(
  sd_repo_root(), "R", "panels", "CP24", "MANIFEST_RECONSTRUCTED.R"
))

analysis_dir <- sd_analysis_dir("CP24")
rdata_dir <- file.path(analysis_dir, "11_RData")
out_dir <- file.path(analysis_dir, "04_metadata_merge_age_QC")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

load(file.path(
  rdata_dir,
  "SDY2583_CP24_STEP2_feature_extraction_RECONSTRUCTED.RData"
))

metadata_merge <- rp_merge_metadata(
  "CP24",
  feature_table,
  unique(CP24_MANIFEST$module_map$feature)
)
analysis_data <- metadata_merge$data
metadata <- metadata_merge$metadata

# Archive-compatible aliases. Internal portable names remain available.
analysis_data$age_years <- analysis_data$age_raw
analysis_data$age_years_clean <- analysis_data$age_for_model
analysis_data$age_group <- as.character(analysis_data$age_group_for_model)
if ("result_file_name" %in% names(analysis_data) &&
    !("file_name" %in% names(analysis_data))) {
  analysis_data$file_name <- analysis_data$result_file_name
}

summary <- tibble::tibble(
  n_rows = nrow(analysis_data),
  n_valid_age = sum(!is.na(analysis_data$age_for_model)),
  n_invalid_age = sum(is.na(analysis_data$age_for_model)),
  n_complete_age_sex = sum(
    !is.na(analysis_data$age_for_model) & !is.na(analysis_data$sex)
  ),
  n_healthy = sum(
    analysis_data$disease_group == "Healthy control",
    na.rm = TRUE
  ),
  n_cancer = sum(
    analysis_data$disease_group == "Cancer patient",
    na.rm = TRUE
  ),
  metadata_source = metadata_merge$source_file
)

readr::write_csv(
  analysis_data,
  file.path(
    out_dir,
    "SDY2583_CP24_analysis_data_with_metadata_CLEANAGE_RECONSTRUCTED.csv"
  )
)
readr::write_csv(
  summary,
  file.path(out_dir, "SDY2583_CP24_metadata_age_QC_summary_RECONSTRUCTED.csv")
)

save(
  analysis_data,
  metadata,
  summary,
  file = file.path(
    rdata_dir,
    "SDY2583_CP24_STEP3A_metadata_merge_age_QC_RECONSTRUCTED.RData"
  )
)
print(summary)
