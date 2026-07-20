# CP28 Step 3A reconstructed metadata merge and archive-compatible age QC.

rm(list = ls())
source(file.path(
  Sys.getenv("SDY2583_REPO_ROOT", unset = "."),
  "R", "shared", "reconstructed_panel_framework.R"
))
rp_install_and_load(c("dplyr", "readr", "stringr", "tibble", "purrr"))
source(file.path(
  sd_repo_root(), "R", "panels", "CP28", "MANIFEST_RECONSTRUCTED.R"
))

analysis_dir <- sd_analysis_dir("CP28")
rdata_dir <- file.path(analysis_dir, "11_RData")
out_dir <- file.path(analysis_dir, "04_metadata_merge_age_QC")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

load(file.path(
  rdata_dir,
  "SDY2583_CP28_STEP2_feature_extraction_RECONSTRUCTED.RData"
))
metadata_merge <- rp_merge_metadata(
  "CP28",
  feature_table,
  unique(c(
    CP28_MANIFEST$module_map$feature,
    CP28_MANIFEST$threshold_targets
  ))
)
analysis_data <- metadata_merge$data
metadata <- metadata_merge$metadata

analysis_data <- analysis_data |>
  dplyr::mutate(
    age_group_for_model = factor(
      dplyr::case_when(
        is.na(age_for_model) ~ NA_character_,
        age_for_model < 40 ~ "Young <40",
        age_for_model < 60 ~ "Middle 40-59",
        TRUE ~ "Older 60+"
      ),
      levels = c("Young <40", "Middle 40-59", "Older 60+")
    ),
    model_ready = !is.na(disease_group) & !is.na(age_for_model) & !is.na(sex)
  )

summary <- tibble::tibble(
  n_rows = nrow(analysis_data),
  n_valid_age = sum(!is.na(analysis_data$age_for_model)),
  n_model_ready = sum(analysis_data$model_ready),
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
  file.path(out_dir, "SDY2583_CP28_analysis_data_with_metadata_RECONSTRUCTED.csv")
)
readr::write_csv(
  summary,
  file.path(out_dir, "SDY2583_CP28_metadata_age_QC_summary_RECONSTRUCTED.csv")
)
save(
  analysis_data,
  metadata,
  summary,
  file = file.path(
    rdata_dir,
    "SDY2583_CP28_STEP3A_metadata_merge_age_QC_RECONSTRUCTED.RData"
  )
)
print(summary)
