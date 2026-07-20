# CP25 Step 3A reconstructed metadata merge and age QC.

rm(list = ls())
source(file.path(
  Sys.getenv("SDY2583_REPO_ROOT", unset = "."),
  "R", "shared", "reconstructed_panel_framework.R"
))
rp_install_and_load(c("dplyr", "readr", "stringr", "tibble", "purrr"))
source(file.path(
  sd_repo_root(), "R", "panels", "CP25", "MANIFEST_RECONSTRUCTED.R"
))

analysis_dir <- sd_analysis_dir("CP25")
rdata_dir <- file.path(analysis_dir, "11_RData")
out_dir <- file.path(analysis_dir, "04_metadata_merge")
age_dir <- file.path(analysis_dir, "05_age_QC")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(age_dir, recursive = TRUE, showWarnings = FALSE)

load(file.path(
  rdata_dir,
  "SDY2583_CP25_STEP2_feature_extraction_RECONSTRUCTED.RData"
))

metadata_merge <- rp_merge_metadata(
  "CP25",
  feature_table,
  CP25_MANIFEST$module_map$feature
)
analysis_data <- metadata_merge$data
metadata <- metadata_merge$metadata

# Match the archived human-readable age-group labels.
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
    model_ready = !is.na(disease_group) & !is.na(age_for_model) & !is.na(sex),
    primary_gate_event_qc = !is.na(n_cd3_cd4_primary) &
      n_cd3_cd4_primary >= CP25_MANIFEST$benchmarks$event_qc
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

age_summary <- analysis_data |>
  dplyr::group_by(disease_group) |>
  dplyr::summarise(
    n_total = dplyr::n(),
    n_valid_age = sum(!is.na(age_for_model)),
    mean_age = mean(age_for_model, na.rm = TRUE),
    sd_age = stats::sd(age_for_model, na.rm = TRUE),
    median_age = stats::median(age_for_model, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(
  analysis_data,
  file.path(out_dir, "SDY2583_CP25_analysis_data_with_metadata_RECONSTRUCTED.csv")
)
readr::write_csv(
  summary,
  file.path(out_dir, "SDY2583_CP25_metadata_merge_summary_RECONSTRUCTED.csv")
)
readr::write_csv(
  age_summary,
  file.path(age_dir, "SDY2583_CP25_age_QC_summary_RECONSTRUCTED.csv")
)

save(
  analysis_data,
  metadata,
  summary,
  age_summary,
  file = file.path(
    rdata_dir,
    "SDY2583_CP25_STEP3A_metadata_merge_age_QC_RECONSTRUCTED.RData"
  )
)
print(summary)
