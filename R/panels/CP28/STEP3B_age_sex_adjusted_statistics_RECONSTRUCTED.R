# CP28 Step 3B reconstructed archive-compatible feature models.
# Duplicate endpoints that belong to two modules are modeled once per module row.

rm(list = ls())
source(file.path(
  Sys.getenv("SDY2583_REPO_ROOT", unset = "."),
  "R", "shared", "reconstructed_panel_framework.R"
))
rp_install_and_load(c("dplyr", "readr", "tibble", "purrr"))
source(file.path(
  sd_repo_root(), "R", "panels", "CP28", "MANIFEST_RECONSTRUCTED.R"
))

analysis_dir <- sd_analysis_dir("CP28")
rdata_dir <- file.path(analysis_dir, "11_RData")
out_dir <- file.path(analysis_dir, "06_statistics_age_sex_adjusted")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

load(file.path(
  rdata_dir,
  "SDY2583_CP28_STEP3A_metadata_merge_age_QC_RECONSTRUCTED.RData"
))
module_map <- CP28_MANIFEST$module_map

format_family <- function(result) {
  result |>
    dplyr::group_by(module) |>
    dplyr::mutate(FDR_within_module = p.adjust(p_value, method = "BH")) |>
    dplyr::ungroup() |>
    dplyr::mutate(FDR_all = p.adjust(p_value, method = "BH")) |>
    dplyr::select(
      module, feature, n_model, n_healthy, n_cancer,
      healthy_mean, cancer_mean, beta_cancer_vs_healthy,
      ci_low, ci_high, p_value, direction,
      FDR_within_module, FDR_all
    )
}

run_main_family <- function(binary_sex = FALSE) {
  purrr::pmap_dfr(module_map, function(module, feature) {
    rp_fit_one(
      analysis_data,
      feature,
      subset_label = ifelse(binary_sex, "binary_sex_only", "main"),
      binary_sex = binary_sex,
      event_col = NULL,
      min_events = 0
    ) |>
      dplyr::mutate(module = module, .before = feature)
  }) |>
    format_family()
}

run_event_family <- function() {
  purrr::pmap_dfr(module_map, function(module, feature) {
    event_column <- unname(CP28_MANIFEST$feature_event_qc[[module]])
    if (is.na(event_column) || !nzchar(event_column)) event_column <- NULL
    rp_fit_one(
      analysis_data,
      feature,
      subset_label = "event_QC_ge300",
      binary_sex = FALSE,
      event_col = event_column,
      min_events = ifelse(is.null(event_column), 0, CP28_MANIFEST$event_qc)
    ) |>
      dplyr::mutate(module = module, .before = feature)
  }) |>
    format_family()
}

main_statistics <- run_main_family(FALSE)
binary_sex_statistics <- run_main_family(TRUE)
event_qc_statistics <- run_event_family()

direction_sensitivity <- main_statistics |>
  dplyr::select(
    module, feature,
    beta_main = beta_cancer_vs_healthy,
    direction_main = direction,
    FDR_all_main = FDR_all
  ) |>
  dplyr::left_join(
    binary_sex_statistics |>
      dplyr::select(
        module, feature,
        beta_binary = beta_cancer_vs_healthy,
        direction_binary = direction,
        FDR_all_binary = FDR_all
      ),
    by = c("module", "feature")
  ) |>
  dplyr::left_join(
    event_qc_statistics |>
      dplyr::select(
        module, feature,
        beta_event_QC = beta_cancer_vs_healthy,
        direction_event_QC = direction,
        FDR_all_event_QC = FDR_all
      ),
    by = c("module", "feature")
  ) |>
  dplyr::mutate(
    direction_preserved_binary = direction_main == direction_binary,
    direction_preserved_event_QC = direction_main == direction_event_QC,
    direction_preserved_all = direction_preserved_binary &
      direction_preserved_event_QC
  )

readr::write_csv(
  module_map,
  file.path(out_dir, "SDY2583_CP28_module_feature_map_STEP3B_RECONSTRUCTED.csv")
)
readr::write_csv(
  main_statistics,
  file.path(out_dir, "SDY2583_CP28_age_sex_adjusted_statistics_STEP3B_RECONSTRUCTED.csv")
)
readr::write_csv(
  binary_sex_statistics,
  file.path(out_dir, "SDY2583_CP28_binary_sex_sensitivity_statistics_STEP3B_RECONSTRUCTED.csv")
)
readr::write_csv(
  event_qc_statistics,
  file.path(out_dir, "SDY2583_CP28_event_QC_sensitivity_statistics_STEP3B_RECONSTRUCTED.csv")
)
readr::write_csv(
  direction_sensitivity,
  file.path(out_dir, "SDY2583_CP28_direction_sensitivity_STEP3B_RECONSTRUCTED.csv")
)

save(
  main_statistics,
  binary_sex_statistics,
  event_qc_statistics,
  direction_sensitivity,
  module_map,
  file = file.path(
    rdata_dir,
    "SDY2583_CP28_STEP3B_statistics_RECONSTRUCTED.RData"
  )
)
cat("CP28 modeled module rows:", nrow(main_statistics), "\n")
