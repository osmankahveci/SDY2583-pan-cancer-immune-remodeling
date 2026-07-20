# CP25 Step 5 reconstructed targeted threshold-sensitivity analysis.
# The exact archived multiplicity family contains 64 features and 9 scores.

rm(list = ls())
source(file.path(
  Sys.getenv("SDY2583_REPO_ROOT", unset = "."),
  "R", "shared", "reconstructed_panel_framework.R"
))
rp_install_and_load(
  c("dplyr", "readr", "tibble", "purrr", "tidyr"),
  "flowCore"
)
source(file.path(
  sd_repo_root(), "R", "panels", "CP25", "MANIFEST_RECONSTRUCTED.R"
))

analysis_dir <- sd_analysis_dir("CP25")
rdata_dir <- file.path(analysis_dir, "11_RData")
out_dir <- file.path(analysis_dir, "12_threshold_sensitivity")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

load(file.path(
  rdata_dir,
  "SDY2583_CP25_STEP2_feature_extraction_RECONSTRUCTED.RData"
))
load(file.path(
  rdata_dir,
  "SDY2583_CP25_STEP3A_metadata_merge_age_QC_RECONSTRUCTED.RData"
))
metadata_for_merge <- metadata

main_thresholds <- CP25_MANIFEST$thresholds
permissive_thresholds <- main_thresholds - 0.2
stringent_thresholds <- main_thresholds + 0.2

# Low-gate thresholds move in the opposite direction.
permissive_thresholds[c("DUMP_LOW", "IL7RA")] <-
  main_thresholds[c("DUMP_LOW", "IL7RA")] + 0.2
stringent_thresholds[c("DUMP_LOW", "IL7RA")] <-
  main_thresholds[c("DUMP_LOW", "IL7RA")] - 0.2

threshold_sets <- list(
  permissive = permissive_thresholds,
  main = main_thresholds,
  stringent = stringent_thresholds
)

score_names <- c(
  names(CP25_MANIFEST$score_definitions),
  CP25_MANIFEST$integrated_score
)
feature_targets <- CP25_MANIFEST$threshold_target_features
targets <- c(feature_targets, score_names)
stopifnot(
  length(feature_targets) == 64L,
  length(score_names) == 9L,
  length(targets) == 73L,
  length(unique(targets)) == 73L
)

run_one_set <- function(thresholds, label) {
  feature_data <- purrr::map_dfr(
    fcs_files,
    function(path) extract_one(path, thresholds)
  )
  data_with_metadata <- feature_data |>
    dplyr::left_join(metadata_for_merge, by = "subject_id")
  scored <- rp_build_scores(
    data_with_metadata,
    CP25_MANIFEST$score_definitions,
    CP25_MANIFEST$integrated_score
  )

  statistics <- purrr::map_dfr(targets, function(target) {
    rp_fit_one(
      scored,
      target,
      subset_label = label,
      binary_sex = FALSE,
      event_col = CP25_MANIFEST$primary_event_col,
      min_events = 0
    )
  }) |>
    dplyr::mutate(
      threshold_set = label,
      targeted_FDR = p.adjust(p_value, method = "BH")
    ) |>
    dplyr::select(
      threshold_set,
      feature,
      n_model,
      healthy_mean,
      cancer_mean,
      beta_cancer_vs_healthy,
      ci_low,
      ci_high,
      p_value,
      targeted_FDR,
      direction
    )

  list(data = scored, statistics = statistics)
}

results <- purrr::imap(threshold_sets, run_one_set)
statistics_long <- dplyr::bind_rows(lapply(results, `[[`, "statistics"))
scored_data_long <- dplyr::bind_rows(
  purrr::imap(results, function(result, label) {
    result$data$threshold_set <- label
    result$data
  })
)

summary_wide <- statistics_long |>
  dplyr::select(
    threshold_set,
    feature,
    beta_cancer_vs_healthy,
    targeted_FDR,
    direction
  ) |>
  tidyr::pivot_wider(
    names_from = threshold_set,
    values_from = c(beta_cancer_vs_healthy, targeted_FDR, direction)
  ) |>
  dplyr::mutate(
    direction_preserved = direction_main == direction_permissive &
      direction_main == direction_stringent,
    FDR_lt_0p05_all_sets = targeted_FDR_main < 0.05 &
      targeted_FDR_permissive < 0.05 &
      targeted_FDR_stringent < 0.05,
    robust_interpretation = dplyr::case_when(
      direction_preserved & FDR_lt_0p05_all_sets ~
        "Direction and FDR-significance preserved across threshold sets",
      direction_preserved ~
        "Direction preserved; FDR-significance not preserved in all threshold sets",
      TRUE ~ "Direction not preserved across threshold sets"
    )
  )

threshold_table <- purrr::imap_dfr(
  threshold_sets,
  function(thresholds, label) {
    tibble::tibble(
      threshold_set = label,
      marker = names(thresholds),
      threshold = as.numeric(thresholds)
    )
  }
)

readr::write_csv(
  threshold_table,
  file.path(
    out_dir,
    "SDY2583_CP25_threshold_sensitivity_threshold_sets_RECONSTRUCTED.csv"
  )
)
readr::write_csv(
  scored_data_long,
  file.path(
    out_dir,
    "SDY2583_CP25_threshold_sensitivity_scored_data_LONG_RECONSTRUCTED.csv"
  )
)
readr::write_csv(
  statistics_long,
  file.path(
    out_dir,
    "SDY2583_CP25_threshold_sensitivity_statistics_LONG_RECONSTRUCTED.csv"
  )
)
readr::write_csv(
  summary_wide,
  file.path(
    out_dir,
    "SDY2583_CP25_threshold_sensitivity_summary_WIDE_RECONSTRUCTED.csv"
  )
)

save(
  results,
  scored_data_long,
  statistics_long,
  summary_wide,
  targets,
  threshold_sets,
  file = file.path(
    rdata_dir,
    "SDY2583_CP25_STEP5_threshold_sensitivity_RECONSTRUCTED.RData"
  )
)
print(table(summary_wide$robust_interpretation))
