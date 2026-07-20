# CP28 Step 5 reconstructed exact 80-outcome threshold sensitivity.

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
  sd_repo_root(), "R", "panels", "CP28", "MANIFEST_RECONSTRUCTED.R"
))

analysis_dir <- sd_analysis_dir("CP28")
rdata_dir <- file.path(analysis_dir, "11_RData")
out_dir <- file.path(analysis_dir, "12_threshold_sensitivity")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

load(file.path(
  rdata_dir,
  "SDY2583_CP28_STEP2_feature_extraction_RECONSTRUCTED.RData"
))
load(file.path(
  rdata_dir,
  "SDY2583_CP28_STEP3A_metadata_merge_age_QC_RECONSTRUCTED.RData"
))
metadata_for_merge <- metadata

main_thresholds <- CP28_MANIFEST$thresholds
permissive_thresholds <- main_thresholds - 0.2
stringent_thresholds <- main_thresholds + 0.2
permissive_thresholds["DUMP_LOW"] <- main_thresholds["DUMP_LOW"] + 0.2
stringent_thresholds["DUMP_LOW"] <- main_thresholds["DUMP_LOW"] - 0.2
threshold_sets <- list(
  main = main_thresholds,
  permissive = permissive_thresholds,
  stringent = stringent_thresholds
)

targets <- CP28_MANIFEST$threshold_targets
stopifnot(length(targets) == 80L, length(unique(targets)) == 80L)

run_one_set <- function(thresholds, threshold_set) {
  feature_data <- purrr::map_dfr(
    fcs_files,
    function(path) extract_one(path, thresholds)
  )
  data_with_metadata <- feature_data |>
    dplyr::left_join(metadata_for_merge, by = "subject_id")
  scored <- rp_build_scores(
    data_with_metadata,
    CP28_MANIFEST$score_definitions,
    CP28_MANIFEST$integrated_score
  )
  scored$threshold_set <- threshold_set

  statistics <- purrr::map_dfr(targets, function(feature_name) {
    rp_fit_one(
      scored,
      feature_name,
      subset_label = threshold_set,
      binary_sex = FALSE,
      event_col = NULL,
      min_events = 0
    )
  }) |>
    dplyr::mutate(
      threshold_set = threshold_set,
      targeted_FDR = p.adjust(p_value, method = "BH")
    ) |>
    dplyr::select(
      threshold_set, feature, n_model,
      healthy_mean, cancer_mean,
      beta_cancer_vs_healthy, ci_low, ci_high,
      p_value, targeted_FDR, direction
    )

  list(data = scored, statistics = statistics)
}

results <- purrr::imap(threshold_sets, run_one_set)
scored_data_long <- dplyr::bind_rows(lapply(results, `[[`, "data"))
statistics_long <- dplyr::bind_rows(lapply(results, `[[`, "statistics"))

summary_wide <- statistics_long |>
  dplyr::select(
    threshold_set, feature,
    beta_cancer_vs_healthy, targeted_FDR, direction
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

selected_features_long <- scored_data_long |>
  dplyr::select(
    subject_id, threshold_set,
    dplyr::any_of(c(
      "disease_group", "age_for_model", "sex", "sex_binary",
      "model_ready", "n_dump_low", "n_cd3_pos_dump_low",
      "n_cd3_cd8_pos", "n_cd3neg_cd56pos_nk_like"
    )),
    dplyr::all_of(targets)
  )

threshold_table <- purrr::imap_dfr(
  threshold_sets,
  function(thresholds, threshold_set) {
    tibble::tibble(
      threshold_set = threshold_set,
      marker = names(thresholds),
      threshold = as.numeric(thresholds)
    )
  }
)

readr::write_csv(
  threshold_table,
  file.path(out_dir, "SDY2583_CP28_threshold_sensitivity_threshold_sets_RECONSTRUCTED.csv")
)
readr::write_csv(
  selected_features_long,
  file.path(out_dir, "SDY2583_CP28_threshold_sensitivity_selected_features_LONG_RECONSTRUCTED.csv")
)
readr::write_csv(
  scored_data_long,
  file.path(out_dir, "SDY2583_CP28_threshold_sensitivity_scored_data_LONG_RECONSTRUCTED.csv")
)
readr::write_csv(
  statistics_long,
  file.path(out_dir, "SDY2583_CP28_threshold_sensitivity_statistics_LONG_RECONSTRUCTED.csv")
)
readr::write_csv(
  summary_wide,
  file.path(out_dir, "SDY2583_CP28_threshold_sensitivity_summary_WIDE_RECONSTRUCTED.csv")
)

save(
  results,
  scored_data_long,
  selected_features_long,
  statistics_long,
  summary_wide,
  targets,
  threshold_sets,
  file = file.path(
    rdata_dir,
    "SDY2583_CP28_STEP5_threshold_sensitivity_RECONSTRUCTED.RData"
  )
)
print(table(summary_wide$robust_interpretation))
