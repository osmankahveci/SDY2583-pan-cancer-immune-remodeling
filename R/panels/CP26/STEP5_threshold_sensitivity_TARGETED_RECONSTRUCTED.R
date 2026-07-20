# CP26 Step 5 reconstructed targeted threshold sensitivity (20 outcomes).
# Threshold-set scores are standardized independently within each all-850 set.

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
  sd_repo_root(), "R", "panels", "CP26", "MANIFEST_RECONSTRUCTED.R"
))

analysis_dir <- sd_analysis_dir("CP26")
rdata_dir <- file.path(analysis_dir, "11_RData")
out_dir <- file.path(analysis_dir, "12_threshold_sensitivity")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

load(file.path(
  rdata_dir,
  "SDY2583_CP26_STEP2_feature_extraction_RECONSTRUCTED.RData"
))
load(file.path(
  rdata_dir,
  "SDY2583_CP26_STEP3A_metadata_merge_age_QC_RECONSTRUCTED.RData"
))
metadata_for_merge <- metadata

main_thresholds <- CP26_MANIFEST$thresholds
permissive_thresholds <- main_thresholds - 0.2
stringent_thresholds <- main_thresholds + 0.2
permissive_thresholds["DUMP_LOW"] <- main_thresholds["DUMP_LOW"] + 0.2
stringent_thresholds["DUMP_LOW"] <- main_thresholds["DUMP_LOW"] - 0.2
threshold_sets <- list(
  main = main_thresholds,
  permissive = permissive_thresholds,
  stringent = stringent_thresholds
)

targets <- c(
  CP26_MANIFEST$targeted_outcomes,
  names(CP26_MANIFEST$score_definitions),
  CP26_MANIFEST$integrated_score
)
stopifnot(length(targets) == 20L, length(unique(targets)) == 20L)

run_one_set <- function(thresholds, label) {
  feature_data <- purrr::map_dfr(
    fcs_files,
    function(path) extract_one(path, thresholds)
  )
  data_with_metadata <- feature_data |>
    dplyr::left_join(metadata_for_merge, by = "subject_id")
  scored <- cp26_build_official_scores(data_with_metadata)
  scored$threshold_set <- label

  statistics <- purrr::map_dfr(targets, function(target) {
    rp_fit_one(
      scored,
      target,
      subset_label = label,
      binary_sex = FALSE,
      event_col = CP26_MANIFEST$primary_event_col,
      min_events = 0
    )
  }) |>
    dplyr::mutate(
      threshold_set = label,
      FDR_targeted = p.adjust(p_value, method = "BH")
    ) |>
    dplyr::transmute(
      threshold_set,
      variable = feature,
      n_model,
      healthy_mean,
      cancer_mean,
      beta_cancer_vs_healthy,
      p_value,
      FDR_targeted,
      direction
    )

  list(data = scored, statistics = statistics)
}

results <- purrr::imap(threshold_sets, run_one_set)
threshold_analysis_scored <- dplyr::bind_rows(lapply(results, `[[`, "data"))
statistics_long <- dplyr::bind_rows(lapply(results, `[[`, "statistics"))

robustness <- statistics_long |>
  tidyr::pivot_wider(
    names_from = threshold_set,
    values_from = c(
      n_model,
      healthy_mean,
      cancer_mean,
      beta_cancer_vs_healthy,
      p_value,
      FDR_targeted,
      direction
    ),
    names_sep = "_"
  ) |>
  dplyr::mutate(
    direction_preserved_all_sets = direction_main == direction_permissive &
      direction_main == direction_stringent,
    FDR_lt_0p05_all_sets = FDR_targeted_main < 0.05 &
      FDR_targeted_permissive < 0.05 &
      FDR_targeted_stringent < 0.05,
    robust_class = dplyr::case_when(
      direction_preserved_all_sets & FDR_lt_0p05_all_sets ~
        "direction_and_FDR_preserved_all_sets",
      direction_preserved_all_sets ~
        "direction_preserved_FDR_not_all_sets",
      TRUE ~ "direction_not_preserved"
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
robustness_counts <- robustness |>
  dplyr::count(robust_class, name = "n_variables")

readr::write_csv(
  threshold_table,
  file.path(out_dir, "SDY2583_CP26_threshold_sets_STEP5_RECONSTRUCTED.csv")
)
readr::write_csv(
  threshold_analysis_scored,
  file.path(out_dir, "SDY2583_CP26_threshold_analysis_scored_STEP5_RECONSTRUCTED.csv")
)
readr::write_csv(
  statistics_long,
  file.path(
    out_dir,
    "SDY2583_CP26_threshold_sensitivity_statistics_STEP5_RECONSTRUCTED.csv"
  )
)
readr::write_csv(
  robustness,
  file.path(
    out_dir,
    "SDY2583_CP26_threshold_sensitivity_robustness_STEP5_RECONSTRUCTED.csv"
  )
)
readr::write_csv(
  robustness_counts,
  file.path(
    out_dir,
    "SDY2583_CP26_threshold_sensitivity_robustness_counts_STEP5_RECONSTRUCTED.csv"
  )
)

save(
  results,
  threshold_analysis_scored,
  statistics_long,
  robustness,
  robustness_counts,
  targets,
  threshold_sets,
  file = file.path(
    rdata_dir,
    "SDY2583_CP26_STEP5_threshold_sensitivity_RECONSTRUCTED.RData"
  )
)
print(robustness_counts)
