# CP24 reconstructed targeted threshold-sensitivity extension.
#
# This ±0.2 transformed-unit analysis supplements, but does not replace, the
# archived official event-count sensitivity. Scores use the exact official CP24
# hierarchy defined in MANIFEST_RECONSTRUCTED.R.

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
  sd_repo_root(), "R", "panels", "CP24", "MANIFEST_RECONSTRUCTED.R"
))

analysis_dir <- sd_analysis_dir("CP24")
rdata_dir <- file.path(analysis_dir, "11_RData")
out_dir <- file.path(analysis_dir, "12_threshold_sensitivity_RECONSTRUCTED")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

load(file.path(
  rdata_dir,
  "SDY2583_CP24_STEP2_feature_extraction_RECONSTRUCTED.RData"
))
load(file.path(
  rdata_dir,
  "SDY2583_CP24_STEP3A_metadata_merge_age_QC_RECONSTRUCTED.RData"
))
metadata_for_merge <- metadata

main_thresholds <- CP24_MANIFEST$thresholds
threshold_sets <- list(
  permissive = main_thresholds - 0.2,
  main = main_thresholds,
  stringent = main_thresholds + 0.2
)

targets <- unique(c(
  "pct_cd3_pos_total",
  "pct_cd3_cd8_pos_total",
  "pct_cd8_within_cd3",
  "median_CD62L_in_CD3CD8",
  "median_CD27_in_CD3CD8",
  "pct_naive_like",
  "pct_temra_like",
  "pct_temra_cd27neg",
  "pct_cd62lneg_cd27neg",
  "pct_cd57_cx3cr1_cd95_pos",
  "pct_pd1_temra_cd27neg",
  "pct_pd1high_temra_cd27neg",
  "pct_pd1_cd57_cx3cr1_cd95_pos",
  CP24_MANIFEST$score_names
))

run_one_set <- function(thresholds, label) {
  feature_data <- purrr::map_dfr(
    fcs_files,
    function(path) extract_one(path, thresholds)
  )
  data_with_metadata <- feature_data |>
    dplyr::left_join(metadata_for_merge, by = "subject_id")
  scored <- cp24_build_official_scores(data_with_metadata)

  statistics <- purrr::map_dfr(targets, function(target) {
    rp_fit_one(
      scored,
      target,
      subset_label = label,
      binary_sex = FALSE,
      event_col = CP24_MANIFEST$primary_event_col,
      min_events = 0
    )
  }) |>
    dplyr::mutate(
      targeted_fdr = p.adjust(p_value, method = "BH"),
      threshold_set = label
    )

  list(data = scored, statistics = statistics)
}

results <- purrr::imap(threshold_sets, run_one_set)
statistics_long <- dplyr::bind_rows(lapply(results, `[[`, "statistics"))

summary_wide <- statistics_long |>
  dplyr::select(
    threshold_set,
    feature,
    beta_cancer_vs_healthy,
    p_value,
    targeted_fdr,
    direction
  ) |>
  tidyr::pivot_wider(
    names_from = threshold_set,
    values_from = c(
      beta_cancer_vs_healthy,
      p_value,
      targeted_fdr,
      direction
    )
  ) |>
  dplyr::mutate(
    direction_preserved = direction_main == direction_permissive &
      direction_main == direction_stringent,
    fdr_preserved = targeted_fdr_main < 0.05 &
      targeted_fdr_permissive < 0.05 &
      targeted_fdr_stringent < 0.05,
    robust_class = dplyr::case_when(
      direction_preserved & fdr_preserved ~
        "direction_and_FDR_preserved",
      direction_preserved ~
        "direction_preserved_FDR_not_preserved",
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

readr::write_csv(
  threshold_table,
  file.path(out_dir, "SDY2583_CP24_threshold_sets_RECONSTRUCTED.csv")
)
readr::write_csv(
  statistics_long,
  file.path(
    out_dir,
    "SDY2583_CP24_threshold_sensitivity_statistics_LONG_RECONSTRUCTED.csv"
  )
)
readr::write_csv(
  summary_wide,
  file.path(
    out_dir,
    "SDY2583_CP24_threshold_sensitivity_summary_WIDE_RECONSTRUCTED.csv"
  )
)

save(
  results,
  statistics_long,
  summary_wide,
  targets,
  threshold_sets,
  file = file.path(
    rdata_dir,
    "SDY2583_CP24_STEP5_threshold_sensitivity_RECONSTRUCTED.RData"
  )
)

print(table(summary_wide$robust_class))
