# CP26 fixed-benchmark and optional archive-table validation.
source(file.path(Sys.getenv("SDY2583_REPO_ROOT", unset = "."), "R", "shared", "reconstructed_panel_framework.R"))
rp_install_and_load(c("dplyr", "readr", "tibble", "purrr"))
source(file.path(sd_repo_root(), "R", "panels", "CP26", "MANIFEST_RECONSTRUCTED.R"))

analysis_dir <- sd_analysis_dir("CP26")
out_dir <- file.path(analysis_dir, "validation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
tolerance <- suppressWarnings(as.numeric(Sys.getenv("SDY2583_VALIDATION_TOLERANCE", unset = "1e-6")))
if (!is.finite(tolerance)) tolerance <- 1e-6
near <- function(x, y) is.finite(x) && abs(x - y) <= tolerance * max(1, abs(y))
read_required <- function(path) {
  if (!file.exists(path)) stop("Missing CP26 output: ", path)
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
}

feature_summary <- read_required(file.path(analysis_dir, "02_feature_extraction", "SDY2583_CP26_feature_extraction_summary_STEP2_NK_panel_RECONSTRUCTED.csv"))
metadata_summary <- read_required(file.path(analysis_dir, "04_metadata_merge_age_QC", "SDY2583_CP26_metadata_age_QC_summary_RECONSTRUCTED.csv"))
main_statistics <- read_required(file.path(analysis_dir, "05_statistics_age_sex_adjusted", "SDY2583_CP26_age_sex_adjusted_statistics_STEP3B_RECONSTRUCTED.csv"))
binary_statistics <- read_required(file.path(analysis_dir, "05_statistics_age_sex_adjusted", "SDY2583_CP26_binary_sex_sensitivity_statistics_STEP3B_RECONSTRUCTED.csv"))
event_statistics <- read_required(file.path(analysis_dir, "05_statistics_age_sex_adjusted", "SDY2583_CP26_event_QC_sensitivity_statistics_STEP3B_RECONSTRUCTED.csv"))
score_data <- read_required(file.path(analysis_dir, "06_composite_scores", "SDY2583_CP26_analysis_data_with_composite_scores_STEP4_RECONSTRUCTED.csv"))
score_map <- read_required(file.path(analysis_dir, "06_composite_scores", "SDY2583_CP26_score_feature_map_STEP4_RECONSTRUCTED.csv"))
score_statistics <- read_required(file.path(analysis_dir, "06_composite_scores", "SDY2583_CP26_composite_score_statistics_STEP4_RECONSTRUCTED.csv"))
matched_5_data <- read_required(file.path(analysis_dir, "07_age_sensitivity_matching_interaction", "SDY2583_CP26_same_sex_age_matched_5y_STEP4B_RECONSTRUCTED.csv"))
matched_10_data <- read_required(file.path(analysis_dir, "07_age_sensitivity_matching_interaction", "SDY2583_CP26_same_sex_age_matched_10y_STEP4B_RECONSTRUCTED.csv"))
matched_5_results <- read_required(file.path(analysis_dir, "07_age_sensitivity_matching_interaction", "SDY2583_CP26_matched_5y_results_STEP4B_RECONSTRUCTED.csv"))
matched_10_results <- read_required(file.path(analysis_dir, "07_age_sensitivity_matching_interaction", "SDY2583_CP26_matched_10y_results_STEP4B_RECONSTRUCTED.csv"))
matching_summary <- read_required(file.path(analysis_dir, "07_age_sensitivity_matching_interaction", "SDY2583_CP26_matching_summary_STEP4B_RECONSTRUCTED.csv"))
interactions <- read_required(file.path(analysis_dir, "07_age_sensitivity_matching_interaction", "SDY2583_CP26_disease_by_age_interaction_results_STEP4B_RECONSTRUCTED.csv"))
threshold_scored <- read_required(file.path(analysis_dir, "12_threshold_sensitivity", "SDY2583_CP26_threshold_analysis_scored_STEP5_RECONSTRUCTED.csv"))
threshold_statistics <- read_required(file.path(analysis_dir, "12_threshold_sensitivity", "SDY2583_CP26_threshold_sensitivity_statistics_STEP5_RECONSTRUCTED.csv"))
threshold_robustness <- read_required(file.path(analysis_dir, "12_threshold_sensitivity", "SDY2583_CP26_threshold_sensitivity_robustness_STEP5_RECONSTRUCTED.csv"))
threshold_counts <- read_required(file.path(analysis_dir, "12_threshold_sensitivity", "SDY2583_CP26_threshold_sensitivity_robustness_counts_STEP5_RECONSTRUCTED.csv"))

count_class <- function(label) {
  value <- threshold_counts$n_variables[threshold_counts$robust_class == label]
  if (length(value)) as.integer(value[1]) else 0L
}
integrated_beta <- score_statistics$beta_cancer_vs_healthy[score_statistics$score == CP26_MANIFEST$integrated_score][1]

checks <- tibble::tribble(
  ~check, ~observed, ~expected, ~pass,
  "FCS files", feature_summary$n_files[1], 850, feature_summary$n_files[1] == 850,
  "successful feature extraction", feature_summary$n_success[1], 850, feature_summary$n_success[1] == 850,
  "remaining feature failures", feature_summary$n_failed[1], 0, feature_summary$n_failed[1] == 0,
  "median total events", feature_summary$median_total_events[1], 299279.5, near(feature_summary$median_total_events[1], 299279.5),
  "median CD45 dump-low events", feature_summary$median_cd45_dump_low_events[1], 3839, near(feature_summary$median_cd45_dump_low_events[1], 3839),
  "median NK-like events", feature_summary$median_nk_like_events[1], 3505, near(feature_summary$median_nk_like_events[1], 3505),
  "valid ages", metadata_summary$n_valid_age[1], 832, metadata_summary$n_valid_age[1] == 832,
  "main feature models", nrow(main_statistics), 15, nrow(main_statistics) == 15,
  "binary feature models", nrow(binary_statistics), 15, nrow(binary_statistics) == 15,
  "event-QC feature models", nrow(event_statistics), 15, nrow(event_statistics) == 15,
  "main score rows", nrow(score_data), 832, nrow(score_data) == 832,
  "score-map rows", nrow(score_map), 20, nrow(score_map) == 20,
  "integrated raw components", sum(score_map$score == CP26_MANIFEST$integrated_score), 10, sum(score_map$score == CP26_MANIFEST$integrated_score) == 10,
  "score models", nrow(score_statistics), 5, nrow(score_statistics) == 5,
  "integrated beta", integrated_beta, 0.364755491769086, near(integrated_beta, 0.364755491769086),
  "5-year matched rows", nrow(matched_5_data), 526, nrow(matched_5_data) == 526,
  "10-year matched rows", nrow(matched_10_data), 544, nrow(matched_10_data) == 544,
  "5-year matched pairs", dplyr::n_distinct(matched_5_data$matched_pair_id), 263, dplyr::n_distinct(matched_5_data$matched_pair_id) == 263,
  "10-year matched pairs", dplyr::n_distinct(matched_10_data$matched_pair_id), 272, dplyr::n_distinct(matched_10_data$matched_pair_id) == 272,
  "5-year matched models", nrow(matched_5_results), 20, nrow(matched_5_results) == 20,
  "10-year matched models", nrow(matched_10_results), 20, nrow(matched_10_results) == 20,
  "matching summary rows", nrow(matching_summary), 2, nrow(matching_summary) == 2,
  "interaction models", nrow(interactions), 20, nrow(interactions) == 20,
  "threshold scored rows", nrow(threshold_scored), 2550, nrow(threshold_scored) == 2550,
  "threshold rows per set", min(table(threshold_scored$threshold_set)), 850, all(table(threshold_scored$threshold_set) == 850),
  "threshold models", nrow(threshold_statistics), 60, nrow(threshold_statistics) == 60,
  "threshold variables", nrow(threshold_robustness), 20, nrow(threshold_robustness) == 20,
  "threshold direction+FDR robust", count_class("direction_and_FDR_preserved_all_sets"), 10, count_class("direction_and_FDR_preserved_all_sets") == 10,
  "threshold direction-only", count_class("direction_preserved_FDR_not_all_sets"), 7, count_class("direction_preserved_FDR_not_all_sets") == 7,
  "threshold direction failures", count_class("direction_not_preserved"), 3, count_class("direction_not_preserved") == 3
)

reference_dir <- path.expand(Sys.getenv("SDY2583_CP26_REFERENCE_DIR", unset = ""))
if (nzchar(reference_dir) && dir.exists(reference_dir)) {
  comparison_plan <- tibble::tribble(
    ~generated, ~reference_name, ~keys,
    file.path(analysis_dir, "02_feature_extraction", "SDY2583_CP26_FULL_850_feature_table_STEP2_NK_panel_RECONSTRUCTED.csv"), "SDY2583_CP26_FULL_850_feature_table_STEP2_NK_panel_REPAIRED.csv", "subject_id",
    file.path(analysis_dir, "04_metadata_merge_age_QC", "SDY2583_CP26_analysis_data_with_metadata_RECONSTRUCTED.csv"), "SDY2583_CP26_analysis_data_STEP3A_metadata_age_QC.csv", "subject_id",
    file.path(analysis_dir, "05_statistics_age_sex_adjusted", "SDY2583_CP26_age_sex_adjusted_statistics_STEP3B_RECONSTRUCTED.csv"), "SDY2583_CP26_age_sex_adjusted_statistics_STEP3B.csv", "module,feature",
    file.path(analysis_dir, "06_composite_scores", "SDY2583_CP26_analysis_data_with_composite_scores_STEP4_RECONSTRUCTED.csv"), "SDY2583_CP26_analysis_data_with_composite_scores_STEP4.csv", "subject_id",
    file.path(analysis_dir, "06_composite_scores", "SDY2583_CP26_score_feature_map_STEP4_RECONSTRUCTED.csv"), "SDY2583_CP26_score_feature_map_STEP4.csv", "score,feature",
    file.path(analysis_dir, "06_composite_scores", "SDY2583_CP26_composite_score_statistics_STEP4_RECONSTRUCTED.csv"), "SDY2583_CP26_composite_score_statistics_STEP4.csv", "score",
    file.path(analysis_dir, "07_age_sensitivity_matching_interaction", "SDY2583_CP26_same_sex_age_matched_5y_STEP4B_RECONSTRUCTED.csv"), "SDY2583_CP26_same_sex_age_matched_5y_STEP4B.csv", "matched_pair_id,subject_id",
    file.path(analysis_dir, "07_age_sensitivity_matching_interaction", "SDY2583_CP26_same_sex_age_matched_10y_STEP4B_RECONSTRUCTED.csv"), "SDY2583_CP26_same_sex_age_matched_10y_STEP4B.csv", "matched_pair_id,subject_id",
    file.path(analysis_dir, "07_age_sensitivity_matching_interaction", "SDY2583_CP26_matched_5y_results_STEP4B_RECONSTRUCTED.csv"), "SDY2583_CP26_matched_5y_results_STEP4B.csv", "variable",
    file.path(analysis_dir, "07_age_sensitivity_matching_interaction", "SDY2583_CP26_matched_10y_results_STEP4B_RECONSTRUCTED.csv"), "SDY2583_CP26_matched_10y_results_STEP4B.csv", "variable",
    file.path(analysis_dir, "07_age_sensitivity_matching_interaction", "SDY2583_CP26_disease_by_age_interaction_results_STEP4B_RECONSTRUCTED.csv"), "SDY2583_CP26_disease_by_age_interaction_results_STEP4B.csv", "variable",
    file.path(analysis_dir, "12_threshold_sensitivity", "SDY2583_CP26_threshold_analysis_scored_STEP5_RECONSTRUCTED.csv"), "SDY2583_CP26_threshold_analysis_scored_STEP5.csv", "threshold_set,subject_id",
    file.path(analysis_dir, "12_threshold_sensitivity", "SDY2583_CP26_threshold_sensitivity_statistics_STEP5_RECONSTRUCTED.csv"), "SDY2583_CP26_threshold_sensitivity_statistics_STEP5.csv", "threshold_set,variable",
    file.path(analysis_dir, "12_threshold_sensitivity", "SDY2583_CP26_threshold_sensitivity_robustness_STEP5_RECONSTRUCTED.csv"), "SDY2583_CP26_threshold_sensitivity_robustness_STEP5.csv", "variable"
  )
  references <- list.files(reference_dir, "\\.csv$", full.names = TRUE, recursive = TRUE)
  archive_checks <- purrr::pmap_dfr(comparison_plan, function(generated, reference_name, keys) {
    hit <- references[basename(references) == reference_name]
    comparison <- if (length(hit)) {
      rp_compare_csv(generated, hit[1], strsplit(keys, ",", fixed = TRUE)[[1]], tolerance)
    } else tibble::tibble(pass = FALSE, detail = "reference missing")
    tibble::tibble(check = paste0("archive table: ", reference_name), observed = comparison$detail[1], expected = "matched", pass = comparison$pass[1])
  })
  checks <- dplyr::bind_rows(checks, archive_checks)
}

report_file <- file.path(out_dir, "SDY2583_CP26_reconstruction_validation_report.csv")
readr::write_csv(checks, report_file)
print(checks, n = nrow(checks))
if (any(!checks$pass | is.na(checks$pass))) stop("CP26 validation failed. Review ", report_file)
message("CP26 validation passed.")
