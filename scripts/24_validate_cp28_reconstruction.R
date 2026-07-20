# CP28 fixed-benchmark and optional archive-table validation.
source(file.path(Sys.getenv("SDY2583_REPO_ROOT", unset = "."), "R", "shared", "reconstructed_panel_framework.R"))
rp_install_and_load(c("dplyr", "readr", "tibble", "purrr"))
source(file.path(sd_repo_root(), "R", "panels", "CP28", "MANIFEST_RECONSTRUCTED.R"))

analysis_dir <- sd_analysis_dir("CP28")
out_dir <- file.path(analysis_dir, "validation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
tolerance <- suppressWarnings(as.numeric(Sys.getenv("SDY2583_VALIDATION_TOLERANCE", unset = "1e-6")))
if (!is.finite(tolerance)) tolerance <- 1e-6
near <- function(x, y) is.finite(x) && abs(x - y) <= tolerance * max(1, abs(y))
read_required <- function(path) {
  if (!file.exists(path)) stop("Missing CP28 output: ", path)
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
}

feature_summary <- read_required(file.path(analysis_dir, "02_feature_extraction", "SDY2583_CP28_feature_extraction_summary_STEP2_TNK_interface_RECONSTRUCTED.csv"))
metadata_summary <- read_required(file.path(analysis_dir, "04_metadata_merge_age_QC", "SDY2583_CP28_metadata_age_QC_summary_RECONSTRUCTED.csv"))
main_statistics <- read_required(file.path(analysis_dir, "06_statistics_age_sex_adjusted", "SDY2583_CP28_age_sex_adjusted_statistics_STEP3B_RECONSTRUCTED.csv"))
binary_statistics <- read_required(file.path(analysis_dir, "06_statistics_age_sex_adjusted", "SDY2583_CP28_binary_sex_sensitivity_statistics_STEP3B_RECONSTRUCTED.csv"))
event_statistics <- read_required(file.path(analysis_dir, "06_statistics_age_sex_adjusted", "SDY2583_CP28_event_QC_sensitivity_statistics_STEP3B_RECONSTRUCTED.csv"))
score_data <- read_required(file.path(analysis_dir, "08_composite_scores", "SDY2583_CP28_analysis_data_with_composite_scores_RECONSTRUCTED.csv"))
score_map <- read_required(file.path(analysis_dir, "08_composite_scores", "SDY2583_CP28_composite_score_feature_sets_RECONSTRUCTED.csv"))
score_statistics <- read_required(file.path(analysis_dir, "08_composite_scores", "SDY2583_CP28_composite_score_statistics_RECONSTRUCTED.csv"))
event_score_statistics <- read_required(file.path(analysis_dir, "08_composite_scores", "SDY2583_CP28_composite_score_event_QC_sensitivity_statistics_RECONSTRUCTED.csv"))
age_stratified <- read_required(file.path(analysis_dir, "10_age_sensitivity", "SDY2583_CP28_age_stratified_statistics_RECONSTRUCTED.csv"))
pairs_5 <- read_required(file.path(analysis_dir, "10_age_sensitivity", "SDY2583_CP28_same_sex_age_matched_pairs_5yr_RECONSTRUCTED.csv"))
pairs_10 <- read_required(file.path(analysis_dir, "10_age_sensitivity", "SDY2583_CP28_same_sex_age_matched_pairs_10yr_RECONSTRUCTED.csv"))
matched_statistics <- read_required(file.path(analysis_dir, "10_age_sensitivity", "SDY2583_CP28_age_matched_sensitivity_statistics_RECONSTRUCTED.csv"))
matching_summary <- read_required(file.path(analysis_dir, "10_age_sensitivity", "SDY2583_CP28_age_matched_sensitivity_summary_RECONSTRUCTED.csv"))
interactions <- read_required(file.path(analysis_dir, "10_age_sensitivity", "SDY2583_CP28_disease_by_age_interaction_statistics_RECONSTRUCTED.csv"))
threshold_scored <- read_required(file.path(analysis_dir, "12_threshold_sensitivity", "SDY2583_CP28_threshold_sensitivity_scored_data_LONG_RECONSTRUCTED.csv"))
threshold_statistics <- read_required(file.path(analysis_dir, "12_threshold_sensitivity", "SDY2583_CP28_threshold_sensitivity_statistics_LONG_RECONSTRUCTED.csv"))
threshold_summary <- read_required(file.path(analysis_dir, "12_threshold_sensitivity", "SDY2583_CP28_threshold_sensitivity_summary_WIDE_RECONSTRUCTED.csv"))

integrated_beta <- score_statistics$beta_cancer_vs_healthy[score_statistics$feature == CP28_MANIFEST$integrated_score][1]
count_class <- function(label) sum(threshold_summary$robust_interpretation == label, na.rm = TRUE)
module_sizes <- table(main_statistics$module)
event_n <- event_statistics |>
  dplyr::group_by(module) |>
  dplyr::summarise(n_model = dplyr::first(n_model), .groups = "drop")
get_event_n <- function(module) event_n$n_model[event_n$module == module][1]

checks <- tibble::tribble(
  ~check, ~observed, ~expected, ~pass,
  "FCS files", feature_summary$n_files[1], 850, feature_summary$n_files[1] == 850,
  "successful feature extraction", feature_summary$n_success[1], 850, feature_summary$n_success[1] == 850,
  "median total events", feature_summary$median_total_events[1], 309700.5, near(feature_summary$median_total_events[1], 309700.5),
  "median dump-low events", feature_summary$median_dump_low_events[1], 77505, near(feature_summary$median_dump_low_events[1], 77505),
  "median CD3 dump-low events", feature_summary$median_cd3_dump_low_events[1], 25096, near(feature_summary$median_cd3_dump_low_events[1], 25096),
  "median CD3+CD8+ events", feature_summary$median_cd3_cd8_events[1], 7440, near(feature_summary$median_cd3_cd8_events[1], 7440),
  "median CD3-negative CD56-positive events", feature_summary$median_cd3neg_cd56pos_events[1], 3066.5, near(feature_summary$median_cd3neg_cd56pos_events[1], 3066.5),
  "valid ages", metadata_summary$n_valid_age[1], 832, metadata_summary$n_valid_age[1] == 832,
  "main module rows", nrow(main_statistics), 65, nrow(main_statistics) == 65,
  "main unique features", dplyr::n_distinct(main_statistics$feature), 62, dplyr::n_distinct(main_statistics$feature) == 62,
  "binary module rows", nrow(binary_statistics), 65, nrow(binary_statistics) == 65,
  "event-QC module rows", nrow(event_statistics), 65, nrow(event_statistics) == 65,
  "composition event-QC N", get_event_n("composition_T_NK_interface"), 832, get_event_n("composition_T_NK_interface") == 832,
  "CD8 differentiation event-QC N", get_event_n("cd8_differentiation_CD45RA_CD27"), 823, get_event_n("cd8_differentiation_CD45RA_CD27") == 823,
  "CD8 effector event-QC N", get_event_n("cd8_effector_NK_like_axis"), 823, get_event_n("cd8_effector_NK_like_axis") == 823,
  "CD3 innate event-QC N", get_event_n("cd3_innate_like_T_cell_axis"), 829, get_event_n("cd3_innate_like_T_cell_axis") == 829,
  "NK-like event-QC N", get_event_n("cd3neg_CD56pos_NK_like_axis"), 810, get_event_n("cd3neg_CD56pos_NK_like_axis") == 810,
  "subject-level score rows", nrow(score_data), 850, nrow(score_data) == 850,
  "score-map rows", nrow(score_map), 48, nrow(score_map) == 48,
  "primary score models", nrow(score_statistics), 9, nrow(score_statistics) == 9,
  "integrated score beta", integrated_beta, 0.256754715641472, near(integrated_beta, 0.256754715641472),
  "score event-QC models", nrow(event_score_statistics), 9, nrow(event_score_statistics) == 9,
  "integrated event-QC N", event_score_statistics$n_model[event_score_statistics$feature == CP28_MANIFEST$integrated_score][1], 808, event_score_statistics$n_model[event_score_statistics$feature == CP28_MANIFEST$integrated_score][1] == 808,
  "age-stratified models", nrow(age_stratified), 156, nrow(age_stratified) == 156,
  "age target family", dplyr::n_distinct(age_stratified$feature), 52, dplyr::n_distinct(age_stratified$feature) == 52,
  "5-year matched rows", nrow(pairs_5), 530, nrow(pairs_5) == 530,
  "10-year matched rows", nrow(pairs_10), 546, nrow(pairs_10) == 546,
  "5-year matched pairs", dplyr::n_distinct(pairs_5$pair_id), 265, dplyr::n_distinct(pairs_5$pair_id) == 265,
  "10-year matched pairs", dplyr::n_distinct(pairs_10$pair_id), 273, dplyr::n_distinct(pairs_10$pair_id) == 273,
  "matched models", nrow(matched_statistics), 104, nrow(matched_statistics) == 104,
  "matching summary rows", nrow(matching_summary), 2, nrow(matching_summary) == 2,
  "interaction models", nrow(interactions), 52, nrow(interactions) == 52,
  "threshold scored rows", nrow(threshold_scored), 2550, nrow(threshold_scored) == 2550,
  "threshold rows per set", min(table(threshold_scored$threshold_set)), 850, all(table(threshold_scored$threshold_set) == 850),
  "threshold models", nrow(threshold_statistics), 240, nrow(threshold_statistics) == 240,
  "threshold variables", nrow(threshold_summary), 80, nrow(threshold_summary) == 80,
  "threshold direction+FDR robust", count_class("Direction and FDR-significance preserved across threshold sets"), 35, count_class("Direction and FDR-significance preserved across threshold sets") == 35,
  "threshold direction-only", count_class("Direction preserved; FDR-significance not preserved in all threshold sets"), 34, count_class("Direction preserved; FDR-significance not preserved in all threshold sets") == 34,
  "threshold direction failures", count_class("Direction not preserved across threshold sets"), 11, count_class("Direction not preserved across threshold sets") == 11
)

reference_dir <- path.expand(Sys.getenv("SDY2583_CP28_REFERENCE_DIR", unset = ""))
if (nzchar(reference_dir) && dir.exists(reference_dir)) {
  comparison_plan <- tibble::tribble(
    ~generated, ~reference_name, ~keys,
    file.path(analysis_dir, "02_feature_extraction", "SDY2583_CP28_FULL_850_feature_table_STEP2_TNK_interface_RECONSTRUCTED.csv"), "SDY2583_CP28_FULL_850_feature_table_STEP2_TNK_interface.csv", "subject_id",
    file.path(analysis_dir, "04_metadata_merge_age_QC", "SDY2583_CP28_analysis_data_with_metadata_RECONSTRUCTED.csv"), "SDY2583_CP28_analysis_data_STEP3A_metadata_merge_age_QC.csv", "subject_id",
    file.path(analysis_dir, "06_statistics_age_sex_adjusted", "SDY2583_CP28_age_sex_adjusted_statistics_STEP3B_RECONSTRUCTED.csv"), "SDY2583_CP28_age_sex_adjusted_statistics_STEP3B.csv", "module,feature",
    file.path(analysis_dir, "06_statistics_age_sex_adjusted", "SDY2583_CP28_binary_sex_sensitivity_statistics_STEP3B_RECONSTRUCTED.csv"), "SDY2583_CP28_binary_sex_sensitivity_statistics_STEP3B.csv", "module,feature",
    file.path(analysis_dir, "06_statistics_age_sex_adjusted", "SDY2583_CP28_event_QC_sensitivity_statistics_STEP3B_RECONSTRUCTED.csv"), "SDY2583_CP28_event_QC_sensitivity_statistics_STEP3B.csv", "module,feature",
    file.path(analysis_dir, "08_composite_scores", "SDY2583_CP28_analysis_data_with_composite_scores_RECONSTRUCTED.csv"), "SDY2583_CP28_analysis_data_with_composite_scores.csv", "subject_id",
    file.path(analysis_dir, "08_composite_scores", "SDY2583_CP28_composite_score_feature_sets_RECONSTRUCTED.csv"), "SDY2583_CP28_composite_score_feature_sets.csv", "composite_score,feature",
    file.path(analysis_dir, "08_composite_scores", "SDY2583_CP28_composite_score_statistics_RECONSTRUCTED.csv"), "SDY2583_CP28_composite_score_statistics.csv", "feature",
    file.path(analysis_dir, "08_composite_scores", "SDY2583_CP28_composite_score_event_QC_sensitivity_statistics_RECONSTRUCTED.csv"), "SDY2583_CP28_composite_score_event_QC_sensitivity_statistics.csv", "feature",
    file.path(analysis_dir, "10_age_sensitivity", "SDY2583_CP28_age_stratified_statistics_RECONSTRUCTED.csv"), "SDY2583_CP28_age_stratified_statistics.csv", "age_group_for_model,feature",
    file.path(analysis_dir, "10_age_sensitivity", "SDY2583_CP28_same_sex_age_matched_pairs_5yr_RECONSTRUCTED.csv"), "SDY2583_CP28_same_sex_age_matched_pairs_5yr.csv", "pair_id,subject_id",
    file.path(analysis_dir, "10_age_sensitivity", "SDY2583_CP28_same_sex_age_matched_pairs_10yr_RECONSTRUCTED.csv"), "SDY2583_CP28_same_sex_age_matched_pairs_10yr.csv", "pair_id,subject_id",
    file.path(analysis_dir, "10_age_sensitivity", "SDY2583_CP28_age_matched_sensitivity_statistics_RECONSTRUCTED.csv"), "SDY2583_CP28_age_matched_sensitivity_statistics.csv", "caliper_years,feature",
    file.path(analysis_dir, "10_age_sensitivity", "SDY2583_CP28_disease_by_age_interaction_statistics_RECONSTRUCTED.csv"), "SDY2583_CP28_disease_by_age_interaction_statistics.csv", "feature",
    file.path(analysis_dir, "12_threshold_sensitivity", "SDY2583_CP28_threshold_sensitivity_scored_data_LONG_RECONSTRUCTED.csv"), "SDY2583_CP28_threshold_sensitivity_scored_data_LONG.csv", "threshold_set,subject_id",
    file.path(analysis_dir, "12_threshold_sensitivity", "SDY2583_CP28_threshold_sensitivity_statistics_LONG_RECONSTRUCTED.csv"), "SDY2583_CP28_threshold_sensitivity_statistics_LONG.csv", "threshold_set,feature",
    file.path(analysis_dir, "12_threshold_sensitivity", "SDY2583_CP28_threshold_sensitivity_summary_WIDE_RECONSTRUCTED.csv"), "SDY2583_CP28_threshold_sensitivity_summary_WIDE.csv", "feature"
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

report_file <- file.path(out_dir, "SDY2583_CP28_reconstruction_validation_report.csv")
readr::write_csv(checks, report_file)
print(checks, n = nrow(checks))
if (any(!checks$pass | is.na(checks$pass))) stop("CP28 validation failed. Review ", report_file)
message("CP28 validation passed.")
