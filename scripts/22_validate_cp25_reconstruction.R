# Validate reconstructed CP25 outputs against fixed archive benchmarks and
# optional archived CSV tables.

source(file.path(
  Sys.getenv("SDY2583_REPO_ROOT", unset = "."),
  "R", "shared", "reconstructed_panel_framework.R"
))
rp_install_and_load(c("dplyr", "readr", "tibble", "purrr"))
source(file.path(
  sd_repo_root(), "R", "panels", "CP25", "MANIFEST_RECONSTRUCTED.R"
))

analysis_dir <- sd_analysis_dir("CP25")
out_dir <- file.path(analysis_dir, "validation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

tolerance <- suppressWarnings(as.numeric(Sys.getenv(
  "SDY2583_VALIDATION_TOLERANCE",
  unset = "1e-6"
)))
if (!is.finite(tolerance)) tolerance <- 1e-6
near <- function(x, y, tol = tolerance) {
  is.finite(x) && abs(x - y) <= tol * max(1, abs(y))
}
read_required <- function(path) {
  if (!file.exists(path)) stop("Missing CP25 output: ", path)
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
}

feature_summary <- read_required(file.path(
  analysis_dir, "02_feature_extraction",
  "SDY2583_CP25_feature_extraction_summary_STEP2_RECONSTRUCTED.csv"
))
metadata_summary <- read_required(file.path(
  analysis_dir, "04_metadata_merge",
  "SDY2583_CP25_metadata_merge_summary_RECONSTRUCTED.csv"
))
main_statistics <- read_required(file.path(
  analysis_dir, "06_statistics",
  "SDY2583_CP25_FULL850_age_sex_adjusted_statistics_MAIN_RECONSTRUCTED.csv"
))
binary_statistics <- read_required(file.path(
  analysis_dir, "07_sensitivity",
  "SDY2583_CP25_binary_sex_sensitivity_statistics_RECONSTRUCTED.csv"
))
event_statistics <- read_required(file.path(
  analysis_dir, "07_sensitivity",
  "SDY2583_CP25_event_QC_ge300_sensitivity_statistics_RECONSTRUCTED.csv"
))
score_data <- read_required(file.path(
  analysis_dir, "08_composite_scores",
  "SDY2583_CP25_analysis_data_with_composite_scores_RECONSTRUCTED.csv"
))
score_statistics <- read_required(file.path(
  analysis_dir, "08_composite_scores",
  "SDY2583_CP25_composite_score_statistics_RECONSTRUCTED.csv"
))
age_stratified <- read_required(file.path(
  analysis_dir, "10_age_sensitivity",
  "SDY2583_CP25_age_stratified_statistics_RECONSTRUCTED.csv"
))
matched_statistics <- read_required(file.path(
  analysis_dir, "10_age_sensitivity",
  "SDY2583_CP25_age_matched_sensitivity_statistics_RECONSTRUCTED.csv"
))
match_summary <- read_required(file.path(
  analysis_dir, "10_age_sensitivity",
  "SDY2583_CP25_age_matched_sensitivity_summary_RECONSTRUCTED.csv"
))
pairs_5 <- read_required(file.path(
  analysis_dir, "10_age_sensitivity",
  "SDY2583_CP25_same_sex_age_matched_pairs_5yr_RECONSTRUCTED.csv"
))
pairs_10 <- read_required(file.path(
  analysis_dir, "10_age_sensitivity",
  "SDY2583_CP25_same_sex_age_matched_pairs_10yr_RECONSTRUCTED.csv"
))
interactions <- read_required(file.path(
  analysis_dir, "10_age_sensitivity",
  "SDY2583_CP25_disease_by_age_interaction_statistics_RECONSTRUCTED.csv"
))
threshold_statistics <- read_required(file.path(
  analysis_dir, "12_threshold_sensitivity",
  "SDY2583_CP25_threshold_sensitivity_statistics_LONG_RECONSTRUCTED.csv"
))
threshold_summary <- read_required(file.path(
  analysis_dir, "12_threshold_sensitivity",
  "SDY2583_CP25_threshold_sensitivity_summary_WIDE_RECONSTRUCTED.csv"
))

integrated_beta <- score_statistics$beta_cancer_vs_healthy[
  score_statistics$score == CP25_MANIFEST$integrated_score
][1]
robust_counts <- table(threshold_summary$robust_interpretation)
count_class <- function(label) {
  value <- unname(robust_counts[label])
  if (length(value) == 0L || is.na(value)) 0L else as.integer(value)
}

checks <- tibble::tribble(
  ~check, ~observed, ~expected, ~pass,
  "FCS files", as.character(feature_summary$n_files[1]), "850", feature_summary$n_files[1] == 850,
  "successful feature extraction", as.character(feature_summary$n_success[1]), "850", feature_summary$n_success[1] == 850,
  "median total events", as.character(feature_summary$median_total_events[1]), "302748", near(feature_summary$median_total_events[1], 302748),
  "median raw CD3+CD4+ events", as.character(feature_summary$median_raw_cd3cd4_events[1]), "26488", near(feature_summary$median_raw_cd3cd4_events[1], 26488),
  "median primary-gate events", as.character(feature_summary$median_primary_events[1]), "16027", near(feature_summary$median_primary_events[1], 16027),
  "valid age/model-ready subjects", as.character(metadata_summary$n_model_ready[1]), "832", metadata_summary$n_model_ready[1] == 832,
  "main feature models", as.character(nrow(main_statistics)), "77", nrow(main_statistics) == 77,
  "binary-sex feature models", as.character(nrow(binary_statistics)), "77", nrow(binary_statistics) == 77,
  "event-QC feature models", as.character(nrow(event_statistics)), "77", nrow(event_statistics) == 77,
  "subject-level score rows", as.character(nrow(score_data)), "850", nrow(score_data) == 850,
  "composite score models", as.character(nrow(score_statistics)), "9", nrow(score_statistics) == 9,
  "all composite scores FDR-significant", as.character(sum(score_statistics$fdr < 0.05, na.rm = TRUE)), "9", sum(score_statistics$fdr < 0.05, na.rm = TRUE) == 9,
  "integrated score beta", as.character(integrated_beta), "0.4392101623", near(integrated_beta, 0.4392101623),
  "age-stratified models", as.character(nrow(age_stratified)), "258", nrow(age_stratified) == 258,
  "matched models", as.character(nrow(matched_statistics)), "172", nrow(matched_statistics) == 172,
  "5-year pairs", as.character(nrow(pairs_5) / 2), "265", nrow(pairs_5) == 530 && dplyr::n_distinct(pairs_5$pair_id) == 265,
  "10-year pairs", as.character(nrow(pairs_10) / 2), "273", nrow(pairs_10) == 546 && dplyr::n_distinct(pairs_10$pair_id) == 273,
  "5-year summary pairs", as.character(match_summary$n_pairs[match_summary$caliper_years == 5][1]), "265", match_summary$n_pairs[match_summary$caliper_years == 5][1] == 265,
  "10-year summary pairs", as.character(match_summary$n_pairs[match_summary$caliper_years == 10][1]), "273", match_summary$n_pairs[match_summary$caliper_years == 10][1] == 273,
  "interaction models", as.character(nrow(interactions)), "86", nrow(interactions) == 86,
  "threshold models", as.character(nrow(threshold_statistics)), "219", nrow(threshold_statistics) == 219,
  "threshold outcomes", as.character(nrow(threshold_summary)), "73", nrow(threshold_summary) == 73,
  "threshold target feature set", as.character(sum(CP25_MANIFEST$threshold_target_features %in% threshold_summary$feature)), "64", all(CP25_MANIFEST$threshold_target_features %in% threshold_summary$feature),
  "threshold direction+FDR robust", as.character(count_class("Direction and FDR-significance preserved across threshold sets")), "48", count_class("Direction and FDR-significance preserved across threshold sets") == 48,
  "threshold direction-only", as.character(count_class("Direction preserved; FDR-significance not preserved in all threshold sets")), "23", count_class("Direction preserved; FDR-significance not preserved in all threshold sets") == 23,
  "threshold direction failures", as.character(count_class("Direction not preserved across threshold sets")), "2", count_class("Direction not preserved across threshold sets") == 2
)

reference_dir <- path.expand(Sys.getenv("SDY2583_CP25_REFERENCE_DIR", unset = ""))
if (nzchar(reference_dir) && dir.exists(reference_dir)) {
  plan <- tibble::tribble(
    ~generated, ~reference_name, ~keys,
    file.path(analysis_dir, "04_metadata_merge", "SDY2583_CP25_analysis_data_with_metadata_RECONSTRUCTED.csv"), "SDY2583_CP25_analysis_data_STEP3A_metadata_merged.csv", "subject_id",
    file.path(analysis_dir, "06_statistics", "SDY2583_CP25_FULL850_age_sex_adjusted_statistics_MAIN_RECONSTRUCTED.csv"), "SDY2583_CP25_FULL850_age_sex_adjusted_statistics_MAIN.csv", "module,feature",
    file.path(analysis_dir, "07_sensitivity", "SDY2583_CP25_binary_sex_sensitivity_statistics_RECONSTRUCTED.csv"), "SDY2583_CP25_binary_sex_sensitivity_statistics.csv", "module,feature",
    file.path(analysis_dir, "07_sensitivity", "SDY2583_CP25_event_QC_ge300_sensitivity_statistics_RECONSTRUCTED.csv"), "SDY2583_CP25_event_QC_ge300_sensitivity_statistics.csv", "module,feature",
    file.path(analysis_dir, "08_composite_scores", "SDY2583_CP25_analysis_data_with_composite_scores_RECONSTRUCTED.csv"), "SDY2583_CP25_analysis_data_with_composite_scores.csv", "subject_id",
    file.path(analysis_dir, "08_composite_scores", "SDY2583_CP25_composite_score_statistics_RECONSTRUCTED.csv"), "SDY2583_CP25_composite_score_statistics.csv", "score,feature",
    file.path(analysis_dir, "10_age_sensitivity", "SDY2583_CP25_age_stratified_statistics_RECONSTRUCTED.csv"), "SDY2583_CP25_age_stratified_statistics.csv", "age_group_for_model,feature",
    file.path(analysis_dir, "10_age_sensitivity", "SDY2583_CP25_same_sex_age_matched_pairs_5yr_RECONSTRUCTED.csv"), "SDY2583_CP25_same_sex_age_matched_pairs_5yr.csv", "pair_id,subject_id",
    file.path(analysis_dir, "10_age_sensitivity", "SDY2583_CP25_same_sex_age_matched_pairs_10yr_RECONSTRUCTED.csv"), "SDY2583_CP25_same_sex_age_matched_pairs_10yr.csv", "pair_id,subject_id",
    file.path(analysis_dir, "10_age_sensitivity", "SDY2583_CP25_age_matched_sensitivity_statistics_RECONSTRUCTED.csv"), "SDY2583_CP25_age_matched_sensitivity_statistics.csv", "caliper_years,feature",
    file.path(analysis_dir, "10_age_sensitivity", "SDY2583_CP25_disease_by_age_interaction_statistics_RECONSTRUCTED.csv"), "SDY2583_CP25_disease_by_age_interaction_statistics.csv", "feature",
    file.path(analysis_dir, "12_threshold_sensitivity", "SDY2583_CP25_threshold_sensitivity_statistics_LONG_RECONSTRUCTED.csv"), "SDY2583_CP25_threshold_sensitivity_statistics_LONG.csv", "threshold_set,feature",
    file.path(analysis_dir, "12_threshold_sensitivity", "SDY2583_CP25_threshold_sensitivity_summary_WIDE_RECONSTRUCTED.csv"), "SDY2583_CP25_threshold_sensitivity_summary_WIDE.csv", "feature"
  )
  references <- list.files(reference_dir, "\\.csv$", full.names = TRUE, recursive = TRUE)
  archive_checks <- purrr::pmap_dfr(plan, function(generated, reference_name, keys) {
    hit <- references[basename(references) == reference_name]
    comparison <- if (length(hit)) {
      rp_compare_csv(
        generated,
        hit[1],
        strsplit(keys, ",", fixed = TRUE)[[1]],
        tolerance = tolerance
      )
    } else {
      tibble::tibble(pass = FALSE, detail = "reference missing")
    }
    tibble::tibble(
      check = paste0("archive table: ", reference_name),
      observed = comparison$detail[1],
      expected = "matched",
      pass = comparison$pass[1]
    )
  })
  checks <- dplyr::bind_rows(checks, archive_checks)
}

report_file <- file.path(
  out_dir,
  "SDY2583_CP25_reconstruction_validation_report.csv"
)
readr::write_csv(checks, report_file)
print(checks, n = nrow(checks))
if (any(!checks$pass | is.na(checks$pass))) {
  stop("CP25 validation failed. Review ", report_file)
}
message("CP25 validation passed.")
