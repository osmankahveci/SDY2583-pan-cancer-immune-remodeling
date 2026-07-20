# ============================================================
# Validate the complete reconstructed CP7 pipeline.
#
# Fixed checks run without the archive. When SDY2583_CP7_REFERENCE_DIR points
# to the archived CP7 analysis tree, exact table-level comparisons are added.
# ============================================================

source(file.path(
  Sys.getenv("SDY2583_REPO_ROOT", unset = "."),
  "R", "shared", "reconstructed_panel_framework.R"
))
rp_install_and_load(c("dplyr", "readr", "tibble", "purrr", "stringr"))
source(file.path(sd_repo_root(), "R", "panels", "CP7", "MANIFEST_RECONSTRUCTED.R"))

analysis_dir <- sd_analysis_dir("CP7")
validation_dir <- file.path(analysis_dir, "validation")
dir.create(validation_dir, recursive = TRUE, showWarnings = FALSE)
reference_dir <- path.expand(Sys.getenv("SDY2583_CP7_REFERENCE_DIR", unset = ""))
strict <- tolower(Sys.getenv("SDY2583_VALIDATION_STRICT", unset = "true")) %in%
  c("true", "1", "yes", "y")
tolerance <- suppressWarnings(as.numeric(
  Sys.getenv("SDY2583_VALIDATION_TOLERANCE", unset = "1e-6")
))
if (!is.finite(tolerance) || tolerance <= 0) tolerance <- 1e-6

read_required <- function(path) {
  if (!file.exists(path)) stop("Missing reconstructed output: ", path)
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
}
near <- function(observed, expected) {
  length(observed) == 1L && length(expected) == 1L &&
    is.finite(observed) && is.finite(expected) &&
    abs(observed - expected) <= tolerance * max(1, abs(expected))
}

report_rows <- list()
add_check <- function(check, observed, expected, pass, detail = "") {
  report_rows[[length(report_rows) + 1L]] <<- tibble::tibble(
    check = check,
    observed = as.character(observed),
    expected = as.character(expected),
    pass = isTRUE(pass),
    detail = detail
  )
}

step1 <- read_required(file.path(
  analysis_dir, "01_channel_marker_QC", "SDY2583_CP7_STEP1_QC_summary.csv"
))
step2 <- read_required(file.path(
  analysis_dir, "02_feature_extraction",
  "SDY2583_CP7_FULL_850_STEP2_QC_summary_RECONSTRUCTED.csv"
))
metadata <- read_required(file.path(
  analysis_dir, "04_metadata_merge",
  "SDY2583_CP7_metadata_merge_summary_RECONSTRUCTED.csv"
))
main_statistics <- read_required(file.path(
  analysis_dir, "06_statistics",
  "SDY2583_CP7_FULL850_age_sex_adjusted_statistics_MAIN_RECONSTRUCTED.csv"
))
composite_data <- read_required(file.path(
  analysis_dir, "08_composite_scores",
  "SDY2583_CP7_analysis_data_with_composite_scores_RECONSTRUCTED.csv"
))
composite_statistics <- read_required(file.path(
  analysis_dir, "08_composite_scores",
  "SDY2583_CP7_composite_score_statistics_RECONSTRUCTED.csv"
))
match_summary <- read_required(file.path(
  analysis_dir, "10_age_sensitivity",
  "SDY2583_CP7_age_matched_sensitivity_summary_RECONSTRUCTED.csv"
))
age_feature_map <- read_required(file.path(
  analysis_dir, "10_age_sensitivity",
  "SDY2583_CP7_age_sensitivity_feature_map_RECONSTRUCTED.csv"
))
threshold_statistics <- read_required(file.path(
  analysis_dir, "12_threshold_sensitivity",
  "SDY2583_CP7_threshold_sensitivity_statistics_LONG_RECONSTRUCTED.csv"
))
threshold_summary <- read_required(file.path(
  analysis_dir, "12_threshold_sensitivity",
  "SDY2583_CP7_threshold_sensitivity_summary_WIDE_RECONSTRUCTED.csv"
))

get_main <- function(feature, column = "beta_cancer_vs_healthy") {
  row <- main_statistics[main_statistics$feature == feature, , drop = FALSE]
  if (nrow(row) != 1L || !(column %in% names(row))) return(NA_real_)
  as.numeric(row[[column]][1])
}
get_score <- function(score, column = "beta_cancer_vs_healthy") {
  row <- composite_statistics[composite_statistics$score == score, , drop = FALSE]
  if (nrow(row) != 1L || !(column %in% names(row))) return(NA_real_)
  as.numeric(row[[column]][1])
}
get_pairs <- function(caliper) {
  row <- match_summary[match_summary$caliper_years == caliper, , drop = FALSE]
  if (nrow(row) != 1L) return(NA_integer_)
  as.integer(row$n_pairs[1])
}

add_check("FCS files", step1$n_fcs_files[1], 850, step1$n_fcs_files[1] == 850)
add_check(
  "Unique subjects", step1$n_unique_subjects[1], 850,
  step1$n_unique_subjects[1] == 850
)
add_check(
  "Successful Step 2 extractions", step2$feature_ok[1], 850,
  step2$feature_ok[1] == 850
)
add_check(
  "Failed Step 2 extractions", step2$feature_failed[1], 0,
  step2$feature_failed[1] == 0
)
add_check(
  "Median total events", step2$median_total_events[1], 301989.5,
  near(step2$median_total_events[1], 301989.5)
)
add_check(
  "Median CD3+CD8+ events", step2$median_cd3cd8_events[1], 12411.5,
  near(step2$median_cd3cd8_events[1], 12411.5)
)
add_check(
  "Minimum CD3+CD8+ events", step2$min_cd3cd8_events[1], 407,
  step2$min_cd3cd8_events[1] == 407
)
add_check(
  "Maximum CD3+CD8+ events", step2$max_cd3cd8_events[1], 59837,
  step2$max_cd3cd8_events[1] == 59837
)
add_check(
  "Valid ages", metadata$n_valid_age[1], 832,
  metadata$n_valid_age[1] == 832
)
add_check("Main feature count", nrow(main_statistics), 55, nrow(main_statistics) == 55)
add_check(
  "Primary composite subject count", nrow(composite_data), 832,
  nrow(composite_data) == 832,
  "Primary scores must be standardized in the model-ready set."
)
add_check(
  "Composite score count", nrow(composite_statistics), 6,
  nrow(composite_statistics) == 6
)
add_check(
  "Age-sensitivity outcome count", nrow(age_feature_map), 32,
  nrow(age_feature_map) == 32
)
add_check(
  "Five-year matched pairs", get_pairs(5),
  CP7_MANIFEST$benchmarks$matched_pairs_5y,
  get_pairs(5) == CP7_MANIFEST$benchmarks$matched_pairs_5y
)
add_check(
  "Ten-year matched pairs", get_pairs(10),
  CP7_MANIFEST$benchmarks$matched_pairs_10y,
  get_pairs(10) == CP7_MANIFEST$benchmarks$matched_pairs_10y
)
add_check(
  "Threshold target count", nrow(threshold_summary),
  CP7_MANIFEST$benchmarks$threshold_target_count,
  nrow(threshold_summary) == CP7_MANIFEST$benchmarks$threshold_target_count
)
add_check(
  "Threshold long-table row count", nrow(threshold_statistics), 33,
  nrow(threshold_statistics) == 33
)
add_check(
  "Threshold directions preserved", sum(threshold_summary$direction_preserved), 11,
  all(threshold_summary$direction_preserved)
)
add_check(
  "Threshold FDR preserved", sum(threshold_summary$FDR_lt_0p05_all_sets), 11,
  all(threshold_summary$FDR_lt_0p05_all_sets)
)

anchor_expectations <- tibble::tribble(
  ~label, ~feature, ~column, ~expected,
  "CD3 beta", "pct_cd3_pos_total", "beta_cancer_vs_healthy", -4.71800773276706,
  "CD3 global FDR", "pct_cd3_pos_total", "fdr_all", 4.54270592950129e-17,
  "CD8 within CD3 beta", "pct_cd8_within_cd3", "beta_cancer_vs_healthy", 4.28759233054575,
  "TEMRA beta", "pct_temra_like", "beta_cancer_vs_healthy", 4.36723790655398,
  "CD39-positive beta", "pct_cd39_pos", "beta_cancer_vs_healthy", 3.23434446149965
)
for (i in seq_len(nrow(anchor_expectations))) {
  observed <- get_main(
    anchor_expectations$feature[i],
    anchor_expectations$column[i]
  )
  add_check(
    anchor_expectations$label[i], observed,
    anchor_expectations$expected[i],
    near(observed, anchor_expectations$expected[i])
  )
}

integrated_beta <- get_score(CP7_MANIFEST$integrated_score)
add_check(
  "Integrated primary composite beta", integrated_beta,
  CP7_MANIFEST$benchmarks$integrated_beta,
  near(integrated_beta, CP7_MANIFEST$benchmarks$integrated_beta),
  "The 832-subject primary-score standardization is required."
)

find_reference <- function(basename_expected) {
  if (!nzchar(reference_dir) || !dir.exists(reference_dir)) return(NA_character_)
  hits <- list.files(
    reference_dir,
    pattern = paste0("^", basename_expected, "$"),
    recursive = TRUE,
    full.names = TRUE
  )
  if (length(hits) == 0L) NA_character_ else hits[1]
}
compare_archive <- function(label, generated, reference_basename, keys) {
  reference <- find_reference(reference_basename)
  if (is.na(reference)) {
    add_check(
      paste0("Archive table: ", label), "reference not found",
      reference_basename, FALSE,
      "The reference directory was provided but this table was not found."
    )
    return(invisible(NULL))
  }
  result <- rp_compare_csv(generated, reference, keys, tolerance)
  add_check(
    paste0("Archive table: ", label), result$detail[1],
    "matched within tolerance", result$pass[1],
    paste("Generated:", generated, "Reference:", reference)
  )
}

if (nzchar(reference_dir) && dir.exists(reference_dir)) {
  table_checks <- list(
    list(
      "Step 2 feature table",
      file.path(analysis_dir, "02_feature_extraction", "SDY2583_CP7_FULL_850_feature_table_STEP2_RECONSTRUCTED.csv"),
      "SDY2583_CP7_FULL_850_feature_table_STEP2.csv", c("subject_id")
    ),
    list(
      "Main adjusted statistics",
      file.path(analysis_dir, "06_statistics", "SDY2583_CP7_FULL850_age_sex_adjusted_statistics_MAIN_RECONSTRUCTED.csv"),
      "SDY2583_CP7_FULL850_age_sex_adjusted_statistics_MAIN.csv", c("module", "feature")
    ),
    list(
      "Composite subject scores",
      file.path(analysis_dir, "08_composite_scores", "SDY2583_CP7_analysis_data_with_composite_scores_RECONSTRUCTED.csv"),
      "SDY2583_CP7_analysis_data_with_composite_scores.csv", c("subject_id")
    ),
    list(
      "Composite statistics",
      file.path(analysis_dir, "08_composite_scores", "SDY2583_CP7_composite_score_statistics_RECONSTRUCTED.csv"),
      "SDY2583_CP7_composite_score_statistics.csv", c("score")
    ),
    list(
      "Age feature map",
      file.path(analysis_dir, "10_age_sensitivity", "SDY2583_CP7_age_sensitivity_feature_map_RECONSTRUCTED.csv"),
      "SDY2583_CP7_age_sensitivity_feature_map.csv", c("module", "feature")
    ),
    list(
      "Age-stratified statistics",
      file.path(analysis_dir, "10_age_sensitivity", "SDY2583_CP7_age_stratified_statistics_RECONSTRUCTED.csv"),
      "SDY2583_CP7_age_stratified_statistics.csv", c("subset_label", "module", "feature")
    ),
    list(
      "Matched summary",
      file.path(analysis_dir, "10_age_sensitivity", "SDY2583_CP7_age_matched_sensitivity_summary_RECONSTRUCTED.csv"),
      "SDY2583_CP7_age_matched_sensitivity_summary.csv", c("subset_label")
    ),
    list(
      "Matched statistics",
      file.path(analysis_dir, "10_age_sensitivity", "SDY2583_CP7_age_matched_sensitivity_statistics_RECONSTRUCTED.csv"),
      "SDY2583_CP7_age_matched_sensitivity_statistics.csv", c("subset_label", "module", "feature")
    ),
    list(
      "Matched 5-year subjects",
      file.path(analysis_dir, "10_age_sensitivity", "SDY2583_CP7_age_matched_dataset_caliper_5y_RECONSTRUCTED.csv"),
      "SDY2583_CP7_age_matched_dataset_caliper_5y.csv", c("matched_pair_id", "disease_group")
    ),
    list(
      "Matched 10-year subjects",
      file.path(analysis_dir, "10_age_sensitivity", "SDY2583_CP7_age_matched_dataset_caliper_10y_RECONSTRUCTED.csv"),
      "SDY2583_CP7_age_matched_dataset_caliper_10y.csv", c("matched_pair_id", "disease_group")
    ),
    list(
      "Disease-by-age interactions",
      file.path(analysis_dir, "10_age_sensitivity", "SDY2583_CP7_disease_by_age_interaction_statistics_RECONSTRUCTED.csv"),
      "SDY2583_CP7_disease_by_age_interaction_statistics.csv", c("module", "feature")
    ),
    list(
      "Threshold statistics",
      file.path(analysis_dir, "12_threshold_sensitivity", "SDY2583_CP7_threshold_sensitivity_statistics_LONG_RECONSTRUCTED.csv"),
      "SDY2583_CP7_threshold_sensitivity_statistics_LONG.csv", c("threshold_set", "feature")
    ),
    list(
      "Threshold summary",
      file.path(analysis_dir, "12_threshold_sensitivity", "SDY2583_CP7_threshold_sensitivity_summary_WIDE_RECONSTRUCTED.csv"),
      "SDY2583_CP7_threshold_sensitivity_summary_WIDE.csv", c("feature")
    ),
    list(
      "Threshold subject scores",
      file.path(analysis_dir, "12_threshold_sensitivity", "SDY2583_CP7_threshold_sensitivity_scored_data_LONG_RECONSTRUCTED.csv"),
      "SDY2583_CP7_threshold_sensitivity_scored_data_LONG.csv", c("threshold_set", "subject_id")
    )
  )
  purrr::walk(table_checks, ~do.call(compare_archive, .x))
}

validation_report <- dplyr::bind_rows(report_rows)
readr::write_csv(
  validation_report,
  file.path(validation_dir, "SDY2583_CP7_reconstruction_validation_report.csv")
)
print(validation_report, n = nrow(validation_report))
cat("\nPassed:", sum(validation_report$pass), "of", nrow(validation_report), "checks.\n")
if (strict && any(!validation_report$pass)) {
  stop("CP7 reconstruction validation failed. Review the validation report.")
}
