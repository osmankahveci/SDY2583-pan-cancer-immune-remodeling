# ============================================================
# Validate the reconstructed CP8 pipeline against fixed archive benchmarks
# and, when SDY2583_CP8_REFERENCE_DIR is configured, archived CSV tables.
# ============================================================

source(file.path(
  Sys.getenv("SDY2583_REPO_ROOT", unset = "."),
  "R", "shared", "reconstructed_panel_framework.R"
))
rp_install_and_load(c("dplyr", "readr", "tibble", "purrr", "stringr"))
source(file.path(sd_repo_root(), "R", "panels", "CP8", "MANIFEST_RECONSTRUCTED.R"))

analysis_dir <- sd_analysis_dir("CP8")
out_dir <- file.path(analysis_dir, "validation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

strict <- tolower(
  Sys.getenv("SDY2583_VALIDATION_STRICT", unset = "true")
) %in% c("true", "1", "yes")
tolerance <- suppressWarnings(as.numeric(
  Sys.getenv("SDY2583_VALIDATION_TOLERANCE", unset = "1e-6")
))
if (!is.finite(tolerance)) tolerance <- 1e-6

read_required <- function(path) {
  if (!file.exists(path)) stop("Missing reconstructed CP8 output: ", path)
  readr::read_csv(path, show_col_types = FALSE)
}

near <- function(x, y, tol = tolerance) {
  length(x) > 0L && is.finite(x[1]) &&
    abs(x[1] - y) <= tol * max(1, abs(y))
}

step2 <- read_required(file.path(
  analysis_dir, "02_feature_extraction",
  "SDY2583_CP8_feature_extraction_summary_STEP2_RECONSTRUCTED.csv"
))
metadata_summary <- read_required(file.path(
  analysis_dir, "04_metadata_merge",
  "SDY2583_CP8_metadata_merge_summary_RECONSTRUCTED.csv"
))
main_statistics <- read_required(file.path(
  analysis_dir, "06_statistics",
  "SDY2583_CP8_FULL850_age_sex_adjusted_statistics_MAIN_RECONSTRUCTED.csv"
))
binary_statistics <- read_required(file.path(
  analysis_dir, "07_sensitivity",
  "SDY2583_CP8_binary_sex_sensitivity_statistics_RECONSTRUCTED.csv"
))
scored_data <- read_required(file.path(
  analysis_dir, "08_composite_scores",
  "SDY2583_CP8_analysis_data_with_composite_scores_RECONSTRUCTED.csv"
))
score_statistics <- read_required(file.path(
  analysis_dir, "08_composite_scores",
  "SDY2583_CP8_composite_score_statistics_RECONSTRUCTED.csv"
))
age_stratified <- read_required(file.path(
  analysis_dir, "10_age_sensitivity",
  "SDY2583_CP8_age_stratified_statistics_RECONSTRUCTED.csv"
))
match_summary <- read_required(file.path(
  analysis_dir, "10_age_sensitivity",
  "SDY2583_CP8_age_matched_sensitivity_summary_RECONSTRUCTED.csv"
))
matched_statistics <- read_required(file.path(
  analysis_dir, "10_age_sensitivity",
  "SDY2583_CP8_age_matched_sensitivity_statistics_RECONSTRUCTED.csv"
))
interaction_statistics <- read_required(file.path(
  analysis_dir, "10_age_sensitivity",
  "SDY2583_CP8_disease_by_age_interaction_statistics_RECONSTRUCTED.csv"
))
pairs_5y <- read_required(file.path(
  analysis_dir, "10_age_sensitivity",
  "SDY2583_CP8_same_sex_age_matched_pairs_5yr_RECONSTRUCTED.csv"
))
pairs_10y <- read_required(file.path(
  analysis_dir, "10_age_sensitivity",
  "SDY2583_CP8_same_sex_age_matched_pairs_10yr_RECONSTRUCTED.csv"
))
threshold_statistics <- read_required(file.path(
  analysis_dir, "12_threshold_sensitivity",
  "SDY2583_CP8_threshold_sensitivity_statistics_LONG_RECONSTRUCTED.csv"
))
threshold_summary <- read_required(file.path(
  analysis_dir, "12_threshold_sensitivity",
  "SDY2583_CP8_threshold_sensitivity_summary_WIDE_RECONSTRUCTED.csv"
))
threshold_scored <- read_required(file.path(
  analysis_dir, "12_threshold_sensitivity",
  "SDY2583_CP8_threshold_sensitivity_scored_data_LONG_RECONSTRUCTED.csv"
))

score_key <- if ("score" %in% names(score_statistics)) "score" else "feature"
integrated_beta <- score_statistics$beta_cancer_vs_healthy[
  score_statistics[[score_key]] == CP8_MANIFEST$integrated_score
][1]

checks <- tibble::tribble(
  ~check, ~observed, ~expected, ~pass,
  "FCS files", step2$n_files[1], 850, step2$n_files[1] == 850,
  "successful extractions", step2$n_success[1], 850, step2$n_success[1] == 850,
  "failed extractions", step2$n_failed[1], 0, step2$n_failed[1] == 0,
  "median total events", step2$median_total_events[1],
    CP8_MANIFEST$benchmarks$median_total_events,
    near(step2$median_total_events[1], CP8_MANIFEST$benchmarks$median_total_events),
  "median CD3CD4 events", step2$median_cd3_cd4_events[1],
    CP8_MANIFEST$benchmarks$median_primary_events,
    near(step2$median_cd3_cd4_events[1], CP8_MANIFEST$benchmarks$median_primary_events),
  "valid ages", metadata_summary$n_valid_age[1],
    CP8_MANIFEST$benchmarks$n_model,
    metadata_summary$n_valid_age[1] == CP8_MANIFEST$benchmarks$n_model,
  "main feature count", nrow(main_statistics),
    CP8_MANIFEST$benchmarks$main_feature_count,
    nrow(main_statistics) == CP8_MANIFEST$benchmarks$main_feature_count,
  "binary-sex feature count", nrow(binary_statistics),
    CP8_MANIFEST$benchmarks$main_feature_count,
    nrow(binary_statistics) == CP8_MANIFEST$benchmarks$main_feature_count,
  "subject-level scored rows", nrow(scored_data), 850, nrow(scored_data) == 850,
  "composite score count", nrow(score_statistics), 6, nrow(score_statistics) == 6,
  "integrated score beta", integrated_beta,
    CP8_MANIFEST$benchmarks$integrated_beta,
    near(integrated_beta, CP8_MANIFEST$benchmarks$integrated_beta),
  "age-stratified model count", nrow(age_stratified),
    3L * CP8_MANIFEST$benchmarks$age_outcome_count,
    nrow(age_stratified) == 3L * CP8_MANIFEST$benchmarks$age_outcome_count,
  "five-year matched pairs", match_summary$n_pairs[match_summary$caliper_years == 5][1],
    CP8_MANIFEST$benchmarks$matched_pairs_5y,
    match_summary$n_pairs[match_summary$caliper_years == 5][1] ==
      CP8_MANIFEST$benchmarks$matched_pairs_5y,
  "ten-year matched pairs", match_summary$n_pairs[match_summary$caliper_years == 10][1],
    CP8_MANIFEST$benchmarks$matched_pairs_10y,
    match_summary$n_pairs[match_summary$caliper_years == 10][1] ==
      CP8_MANIFEST$benchmarks$matched_pairs_10y,
  "five-year pair rows", nrow(pairs_5y),
    2L * CP8_MANIFEST$benchmarks$matched_pairs_5y,
    nrow(pairs_5y) == 2L * CP8_MANIFEST$benchmarks$matched_pairs_5y,
  "ten-year pair rows", nrow(pairs_10y),
    2L * CP8_MANIFEST$benchmarks$matched_pairs_10y,
    nrow(pairs_10y) == 2L * CP8_MANIFEST$benchmarks$matched_pairs_10y,
  "matched model count", nrow(matched_statistics),
    2L * CP8_MANIFEST$benchmarks$age_outcome_count,
    nrow(matched_statistics) == 2L * CP8_MANIFEST$benchmarks$age_outcome_count,
  "interaction model count", nrow(interaction_statistics),
    CP8_MANIFEST$benchmarks$age_outcome_count,
    nrow(interaction_statistics) == CP8_MANIFEST$benchmarks$age_outcome_count,
  "threshold model count", nrow(threshold_statistics),
    3L * CP8_MANIFEST$benchmarks$threshold_outcome_count,
    nrow(threshold_statistics) == 3L * CP8_MANIFEST$benchmarks$threshold_outcome_count,
  "threshold summary count", nrow(threshold_summary),
    CP8_MANIFEST$benchmarks$threshold_outcome_count,
    nrow(threshold_summary) == CP8_MANIFEST$benchmarks$threshold_outcome_count,
  "threshold scored rows", nrow(threshold_scored), 2550,
    nrow(threshold_scored) == 2550,
  "threshold direction robustness",
    sum(threshold_summary$direction_preserved, na.rm = TRUE),
    CP8_MANIFEST$benchmarks$threshold_outcome_count,
    sum(threshold_summary$direction_preserved, na.rm = TRUE) ==
      CP8_MANIFEST$benchmarks$threshold_outcome_count,
  "threshold FDR robustness",
    sum(threshold_summary$FDR_lt_0p05_all_sets, na.rm = TRUE),
    CP8_MANIFEST$benchmarks$threshold_outcome_count,
    sum(threshold_summary$FDR_lt_0p05_all_sets, na.rm = TRUE) ==
      CP8_MANIFEST$benchmarks$threshold_outcome_count
)

reference_dir <- path.expand(
  Sys.getenv("SDY2583_CP8_REFERENCE_DIR", unset = "")
)

compare_reference <- function(generated, reference_basename, keys) {
  all_reference_files <- list.files(
    reference_dir,
    full.names = TRUE,
    recursive = TRUE
  )
  hits <- all_reference_files[basename(all_reference_files) == reference_basename]
  if (length(hits) == 0L) {
    return(tibble::tibble(pass = FALSE, detail = "reference missing"))
  }
  rp_compare_csv(
    generated = generated,
    reference = hits[1],
    keys = keys,
    tolerance = tolerance
  )
}

if (nzchar(reference_dir) && dir.exists(reference_dir)) {
  comparisons <- list(
    list(
      generated = file.path(
        analysis_dir, "02_feature_extraction",
        "SDY2583_CP8_FULL_850_feature_table_STEP2_RECONSTRUCTED.csv"
      ),
      reference = "SDY2583_CP8_FULL_850_feature_table_STEP2.csv",
      keys = "subject_id"
    ),
    list(
      generated = file.path(
        analysis_dir, "06_statistics",
        "SDY2583_CP8_FULL850_age_sex_adjusted_statistics_MAIN_RECONSTRUCTED.csv"
      ),
      reference = "SDY2583_CP8_FULL850_age_sex_adjusted_statistics_MAIN.csv",
      keys = c("module", "feature")
    ),
    list(
      generated = file.path(
        analysis_dir, "07_sensitivity",
        "SDY2583_CP8_binary_sex_sensitivity_statistics_RECONSTRUCTED.csv"
      ),
      reference = "SDY2583_CP8_binary_sex_sensitivity_statistics.csv",
      keys = c("module", "feature")
    ),
    list(
      generated = file.path(
        analysis_dir, "08_composite_scores",
        "SDY2583_CP8_analysis_data_with_composite_scores_RECONSTRUCTED.csv"
      ),
      reference = "SDY2583_CP8_analysis_data_with_composite_scores.csv",
      keys = "subject_id"
    ),
    list(
      generated = file.path(
        analysis_dir, "08_composite_scores",
        "SDY2583_CP8_composite_score_statistics_RECONSTRUCTED.csv"
      ),
      reference = "SDY2583_CP8_composite_score_statistics.csv",
      keys = character()
    ),
    list(
      generated = file.path(
        analysis_dir, "10_age_sensitivity",
        "SDY2583_CP8_age_stratified_statistics_RECONSTRUCTED.csv"
      ),
      reference = "SDY2583_CP8_age_stratified_statistics.csv",
      keys = c("age_group_for_model", "feature")
    ),
    list(
      generated = file.path(
        analysis_dir, "10_age_sensitivity",
        "SDY2583_CP8_same_sex_age_matched_pairs_5yr_RECONSTRUCTED.csv"
      ),
      reference = "SDY2583_CP8_same_sex_age_matched_pairs_5yr.csv",
      keys = c("pair_id", "subject_id")
    ),
    list(
      generated = file.path(
        analysis_dir, "10_age_sensitivity",
        "SDY2583_CP8_same_sex_age_matched_pairs_10yr_RECONSTRUCTED.csv"
      ),
      reference = "SDY2583_CP8_same_sex_age_matched_pairs_10yr.csv",
      keys = c("pair_id", "subject_id")
    ),
    list(
      generated = file.path(
        analysis_dir, "10_age_sensitivity",
        "SDY2583_CP8_age_matched_sensitivity_statistics_RECONSTRUCTED.csv"
      ),
      reference = "SDY2583_CP8_age_matched_sensitivity_statistics.csv",
      keys = c("caliper_years", "feature")
    ),
    list(
      generated = file.path(
        analysis_dir, "10_age_sensitivity",
        "SDY2583_CP8_disease_by_age_interaction_statistics_RECONSTRUCTED.csv"
      ),
      reference = "SDY2583_CP8_disease_by_age_interaction_statistics.csv",
      keys = "feature"
    ),
    list(
      generated = file.path(
        analysis_dir, "12_threshold_sensitivity",
        "SDY2583_CP8_threshold_sensitivity_scored_data_LONG_RECONSTRUCTED.csv"
      ),
      reference = "SDY2583_CP8_threshold_sensitivity_scored_data_LONG.csv",
      keys = c("threshold_set", "subject_id")
    ),
    list(
      generated = file.path(
        analysis_dir, "12_threshold_sensitivity",
        "SDY2583_CP8_threshold_sensitivity_statistics_LONG_RECONSTRUCTED.csv"
      ),
      reference = "SDY2583_CP8_threshold_sensitivity_statistics_LONG.csv",
      keys = c("threshold_set", "feature")
    ),
    list(
      generated = file.path(
        analysis_dir, "12_threshold_sensitivity",
        "SDY2583_CP8_threshold_sensitivity_summary_WIDE_RECONSTRUCTED.csv"
      ),
      reference = "SDY2583_CP8_threshold_sensitivity_summary_WIDE.csv",
      keys = "feature"
    )
  )

  for (comparison in comparisons) {
    result <- compare_reference(
      comparison$generated,
      comparison$reference,
      comparison$keys
    )
    checks <- dplyr::bind_rows(
      checks,
      tibble::tibble(
        check = comparison$reference,
        observed = result$detail[1],
        expected = "table match",
        pass = result$pass[1]
      )
    )
  }
}

readr::write_csv(
  checks,
  file.path(out_dir, "SDY2583_CP8_reconstruction_validation_report.csv")
)
print(checks, n = nrow(checks))

if (strict && any(!checks$pass)) {
  stop("CP8 validation failed.")
}
