# ============================================================
# Validate the reconstructed CP7 pipeline against archived benchmarks.
#
# Usage from the repository root:
#   source("config/paths.R")            # optional local path configuration
#   source("scripts/20_validate_cp7_reconstruction.R")
#
# Optional deep validation:
#   Set SDY2583_CP7_REFERENCE_DIR to the local folder containing the archived
#   CP7 CSV outputs. The script then compares generated reconstructed tables
#   with their archived counterparts in addition to checking fixed benchmarks.
# ============================================================

source(file.path(Sys.getenv("SDY2583_REPO_ROOT", unset = "."), "R", "shared", "bootstrap.R"))

for (p in c("dplyr", "readr", "tibble", "purrr", "stringr")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(purrr)
  library(stringr)
})

analysis_dir <- sd_analysis_dir("CP7")
validation_dir <- file.path(analysis_dir, "validation")
dir.create(validation_dir, recursive = TRUE, showWarnings = FALSE)

reference_dir <- path.expand(Sys.getenv("SDY2583_CP7_REFERENCE_DIR", unset = ""))
strict <- tolower(Sys.getenv("SDY2583_VALIDATION_STRICT", unset = "true")) %in% c("true", "1", "yes", "y")
numeric_tolerance <- suppressWarnings(as.numeric(Sys.getenv("SDY2583_VALIDATION_TOLERANCE", unset = "1e-6")))
if (!is.finite(numeric_tolerance) || numeric_tolerance <= 0) numeric_tolerance <- 1e-6

report_rows <- list()
add_check <- function(check, observed, expected, pass, detail = "") {
  report_rows[[length(report_rows) + 1L]] <<- tibble(
    check = check,
    observed = as.character(observed),
    expected = as.character(expected),
    pass = isTRUE(pass),
    detail = detail
  )
}

near <- function(observed, expected, tolerance = numeric_tolerance) {
  if (length(observed) != 1L || length(expected) != 1L) return(FALSE)
  if (is.na(observed) || is.na(expected)) return(FALSE)
  abs(observed - expected) <= tolerance * max(1, abs(expected))
}

read_required <- function(path) {
  if (!file.exists(path)) stop("Required reconstructed output is missing: ", path)
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
}

# ------------------------------------------------------------
# 1. Fixed, non-sensitive archive benchmarks
# ------------------------------------------------------------

step1_summary_path <- file.path(analysis_dir, "01_channel_marker_QC", "SDY2583_CP7_STEP1_QC_summary.csv")
step2_summary_path <- file.path(analysis_dir, "02_feature_extraction", "SDY2583_CP7_FULL_850_STEP2_QC_summary_RECONSTRUCTED.csv")
metadata_summary_path <- file.path(analysis_dir, "04_metadata_merge", "SDY2583_CP7_metadata_merge_summary_RECONSTRUCTED.csv")
main_statistics_path <- file.path(analysis_dir, "06_statistics", "SDY2583_CP7_FULL850_age_sex_adjusted_statistics_MAIN_RECONSTRUCTED.csv")
composite_statistics_path <- file.path(analysis_dir, "08_composite_scores", "SDY2583_CP7_composite_score_statistics_RECONSTRUCTED.csv")

step1 <- read_required(step1_summary_path)
step2 <- read_required(step2_summary_path)
metadata <- read_required(metadata_summary_path)
main_stats <- read_required(main_statistics_path)
composite_stats <- read_required(composite_statistics_path)

add_check("CP7 FCS files", step1$n_fcs_files[[1]], 850, step1$n_fcs_files[[1]] == 850)
add_check("CP7 unique subjects in Step 1", step1$n_unique_subjects[[1]], 850, step1$n_unique_subjects[[1]] == 850)
add_check("CP7 Step 2 successful files", step2$feature_ok[[1]], 850, step2$feature_ok[[1]] == 850)
add_check("CP7 Step 2 failed files", step2$feature_failed[[1]], 0, step2$feature_failed[[1]] == 0)
add_check("Median total events", step2$median_total_events[[1]], 301989.5, near(step2$median_total_events[[1]], 301989.5))
add_check("Median CD3+CD8+ events", step2$median_cd3cd8_events[[1]], 12411.5, near(step2$median_cd3cd8_events[[1]], 12411.5))
add_check("Minimum CD3+CD8+ events", step2$min_cd3cd8_events[[1]], 407, step2$min_cd3cd8_events[[1]] == 407)
add_check("Maximum CD3+CD8+ events", step2$max_cd3cd8_events[[1]], 59837, step2$max_cd3cd8_events[[1]] == 59837)
add_check("Mapped disease groups", metadata$n_disease_mapped[[1]], 850, metadata$n_disease_mapped[[1]] == 850)
add_check("Valid ages", metadata$n_valid_age[[1]], 832, metadata$n_valid_age[[1]] == 832)
add_check("Overall healthy controls", metadata$n_healthy[[1]], 408, metadata$n_healthy[[1]] == 408)
add_check("Overall cancer patients", metadata$n_cancer[[1]], 442, metadata$n_cancer[[1]] == 442)
add_check("Main-model feature count", nrow(main_stats), 55, nrow(main_stats) == 55)

anchor_expectations <- tribble(
  ~feature, ~column, ~expected,
  "pct_cd3_pos_total", "beta_cancer_vs_healthy", -4.71800773276706,
  "pct_cd3_pos_total", "fdr_all", 4.54270592950129e-17,
  "pct_cd8_within_cd3", "beta_cancer_vs_healthy", 4.28759233054575,
  "pct_temra_like", "beta_cancer_vs_healthy", 4.36723790655398,
  "pct_cd39_pos", "beta_cancer_vs_healthy", 3.23434446149965
)

for (i in seq_len(nrow(anchor_expectations))) {
  feature <- anchor_expectations$feature[[i]]
  column <- anchor_expectations$column[[i]]
  expected <- anchor_expectations$expected[[i]]
  row <- main_stats %>% filter(.data$feature == feature)
  observed <- if (nrow(row) == 1L && column %in% names(row)) row[[column]][[1]] else NA_real_
  add_check(
    paste(feature, column), observed, expected, near(observed, expected),
    "Archived age- and sex-adjusted CP7 anchor"
  )
}

integrated <- composite_stats %>%
  filter(score == "CP7_integrated_checkpoint_remodeling_score")
integrated_beta <- if (nrow(integrated) == 1L) integrated$beta_cancer_vs_healthy[[1]] else NA_real_
add_check(
  "Integrated CP7 composite beta",
  integrated_beta,
  0.351119363746345,
  near(integrated_beta, 0.351119363746345),
  "Archived age- and sex-adjusted composite anchor"
)

# ------------------------------------------------------------
# 2. Optional table-level comparison with the local archive
# ------------------------------------------------------------

find_reference <- function(expected_basename) {
  if (!nzchar(reference_dir) || !dir.exists(reference_dir)) return(NA_character_)
  hits <- list.files(reference_dir, pattern = paste0("^", stringr::fixed(expected_basename), "$"),
                     full.names = TRUE, recursive = TRUE)
  if (length(hits) == 0L) NA_character_ else hits[[1L]]
}

compare_tables <- function(generated_path, reference_basename, key_columns) {
  reference_path <- find_reference(reference_basename)
  label <- paste0("Table comparison: ", reference_basename)
  if (is.na(reference_path)) {
    add_check(label, "not found", "reference file present", FALSE,
              "Set SDY2583_CP7_REFERENCE_DIR to the archived CP7 output root.")
    return(invisible(NULL))
  }

  generated <- read_required(generated_path)
  reference <- readr::read_csv(reference_path, show_col_types = FALSE, progress = FALSE)
  missing_keys <- setdiff(key_columns, intersect(names(generated), names(reference)))
  if (length(missing_keys) > 0L) {
    add_check(label, paste(missing_keys, collapse = ", "), "all key columns present", FALSE,
              "Key columns missing from one or both tables")
    return(invisible(NULL))
  }

  generated <- generated %>% arrange(across(all_of(key_columns)))
  reference <- reference %>% arrange(across(all_of(key_columns)))
  common_columns <- intersect(names(reference), names(generated))

  if (nrow(generated) != nrow(reference)) {
    add_check(label, nrow(generated), nrow(reference), FALSE, "Row-count mismatch")
    return(invisible(NULL))
  }

  failures <- character()
  for (column in common_columns) {
    g <- generated[[column]]
    r <- reference[[column]]
    if (is.numeric(g) && is.numeric(r)) {
      difference <- abs(g - r)
      limit <- numeric_tolerance * pmax(1, abs(r))
      bad <- which(!(is.na(g) & is.na(r)) & (is.na(g) != is.na(r) | (!is.na(difference) & difference > limit)))
      if (length(bad) > 0L) failures <- c(failures, paste0(column, "[", length(bad), "]"))
    } else {
      g_chr <- ifelse(is.na(g), "<NA>", as.character(g))
      r_chr <- ifelse(is.na(r), "<NA>", as.character(r))
      bad <- which(g_chr != r_chr)
      if (length(bad) > 0L) failures <- c(failures, paste0(column, "[", length(bad), "]"))
    }
  }

  add_check(
    label,
    ifelse(length(failures) == 0L, "matched", paste(failures, collapse = "; ")),
    "matched within tolerance",
    length(failures) == 0L,
    paste("Generated:", generated_path, "Reference:", reference_path)
  )
}

if (nzchar(reference_dir) && dir.exists(reference_dir)) {
  compare_tables(
    file.path(analysis_dir, "02_feature_extraction", "SDY2583_CP7_FULL_850_feature_table_STEP2_RECONSTRUCTED.csv"),
    "SDY2583_CP7_FULL_850_feature_table_STEP2.csv",
    c("subject_id")
  )
  compare_tables(
    main_statistics_path,
    "SDY2583_CP7_FULL850_age_sex_adjusted_statistics_MAIN.csv",
    c("module", "feature")
  )
  compare_tables(
    file.path(analysis_dir, "07_sensitivity", "SDY2583_CP7_binary_sex_only_sensitivity_statistics_RECONSTRUCTED.csv"),
    "SDY2583_CP7_binary_sex_only_sensitivity_statistics.csv",
    c("module", "feature")
  )
  compare_tables(
    file.path(analysis_dir, "07_sensitivity", "SDY2583_CP7_event_threshold_sensitivity_statistics_RECONSTRUCTED.csv"),
    "SDY2583_CP7_event_threshold_sensitivity_statistics.csv",
    c("subset_label", "module", "feature")
  )
  compare_tables(
    composite_statistics_path,
    "SDY2583_CP7_composite_score_statistics.csv",
    c("score")
  )
}

validation_report <- bind_rows(report_rows)
readr::write_csv(validation_report, file.path(validation_dir, "SDY2583_CP7_reconstruction_validation_report.csv"))

print(validation_report, n = nrow(validation_report))
failed <- validation_report %>% filter(!pass)
cat("\nPassed:", sum(validation_report$pass), "of", nrow(validation_report), "checks.\n")

if (nrow(failed) > 0L && strict) {
  stop("CP7 reconstruction validation failed. Review the validation report under outputs/CP7/validation/.")
}
