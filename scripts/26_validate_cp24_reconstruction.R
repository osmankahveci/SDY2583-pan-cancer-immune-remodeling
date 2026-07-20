# Validate reconstructed CP24 full-cohort outputs against fixed official
# clean-age benchmarks and, when configured, archived reference tables.

source(file.path(
  Sys.getenv("SDY2583_REPO_ROOT", unset = "."),
  "R", "shared", "reconstructed_panel_framework.R"
))
rp_install_and_load(c("dplyr", "readr", "tibble", "purrr", "stringr"))
source(file.path(
  sd_repo_root(), "R", "panels", "CP24", "MANIFEST_RECONSTRUCTED.R"
))

analysis_dir <- sd_analysis_dir("CP24")
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
  if (!file.exists(path)) stop("Missing CP24 output: ", path)
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
}

compare_archive_csv <- function(generated, reference, keys, tol = tolerance) {
  if (!file.exists(generated) || !file.exists(reference)) {
    return(tibble::tibble(pass = FALSE, detail = "missing file"))
  }

  generated_data <- readr::read_csv(
    generated,
    show_col_types = FALSE,
    progress = FALSE
  )
  reference_data <- readr::read_csv(
    reference,
    show_col_types = FALSE,
    progress = FALSE
  )

  ignored <- c(
    "file_path", "selected_fcs_dir", "analysis_dir", "out_dir", "rdata_dir",
    "metadata_source", "source_file", "source_path", "reference_file"
  )
  common <- setdiff(
    intersect(names(generated_data), names(reference_data)),
    ignored
  )
  keys <- intersect(keys, common)

  if (length(keys) > 0L) {
    generated_data <- generated_data |>
      dplyr::arrange(dplyr::across(dplyr::all_of(keys)))
    reference_data <- reference_data |>
      dplyr::arrange(dplyr::across(dplyr::all_of(keys)))
  }

  if (nrow(generated_data) != nrow(reference_data)) {
    return(tibble::tibble(
      pass = FALSE,
      detail = paste(
        "row mismatch",
        nrow(generated_data),
        nrow(reference_data)
      )
    ))
  }

  failures <- character()
  for (column in common) {
    generated_value <- generated_data[[column]]
    reference_value <- reference_data[[column]]

    if (is.numeric(generated_value) && is.numeric(reference_value)) {
      bad <- !(is.na(generated_value) & is.na(reference_value)) &
        (
          is.na(generated_value) != is.na(reference_value) |
            abs(generated_value - reference_value) >
              tol * pmax(1, abs(reference_value))
        )
    } else {
      bad <- ifelse(
        is.na(generated_value),
        "<NA>",
        as.character(generated_value)
      ) != ifelse(
        is.na(reference_value),
        "<NA>",
        as.character(reference_value)
      )
    }

    if (any(bad, na.rm = TRUE)) {
      failures <- c(
        failures,
        paste0(column, "[", sum(bad, na.rm = TRUE), "]")
      )
    }
  }

  tibble::tibble(
    pass = length(failures) == 0L,
    detail = ifelse(
      length(failures) == 0L,
      "matched",
      paste(failures, collapse = "; ")
    )
  )
}

feature_summary <- read_required(file.path(
  analysis_dir,
  "02_feature_extraction",
  "SDY2583_CP24_feature_extraction_summary_RECONSTRUCTED.csv"
))
metadata_summary <- read_required(file.path(
  analysis_dir,
  "04_metadata_merge_age_QC",
  "SDY2583_CP24_metadata_age_QC_summary_RECONSTRUCTED.csv"
))
main_statistics <- read_required(file.path(
  analysis_dir,
  "09_statistics_RECONSTRUCTED",
  "SDY2583_CP24_CLEANAGE_model_stats_RECONSTRUCTED.csv"
))
event_statistics <- read_required(file.path(
  analysis_dir,
  "09_statistics_RECONSTRUCTED",
  "SDY2583_CP24_event_count_sensitivity_all_outcomes_RECONSTRUCTED.csv"
))
direction_summary <- read_required(file.path(
  analysis_dir,
  "09_statistics_RECONSTRUCTED",
  "SDY2583_CP24_event_count_sensitivity_direction_summary_RECONSTRUCTED.csv"
))
score_data <- read_required(file.path(
  analysis_dir,
  "08_composite_scores_RECONSTRUCTED",
  "SDY2583_CP24_score_dataset_RECONSTRUCTED.csv"
))
score_statistics <- read_required(file.path(
  analysis_dir,
  "08_composite_scores_RECONSTRUCTED",
  "SDY2583_CP24_composite_score_results_CLEANAGE_RECONSTRUCTED.csv"
))
event_score_statistics <- read_required(file.path(
  analysis_dir,
  "08_composite_scores_RECONSTRUCTED",
  "SDY2583_CP24_composite_event_sensitivity_RECONSTRUCTED.csv"
))

integrated_beta <- score_statistics$beta_cancer[
  score_statistics$score == CP24_MANIFEST$integrated_score
][1]
core_beta <- score_statistics$beta_cancer[
  score_statistics$score == "CP24_core_differentiation_score"
][1]
pd1_score_fdr <- score_statistics$p_FDR_scores[
  score_statistics$score == "CP24_PD1_associated_terminal_score"
][1]

model_n_by_set <- event_statistics |>
  dplyr::group_by(analysis_set) |>
  dplyr::summarise(
    minimum = min(n_model, na.rm = TRUE),
    maximum = max(n_model, na.rm = TRUE),
    .groups = "drop"
  )
score_model_n_by_set <- event_score_statistics |>
  dplyr::group_by(analysis_set) |>
  dplyr::summarise(
    minimum = min(n_model, na.rm = TRUE),
    maximum = max(n_model, na.rm = TRUE),
    .groups = "drop"
  )

get_set_n <- function(table, label) {
  value <- table$minimum[table$analysis_set == label]
  if (length(value) == 0L) return(NA_real_)
  value[1]
}

checks <- tibble::tribble(
  ~check, ~observed, ~expected, ~pass,
  "FCS files", as.character(feature_summary$n_files[1]), "850", feature_summary$n_files[1] == 850,
  "unique subjects", as.character(feature_summary$n_subjects[1]), "850", feature_summary$n_subjects[1] == 850,
  "successful feature extractions", as.character(feature_summary$n_success[1]), "850", feature_summary$n_success[1] == 850,
  "median total events", as.character(feature_summary$median_total_events[1]), "302822", near(feature_summary$median_total_events[1], 302822),
  "median CD3+CD8+ events", as.character(feature_summary$median_cd3_cd8_events[1]), "11836", near(feature_summary$median_cd3_cd8_events[1], 11836),
  "subjects with >=500 CD3+CD8+ events", as.character(feature_summary$n_ge500[1]), "849", feature_summary$n_ge500[1] == 849,
  "subjects with >=1000 CD3+CD8+ events", as.character(feature_summary$n_ge1000[1]), "844", feature_summary$n_ge1000[1] == 844,
  "metadata rows", as.character(metadata_summary$n_rows[1]), "850", metadata_summary$n_rows[1] == 850,
  "valid clean ages", as.character(metadata_summary$n_valid_age[1]), "832", metadata_summary$n_valid_age[1] == 832,
  "complete age and sex", as.character(metadata_summary$n_complete_age_sex[1]), "828", metadata_summary$n_complete_age_sex[1] == 828,
  "healthy subjects", as.character(metadata_summary$n_healthy[1]), "408", metadata_summary$n_healthy[1] == 408,
  "cancer subjects", as.character(metadata_summary$n_cancer[1]), "442", metadata_summary$n_cancer[1] == 442,
  "official module rows", as.character(nrow(CP24_MANIFEST$module_map)), "55", nrow(CP24_MANIFEST$module_map) == 55,
  "official unique outcomes", as.character(dplyr::n_distinct(CP24_MANIFEST$module_map$feature)), "54", dplyr::n_distinct(CP24_MANIFEST$module_map$feature) == 54,
  "main statistic rows", as.character(nrow(main_statistics)), "55", nrow(main_statistics) == 55,
  "main globally significant rows", as.character(sum(main_statistics$p_FDR_all < 0.05, na.rm = TRUE)), "16", sum(main_statistics$p_FDR_all < 0.05, na.rm = TRUE) == 16,
  "main module-significant rows", as.character(sum(main_statistics$p_FDR_module < 0.05, na.rm = TRUE)), "17", sum(main_statistics$p_FDR_module < 0.05, na.rm = TRUE) == 17,
  "event sensitivity rows", as.character(nrow(event_statistics)), "165", nrow(event_statistics) == 165,
  "all-event reported model N", as.character(get_set_n(model_n_by_set, "all_850_clean_age")), "832", get_set_n(model_n_by_set, "all_850_clean_age") == 832,
  ">=500 reported model N", as.character(get_set_n(model_n_by_set, "cd3cd8_ge500_clean_age")), "831", get_set_n(model_n_by_set, "cd3cd8_ge500_clean_age") == 831,
  ">=1000 reported model N", as.character(get_set_n(model_n_by_set, "cd3cd8_ge1000_clean_age")), "826", get_set_n(model_n_by_set, "cd3cd8_ge1000_clean_age") == 826,
  "event directions preserved", as.character(sum(direction_summary$direction_preserved, na.rm = TRUE)), "55", sum(direction_summary$direction_preserved, na.rm = TRUE) == 55,
  "subject-level score rows", as.character(nrow(score_data)), "850", nrow(score_data) == 850,
  "official score columns present", as.character(sum(CP24_MANIFEST$score_names %in% names(score_data))), "7", all(CP24_MANIFEST$score_names %in% names(score_data)),
  "primary score models", as.character(nrow(score_statistics)), "7", nrow(score_statistics) == 7,
  "significant primary scores", as.character(sum(score_statistics$p_FDR_scores < 0.05, na.rm = TRUE)), "6", sum(score_statistics$p_FDR_scores < 0.05, na.rm = TRUE) == 6,
  "PD-1 terminal score nonsignificant", as.character(pd1_score_fdr), ">=0.05", is.finite(pd1_score_fdr) && pd1_score_fdr >= 0.05,
  "core differentiation beta", as.character(core_beta), "0.167915840013875", near(core_beta, 0.167915840013875),
  "integrated remodeling beta", as.character(integrated_beta), "0.113850267291264", near(integrated_beta, 0.113850267291264),
  "score sensitivity rows", as.character(nrow(event_score_statistics)), "21", nrow(event_score_statistics) == 21,
  "all-event score model N", as.character(get_set_n(score_model_n_by_set, "all_850_clean_age")), "832", get_set_n(score_model_n_by_set, "all_850_clean_age") == 832,
  ">=500 score model N", as.character(get_set_n(score_model_n_by_set, "cd3cd8_ge500_clean_age")), "831", get_set_n(score_model_n_by_set, "cd3cd8_ge500_clean_age") == 831,
  ">=1000 score model N", as.character(get_set_n(score_model_n_by_set, "cd3cd8_ge1000_clean_age")), "826", get_set_n(score_model_n_by_set, "cd3cd8_ge1000_clean_age") == 826,
  "quantitative figures", as.character(length(list.files(file.path(analysis_dir, "10_figures_RECONSTRUCTED"), "\\.pdf$", recursive = FALSE))), ">=4", length(list.files(file.path(analysis_dir, "10_figures_RECONSTRUCTED"), "\\.pdf$", recursive = FALSE)) >= 4
)

reference_dir <- path.expand(Sys.getenv("SDY2583_CP24_REFERENCE_DIR", unset = ""))
if (nzchar(reference_dir) && dir.exists(reference_dir)) {
  comparison_plan <- tibble::tribble(
    ~generated, ~reference_name, ~keys,
    file.path(analysis_dir, "04_metadata_merge_age_QC", "SDY2583_CP24_analysis_data_with_metadata_CLEANAGE_RECONSTRUCTED.csv"), "SDY2583_CP24_FULL_850_FINAL_feature_table.csv", "subject_id",
    file.path(analysis_dir, "09_statistics_RECONSTRUCTED", "SDY2583_CP24_CLEANAGE_model_stats_RECONSTRUCTED.csv"), "SDY2583_CP24_FULL_850_FINAL_all_outcome_statistics.csv", "module,outcome",
    file.path(analysis_dir, "09_statistics_RECONSTRUCTED", "SDY2583_CP24_event_count_sensitivity_all_outcomes_RECONSTRUCTED.csv"), "SDY2583_CP24_FULL_850_FINAL_sensitivity_all_outcomes.csv", "analysis_set,module,outcome",
    file.path(analysis_dir, "08_composite_scores_RECONSTRUCTED", "SDY2583_CP24_score_dataset_RECONSTRUCTED.csv"), "SDY2583_CP24_FULL_850_STEP5_score_dataset.csv", "subject_id",
    file.path(analysis_dir, "08_composite_scores_RECONSTRUCTED", "SDY2583_CP24_composite_score_results_CLEANAGE_RECONSTRUCTED.csv"), "SDY2583_CP24_FULL_850_FINAL_composite_score_results.csv", "score",
    file.path(analysis_dir, "08_composite_scores_RECONSTRUCTED", "SDY2583_CP24_composite_event_sensitivity_RECONSTRUCTED.csv"), "SDY2583_CP24_FULL_850_FINAL_score_sensitivity_results.csv", "analysis_set,score"
  )

  all_reference_csv <- list.files(
    reference_dir,
    "\\.csv$",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )

  archive_checks <- purrr::pmap_dfr(
    comparison_plan,
    function(generated, reference_name, keys) {
      hit <- all_reference_csv[basename(all_reference_csv) == reference_name]
      if (length(hit) == 0L) {
        return(tibble::tibble(
          check = paste0("archive table: ", reference_name),
          observed = "reference missing",
          expected = "matched",
          pass = FALSE
        ))
      }

      comparison <- compare_archive_csv(
        generated,
        hit[1],
        strsplit(keys, ",", fixed = TRUE)[[1]],
        tol = tolerance
      )
      tibble::tibble(
        check = paste0("archive table: ", reference_name),
        observed = comparison$detail[1],
        expected = "matched",
        pass = comparison$pass[1]
      )
    }
  )

  checks <- dplyr::bind_rows(checks, archive_checks)
}

report_file <- file.path(
  out_dir,
  "SDY2583_CP24_reconstruction_validation_report.csv"
)
readr::write_csv(checks, report_file)
print(checks, n = nrow(checks))

if (any(!checks$pass | is.na(checks$pass))) {
  stop("CP24 validation failed. Review ", report_file)
}
message("CP24 validation passed.")
