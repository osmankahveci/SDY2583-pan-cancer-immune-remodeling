# CP22-specific fixed-benchmark and optional archive-table validation.
source(file.path(Sys.getenv("SDY2583_REPO_ROOT", unset = "."), "R", "shared", "reconstructed_panel_framework.R"))
rp_install_and_load(c("dplyr", "readr", "tibble", "purrr", "stringr"))

a <- sd_analysis_dir("CP22")
out <- file.path(a, "validation")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
tol <- suppressWarnings(as.numeric(Sys.getenv("SDY2583_VALIDATION_TOLERANCE", unset = "1e-6")))
if (!is.finite(tol)) tol <- 1e-6
near <- function(x, y, tolerance = tol) is.finite(x) && abs(x - y) <= tolerance * max(1, abs(y))
read_required <- function(path) {
  if (!file.exists(path)) stop("Missing CP22 output: ", path)
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
}

inventory <- read_required(file.path(a, "01_fcs_inventory", "SDY2583_CP22_FCS_inventory_summary_STEP1.csv"))
step2 <- read_required(file.path(a, "02_feature_extraction", "SDY2583_CP22_feature_extraction_summary_STEP2_Bcell_humoral.csv"))
meta <- read_required(file.path(a, "03_metadata_merge_age_QC", "SDY2583_CP22_metadata_merge_summary_STEP3A.csv"))
model_summary <- read_required(file.path(a, "04_age_sex_adjusted_statistics", "SDY2583_CP22_overall_model_summary_STEP3B.csv"))
primary <- read_required(file.path(a, "04_age_sex_adjusted_statistics", "SDY2583_CP22_primary_age_sex_adjusted_results_STEP3B.csv"))
score_data <- read_required(file.path(a, "05_composite_scores", "SDY2583_CP22_analysis_data_with_composite_scores_STEP4.csv"))
primary_scores <- read_required(file.path(a, "05_composite_scores", "SDY2583_CP22_primary_composite_results_STEP4.csv"))
binary_scores <- read_required(file.path(a, "05_composite_scores", "SDY2583_CP22_binary_sex_composite_results_STEP4.csv"))
event_scores <- read_required(file.path(a, "05_composite_scores", "SDY2583_CP22_event_QC_composite_results_STEP4.csv"))
standard_scores <- read_required(file.path(a, "05_composite_scores", "SDY2583_CP22_standard_dump_composite_results_STEP4.csv"))
matched_summary <- read_required(file.path(a, "06_age_sensitivity_caliper_interaction", "SDY2583_CP22_matched_dataset_summary_STEP4B.csv"))
matched_5 <- read_required(file.path(a, "06_age_sensitivity_caliper_interaction", "SDY2583_CP22_same_sex_age_matched_5y_results_STEP4B.csv"))
matched_10 <- read_required(file.path(a, "06_age_sensitivity_caliper_interaction", "SDY2583_CP22_same_sex_age_matched_10y_results_STEP4B.csv"))
interactions <- read_required(file.path(a, "06_age_sensitivity_caliper_interaction", "SDY2583_CP22_disease_by_age_interaction_results_STEP4B.csv"))
threshold_scores <- read_required(file.path(a, "07_threshold_sensitivity", "SDY2583_CP22_threshold_analysis_data_with_scores_STEP5.csv"))
threshold_models <- read_required(file.path(a, "07_threshold_sensitivity", "SDY2583_CP22_threshold_model_results_STEP5.csv"))
threshold_robustness <- read_required(file.path(a, "07_threshold_sensitivity", "SDY2583_CP22_threshold_robustness_summary_STEP5.csv"))

integrated_name <- "CP22_integrated_humoral_B_cell_remodeling_score"
integrated_beta <- primary_scores$beta_cancer_vs_healthy[primary_scores$composite_score == integrated_name][1]
robust_n <- sum(threshold_robustness$robustness_class == "direction_and_global_FDR_preserved_all_sets", na.rm = TRUE)
module_only_n <- sum(threshold_robustness$robustness_class == "direction_and_module_FDR_preserved_all_sets", na.rm = TRUE)
direction_only_n <- sum(threshold_robustness$robustness_class == "direction_preserved_FDR_not_all_sets", na.rm = TRUE)
direction_fail_n <- sum(threshold_robustness$robustness_class == "direction_not_preserved", na.rm = TRUE)
all_composite_rows <- nrow(primary_scores) + nrow(binary_scores) + nrow(event_scores) + nrow(standard_scores)

checks <- tibble::tribble(
  ~check, ~observed, ~expected, ~pass,
  "FCS files", as.character(inventory$n_fcs_files[1]), "850", inventory$n_fcs_files[1] == 850,
  "unique subjects", as.character(inventory$n_unique_subjects[1]), "850", inventory$n_unique_subjects[1] == 850,
  "read-success files", as.character(inventory$n_read_ok[1]), "850", inventory$n_read_ok[1] == 850,
  "files with spillover", as.character(inventory$n_with_spillover_keyword[1]), "850", inventory$n_with_spillover_keyword[1] == 850,
  "median total events", as.character(inventory$median_events[1]), "303227.5", near(inventory$median_events[1], 303227.5),
  "channel mismatch files", as.character(inventory$n_channel_mismatch[1]), "2", inventory$n_channel_mismatch[1] == 2,
  "marker mismatch files", as.character(inventory$n_marker_mismatch[1]), "9", inventory$n_marker_mismatch[1] == 9,
  "successful feature extractions", as.character(step2$n_feature_ok[1]), "850", step2$n_feature_ok[1] == 850,
  "compensation applied", as.character(step2$n_compensation_applied[1]), "850", step2$n_compensation_applied[1] == 850,
  "transformation applied", as.character(step2$n_transform_applied[1]), "850", step2$n_transform_applied[1] == 850,
  "median CD19 B-cell events", as.character(step2$median_cd19_b_events[1]), "2666.5", near(step2$median_cd19_b_events[1], 2666.5),
  "minimum CD19 B-cell events", as.character(step2$min_cd19_b_events[1]), "6", step2$min_cd19_b_events[1] == 6,
  "model-ready subjects", as.character(meta$n_model_ready_age_sex[1]), "832", meta$n_model_ready_age_sex[1] == 832,
  "primary modeled features", as.character(nrow(primary)), "48", nrow(primary) == 48,
  "primary globally significant features", as.character(model_summary$n_primary_global_FDR_lt_0p05[1]), "35", model_summary$n_primary_global_FDR_lt_0p05[1] == 35,
  "feature directions preserved across sensitivities", as.character(model_summary$n_direction_preserved_all_sensitivities[1]), "44", model_summary$n_direction_preserved_all_sensitivities[1] == 44,
  "feature FDR preserved across sensitivities", as.character(model_summary$n_global_FDR_preserved_all_sensitivities[1]), "31", model_summary$n_global_FDR_preserved_all_sensitivities[1] == 31,
  "subject-level composite rows", as.character(nrow(score_data)), "850", nrow(score_data) == 850,
  "primary composite models", as.character(nrow(primary_scores)), "7", nrow(primary_scores) == 7,
  "all composite model rows", as.character(all_composite_rows), "28", all_composite_rows == 28,
  "integrated composite beta", as.character(integrated_beta), "0.4072877744339053", near(integrated_beta, 0.4072877744339053),
  "5-year matched pairs", as.character(matched_summary$n_pairs[matched_summary$caliper_years == 5][1]), "265", matched_summary$n_pairs[matched_summary$caliper_years == 5][1] == 265,
  "10-year matched pairs", as.character(matched_summary$n_pairs[matched_summary$caliper_years == 10][1]), "273", matched_summary$n_pairs[matched_summary$caliper_years == 10][1] == 273,
  "5-year matched composite models", as.character(nrow(matched_5)), "7", nrow(matched_5) == 7,
  "10-year matched composite models", as.character(nrow(matched_10)), "7", nrow(matched_10) == 7,
  "disease-by-age interaction models", as.character(nrow(interactions)), "7", nrow(interactions) == 7,
  "threshold scored rows", as.character(nrow(threshold_scores)), "2550", nrow(threshold_scores) == 2550,
  "threshold model rows", as.character(nrow(threshold_models)), "165", nrow(threshold_models) == 165,
  "threshold variables", as.character(dplyr::n_distinct(threshold_models$variable)), "55", dplyr::n_distinct(threshold_models$variable) == 55,
  "threshold globally robust variables", as.character(robust_n), "39", robust_n == 39,
  "threshold module-only robust variables", as.character(module_only_n), "0", module_only_n == 0,
  "threshold direction-only variables", as.character(direction_only_n), "6", direction_only_n == 6,
  "threshold direction failures", as.character(direction_fail_n), "10", direction_fail_n == 10
)

reference_dir <- path.expand(Sys.getenv("SDY2583_CP22_REFERENCE_DIR", unset = ""))
if (nzchar(reference_dir) && dir.exists(reference_dir)) {
  comparison_plan <- tibble::tribble(
    ~generated, ~reference_name, ~keys,
    file.path(a, "02_feature_extraction", "SDY2583_CP22_FULL_850_feature_table_STEP2_Bcell_humoral.csv"), "SDY2583_CP22_FULL_850_feature_table_STEP2_Bcell_humoral.csv", "subject_id",
    file.path(a, "03_metadata_merge_age_QC", "SDY2583_CP22_analysis_data_STEP3A_metadata_age_QC.csv"), "SDY2583_CP22_analysis_data_STEP3A_metadata_age_QC.csv", "subject_id",
    file.path(a, "04_age_sex_adjusted_statistics", "SDY2583_CP22_primary_age_sex_adjusted_results_STEP3B.csv"), "SDY2583_CP22_primary_age_sex_adjusted_results_STEP3B.csv", "feature",
    file.path(a, "04_age_sex_adjusted_statistics", "SDY2583_CP22_binary_sex_sensitivity_results_STEP3B.csv"), "SDY2583_CP22_binary_sex_sensitivity_results_STEP3B.csv", "feature",
    file.path(a, "04_age_sex_adjusted_statistics", "SDY2583_CP22_event_QC_sensitivity_results_STEP3B.csv"), "SDY2583_CP22_event_QC_sensitivity_results_STEP3B.csv", "feature",
    file.path(a, "04_age_sex_adjusted_statistics", "SDY2583_CP22_standard_dump_sensitivity_results_STEP3B.csv"), "SDY2583_CP22_standard_dump_sensitivity_results_STEP3B.csv", "feature",
    file.path(a, "05_composite_scores", "SDY2583_CP22_analysis_data_with_composite_scores_STEP4.csv"), "SDY2583_CP22_analysis_data_with_composite_scores_STEP4.csv", "subject_id",
    file.path(a, "05_composite_scores", "SDY2583_CP22_primary_composite_results_STEP4.csv"), "SDY2583_CP22_primary_composite_results_STEP4.csv", "composite_score",
    file.path(a, "05_composite_scores", "SDY2583_CP22_binary_sex_composite_results_STEP4.csv"), "SDY2583_CP22_binary_sex_composite_results_STEP4.csv", "composite_score",
    file.path(a, "05_composite_scores", "SDY2583_CP22_event_QC_composite_results_STEP4.csv"), "SDY2583_CP22_event_QC_composite_results_STEP4.csv", "composite_score",
    file.path(a, "05_composite_scores", "SDY2583_CP22_standard_dump_composite_results_STEP4.csv"), "SDY2583_CP22_standard_dump_composite_results_STEP4.csv", "composite_score",
    file.path(a, "06_age_sensitivity_caliper_interaction", "SDY2583_CP22_same_sex_age_matched_5y_results_STEP4B.csv"), "SDY2583_CP22_same_sex_age_matched_5y_results_STEP4B.csv", "composite_score",
    file.path(a, "06_age_sensitivity_caliper_interaction", "SDY2583_CP22_same_sex_age_matched_10y_results_STEP4B.csv"), "SDY2583_CP22_same_sex_age_matched_10y_results_STEP4B.csv", "composite_score",
    file.path(a, "06_age_sensitivity_caliper_interaction", "SDY2583_CP22_disease_by_age_interaction_results_STEP4B.csv"), "SDY2583_CP22_disease_by_age_interaction_results_STEP4B.csv", "composite_score",
    file.path(a, "07_threshold_sensitivity", "SDY2583_CP22_threshold_analysis_data_with_scores_STEP5.csv"), "SDY2583_CP22_threshold_analysis_data_with_scores_STEP5.csv", "threshold_set,subject_id",
    file.path(a, "07_threshold_sensitivity", "SDY2583_CP22_threshold_model_results_STEP5.csv"), "SDY2583_CP22_threshold_model_results_STEP5.csv", "variable,threshold_set",
    file.path(a, "07_threshold_sensitivity", "SDY2583_CP22_threshold_robustness_summary_STEP5.csv"), "SDY2583_CP22_threshold_robustness_summary_STEP5.csv", "variable"
  )
  all_reference_csv <- list.files(reference_dir, "\\.csv$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  archive_checks <- purrr::pmap_dfr(comparison_plan, function(generated, reference_name, keys) {
    hit <- all_reference_csv[basename(all_reference_csv) == reference_name]
    if (length(hit) == 0L) return(tibble::tibble(check = paste0("archive table: ", reference_name), observed = "reference missing", expected = "matched", pass = FALSE))
    cmp <- rp_compare_csv(generated, hit[1], strsplit(keys, ",", fixed = TRUE)[[1]], tolerance = tol)
    tibble::tibble(check = paste0("archive table: ", reference_name), observed = cmp$detail[1], expected = "matched", pass = cmp$pass[1])
  })
  checks <- dplyr::bind_rows(checks, archive_checks)
}

report_file <- file.path(out, "SDY2583_CP22_recovered_pipeline_validation_report.csv")
readr::write_csv(checks, report_file)
print(checks, n = nrow(checks))
if (any(!checks$pass | is.na(checks$pass))) stop("CP22 validation failed. Review ", report_file)
message("CP22 validation passed.")
