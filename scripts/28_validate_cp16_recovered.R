# CP16-specific fixed-benchmark and optional archive-table validation.
source(file.path(Sys.getenv("SDY2583_REPO_ROOT", unset = "."), "R", "shared", "reconstructed_panel_framework.R"))
rp_install_and_load(c("dplyr", "readr", "tibble", "purrr", "stringr"))

a <- sd_analysis_dir("CP16")
out <- file.path(a, "validation")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
tol <- suppressWarnings(as.numeric(Sys.getenv("SDY2583_VALIDATION_TOLERANCE", unset = "1e-6")))
if (!is.finite(tol)) tol <- 1e-6
near <- function(x, y, tolerance = tol) is.finite(x) && abs(x - y) <= tolerance * max(1, abs(y))
read_required <- function(path) {
  if (!file.exists(path)) stop("Missing CP16 output: ", path)
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
}

inventory <- read_required(file.path(a, "01_fcs_inventory_marker_QC", "SDY2583_CP16_fcs_inventory_summary_STEP1.csv"))
step2 <- read_required(file.path(a, "02_feature_extraction", "SDY2583_CP16_feature_extraction_summary_STEP2.csv"))
meta <- read_required(file.path(a, "03_metadata_merge_age_QC", "SDY2583_CP16_merge_summary_STEP3A.csv"))
model_summary <- read_required(file.path(a, "04_age_sex_adjusted_models", "SDY2583_CP16_model_subject_summary_STEP3B.csv"))
primary <- read_required(file.path(a, "04_age_sex_adjusted_models", "SDY2583_CP16_primary_age_sex_adjusted_results_STEP3B.csv"))
score_data <- read_required(file.path(a, "05_composite_scores", "SDY2583_CP16_scores_data_STEP4.csv"))
primary_scores <- read_required(file.path(a, "05_composite_scores", "SDY2583_CP16_primary_composite_results_STEP4.csv"))
all_scores <- read_required(file.path(a, "05_composite_scores", "SDY2583_CP16_all_composite_model_results_STEP4.csv"))
matched_summary <- read_required(file.path(a, "06_age_sensitivity_caliper_interaction", "SDY2583_CP16_matched_dataset_summary_STEP4B.csv"))
matched_5 <- read_required(file.path(a, "06_age_sensitivity_caliper_interaction", "SDY2583_CP16_same_sex_age_matched_5y_results_STEP4B.csv"))
matched_10 <- read_required(file.path(a, "06_age_sensitivity_caliper_interaction", "SDY2583_CP16_same_sex_age_matched_10y_results_STEP4B.csv"))
interactions <- read_required(file.path(a, "06_age_sensitivity_caliper_interaction", "SDY2583_CP16_disease_by_age_interaction_results_STEP4B.csv"))
threshold_scores <- read_required(file.path(a, "07_threshold_sensitivity", "SDY2583_CP16_threshold_scores_data_STEP5.csv"))
threshold_models <- read_required(file.path(a, "07_threshold_sensitivity", "SDY2583_CP16_threshold_model_results_STEP5.csv"))
threshold_robustness <- read_required(file.path(a, "07_threshold_sensitivity", "SDY2583_CP16_threshold_robustness_summary_STEP5.csv"))

integrated_name <- "CP16_integrated_APC_DC_myeloid_remodeling_score"
integrated_beta <- primary_scores$beta_cancer_vs_healthy[primary_scores$composite_score == integrated_name][1]
robust_n <- sum(threshold_robustness$robustness_class == "direction_and_global_FDR_preserved_all_sets", na.rm = TRUE)
direction_only_n <- sum(threshold_robustness$robustness_class == "direction_preserved_FDR_not_all_sets", na.rm = TRUE)
direction_fail_n <- sum(threshold_robustness$robustness_class == "direction_not_preserved", na.rm = TRUE)

checks <- tibble::tribble(
  ~check, ~observed, ~expected, ~pass,
  "FCS files", as.character(inventory$n_fcs_files[1]), "847", inventory$n_fcs_files[1] == 847,
  "unique subjects", as.character(inventory$n_unique_subjects[1]), "847", inventory$n_unique_subjects[1] == 847,
  "read-success files", as.character(inventory$n_read_ok[1]), "847", inventory$n_read_ok[1] == 847,
  "files with spillover", as.character(inventory$n_with_spillover_keyword[1]), "847", inventory$n_with_spillover_keyword[1] == 847,
  "modal channel count", as.character(inventory$modal_n_channels[1]), "19", inventory$modal_n_channels[1] == 19,
  "median total events", as.character(inventory$median_events[1]), "290935", near(inventory$median_events[1], 290935),
  "channel-order mismatch files", as.character(inventory$n_channel_order_mismatch_files[1]), "2", inventory$n_channel_order_mismatch_files[1] == 2,
  "marker-order mismatch files", as.character(inventory$n_marker_order_mismatch_files[1]), "2", inventory$n_marker_order_mismatch_files[1] == 2,
  "successful feature extractions", as.character(step2$n_feature_ok[1]), "847", step2$n_feature_ok[1] == 847,
  "median CD45 dump-low events", as.character(step2$median_cd45_dump_low_events[1]), "427", near(step2$median_cd45_dump_low_events[1], 427),
  "median HLA-DR APC-core events", as.character(step2$median_hladr_apc_core_events[1]), "357", near(step2$median_hladr_apc_core_events[1], 357),
  "minimum HLA-DR APC-core events", as.character(step2$min_hladr_apc_core_events[1]), "9", step2$min_hladr_apc_core_events[1] == 9,
  "model-ready subjects", as.character(meta$n_model_ready_age_sex[1]), "829", meta$n_model_ready_age_sex[1] == 829,
  "primary modeled features", as.character(nrow(primary)), "54", nrow(primary) == 54,
  "primary globally significant features", as.character(model_summary$n_primary_global_FDR_lt_0p05[1]), "40", model_summary$n_primary_global_FDR_lt_0p05[1] == 40,
  "feature directions preserved across sensitivities", as.character(model_summary$n_direction_preserved_all_sensitivities[1]), "51", model_summary$n_direction_preserved_all_sensitivities[1] == 51,
  "feature FDR preserved across sensitivities", as.character(model_summary$n_FDR_preserved_all_sensitivities[1]), "36", model_summary$n_FDR_preserved_all_sensitivities[1] == 36,
  "subject-level composite rows", as.character(nrow(score_data)), "847", nrow(score_data) == 847,
  "primary composite models", as.character(nrow(primary_scores)), "6", nrow(primary_scores) == 6,
  "all composite model rows", as.character(nrow(all_scores)), "30", nrow(all_scores) == 30,
  "integrated composite beta", as.character(integrated_beta), "0.447860154036832", near(integrated_beta, 0.447860154036832),
  "5-year matched pairs", as.character(matched_summary$n_pairs[matched_summary$caliper_years == 5][1]), "264", matched_summary$n_pairs[matched_summary$caliper_years == 5][1] == 264,
  "10-year matched pairs", as.character(matched_summary$n_pairs[matched_summary$caliper_years == 10][1]), "272", matched_summary$n_pairs[matched_summary$caliper_years == 10][1] == 272,
  "5-year matched composite models", as.character(nrow(matched_5)), "6", nrow(matched_5) == 6,
  "10-year matched composite models", as.character(nrow(matched_10)), "6", nrow(matched_10) == 6,
  "disease-by-age interaction models", as.character(nrow(interactions)), "6", nrow(interactions) == 6,
  "threshold scored rows", as.character(nrow(threshold_scores)), "2541", nrow(threshold_scores) == 2541,
  "threshold model rows", as.character(nrow(threshold_models)), "180", nrow(threshold_models) == 180,
  "threshold variables", as.character(dplyr::n_distinct(threshold_models$variable)), "60", dplyr::n_distinct(threshold_models$variable) == 60,
  "threshold robust variables", as.character(robust_n), "33", robust_n == 33,
  "threshold direction-only variables", as.character(direction_only_n), "14", direction_only_n == 14,
  "threshold direction failures", as.character(direction_fail_n), "13", direction_fail_n == 13
)

reference_dir <- path.expand(Sys.getenv("SDY2583_CP16_REFERENCE_DIR", unset = ""))
if (nzchar(reference_dir) && dir.exists(reference_dir)) {
  comparison_plan <- tibble::tribble(
    ~generated, ~reference_name, ~keys,
    file.path(a, "02_feature_extraction", "SDY2583_CP16_features_STEP2.csv"), "SDY2583_CP16_features_STEP2.csv", "subject_id",
    file.path(a, "03_metadata_merge_age_QC", "SDY2583_CP16_analysis_data_STEP3A.csv"), "SDY2583_CP16_analysis_data_STEP3A.csv", "subject_id",
    file.path(a, "04_age_sex_adjusted_models", "SDY2583_CP16_primary_age_sex_adjusted_results_STEP3B.csv"), "SDY2583_CP16_primary_age_sex_adjusted_results_STEP3B.csv", "feature",
    file.path(a, "04_age_sex_adjusted_models", "SDY2583_CP16_binary_sex_sensitivity_results_STEP3B.csv"), "SDY2583_CP16_binary_sex_sensitivity_results_STEP3B.csv", "feature",
    file.path(a, "04_age_sex_adjusted_models", "SDY2583_CP16_CD45_dump_low_event_QC_results_STEP3B.csv"), "SDY2583_CP16_CD45_dump_low_event_QC_results_STEP3B.csv", "feature",
    file.path(a, "04_age_sex_adjusted_models", "SDY2583_CP16_HLADR_APC_core_event_QC_results_STEP3B.csv"), "SDY2583_CP16_HLADR_APC_core_event_QC_results_STEP3B.csv", "feature",
    file.path(a, "04_age_sex_adjusted_models", "SDY2583_CP16_technical_QC_sensitivity_results_STEP3B.csv"), "SDY2583_CP16_technical_QC_sensitivity_results_STEP3B.csv", "feature",
    file.path(a, "05_composite_scores", "SDY2583_CP16_scores_data_STEP4.csv"), "SDY2583_CP16_scores_data_STEP4.csv", "subject_id",
    file.path(a, "05_composite_scores", "SDY2583_CP16_all_composite_model_results_STEP4.csv"), "SDY2583_CP16_all_composite_model_results_STEP4.csv", "model_label,composite_score",
    file.path(a, "06_age_sensitivity_caliper_interaction", "SDY2583_CP16_same_sex_age_matched_pairs_5y_STEP4B.csv"), "SDY2583_CP16_same_sex_age_matched_pairs_5y_STEP4B.csv", "pair_id",
    file.path(a, "06_age_sensitivity_caliper_interaction", "SDY2583_CP16_same_sex_age_matched_pairs_10y_STEP4B.csv"), "SDY2583_CP16_same_sex_age_matched_pairs_10y_STEP4B.csv", "pair_id",
    file.path(a, "06_age_sensitivity_caliper_interaction", "SDY2583_CP16_same_sex_age_matched_5y_results_STEP4B.csv"), "SDY2583_CP16_same_sex_age_matched_5y_results_STEP4B.csv", "composite_score",
    file.path(a, "06_age_sensitivity_caliper_interaction", "SDY2583_CP16_same_sex_age_matched_10y_results_STEP4B.csv"), "SDY2583_CP16_same_sex_age_matched_10y_results_STEP4B.csv", "composite_score",
    file.path(a, "06_age_sensitivity_caliper_interaction", "SDY2583_CP16_disease_by_age_interaction_results_STEP4B.csv"), "SDY2583_CP16_disease_by_age_interaction_results_STEP4B.csv", "composite_score",
    file.path(a, "07_threshold_sensitivity", "SDY2583_CP16_threshold_scores_data_STEP5.csv"), "SDY2583_CP16_threshold_scores_data_STEP5.csv", "threshold_set,subject_id",
    file.path(a, "07_threshold_sensitivity", "SDY2583_CP16_threshold_model_results_STEP5.csv"), "SDY2583_CP16_threshold_model_results_STEP5.csv", "variable,threshold_set",
    file.path(a, "07_threshold_sensitivity", "SDY2583_CP16_threshold_robustness_summary_STEP5.csv"), "SDY2583_CP16_threshold_robustness_summary_STEP5.csv", "variable"
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

report_file <- file.path(out, "SDY2583_CP16_recovered_pipeline_validation_report.csv")
readr::write_csv(checks, report_file)
print(checks, n = nrow(checks))
if (any(!checks$pass | is.na(checks$pass))) stop("CP16 validation failed. Review ", report_file)
message("CP16 validation passed.")
