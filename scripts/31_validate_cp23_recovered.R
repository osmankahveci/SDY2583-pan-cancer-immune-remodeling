# CP23-specific fixed-benchmark and optional archive-table validation.
source(file.path(Sys.getenv("SDY2583_REPO_ROOT", unset = "."), "R", "shared", "reconstructed_panel_framework.R"))
rp_install_and_load(c("dplyr", "readr", "tibble", "purrr", "stringr"))

a <- sd_analysis_dir("CP23")
out <- file.path(a, "validation")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

tol <- suppressWarnings(as.numeric(Sys.getenv("SDY2583_VALIDATION_TOLERANCE", unset = "1e-6")))
if (!is.finite(tol)) tol <- 1e-6
near <- function(x, y, tolerance = tol) is.finite(x) && abs(x - y) <= tolerance * max(1, abs(y))

read_required <- function(path) {
  if (!file.exists(path)) stop("Missing CP23 output: ", path)
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
}

# Environment-dependent absolute paths are intentionally excluded; all scientific
# and identifier columns shared by generated and archived tables remain compared.
compare_archive_csv <- function(generated, reference, keys, tolerance = tol) {
  if (!file.exists(generated) || !file.exists(reference)) {
    return(tibble::tibble(pass = FALSE, detail = "missing file"))
  }
  g <- readr::read_csv(generated, show_col_types = FALSE, progress = FALSE)
  r <- readr::read_csv(reference, show_col_types = FALSE, progress = FALSE)
  ignored <- c(
    "file_path", "selected_fcs_dir", "analysis_dir", "out_dir", "rdata_dir",
    "source_file", "source_path", "reference_file"
  )
  common <- setdiff(intersect(names(g), names(r)), ignored)
  keys <- intersect(keys, common)
  if (length(keys) > 0L) {
    g <- g |> dplyr::arrange(dplyr::across(dplyr::all_of(keys)))
    r <- r |> dplyr::arrange(dplyr::across(dplyr::all_of(keys)))
  }
  if (nrow(g) != nrow(r)) {
    return(tibble::tibble(pass = FALSE, detail = paste("row mismatch", nrow(g), nrow(r))))
  }
  failures <- character()
  for (col in common) {
    gv <- g[[col]]
    rv <- r[[col]]
    if (is.numeric(gv) && is.numeric(rv)) {
      bad <- !(is.na(gv) & is.na(rv)) &
        (is.na(gv) != is.na(rv) | abs(gv - rv) > tolerance * pmax(1, abs(rv)))
    } else {
      bad <- ifelse(is.na(gv), "<NA>", as.character(gv)) !=
        ifelse(is.na(rv), "<NA>", as.character(rv))
    }
    if (any(bad, na.rm = TRUE)) {
      failures <- c(failures, paste0(col, "[", sum(bad, na.rm = TRUE), "]"))
    }
  }
  tibble::tibble(
    pass = length(failures) == 0L,
    detail = ifelse(length(failures) == 0L, "matched", paste(failures, collapse = "; "))
  )
}

inventory <- read_required(file.path(
  a, "01_fcs_inventory_marker_QC", "SDY2583_CP23_fcs_inventory_summary_STEP1.csv"
))
step2 <- read_required(file.path(
  a, "02_feature_extraction", "SDY2583_CP23_feature_extraction_summary_STEP2.csv"
))
meta <- read_required(file.path(
  a, "03_metadata_merge_age_QC", "SDY2583_CP23_merge_summary_STEP3A.csv"
))
feature_models <- read_required(file.path(
  a, "04_age_sex_adjusted_models", "SDY2583_CP23_all_feature_model_results_STEP3B.csv"
))
feature_robustness <- read_required(file.path(
  a, "04_age_sex_adjusted_models", "SDY2583_CP23_sensitivity_preservation_summary_STEP3B.csv"
))
score_data <- read_required(file.path(
  a, "05_composite_scores", "SDY2583_CP23_score_data_STEP4.csv"
))
composite_models <- read_required(file.path(
  a, "05_composite_scores", "SDY2583_CP23_all_composite_model_results_STEP4.csv"
))
matched_5_data <- read_required(file.path(
  a, "06_age_sensitivity_caliper_interaction", "SDY2583_CP23_same_sex_age_matched_5y_dataset_STEP4B.csv"
))
matched_10_data <- read_required(file.path(
  a, "06_age_sensitivity_caliper_interaction", "SDY2583_CP23_same_sex_age_matched_10y_dataset_STEP4B.csv"
))
matched_5 <- read_required(file.path(
  a, "06_age_sensitivity_caliper_interaction", "SDY2583_CP23_same_sex_age_matched_5y_results_STEP4B.csv"
))
matched_10 <- read_required(file.path(
  a, "06_age_sensitivity_caliper_interaction", "SDY2583_CP23_same_sex_age_matched_10y_results_STEP4B.csv"
))
interactions <- read_required(file.path(
  a, "06_age_sensitivity_caliper_interaction", "SDY2583_CP23_disease_by_age_interaction_results_STEP4B.csv"
))
threshold_scores <- read_required(file.path(
  a, "07_threshold_sensitivity", "SDY2583_CP23_threshold_score_data_STEP5.csv"
))
threshold_models <- read_required(file.path(
  a, "07_threshold_sensitivity", "SDY2583_CP23_threshold_model_results_STEP5.csv"
))
threshold_robustness <- read_required(file.path(
  a, "07_threshold_sensitivity", "SDY2583_CP23_threshold_wide_robustness_STEP5.csv"
))

integrated_name <- "CP23_integrated_monocyte_macrophage_like_myeloid_remodeling_score"
primary_composites <- composite_models |>
  dplyr::filter(model_label == "primary_age_sex_adjusted")
integrated_beta <- primary_composites$beta_cancer_vs_healthy[
  primary_composites$composite_score == integrated_name
][1]

feature_robust_n <- sum(
  feature_robustness$robustness_class ==
    "direction_and_global_FDR_preserved_all_sensitivities",
  na.rm = TRUE
)
feature_direction_only_n <- sum(
  feature_robustness$robustness_class ==
    "direction_preserved_FDR_not_all_sensitivities",
  na.rm = TRUE
)
feature_direction_fail_n <- sum(
  feature_robustness$robustness_class == "direction_not_preserved",
  na.rm = TRUE
)

threshold_robust_n <- sum(
  threshold_robustness$robustness_class ==
    "direction_and_global_FDR_preserved_all_sets",
  na.rm = TRUE
)
threshold_direction_only_n <- sum(
  threshold_robustness$robustness_class ==
    "direction_preserved_FDR_not_all_sets",
  na.rm = TRUE
)
threshold_direction_fail_n <- sum(
  threshold_robustness$robustness_class == "direction_not_preserved",
  na.rm = TRUE
)

checks <- tibble::tribble(
  ~check, ~observed, ~expected, ~pass,
  "FCS files", as.character(inventory$n_fcs_files[1]), "850", inventory$n_fcs_files[1] == 850,
  "unique subjects", as.character(inventory$n_unique_subjects[1]), "850", inventory$n_unique_subjects[1] == 850,
  "read-success files", as.character(inventory$n_read_ok[1]), "850", inventory$n_read_ok[1] == 850,
  "files with spillover", as.character(inventory$n_with_spillover_keyword[1]), "850", inventory$n_with_spillover_keyword[1] == 850,
  "median total events", as.character(inventory$median_events[1]), "280931", near(inventory$median_events[1], 280931),
  "channel-order mismatch files", as.character(inventory$n_channel_order_mismatch_files[1]), "2", inventory$n_channel_order_mismatch_files[1] == 2,
  "marker-order mismatch files", as.character(inventory$n_marker_order_mismatch_files[1]), "2", inventory$n_marker_order_mismatch_files[1] == 2,
  "marker-set mismatch files", as.character(inventory$n_marker_set_mismatch_files[1]), "0", inventory$n_marker_set_mismatch_files[1] == 0,
  "successful feature extractions", as.character(step2$n_feature_ok[1]), "850", step2$n_feature_ok[1] == 850,
  "median CD45 dump-low events", as.character(step2$median_cd45_dump_low_events[1]), "15173.5", near(step2$median_cd45_dump_low_events[1], 15173.5),
  "minimum CD45 dump-low events", as.character(step2$min_cd45_dump_low_events[1]), "22", step2$min_cd45_dump_low_events[1] == 22,
  "median CD33+HLA-DR+ myeloid-like events", as.character(step2$median_cd33_hladr_myeloid_like_events[1]), "2818", near(step2$median_cd33_hladr_myeloid_like_events[1], 2818),
  "minimum CD33+HLA-DR+ myeloid-like events", as.character(step2$min_cd33_hladr_myeloid_like_events[1]), "4", step2$min_cd33_hladr_myeloid_like_events[1] == 4,
  "model-ready age/sex subjects", as.character(meta$n_model_ready_age_sex[1]), "832", meta$n_model_ready_age_sex[1] == 832,
  "binary-sex subjects", as.character(meta$n_with_binary_sex[1]), "835", meta$n_with_binary_sex[1] == 835,
  "feature model rows", as.character(nrow(feature_models)), "345", nrow(feature_models) == 345,
  "modeled features", as.character(dplyr::n_distinct(feature_models$feature)), "69", dplyr::n_distinct(feature_models$feature) == 69,
  "feature robust across all sensitivities", as.character(feature_robust_n), "15", feature_robust_n == 15,
  "feature direction-only robustness", as.character(feature_direction_only_n), "48", feature_direction_only_n == 48,
  "feature direction failures", as.character(feature_direction_fail_n), "6", feature_direction_fail_n == 6,
  "subject-level composite rows", as.character(nrow(score_data)), "850", nrow(score_data) == 850,
  "subjects with integrated composite", as.character(sum(!is.na(score_data[[integrated_name]]))), "850", sum(!is.na(score_data[[integrated_name]])) == 850,
  "composite model rows", as.character(nrow(composite_models)), "25", nrow(composite_models) == 25,
  "primary composite models", as.character(nrow(primary_composites)), "5", nrow(primary_composites) == 5,
  "integrated composite beta", as.character(integrated_beta), "0.2476159779193443", near(integrated_beta, 0.2476159779193443),
  "5-year matched rows", as.character(nrow(matched_5_data)), "530", nrow(matched_5_data) == 530,
  "5-year matched pairs", as.character(dplyr::n_distinct(matched_5_data$matched_pair_id)), "265", dplyr::n_distinct(matched_5_data$matched_pair_id) == 265,
  "10-year matched rows", as.character(nrow(matched_10_data)), "546", nrow(matched_10_data) == 546,
  "10-year matched pairs", as.character(dplyr::n_distinct(matched_10_data$matched_pair_id)), "273", dplyr::n_distinct(matched_10_data$matched_pair_id) == 273,
  "5-year matched composite models", as.character(nrow(matched_5)), "5", nrow(matched_5) == 5,
  "10-year matched composite models", as.character(nrow(matched_10)), "5", nrow(matched_10) == 5,
  "disease-by-age interaction models", as.character(nrow(interactions)), "5", nrow(interactions) == 5,
  "threshold scored rows", as.character(nrow(threshold_scores)), "2550", nrow(threshold_scores) == 2550,
  "threshold subjects per set", as.character(min(table(threshold_scores$threshold_set))), "850", all(table(threshold_scores$threshold_set) == 850),
  "threshold model rows", as.character(nrow(threshold_models)), "60", nrow(threshold_models) == 60,
  "threshold variables", as.character(dplyr::n_distinct(threshold_models$variable)), "20", dplyr::n_distinct(threshold_models$variable) == 20,
  "threshold composite variables", as.character(dplyr::n_distinct(threshold_models$variable[threshold_models$variable_type == "composite"])), "5", dplyr::n_distinct(threshold_models$variable[threshold_models$variable_type == "composite"]) == 5,
  "threshold feature variables", as.character(dplyr::n_distinct(threshold_models$variable[threshold_models$variable_type == "feature"])), "15", dplyr::n_distinct(threshold_models$variable[threshold_models$variable_type == "feature"]) == 15,
  "threshold globally robust variables", as.character(threshold_robust_n), "16", threshold_robust_n == 16,
  "threshold direction-only variables", as.character(threshold_direction_only_n), "4", threshold_direction_only_n == 4,
  "threshold direction failures", as.character(threshold_direction_fail_n), "0", threshold_direction_fail_n == 0
)

reference_dir <- path.expand(Sys.getenv("SDY2583_CP23_REFERENCE_DIR", unset = ""))
if (nzchar(reference_dir) && dir.exists(reference_dir)) {
  comparison_plan <- tibble::tribble(
    ~generated, ~reference_name, ~keys,
    file.path(a, "02_feature_extraction", "SDY2583_CP23_features_STEP2.csv"), "SDY2583_CP23_features_STEP2.csv", "subject_id",
    file.path(a, "03_metadata_merge_age_QC", "SDY2583_CP23_analysis_data_STEP3A.csv"), "SDY2583_CP23_analysis_data_STEP3A.csv", "subject_id",
    file.path(a, "04_age_sex_adjusted_models", "SDY2583_CP23_all_feature_model_results_STEP3B.csv"), "SDY2583_CP23_all_feature_model_results_STEP3B.csv", "model_label,feature",
    file.path(a, "04_age_sex_adjusted_models", "SDY2583_CP23_sensitivity_preservation_summary_STEP3B.csv"), "SDY2583_CP23_sensitivity_preservation_summary_STEP3B.csv", "feature",
    file.path(a, "05_composite_scores", "SDY2583_CP23_score_data_STEP4.csv"), "SDY2583_CP23_score_data_STEP4.csv", "subject_id",
    file.path(a, "05_composite_scores", "SDY2583_CP23_all_composite_model_results_STEP4.csv"), "SDY2583_CP23_all_composite_model_results_STEP4.csv", "model_label,composite_score",
    file.path(a, "06_age_sensitivity_caliper_interaction", "SDY2583_CP23_same_sex_age_matched_5y_dataset_STEP4B.csv"), "SDY2583_CP23_same_sex_age_matched_5y_dataset_STEP4B.csv", "subject_id",
    file.path(a, "06_age_sensitivity_caliper_interaction", "SDY2583_CP23_same_sex_age_matched_10y_dataset_STEP4B.csv"), "SDY2583_CP23_same_sex_age_matched_10y_dataset_STEP4B.csv", "subject_id",
    file.path(a, "06_age_sensitivity_caliper_interaction", "SDY2583_CP23_same_sex_age_matched_5y_results_STEP4B.csv"), "SDY2583_CP23_same_sex_age_matched_5y_results_STEP4B.csv", "composite_score",
    file.path(a, "06_age_sensitivity_caliper_interaction", "SDY2583_CP23_same_sex_age_matched_10y_results_STEP4B.csv"), "SDY2583_CP23_same_sex_age_matched_10y_results_STEP4B.csv", "composite_score",
    file.path(a, "06_age_sensitivity_caliper_interaction", "SDY2583_CP23_disease_by_age_interaction_results_STEP4B.csv"), "SDY2583_CP23_disease_by_age_interaction_results_STEP4B.csv", "composite_score",
    file.path(a, "07_threshold_sensitivity", "SDY2583_CP23_threshold_score_data_STEP5.csv"), "SDY2583_CP23_threshold_score_data_STEP5.csv", "threshold_set,subject_id",
    file.path(a, "07_threshold_sensitivity", "SDY2583_CP23_threshold_model_results_STEP5.csv"), "SDY2583_CP23_threshold_model_results_STEP5.csv", "threshold_set,variable",
    file.path(a, "07_threshold_sensitivity", "SDY2583_CP23_threshold_wide_robustness_STEP5.csv"), "SDY2583_CP23_threshold_wide_robustness_STEP5.csv", "variable"
  )
  all_reference_csv <- list.files(
    reference_dir, "\\.csv$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE
  )
  archive_checks <- purrr::pmap_dfr(comparison_plan, function(generated, reference_name, keys) {
    hit <- all_reference_csv[basename(all_reference_csv) == reference_name]
    if (length(hit) == 0L) {
      return(tibble::tibble(
        check = paste0("archive table: ", reference_name),
        observed = "reference missing",
        expected = "matched",
        pass = FALSE
      ))
    }
    cmp <- compare_archive_csv(
      generated, hit[1], strsplit(keys, ",", fixed = TRUE)[[1]], tolerance = tol
    )
    tibble::tibble(
      check = paste0("archive table: ", reference_name),
      observed = cmp$detail[1],
      expected = "matched",
      pass = cmp$pass[1]
    )
  })
  checks <- dplyr::bind_rows(checks, archive_checks)
}

report_file <- file.path(out, "SDY2583_CP23_recovered_pipeline_validation_report.csv")
readr::write_csv(checks, report_file)
print(checks, n = nrow(checks))
if (any(!checks$pass | is.na(checks$pass))) {
  stop("CP23 validation failed. Review ", report_file)
}
message("CP23 validation passed.")
