# ============================================================
# SDY2583 CP22
# STEP 3B SAFE: Age/sex-adjusted cancer-vs-healthy feature models
# B-cell / humoral remodeling panel
#
# Primary model:
#   feature ~ disease_group + age_for_model + sex
#
# Sensitivities:
#   1) binary-sex only: Female/Male
#   2) B-cell event QC: n_cd19_b >= 300
#   3) standard dump marker only: Viability_CD3_CD7_CD13
# ============================================================

rm(list = ls())
source(file.path(Sys.getenv("SDY2583_REPO_ROOT", unset = "."), "R", "shared", "bootstrap.R"))

# ------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------

cran_pkgs <- c("dplyr", "readr", "stringr", "tibble", "broom", "purrr")

for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(broom)
  library(purrr)
})

# ------------------------------------------------------------
# 2. Paths
# ------------------------------------------------------------

analysis_dir <- sd_analysis_dir("CP22")

step3a_rdata <- file.path(
  analysis_dir,
  "11_RData",
  "SDY2583_CP22_STEP3A_metadata_merge_age_QC.RData"
)

if (!file.exists(step3a_rdata)) {
  stop("CP22 Step 3A RData bulunamadı: ", step3a_rdata)
}

out_dir <- file.path(analysis_dir, "04_age_sex_adjusted_statistics")
rdata_dir <- file.path(analysis_dir, "11_RData")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rdata_dir, recursive = TRUE, showWarnings = FALSE)

load(step3a_rdata)

if (!exists("cp22_analysis_data")) {
  stop("Step 3A RData içinde cp22_analysis_data bulunamadı.")
}

dat <- cp22_analysis_data

# ------------------------------------------------------------
# 3. Feature dictionary
# ------------------------------------------------------------

feature_dictionary <- tibble::tribble(
  ~feature, ~module, ~interpretation_direction_if_higher_in_cancer,

  "pct_cd19_b_total", "B_cell_composition", "higher total CD19+ B-cell representation",
  "pct_cd19_b_within_dump_low", "B_cell_composition", "higher CD19+ B-cell representation within dump-low compartment",

  "pct_naive_like_within_b", "core_B_cell_states", "higher naive-like B-cell fraction",
  "pct_unswitched_memory_like_within_b", "core_B_cell_states", "higher unswitched memory-like B-cell fraction",
  "pct_switched_memory_like_within_b", "core_B_cell_states", "higher switched memory-like B-cell fraction",
  "pct_double_negative_like_within_b", "core_B_cell_states", "higher double-negative-like B-cell fraction",
  "pct_cd27_pos_within_b", "core_B_cell_states", "higher CD27+ memory-like B-cell fraction",
  "pct_igd_pos_within_b", "core_B_cell_states", "higher IgD+ B-cell fraction",

  "pct_igm_pos_within_b", "immunoglobulin_isotype_B_cells", "higher IgM+ B-cell fraction",
  "pct_iga_pos_within_b", "immunoglobulin_isotype_B_cells", "higher IgA+ B-cell fraction",
  "pct_igg_pos_within_b", "immunoglobulin_isotype_B_cells", "higher IgG+ B-cell fraction",
  "pct_igm_unswitched_memory_like_within_b", "immunoglobulin_isotype_B_cells", "higher IgM+ unswitched memory-like B-cell fraction",
  "pct_iga_switched_memory_like_within_b", "immunoglobulin_isotype_B_cells", "higher IgA+ switched memory-like B-cell fraction",
  "pct_igg_switched_memory_like_within_b", "immunoglobulin_isotype_B_cells", "higher IgG+ switched memory-like B-cell fraction",
  "pct_iga_igg_double_negative_switched_like_within_b", "immunoglobulin_isotype_B_cells", "higher IgA-IgG- switched-memory-like fraction",
  "pct_iga_within_switched_memory_like", "immunoglobulin_isotype_B_cells", "higher IgA fraction within switched memory-like B cells",
  "pct_igg_within_switched_memory_like", "immunoglobulin_isotype_B_cells", "higher IgG fraction within switched memory-like B cells",
  "pct_iga_igg_double_negative_within_switched_memory_like", "immunoglobulin_isotype_B_cells", "higher IgA-IgG- fraction within switched memory-like B cells",
  "pct_igm_within_unswitched_memory_like", "immunoglobulin_isotype_B_cells", "higher IgM fraction within unswitched memory-like B cells",

  "pct_cd38_pos_within_b", "plasmablast_plasma_cell_axis", "higher CD38+ activated/plasmablast-like B-cell fraction",
  "pct_cd38high_within_b", "plasmablast_plasma_cell_axis", "higher CD38-high B-cell fraction",
  "pct_cd138_pos_within_b", "plasmablast_plasma_cell_axis", "higher CD138+ plasma-cell-like B-cell fraction",
  "pct_plasmablast_like_within_b", "plasmablast_plasma_cell_axis", "higher CD38-high CD27+ plasmablast-like B-cell fraction",
  "pct_class_switched_plasmablast_like_within_b", "plasmablast_plasma_cell_axis", "higher class-switched plasmablast-like B-cell fraction",
  "pct_iga_plasmablast_like_within_b", "plasmablast_plasma_cell_axis", "higher IgA+ plasmablast-like B-cell fraction",
  "pct_igg_plasmablast_like_within_b", "plasmablast_plasma_cell_axis", "higher IgG+ plasmablast-like B-cell fraction",
  "pct_plasma_cell_like_within_b", "plasmablast_plasma_cell_axis", "higher CD138+ plasma-cell-like B-cell fraction",
  "pct_cd38_cd138_plasma_cell_like_within_b", "plasmablast_plasma_cell_axis", "higher CD38+CD138+ plasma-cell-like B-cell fraction",

  "pct_cd24_pos_within_b", "transitional_regulatory_like_B_cells", "higher CD24+ B-cell fraction",
  "pct_cd10_pos_within_b", "transitional_regulatory_like_B_cells", "higher CD10+ immature/transitional-like B-cell fraction",
  "pct_transitional_like_within_b", "transitional_regulatory_like_B_cells", "higher CD24+CD38+ transitional-like B-cell fraction",
  "pct_immature_transitional_like_within_b", "transitional_regulatory_like_B_cells", "higher CD10+ immature/transitional-like B-cell fraction",
  "pct_cd10_cd24_cd38_transitional_like_within_b", "transitional_regulatory_like_B_cells", "higher CD10+CD24+CD38+ transitional-like B-cell fraction",
  "pct_cd39_pos_within_b", "transitional_regulatory_like_B_cells", "higher CD39+ regulatory-like B-cell fraction",
  "pct_cd39_regulatory_like_within_b", "transitional_regulatory_like_B_cells", "higher CD39+ regulatory-like B-cell fraction",
  "pct_cd39_cd24_regulatory_like_within_b", "transitional_regulatory_like_B_cells", "higher CD39+CD24+ regulatory-like B-cell fraction",
  "pct_cd39_cd24_cd38_regulatory_transitional_like_within_b", "transitional_regulatory_like_B_cells", "higher CD39+CD24+CD38+ regulatory/transitional-like B-cell fraction",

  "median_CD19_in_B", "B_cell_marker_medians", "higher CD19 median intensity within B cells",
  "median_CD27_in_B", "B_cell_marker_medians", "higher CD27 median intensity within B cells",
  "median_IgD_in_B", "B_cell_marker_medians", "higher IgD median intensity within B cells",
  "median_IgM_in_B", "B_cell_marker_medians", "higher IgM median intensity within B cells",
  "median_IgA_in_B", "B_cell_marker_medians", "higher IgA median intensity within B cells",
  "median_IgG_in_B", "B_cell_marker_medians", "higher IgG median intensity within B cells",
  "median_CD38_in_B", "B_cell_marker_medians", "higher CD38 median intensity within B cells",
  "median_CD138_in_B", "B_cell_marker_medians", "higher CD138 median intensity within B cells",
  "median_CD24_in_B", "B_cell_marker_medians", "higher CD24 median intensity within B cells",
  "median_CD10_in_B", "B_cell_marker_medians", "higher CD10 median intensity within B cells",
  "median_CD39_in_B", "B_cell_marker_medians", "higher CD39 median intensity within B cells"
)

# Keep only features that exist in current data
feature_dictionary <- feature_dictionary %>%
  dplyr::filter(feature %in% names(dat))

target_features <- feature_dictionary$feature

if (length(target_features) == 0) {
  stop("No target features found in cp22_analysis_data.")
}

# ------------------------------------------------------------
# 4. Helper functions
# ------------------------------------------------------------

safe_mean <- function(x) {
  if (sum(!is.na(x)) == 0) return(NA_real_)
  mean(x, na.rm = TRUE)
}

safe_sd <- function(x) {
  if (sum(!is.na(x)) <= 1) return(NA_real_)
  sd(x, na.rm = TRUE)
}

run_lm_feature <- function(df, feature_name, model_label = "primary") {

  model_df <- df %>%
    dplyr::transmute(
      value = suppressWarnings(as.numeric(.data[[feature_name]])),
      disease_group = disease_group,
      age_for_model = suppressWarnings(as.numeric(age_for_model)),
      sex = sex
    ) %>%
    dplyr::filter(
      !is.na(value),
      !is.na(disease_group),
      !is.na(age_for_model),
      !is.na(sex)
    ) %>%
    dplyr::mutate(
      disease_group = factor(
        as.character(disease_group),
        levels = c("Healthy control", "Cancer patient")
      ),
      sex = droplevels(factor(as.character(sex)))
    )

  n_model <- nrow(model_df)
  n_healthy <- sum(model_df$disease_group == "Healthy control")
  n_cancer <- sum(model_df$disease_group == "Cancer patient")

  if (n_model < 50 || n_healthy < 20 || n_cancer < 20) {
    return(NULL)
  }

  if (length(unique(model_df$disease_group)) < 2) {
    return(NULL)
  }

  fit <- tryCatch(
    stats::lm(value ~ disease_group + age_for_model + sex, data = model_df),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    return(NULL)
  }

  tt <- tryCatch(
    broom::tidy(fit, conf.int = TRUE),
    error = function(e) NULL
  )

  if (is.null(tt)) {
    return(NULL)
  }

  disease_term <- "disease_groupCancer patient"

  if (!(disease_term %in% tt$term)) {
    return(NULL)
  }

  out <- tt %>%
    dplyr::filter(term == disease_term)

  desc <- model_df %>%
    dplyr::group_by(disease_group) %>%
    dplyr::summarise(
      n = dplyr::n(),
      mean = mean(value, na.rm = TRUE),
      sd = sd(value, na.rm = TRUE),
      median = median(value, na.rm = TRUE),
      q1 = as.numeric(quantile(value, 0.25, na.rm = TRUE)),
      q3 = as.numeric(quantile(value, 0.75, na.rm = TRUE)),
      .groups = "drop"
    )

  healthy_mean <- desc$mean[desc$disease_group == "Healthy control"]
  cancer_mean <- desc$mean[desc$disease_group == "Cancer patient"]
  healthy_median <- desc$median[desc$disease_group == "Healthy control"]
  cancer_median <- desc$median[desc$disease_group == "Cancer patient"]

  if (length(healthy_mean) == 0) healthy_mean <- NA_real_
  if (length(cancer_mean) == 0) cancer_mean <- NA_real_
  if (length(healthy_median) == 0) healthy_median <- NA_real_
  if (length(cancer_median) == 0) cancer_median <- NA_real_

  tibble::tibble(
    model_label = model_label,
    feature = feature_name,
    n_model = n_model,
    n_healthy = n_healthy,
    n_cancer = n_cancer,
    healthy_mean = healthy_mean,
    cancer_mean = cancer_mean,
    healthy_median = healthy_median,
    cancer_median = cancer_median,
    beta_cancer_vs_healthy = out$estimate[1],
    conf_low = out$conf.low[1],
    conf_high = out$conf.high[1],
    p_value = out$p.value[1],
    direction = dplyr::case_when(
      out$estimate[1] > 0 ~ "higher_in_cancer",
      out$estimate[1] < 0 ~ "lower_in_cancer",
      TRUE ~ "no_direction"
    )
  )
}

add_fdr_and_dictionary <- function(results_df) {

  if (is.null(results_df) || nrow(results_df) == 0) {
    return(results_df)
  }

  results_df %>%
    dplyr::left_join(feature_dictionary, by = "feature") %>%
    dplyr::group_by(module) %>%
    dplyr::mutate(FDR_within_module = p.adjust(p_value, method = "BH")) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(FDR_global = p.adjust(p_value, method = "BH")) %>%
    dplyr::arrange(FDR_global, p_value)
}

# ------------------------------------------------------------
# 5. Primary age/sex-adjusted models
# ------------------------------------------------------------

primary_df <- dat %>%
  dplyr::filter(
    feature_ok == TRUE,
    !is.na(disease_group),
    !is.na(age_for_model),
    !is.na(sex)
  )

primary_results <- dplyr::bind_rows(
  lapply(target_features, function(ff) {
    run_lm_feature(primary_df, ff, model_label = "primary_age_sex_adjusted")
  })
) %>%
  add_fdr_and_dictionary()

# ------------------------------------------------------------
# 6. Binary-sex sensitivity
# ------------------------------------------------------------

binary_sex_df <- primary_df %>%
  dplyr::filter(as.character(sex) %in% c("Female", "Male")) %>%
  dplyr::mutate(sex = factor(as.character(sex), levels = c("Female", "Male")))

binary_sex_results <- dplyr::bind_rows(
  lapply(target_features, function(ff) {
    run_lm_feature(binary_sex_df, ff, model_label = "binary_sex_sensitivity")
  })
) %>%
  add_fdr_and_dictionary()

# ------------------------------------------------------------
# 7. B-cell event-QC sensitivity: n_cd19_b >= 300
# ------------------------------------------------------------

event_qc_df <- primary_df %>%
  dplyr::filter(!is.na(n_cd19_b), n_cd19_b >= 300)

event_qc_results <- dplyr::bind_rows(
  lapply(target_features, function(ff) {
    run_lm_feature(event_qc_df, ff, model_label = "event_QC_n_cd19_b_ge_300")
  })
) %>%
  add_fdr_and_dictionary()

# ------------------------------------------------------------
# 8. Standard dump marker sensitivity
# ------------------------------------------------------------

standard_dump_df <- primary_df %>%
  dplyr::filter(dump_marker_variant == "Viability_CD3_CD7_CD13")

standard_dump_results <- dplyr::bind_rows(
  lapply(target_features, function(ff) {
    run_lm_feature(standard_dump_df, ff, model_label = "standard_dump_marker_only")
  })
) %>%
  add_fdr_and_dictionary()

# ------------------------------------------------------------
# 9. Direction preservation summary
# ------------------------------------------------------------

direction_sensitivity_summary <- primary_results %>%
  dplyr::select(
    feature,
    module,
    primary_beta = beta_cancer_vs_healthy,
    primary_direction = direction,
    primary_FDR_global = FDR_global,
    primary_FDR_within_module = FDR_within_module
  ) %>%
  dplyr::left_join(
    binary_sex_results %>%
      dplyr::select(
        feature,
        binary_beta = beta_cancer_vs_healthy,
        binary_direction = direction,
        binary_FDR_global = FDR_global,
        binary_FDR_within_module = FDR_within_module
      ),
    by = "feature"
  ) %>%
  dplyr::left_join(
    event_qc_results %>%
      dplyr::select(
        feature,
        event_QC_beta = beta_cancer_vs_healthy,
        event_QC_direction = direction,
        event_QC_FDR_global = FDR_global,
        event_QC_FDR_within_module = FDR_within_module
      ),
    by = "feature"
  ) %>%
  dplyr::left_join(
    standard_dump_results %>%
      dplyr::select(
        feature,
        standard_dump_beta = beta_cancer_vs_healthy,
        standard_dump_direction = direction,
        standard_dump_FDR_global = FDR_global,
        standard_dump_FDR_within_module = FDR_within_module
      ),
    by = "feature"
  ) %>%
  dplyr::mutate(
    direction_preserved_binary_sex = primary_direction == binary_direction,
    direction_preserved_event_QC = primary_direction == event_QC_direction,
    direction_preserved_standard_dump = primary_direction == standard_dump_direction,
    direction_preserved_all_sensitivities =
      direction_preserved_binary_sex == TRUE &
      direction_preserved_event_QC == TRUE &
      direction_preserved_standard_dump == TRUE,
    FDR_global_preserved_all_sensitivities =
      primary_FDR_global < 0.05 &
      binary_FDR_global < 0.05 &
      event_QC_FDR_global < 0.05 &
      standard_dump_FDR_global < 0.05
  ) %>%
  dplyr::arrange(primary_FDR_global, primary_FDR_within_module)

# ------------------------------------------------------------
# 10. Summaries
# ------------------------------------------------------------

primary_significance_summary <- primary_results %>%
  dplyr::group_by(module) %>%
  dplyr::summarise(
    n_features = dplyr::n(),
    n_nominal_p_lt_0p05 = sum(p_value < 0.05, na.rm = TRUE),
    n_module_FDR_lt_0p05 = sum(FDR_within_module < 0.05, na.rm = TRUE),
    n_global_FDR_lt_0p05 = sum(FDR_global < 0.05, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(desc(n_global_FDR_lt_0p05), module)

overall_model_summary <- tibble::tibble(
  n_primary_model_subjects = nrow(primary_df),
  n_primary_healthy = sum(primary_df$disease_group == "Healthy control"),
  n_primary_cancer = sum(primary_df$disease_group == "Cancer patient"),
  n_binary_sex_subjects = nrow(binary_sex_df),
  n_event_QC_subjects = nrow(event_qc_df),
  n_standard_dump_subjects = nrow(standard_dump_df),
  n_features_tested = nrow(primary_results),
  n_primary_global_FDR_lt_0p05 = sum(primary_results$FDR_global < 0.05, na.rm = TRUE),
  n_primary_module_FDR_lt_0p05 = sum(primary_results$FDR_within_module < 0.05, na.rm = TRUE),
  n_direction_preserved_all_sensitivities = sum(direction_sensitivity_summary$direction_preserved_all_sensitivities == TRUE, na.rm = TRUE),
  n_global_FDR_preserved_all_sensitivities = sum(direction_sensitivity_summary$FDR_global_preserved_all_sensitivities == TRUE, na.rm = TRUE)
)

top_primary_results <- primary_results %>%
  dplyr::arrange(FDR_global, p_value) %>%
  dplyr::slice_head(n = 30)

top_robust_results <- direction_sensitivity_summary %>%
  dplyr::filter(direction_preserved_all_sensitivities == TRUE) %>%
  dplyr::arrange(primary_FDR_global, primary_FDR_within_module) %>%
  dplyr::slice_head(n = 30)

# ------------------------------------------------------------
# 11. Save outputs
# ------------------------------------------------------------

readr::write_csv(
  feature_dictionary,
  file.path(out_dir, "SDY2583_CP22_feature_dictionary_STEP3B.csv")
)

readr::write_csv(
  primary_results,
  file.path(out_dir, "SDY2583_CP22_primary_age_sex_adjusted_results_STEP3B.csv")
)

readr::write_csv(
  binary_sex_results,
  file.path(out_dir, "SDY2583_CP22_binary_sex_sensitivity_results_STEP3B.csv")
)

readr::write_csv(
  event_qc_results,
  file.path(out_dir, "SDY2583_CP22_event_QC_sensitivity_results_STEP3B.csv")
)

readr::write_csv(
  standard_dump_results,
  file.path(out_dir, "SDY2583_CP22_standard_dump_sensitivity_results_STEP3B.csv")
)

readr::write_csv(
  direction_sensitivity_summary,
  file.path(out_dir, "SDY2583_CP22_direction_sensitivity_summary_STEP3B.csv")
)

readr::write_csv(
  primary_significance_summary,
  file.path(out_dir, "SDY2583_CP22_primary_significance_summary_by_module_STEP3B.csv")
)

readr::write_csv(
  overall_model_summary,
  file.path(out_dir, "SDY2583_CP22_overall_model_summary_STEP3B.csv")
)

readr::write_csv(
  top_primary_results,
  file.path(out_dir, "SDY2583_CP22_top_primary_results_STEP3B.csv")
)

readr::write_csv(
  top_robust_results,
  file.path(out_dir, "SDY2583_CP22_top_direction_robust_results_STEP3B.csv")
)

save(
  feature_dictionary,
  primary_results,
  binary_sex_results,
  event_qc_results,
  standard_dump_results,
  direction_sensitivity_summary,
  primary_significance_summary,
  overall_model_summary,
  top_primary_results,
  top_robust_results,
  primary_df,
  binary_sex_df,
  event_qc_df,
  standard_dump_df,
  file = file.path(
    rdata_dir,
    "SDY2583_CP22_STEP3B_age_sex_adjusted_statistics.RData"
  )
)

# ------------------------------------------------------------
# 12. Console output
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SDY2583 CP22 STEP 3B COMPLETE: AGE/SEX-ADJUSTED STATISTICS\n")
cat("============================================================\n")

cat("\nOverall model summary:\n")
print(as.data.frame(overall_model_summary), row.names = FALSE)

cat("\nPrimary significance summary by module:\n")
print(as.data.frame(primary_significance_summary), row.names = FALSE)

cat("\nTop primary age/sex-adjusted results:\n")
print(as.data.frame(top_primary_results), row.names = FALSE)

cat("\nTop direction-robust sensitivity results:\n")
print(as.data.frame(top_robust_results), row.names = FALSE)

cat("\nFiles saved in:\n")
print(out_dir)

cat("============================================================\n")
