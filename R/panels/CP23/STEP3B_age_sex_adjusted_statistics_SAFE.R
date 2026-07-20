# ============================================================
# SDY2583 CP23
# STEP 3B SAFE: Age/sex-adjusted feature models
#
# Purpose:
#   Test CP23 threshold-defined monocyte/macrophage-like myeloid
#   features for cancer-vs-healthy differences using age/sex-
#   adjusted linear models.
#
# Primary model:
#   feature ~ disease_group + age_for_model + sex
#
# Sensitivity models:
#   1) Binary-sex sensitivity:
#        feature ~ disease_group + age_for_model + sex_binary
#   2) CD45 event-count QC:
#        n_cd45_dump_low >= 1000
#   3) CD33+HLA-DR+ myeloid/APC-like event-count QC:
#        n_cd33_hladr_myeloid_like >= 300
#   4) Technical QC:
#        Exclude channel/marker-order mismatch files
#
# Output:
#   outputs/CP23/04_age_sex_adjusted_models
#   outputs/CP23/11_RData
# ============================================================

rm(list = ls())
source(file.path(Sys.getenv("SDY2583_REPO_ROOT", unset = "."), "R", "shared", "bootstrap.R"))

# ------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------

cran_pkgs <- c("dplyr", "readr", "stringr", "tibble", "broom", "purrr", "tidyr")

for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(broom)
  library(purrr)
  library(tidyr)
})

# Namespace safety
filter <- dplyr::filter
select <- dplyr::select
mutate <- dplyr::mutate
arrange <- dplyr::arrange
summarise <- dplyr::summarise
group_by <- dplyr::group_by
ungroup <- dplyr::ungroup
count <- dplyr::count
distinct <- dplyr::distinct
case_when <- dplyr::case_when
bind_rows <- dplyr::bind_rows
left_join <- dplyr::left_join
n_distinct <- dplyr::n_distinct
n <- dplyr::n

# ------------------------------------------------------------
# 2. Paths and load Step 3A
# ------------------------------------------------------------

analysis_dir <- sd_analysis_dir("CP23")
rdata_dir <- file.path(analysis_dir, "11_RData")
out_dir <- file.path(analysis_dir, "04_age_sex_adjusted_models")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

step3a_rdata <- file.path(rdata_dir, "SDY2583_CP23_STEP3A_metadata_merge_age_QC.RData")

if (!file.exists(step3a_rdata)) {
  stop("Step 3A RData bulunamadı: ", step3a_rdata)
}

load(step3a_rdata)

# Reset paths after RData load.
analysis_dir <- sd_analysis_dir("CP23")
rdata_dir <- file.path(analysis_dir, "11_RData")
out_dir <- file.path(analysis_dir, "04_age_sex_adjusted_models")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!exists("cp23_analysis_data")) stop("cp23_analysis_data bulunamadı.")
if (!exists("feature_dictionary_full")) stop("feature_dictionary_full bulunamadı.")

# ------------------------------------------------------------
# 3. Feature list
# ------------------------------------------------------------

feature_dictionary <- feature_dictionary_full %>%
  filter(feature %in% names(cp23_analysis_data)) %>%
  distinct(feature, .keep_all = TRUE)

feature_cols <- feature_dictionary$feature

if (length(feature_cols) == 0) {
  stop("Test edilecek CP23 feature bulunamadı.")
}

# ------------------------------------------------------------
# 4. Model helper
# ------------------------------------------------------------

run_feature_model <- function(df, feature_name, model_label, sex_var = "sex") {

  if (!(feature_name %in% names(df))) return(NULL)

  work <- df %>%
    transmute(
      value = suppressWarnings(as.numeric(.data[[feature_name]])),
      disease_group = disease_group,
      age_for_model = age_for_model,
      sex_model = .data[[sex_var]]
    ) %>%
    filter(
      !is.na(value),
      !is.na(disease_group),
      !is.na(age_for_model),
      !is.na(sex_model)
    ) %>%
    mutate(
      disease_group = factor(
        as.character(disease_group),
        levels = c("Healthy control", "Cancer patient")
      ),
      sex_model = factor(sex_model)
    )

  if (nrow(work) < 30) return(NULL)
  if (length(unique(work$disease_group)) < 2) return(NULL)

  n_healthy <- sum(work$disease_group == "Healthy control")
  n_cancer <- sum(work$disease_group == "Cancer patient")

  if (n_healthy < 10 || n_cancer < 10) return(NULL)

  fit <- tryCatch(
    lm(value ~ disease_group + age_for_model + sex_model, data = work),
    error = function(e) NULL
  )

  if (is.null(fit)) return(NULL)

  tt <- tryCatch(broom::tidy(fit, conf.int = TRUE), error = function(e) NULL)
  if (is.null(tt)) return(NULL)

  disease_term <- "disease_groupCancer patient"
  if (!(disease_term %in% tt$term)) return(NULL)

  term_row <- tt %>% filter(term == disease_term)

  tibble(
    model_label = model_label,
    feature = feature_name,
    n_model = nrow(work),
    n_healthy = n_healthy,
    n_cancer = n_cancer,
    healthy_mean = mean(work$value[work$disease_group == "Healthy control"], na.rm = TRUE),
    cancer_mean = mean(work$value[work$disease_group == "Cancer patient"], na.rm = TRUE),
    healthy_median = median(work$value[work$disease_group == "Healthy control"], na.rm = TRUE),
    cancer_median = median(work$value[work$disease_group == "Cancer patient"], na.rm = TRUE),
    beta_cancer_vs_healthy = term_row$estimate[1],
    conf_low = term_row$conf.low[1],
    conf_high = term_row$conf.high[1],
    p_value = term_row$p.value[1],
    direction = case_when(
      term_row$estimate[1] > 0 ~ "higher_in_cancer",
      term_row$estimate[1] < 0 ~ "lower_in_cancer",
      TRUE ~ "no_direction"
    )
  )
}

add_fdr_and_dictionary <- function(res_df) {
  if (is.null(res_df) || nrow(res_df) == 0) return(tibble())

  res_df %>%
    left_join(feature_dictionary, by = "feature") %>%
    group_by(module) %>%
    mutate(FDR_within_module = p.adjust(p_value, method = "BH")) %>%
    ungroup() %>%
    mutate(FDR_global = p.adjust(p_value, method = "BH")) %>%
    arrange(FDR_global, p_value)
}

# ------------------------------------------------------------
# 5. Model datasets
# ------------------------------------------------------------

primary_df <- cp23_analysis_data %>%
  filter(model_ready_age_sex == TRUE)

binary_sex_df <- cp23_analysis_data %>%
  filter(model_ready_binary_sex == TRUE)

cd45_event_qc_df <- cp23_analysis_data %>%
  filter(
    model_ready_age_sex == TRUE,
    !is.na(n_cd45_dump_low),
    n_cd45_dump_low >= 1000
  )

hladr_event_qc_df <- cp23_analysis_data %>%
  filter(
    model_ready_age_sex == TRUE,
    !is.na(n_cd33_hladr_myeloid_like),
    n_cd33_hladr_myeloid_like >= 300
  )

technical_qc_df <- cp23_analysis_data %>%
  filter(
    model_ready_age_sex == TRUE,
    channel_order_mismatch_file != TRUE,
    marker_order_mismatch_file != TRUE,
    marker_set_mismatch_file != TRUE
  )

# ------------------------------------------------------------
# 6. Run primary and sensitivity models
# ------------------------------------------------------------

primary_results <- bind_rows(lapply(feature_cols, function(ff) {
  run_feature_model(primary_df, ff, "primary_age_sex_adjusted", sex_var = "sex")
})) %>%
  add_fdr_and_dictionary()

binary_sex_results <- bind_rows(lapply(feature_cols, function(ff) {
  run_feature_model(binary_sex_df, ff, "binary_sex_sensitivity", sex_var = "sex_binary")
})) %>%
  add_fdr_and_dictionary()

cd45_event_qc_results <- bind_rows(lapply(feature_cols, function(ff) {
  run_feature_model(cd45_event_qc_df, ff, "cd45_event_QC_sensitivity", sex_var = "sex")
})) %>%
  add_fdr_and_dictionary()

hladr_event_qc_results <- bind_rows(lapply(feature_cols, function(ff) {
  run_feature_model(hladr_event_qc_df, ff, "hladr_myeloid_event_QC_sensitivity", sex_var = "sex")
})) %>%
  add_fdr_and_dictionary()

technical_qc_results <- bind_rows(lapply(feature_cols, function(ff) {
  run_feature_model(technical_qc_df, ff, "technical_QC_sensitivity", sex_var = "sex")
})) %>%
  add_fdr_and_dictionary()

all_model_results <- bind_rows(
  primary_results,
  binary_sex_results,
  cd45_event_qc_results,
  hladr_event_qc_results,
  technical_qc_results
)

# ------------------------------------------------------------
# 7. Sensitivity-preservation table
# ------------------------------------------------------------

primary_for_join <- primary_results %>%
  select(
    feature,
    module,
    interpretation,
    primary_beta = beta_cancer_vs_healthy,
    primary_direction = direction,
    primary_FDR_global = FDR_global,
    primary_FDR_within_module = FDR_within_module
  )

binary_for_join <- binary_sex_results %>%
  select(
    feature,
    binary_sex_beta = beta_cancer_vs_healthy,
    binary_sex_direction = direction,
    binary_sex_FDR_global = FDR_global
  )

cd45_for_join <- cd45_event_qc_results %>%
  select(
    feature,
    cd45_event_QC_beta = beta_cancer_vs_healthy,
    cd45_event_QC_direction = direction,
    cd45_event_QC_FDR_global = FDR_global
  )

hladr_for_join <- hladr_event_qc_results %>%
  select(
    feature,
    hladr_event_QC_beta = beta_cancer_vs_healthy,
    hladr_event_QC_direction = direction,
    hladr_event_QC_FDR_global = FDR_global
  )

technical_for_join <- technical_qc_results %>%
  select(
    feature,
    technical_QC_beta = beta_cancer_vs_healthy,
    technical_QC_direction = direction,
    technical_QC_FDR_global = FDR_global
  )

sensitivity_preservation <- primary_for_join %>%
  left_join(binary_for_join, by = "feature") %>%
  left_join(cd45_for_join, by = "feature") %>%
  left_join(hladr_for_join, by = "feature") %>%
  left_join(technical_for_join, by = "feature") %>%
  mutate(
    direction_preserved_all_sensitivities =
      primary_direction == binary_sex_direction &
      primary_direction == cd45_event_QC_direction &
      primary_direction == hladr_event_QC_direction &
      primary_direction == technical_QC_direction,
    FDR_global_preserved_all_sensitivities =
      primary_FDR_global < 0.05 &
      binary_sex_FDR_global < 0.05 &
      cd45_event_QC_FDR_global < 0.05 &
      hladr_event_QC_FDR_global < 0.05 &
      technical_QC_FDR_global < 0.05,
    robustness_class = case_when(
      direction_preserved_all_sensitivities == TRUE &
        FDR_global_preserved_all_sensitivities == TRUE ~ "direction_and_global_FDR_preserved_all_sensitivities",
      direction_preserved_all_sensitivities == TRUE &
        FDR_global_preserved_all_sensitivities == FALSE ~ "direction_preserved_FDR_not_all_sensitivities",
      TRUE ~ "direction_not_preserved"
    )
  ) %>%
  arrange(primary_FDR_global, desc(FDR_global_preserved_all_sensitivities), feature)

robust_primary_results <- sensitivity_preservation %>%
  filter(
    direction_preserved_all_sensitivities == TRUE,
    FDR_global_preserved_all_sensitivities == TRUE
  ) %>%
  arrange(primary_FDR_global)

direction_preserved_fdr_not_all <- sensitivity_preservation %>%
  filter(
    direction_preserved_all_sensitivities == TRUE,
    FDR_global_preserved_all_sensitivities == FALSE
  ) %>%
  arrange(primary_FDR_global)

direction_not_preserved <- sensitivity_preservation %>%
  filter(direction_preserved_all_sensitivities != TRUE) %>%
  arrange(primary_FDR_global)

# ------------------------------------------------------------
# 8. Summaries
# ------------------------------------------------------------

model_subject_summary <- tibble(
  n_primary_model_subjects = nrow(primary_df),
  n_primary_healthy = sum(primary_df$disease_group == "Healthy control", na.rm = TRUE),
  n_primary_cancer = sum(primary_df$disease_group == "Cancer patient", na.rm = TRUE),
  n_binary_sex_subjects = nrow(binary_sex_df),
  n_cd45_event_QC_subjects = nrow(cd45_event_qc_df),
  n_hladr_event_QC_subjects = nrow(hladr_event_qc_df),
  n_technical_QC_subjects = nrow(technical_qc_df),
  n_features_tested = length(feature_cols),
  n_primary_global_FDR_lt_0p05 = sum(primary_results$FDR_global < 0.05, na.rm = TRUE),
  n_primary_module_FDR_lt_0p05 = sum(primary_results$FDR_within_module < 0.05, na.rm = TRUE),
  n_direction_preserved_all_sensitivities = sum(sensitivity_preservation$direction_preserved_all_sensitivities == TRUE, na.rm = TRUE),
  n_FDR_preserved_all_sensitivities = sum(sensitivity_preservation$FDR_global_preserved_all_sensitivities == TRUE, na.rm = TRUE)
)

module_summary <- primary_results %>%
  group_by(module) %>%
  summarise(
    n_features = n(),
    n_nominal_p_lt_0p05 = sum(p_value < 0.05, na.rm = TRUE),
    n_module_FDR_lt_0p05 = sum(FDR_within_module < 0.05, na.rm = TRUE),
    n_global_FDR_lt_0p05 = sum(FDR_global < 0.05, na.rm = TRUE),
    min_p = min(p_value, na.rm = TRUE),
    min_FDR_global = min(FDR_global, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(min_FDR_global)

robustness_count_summary <- sensitivity_preservation %>%
  count(robustness_class, name = "n") %>%
  arrange(desc(n))

# ------------------------------------------------------------
# 9. Save outputs
# ------------------------------------------------------------

write_csv(primary_results, file.path(out_dir, "SDY2583_CP23_primary_age_sex_adjusted_feature_results_STEP3B.csv"))
write_csv(binary_sex_results, file.path(out_dir, "SDY2583_CP23_binary_sex_sensitivity_feature_results_STEP3B.csv"))
write_csv(cd45_event_qc_results, file.path(out_dir, "SDY2583_CP23_cd45_event_QC_sensitivity_feature_results_STEP3B.csv"))
write_csv(hladr_event_qc_results, file.path(out_dir, "SDY2583_CP23_hladr_myeloid_event_QC_sensitivity_feature_results_STEP3B.csv"))
write_csv(technical_qc_results, file.path(out_dir, "SDY2583_CP23_technical_QC_sensitivity_feature_results_STEP3B.csv"))
write_csv(all_model_results, file.path(out_dir, "SDY2583_CP23_all_feature_model_results_STEP3B.csv"))
write_csv(sensitivity_preservation, file.path(out_dir, "SDY2583_CP23_sensitivity_preservation_summary_STEP3B.csv"))
write_csv(robust_primary_results, file.path(out_dir, "SDY2583_CP23_robust_primary_results_preserved_all_sensitivities_STEP3B.csv"))
write_csv(direction_preserved_fdr_not_all, file.path(out_dir, "SDY2583_CP23_direction_preserved_FDR_not_all_sensitivities_STEP3B.csv"))
write_csv(direction_not_preserved, file.path(out_dir, "SDY2583_CP23_direction_not_preserved_STEP3B.csv"))
write_csv(model_subject_summary, file.path(out_dir, "SDY2583_CP23_model_subject_summary_STEP3B.csv"))
write_csv(module_summary, file.path(out_dir, "SDY2583_CP23_module_summary_STEP3B.csv"))
write_csv(robustness_count_summary, file.path(out_dir, "SDY2583_CP23_robustness_count_summary_STEP3B.csv"))

save(
  feature_cols,
  feature_dictionary,
  primary_df,
  binary_sex_df,
  cd45_event_qc_df,
  hladr_event_qc_df,
  technical_qc_df,
  primary_results,
  binary_sex_results,
  cd45_event_qc_results,
  hladr_event_qc_results,
  technical_qc_results,
  all_model_results,
  sensitivity_preservation,
  robust_primary_results,
  direction_preserved_fdr_not_all,
  direction_not_preserved,
  model_subject_summary,
  module_summary,
  robustness_count_summary,
  file = file.path(rdata_dir, "SDY2583_CP23_STEP3B_age_sex_adjusted_statistics.RData")
)

# ------------------------------------------------------------
# 10. Console output
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SDY2583 CP23 STEP 3B COMPLETE: AGE/SEX-ADJUSTED MODELS\n")
cat("============================================================\n")

cat("\nModel subject summary:\n")
print(as.data.frame(model_subject_summary), row.names = FALSE)

cat("\nModule summary:\n")
print(as.data.frame(module_summary), row.names = FALSE)

cat("\nRobustness count summary:\n")
print(as.data.frame(robustness_count_summary), row.names = FALSE)

cat("\nTop primary age/sex-adjusted results:\n")
print(as.data.frame(primary_results %>% slice_head(n = 45)), row.names = FALSE)

cat("\nRobust primary results preserved across all sensitivities:\n")
print(as.data.frame(robust_primary_results %>% slice_head(n = 60)), row.names = FALSE)

cat("\nDirection preserved but FDR not preserved across all sensitivities:\n")
print(as.data.frame(direction_preserved_fdr_not_all %>% slice_head(n = 60)), row.names = FALSE)

cat("\nDirection not preserved across sensitivities:\n")
print(as.data.frame(direction_not_preserved %>% slice_head(n = 60)), row.names = FALSE)

cat("\nOutputs saved in:\n")
print(out_dir)

cat("============================================================\n")
