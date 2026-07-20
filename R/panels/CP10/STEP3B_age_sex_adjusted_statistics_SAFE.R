# ============================================================
# SDY2583 CP10
# STEP 3B SAFE: Age/sex-adjusted cancer vs healthy feature models
#
# Primary model:
#   feature ~ disease_group + age_for_model + sex
#
# Sensitivity models:
#   1) binary-sex sensitivity
#   2) CD45 event-count QC sensitivity: n_cd45_viable >= 1000
#   3) technical sensitivity: exclude channel/marker mismatch files
#
# Output:
# outputs/CP10/04_age_sex_adjusted_models
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
# 2. Paths
# ------------------------------------------------------------

analysis_dir <- sd_analysis_dir("CP10")
rdata_dir <- file.path(analysis_dir, "11_RData")
out_dir <- file.path(analysis_dir, "04_age_sex_adjusted_models")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rdata_dir, recursive = TRUE, showWarnings = FALSE)

step3a_rdata <- file.path(rdata_dir, "SDY2583_CP10_STEP3A_metadata_merge_age_QC.RData")

if (!file.exists(step3a_rdata)) {
  stop("Step 3A RData bulunamadı: ", step3a_rdata)
}

load(step3a_rdata)

if (!exists("cp10_analysis_data")) {
  stop("Step 3A RData içinde cp10_analysis_data bulunamadı.")
}

# ------------------------------------------------------------
# 3. Feature dictionary and target feature list
# ------------------------------------------------------------

if (!exists("feature_dictionary")) {
  stop("feature_dictionary bulunamadı. Step 2/3A RData içinde olmalı.")
}

target_features <- feature_dictionary$feature
target_features <- target_features[target_features %in% names(cp10_analysis_data)]

# Add selected marker medians that are not in feature_dictionary
additional_features <- c(
  "median_HLA_DR_in_CD45",
  "median_CD11c_in_CD45",
  "median_CD13_in_CD45",
  "median_CD66b_in_CD45",
  "median_CCR3_in_CD45",
  "median_CD123_in_CD45",
  "median_CD14_in_CD45"
)

additional_features <- additional_features[additional_features %in% names(cp10_analysis_data)]

target_features <- unique(c(target_features, additional_features))

additional_dictionary <- tibble(
  feature = additional_features,
  module = "marker_medians",
  interpretation = additional_features
)

feature_dictionary_full <- bind_rows(
  feature_dictionary,
  additional_dictionary
) %>%
  distinct(feature, .keep_all = TRUE) %>%
  filter(feature %in% target_features)

# ------------------------------------------------------------
# 4. Modeling helper
# ------------------------------------------------------------

run_lm_feature <- function(df, feature_name, model_label) {

  model_df <- df %>%
    transmute(
      value = suppressWarnings(as.numeric(.data[[feature_name]])),
      disease_group = disease_group,
      age_for_model = suppressWarnings(as.numeric(age_for_model)),
      sex = sex
    ) %>%
    filter(
      !is.na(value),
      !is.na(disease_group),
      !is.na(age_for_model),
      !is.na(sex)
    ) %>%
    mutate(
      disease_group = factor(as.character(disease_group), levels = c("Healthy control", "Cancer patient")),
      sex = droplevels(factor(as.character(sex)))
    )

  n_model <- nrow(model_df)
  n_healthy <- sum(model_df$disease_group == "Healthy control")
  n_cancer <- sum(model_df$disease_group == "Cancer patient")

  if (n_model < 50 || n_healthy < 20 || n_cancer < 20) {
    return(NULL)
  }

  fit <- tryCatch(
    lm(value ~ disease_group + age_for_model + sex, data = model_df),
    error = function(e) NULL
  )

  if (is.null(fit)) return(NULL)

  tt <- tryCatch(
    broom::tidy(fit, conf.int = TRUE),
    error = function(e) NULL
  )

  if (is.null(tt)) return(NULL)

  term <- "disease_groupCancer patient"
  if (!(term %in% tt$term)) return(NULL)

  out <- tt %>% filter(term == !!term)

  desc <- model_df %>%
    group_by(disease_group) %>%
    summarise(
      n = n(),
      mean = mean(value, na.rm = TRUE),
      median = median(value, na.rm = TRUE),
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

  tibble(
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
    direction = case_when(
      out$estimate[1] > 0 ~ "higher_in_cancer",
      out$estimate[1] < 0 ~ "lower_in_cancer",
      TRUE ~ "no_direction"
    )
  )
}

run_model_set <- function(df, model_label) {
  bind_rows(
    lapply(target_features, function(ff) {
      run_lm_feature(df, ff, model_label)
    })
  ) %>%
    left_join(feature_dictionary_full, by = "feature") %>%
    group_by(module) %>%
    mutate(FDR_within_module = p.adjust(p_value, method = "BH")) %>%
    ungroup() %>%
    mutate(FDR_global = p.adjust(p_value, method = "BH")) %>%
    arrange(FDR_global, p_value)
}

# ------------------------------------------------------------
# 5. Run models
# ------------------------------------------------------------

primary_df <- cp10_analysis_data %>%
  filter(model_ready_age_sex == TRUE)

binary_sex_df <- cp10_analysis_data %>%
  filter(
    model_ready_age_sex == TRUE,
    !is.na(sex_binary)
  ) %>%
  mutate(sex = sex_binary)

event_qc_df <- cp10_analysis_data %>%
  filter(
    model_ready_age_sex == TRUE,
    !is.na(n_cd45_viable),
    n_cd45_viable >= 1000
  )

technical_qc_df <- cp10_analysis_data %>%
  filter(
    model_ready_age_sex == TRUE,
    channel_order_mismatch_file == FALSE,
    marker_mismatch_file == FALSE
  )

primary_results <- run_model_set(primary_df, "primary_age_sex_adjusted")
binary_sex_results <- run_model_set(binary_sex_df, "binary_sex_sensitivity")
event_qc_results <- run_model_set(event_qc_df, "CD45_event_QC_ge_1000")
technical_qc_results <- run_model_set(technical_qc_df, "technical_QC_exclude_mismatch_files")

all_model_results <- bind_rows(
  primary_results,
  binary_sex_results,
  event_qc_results,
  technical_qc_results
)

# ------------------------------------------------------------
# 6. Sensitivity summary
# ------------------------------------------------------------

get_result_value <- function(df, feature_name, model_label, col_name) {
  x <- df %>%
    filter(feature == feature_name, model_label == !!model_label) %>%
    pull(!!rlang::sym(col_name))
  if (length(x) == 0) return(NA)
  x[1]
}

make_sensitivity_one <- function(feature_name) {

  primary_beta <- as.numeric(get_result_value(all_model_results, feature_name, "primary_age_sex_adjusted", "beta_cancer_vs_healthy"))
  primary_direction <- as.character(get_result_value(all_model_results, feature_name, "primary_age_sex_adjusted", "direction"))
  primary_FDR <- as.numeric(get_result_value(all_model_results, feature_name, "primary_age_sex_adjusted", "FDR_global"))

  binary_beta <- as.numeric(get_result_value(all_model_results, feature_name, "binary_sex_sensitivity", "beta_cancer_vs_healthy"))
  binary_direction <- as.character(get_result_value(all_model_results, feature_name, "binary_sex_sensitivity", "direction"))
  binary_FDR <- as.numeric(get_result_value(all_model_results, feature_name, "binary_sex_sensitivity", "FDR_global"))

  event_beta <- as.numeric(get_result_value(all_model_results, feature_name, "CD45_event_QC_ge_1000", "beta_cancer_vs_healthy"))
  event_direction <- as.character(get_result_value(all_model_results, feature_name, "CD45_event_QC_ge_1000", "direction"))
  event_FDR <- as.numeric(get_result_value(all_model_results, feature_name, "CD45_event_QC_ge_1000", "FDR_global"))

  tech_beta <- as.numeric(get_result_value(all_model_results, feature_name, "technical_QC_exclude_mismatch_files", "beta_cancer_vs_healthy"))
  tech_direction <- as.character(get_result_value(all_model_results, feature_name, "technical_QC_exclude_mismatch_files", "direction"))
  tech_FDR <- as.numeric(get_result_value(all_model_results, feature_name, "technical_QC_exclude_mismatch_files", "FDR_global"))

  directions <- c(primary_direction, binary_direction, event_direction, tech_direction)
  direction_preserved <- all(!is.na(directions)) && length(unique(directions)) == 1

  tibble(
    feature = feature_name,
    module = feature_dictionary_full$module[match(feature_name, feature_dictionary_full$feature)],
    interpretation = feature_dictionary_full$interpretation[match(feature_name, feature_dictionary_full$feature)],
    primary_beta = primary_beta,
    primary_direction = primary_direction,
    primary_FDR_global = primary_FDR,
    binary_sex_beta = binary_beta,
    binary_sex_direction = binary_direction,
    binary_sex_FDR_global = binary_FDR,
    event_QC_beta = event_beta,
    event_QC_direction = event_direction,
    event_QC_FDR_global = event_FDR,
    technical_QC_beta = tech_beta,
    technical_QC_direction = tech_direction,
    technical_QC_FDR_global = tech_FDR,
    direction_preserved_all_sensitivities = direction_preserved,
    FDR_global_preserved_all_sensitivities =
      direction_preserved &&
      primary_FDR < 0.05 &&
      binary_FDR < 0.05 &&
      event_FDR < 0.05 &&
      tech_FDR < 0.05
  )
}

sensitivity_summary <- bind_rows(lapply(target_features, make_sensitivity_one)) %>%
  arrange(primary_FDR_global)

# ------------------------------------------------------------
# 7. Summaries
# ------------------------------------------------------------

model_subject_summary <- tibble(
  n_primary_model_subjects = nrow(primary_df),
  n_primary_healthy = sum(primary_df$disease_group == "Healthy control"),
  n_primary_cancer = sum(primary_df$disease_group == "Cancer patient"),
  n_binary_sex_subjects = nrow(binary_sex_df),
  n_event_QC_subjects = nrow(event_qc_df),
  n_technical_QC_subjects = nrow(technical_qc_df),
  n_features_tested = length(target_features),
  n_primary_global_FDR_lt_0p05 = sum(primary_results$FDR_global < 0.05, na.rm = TRUE),
  n_primary_module_FDR_lt_0p05 = sum(primary_results$FDR_within_module < 0.05, na.rm = TRUE),
  n_direction_preserved_all_sensitivities = sum(sensitivity_summary$direction_preserved_all_sensitivities == TRUE, na.rm = TRUE),
  n_FDR_preserved_all_sensitivities = sum(sensitivity_summary$FDR_global_preserved_all_sensitivities == TRUE, na.rm = TRUE)
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

top_primary_results <- primary_results %>%
  arrange(FDR_global, p_value)

robust_primary_results <- sensitivity_summary %>%
  filter(FDR_global_preserved_all_sensitivities == TRUE) %>%
  arrange(primary_FDR_global)

direction_preserved_not_FDR <- sensitivity_summary %>%
  filter(
    direction_preserved_all_sensitivities == TRUE,
    FDR_global_preserved_all_sensitivities == FALSE
  ) %>%
  arrange(primary_FDR_global)

direction_not_preserved <- sensitivity_summary %>%
  filter(direction_preserved_all_sensitivities == FALSE) %>%
  arrange(primary_FDR_global)

# ------------------------------------------------------------
# 8. Save outputs
# ------------------------------------------------------------

write_csv(feature_dictionary_full, file.path(out_dir, "SDY2583_CP10_feature_dictionary_full_STEP3B.csv"))
write_csv(primary_results, file.path(out_dir, "SDY2583_CP10_primary_age_sex_adjusted_results_STEP3B.csv"))
write_csv(binary_sex_results, file.path(out_dir, "SDY2583_CP10_binary_sex_sensitivity_results_STEP3B.csv"))
write_csv(event_qc_results, file.path(out_dir, "SDY2583_CP10_CD45_event_QC_sensitivity_results_STEP3B.csv"))
write_csv(technical_qc_results, file.path(out_dir, "SDY2583_CP10_technical_QC_sensitivity_results_STEP3B.csv"))
write_csv(all_model_results, file.path(out_dir, "SDY2583_CP10_all_model_results_STEP3B.csv"))
write_csv(sensitivity_summary, file.path(out_dir, "SDY2583_CP10_sensitivity_summary_STEP3B.csv"))
write_csv(model_subject_summary, file.path(out_dir, "SDY2583_CP10_model_subject_summary_STEP3B.csv"))
write_csv(module_summary, file.path(out_dir, "SDY2583_CP10_module_summary_STEP3B.csv"))
write_csv(robust_primary_results, file.path(out_dir, "SDY2583_CP10_robust_primary_results_STEP3B.csv"))
write_csv(direction_preserved_not_FDR, file.path(out_dir, "SDY2583_CP10_direction_preserved_not_FDR_STEP3B.csv"))
write_csv(direction_not_preserved, file.path(out_dir, "SDY2583_CP10_direction_not_preserved_STEP3B.csv"))

save(
  cp10_analysis_data,
  feature_dictionary_full,
  target_features,
  primary_results,
  binary_sex_results,
  event_qc_results,
  technical_qc_results,
  all_model_results,
  sensitivity_summary,
  model_subject_summary,
  module_summary,
  robust_primary_results,
  direction_preserved_not_FDR,
  direction_not_preserved,
  file = file.path(rdata_dir, "SDY2583_CP10_STEP3B_age_sex_adjusted_statistics.RData")
)

# ------------------------------------------------------------
# 9. Console output
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SDY2583 CP10 STEP 3B COMPLETE: AGE/SEX-ADJUSTED MODELS\n")
cat("============================================================\n")

cat("\nModel subject summary:\n")
print(as.data.frame(model_subject_summary), row.names = FALSE)

cat("\nModule summary:\n")
print(as.data.frame(module_summary), row.names = FALSE)

cat("\nTop primary age/sex-adjusted results:\n")
print(as.data.frame(top_primary_results %>% slice_head(n = 40)), row.names = FALSE)

cat("\nRobust primary results preserved across sensitivities:\n")
print(as.data.frame(robust_primary_results %>% slice_head(n = 40)), row.names = FALSE)

cat("\nDirection preserved but FDR not preserved across all sensitivities:\n")
print(as.data.frame(direction_preserved_not_FDR %>% slice_head(n = 40)), row.names = FALSE)

cat("\nDirection not preserved across sensitivities:\n")
print(as.data.frame(direction_not_preserved), row.names = FALSE)

cat("\nOutputs saved in:\n")
print(out_dir)

cat("============================================================\n")
