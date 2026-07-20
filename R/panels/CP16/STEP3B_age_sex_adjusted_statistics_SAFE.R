# ============================================================
# SDY2583 CP16
# STEP 3B SAFE: Age/sex-adjusted cancer vs healthy feature models
#
# Primary model:
#   feature ~ disease_group + age_for_model + sex
#
# Sensitivity models:
#   1) binary-sex sensitivity
#   2) CD45 dump-low event QC: n_cd45_dump_low >= 300
#   3) HLA-DR APC-like core event QC: n_hladr_apc_core >= 300
#   4) technical QC: exclude channel/marker-order mismatch files
#
# CP16-specific caution:
#   - CP16 APC/DC-like denominators are relatively small.
#   - Event-count sensitivity is essential.
#   - Use phenotype-like terms only:
#       APC-like myeloid, monocyte-like, cDC1-like, cDC2-like, pDC-like.
#
# Output:
#   outputs/CP16/04_age_sex_adjusted_models
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

analysis_dir <- sd_analysis_dir("CP16")
rdata_dir <- file.path(analysis_dir, "11_RData")
out_dir <- file.path(analysis_dir, "04_age_sex_adjusted_models")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rdata_dir, recursive = TRUE, showWarnings = FALSE)

step3a_rdata <- file.path(rdata_dir, "SDY2583_CP16_STEP3A_metadata_merge_age_QC.RData")

if (!file.exists(step3a_rdata)) {
  stop("Step 3A RData bulunamadı: ", step3a_rdata)
}

load(step3a_rdata)

# Reset paths after loading, because loaded RData may contain older out_dir.
analysis_dir <- sd_analysis_dir("CP16")
rdata_dir <- file.path(analysis_dir, "11_RData")
out_dir <- file.path(analysis_dir, "04_age_sex_adjusted_models")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!exists("cp16_analysis_data")) {
  stop("Step 3A RData içinde cp16_analysis_data bulunamadı.")
}

if (!exists("feature_dictionary")) {
  stop("feature_dictionary bulunamadı.")
}

# ------------------------------------------------------------
# 3. Feature list
# ------------------------------------------------------------

target_features <- feature_dictionary$feature
target_features <- target_features[target_features %in% names(cp16_analysis_data)]

# Add marker median features as exploratory marker-intensity endpoints.
additional_features <- c(
  "median_HLA_DR_in_cd45_dump_low",
  "median_CD11c_in_cd45_dump_low",
  "median_CD14_in_cd45_dump_low",
  "median_CD16_in_cd45_dump_low",
  "median_CD1c_in_cd45_dump_low",
  "median_CD123_in_cd45_dump_low",
  "median_CD141_in_cd45_dump_low",
  "median_CLEC9A_in_cd45_dump_low",
  "median_FceRI_in_cd45_dump_low",
  "median_CD13_in_cd45_dump_low",
  "median_HLA_DR_in_cd14_mono_like",
  "median_CD16_in_cd14_mono_like",
  "median_CD1c_in_hladr_core",
  "median_FceRI_in_cd1c_apc_like",
  "median_CLEC9A_in_cd141_pos_hladr",
  "median_CD123_in_hladr_core",
  "median_CD13_in_hladr_core"
)

additional_features <- additional_features[additional_features %in% names(cp16_analysis_data)]

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
  filter(feature %in% unique(c(target_features, additional_features)))

target_features <- unique(c(target_features, additional_features))

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
# 5. Model datasets
# ------------------------------------------------------------

primary_df <- cp16_analysis_data %>%
  filter(model_ready_age_sex == TRUE)

binary_sex_df <- cp16_analysis_data %>%
  filter(
    model_ready_age_sex == TRUE,
    !is.na(sex_binary)
  ) %>%
  mutate(sex = sex_binary)

cd45_event_qc_df <- cp16_analysis_data %>%
  filter(
    model_ready_age_sex == TRUE,
    !is.na(n_cd45_dump_low),
    n_cd45_dump_low >= 300
  )

hladr_event_qc_df <- cp16_analysis_data %>%
  filter(
    model_ready_age_sex == TRUE,
    !is.na(n_hladr_apc_core),
    n_hladr_apc_core >= 300
  )

technical_qc_df <- cp16_analysis_data %>%
  filter(
    model_ready_age_sex == TRUE,
    channel_order_mismatch_file == FALSE,
    marker_order_mismatch_file == FALSE,
    marker_set_mismatch_file == FALSE
  )

# ------------------------------------------------------------
# 6. Run models
# ------------------------------------------------------------

primary_results <- run_model_set(primary_df, "primary_age_sex_adjusted")
binary_sex_results <- run_model_set(binary_sex_df, "binary_sex_sensitivity")
cd45_event_qc_results <- run_model_set(cd45_event_qc_df, "CD45_dump_low_event_QC_ge_300")
hladr_event_qc_results <- run_model_set(hladr_event_qc_df, "HLADR_APC_core_event_QC_ge_300")
technical_qc_results <- run_model_set(technical_qc_df, "technical_QC_exclude_mismatch_files")

all_model_results <- bind_rows(
  primary_results,
  binary_sex_results,
  cd45_event_qc_results,
  hladr_event_qc_results,
  technical_qc_results
)

# ------------------------------------------------------------
# 7. Sensitivity summary
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

  cd45_beta <- as.numeric(get_result_value(all_model_results, feature_name, "CD45_dump_low_event_QC_ge_300", "beta_cancer_vs_healthy"))
  cd45_direction <- as.character(get_result_value(all_model_results, feature_name, "CD45_dump_low_event_QC_ge_300", "direction"))
  cd45_FDR <- as.numeric(get_result_value(all_model_results, feature_name, "CD45_dump_low_event_QC_ge_300", "FDR_global"))

  hladr_beta <- as.numeric(get_result_value(all_model_results, feature_name, "HLADR_APC_core_event_QC_ge_300", "beta_cancer_vs_healthy"))
  hladr_direction <- as.character(get_result_value(all_model_results, feature_name, "HLADR_APC_core_event_QC_ge_300", "direction"))
  hladr_FDR <- as.numeric(get_result_value(all_model_results, feature_name, "HLADR_APC_core_event_QC_ge_300", "FDR_global"))

  tech_beta <- as.numeric(get_result_value(all_model_results, feature_name, "technical_QC_exclude_mismatch_files", "beta_cancer_vs_healthy"))
  tech_direction <- as.character(get_result_value(all_model_results, feature_name, "technical_QC_exclude_mismatch_files", "direction"))
  tech_FDR <- as.numeric(get_result_value(all_model_results, feature_name, "technical_QC_exclude_mismatch_files", "FDR_global"))

  directions <- c(primary_direction, binary_direction, cd45_direction, hladr_direction, tech_direction)
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
    cd45_event_QC_beta = cd45_beta,
    cd45_event_QC_direction = cd45_direction,
    cd45_event_QC_FDR_global = cd45_FDR,
    hladr_event_QC_beta = hladr_beta,
    hladr_event_QC_direction = hladr_direction,
    hladr_event_QC_FDR_global = hladr_FDR,
    technical_QC_beta = tech_beta,
    technical_QC_direction = tech_direction,
    technical_QC_FDR_global = tech_FDR,
    direction_preserved_all_sensitivities = direction_preserved,
    FDR_global_preserved_all_sensitivities =
      direction_preserved &&
      primary_FDR < 0.05 &&
      binary_FDR < 0.05 &&
      cd45_FDR < 0.05 &&
      hladr_FDR < 0.05 &&
      tech_FDR < 0.05
  )
}

sensitivity_summary <- bind_rows(lapply(target_features, make_sensitivity_one)) %>%
  arrange(primary_FDR_global)

# ------------------------------------------------------------
# 8. Summaries
# ------------------------------------------------------------

model_subject_summary <- tibble(
  n_primary_model_subjects = nrow(primary_df),
  n_primary_healthy = sum(primary_df$disease_group == "Healthy control"),
  n_primary_cancer = sum(primary_df$disease_group == "Cancer patient"),
  n_binary_sex_subjects = nrow(binary_sex_df),
  n_cd45_event_QC_subjects = nrow(cd45_event_qc_df),
  n_hladr_event_QC_subjects = nrow(hladr_event_qc_df),
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
# 9. Save outputs
# ------------------------------------------------------------

write_csv(feature_dictionary_full, file.path(out_dir, "SDY2583_CP16_feature_dictionary_full_STEP3B.csv"))
write_csv(primary_results, file.path(out_dir, "SDY2583_CP16_primary_age_sex_adjusted_results_STEP3B.csv"))
write_csv(binary_sex_results, file.path(out_dir, "SDY2583_CP16_binary_sex_sensitivity_results_STEP3B.csv"))
write_csv(cd45_event_qc_results, file.path(out_dir, "SDY2583_CP16_CD45_dump_low_event_QC_results_STEP3B.csv"))
write_csv(hladr_event_qc_results, file.path(out_dir, "SDY2583_CP16_HLADR_APC_core_event_QC_results_STEP3B.csv"))
write_csv(technical_qc_results, file.path(out_dir, "SDY2583_CP16_technical_QC_sensitivity_results_STEP3B.csv"))
write_csv(all_model_results, file.path(out_dir, "SDY2583_CP16_all_model_results_STEP3B.csv"))
write_csv(sensitivity_summary, file.path(out_dir, "SDY2583_CP16_sensitivity_summary_STEP3B.csv"))
write_csv(model_subject_summary, file.path(out_dir, "SDY2583_CP16_model_subject_summary_STEP3B.csv"))
write_csv(module_summary, file.path(out_dir, "SDY2583_CP16_module_summary_STEP3B.csv"))
write_csv(robust_primary_results, file.path(out_dir, "SDY2583_CP16_robust_primary_results_STEP3B.csv"))
write_csv(direction_preserved_not_FDR, file.path(out_dir, "SDY2583_CP16_direction_preserved_not_FDR_STEP3B.csv"))
write_csv(direction_not_preserved, file.path(out_dir, "SDY2583_CP16_direction_not_preserved_STEP3B.csv"))

save(
  cp16_analysis_data,
  feature_dictionary_full,
  target_features,
  primary_results,
  binary_sex_results,
  cd45_event_qc_results,
  hladr_event_qc_results,
  technical_qc_results,
  all_model_results,
  sensitivity_summary,
  model_subject_summary,
  module_summary,
  robust_primary_results,
  direction_preserved_not_FDR,
  direction_not_preserved,
  file = file.path(rdata_dir, "SDY2583_CP16_STEP3B_age_sex_adjusted_statistics.RData")
)

# ------------------------------------------------------------
# 10. Console output
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SDY2583 CP16 STEP 3B COMPLETE: AGE/SEX-ADJUSTED MODELS\n")
cat("============================================================\n")

cat("\nModel subject summary:\n")
print(as.data.frame(model_subject_summary), row.names = FALSE)

cat("\nModule summary:\n")
print(as.data.frame(module_summary), row.names = FALSE)

cat("\nTop primary age/sex-adjusted results:\n")
print(as.data.frame(top_primary_results %>% slice_head(n = 45)), row.names = FALSE)

cat("\nRobust primary results preserved across all sensitivities:\n")
print(as.data.frame(robust_primary_results %>% slice_head(n = 45)), row.names = FALSE)

cat("\nDirection preserved but FDR not preserved across all sensitivities:\n")
print(as.data.frame(direction_preserved_not_FDR %>% slice_head(n = 45)), row.names = FALSE)

cat("\nDirection not preserved across sensitivities:\n")
print(as.data.frame(direction_not_preserved), row.names = FALSE)

cat("\nOutputs saved in:\n")
print(out_dir)

cat("============================================================\n")
