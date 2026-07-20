# ============================================================
# SDY2583 CP23
# STEP 4 SAFE v2: Composite scores
#
# FIX in v2:
#   Step 3B RData does not necessarily contain cp23_analysis_data.
#   Therefore this script first loads Step 3A RData to recover
#   cp23_analysis_data, then loads Step 3B RData for model results.
#
# Composite orientation:
#   Higher composite score = stronger cancer-associated CP23
#   monocyte/macrophage-like myeloid remodeling pattern.
# ============================================================

rm(list = ls())
source(file.path(Sys.getenv("SDY2583_REPO_ROOT", unset = "."), "R", "shared", "bootstrap.R"))

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

analysis_dir <- sd_analysis_dir("CP23")
rdata_dir <- file.path(analysis_dir, "11_RData")
out_dir <- file.path(analysis_dir, "05_composite_scores")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

step3a_rdata <- file.path(rdata_dir, "SDY2583_CP23_STEP3A_metadata_merge_age_QC.RData")
step3b_rdata <- file.path(rdata_dir, "SDY2583_CP23_STEP3B_age_sex_adjusted_statistics.RData")

if (!file.exists(step3a_rdata)) {
  stop("Step 3A RData bulunamadı: ", step3a_rdata)
}

if (!file.exists(step3b_rdata)) {
  stop("Step 3B RData bulunamadı: ", step3b_rdata)
}

# Load Step 3A first: contains cp23_analysis_data.
load(step3a_rdata)

if (!exists("cp23_analysis_data")) {
  stop("Step 3A RData içinde cp23_analysis_data bulunamadı.")
}

# Preserve Step 3A objects before loading Step 3B.
cp23_analysis_data_STEP3A <- cp23_analysis_data
feature_dictionary_full_STEP3A <- if (exists("feature_dictionary_full")) feature_dictionary_full else NULL

# Load Step 3B: contains primary_results and sensitivity results.
load(step3b_rdata)

# Restore objects needed from Step 3A if Step 3B did not carry them.
if (!exists("cp23_analysis_data")) {
  cp23_analysis_data <- cp23_analysis_data_STEP3A
}

if (!exists("feature_dictionary_full") && !is.null(feature_dictionary_full_STEP3A)) {
  feature_dictionary_full <- feature_dictionary_full_STEP3A
}

# Reset paths after RData loads.
analysis_dir <- sd_analysis_dir("CP23")
rdata_dir <- file.path(analysis_dir, "11_RData")
out_dir <- file.path(analysis_dir, "05_composite_scores")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!exists("cp23_analysis_data")) stop("cp23_analysis_data bulunamadı.")
if (!exists("primary_results")) stop("primary_results bulunamadı.")

# ------------------------------------------------------------
# Composite definitions
# ------------------------------------------------------------

composite_definitions <- tribble(
  ~composite_score, ~feature, ~component_direction, ~component_axis,

  "CP23_CD45_dump_low_myeloid_enrichment_score",
  "pct_cd45_dump_low_within_total",
  1,
  "CD45_dump_low_composition",

  "CP23_CD9_CD84_activation_attenuation_score",
  "pct_cd84_pos_within_cd14_mono_like",
  -1,
  "CD9_CD84_activation_attenuation",

  "CP23_CD9_CD84_activation_attenuation_score",
  "median_CD84_in_cd45_dump_low",
  -1,
  "CD9_CD84_activation_attenuation",

  "CP23_CD9_CD84_activation_attenuation_score",
  "pct_cd9_cd84_pos_within_cd14_mono_like",
  -1,
  "CD9_CD84_activation_attenuation",

  "CP23_CD9_CD84_activation_attenuation_score",
  "median_CD9_in_cd14_mono_like",
  -1,
  "CD9_CD84_activation_attenuation",

  "CP23_CD9_CD84_activation_attenuation_score",
  "pct_cd9_pos_within_cd14_mono_like",
  -1,
  "CD9_CD84_activation_attenuation",

  "CP23_CD9_CD84_activation_attenuation_score",
  "pct_cd14_cd9_cd84_pos_within_cd45_dump_low",
  -1,
  "CD9_CD84_activation_attenuation",

  "CP23_FcERI_myeloid_attenuation_score",
  "median_FceRI_in_cd14_mono_like",
  -1,
  "FcERI_myeloid_attenuation",

  "CP23_FcERI_myeloid_attenuation_score",
  "median_FceRI_in_cd45_dump_low",
  -1,
  "FcERI_myeloid_attenuation",

  "CP23_FcERI_myeloid_attenuation_score",
  "pct_fceri_pos_within_cd45_dump_low",
  -1,
  "FcERI_myeloid_attenuation",

  "CP23_FcERI_myeloid_attenuation_score",
  "pct_fceri_pos_within_cd14_mono_like",
  -1,
  "FcERI_myeloid_attenuation",

  "CP23_FcERI_myeloid_attenuation_score",
  "pct_fceri_hladr_apc_like_within_hladr_pos",
  -1,
  "FcERI_myeloid_attenuation",

  "CP23_CD33_HLA_DR_myeloid_repatterning_score",
  "pct_cd33_hladr_myeloid_like_within_hladr_pos",
  -1,
  "CD33_HLA_DR_myeloid_repatterning",

  "CP23_CD33_HLA_DR_myeloid_repatterning_score",
  "pct_cd33_pos_within_cd45_dump_low",
  -1,
  "CD33_HLA_DR_myeloid_repatterning",

  "CP23_CD33_HLA_DR_myeloid_repatterning_score",
  "median_CD206_in_cd33_hladr_myeloid_like",
  1,
  "CD33_HLA_DR_myeloid_repatterning",

  "CP23_integrated_monocyte_macrophage_like_myeloid_remodeling_score",
  "pct_cd45_dump_low_within_total",
  1,
  "integrated",

  "CP23_integrated_monocyte_macrophage_like_myeloid_remodeling_score",
  "pct_cd84_pos_within_cd14_mono_like",
  -1,
  "integrated",

  "CP23_integrated_monocyte_macrophage_like_myeloid_remodeling_score",
  "median_CD84_in_cd45_dump_low",
  -1,
  "integrated",

  "CP23_integrated_monocyte_macrophage_like_myeloid_remodeling_score",
  "pct_cd9_cd84_pos_within_cd14_mono_like",
  -1,
  "integrated",

  "CP23_integrated_monocyte_macrophage_like_myeloid_remodeling_score",
  "median_CD9_in_cd14_mono_like",
  -1,
  "integrated",

  "CP23_integrated_monocyte_macrophage_like_myeloid_remodeling_score",
  "pct_cd9_pos_within_cd14_mono_like",
  -1,
  "integrated",

  "CP23_integrated_monocyte_macrophage_like_myeloid_remodeling_score",
  "pct_cd14_cd9_cd84_pos_within_cd45_dump_low",
  -1,
  "integrated",

  "CP23_integrated_monocyte_macrophage_like_myeloid_remodeling_score",
  "median_FceRI_in_cd14_mono_like",
  -1,
  "integrated",

  "CP23_integrated_monocyte_macrophage_like_myeloid_remodeling_score",
  "median_FceRI_in_cd45_dump_low",
  -1,
  "integrated",

  "CP23_integrated_monocyte_macrophage_like_myeloid_remodeling_score",
  "pct_fceri_pos_within_cd45_dump_low",
  -1,
  "integrated",

  "CP23_integrated_monocyte_macrophage_like_myeloid_remodeling_score",
  "pct_fceri_pos_within_cd14_mono_like",
  -1,
  "integrated",

  "CP23_integrated_monocyte_macrophage_like_myeloid_remodeling_score",
  "pct_fceri_hladr_apc_like_within_hladr_pos",
  -1,
  "integrated",

  "CP23_integrated_monocyte_macrophage_like_myeloid_remodeling_score",
  "pct_cd33_hladr_myeloid_like_within_hladr_pos",
  -1,
  "integrated",

  "CP23_integrated_monocyte_macrophage_like_myeloid_remodeling_score",
  "pct_cd33_pos_within_cd45_dump_low",
  -1,
  "integrated",

  "CP23_integrated_monocyte_macrophage_like_myeloid_remodeling_score",
  "median_CD206_in_cd33_hladr_myeloid_like",
  1,
  "integrated"
)

missing_composite_components <- composite_definitions %>%
  filter(!(feature %in% names(cp23_analysis_data)))

composite_definitions_available <- composite_definitions %>%
  filter(feature %in% names(cp23_analysis_data))

if (nrow(composite_definitions_available) == 0) {
  stop("Composite hesaplanacak uygun component feature bulunamadı.")
}

feature_annotation <- primary_results %>%
  select(feature, module, interpretation, beta_cancer_vs_healthy, FDR_global, FDR_within_module) %>%
  distinct(feature, .keep_all = TRUE)

composite_definitions_annotated <- composite_definitions_available %>%
  left_join(feature_annotation, by = "feature")

# ------------------------------------------------------------
# Build oriented z-score composites
# ------------------------------------------------------------

z_safe <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  mu <- mean(x, na.rm = TRUE)
  sig <- stats::sd(x, na.rm = TRUE)
  if (is.na(sig) || sig == 0) return(rep(NA_real_, length(x)))
  (x - mu) / sig
}

cp23_score_data <- cp23_analysis_data

for (i in seq_len(nrow(composite_definitions_available))) {
  ff <- composite_definitions_available$feature[i]
  direction <- composite_definitions_available$component_direction[i]
  zname <- paste0("z_oriented__", ff)

  if (!(zname %in% names(cp23_score_data))) {
    cp23_score_data[[zname]] <- z_safe(cp23_score_data[[ff]]) * direction
  }
}

score_names <- unique(composite_definitions_available$composite_score)

for (ss in score_names) {
  comp_features <- composite_definitions_available %>%
    filter(composite_score == ss) %>%
    pull(feature)

  zcols <- paste0("z_oriented__", comp_features)
  zcols <- zcols[zcols %in% names(cp23_score_data)]

  cp23_score_data[[ss]] <- rowMeans(cp23_score_data[, zcols, drop = FALSE], na.rm = TRUE)

  all_na <- apply(cp23_score_data[, zcols, drop = FALSE], 1, function(x) all(is.na(x)))
  cp23_score_data[[ss]][all_na] <- NA_real_
}

score_label_map <- tibble(
  composite_score = c(
    "CP23_CD45_dump_low_myeloid_enrichment_score",
    "CP23_CD9_CD84_activation_attenuation_score",
    "CP23_FcERI_myeloid_attenuation_score",
    "CP23_CD33_HLA_DR_myeloid_repatterning_score",
    "CP23_integrated_monocyte_macrophage_like_myeloid_remodeling_score"
  ),
  score_label = c(
    "CD45+ dump-low myeloid enrichment",
    "CD9/CD84 activation attenuation",
    "FcERI-associated myeloid attenuation",
    "CD33/HLA-DR myeloid repatterning",
    "Integrated CP23 monocyte/macrophage-like myeloid remodeling"
  )
)

composite_score_dictionary <- tibble(composite_score = score_names) %>%
  left_join(score_label_map, by = "composite_score") %>%
  mutate(score_label = ifelse(is.na(score_label), composite_score, score_label))

# ------------------------------------------------------------
# Model helper
# ------------------------------------------------------------

run_score_model <- function(df, score_name, model_label, sex_var = "sex") {

  if (!(score_name %in% names(df))) return(NULL)

  work <- df %>%
    transmute(
      value = suppressWarnings(as.numeric(.data[[score_name]])),
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
    composite_score = score_name,
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

add_score_fdr <- function(res_df) {
  if (is.null(res_df) || nrow(res_df) == 0) return(tibble())

  res_df %>%
    left_join(composite_score_dictionary, by = "composite_score") %>%
    mutate(FDR_global = p.adjust(p_value, method = "BH")) %>%
    arrange(FDR_global, p_value)
}

primary_df <- cp23_score_data %>% filter(model_ready_age_sex == TRUE)
binary_sex_df <- cp23_score_data %>% filter(model_ready_binary_sex == TRUE)

cd45_event_qc_df <- cp23_score_data %>%
  filter(model_ready_age_sex == TRUE, !is.na(n_cd45_dump_low), n_cd45_dump_low >= 1000)

hladr_event_qc_df <- cp23_score_data %>%
  filter(model_ready_age_sex == TRUE, !is.na(n_cd33_hladr_myeloid_like), n_cd33_hladr_myeloid_like >= 300)

technical_qc_df <- cp23_score_data %>%
  filter(
    model_ready_age_sex == TRUE,
    channel_order_mismatch_file != TRUE,
    marker_order_mismatch_file != TRUE,
    marker_set_mismatch_file != TRUE
  )

primary_composite_results <- bind_rows(lapply(score_names, function(ss) {
  run_score_model(primary_df, ss, "primary_age_sex_adjusted", sex_var = "sex")
})) %>% add_score_fdr()

binary_sex_composite_results <- bind_rows(lapply(score_names, function(ss) {
  run_score_model(binary_sex_df, ss, "binary_sex_sensitivity", sex_var = "sex_binary")
})) %>% add_score_fdr()

cd45_event_qc_composite_results <- bind_rows(lapply(score_names, function(ss) {
  run_score_model(cd45_event_qc_df, ss, "cd45_event_QC_sensitivity", sex_var = "sex")
})) %>% add_score_fdr()

hladr_event_qc_composite_results <- bind_rows(lapply(score_names, function(ss) {
  run_score_model(hladr_event_qc_df, ss, "hladr_myeloid_event_QC_sensitivity", sex_var = "sex")
})) %>% add_score_fdr()

technical_qc_composite_results <- bind_rows(lapply(score_names, function(ss) {
  run_score_model(technical_qc_df, ss, "technical_QC_sensitivity", sex_var = "sex")
})) %>% add_score_fdr()

all_composite_model_results <- bind_rows(
  primary_composite_results,
  binary_sex_composite_results,
  cd45_event_qc_composite_results,
  hladr_event_qc_composite_results,
  technical_qc_composite_results
)

primary_join <- primary_composite_results %>%
  select(composite_score, score_label, primary_beta = beta_cancer_vs_healthy,
         primary_direction = direction, primary_FDR_global = FDR_global)

binary_join <- binary_sex_composite_results %>%
  select(composite_score, binary_sex_beta = beta_cancer_vs_healthy,
         binary_sex_direction = direction, binary_sex_FDR_global = FDR_global)

cd45_join <- cd45_event_qc_composite_results %>%
  select(composite_score, cd45_event_QC_beta = beta_cancer_vs_healthy,
         cd45_event_QC_direction = direction, cd45_event_QC_FDR_global = FDR_global)

hladr_join <- hladr_event_qc_composite_results %>%
  select(composite_score, hladr_event_QC_beta = beta_cancer_vs_healthy,
         hladr_event_QC_direction = direction, hladr_event_QC_FDR_global = FDR_global)

technical_join <- technical_qc_composite_results %>%
  select(composite_score, technical_QC_beta = beta_cancer_vs_healthy,
         technical_QC_direction = direction, technical_QC_FDR_global = FDR_global)

composite_sensitivity_preservation <- primary_join %>%
  left_join(binary_join, by = "composite_score") %>%
  left_join(cd45_join, by = "composite_score") %>%
  left_join(hladr_join, by = "composite_score") %>%
  left_join(technical_join, by = "composite_score") %>%
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
  arrange(primary_FDR_global, composite_score)

robust_composite_results <- composite_sensitivity_preservation %>%
  filter(direction_preserved_all_sensitivities == TRUE,
         FDR_global_preserved_all_sensitivities == TRUE) %>%
  arrange(primary_FDR_global)

composite_model_subject_summary <- tibble(
  n_primary_model_subjects = nrow(primary_df),
  n_primary_healthy = sum(primary_df$disease_group == "Healthy control", na.rm = TRUE),
  n_primary_cancer = sum(primary_df$disease_group == "Cancer patient", na.rm = TRUE),
  n_binary_sex_subjects = nrow(binary_sex_df),
  n_cd45_event_QC_subjects = nrow(cd45_event_qc_df),
  n_hladr_event_QC_subjects = nrow(hladr_event_qc_df),
  n_technical_QC_subjects = nrow(technical_qc_df),
  n_composite_scores = length(score_names),
  n_primary_FDR_lt_0p05 = sum(primary_composite_results$FDR_global < 0.05, na.rm = TRUE),
  n_direction_preserved_all_sensitivities = sum(composite_sensitivity_preservation$direction_preserved_all_sensitivities == TRUE, na.rm = TRUE),
  n_FDR_preserved_all_sensitivities = sum(composite_sensitivity_preservation$FDR_global_preserved_all_sensitivities == TRUE, na.rm = TRUE)
)

component_summary_by_composite <- composite_definitions_annotated %>%
  group_by(composite_score, component_axis) %>%
  summarise(
    n_components = n(),
    n_components_primary_global_FDR_lt_0p05 = sum(FDR_global < 0.05, na.rm = TRUE),
    n_components_direction_matches_composite =
      sum((component_direction == 1 & beta_cancer_vs_healthy > 0) |
            (component_direction == -1 & beta_cancer_vs_healthy < 0),
          na.rm = TRUE),
    min_component_FDR_global = min(FDR_global, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(composite_score)

score_availability_summary <- cp23_score_data %>%
  summarise(
    across(all_of(score_names), ~ sum(!is.na(.x)), .names = "n_nonmissing_{.col}")
  )

write_csv(cp23_score_data, file.path(out_dir, "SDY2583_CP23_score_data_STEP4.csv"))
write_csv(composite_definitions_annotated, file.path(out_dir, "SDY2583_CP23_composite_definitions_STEP4.csv"))
write_csv(missing_composite_components, file.path(out_dir, "SDY2583_CP23_missing_composite_components_STEP4.csv"))
write_csv(composite_score_dictionary, file.path(out_dir, "SDY2583_CP23_composite_score_dictionary_STEP4.csv"))
write_csv(primary_composite_results, file.path(out_dir, "SDY2583_CP23_primary_composite_results_STEP4.csv"))
write_csv(binary_sex_composite_results, file.path(out_dir, "SDY2583_CP23_binary_sex_composite_results_STEP4.csv"))
write_csv(cd45_event_qc_composite_results, file.path(out_dir, "SDY2583_CP23_cd45_event_QC_composite_results_STEP4.csv"))
write_csv(hladr_event_qc_composite_results, file.path(out_dir, "SDY2583_CP23_hladr_myeloid_event_QC_composite_results_STEP4.csv"))
write_csv(technical_qc_composite_results, file.path(out_dir, "SDY2583_CP23_technical_QC_composite_results_STEP4.csv"))
write_csv(all_composite_model_results, file.path(out_dir, "SDY2583_CP23_all_composite_model_results_STEP4.csv"))
write_csv(composite_sensitivity_preservation, file.path(out_dir, "SDY2583_CP23_composite_sensitivity_preservation_STEP4.csv"))
write_csv(robust_composite_results, file.path(out_dir, "SDY2583_CP23_robust_composite_results_STEP4.csv"))
write_csv(composite_model_subject_summary, file.path(out_dir, "SDY2583_CP23_composite_model_subject_summary_STEP4.csv"))
write_csv(component_summary_by_composite, file.path(out_dir, "SDY2583_CP23_component_summary_by_composite_STEP4.csv"))
write_csv(score_availability_summary, file.path(out_dir, "SDY2583_CP23_score_availability_summary_STEP4.csv"))

save(
  cp23_score_data,
  composite_definitions,
  composite_definitions_available,
  composite_definitions_annotated,
  missing_composite_components,
  composite_score_dictionary,
  score_names,
  primary_composite_results,
  binary_sex_composite_results,
  cd45_event_qc_composite_results,
  hladr_event_qc_composite_results,
  technical_qc_composite_results,
  all_composite_model_results,
  composite_sensitivity_preservation,
  robust_composite_results,
  composite_model_subject_summary,
  component_summary_by_composite,
  score_availability_summary,
  file = file.path(rdata_dir, "SDY2583_CP23_STEP4_composite_scores.RData")
)

cat("\n============================================================\n")
cat("SDY2583 CP23 STEP 4 COMPLETE: COMPOSITE SCORES\n")
cat("============================================================\n")

cat("\nComposite model subject summary:\n")
print(as.data.frame(composite_model_subject_summary), row.names = FALSE)

cat("\nComposite definitions:\n")
print(as.data.frame(composite_definitions_annotated), row.names = FALSE)

cat("\nComponent summary by composite:\n")
print(as.data.frame(component_summary_by_composite), row.names = FALSE)

cat("\nScore availability summary:\n")
print(as.data.frame(score_availability_summary), row.names = FALSE)

cat("\nPrimary composite results:\n")
print(as.data.frame(primary_composite_results), row.names = FALSE)

cat("\nComposite sensitivity preservation:\n")
print(as.data.frame(composite_sensitivity_preservation), row.names = FALSE)

cat("\nRobust composite results:\n")
print(as.data.frame(robust_composite_results), row.names = FALSE)

cat("\nMissing composite component features:\n")
print(as.data.frame(missing_composite_components), row.names = FALSE)

cat("\nOutputs saved in:\n")
print(out_dir)

cat("============================================================\n")
