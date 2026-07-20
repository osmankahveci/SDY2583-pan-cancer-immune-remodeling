# ============================================================
# SDY2583 CP16
# STEP 4 SAFE: Composite scores for APC/DC-like myeloid remodeling
#
# Input:
#   outputs/CP16/11_RData/
#     SDY2583_CP16_STEP3B_age_sex_adjusted_statistics.RData
#
# Composite scores:
#   1) CP16_CD1c_cDC2_like_depletion_score
#   2) CP16_CD123_pDC_like_depletion_score
#   3) CP16_CD141_CLEC9A_cDC1_like_depletion_score
#   4) CP16_monocyte_HLA_DR_low_remodeling_score
#   5) CP16_APC_core_attenuation_score
#   6) CP16_integrated_APC_DC_myeloid_remodeling_score
#
# Higher composite score = stronger cancer-associated remodeling direction.
#
# CP16-specific interpretation:
#   - Use phenotype-like terminology only.
#   - Do not claim definitive cDC1, cDC2, pDC, MDSC, or mature DC identity.
# ============================================================

rm(list = ls())
source(file.path(Sys.getenv("SDY2583_REPO_ROOT", unset = "."), "R", "shared", "bootstrap.R"))

# ------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------

cran_pkgs <- c("dplyr", "readr", "stringr", "tibble", "broom", "tidyr", "purrr")

for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(broom)
  library(tidyr)
  library(purrr)
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
out_dir <- file.path(analysis_dir, "05_composite_scores")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rdata_dir, recursive = TRUE, showWarnings = FALSE)

step3b_rdata <- file.path(rdata_dir, "SDY2583_CP16_STEP3B_age_sex_adjusted_statistics.RData")

if (!file.exists(step3b_rdata)) {
  stop("Step 3B RData bulunamadı: ", step3b_rdata)
}

load(step3b_rdata)

# Reset paths after RData load.
analysis_dir <- sd_analysis_dir("CP16")
rdata_dir <- file.path(analysis_dir, "11_RData")
out_dir <- file.path(analysis_dir, "05_composite_scores")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!exists("cp16_analysis_data")) stop("cp16_analysis_data bulunamadı.")
if (!exists("primary_results")) stop("primary_results bulunamadı.")
if (!exists("feature_dictionary_full")) stop("feature_dictionary_full bulunamadı.")

# ------------------------------------------------------------
# 3. Composite definitions
# ------------------------------------------------------------
# component_direction:
#   +1 = feature is expected to be higher in cancer
#   -1 = feature is expected to be lower in cancer, so depletion is scored higher

composite_definitions_raw <- tibble::tribble(
  ~composite_score, ~feature, ~component_direction, ~component_axis,

  # ----------------------------------------------------------
  # CD1c/cDC2-like depletion
  # ----------------------------------------------------------
  "CP16_CD1c_cDC2_like_depletion_score", "pct_cd1c_fceri_apc_like_within_cd45_dump_low", -1, "CD1c+FcERI+ APC-like depletion",
  "CP16_CD1c_cDC2_like_depletion_score", "pct_cd1c_apc_like_within_cd45_dump_low", -1, "CD1c+ APC-like depletion",
  "CP16_CD1c_cDC2_like_depletion_score", "pct_cd1c_cd14low_cdc2_like_within_cd45_dump_low", -1, "CD1c+CD14-low cDC2-like depletion",
  "CP16_CD1c_cDC2_like_depletion_score", "pct_cd1c_apc_like_within_hladr_core", -1, "CD1c+ within HLA-DR+ APC-like core depletion",
  "CP16_CD1c_cDC2_like_depletion_score", "pct_cd1c_cd14low_cdc2_like_within_hladr_core", -1, "CD1c+CD14-low within HLA-DR+ core depletion",
  "CP16_CD1c_cDC2_like_depletion_score", "pct_cd1c_fceri_apc_like_within_cd1c_apc_like", -1, "FcERI+ fraction within CD1c+ APC-like depletion",
  "CP16_CD1c_cDC2_like_depletion_score", "pct_fceri_pos_hladr_within_cd45_dump_low", -1, "FcERI+HLA-DR+ APC-like depletion",

  # ----------------------------------------------------------
  # CD123/pDC-like depletion
  # ----------------------------------------------------------
  "CP16_CD123_pDC_like_depletion_score", "pct_cd123_hladr_dc_like_within_cd45_dump_low", -1, "CD123+HLA-DR+ DC-like depletion",
  "CP16_CD123_pDC_like_depletion_score", "pct_cd123_hladr_dc_like_within_hladr_core", -1, "CD123+ fraction within HLA-DR+ core depletion",
  "CP16_CD123_pDC_like_depletion_score", "pct_cd123_pdc_like_within_cd45_dump_low", -1, "CD123+ pDC-like depletion",
  "CP16_CD123_pDC_like_depletion_score", "pct_cd123_pdc_like_within_hladr_core", -1, "pDC-like fraction within HLA-DR+ core depletion",
  "CP16_CD123_pDC_like_depletion_score", "pct_cd123_cd11c_mixed_apc_like_within_cd45_dump_low", -1, "CD123+CD11c+ mixed APC-like depletion",

  # ----------------------------------------------------------
  # CD141/CLEC9A cDC1-like depletion
  # ----------------------------------------------------------
  "CP16_CD141_CLEC9A_cDC1_like_depletion_score", "pct_cd141_pos_hladr_within_cd45_dump_low", -1, "CD141+HLA-DR+ depletion",
  "CP16_CD141_CLEC9A_cDC1_like_depletion_score", "pct_cd141_clec9a_cdc1_like_within_cd45_dump_low", -1, "CD141+CLEC9A+ cDC1-like depletion",
  "CP16_CD141_CLEC9A_cDC1_like_depletion_score", "pct_cd141_clec9a_cdc1_like_within_hladr_core", -1, "CD141+CLEC9A+ fraction within HLA-DR+ core depletion",

  # ----------------------------------------------------------
  # Monocyte-like HLA-DR-low remodeling
  # ----------------------------------------------------------
  "CP16_monocyte_HLA_DR_low_remodeling_score", "pct_cd14_hladr_low_within_cd45_dump_low", 1, "CD14+HLA-DR-low monocyte-like enrichment",
  "CP16_monocyte_HLA_DR_low_remodeling_score", "pct_hladr_low_within_cd14_mono_like", 1, "HLA-DR-low fraction within CD14+ monocyte-like enrichment",
  "CP16_monocyte_HLA_DR_low_remodeling_score", "pct_hladr_pos_within_cd14_mono_like", -1, "HLA-DR+ fraction within CD14+ monocyte-like depletion",
  "CP16_monocyte_HLA_DR_low_remodeling_score", "ratio_cd11c_hladr_apc_to_cd14_mono_like", -1, "lower CD11c+HLA-DR+ APC-like to CD14+ monocyte-like ratio",
  "CP16_monocyte_HLA_DR_low_remodeling_score", "ratio_cd14_mono_like_to_cd123_pdc_like", 1, "higher CD14+ monocyte-like to CD123+ pDC-like ratio",
  "CP16_monocyte_HLA_DR_low_remodeling_score", "pct_cd14_mono_like_within_cd45_dump_low", 1, "CD14+ monocyte-like enrichment",

  # ----------------------------------------------------------
  # APC core attenuation
  # ----------------------------------------------------------
  "CP16_APC_core_attenuation_score", "pct_hladr_apc_core_within_cd45_dump_low", -1, "HLA-DR+ APC-like core attenuation",
  "CP16_APC_core_attenuation_score", "pct_cd11c_hladr_apc_like_within_cd45_dump_low", -1, "CD11c+HLA-DR+ APC-like attenuation",
  "CP16_APC_core_attenuation_score", "pct_cd11c_hladr_apc_like_within_hladr_core", -1, "CD11c+ fraction within HLA-DR+ core attenuation",
  "CP16_APC_core_attenuation_score", "pct_cd13_cd11c_hladr_apc_like_within_cd45_dump_low", -1, "CD13+CD11c+HLA-DR+ APC-like attenuation",

  # ----------------------------------------------------------
  # Integrated APC/DC-like myeloid remodeling
  # ----------------------------------------------------------
  "CP16_integrated_APC_DC_myeloid_remodeling_score", "pct_cd1c_fceri_apc_like_within_cd45_dump_low", -1, "CD1c+FcERI+ APC-like depletion",
  "CP16_integrated_APC_DC_myeloid_remodeling_score", "pct_cd1c_apc_like_within_cd45_dump_low", -1, "CD1c+ APC-like depletion",
  "CP16_integrated_APC_DC_myeloid_remodeling_score", "pct_cd1c_cd14low_cdc2_like_within_cd45_dump_low", -1, "CD1c+CD14-low cDC2-like depletion",
  "CP16_integrated_APC_DC_myeloid_remodeling_score", "pct_cd123_hladr_dc_like_within_cd45_dump_low", -1, "CD123+HLA-DR+ DC-like depletion",
  "CP16_integrated_APC_DC_myeloid_remodeling_score", "pct_cd123_pdc_like_within_cd45_dump_low", -1, "CD123+ pDC-like depletion",
  "CP16_integrated_APC_DC_myeloid_remodeling_score", "pct_cd141_clec9a_cdc1_like_within_cd45_dump_low", -1, "CD141+CLEC9A+ cDC1-like depletion",
  "CP16_integrated_APC_DC_myeloid_remodeling_score", "pct_cd141_pos_hladr_within_cd45_dump_low", -1, "CD141+HLA-DR+ depletion",
  "CP16_integrated_APC_DC_myeloid_remodeling_score", "pct_cd14_hladr_low_within_cd45_dump_low", 1, "CD14+HLA-DR-low monocyte-like enrichment",
  "CP16_integrated_APC_DC_myeloid_remodeling_score", "pct_hladr_low_within_cd14_mono_like", 1, "HLA-DR-low within CD14+ monocyte-like enrichment",
  "CP16_integrated_APC_DC_myeloid_remodeling_score", "pct_hladr_pos_within_cd14_mono_like", -1, "HLA-DR+ within CD14+ monocyte-like depletion",
  "CP16_integrated_APC_DC_myeloid_remodeling_score", "pct_hladr_apc_core_within_cd45_dump_low", -1, "HLA-DR+ APC-like core attenuation",
  "CP16_integrated_APC_DC_myeloid_remodeling_score", "ratio_cd11c_hladr_apc_to_cd14_mono_like", -1, "lower APC-like to monocyte-like balance"
)

# Keep only features that exist in analysis data.
composite_definitions <- composite_definitions_raw %>%
  filter(feature %in% names(cp16_analysis_data)) %>%
  left_join(
    feature_dictionary_full %>%
      select(feature, module, interpretation),
    by = "feature"
  ) %>%
  arrange(composite_score, feature)

missing_component_features <- composite_definitions_raw %>%
  filter(!(feature %in% names(cp16_analysis_data)))

# ------------------------------------------------------------
# 4. Build composite scores
# ------------------------------------------------------------

zscore_safe <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  s <- sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - m) / s
}

cp16_scores_data <- cp16_analysis_data

score_cols <- unique(composite_definitions$composite_score)

component_availability_rows <- list()

for (sc in score_cols) {

  defs <- composite_definitions %>%
    filter(composite_score == sc)

  component_mat <- matrix(NA_real_, nrow = nrow(cp16_scores_data), ncol = nrow(defs))
  colnames(component_mat) <- defs$feature

  for (j in seq_len(nrow(defs))) {
    ff <- defs$feature[j]
    dd <- defs$component_direction[j]
    component_mat[, j] <- zscore_safe(cp16_scores_data[[ff]]) * dd
  }

  n_available_components <- rowSums(!is.na(component_mat))
  min_required <- ceiling(nrow(defs) / 2)
  score_value <- rowMeans(component_mat, na.rm = TRUE)
  score_value[n_available_components < min_required] <- NA_real_

  cp16_scores_data[[sc]] <- score_value
  cp16_scores_data[[paste0(sc, "_n_components_available")]] <- n_available_components

  component_availability_rows[[sc]] <- tibble(
    composite_score = sc,
    n_defined_components = nrow(defs),
    min_required_components = min_required,
    n_subjects_with_score = sum(!is.na(score_value)),
    n_subjects_missing_score = sum(is.na(score_value)),
    median_components_available = median(n_available_components, na.rm = TRUE),
    min_components_available = min(n_available_components, na.rm = TRUE),
    max_components_available = max(n_available_components, na.rm = TRUE)
  )
}

component_availability_summary <- bind_rows(component_availability_rows)

composite_score_cols <- score_cols

# ------------------------------------------------------------
# 5. Model helper
# ------------------------------------------------------------

run_score_model <- function(df, score_name, model_label) {

  model_df <- df %>%
    transmute(
      value = suppressWarnings(as.numeric(.data[[score_name]])),
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

  if (n_model < 50 || n_healthy < 20 || n_cancer < 20) return(NULL)

  fit <- tryCatch(
    lm(value ~ disease_group + age_for_model + sex, data = model_df),
    error = function(e) NULL
  )

  if (is.null(fit)) return(NULL)

  tt <- tryCatch(broom::tidy(fit, conf.int = TRUE), error = function(e) NULL)
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

  tibble(
    model_label = model_label,
    composite_score = score_name,
    n_model = n_model,
    n_healthy = n_healthy,
    n_cancer = n_cancer,
    healthy_mean = desc$mean[desc$disease_group == "Healthy control"][1],
    cancer_mean = desc$mean[desc$disease_group == "Cancer patient"][1],
    healthy_median = desc$median[desc$disease_group == "Healthy control"][1],
    cancer_median = desc$median[desc$disease_group == "Cancer patient"][1],
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

run_score_model_set <- function(df, model_label) {
  bind_rows(lapply(composite_score_cols, function(sc) {
    run_score_model(df, sc, model_label)
  })) %>%
    mutate(FDR_global = p.adjust(p_value, method = "BH")) %>%
    arrange(FDR_global, p_value)
}

# ------------------------------------------------------------
# 6. Model datasets
# ------------------------------------------------------------

primary_df <- cp16_scores_data %>%
  filter(model_ready_age_sex == TRUE)

binary_sex_df <- cp16_scores_data %>%
  filter(
    model_ready_age_sex == TRUE,
    !is.na(sex_binary)
  ) %>%
  mutate(sex = sex_binary)

cd45_event_qc_df <- cp16_scores_data %>%
  filter(
    model_ready_age_sex == TRUE,
    !is.na(n_cd45_dump_low),
    n_cd45_dump_low >= 300
  )

hladr_event_qc_df <- cp16_scores_data %>%
  filter(
    model_ready_age_sex == TRUE,
    !is.na(n_hladr_apc_core),
    n_hladr_apc_core >= 300
  )

technical_qc_df <- cp16_scores_data %>%
  filter(
    model_ready_age_sex == TRUE,
    channel_order_mismatch_file == FALSE,
    marker_order_mismatch_file == FALSE,
    marker_set_mismatch_file == FALSE
  )

# ------------------------------------------------------------
# 7. Run models
# ------------------------------------------------------------

primary_composite_results <- run_score_model_set(primary_df, "primary_age_sex_adjusted")
binary_sex_composite_results <- run_score_model_set(binary_sex_df, "binary_sex_sensitivity")
cd45_event_qc_composite_results <- run_score_model_set(cd45_event_qc_df, "CD45_dump_low_event_QC_ge_300")
hladr_event_qc_composite_results <- run_score_model_set(hladr_event_qc_df, "HLADR_APC_core_event_QC_ge_300")
technical_qc_composite_results <- run_score_model_set(technical_qc_df, "technical_QC_exclude_mismatch_files")

all_composite_model_results <- bind_rows(
  primary_composite_results,
  binary_sex_composite_results,
  cd45_event_qc_composite_results,
  hladr_event_qc_composite_results,
  technical_qc_composite_results
)

# ------------------------------------------------------------
# 8. Composite sensitivity summary
# ------------------------------------------------------------

get_val <- function(df, score_name, model_label, col_name) {
  x <- df %>%
    filter(composite_score == score_name, model_label == !!model_label) %>%
    pull(!!rlang::sym(col_name))
  if (length(x) == 0) return(NA)
  x[1]
}

make_composite_sensitivity <- function(score_name) {

  primary_beta <- as.numeric(get_val(all_composite_model_results, score_name, "primary_age_sex_adjusted", "beta_cancer_vs_healthy"))
  primary_direction <- as.character(get_val(all_composite_model_results, score_name, "primary_age_sex_adjusted", "direction"))
  primary_FDR <- as.numeric(get_val(all_composite_model_results, score_name, "primary_age_sex_adjusted", "FDR_global"))

  binary_beta <- as.numeric(get_val(all_composite_model_results, score_name, "binary_sex_sensitivity", "beta_cancer_vs_healthy"))
  binary_direction <- as.character(get_val(all_composite_model_results, score_name, "binary_sex_sensitivity", "direction"))
  binary_FDR <- as.numeric(get_val(all_composite_model_results, score_name, "binary_sex_sensitivity", "FDR_global"))

  cd45_beta <- as.numeric(get_val(all_composite_model_results, score_name, "CD45_dump_low_event_QC_ge_300", "beta_cancer_vs_healthy"))
  cd45_direction <- as.character(get_val(all_composite_model_results, score_name, "CD45_dump_low_event_QC_ge_300", "direction"))
  cd45_FDR <- as.numeric(get_val(all_composite_model_results, score_name, "CD45_dump_low_event_QC_ge_300", "FDR_global"))

  hladr_beta <- as.numeric(get_val(all_composite_model_results, score_name, "HLADR_APC_core_event_QC_ge_300", "beta_cancer_vs_healthy"))
  hladr_direction <- as.character(get_val(all_composite_model_results, score_name, "HLADR_APC_core_event_QC_ge_300", "direction"))
  hladr_FDR <- as.numeric(get_val(all_composite_model_results, score_name, "HLADR_APC_core_event_QC_ge_300", "FDR_global"))

  tech_beta <- as.numeric(get_val(all_composite_model_results, score_name, "technical_QC_exclude_mismatch_files", "beta_cancer_vs_healthy"))
  tech_direction <- as.character(get_val(all_composite_model_results, score_name, "technical_QC_exclude_mismatch_files", "direction"))
  tech_FDR <- as.numeric(get_val(all_composite_model_results, score_name, "technical_QC_exclude_mismatch_files", "FDR_global"))

  directions <- c(primary_direction, binary_direction, cd45_direction, hladr_direction, tech_direction)
  direction_preserved <- all(!is.na(directions)) && length(unique(directions)) == 1

  tibble(
    composite_score = score_name,
    n_components = sum(composite_definitions$composite_score == score_name),
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

composite_sensitivity_summary <- bind_rows(lapply(composite_score_cols, make_composite_sensitivity)) %>%
  arrange(primary_FDR_global)

# ------------------------------------------------------------
# 9. Summaries
# ------------------------------------------------------------

component_summary_by_composite <- composite_definitions %>%
  group_by(composite_score) %>%
  summarise(
    n_components = n(),
    n_components_primary_global_FDR_lt_0p05 =
      sum(primary_results$FDR_global[match(feature, primary_results$feature)] < 0.05, na.rm = TRUE),
    n_components_direction_matches_composite =
      sum(
        case_when(
          component_direction == 1 ~ primary_results$direction[match(feature, primary_results$feature)] == "higher_in_cancer",
          component_direction == -1 ~ primary_results$direction[match(feature, primary_results$feature)] == "lower_in_cancer",
          TRUE ~ FALSE
        ),
        na.rm = TRUE
      ),
    min_component_FDR_global =
      min(primary_results$FDR_global[match(feature, primary_results$feature)], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(composite_score)

composite_model_subject_summary <- tibble(
  n_primary_model_subjects = nrow(primary_df),
  n_primary_healthy = sum(primary_df$disease_group == "Healthy control"),
  n_primary_cancer = sum(primary_df$disease_group == "Cancer patient"),
  n_binary_sex_subjects = nrow(binary_sex_df),
  n_cd45_event_QC_subjects = nrow(cd45_event_qc_df),
  n_hladr_event_QC_subjects = nrow(hladr_event_qc_df),
  n_technical_QC_subjects = nrow(technical_qc_df),
  n_composite_scores = length(composite_score_cols),
  n_primary_FDR_lt_0p05 = sum(primary_composite_results$FDR_global < 0.05, na.rm = TRUE),
  n_direction_preserved_all_sensitivities = sum(composite_sensitivity_summary$direction_preserved_all_sensitivities == TRUE, na.rm = TRUE),
  n_FDR_preserved_all_sensitivities = sum(composite_sensitivity_summary$FDR_global_preserved_all_sensitivities == TRUE, na.rm = TRUE)
)

# Integration object for later ALL9 matrix update
cp16_scores_for_integration <- cp16_scores_data %>%
  select(subject_id, all_of(composite_score_cols))

# ------------------------------------------------------------
# 10. Save outputs
# ------------------------------------------------------------

write_csv(composite_definitions, file.path(out_dir, "SDY2583_CP16_composite_definitions_STEP4.csv"))
write_csv(missing_component_features, file.path(out_dir, "SDY2583_CP16_missing_composite_component_features_STEP4.csv"))
write_csv(component_availability_summary, file.path(out_dir, "SDY2583_CP16_component_availability_summary_STEP4.csv"))
write_csv(component_summary_by_composite, file.path(out_dir, "SDY2583_CP16_component_summary_by_composite_STEP4.csv"))
write_csv(cp16_scores_data, file.path(out_dir, "SDY2583_CP16_scores_data_STEP4.csv"))
write_csv(primary_composite_results, file.path(out_dir, "SDY2583_CP16_primary_composite_results_STEP4.csv"))
write_csv(binary_sex_composite_results, file.path(out_dir, "SDY2583_CP16_binary_sex_composite_results_STEP4.csv"))
write_csv(cd45_event_qc_composite_results, file.path(out_dir, "SDY2583_CP16_CD45_event_QC_composite_results_STEP4.csv"))
write_csv(hladr_event_qc_composite_results, file.path(out_dir, "SDY2583_CP16_HLADR_event_QC_composite_results_STEP4.csv"))
write_csv(technical_qc_composite_results, file.path(out_dir, "SDY2583_CP16_technical_QC_composite_results_STEP4.csv"))
write_csv(all_composite_model_results, file.path(out_dir, "SDY2583_CP16_all_composite_model_results_STEP4.csv"))
write_csv(composite_sensitivity_summary, file.path(out_dir, "SDY2583_CP16_composite_sensitivity_summary_STEP4.csv"))
write_csv(composite_model_subject_summary, file.path(out_dir, "SDY2583_CP16_composite_model_subject_summary_STEP4.csv"))
write_csv(cp16_scores_for_integration, file.path(out_dir, "SDY2583_CP16_scores_for_integration_STEP4.csv"))

save(
  cp16_scores_data,
  cp16_scores_for_integration,
  composite_score_cols,
  composite_definitions,
  composite_definitions_raw,
  missing_component_features,
  component_availability_summary,
  component_summary_by_composite,
  primary_composite_results,
  binary_sex_composite_results,
  cd45_event_qc_composite_results,
  hladr_event_qc_composite_results,
  technical_qc_composite_results,
  all_composite_model_results,
  composite_sensitivity_summary,
  composite_model_subject_summary,
  file = file.path(rdata_dir, "SDY2583_CP16_STEP4_composite_scores.RData")
)

save(
  cp16_scores_for_integration,
  composite_score_cols,
  file = file.path(rdata_dir, "SDY2583_CP16_STEP4_composite_scores_for_integration.RData")
)

# ------------------------------------------------------------
# 11. Console output
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SDY2583 CP16 STEP 4 COMPLETE: COMPOSITE SCORES\n")
cat("============================================================\n")

cat("\nComposite model subject summary:\n")
print(as.data.frame(composite_model_subject_summary), row.names = FALSE)

cat("\nComposite definitions:\n")
print(as.data.frame(composite_definitions), row.names = FALSE)

cat("\nComponent summary by composite:\n")
print(as.data.frame(component_summary_by_composite), row.names = FALSE)

cat("\nPrimary composite results:\n")
print(as.data.frame(primary_composite_results), row.names = FALSE)

cat("\nComposite sensitivity summary:\n")
print(as.data.frame(composite_sensitivity_summary), row.names = FALSE)

cat("\nMissing composite component features:\n")
print(as.data.frame(missing_component_features), row.names = FALSE)

cat("\nOutputs saved in:\n")
print(out_dir)

cat("============================================================\n")
