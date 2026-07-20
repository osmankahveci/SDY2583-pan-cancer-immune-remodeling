# ============================================================
# SDY2583 CP10
# STEP 4 SAFE: Composite scores for myeloid / granulocytic remodeling
#
# Important interpretation:
# - CP10 composites are phenotype-based.
# - They do NOT claim definitive neutrophils, eosinophils, DCs, or MDSCs.
# - HLA-DR-related CP10 results are interpreted mainly as expansion of
#   CD14+ monocyte-like / HLA-DR+ and HLA-DR-low phenotypes within CD45+,
#   not as definitive per-monocyte HLA-DR suppression unless the relevant
#   within-CD14 or median-HLA-DR feature is robust.
#
# Output:
# outputs/CP10/05_composite_scores
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
out_dir <- file.path(analysis_dir, "05_composite_scores")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rdata_dir, recursive = TRUE, showWarnings = FALSE)

step3a_rdata <- file.path(rdata_dir, "SDY2583_CP10_STEP3A_metadata_merge_age_QC.RData")
step3b_rdata <- file.path(rdata_dir, "SDY2583_CP10_STEP3B_age_sex_adjusted_statistics.RData")

if (!file.exists(step3a_rdata)) stop("Step 3A RData bulunamadı: ", step3a_rdata)
if (!file.exists(step3b_rdata)) stop("Step 3B RData bulunamadı: ", step3b_rdata)

load(step3a_rdata)
load(step3b_rdata)

if (!exists("cp10_analysis_data")) stop("cp10_analysis_data bulunamadı.")
if (!exists("primary_results")) stop("primary_results bulunamadı.")

# ------------------------------------------------------------
# 3. Composite definitions
# ------------------------------------------------------------

# component_direction:
#   +1 = higher feature value increases the composite
#   -1 = lower feature value increases the composite after sign reversal
#
# The chosen directions make higher composite score reflect the
# cancer-associated remodeling direction observed in Step 3B.

composite_definitions <- tibble::tribble(
  ~composite_score, ~feature, ~component_direction, ~module, ~interpretation,

  # Granulocyte-like enrichment
  "CP10_granulocyte_like_enrichment_score", "pct_cd13_cd66b_gran_like_within_cd45",  1, "granulocyte_like_axis", "CD13+CD66b+ granulocyte-like enrichment",
  "CP10_granulocyte_like_enrichment_score", "pct_cd66b_gran_like_within_cd45",       1, "granulocyte_like_axis", "CD66b+ granulocyte-like enrichment",
  "CP10_granulocyte_like_enrichment_score", "pct_cd13_pos_within_cd66b_gran_like",  1, "granulocyte_like_axis", "CD13 enrichment within CD66b+ granulocyte-like cells",
  "CP10_granulocyte_like_enrichment_score", "median_CD13_in_CD45",                  1, "granulocyte_like_axis", "CD13 median intensity in CD45+ leukocytes",
  "CP10_granulocyte_like_enrichment_score", "median_CD66b_in_CD45",                 1, "granulocyte_like_axis", "CD66b median intensity in CD45+ leukocytes",
  "CP10_granulocyte_like_enrichment_score", "median_CD13_in_CD66b_gran_like",       1, "granulocyte_like_axis", "CD13 median intensity in CD66b+ granulocyte-like cells",
  "CP10_granulocyte_like_enrichment_score", "median_CD66b_in_CD66b_gran_like",      1, "granulocyte_like_axis", "CD66b median intensity in CD66b+ granulocyte-like cells",

  # Myeloid/granulocytic-to-lymphoid balance
  "CP10_myeloid_granulocytic_to_lymphoid_balance_score", "pct_myeloid_granulocytic_like_any_within_cd45",  1, "integrated_myeloid_granulocytic_balance", "Union myeloid/granulocytic-like fraction within CD45+ leukocytes",
  "CP10_myeloid_granulocytic_to_lymphoid_balance_score", "ratio_myeloid_granulocytic_to_lymphoid_like",    1, "integrated_myeloid_granulocytic_balance", "Myeloid/granulocytic-to-lymphoid-like ratio",
  "CP10_myeloid_granulocytic_to_lymphoid_balance_score", "ratio_cd66b_gran_like_to_t_like",                1, "integrated_myeloid_granulocytic_balance", "CD66b+ granulocyte-like to CD3+ T-like ratio",
  "CP10_myeloid_granulocytic_to_lymphoid_balance_score", "ratio_cd14_mono_like_to_t_like",                 1, "integrated_myeloid_granulocytic_balance", "CD14+ monocyte-like to CD3+ T-like ratio",
  "CP10_myeloid_granulocytic_to_lymphoid_balance_score", "pct_lymphoid_like_any_within_cd45",             -1, "integrated_myeloid_granulocytic_balance", "Lower CD3/CD19/CD56-defined lymphoid-like fraction",
  "CP10_myeloid_granulocytic_to_lymphoid_balance_score", "pct_t_like_within_cd45",                        -1, "integrated_myeloid_granulocytic_balance", "Lower CD3+ T-cell-like fraction",
  "CP10_myeloid_granulocytic_to_lymphoid_balance_score", "pct_b_like_within_cd45",                        -1, "integrated_myeloid_granulocytic_balance", "Lower CD19+ B-cell-like fraction",

  # Monocyte-like enrichment / HLA-DR phenotype redistribution
  "CP10_monocyte_like_HLA_DR_phenotype_remodeling_score", "pct_cd14_mono_like_within_cd45",                1, "monocyte_HLA_DR_axis", "CD14+ monocyte-like enrichment",
  "CP10_monocyte_like_HLA_DR_phenotype_remodeling_score", "median_CD14_in_CD45",                           1, "monocyte_HLA_DR_axis", "CD14 median intensity in CD45+ leukocytes",
  "CP10_monocyte_like_HLA_DR_phenotype_remodeling_score", "pct_cd14_hladr_low_within_cd45",                1, "monocyte_HLA_DR_axis", "CD14+HLA-DR-low monocyte-like enrichment within CD45+ leukocytes",
  "CP10_monocyte_like_HLA_DR_phenotype_remodeling_score", "pct_cd14_hladr_pos_within_cd45",                1, "monocyte_HLA_DR_axis", "CD14+HLA-DR+ monocyte-like enrichment within CD45+ leukocytes",
  "CP10_monocyte_like_HLA_DR_phenotype_remodeling_score", "pct_cd14_cd11c_hladr_pos_within_cd45",          1, "monocyte_HLA_DR_axis", "CD14+CD11c+HLA-DR+ myeloid APC-like enrichment",

  # APC-like myeloid enrichment
  "CP10_APC_like_myeloid_enrichment_score", "pct_cd11c_hladr_apc_like_within_cd45",        1, "APC_like_myeloid_axis", "CD11c+HLA-DR+ APC-like myeloid enrichment",
  "CP10_APC_like_myeloid_enrichment_score", "pct_cd14_cd11c_hladr_pos_within_cd45",        1, "APC_like_myeloid_axis", "CD14+CD11c+HLA-DR+ myeloid APC-like enrichment",
  "CP10_APC_like_myeloid_enrichment_score", "median_CD11c_in_CD45",                        1, "APC_like_myeloid_axis", "CD11c median intensity in CD45+ leukocytes",

  # CCR3/CD193 eosinophil-like granulocytic enrichment
  "CP10_CCR3_eosinophil_like_granulocytic_enrichment_score", "pct_ccr3_pos_within_cd45",                         1, "CCR3_eosinophil_like_axis", "CCR3/CD193+ fraction within CD45+ leukocytes",
  "CP10_CCR3_eosinophil_like_granulocytic_enrichment_score", "pct_ccr3_cd66b_eosinophil_like_within_cd45",        1, "CCR3_eosinophil_like_axis", "CCR3+CD66b+ eosinophil-like granulocytic enrichment",
  "CP10_CCR3_eosinophil_like_granulocytic_enrichment_score", "pct_ccr3_cd13_cd66b_gran_like_within_cd45",         1, "CCR3_eosinophil_like_axis", "CCR3+CD13+CD66b+ granulocyte-like enrichment",
  "CP10_CCR3_eosinophil_like_granulocytic_enrichment_score", "pct_ccr3_pos_within_cd66b_gran_like",               1, "CCR3_eosinophil_like_axis", "CCR3+ fraction within CD66b+ granulocyte-like cells",
  "CP10_CCR3_eosinophil_like_granulocytic_enrichment_score", "median_CCR3_in_CD45",                               1, "CCR3_eosinophil_like_axis", "CCR3 median intensity in CD45+ leukocytes",
  "CP10_CCR3_eosinophil_like_granulocytic_enrichment_score", "median_CCR3_in_CD66b_gran_like",                    1, "CCR3_eosinophil_like_axis", "CCR3 median intensity in CD66b+ granulocyte-like cells",

  # Integrated score
  "CP10_integrated_myeloid_granulocytic_remodeling_score", "pct_cd13_cd66b_gran_like_within_cd45",                1, "integrated", "Granulocyte-like enrichment",
  "CP10_integrated_myeloid_granulocytic_remodeling_score", "pct_cd66b_gran_like_within_cd45",                     1, "integrated", "CD66b+ granulocyte-like enrichment",
  "CP10_integrated_myeloid_granulocytic_remodeling_score", "pct_cd14_mono_like_within_cd45",                      1, "integrated", "CD14+ monocyte-like enrichment",
  "CP10_integrated_myeloid_granulocytic_remodeling_score", "pct_myeloid_granulocytic_like_any_within_cd45",       1, "integrated", "Union myeloid/granulocytic-like enrichment",
  "CP10_integrated_myeloid_granulocytic_remodeling_score", "ratio_myeloid_granulocytic_to_lymphoid_like",         1, "integrated", "Myeloid/granulocytic-to-lymphoid-like balance",
  "CP10_integrated_myeloid_granulocytic_remodeling_score", "ratio_cd66b_gran_like_to_t_like",                     1, "integrated", "Granulocyte-like to T-like balance",
  "CP10_integrated_myeloid_granulocytic_remodeling_score", "pct_cd11c_hladr_apc_like_within_cd45",                1, "integrated", "APC-like myeloid enrichment",
  "CP10_integrated_myeloid_granulocytic_remodeling_score", "pct_ccr3_cd66b_eosinophil_like_within_cd45",          1, "integrated", "CCR3+CD66b+ eosinophil-like granulocytic enrichment",
  "CP10_integrated_myeloid_granulocytic_remodeling_score", "pct_lymphoid_like_any_within_cd45",                  -1, "integrated", "Lower lymphoid-like fraction",
  "CP10_integrated_myeloid_granulocytic_remodeling_score", "pct_t_like_within_cd45",                             -1, "integrated", "Lower T-cell-like fraction"
)

# Keep only components present in data.
missing_components <- composite_definitions %>%
  filter(!(feature %in% names(cp10_analysis_data)))

if (nrow(missing_components) > 0) {
  cat("\nWarning: missing composite components will be dropped:\n")
  print(as.data.frame(missing_components), row.names = FALSE)
}

composite_definitions <- composite_definitions %>%
  filter(feature %in% names(cp10_analysis_data))

# ------------------------------------------------------------
# 4. Composite-score construction helper
# ------------------------------------------------------------

zscore_safe <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  s <- sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - m) / s
}

cp10_scores_data <- cp10_analysis_data

component_score_long <- list()

for (score_name in unique(composite_definitions$composite_score)) {

  defs <- composite_definitions %>%
    filter(composite_score == score_name)

  component_mat <- matrix(
    NA_real_,
    nrow = nrow(cp10_scores_data),
    ncol = nrow(defs)
  )

  colnames(component_mat) <- defs$feature

  for (j in seq_len(nrow(defs))) {
    ff <- defs$feature[j]
    dd <- defs$component_direction[j]
    component_mat[, j] <- zscore_safe(cp10_scores_data[[ff]]) * dd

    component_score_long[[length(component_score_long) + 1]] <- tibble(
      subject_id = cp10_scores_data$subject_id,
      composite_score = score_name,
      feature = ff,
      component_direction = dd,
      oriented_z = component_mat[, j]
    )
  }

  cp10_scores_data[[score_name]] <- rowMeans(component_mat, na.rm = TRUE)
}

component_score_long <- bind_rows(component_score_long)

composite_score_cols <- unique(composite_definitions$composite_score)

# ------------------------------------------------------------
# 5. Model helpers
# ------------------------------------------------------------

run_lm_score <- function(df, score_name, model_label) {

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

run_composite_model_set <- function(df, model_label) {
  bind_rows(lapply(composite_score_cols, function(sc) {
    run_lm_score(df, sc, model_label)
  })) %>%
    mutate(FDR_global = p.adjust(p_value, method = "BH")) %>%
    arrange(FDR_global, p_value)
}

# ------------------------------------------------------------
# 6. Run primary and sensitivity models
# ------------------------------------------------------------

primary_df <- cp10_scores_data %>%
  filter(model_ready_age_sex == TRUE)

binary_sex_df <- cp10_scores_data %>%
  filter(
    model_ready_age_sex == TRUE,
    !is.na(sex_binary)
  ) %>%
  mutate(sex = sex_binary)

event_qc_df <- cp10_scores_data %>%
  filter(
    model_ready_age_sex == TRUE,
    !is.na(n_cd45_viable),
    n_cd45_viable >= 1000
  )

technical_qc_df <- cp10_scores_data %>%
  filter(
    model_ready_age_sex == TRUE,
    channel_order_mismatch_file == FALSE,
    marker_mismatch_file == FALSE
  )

primary_composite_results <- run_composite_model_set(primary_df, "primary_age_sex_adjusted")
binary_sex_composite_results <- run_composite_model_set(binary_sex_df, "binary_sex_sensitivity")
event_qc_composite_results <- run_composite_model_set(event_qc_df, "CD45_event_QC_ge_1000")
technical_qc_composite_results <- run_composite_model_set(technical_qc_df, "technical_QC_exclude_mismatch_files")

all_composite_model_results <- bind_rows(
  primary_composite_results,
  binary_sex_composite_results,
  event_qc_composite_results,
  technical_qc_composite_results
)

# ------------------------------------------------------------
# 7. Sensitivity summary
# ------------------------------------------------------------

get_val <- function(df, score_name, model_label, col_name) {
  x <- df %>%
    filter(composite_score == score_name, model_label == !!model_label) %>%
    pull(!!rlang::sym(col_name))
  if (length(x) == 0) return(NA)
  x[1]
}

make_sensitivity_one <- function(score_name) {

  primary_beta <- as.numeric(get_val(all_composite_model_results, score_name, "primary_age_sex_adjusted", "beta_cancer_vs_healthy"))
  primary_direction <- as.character(get_val(all_composite_model_results, score_name, "primary_age_sex_adjusted", "direction"))
  primary_FDR <- as.numeric(get_val(all_composite_model_results, score_name, "primary_age_sex_adjusted", "FDR_global"))

  binary_beta <- as.numeric(get_val(all_composite_model_results, score_name, "binary_sex_sensitivity", "beta_cancer_vs_healthy"))
  binary_direction <- as.character(get_val(all_composite_model_results, score_name, "binary_sex_sensitivity", "direction"))
  binary_FDR <- as.numeric(get_val(all_composite_model_results, score_name, "binary_sex_sensitivity", "FDR_global"))

  event_beta <- as.numeric(get_val(all_composite_model_results, score_name, "CD45_event_QC_ge_1000", "beta_cancer_vs_healthy"))
  event_direction <- as.character(get_val(all_composite_model_results, score_name, "CD45_event_QC_ge_1000", "direction"))
  event_FDR <- as.numeric(get_val(all_composite_model_results, score_name, "CD45_event_QC_ge_1000", "FDR_global"))

  tech_beta <- as.numeric(get_val(all_composite_model_results, score_name, "technical_QC_exclude_mismatch_files", "beta_cancer_vs_healthy"))
  tech_direction <- as.character(get_val(all_composite_model_results, score_name, "technical_QC_exclude_mismatch_files", "direction"))
  tech_FDR <- as.numeric(get_val(all_composite_model_results, score_name, "technical_QC_exclude_mismatch_files", "FDR_global"))

  directions <- c(primary_direction, binary_direction, event_direction, tech_direction)
  direction_preserved <- all(!is.na(directions)) && length(unique(directions)) == 1

  tibble(
    composite_score = score_name,
    n_components = nrow(composite_definitions %>% filter(composite_score == score_name)),
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

composite_sensitivity_summary <- bind_rows(lapply(composite_score_cols, make_sensitivity_one)) %>%
  arrange(primary_FDR_global)

# ------------------------------------------------------------
# 8. Component-level summary
# ------------------------------------------------------------

component_level_primary_results <- composite_definitions %>%
  left_join(
    primary_results %>%
      select(
        feature,
        beta_cancer_vs_healthy,
        p_value,
        FDR_global,
        FDR_within_module,
        direction
      ),
    by = "feature"
  ) %>%
  arrange(composite_score, FDR_global)

component_summary_by_composite <- component_level_primary_results %>%
  group_by(composite_score) %>%
  summarise(
    n_components = n(),
    n_components_primary_global_FDR_lt_0p05 = sum(FDR_global < 0.05, na.rm = TRUE),
    n_components_direction_matches_composite =
      sum(
        (component_direction == 1 & direction == "higher_in_cancer") |
          (component_direction == -1 & direction == "lower_in_cancer"),
        na.rm = TRUE
      ),
    min_component_FDR_global = min(FDR_global, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(min_component_FDR_global)

composite_model_subject_summary <- tibble(
  n_primary_model_subjects = nrow(primary_df),
  n_primary_healthy = sum(primary_df$disease_group == "Healthy control"),
  n_primary_cancer = sum(primary_df$disease_group == "Cancer patient"),
  n_binary_sex_subjects = nrow(binary_sex_df),
  n_event_QC_subjects = nrow(event_qc_df),
  n_technical_QC_subjects = nrow(technical_qc_df),
  n_composite_scores = length(composite_score_cols),
  n_primary_FDR_lt_0p05 = sum(primary_composite_results$FDR_global < 0.05, na.rm = TRUE),
  n_direction_preserved_all_sensitivities = sum(composite_sensitivity_summary$direction_preserved_all_sensitivities == TRUE, na.rm = TRUE),
  n_FDR_preserved_all_sensitivities = sum(composite_sensitivity_summary$FDR_global_preserved_all_sensitivities == TRUE, na.rm = TRUE)
)

# ------------------------------------------------------------
# 9. Save outputs
# ------------------------------------------------------------

write_csv(cp10_scores_data, file.path(out_dir, "SDY2583_CP10_composite_scores_data_STEP4.csv"))
write_csv(composite_definitions, file.path(out_dir, "SDY2583_CP10_composite_definitions_STEP4.csv"))
write_csv(component_score_long, file.path(out_dir, "SDY2583_CP10_composite_component_oriented_z_STEP4.csv"))
write_csv(primary_composite_results, file.path(out_dir, "SDY2583_CP10_primary_composite_results_STEP4.csv"))
write_csv(binary_sex_composite_results, file.path(out_dir, "SDY2583_CP10_binary_sex_composite_results_STEP4.csv"))
write_csv(event_qc_composite_results, file.path(out_dir, "SDY2583_CP10_CD45_event_QC_composite_results_STEP4.csv"))
write_csv(technical_qc_composite_results, file.path(out_dir, "SDY2583_CP10_technical_QC_composite_results_STEP4.csv"))
write_csv(all_composite_model_results, file.path(out_dir, "SDY2583_CP10_all_composite_model_results_STEP4.csv"))
write_csv(composite_sensitivity_summary, file.path(out_dir, "SDY2583_CP10_composite_sensitivity_summary_STEP4.csv"))
write_csv(component_level_primary_results, file.path(out_dir, "SDY2583_CP10_component_level_primary_results_STEP4.csv"))
write_csv(component_summary_by_composite, file.path(out_dir, "SDY2583_CP10_component_summary_by_composite_STEP4.csv"))
write_csv(composite_model_subject_summary, file.path(out_dir, "SDY2583_CP10_composite_model_subject_summary_STEP4.csv"))

save(
  cp10_scores_data,
  composite_score_cols,
  composite_definitions,
  component_score_long,
  primary_composite_results,
  binary_sex_composite_results,
  event_qc_composite_results,
  technical_qc_composite_results,
  all_composite_model_results,
  composite_sensitivity_summary,
  component_level_primary_results,
  component_summary_by_composite,
  composite_model_subject_summary,
  file = file.path(rdata_dir, "SDY2583_CP10_STEP4_composite_scores.RData")
)

# Save lean object for integrated clinical matrix later.
cp10_scores_for_integration <- cp10_scores_data %>%
  select(subject_id, all_of(composite_score_cols))

save(
  cp10_scores_for_integration,
  composite_score_cols,
  composite_definitions,
  file = file.path(rdata_dir, "SDY2583_CP10_STEP4_composite_scores_for_integration.RData")
)

# ------------------------------------------------------------
# 10. Console output
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SDY2583 CP10 STEP 4 COMPLETE: COMPOSITE SCORES\n")
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

cat("\nOutputs saved in:\n")
print(out_dir)

cat("============================================================\n")
