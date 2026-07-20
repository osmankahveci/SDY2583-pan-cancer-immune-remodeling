# ============================================================
# SDY2583 CP10
# STEP 7A SAFE: Main manuscript-ready figures
#
# Figures produced:
#   Figure 1: CP10 composite scores, cancer vs healthy
#   Figure 2: Key CP10 feature-level axes, cancer vs healthy
#   Figure 3: CP10 feature-level forest plot
#   Figure 4: CP10 composite threshold sensitivity
#   Figure 5: CP10 age-matched composite effects
#   Figure 6: CP10 cancer-subgroup clinical annotation, breast vs rest
#
# Notes:
# - These are statistical/main-result figures, not FlowJo-style FCS plots.
# - FlowJo-style representative gating figures will be Step 7B.
# - Phenotype-like terminology is used throughout.
# ============================================================

rm(list = ls())
source(file.path(Sys.getenv("SDY2583_REPO_ROOT", unset = "."), "R", "shared", "bootstrap.R"))

# ------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------

cran_pkgs <- c(
  "dplyr", "readr", "stringr", "tibble", "tidyr",
  "ggplot2", "forcats", "scales", "patchwork"
)

for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
  library(forcats)
  library(scales)
  library(patchwork)
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

cp10_analysis_dir <- sd_analysis_dir("CP10")
cp10_rdata_dir <- file.path(cp10_analysis_dir, "11_RData")

integrated_matrix_dir <- sd_integrated_dir()
clinical_rdata_dir <- file.path(integrated_matrix_dir, "RData")

out_dir <- file.path(cp10_analysis_dir, "08_main_figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

step3b_rdata <- file.path(cp10_rdata_dir, "SDY2583_CP10_STEP3B_age_sex_adjusted_statistics.RData")
step4_rdata <- file.path(cp10_rdata_dir, "SDY2583_CP10_STEP4_composite_scores.RData")
step4b_rdata <- file.path(cp10_rdata_dir, "SDY2583_CP10_STEP4B_age_sensitivity_caliper_interaction.RData")
step5_rdata <- file.path(cp10_rdata_dir, "SDY2583_CP10_STEP5_threshold_sensitivity_TARGETED.RData")
step6_rdata <- file.path(clinical_rdata_dir, "SDY2583_CP10_clinical_annotation_STEP6.RData")

required_files <- c(step3b_rdata, step4_rdata, step4b_rdata, step5_rdata, step6_rdata)
missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop("Eksik RData dosyaları var:\n", paste(missing_files, collapse = "\n"))
}

load(step3b_rdata)
load(step4_rdata)
load(step4b_rdata)
load(step5_rdata)
load(step6_rdata)

# Required objects check
# ------------------------------------------------------------
# 2B. Backward-compatible fallback for Step 6 object names
# ------------------------------------------------------------
# Some Step 6 runs saved the cancer-subgroup one-vs-rest table only as CSV
# or with a different in-memory name. This block prevents Step 7A from
# failing when the RData object is absent.

step6_out_dir <- file.path(
  integrated_matrix_dir,
  "06_CP10_clinical_annotation"
)

if (!exists("cancer_subgroup_one_vs_rest_results")) {
  fallback_csv <- file.path(
    step6_out_dir,
    "SDY2583_CP10_cancer_subgroup_one_vs_rest_results_STEP6.csv"
  )

  if (file.exists(fallback_csv)) {
    cancer_subgroup_one_vs_rest_results <- readr::read_csv(
      fallback_csv,
      show_col_types = FALSE
    )
  }
}

if (!exists("cancer_subgroup_one_vs_rest_results")) {
  # Allow the main figures to run without Figure 6 if clinical subgroup table
  # is unavailable. Figure 1-5 do not require this table.
  cancer_subgroup_one_vs_rest_results <- tibble(
    score = character(),
    cancer_subgroup = character(),
    n_model = integer(),
    n_subgroup = integer(),
    n_rest = integer(),
    beta_subgroup_vs_rest = numeric(),
    conf_low = numeric(),
    conf_high = numeric(),
    p_value = numeric(),
    direction = character(),
    FDR_within_subgroup = numeric(),
    FDR_within_score = numeric(),
    FDR_global = numeric()
  )
}

need_objects <- c(
  "cp10_scores_data",
  "primary_composite_results",
  "primary_results",
  "matched_5y_results",
  "matched_10y_results",
  "threshold_model_results",
  "threshold_robustness_summary",
  "cancer_df",
  "cancer_subgroup_one_vs_rest_results"
)

for (oo in need_objects) {
  if (!exists(oo)) stop("Required object missing after load: ", oo)
}

# ------------------------------------------------------------
# 3. General plotting helpers
# ------------------------------------------------------------

theme_cp10 <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      strip.background = element_rect(fill = "grey95", color = "grey70"),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 35, hjust = 1),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = base_size - 1),
      legend.position = "top"
    )
}

save_plot_both <- function(plot_obj, file_stub, width = 9, height = 6) {
  pdf_file <- file.path(out_dir, paste0(file_stub, ".pdf"))
  png_file <- file.path(out_dir, paste0(file_stub, ".png"))

  ggsave(pdf_file, plot_obj, width = width, height = height, device = cairo_pdf)
  ggsave(png_file, plot_obj, width = width, height = height, dpi = 300)

  tibble(
    figure = file_stub,
    pdf_file = pdf_file,
    png_file = png_file
  )
}

fmt_fdr <- function(x) {
  ifelse(
    is.na(x),
    "FDR=NA",
    ifelse(x < 0.001, paste0("FDR=", formatC(x, format = "e", digits = 2)),
           paste0("FDR=", signif(x, 3)))
  )
}

score_labels <- c(
  CP10_granulocyte_like_enrichment_score =
    "Granulocyte-like\nenrichment",
  CP10_integrated_myeloid_granulocytic_remodeling_score =
    "Integrated myeloid/\ngranulocytic remodeling",
  CP10_myeloid_granulocytic_to_lymphoid_balance_score =
    "Myeloid/granulocytic-\nto-lymphoid balance",
  CP10_monocyte_like_HLA_DR_phenotype_remodeling_score =
    "Monocyte-like\nHLA-DR phenotype",
  CP10_APC_like_myeloid_enrichment_score =
    "APC-like myeloid\nenrichment",
  CP10_CCR3_eosinophil_like_granulocytic_enrichment_score =
    "CCR3 eosinophil-like\ngranulocytic"
)

feature_labels <- c(
  pct_cd13_cd66b_gran_like_within_cd45 =
    "CD13+CD66b+\ngranulocyte-like",
  pct_cd66b_gran_like_within_cd45 =
    "CD66b+\ngranulocyte-like",
  pct_cd14_mono_like_within_cd45 =
    "CD14+\nmonocyte-like",
  pct_lymphoid_like_any_within_cd45 =
    "CD3/CD19/CD56-defined\nlymphoid-like",
  pct_t_like_within_cd45 =
    "CD3+\nT-cell-like",
  ratio_myeloid_granulocytic_to_lymphoid_like =
    "Myeloid/granulocytic-\nto-lymphoid ratio",
  ratio_cd66b_gran_like_to_t_like =
    "CD66b granulocyte-like\n/ T-like ratio",
  pct_cd11c_hladr_apc_like_within_cd45 =
    "CD11c+HLA-DR+\nAPC-like myeloid",
  pct_cd14_cd11c_hladr_pos_within_cd45 =
    "CD14+CD11c+HLA-DR+\nmyeloid APC-like",
  pct_ccr3_cd66b_eosinophil_like_within_cd45 =
    "CCR3+CD66b+\neosinophil-like"
)

figure_index <- list()

# ------------------------------------------------------------
# 4. Figure 1: Composite scores, cancer vs healthy
# ------------------------------------------------------------

robust_composite_scores <- c(
  "CP10_granulocyte_like_enrichment_score",
  "CP10_integrated_myeloid_granulocytic_remodeling_score",
  "CP10_myeloid_granulocytic_to_lymphoid_balance_score",
  "CP10_monocyte_like_HLA_DR_phenotype_remodeling_score",
  "CP10_APC_like_myeloid_enrichment_score"
)

composite_plot_data <- cp10_scores_data %>%
  filter(model_ready_age_sex == TRUE) %>%
  select(subject_id, disease_group, all_of(robust_composite_scores)) %>%
  pivot_longer(
    cols = all_of(robust_composite_scores),
    names_to = "composite_score",
    values_to = "score_value"
  ) %>%
  mutate(
    disease_group = factor(as.character(disease_group), levels = c("Healthy control", "Cancer patient")),
    composite_label = factor(
      score_labels[composite_score],
      levels = score_labels[robust_composite_scores]
    )
  )

composite_stats <- primary_composite_results %>%
  filter(composite_score %in% robust_composite_scores) %>%
  transmute(
    composite_score,
    composite_label = factor(score_labels[composite_score], levels = score_labels[robust_composite_scores]),
    label = paste0(
      "beta=", sprintf("%.2f", beta_cancer_vs_healthy),
      "\n", fmt_fdr(FDR_global)
    ),
    y = max(composite_plot_data$score_value[composite_plot_data$composite_score == composite_score], na.rm = TRUE) + 0.18
  )

fig1 <- ggplot(composite_plot_data, aes(x = disease_group, y = score_value)) +
  geom_boxplot(outlier.shape = NA, width = 0.55) +
  geom_jitter(width = 0.12, alpha = 0.20, size = 0.6) +
  facet_wrap(~ composite_label, scales = "free_y", nrow = 1) +
  geom_text(
    data = composite_stats,
    aes(x = 1.5, y = y, label = label),
    inherit.aes = FALSE,
    size = 3.1
  ) +
  labs(
    title = "CP10 composite scores identify myeloid/granulocytic remodeling in cancer",
    subtitle = "Scores are directionally oriented; higher values indicate the cancer-associated remodeling direction",
    x = NULL,
    y = "Composite score"
  ) +
  theme_cp10(base_size = 10)

figure_index[["Figure1_CP10_composite_scores_cancer_vs_healthy"]] <-
  save_plot_both(fig1, "Figure1_CP10_composite_scores_cancer_vs_healthy", width = 13.5, height = 5.2)

# ------------------------------------------------------------
# 5. Figure 2: Key feature-level axes
# ------------------------------------------------------------

key_features <- c(
  "pct_cd13_cd66b_gran_like_within_cd45",
  "pct_cd66b_gran_like_within_cd45",
  "pct_cd14_mono_like_within_cd45",
  "pct_lymphoid_like_any_within_cd45",
  "ratio_myeloid_granulocytic_to_lymphoid_like",
  "pct_cd11c_hladr_apc_like_within_cd45"
)

feature_plot_data <- cp10_scores_data %>%
  filter(model_ready_age_sex == TRUE) %>%
  select(subject_id, disease_group, all_of(key_features)) %>%
  pivot_longer(
    cols = all_of(key_features),
    names_to = "feature",
    values_to = "feature_value"
  ) %>%
  mutate(
    disease_group = factor(as.character(disease_group), levels = c("Healthy control", "Cancer patient")),
    feature_label = factor(feature_labels[feature], levels = feature_labels[key_features])
  )

feature_stats <- primary_results %>%
  filter(feature %in% key_features) %>%
  transmute(
    feature,
    feature_label = factor(feature_labels[feature], levels = feature_labels[key_features]),
    label = paste0(
      "beta=", sprintf("%.2f", beta_cancer_vs_healthy),
      "\n", fmt_fdr(FDR_global)
    ),
    y = sapply(feature, function(ff) max(feature_plot_data$feature_value[feature_plot_data$feature == ff], na.rm = TRUE) * 1.05)
  )

fig2 <- ggplot(feature_plot_data, aes(x = disease_group, y = feature_value)) +
  geom_boxplot(outlier.shape = NA, width = 0.55) +
  geom_jitter(width = 0.12, alpha = 0.20, size = 0.6) +
  facet_wrap(~ feature_label, scales = "free_y", nrow = 2) +
  geom_text(
    data = feature_stats,
    aes(x = 1.5, y = y, label = label),
    inherit.aes = FALSE,
    size = 3
  ) +
  labs(
    title = "Key CP10 feature-level axes",
    subtitle = "Phenotype-like fractions are shown within viable CD45+ leukocytes unless otherwise noted",
    x = NULL,
    y = "Feature value"
  ) +
  theme_cp10(base_size = 10)

figure_index[["Figure2_CP10_key_feature_axes_cancer_vs_healthy"]] <-
  save_plot_both(fig2, "Figure2_CP10_key_feature_axes_cancer_vs_healthy", width = 11.5, height = 7.0)

# ------------------------------------------------------------
# 6. Figure 3: Feature-level forest plot
# ------------------------------------------------------------

forest_features <- c(
  "pct_cd13_cd66b_gran_like_within_cd45",
  "pct_cd66b_gran_like_within_cd45",
  "pct_cd13_pos_within_cd66b_gran_like",
  "pct_lymphoid_like_any_within_cd45",
  "pct_t_like_within_cd45",
  "pct_cd14_mono_like_within_cd45",
  "ratio_myeloid_granulocytic_to_lymphoid_like",
  "ratio_cd66b_gran_like_to_t_like",
  "ratio_cd14_mono_like_to_t_like",
  "pct_cd14_hladr_low_within_cd45",
  "pct_cd14_cd11c_hladr_pos_within_cd45",
  "pct_cd11c_hladr_apc_like_within_cd45",
  "pct_ccr3_cd66b_eosinophil_like_within_cd45"
)

forest_data <- primary_results %>%
  filter(feature %in% forest_features) %>%
  mutate(
    feature_label = ifelse(feature %in% names(feature_labels), feature_labels[feature], feature),
    feature_label = stringr::str_replace_all(feature_label, "\n", " "),
    feature_label = forcats::fct_reorder(feature_label, beta_cancer_vs_healthy),
    significant = FDR_global < 0.05
  )

fig3 <- ggplot(forest_data, aes(x = beta_cancer_vs_healthy, y = feature_label)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbarh(aes(xmin = conf_low, xmax = conf_high), height = 0.22) +
  geom_point(aes(shape = significant), size = 2.3) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), name = "Global FDR < 0.05") +
  labs(
    title = "CP10 age- and sex-adjusted feature-level effects",
    subtitle = "Positive beta indicates higher values in cancer patients",
    x = "Cancer vs healthy beta",
    y = NULL
  ) +
  theme_cp10(base_size = 10) +
  theme(axis.text.x = element_text(angle = 0))

figure_index[["Figure3_CP10_feature_level_forest_plot"]] <-
  save_plot_both(fig3, "Figure3_CP10_feature_level_forest_plot", width = 9.8, height = 6.6)

# ------------------------------------------------------------
# 7. Figure 4: Threshold sensitivity for CP10 composites
# ------------------------------------------------------------

threshold_composite_data <- threshold_model_results %>%
  filter(variable_type == "composite") %>%
  mutate(
    variable_label = ifelse(variable %in% names(score_labels), score_labels[variable], variable),
    variable_label = stringr::str_replace_all(variable_label, "\n", " "),
    threshold_set = factor(threshold_set, levels = c("main", "permissive", "stringent")),
    variable_label = factor(
      variable_label,
      levels = rev(stringr::str_replace_all(score_labels[unique(threshold_model_results$variable[threshold_model_results$variable_type == "composite"])], "\n", " "))
    )
  )

fig4 <- ggplot(threshold_composite_data, aes(x = beta_cancer_vs_healthy, y = variable_label)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbarh(aes(xmin = conf_low, xmax = conf_high), height = 0.18, position = position_dodge(width = 0.6)) +
  geom_point(aes(shape = threshold_set), size = 2.4, position = position_dodge(width = 0.6)) +
  labs(
    title = "CP10 composite-score threshold sensitivity",
    subtitle = "Main, permissive, and stringent threshold sets",
    x = "Cancer vs healthy beta",
    y = NULL,
    shape = "Threshold set"
  ) +
  theme_cp10(base_size = 10) +
  theme(axis.text.x = element_text(angle = 0))

figure_index[["Figure4_CP10_composite_threshold_sensitivity"]] <-
  save_plot_both(fig4, "Figure4_CP10_composite_threshold_sensitivity", width = 10.5, height = 5.8)

# ------------------------------------------------------------
# 8. Figure 5: Age-matched composite effects
# ------------------------------------------------------------

age_matched_data <- bind_rows(
  primary_composite_results %>%
    mutate(model_for_plot = "Primary age/sex-adjusted") %>%
    select(model_for_plot, composite_score, beta_cancer_vs_healthy, conf_low, conf_high, FDR_global),
  matched_5y_results %>%
    mutate(model_for_plot = "Same-sex age matched, 5-year") %>%
    select(model_for_plot, composite_score, beta_cancer_vs_healthy, conf_low, conf_high, FDR_global),
  matched_10y_results %>%
    mutate(model_for_plot = "Same-sex age matched, 10-year") %>%
    select(model_for_plot, composite_score, beta_cancer_vs_healthy, conf_low, conf_high, FDR_global)
) %>%
  filter(composite_score %in% robust_composite_scores) %>%
  mutate(
    composite_label = stringr::str_replace_all(score_labels[composite_score], "\n", " "),
    composite_label = factor(composite_label, levels = rev(stringr::str_replace_all(score_labels[robust_composite_scores], "\n", " "))),
    model_for_plot = factor(
      model_for_plot,
      levels = c("Primary age/sex-adjusted", "Same-sex age matched, 5-year", "Same-sex age matched, 10-year")
    )
  )

fig5 <- ggplot(age_matched_data, aes(x = beta_cancer_vs_healthy, y = composite_label)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbarh(aes(xmin = conf_low, xmax = conf_high), height = 0.18, position = position_dodge(width = 0.65)) +
  geom_point(aes(shape = model_for_plot), size = 2.4, position = position_dodge(width = 0.65)) +
  labs(
    title = "CP10 composite effects remain after same-sex age matching",
    subtitle = "Positive beta indicates higher composite score in cancer patients",
    x = "Cancer vs healthy beta",
    y = NULL,
    shape = "Model"
  ) +
  theme_cp10(base_size = 10) +
  theme(axis.text.x = element_text(angle = 0))

figure_index[["Figure5_CP10_age_matched_composite_effects"]] <-
  save_plot_both(fig5, "Figure5_CP10_age_matched_composite_effects", width = 10.8, height = 5.8)

# ------------------------------------------------------------
# 9. Figure 6: Breast cancer lower CP10 remodeling magnitude
# ------------------------------------------------------------

breast_scores <- c(
  "CP10_integrated_myeloid_granulocytic_remodeling_score",
  "CP10_monocyte_like_HLA_DR_phenotype_remodeling_score",
  "CP10_granulocyte_like_enrichment_score",
  "CP10_myeloid_granulocytic_to_lymphoid_balance_score"
)

breast_plot_data <- cancer_df %>%
  filter(!is.na(cancer_subgroup_model)) %>%
  mutate(
    subgroup_for_plot = ifelse(as.character(cancer_subgroup_model) == "Meme", "Breast cancer", "Other cancers"),
    subgroup_for_plot = factor(subgroup_for_plot, levels = c("Other cancers", "Breast cancer"))
  ) %>%
  select(subject_id, subgroup_for_plot, all_of(breast_scores)) %>%
  pivot_longer(
    cols = all_of(breast_scores),
    names_to = "score",
    values_to = "score_value"
  ) %>%
  mutate(
    score_label = factor(score_labels[score], levels = score_labels[breast_scores])
  )

breast_stats <- cancer_subgroup_one_vs_rest_results %>%
  filter(cancer_subgroup == "Meme", score %in% breast_scores) %>%
  transmute(
    score,
    score_label = factor(score_labels[score], levels = score_labels[breast_scores]),
    label = paste0(
      "beta=", sprintf("%.2f", beta_subgroup_vs_rest),
      "\n", fmt_fdr(FDR_global)
    ),
    y = sapply(score, function(sc) max(breast_plot_data$score_value[breast_plot_data$score == sc], na.rm = TRUE) + 0.18)
  )

# If subgroup statistics are unavailable for any reason, still draw Figure 6
# without annotation labels.
if (nrow(breast_stats) == 0) {
  breast_stats <- tibble(
    score = character(),
    score_label = factor(character(), levels = score_labels[breast_scores]),
    label = character(),
    y = numeric()
  )
}

fig6 <- ggplot(breast_plot_data, aes(x = subgroup_for_plot, y = score_value)) +
  geom_boxplot(outlier.shape = NA, width = 0.55) +
  geom_jitter(width = 0.12, alpha = 0.22, size = 0.65) +
  facet_wrap(~ score_label, scales = "free_y", nrow = 1) +
  geom_text(
    data = breast_stats,
    aes(x = 1.5, y = y, label = label),
    inherit.aes = FALSE,
    size = 3
  ) +
  labs(
    title = "Cancer-subgroup annotation: lower CP10 remodeling magnitude in breast cancer",
    subtitle = "Cancer-only one-vs-rest comparison; shown as secondary clinical annotation",
    x = NULL,
    y = "Composite score"
  ) +
  theme_cp10(base_size = 10)

figure_index[["Figure6_CP10_breast_cancer_lower_myeloid_remodeling"]] <-
  save_plot_both(fig6, "Figure6_CP10_breast_cancer_lower_myeloid_remodeling", width = 12.5, height = 5.2)

# ------------------------------------------------------------
# 10. Figure manifest and compact tables
# ------------------------------------------------------------

figure_manifest <- bind_rows(figure_index)

readr::write_csv(
  figure_manifest,
  file.path(out_dir, "SDY2583_CP10_STEP7A_main_figure_manifest.csv")
)

# Tables used in figure annotations
readr::write_csv(
  primary_composite_results,
  file.path(out_dir, "SDY2583_CP10_STEP7A_primary_composite_results_for_figures.csv")
)

readr::write_csv(
  primary_results,
  file.path(out_dir, "SDY2583_CP10_STEP7A_primary_feature_results_for_figures.csv")
)

readr::write_csv(
  cancer_subgroup_one_vs_rest_results,
  file.path(out_dir, "SDY2583_CP10_STEP7A_cancer_subgroup_one_vs_rest_for_figures.csv")
)

save(
  figure_manifest,
  composite_plot_data,
  feature_plot_data,
  forest_data,
  threshold_composite_data,
  age_matched_data,
  breast_plot_data,
  file = file.path(cp10_rdata_dir, "SDY2583_CP10_STEP7A_main_figures.RData")
)

cat("\n============================================================\n")
cat("SDY2583 CP10 STEP 7A COMPLETE: MAIN FIGURES\n")
cat("============================================================\n")

cat("\nFigures saved in:\n")
print(out_dir)

cat("\nFigure manifest:\n")
print(as.data.frame(figure_manifest), row.names = FALSE)

cat("\nGenerated figures:\n")
cat("1. Figure1_CP10_composite_scores_cancer_vs_healthy\n")
cat("2. Figure2_CP10_key_feature_axes_cancer_vs_healthy\n")
cat("3. Figure3_CP10_feature_level_forest_plot\n")
cat("4. Figure4_CP10_composite_threshold_sensitivity\n")
cat("5. Figure5_CP10_age_matched_composite_effects\n")
cat("6. Figure6_CP10_breast_cancer_lower_myeloid_remodeling\n")

cat("============================================================\n")
