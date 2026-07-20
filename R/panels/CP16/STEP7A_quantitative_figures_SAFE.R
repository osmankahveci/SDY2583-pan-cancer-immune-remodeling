# ============================================================
# SDY2583 CP16
# STEP 7A SAFE: Quantitative figures
#
# Purpose:
#   Generate manuscript/supplement-ready quantitative figures for
#   CP16 APC/DC-like myeloid remodeling.
#
# Inputs:
#   outputs/CP16/11_RData/
#     SDY2583_CP16_STEP4_composite_scores.RData
#     SDY2583_CP16_STEP4B_age_sensitivity_caliper_interaction.RData
#     SDY2583_CP16_STEP5_threshold_sensitivity_TARGETED.RData
#
#   data/derived/clinical_integration/RData/
#     SDY2583_CP16_clinical_annotation_STEP6.RData
#
# Outputs:
#   outputs/CP16/08_quantitative_figures
#
# Notes:
#   - This script generates figures from real analysis outputs only.
#   - No synthetic flow events are generated.
#   - FlowJo-style representative gating is NOT included here.
# ============================================================

rm(list = ls())
source(file.path(Sys.getenv("SDY2583_REPO_ROOT", unset = "."), "R", "shared", "bootstrap.R"))

# ------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------

cran_pkgs <- c(
  "dplyr", "readr", "stringr", "tibble", "tidyr",
  "ggplot2", "forcats", "scales"
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
fig_dir <- file.path(analysis_dir, "08_quantitative_figures")

main_fig_dir <- file.path(fig_dir, "main_candidate_figures")
supp_fig_dir <- file.path(fig_dir, "supplementary_candidate_figures")
table_dir <- file.path(fig_dir, "figure_source_tables")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(main_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

matrix_rdata_dir <- file.path(sd_integrated_dir(), "RData")

step4_rdata <- file.path(rdata_dir, "SDY2583_CP16_STEP4_composite_scores.RData")
step4b_rdata <- file.path(rdata_dir, "SDY2583_CP16_STEP4B_age_sensitivity_caliper_interaction.RData")
step5_rdata <- file.path(rdata_dir, "SDY2583_CP16_STEP5_threshold_sensitivity_TARGETED.RData")
step6_rdata <- file.path(matrix_rdata_dir, "SDY2583_CP16_clinical_annotation_STEP6.RData")

required <- c(step4_rdata, step4b_rdata, step5_rdata, step6_rdata)
missing_required <- required[!file.exists(required)]

if (length(missing_required) > 0) {
  stop("Eksik RData dosyaları:\n", paste(missing_required, collapse = "\n"))
}

# ------------------------------------------------------------
# 3. Load data
# ------------------------------------------------------------

load(step4_rdata)
load(step4b_rdata)
load(step5_rdata)
load(step6_rdata)

# Reset paths after RData loads.
analysis_dir <- sd_analysis_dir("CP16")
fig_dir <- file.path(analysis_dir, "08_quantitative_figures")
main_fig_dir <- file.path(fig_dir, "main_candidate_figures")
supp_fig_dir <- file.path(fig_dir, "supplementary_candidate_figures")
table_dir <- file.path(fig_dir, "figure_source_tables")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(main_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

if (!exists("cp16_scores_data")) stop("cp16_scores_data bulunamadı.")
if (!exists("composite_score_cols")) stop("composite_score_cols bulunamadı.")
if (!exists("primary_composite_results")) stop("primary_composite_results bulunamadı.")
if (!exists("matched_5y_results")) stop("matched_5y_results bulunamadı.")
if (!exists("matched_10y_results")) stop("matched_10y_results bulunamadı.")
if (!exists("interaction_results")) stop("interaction_results bulunamadı.")
if (!exists("threshold_model_results")) stop("threshold_model_results bulunamadı.")
if (!exists("threshold_robustness_summary")) stop("threshold_robustness_summary bulunamadı.")

# ------------------------------------------------------------
# 4. Helper functions
# ------------------------------------------------------------

pretty_score <- function(x) {
  dplyr::case_when(
    x == "CP16_integrated_APC_DC_myeloid_remodeling_score" ~ "Integrated APC/DC-like\nmyeloid remodeling",
    x == "CP16_CD1c_cDC2_like_depletion_score" ~ "CD1c/cDC2-like\ndepletion",
    x == "CP16_CD123_pDC_like_depletion_score" ~ "CD123/pDC-like\ndepletion",
    x == "CP16_CD141_CLEC9A_cDC1_like_depletion_score" ~ "CD141/CLEC9A\ncDC1-like depletion",
    x == "CP16_monocyte_HLA_DR_low_remodeling_score" ~ "Monocyte HLA-DR-low\nremodeling",
    x == "CP16_APC_core_attenuation_score" ~ "APC core\nattenuation",
    TRUE ~ stringr::str_replace_all(x, "_", " ")
  )
}

pretty_score_one_line <- function(x) {
  dplyr::case_when(
    x == "CP16_integrated_APC_DC_myeloid_remodeling_score" ~ "Integrated APC/DC-like myeloid remodeling",
    x == "CP16_CD1c_cDC2_like_depletion_score" ~ "CD1c/cDC2-like depletion",
    x == "CP16_CD123_pDC_like_depletion_score" ~ "CD123/pDC-like depletion",
    x == "CP16_CD141_CLEC9A_cDC1_like_depletion_score" ~ "CD141/CLEC9A cDC1-like depletion",
    x == "CP16_monocyte_HLA_DR_low_remodeling_score" ~ "Monocyte HLA-DR-low remodeling",
    x == "CP16_APC_core_attenuation_score" ~ "APC core attenuation",
    TRUE ~ stringr::str_replace_all(x, "_", " ")
  )
}

format_fdr <- function(x) {
  ifelse(
    is.na(x),
    "FDR = NA",
    ifelse(x < 0.001, paste0("FDR = ", formatC(x, format = "e", digits = 2)),
           paste0("FDR = ", signif(x, 3)))
  )
}

save_plot_both <- function(plot_obj, file_stub, out_dir, width = 7, height = 5) {
  png_file <- file.path(out_dir, paste0(file_stub, ".png"))
  pdf_file <- file.path(out_dir, paste0(file_stub, ".pdf"))

  ggplot2::ggsave(png_file, plot_obj, width = width, height = height, dpi = 320, bg = "white")
  ggplot2::ggsave(pdf_file, plot_obj, width = width, height = height, device = cairo_pdf, bg = "white")

  tibble(
    figure_stub = file_stub,
    png_file = png_file,
    pdf_file = pdf_file
  )
}

base_theme <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey95", colour = "grey70"),
      strip.text = element_text(face = "bold"),
      legend.position = "right"
    )
}

figure_manifest <- list()

# ------------------------------------------------------------
# 5. Source table: composite score long data
# ------------------------------------------------------------

composite_long <- cp16_scores_data %>%
  filter(model_ready_age_sex == TRUE) %>%
  select(subject_id, disease_group, age_for_model, sex, all_of(composite_score_cols)) %>%
  pivot_longer(
    cols = all_of(composite_score_cols),
    names_to = "composite_score",
    values_to = "score_value"
  ) %>%
  left_join(
    primary_composite_results %>%
      select(composite_score, beta_cancer_vs_healthy, conf_low, conf_high, p_value, FDR_global),
    by = "composite_score"
  ) %>%
  mutate(
    score_label = pretty_score(composite_score),
    score_label_one_line = pretty_score_one_line(composite_score),
    disease_group = factor(as.character(disease_group), levels = c("Healthy control", "Cancer patient"))
  )

write_csv(composite_long, file.path(table_dir, "CP16_Figure_source_composite_scores_long.csv"))

# ------------------------------------------------------------
# 6. Figure 1: Composite score distributions
#    One separate figure per composite, plus one multipage PDF
# ------------------------------------------------------------

individual_boxplot_manifest <- list()

for (sc in composite_score_cols) {

  df_sc <- composite_long %>%
    filter(composite_score == sc)

  fdr_lab <- format_fdr(unique(df_sc$FDR_global)[1])
  beta_lab <- paste0("β = ", signif(unique(df_sc$beta_cancer_vs_healthy)[1], 3))

  p <- ggplot(df_sc, aes(x = disease_group, y = score_value)) +
    geom_boxplot(width = 0.55, outlier.shape = NA) +
    geom_jitter(width = 0.16, height = 0, alpha = 0.35, size = 1.1) +
    labs(
      title = pretty_score_one_line(sc),
      subtitle = paste(beta_lab, fdr_lab, sep = "; "),
      x = NULL,
      y = "Composite score"
    ) +
    base_theme(12)

  stub <- paste0("CP16_boxplot_", sc)
  individual_boxplot_manifest[[sc]] <- save_plot_both(p, stub, main_fig_dir, width = 6.2, height = 4.8)
}

figure_manifest[["individual_boxplots"]] <- bind_rows(individual_boxplot_manifest)

# Multipage PDF: one composite per page.
multipage_boxplot_pdf <- file.path(main_fig_dir, "CP16_multipage_composite_boxplots_one_score_per_page.pdf")
grDevices::cairo_pdf(multipage_boxplot_pdf, width = 6.2, height = 4.8)

for (sc in composite_score_cols) {
  df_sc <- composite_long %>% filter(composite_score == sc)
  fdr_lab <- format_fdr(unique(df_sc$FDR_global)[1])
  beta_lab <- paste0("β = ", signif(unique(df_sc$beta_cancer_vs_healthy)[1], 3))

  p <- ggplot(df_sc, aes(x = disease_group, y = score_value)) +
    geom_boxplot(width = 0.55, outlier.shape = NA) +
    geom_jitter(width = 0.16, height = 0, alpha = 0.35, size = 1.1) +
    labs(
      title = pretty_score_one_line(sc),
      subtitle = paste(beta_lab, fdr_lab, sep = "; "),
      x = NULL,
      y = "Composite score"
    ) +
    base_theme(12)

  print(p)
}

dev.off()

figure_manifest[["multipage_boxplot_pdf"]] <- tibble(
  figure_stub = "CP16_multipage_composite_boxplots_one_score_per_page",
  png_file = NA_character_,
  pdf_file = multipage_boxplot_pdf
)

# ------------------------------------------------------------
# 7. Figure 2: Primary composite forest plot
# ------------------------------------------------------------

primary_forest_df <- primary_composite_results %>%
  mutate(
    score_label = pretty_score_one_line(composite_score),
    score_label = forcats::fct_reorder(score_label, beta_cancer_vs_healthy)
  )

write_csv(primary_forest_df, file.path(table_dir, "CP16_Figure_source_primary_composite_forest.csv"))

p_primary_forest <- ggplot(
  primary_forest_df,
  aes(x = beta_cancer_vs_healthy, y = score_label)
) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(xmin = conf_low, xmax = conf_high), size = 0.5) +
  labs(
    title = "CP16 composite scores: age/sex-adjusted cancer effect",
    x = "β for cancer vs healthy control",
    y = NULL
  ) +
  base_theme(12)

figure_manifest[["primary_forest"]] <- save_plot_both(
  p_primary_forest,
  "CP16_primary_composite_forest",
  main_fig_dir,
  width = 7.6,
  height = 4.8
)

# ------------------------------------------------------------
# 8. Figure 3: Age-matching sensitivity forest plot
# ------------------------------------------------------------

matched_forest_df <- bind_rows(
  primary_composite_results %>%
    transmute(
      composite_score,
      model_label = "Primary age/sex-adjusted",
      beta = beta_cancer_vs_healthy,
      conf_low,
      conf_high,
      FDR_global
    ),
  matched_5y_results %>%
    transmute(
      composite_score,
      model_label = "Same-sex age-matched 5-year",
      beta = beta_cancer_vs_healthy,
      conf_low,
      conf_high,
      FDR_global
    ),
  matched_10y_results %>%
    transmute(
      composite_score,
      model_label = "Same-sex age-matched 10-year",
      beta = beta_cancer_vs_healthy,
      conf_low,
      conf_high,
      FDR_global
    )
) %>%
  mutate(
    model_label = factor(
      model_label,
      levels = c(
        "Primary age/sex-adjusted",
        "Same-sex age-matched 5-year",
        "Same-sex age-matched 10-year"
      )
    ),
    score_label = pretty_score_one_line(composite_score),
    score_label = forcats::fct_reorder(score_label, beta, .fun = median)
  )

write_csv(matched_forest_df, file.path(table_dir, "CP16_Figure_source_age_matching_forest.csv"))

p_matched_forest <- ggplot(
  matched_forest_df,
  aes(x = beta, y = score_label, shape = model_label)
) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(
    aes(xmin = conf_low, xmax = conf_high),
    position = position_dodge(width = 0.65),
    size = 0.4
  ) +
  labs(
    title = "CP16 composite scores: age-matching sensitivity",
    x = "β for cancer vs healthy control",
    y = NULL,
    shape = NULL
  ) +
  base_theme(11)

figure_manifest[["age_matching_forest"]] <- save_plot_both(
  p_matched_forest,
  "CP16_age_matching_sensitivity_forest",
  main_fig_dir,
  width = 8.5,
  height = 5.2
)

# ------------------------------------------------------------
# 9. Figure 4: Threshold sensitivity for composite scores
# ------------------------------------------------------------

threshold_composite_df <- threshold_model_results %>%
  filter(variable_type == "composite") %>%
  mutate(
    composite_score = variable,
    score_label = pretty_score_one_line(composite_score),
    threshold_set = factor(threshold_set, levels = c("permissive", "main", "stringent")),
    score_label = forcats::fct_reorder(score_label, beta_cancer_vs_healthy, .fun = median)
  )

write_csv(threshold_composite_df, file.path(table_dir, "CP16_Figure_source_threshold_composite_betas.csv"))

p_threshold <- ggplot(
  threshold_composite_df,
  aes(x = threshold_set, y = beta_cancer_vs_healthy, group = score_label)
) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(alpha = 0.7) +
  geom_point(size = 2) +
  facet_wrap(~ score_label, ncol = 2, scales = "free_y") +
  labs(
    title = "CP16 composite scores: targeted threshold sensitivity",
    x = "Threshold set",
    y = "β for cancer vs healthy control"
  ) +
  base_theme(10) +
  theme(legend.position = "none")

figure_manifest[["threshold_composite"]] <- save_plot_both(
  p_threshold,
  "CP16_threshold_sensitivity_composite_betas",
  main_fig_dir,
  width = 8.2,
  height = 7.2
)

# ------------------------------------------------------------
# 10. Figure 5: Disease-by-age interaction forest
# ------------------------------------------------------------

interaction_plot_df <- interaction_results %>%
  mutate(
    score_label = pretty_score_one_line(composite_score),
    score_label = forcats::fct_reorder(score_label, beta_interaction_cancer_by_age_z)
  )

write_csv(interaction_plot_df, file.path(table_dir, "CP16_Figure_source_disease_by_age_interaction.csv"))

p_interaction <- ggplot(
  interaction_plot_df,
  aes(x = beta_interaction_cancer_by_age_z, y = score_label)
) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(xmin = conf_low, xmax = conf_high), size = 0.5) +
  labs(
    title = "CP16 composite scores: disease-by-age interaction",
    subtitle = "Positive values indicate stronger cancer-associated remodeling at older age",
    x = "β for cancer × age z-score interaction",
    y = NULL
  ) +
  base_theme(11)

figure_manifest[["interaction_forest"]] <- save_plot_both(
  p_interaction,
  "CP16_disease_by_age_interaction_forest",
  supp_fig_dir,
  width = 8.2,
  height = 4.8
)

# ------------------------------------------------------------
# 11. Figure 6: Clinical annotation, tumor subgroup one-vs-rest
# ------------------------------------------------------------

if (exists("cancer_subgroup_one_vs_rest_results") &&
    nrow(cancer_subgroup_one_vs_rest_results) > 0) {

  subgroup_plot_df <- cancer_subgroup_one_vs_rest_results %>%
    mutate(
      score_label = pretty_score_one_line(score),
      cancer_subgroup = factor(cancer_subgroup),
      nominal = p_value < 0.05
    ) %>%
    arrange(p_value)

  write_csv(subgroup_plot_df, file.path(table_dir, "CP16_Figure_source_clinical_subgroup_one_vs_rest.csv"))

  p_subgroup <- ggplot(
    subgroup_plot_df,
    aes(x = beta, y = score_label, shape = nominal)
  ) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_pointrange(
      aes(xmin = conf_low, xmax = conf_high),
      size = 0.35
    ) +
    facet_wrap(~ cancer_subgroup, ncol = 2) +
    labs(
      title = "CP16 clinical annotation: cancer subgroup one-vs-rest",
      subtitle = "Exploratory; primary interpretation should use global FDR",
      x = "β for subgroup vs rest",
      y = NULL,
      shape = "Nominal p < 0.05"
    ) +
    base_theme(9)

  figure_manifest[["clinical_subgroup"]] <- save_plot_both(
    p_subgroup,
    "CP16_clinical_subgroup_one_vs_rest_exploratory",
    supp_fig_dir,
    width = 9.5,
    height = 9.5
  )
}

# ------------------------------------------------------------
# 12. Figure 7: Clinical annotation, active-treatment modalities
# ------------------------------------------------------------

if (exists("active_modality_results") &&
    nrow(active_modality_results) > 0) {

  modality_plot_df <- active_modality_results %>%
    filter(analysis_tier != "descriptive_only") %>%
    mutate(
      score_label = pretty_score_one_line(score),
      exposure = factor(exposure),
      nominal = p_value < 0.05
    )

  write_csv(modality_plot_df, file.path(table_dir, "CP16_Figure_source_active_treatment_modality.csv"))

  if (nrow(modality_plot_df) > 0) {

    p_modality <- ggplot(
      modality_plot_df,
      aes(x = beta, y = score_label, shape = nominal)
    ) +
      geom_vline(xintercept = 0, linetype = "dashed") +
      geom_pointrange(
        aes(xmin = conf_low, xmax = conf_high),
        size = 0.35
      ) +
      facet_wrap(~ exposure, ncol = 2) +
      labs(
        title = "CP16 clinical annotation: active-treatment modality",
        subtitle = "Exploratory; no causal treatment-response interpretation",
        x = "β for exposed vs unexposed within active-treatment patients",
        y = NULL,
        shape = "Nominal p < 0.05"
      ) +
      base_theme(9)

    figure_manifest[["clinical_modality"]] <- save_plot_both(
      p_modality,
      "CP16_active_treatment_modality_exploratory",
      supp_fig_dir,
      width = 9.5,
      height = 7.8
    )
  }
}

# ------------------------------------------------------------
# 13. Compact results table for manuscript writing
# ------------------------------------------------------------

cp16_quantitative_figure_key_results <- primary_composite_results %>%
  transmute(
    composite_score,
    score_label = pretty_score_one_line(composite_score),
    n_model,
    n_healthy,
    n_cancer,
    beta_primary = beta_cancer_vs_healthy,
    conf_low_primary = conf_low,
    conf_high_primary = conf_high,
    p_value_primary = p_value,
    FDR_global_primary = FDR_global
  ) %>%
  left_join(
    matched_direction_summary %>%
      select(
        composite_score,
        matched_5y_beta,
        matched_5y_FDR_global,
        matched_10y_beta,
        matched_10y_FDR_global,
        FDR_global_preserved_both_matched
      ),
    by = "composite_score"
  ) %>%
  left_join(
    threshold_robustness_summary %>%
      filter(variable_type == "composite") %>%
      transmute(
        composite_score = variable,
        threshold_robustness_class = robustness_class,
        threshold_direction_preserved_all_sets = direction_preserved_all_sets,
        threshold_global_FDR_lt_0p05_all_sets = global_FDR_lt_0p05_all_sets
      ),
    by = "composite_score"
  ) %>%
  arrange(FDR_global_primary)

write_csv(
  cp16_quantitative_figure_key_results,
  file.path(table_dir, "CP16_quantitative_figure_key_results.csv")
)

# ------------------------------------------------------------
# 14. Save figure manifest
# ------------------------------------------------------------

figure_manifest_df <- bind_rows(figure_manifest, .id = "figure_group") %>%
  arrange(figure_group, figure_stub)

write_csv(
  figure_manifest_df,
  file.path(fig_dir, "SDY2583_CP16_STEP7A_quantitative_figure_manifest.csv")
)

save(
  composite_long,
  primary_forest_df,
  matched_forest_df,
  threshold_composite_df,
  interaction_plot_df,
  cp16_quantitative_figure_key_results,
  figure_manifest_df,
  file = file.path(rdata_dir, "SDY2583_CP16_STEP7A_quantitative_figures.RData")
)

# ------------------------------------------------------------
# 15. Console output
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SDY2583 CP16 STEP 7A COMPLETE: QUANTITATIVE FIGURES\n")
cat("============================================================\n")

cat("\nFigure output folder:\n")
print(fig_dir)

cat("\nMain candidate figures:\n")
print(main_fig_dir)

cat("\nSupplementary candidate figures:\n")
print(supp_fig_dir)

cat("\nFigure source tables:\n")
print(table_dir)

cat("\nFigure manifest:\n")
print(as.data.frame(figure_manifest_df), row.names = FALSE)

cat("\nKey results table:\n")
print(as.data.frame(cp16_quantitative_figure_key_results), row.names = FALSE)

cat("============================================================\n")
