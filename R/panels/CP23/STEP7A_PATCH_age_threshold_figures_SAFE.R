source(file.path(Sys.getenv("SDY2583_REPO_ROOT", unset = "."), "R", "shared", "bootstrap.R"))

# ============================================================
# SDY2583 CP23
# STEP 7A PATCH SAFE: Add missing age-matching and threshold figures
#
# Why:
#   Step7A completed, but age-matching and threshold figures were skipped
#   because automatic file discovery selected summary files without
#   composite_score/variable columns.
#
# This patch reads the exact Step4B and Step5 result files and creates:
#   1) CP23_age_matching_sensitivity_forest
#   2) CP23_threshold_sensitivity_composite_betas
#
# It also appends/updates:
#   outputs/CP23/08_quantitative_figures/
#      CP23_STEP7A_figure_manifest.csv
#      figure_source_tables/CP23_age_matching_sensitivity_plot_table.csv
#      figure_source_tables/CP23_threshold_sensitivity_composite_plot_table.csv
# ============================================================

cran_pkgs <- c("dplyr", "readr", "tibble", "tidyr", "stringr", "ggplot2", "forcats")
for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(forcats)
})

analysis_dir <- sd_analysis_dir("CP23")
fig_dir <- file.path(analysis_dir, "08_quantitative_figures")
main_dir <- file.path(fig_dir, "main_candidate_figures")
table_dir <- file.path(fig_dir, "figure_source_tables")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(main_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

score_labels <- tibble(
  composite_score = c(
    "CP23_integrated_monocyte_macrophage_like_myeloid_remodeling_score",
    "CP23_CD33_HLA_DR_myeloid_repatterning_score",
    "CP23_CD9_CD84_activation_attenuation_score",
    "CP23_CD45_dump_low_myeloid_enrichment_score",
    "CP23_FcERI_myeloid_attenuation_score"
  ),
  score_label = c(
    "Integrated CP23 monocyte/macrophage-like myeloid remodeling",
    "CD33/HLA-DR myeloid repatterning",
    "CD9/CD84 activation attenuation",
    "CD45+ dump-low myeloid enrichment",
    "FcERI-associated myeloid attenuation"
  )
)

theme_pub <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      strip.background = element_rect(fill = "white", colour = "black"),
      plot.title = element_text(face = "bold", hjust = 0),
      plot.subtitle = element_text(hjust = 0),
      axis.title = element_text(face = "bold"),
      legend.position = "right"
    )
}

label_p <- function(p) {
  ifelse(is.na(p), "NA",
         ifelse(p < 0.001, formatC(p, format = "e", digits = 2),
                sprintf("%.3f", p)))
}

save_plot <- function(p, file_stub, width = 8.5, height = 6.2) {
  png_file <- file.path(main_dir, paste0(file_stub, ".png"))
  pdf_file <- file.path(main_dir, paste0(file_stub, ".pdf"))
  ggsave(png_file, p, width = width, height = height, dpi = 320)
  ggsave(pdf_file, p, width = width, height = height, device = cairo_pdf)
  tibble(figure_stub = file_stub, png_file = png_file, pdf_file = pdf_file)
}

new_manifest_rows <- list()

# ------------------------------------------------------------
# 1. Age-matching sensitivity figure
# ------------------------------------------------------------

matched_file <- file.path(
  analysis_dir,
  "06_age_sensitivity_caliper_interaction",
  "SDY2583_CP23_matched_direction_sensitivity_summary_STEP4B.csv"
)

if (!file.exists(matched_file)) {
  stop("Matched direction sensitivity file not found: ", matched_file)
}

matched_summary <- read_csv(matched_file, show_col_types = FALSE)

if (!("composite_score" %in% names(matched_summary))) {
  stop("Matched file does not contain composite_score: ", matched_file)
}

matched_summary <- matched_summary %>%
  left_join(score_labels, by = "composite_score") %>%
  mutate(score_label = ifelse(is.na(score_label), composite_score, score_label))

required_matched_cols <- c(
  "matched_5y_beta",
  "matched_10y_beta",
  "matched_5y_FDR_global",
  "matched_10y_FDR_global"
)

missing_matched_cols <- setdiff(required_matched_cols, names(matched_summary))
if (length(missing_matched_cols) > 0) {
  stop("Matched file is missing columns: ", paste(missing_matched_cols, collapse = ", "))
}

matched_long <- matched_summary %>%
  select(
    composite_score,
    score_label,
    primary_beta,
    primary_FDR_global,
    matched_5y_beta,
    matched_5y_FDR_global,
    matched_10y_beta,
    matched_10y_FDR_global,
    direction_preserved_both_matched,
    FDR_global_preserved_both_matched
  ) %>%
  pivot_longer(
    cols = c(matched_5y_beta, matched_10y_beta),
    names_to = "matched_model",
    values_to = "beta"
  ) %>%
  mutate(
    FDR_global = ifelse(
      matched_model == "matched_5y_beta",
      matched_5y_FDR_global,
      matched_10y_FDR_global
    ),
    caliper = ifelse(
      matched_model == "matched_5y_beta",
      "Same-sex age-matched 5-year",
      "Same-sex age-matched 10-year"
    ),
    caliper = factor(caliper, levels = c("Same-sex age-matched 5-year", "Same-sex age-matched 10-year")),
    score_label = fct_reorder(score_label, beta)
  )

read_csv_out <- file.path(table_dir, "CP23_age_matching_sensitivity_plot_table.csv")
write_csv(matched_long, read_csv_out)

p_matched <- ggplot(matched_long, aes(x = beta, y = score_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35) +
  geom_point(size = 2.4) +
  geom_text(aes(label = paste0("FDR=", label_p(FDR_global))), nudge_y = 0.22, size = 3) +
  facet_wrap(~ caliper, ncol = 1) +
  labs(
    title = "CP23 age-matching sensitivity",
    subtitle = "Same-sex nearest-neighbor age matching; positive beta indicates higher composite score in cancer.",
    x = "Cancer vs healthy matched beta",
    y = NULL
  ) +
  theme_pub(11)

new_manifest_rows[["age_matching"]] <- save_plot(
  p_matched,
  "CP23_age_matching_sensitivity_forest",
  width = 8.8,
  height = 6.2
)

# ------------------------------------------------------------
# 2. Threshold-sensitivity composite beta figure
# ------------------------------------------------------------

threshold_file <- file.path(
  analysis_dir,
  "07_threshold_sensitivity",
  "SDY2583_CP23_threshold_wide_robustness_STEP5.csv"
)

if (!file.exists(threshold_file)) {
  stop("Threshold wide robustness file not found: ", threshold_file)
}

threshold_wide <- read_csv(threshold_file, show_col_types = FALSE)

if (!("variable" %in% names(threshold_wide))) {
  stop("Threshold file does not contain variable column: ", threshold_file)
}

threshold_comp <- threshold_wide %>%
  filter(variable_type == "composite") %>%
  rename(composite_score = variable) %>%
  filter(composite_score %in% score_labels$composite_score) %>%
  left_join(score_labels, by = "composite_score") %>%
  mutate(score_label = ifelse(is.na(score_label), composite_score, score_label))

required_threshold_cols <- c(
  "beta_cancer_vs_healthy__main",
  "beta_cancer_vs_healthy__permissive",
  "beta_cancer_vs_healthy__stringent",
  "FDR_global__main",
  "FDR_global__permissive",
  "FDR_global__stringent"
)

missing_threshold_cols <- setdiff(required_threshold_cols, names(threshold_comp))
if (length(missing_threshold_cols) > 0) {
  stop("Threshold file is missing columns: ", paste(missing_threshold_cols, collapse = ", "))
}

threshold_long <- threshold_comp %>%
  select(
    composite_score,
    score_label,
    robustness_class,
    direction_preserved_all_sets,
    global_FDR_lt_0p05_all_sets,
    all_of(required_threshold_cols)
  ) %>%
  pivot_longer(
    cols = c(
      beta_cancer_vs_healthy__main,
      beta_cancer_vs_healthy__permissive,
      beta_cancer_vs_healthy__stringent
    ),
    names_to = "threshold_set",
    values_to = "beta"
  ) %>%
  mutate(
    FDR_global = case_when(
      threshold_set == "beta_cancer_vs_healthy__main" ~ FDR_global__main,
      threshold_set == "beta_cancer_vs_healthy__permissive" ~ FDR_global__permissive,
      threshold_set == "beta_cancer_vs_healthy__stringent" ~ FDR_global__stringent,
      TRUE ~ NA_real_
    ),
    threshold_set = case_when(
      threshold_set == "beta_cancer_vs_healthy__permissive" ~ "Permissive",
      threshold_set == "beta_cancer_vs_healthy__main" ~ "Main",
      threshold_set == "beta_cancer_vs_healthy__stringent" ~ "Stringent",
      TRUE ~ threshold_set
    ),
    threshold_set = factor(threshold_set, levels = c("Permissive", "Main", "Stringent")),
    score_label = fct_reorder(score_label, beta)
  )

write_csv(threshold_long, file.path(table_dir, "CP23_threshold_sensitivity_composite_plot_table.csv"))

p_threshold <- ggplot(threshold_long, aes(x = threshold_set, y = beta, group = score_label)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35) +
  geom_line(linewidth = 0.45) +
  geom_point(size = 2) +
  geom_text(aes(label = paste0("FDR=", label_p(FDR_global))), nudge_y = 0.06, size = 2.7) +
  facet_wrap(~ score_label, ncol = 2, scales = "free_y") +
  labs(
    title = "CP23 threshold sensitivity of composite effects",
    subtitle = "Composite beta estimates across permissive, main, and stringent marker thresholds.",
    x = "Threshold set",
    y = "Cancer vs healthy adjusted beta"
  ) +
  theme_pub(10)

new_manifest_rows[["threshold"]] <- save_plot(
  p_threshold,
  "CP23_threshold_sensitivity_composite_betas",
  width = 9.4,
  height = 7.2
)

# ------------------------------------------------------------
# 3. Update manifest
# ------------------------------------------------------------

manifest_file <- file.path(fig_dir, "CP23_STEP7A_figure_manifest.csv")
new_manifest <- bind_rows(new_manifest_rows)

if (file.exists(manifest_file)) {
  old_manifest <- read_csv(manifest_file, show_col_types = FALSE)
  updated_manifest <- old_manifest %>%
    filter(!(figure_stub %in% new_manifest$figure_stub)) %>%
    bind_rows(new_manifest) %>%
    arrange(figure_stub)
} else {
  updated_manifest <- new_manifest
}

write_csv(updated_manifest, manifest_file)

cat("\n============================================================\n")
cat("SDY2583 CP23 STEP 7A PATCH COMPLETE: AGE + THRESHOLD FIGURES\n")
cat("============================================================\n")

cat("\nNew/updated figures:\n")
print(as.data.frame(new_manifest), row.names = FALSE)

cat("\nUpdated figure manifest:\n")
print(as.data.frame(updated_manifest), row.names = FALSE)

cat("\nSource tables written:\n")
print(read_csv_out)
print(file.path(table_dir, "CP23_threshold_sensitivity_composite_plot_table.csv"))

cat("\nMain figure folder:\n")
print(main_dir)

cat("============================================================\n")
