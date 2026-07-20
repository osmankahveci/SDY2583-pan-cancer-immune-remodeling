# ============================================================
# SDY2583 CP22
# STEP 7A SAFE: Quantitative manuscript figures
#
# Creates:
# 1) Composite-score forest plot
# 2) Feature-level robust B-cell remodeling forest plot
# 3) Threshold-sensitivity robustness plot
# 4) CP22 clinical annotation overview plot, if Step 6 exists
#
# Output:
# outputs/CP22/09_main_figures_bracketed
# ============================================================

rm(list = ls())
source(file.path(Sys.getenv("SDY2583_REPO_ROOT", unset = "."), "R", "shared", "bootstrap.R"))

pkgs <- c("dplyr", "readr", "stringr", "ggplot2", "tibble", "forcats")

for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(tibble)
  library(forcats)
})

# Namespace safety: force dplyr verbs in case another package masks them.
filter <- dplyr::filter
select <- dplyr::select
mutate <- dplyr::mutate
arrange <- dplyr::arrange
left_join <- dplyr::left_join
case_when <- dplyr::case_when

# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

analysis_dir <- sd_analysis_dir("CP22")
rdata_dir <- file.path(analysis_dir, "11_RData")
out_dir <- file.path(analysis_dir, "09_main_figures_bracketed")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

step3b_rdata <- file.path(rdata_dir, "SDY2583_CP22_STEP3B_age_sex_adjusted_statistics.RData")
step4_rdata <- file.path(rdata_dir, "SDY2583_CP22_STEP4_composite_scores.RData")
step5_rdata <- file.path(rdata_dir, "SDY2583_CP22_STEP5_threshold_sensitivity_TARGETED.RData")

if (!file.exists(step3b_rdata)) stop("Step 3B RData not found: ", step3b_rdata)
if (!file.exists(step4_rdata)) stop("Step 4 RData not found: ", step4_rdata)
if (!file.exists(step5_rdata)) stop("Step 5 RData not found: ", step5_rdata)

load(step3b_rdata)
load(step4_rdata)
load(step5_rdata)

# ------------------------------------------------------------
# 2. Helpers
# ------------------------------------------------------------

format_fdr <- function(x) {
  ifelse(
    is.na(x),
    "FDR=NA",
    ifelse(x < 0.001, "FDR<0.001", paste0("FDR=", formatC(x, digits = 3, format = "f")))
  )
}

pretty_score <- function(x) {
  x %>%
    str_remove("^CP22_") %>%
    str_remove("_score$") %>%
    str_replace_all("_", " ") %>%
    str_to_sentence()
}

pretty_feature <- function(x) {
  x %>%
    str_replace_all("^pct_", "") %>%
    str_replace_all("_within_b$", " within B") %>%
    str_replace_all("_within_switched_memory_like$", " within switched memory-like") %>%
    str_replace_all("^median_", "median ") %>%
    str_replace_all("_", " ") %>%
    str_to_sentence()
}

first_existing_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

infer_cp22_module <- function(x) {
  dplyr::case_when(
    x %in% c("pct_cd19_b_total", "pct_cd19_b_within_dump_low", "median_CD19_in_B") ~ "B_cell_composition",
    grepl("naive|unswitched|switched_memory|double_negative|igd|cd27", x, ignore.case = TRUE) ~ "core_B_cell_states",
    grepl("igm|iga|igg|IgM|IgA|IgG", x, ignore.case = TRUE) ~ "immunoglobulin_isotype_B_cells",
    grepl("cd38|cd138|plasmablast|plasma", x, ignore.case = TRUE) ~ "plasmablast_plasma_cell_axis",
    grepl("cd10|cd24|cd39|transitional|regulatory", x, ignore.case = TRUE) ~ "transitional_regulatory_like_B_cells",
    TRUE ~ "selected_features"
  )
}

standardize_feature_results <- function(obj) {

  df <- obj

  id_col <- first_existing_col(df, c("feature", "variable", "outcome", "feature_name", "marker_feature"))
  if (is.na(id_col)) {
    cat("\nAvailable columns in primary feature result object:\n")
    print(names(df))
    stop("Feature identifier column bulunamadı. Beklenenlerden biri: feature, variable, outcome, feature_name.")
  }

  beta_col <- first_existing_col(df, c("beta_cancer_vs_healthy", "estimate", "beta", "estimate_cancer_vs_healthy"))
  low_col  <- first_existing_col(df, c("conf_low", "conf.low", "ci_low", "lower_ci"))
  high_col <- first_existing_col(df, c("conf_high", "conf.high", "ci_high", "upper_ci"))
  fdr_col  <- first_existing_col(df, c("FDR_global", "fdr_global", "q_value", "padj", "p_adj_global"))
  p_col    <- first_existing_col(df, c("p_value", "p.value", "p", "P.Value"))

  if (is.na(beta_col) || is.na(low_col) || is.na(high_col)) {
    cat("\nAvailable columns in primary feature result object:\n")
    print(names(df))
    stop("Beta/CI kolonları bulunamadı.")
  }

  out <- df %>%
    mutate(
      feature = as.character(.data[[id_col]]),
      beta_cancer_vs_healthy = suppressWarnings(as.numeric(.data[[beta_col]])),
      conf_low = suppressWarnings(as.numeric(.data[[low_col]])),
      conf_high = suppressWarnings(as.numeric(.data[[high_col]])),
      p_value = if (!is.na(p_col)) suppressWarnings(as.numeric(.data[[p_col]])) else NA_real_,
      FDR_global = if (!is.na(fdr_col)) suppressWarnings(as.numeric(.data[[fdr_col]])) else NA_real_
    )

  if (!("module" %in% names(out))) {
    out$module <- infer_cp22_module(out$feature)
  }

  out
}

# Step 3B object-name compatibility
if (exists("primary_results")) {
  primary_results_norm <- standardize_feature_results(primary_results)
} else if (exists("primary_feature_results")) {
  primary_results_norm <- standardize_feature_results(primary_feature_results)
} else if (exists("age_sex_adjusted_results")) {
  primary_results_norm <- standardize_feature_results(age_sex_adjusted_results)
} else {
  cat("\nObjects loaded from Step 3B RData:\n")
  print(ls())
  stop("Feature-level primary results object bulunamadı: primary_results / primary_feature_results / age_sex_adjusted_results")
}

# ------------------------------------------------------------
# 3. Figure 1: Composite-score forest plot
# ------------------------------------------------------------

composite_plot_df <- primary_composite_results %>%
  mutate(
    score_label = pretty_score(composite_score),
    fdr_label = format_fdr(FDR_global),
    score_label = paste0(score_label, "\n", fdr_label),
    score_label = fct_reorder(score_label, beta_cancer_vs_healthy)
  )

p_composite <- ggplot(
  composite_plot_df,
  aes(
    x = beta_cancer_vs_healthy,
    y = score_label
  )
) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbarh(
    aes(xmin = conf_low, xmax = conf_high),
    height = 0.15
  ) +
  geom_point(size = 2.5) +
  labs(
    title = "CP22 composite humoral B-cell remodeling scores",
    subtitle = "Age- and sex-adjusted cancer vs healthy-control effects",
    x = "Adjusted beta: cancer vs healthy control",
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 9)
  )

ggsave(
  file.path(out_dir, "Figure_CP22_composite_humoral_B_cell_remodeling_forest.png"),
  p_composite,
  width = 9,
  height = 5.5,
  dpi = 300
)

ggsave(
  file.path(out_dir, "Figure_CP22_composite_humoral_B_cell_remodeling_forest.pdf"),
  p_composite,
  width = 9,
  height = 5.5
)

# ------------------------------------------------------------
# 4. Figure 2: Feature-level robust B-cell remodeling
# ------------------------------------------------------------

feature_priority <- c(
  "pct_cd19_b_total",
  "pct_cd19_b_within_dump_low",
  "pct_unswitched_memory_like_within_b",
  "pct_switched_memory_like_within_b",
  "pct_igd_pos_within_b",
  "pct_igm_pos_within_b",
  "pct_igm_unswitched_memory_like_within_b",
  "pct_igg_within_switched_memory_like",
  "pct_iga_within_switched_memory_like",
  "pct_iga_switched_memory_like_within_b",
  "pct_iga_igg_double_negative_switched_like_within_b",
  "pct_class_switched_plasmablast_like_within_b",
  "pct_plasmablast_like_within_b",
  "pct_cd38high_within_b",
  "pct_cd10_cd24_cd38_transitional_like_within_b",
  "pct_cd10_pos_within_b",
  "pct_immature_transitional_like_within_b",
  "pct_cd39_pos_within_b",
  "pct_cd39_regulatory_like_within_b",
  "pct_cd39_cd24_regulatory_like_within_b"
)

feature_plot_df <- primary_results_norm %>%
  filter(feature %in% feature_priority) %>%
  left_join(
    threshold_robustness_summary %>%
      select(variable, robustness_class),
    by = c("feature" = "variable")
  ) %>%
  mutate(
    feature_label = pretty_feature(feature),
    feature_label = paste0(feature_label, "\n", format_fdr(FDR_global)),
    feature_label = factor(feature_label, levels = rev(pretty_feature(feature_priority))),
    module = str_replace_all(module, "_", " "),
    robustness_class = ifelse(is.na(robustness_class), "not tested", robustness_class),
    robust_shape = case_when(
      robustness_class == "direction_and_global_FDR_preserved_all_sets" ~ "Threshold robust",
      robustness_class == "direction_preserved_FDR_not_all_sets" ~ "Direction preserved only",
      TRUE ~ "Threshold-sensitive"
    )
  ) %>%
  arrange(module, FDR_global)

# Preserve only labels present
feature_plot_df <- feature_plot_df %>%
  mutate(
    feature_label = fct_reorder(feature_label, beta_cancer_vs_healthy)
  )

p_feature <- ggplot(
  feature_plot_df,
  aes(
    x = beta_cancer_vs_healthy,
    y = feature_label
  )
) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbarh(
    aes(xmin = conf_low, xmax = conf_high),
    height = 0.12
  ) +
  geom_point(aes(shape = robust_shape), size = 2) +
  facet_grid(
    module ~ .,
    scales = "free_y",
    space = "free_y"
  ) +
  labs(
    title = "CP22 feature-level humoral B-cell remodeling",
    subtitle = "Selected age- and sex-adjusted features; shape indicates threshold robustness",
    x = "Adjusted beta: cancer vs healthy control",
    y = NULL,
    shape = "Threshold sensitivity"
  ) +
  theme_bw(base_size = 9) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text.y = element_text(angle = 0),
    axis.text.y = element_text(size = 7)
  )

ggsave(
  file.path(out_dir, "Figure_CP22_feature_level_humoral_B_cell_remodeling_forest.png"),
  p_feature,
  width = 10,
  height = 12,
  dpi = 300
)

ggsave(
  file.path(out_dir, "Figure_CP22_feature_level_humoral_B_cell_remodeling_forest.pdf"),
  p_feature,
  width = 10,
  height = 12
)

# ------------------------------------------------------------
# 5. Figure 3: Threshold-sensitivity plot
# ------------------------------------------------------------

threshold_selected <- threshold_robustness_summary %>%
  filter(
    variable %in% c(
      "CP22_integrated_humoral_B_cell_remodeling_score",
      "CP22_Ig_isotype_repatterning_score",
      "CP22_core_switched_memory_repatterning_score",
      "CP22_transitional_immature_B_enrichment_score",
      "CP22_plasmablast_activation_enrichment_score",
      "CP22_B_cell_composition_depletion_score",
      "CP22_CD39_regulatory_like_attenuation_score",
      feature_priority
    )
  ) %>%
  arrange(main_FDR_global) %>%
  slice_head(n = 28) %>%
  pull(variable)

threshold_plot_df <- threshold_model_results %>%
  filter(variable %in% threshold_selected) %>%
  mutate(
    variable_label = ifelse(
      str_detect(variable, "^CP22_"),
      pretty_score(variable),
      pretty_feature(variable)
    ),
    variable_label = fct_reorder(variable_label, beta_cancer_vs_healthy),
    threshold_set = factor(threshold_set, levels = c("permissive", "main", "stringent")),
    significant = ifelse(FDR_global < 0.05, "Global FDR < 0.05", "Not global FDR < 0.05")
  )

p_threshold <- ggplot(
  threshold_plot_df,
  aes(
    x = threshold_set,
    y = variable_label
  )
) +
  geom_point(
    aes(
      size = abs(beta_cancer_vs_healthy),
      shape = significant
    )
  ) +
  geom_line(
    aes(group = variable_label),
    linewidth = 0.25
  ) +
  labs(
    title = "CP22 targeted threshold-sensitivity analysis",
    subtitle = "Point size reflects absolute adjusted beta; shape indicates global FDR significance",
    x = "Threshold configuration",
    y = NULL,
    size = "|Beta|",
    shape = "Significance"
  ) +
  theme_bw(base_size = 9) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 7)
  )

ggsave(
  file.path(out_dir, "Figure_CP22_threshold_sensitivity_summary.png"),
  p_threshold,
  width = 9,
  height = 10,
  dpi = 300
)

ggsave(
  file.path(out_dir, "Figure_CP22_threshold_sensitivity_summary.pdf"),
  p_threshold,
  width = 9,
  height = 10
)

# ------------------------------------------------------------
# 6. Figure 4: CP22 clinical annotation overview, if Step 6 exists
# ------------------------------------------------------------

clinical_rdata <- path.expand(
  file.path(sd_integrated_dir(), "RData", "SDY2583_CP22_clinical_annotation_STEP6.RData")
)

if (file.exists(clinical_rdata)) {

  load(clinical_rdata)

  clinical_parts <- list()

  if (exists("therapy_status_binary_results") && nrow(therapy_status_binary_results) > 0) {
    clinical_parts[["therapy_status"]] <- therapy_status_binary_results %>%
      transmute(
        clinical_layer = "Ongoing active treatment",
        score = score,
        beta = beta_ongoing_vs_no_ongoing,
        conf_low = conf_low,
        conf_high = conf_high,
        p_value = p_value,
        FDR_global = FDR_global,
        comparison = "ongoing vs no ongoing"
      )
  }

  if (exists("active_modality_results") && nrow(active_modality_results) > 0) {
    clinical_parts[["modality"]] <- active_modality_results %>%
      transmute(
        clinical_layer = paste0("Modality: ", exposure),
        score = score,
        beta = beta_exposed_vs_unexposed,
        conf_low = conf_low,
        conf_high = conf_high,
        p_value = p_value,
        FDR_global = FDR_global,
        comparison = "exposed vs unexposed"
      )
  }

  if (exists("therapy_line_binary_results") && nrow(therapy_line_binary_results) > 0) {
    clinical_parts[["line"]] <- therapy_line_binary_results %>%
      transmute(
        clinical_layer = "Therapy line",
        score = score,
        beta = beta_later_vs_first,
        conf_low = conf_low,
        conf_high = conf_high,
        p_value = p_value,
        FDR_global = FDR_global,
        comparison = "later vs first"
      )
  }

  if (exists("time_from_start_results") && nrow(time_from_start_results) > 0) {
    clinical_parts[["time"]] <- time_from_start_results %>%
      transmute(
        clinical_layer = "Time from treatment start",
        score = score,
        beta = beta_per_log_time_z,
        conf_low = conf_low,
        conf_high = conf_high,
        p_value = p_value,
        FDR_global = FDR_global,
        comparison = "per log-time z"
      )
  }

  clinical_plot_df <- bind_rows(clinical_parts) %>%
    filter(!is.na(beta)) %>%
    arrange(FDR_global, p_value) %>%
    slice_head(n = 35) %>%
    mutate(
      score_label = pretty_score(score),
      row_label = paste0(clinical_layer, "\n", score_label, "\n", format_fdr(FDR_global)),
      row_label = fct_reorder(row_label, beta),
      significance = ifelse(FDR_global < 0.05, "Global FDR < 0.05", "Not global FDR < 0.05")
    )

  if (nrow(clinical_plot_df) > 0) {

    p_clinical <- ggplot(
      clinical_plot_df,
      aes(
        x = beta,
        y = row_label
      )
    ) +
      geom_vline(xintercept = 0, linetype = "dashed") +
      geom_errorbarh(
        aes(xmin = conf_low, xmax = conf_high),
        height = 0.12
      ) +
      geom_point(aes(shape = significance), size = 2) +
      labs(
        title = "CP22 clinical annotation overview",
        subtitle = "Top-ranked clinical associations among cancer patients; exploratory unless globally FDR-significant",
        x = "Adjusted beta",
        y = NULL,
        shape = "Significance"
      ) +
      theme_bw(base_size = 9) +
      theme(
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 7)
      )

    ggsave(
      file.path(out_dir, "Figure_CP22_clinical_annotation_overview.png"),
      p_clinical,
      width = 10,
      height = 10,
      dpi = 300
    )

    ggsave(
      file.path(out_dir, "Figure_CP22_clinical_annotation_overview.pdf"),
      p_clinical,
      width = 10,
      height = 10
    )
  }
}

# ------------------------------------------------------------
# 7. Save figure-source tables
# ------------------------------------------------------------

readr::write_csv(
  composite_plot_df,
  file.path(out_dir, "FigureSource_CP22_composite_forest.csv")
)

readr::write_csv(
  feature_plot_df,
  file.path(out_dir, "FigureSource_CP22_feature_level_forest.csv")
)

readr::write_csv(
  threshold_plot_df,
  file.path(out_dir, "FigureSource_CP22_threshold_sensitivity.csv")
)

cat("\n============================================================\n")
cat("SDY2583 CP22 STEP 7A COMPLETE: QUANTITATIVE FIGURES\n")
cat("============================================================\n")

cat("\nFigures saved in:\n")
print(out_dir)

cat("\nMain figure files:\n")
print(list.files(out_dir, pattern = "Figure_CP22_.*\\.(png|pdf)$", full.names = FALSE))

cat("============================================================\n")
