# ============================================================
# SDY2583 CP22
# STEP 4 SAFE: Composite score construction and modeling
# B-cell / humoral remodeling panel
#
# Composite scores:
# 1. B-cell composition depletion
# 2. Core switched-memory repatterning
# 3. Immunoglobulin-isotype repatterning
# 4. Plasmablast/activation enrichment
# 5. Transitional/immature B-cell enrichment
# 6. CD39 regulatory-like attenuation
# 7. Integrated humoral B-cell remodeling
# ============================================================

rm(list = ls())
source(file.path(Sys.getenv("SDY2583_REPO_ROOT", unset = "."), "R", "shared", "bootstrap.R"))

# ------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------

cran_pkgs <- c("dplyr", "readr", "stringr", "tibble", "broom", "purrr")

for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(broom)
  library(purrr)
})

# ------------------------------------------------------------
# 2. Paths
# ------------------------------------------------------------

analysis_dir <- sd_analysis_dir("CP22")

step3a_rdata <- file.path(
  analysis_dir,
  "11_RData",
  "SDY2583_CP22_STEP3A_metadata_merge_age_QC.RData"
)

if (!file.exists(step3a_rdata)) {
  stop("CP22 Step 3A RData bulunamadı: ", step3a_rdata)
}

out_dir <- file.path(analysis_dir, "05_composite_scores")
rdata_dir <- file.path(analysis_dir, "11_RData")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rdata_dir, recursive = TRUE, showWarnings = FALSE)

load(step3a_rdata)

if (!exists("cp22_analysis_data")) {
  stop("Step 3A RData içinde cp22_analysis_data bulunamadı.")
}

analysis_df <- cp22_analysis_data

# ------------------------------------------------------------
# 3. Composite score definitions
# direction:
# +1 means higher raw feature contributes to higher composite
# -1 means lower raw feature contributes to higher composite
# Higher composite score is oriented toward the cancer-associated direction
# observed in Step 3B primary models.
# ------------------------------------------------------------

composite_definitions <- tibble::tribble(
  ~composite_score, ~feature, ~direction, ~component_rationale,

  # B-cell composition depletion
  "CP22_B_cell_composition_depletion_score", "pct_cd19_b_total", -1, "lower total CD19+ B-cell representation",
  "CP22_B_cell_composition_depletion_score", "pct_cd19_b_within_dump_low", -1, "lower CD19+ B-cell representation within dump-low compartment",

  # Core switched-memory repatterning
  "CP22_core_switched_memory_repatterning_score", "pct_switched_memory_like_within_b", 1, "higher switched memory-like B-cell fraction",
  "CP22_core_switched_memory_repatterning_score", "pct_unswitched_memory_like_within_b", -1, "lower unswitched memory-like B-cell fraction",
  "CP22_core_switched_memory_repatterning_score", "pct_igd_pos_within_b", -1, "lower IgD+ B-cell fraction",
  "CP22_core_switched_memory_repatterning_score", "pct_igm_unswitched_memory_like_within_b", -1, "lower IgM+ unswitched memory-like fraction",
  "CP22_core_switched_memory_repatterning_score", "pct_double_negative_like_within_b", 1, "higher double-negative-like B-cell fraction",

  # Immunoglobulin-isotype repatterning
  "CP22_Ig_isotype_repatterning_score", "pct_igg_within_switched_memory_like", -1, "lower IgG fraction within switched memory-like B cells",
  "CP22_Ig_isotype_repatterning_score", "pct_iga_within_switched_memory_like", -1, "lower IgA fraction within switched memory-like B cells",
  "CP22_Ig_isotype_repatterning_score", "pct_igm_unswitched_memory_like_within_b", -1, "lower IgM+ unswitched memory-like fraction",
  "CP22_Ig_isotype_repatterning_score", "pct_igm_pos_within_b", -1, "lower IgM+ B-cell fraction",
  "CP22_Ig_isotype_repatterning_score", "pct_iga_switched_memory_like_within_b", 1, "higher IgA+ switched memory-like fraction within B cells",
  "CP22_Ig_isotype_repatterning_score", "pct_iga_igg_double_negative_switched_like_within_b", 1, "higher IgA-IgG- switched-memory-like fraction within B cells",
  "CP22_Ig_isotype_repatterning_score", "pct_iga_igg_double_negative_within_switched_memory_like", 1, "higher IgA-IgG- fraction within switched memory-like B cells",
  "CP22_Ig_isotype_repatterning_score", "median_IgM_in_B", -1, "lower IgM median intensity within B cells",
  "CP22_Ig_isotype_repatterning_score", "median_IgG_in_B", -1, "lower IgG median intensity within B cells",

  # Plasmablast / activation enrichment
  "CP22_plasmablast_activation_enrichment_score", "pct_class_switched_plasmablast_like_within_b", 1, "higher class-switched plasmablast-like fraction",
  "CP22_plasmablast_activation_enrichment_score", "pct_plasmablast_like_within_b", 1, "higher CD38-high CD27+ plasmablast-like fraction",
  "CP22_plasmablast_activation_enrichment_score", "pct_iga_plasmablast_like_within_b", 1, "higher IgA+ plasmablast-like fraction",
  "CP22_plasmablast_activation_enrichment_score", "pct_cd38high_within_b", 1, "higher CD38-high B-cell fraction",
  "CP22_plasmablast_activation_enrichment_score", "pct_cd38_pos_within_b", 1, "higher CD38+ activated/plasmablast-like B-cell fraction",
  "CP22_plasmablast_activation_enrichment_score", "median_CD38_in_B", 1, "higher CD38 median intensity within B cells",

  # Transitional / immature B-cell enrichment
  "CP22_transitional_immature_B_enrichment_score", "pct_cd10_cd24_cd38_transitional_like_within_b", 1, "higher CD10+CD24+CD38+ transitional-like fraction",
  "CP22_transitional_immature_B_enrichment_score", "pct_cd10_pos_within_b", 1, "higher CD10+ immature/transitional-like fraction",
  "CP22_transitional_immature_B_enrichment_score", "pct_immature_transitional_like_within_b", 1, "higher CD10+ immature/transitional-like fraction",
  "CP22_transitional_immature_B_enrichment_score", "median_CD10_in_B", 1, "higher CD10 median intensity within B cells",
  "CP22_transitional_immature_B_enrichment_score", "median_CD24_in_B", 1, "higher CD24 median intensity within B cells",

  # CD39 regulatory-like attenuation
  "CP22_CD39_regulatory_like_attenuation_score", "pct_cd39_pos_within_b", -1, "lower CD39+ B-cell fraction",
  "CP22_CD39_regulatory_like_attenuation_score", "pct_cd39_regulatory_like_within_b", -1, "lower CD39+ regulatory-like B-cell fraction",
  "CP22_CD39_regulatory_like_attenuation_score", "pct_cd39_cd24_regulatory_like_within_b", -1, "lower CD39+CD24+ regulatory-like B-cell fraction",
  "CP22_CD39_regulatory_like_attenuation_score", "pct_cd39_cd24_cd38_regulatory_transitional_like_within_b", -1, "lower CD39+CD24+CD38+ regulatory/transitional-like fraction",
  "CP22_CD39_regulatory_like_attenuation_score", "median_CD39_in_B", -1, "lower CD39 median intensity within B cells",

  # Integrated humoral remodeling
  "CP22_integrated_humoral_B_cell_remodeling_score", "pct_cd19_b_total", -1, "lower total CD19+ B-cell representation",
  "CP22_integrated_humoral_B_cell_remodeling_score", "pct_cd19_b_within_dump_low", -1, "lower CD19+ B-cell representation within dump-low compartment",
  "CP22_integrated_humoral_B_cell_remodeling_score", "pct_switched_memory_like_within_b", 1, "higher switched memory-like B-cell fraction",
  "CP22_integrated_humoral_B_cell_remodeling_score", "pct_unswitched_memory_like_within_b", -1, "lower unswitched memory-like B-cell fraction",
  "CP22_integrated_humoral_B_cell_remodeling_score", "pct_igd_pos_within_b", -1, "lower IgD+ B-cell fraction",
  "CP22_integrated_humoral_B_cell_remodeling_score", "pct_igm_unswitched_memory_like_within_b", -1, "lower IgM+ unswitched memory-like fraction",
  "CP22_integrated_humoral_B_cell_remodeling_score", "pct_igg_within_switched_memory_like", -1, "lower IgG fraction within switched memory-like B cells",
  "CP22_integrated_humoral_B_cell_remodeling_score", "pct_iga_within_switched_memory_like", -1, "lower IgA fraction within switched memory-like B cells",
  "CP22_integrated_humoral_B_cell_remodeling_score", "pct_iga_switched_memory_like_within_b", 1, "higher IgA+ switched memory-like fraction within B cells",
  "CP22_integrated_humoral_B_cell_remodeling_score", "pct_class_switched_plasmablast_like_within_b", 1, "higher class-switched plasmablast-like fraction",
  "CP22_integrated_humoral_B_cell_remodeling_score", "pct_plasmablast_like_within_b", 1, "higher plasmablast-like fraction",
  "CP22_integrated_humoral_B_cell_remodeling_score", "pct_cd38high_within_b", 1, "higher CD38-high B-cell fraction",
  "CP22_integrated_humoral_B_cell_remodeling_score", "pct_cd10_cd24_cd38_transitional_like_within_b", 1, "higher transitional-like B-cell fraction",
  "CP22_integrated_humoral_B_cell_remodeling_score", "pct_cd10_pos_within_b", 1, "higher CD10+ immature/transitional-like fraction",
  "CP22_integrated_humoral_B_cell_remodeling_score", "pct_cd39_pos_within_b", -1, "lower CD39+ B-cell fraction",
  "CP22_integrated_humoral_B_cell_remodeling_score", "pct_cd39_regulatory_like_within_b", -1, "lower CD39+ regulatory-like B-cell fraction"
)

# Keep only features available in analysis_df
composite_definitions <- composite_definitions %>%
  dplyr::filter(feature %in% names(analysis_df))

# ------------------------------------------------------------
# 4. Helper functions
# ------------------------------------------------------------

zscore <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  mu <- mean(x, na.rm = TRUE)
  sig <- stats::sd(x, na.rm = TRUE)
  if (is.na(sig) || sig == 0) {
    return(rep(NA_real_, length(x)))
  }
  (x - mu) / sig
}

make_composite <- function(df, comp_name, defs) {

  defs_comp <- defs %>% dplyr::filter(composite_score == comp_name)

  if (nrow(defs_comp) == 0) {
    return(rep(NA_real_, nrow(df)))
  }

  component_mat <- matrix(NA_real_, nrow = nrow(df), ncol = nrow(defs_comp))

  for (i in seq_len(nrow(defs_comp))) {
    f <- defs_comp$feature[i]
    direction <- defs_comp$direction[i]
    component_mat[, i] <- zscore(df[[f]]) * direction
  }

  out <- rowMeans(component_mat, na.rm = TRUE)
  out[is.nan(out)] <- NA_real_
  out
}

run_lm_score <- function(df, score_name, model_label = "primary") {

  model_df <- df %>%
    dplyr::transmute(
      value = suppressWarnings(as.numeric(.data[[score_name]])),
      disease_group = disease_group,
      age_for_model = suppressWarnings(as.numeric(age_for_model)),
      sex = sex
    ) %>%
    dplyr::filter(
      !is.na(value),
      !is.na(disease_group),
      !is.na(age_for_model),
      !is.na(sex)
    ) %>%
    dplyr::mutate(
      disease_group = factor(
        as.character(disease_group),
        levels = c("Healthy control", "Cancer patient")
      ),
      sex = droplevels(factor(as.character(sex)))
    )

  n_model <- nrow(model_df)
  n_healthy <- sum(model_df$disease_group == "Healthy control")
  n_cancer <- sum(model_df$disease_group == "Cancer patient")

  if (n_model < 50 || n_healthy < 20 || n_cancer < 20) {
    return(NULL)
  }

  fit <- tryCatch(
    stats::lm(value ~ disease_group + age_for_model + sex, data = model_df),
    error = function(e) NULL
  )

  if (is.null(fit)) return(NULL)

  tt <- tryCatch(
    broom::tidy(fit, conf.int = TRUE),
    error = function(e) NULL
  )

  if (is.null(tt)) return(NULL)

  disease_term <- "disease_groupCancer patient"
  if (!(disease_term %in% tt$term)) return(NULL)

  out <- tt %>% dplyr::filter(term == disease_term)

  desc <- model_df %>%
    dplyr::group_by(disease_group) %>%
    dplyr::summarise(
      n = dplyr::n(),
      mean = mean(value, na.rm = TRUE),
      sd = sd(value, na.rm = TRUE),
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

  tibble::tibble(
    model_label = model_label,
    composite_score = score_name,
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
    direction = dplyr::case_when(
      out$estimate[1] > 0 ~ "higher_in_cancer",
      out$estimate[1] < 0 ~ "lower_in_cancer",
      TRUE ~ "no_direction"
    )
  )
}

# ------------------------------------------------------------
# 5. Create composite scores
# ------------------------------------------------------------

composite_names <- unique(composite_definitions$composite_score)

for (cs in composite_names) {
  analysis_df[[cs]] <- make_composite(analysis_df, cs, composite_definitions)
}

composite_availability <- tibble::tibble(
  composite_score = composite_names,
  n_components = vapply(
    composite_names,
    function(cs) sum(composite_definitions$composite_score == cs),
    numeric(1)
  ),
  n_nonmissing = vapply(
    composite_names,
    function(cs) sum(!is.na(analysis_df[[cs]])),
    numeric(1)
  )
)

# ------------------------------------------------------------
# 6. Model datasets
# ------------------------------------------------------------

primary_df <- analysis_df %>%
  dplyr::filter(
    feature_ok == TRUE,
    !is.na(disease_group),
    !is.na(age_for_model),
    !is.na(sex)
  )

binary_sex_df <- primary_df %>%
  dplyr::filter(as.character(sex) %in% c("Female", "Male")) %>%
  dplyr::mutate(sex = factor(as.character(sex), levels = c("Female", "Male")))

event_qc_df <- primary_df %>%
  dplyr::filter(!is.na(n_cd19_b), n_cd19_b >= 300)

standard_dump_df <- primary_df %>%
  dplyr::filter(dump_marker_variant == "Viability_CD3_CD7_CD13")

# ------------------------------------------------------------
# 7. Run composite models
# ------------------------------------------------------------

primary_composite_results <- dplyr::bind_rows(
  lapply(composite_names, function(cs) {
    run_lm_score(primary_df, cs, model_label = "primary_age_sex_adjusted")
  })
) %>%
  dplyr::mutate(FDR_global = p.adjust(p_value, method = "BH")) %>%
  dplyr::arrange(FDR_global, p_value)

binary_sex_composite_results <- dplyr::bind_rows(
  lapply(composite_names, function(cs) {
    run_lm_score(binary_sex_df, cs, model_label = "binary_sex_sensitivity")
  })
) %>%
  dplyr::mutate(FDR_global = p.adjust(p_value, method = "BH")) %>%
  dplyr::arrange(FDR_global, p_value)

event_qc_composite_results <- dplyr::bind_rows(
  lapply(composite_names, function(cs) {
    run_lm_score(event_qc_df, cs, model_label = "event_QC_n_cd19_b_ge_300")
  })
) %>%
  dplyr::mutate(FDR_global = p.adjust(p_value, method = "BH")) %>%
  dplyr::arrange(FDR_global, p_value)

standard_dump_composite_results <- dplyr::bind_rows(
  lapply(composite_names, function(cs) {
    run_lm_score(standard_dump_df, cs, model_label = "standard_dump_marker_only")
  })
) %>%
  dplyr::mutate(FDR_global = p.adjust(p_value, method = "BH")) %>%
  dplyr::arrange(FDR_global, p_value)

# ------------------------------------------------------------
# 8. Sensitivity summary
# ------------------------------------------------------------

composite_direction_sensitivity_summary <- primary_composite_results %>%
  dplyr::select(
    composite_score,
    primary_beta = beta_cancer_vs_healthy,
    primary_direction = direction,
    primary_FDR_global = FDR_global
  ) %>%
  dplyr::left_join(
    binary_sex_composite_results %>%
      dplyr::select(
        composite_score,
        binary_beta = beta_cancer_vs_healthy,
        binary_direction = direction,
        binary_FDR_global = FDR_global
      ),
    by = "composite_score"
  ) %>%
  dplyr::left_join(
    event_qc_composite_results %>%
      dplyr::select(
        composite_score,
        event_QC_beta = beta_cancer_vs_healthy,
        event_QC_direction = direction,
        event_QC_FDR_global = FDR_global
      ),
    by = "composite_score"
  ) %>%
  dplyr::left_join(
    standard_dump_composite_results %>%
      dplyr::select(
        composite_score,
        standard_dump_beta = beta_cancer_vs_healthy,
        standard_dump_direction = direction,
        standard_dump_FDR_global = FDR_global
      ),
    by = "composite_score"
  ) %>%
  dplyr::mutate(
    direction_preserved_binary_sex = primary_direction == binary_direction,
    direction_preserved_event_QC = primary_direction == event_QC_direction,
    direction_preserved_standard_dump = primary_direction == standard_dump_direction,
    direction_preserved_all_sensitivities =
      direction_preserved_binary_sex == TRUE &
      direction_preserved_event_QC == TRUE &
      direction_preserved_standard_dump == TRUE,
    FDR_global_preserved_all_sensitivities =
      primary_FDR_global < 0.05 &
      binary_FDR_global < 0.05 &
      event_QC_FDR_global < 0.05 &
      standard_dump_FDR_global < 0.05
  ) %>%
  dplyr::arrange(primary_FDR_global)

overall_composite_summary <- tibble::tibble(
  n_primary_model_subjects = nrow(primary_df),
  n_primary_healthy = sum(primary_df$disease_group == "Healthy control"),
  n_primary_cancer = sum(primary_df$disease_group == "Cancer patient"),
  n_binary_sex_subjects = nrow(binary_sex_df),
  n_event_QC_subjects = nrow(event_qc_df),
  n_standard_dump_subjects = nrow(standard_dump_df),
  n_composite_scores = length(composite_names),
  n_primary_FDR_lt_0p05 = sum(primary_composite_results$FDR_global < 0.05, na.rm = TRUE),
  n_direction_preserved_all_sensitivities = sum(composite_direction_sensitivity_summary$direction_preserved_all_sensitivities == TRUE, na.rm = TRUE),
  n_FDR_preserved_all_sensitivities = sum(composite_direction_sensitivity_summary$FDR_global_preserved_all_sensitivities == TRUE, na.rm = TRUE)
)

# ------------------------------------------------------------
# 9. Save outputs
# ------------------------------------------------------------

readr::write_csv(
  analysis_df,
  file.path(out_dir, "SDY2583_CP22_analysis_data_with_composite_scores_STEP4.csv")
)

readr::write_csv(
  composite_definitions,
  file.path(out_dir, "SDY2583_CP22_composite_score_definitions_STEP4.csv")
)

readr::write_csv(
  composite_availability,
  file.path(out_dir, "SDY2583_CP22_composite_score_availability_STEP4.csv")
)

readr::write_csv(
  primary_composite_results,
  file.path(out_dir, "SDY2583_CP22_primary_composite_results_STEP4.csv")
)

readr::write_csv(
  binary_sex_composite_results,
  file.path(out_dir, "SDY2583_CP22_binary_sex_composite_results_STEP4.csv")
)

readr::write_csv(
  event_qc_composite_results,
  file.path(out_dir, "SDY2583_CP22_event_QC_composite_results_STEP4.csv")
)

readr::write_csv(
  standard_dump_composite_results,
  file.path(out_dir, "SDY2583_CP22_standard_dump_composite_results_STEP4.csv")
)

readr::write_csv(
  composite_direction_sensitivity_summary,
  file.path(out_dir, "SDY2583_CP22_composite_direction_sensitivity_summary_STEP4.csv")
)

readr::write_csv(
  overall_composite_summary,
  file.path(out_dir, "SDY2583_CP22_overall_composite_summary_STEP4.csv")
)

save(
  analysis_df,
  composite_definitions,
  composite_availability,
  primary_composite_results,
  binary_sex_composite_results,
  event_qc_composite_results,
  standard_dump_composite_results,
  composite_direction_sensitivity_summary,
  overall_composite_summary,
  primary_df,
  binary_sex_df,
  event_qc_df,
  standard_dump_df,
  file = file.path(
    rdata_dir,
    "SDY2583_CP22_STEP4_composite_scores.RData"
  )
)

# Convenience object name for later integration scripts
cp22_scores_data <- analysis_df

save(
  cp22_scores_data,
  analysis_df,
  composite_definitions,
  composite_availability,
  primary_composite_results,
  composite_direction_sensitivity_summary,
  overall_composite_summary,
  file = file.path(
    rdata_dir,
    "SDY2583_CP22_STEP4_composite_scores_for_integration.RData"
  )
)

# ------------------------------------------------------------
# 10. Console output
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SDY2583 CP22 STEP 4 COMPLETE: COMPOSITE SCORES\n")
cat("============================================================\n")

cat("\nComposite availability:\n")
print(as.data.frame(composite_availability), row.names = FALSE)

cat("\nOverall composite summary:\n")
print(as.data.frame(overall_composite_summary), row.names = FALSE)

cat("\nPrimary composite results:\n")
print(as.data.frame(primary_composite_results), row.names = FALSE)

cat("\nComposite sensitivity summary:\n")
print(as.data.frame(composite_direction_sensitivity_summary), row.names = FALSE)

cat("\nFiles saved in:\n")
print(out_dir)

cat("============================================================\n")
