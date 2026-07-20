# ============================================================
# SDY2583 CP10
# STEP 8 SAFE: Manuscript-ready summary tables, figure legends,
# Methods/Results draft text, and supplementary output guide
#
# This step does not re-run analyses.
# It gathers CP10 outputs from Steps 1-7 and writes:
#   1) CP10 manuscript-ready summary tables
#   2) CP10 concise Methods text
#   3) CP10 Results text
#   4) CP10 figure legends
#   5) CP10 supplementary file guide
#
# Output:
#   ~/Desktop/SDY2583_CP10_FULL_850_ANALYSIS/10_manuscript_export
# ============================================================

rm(list = ls())
source(file.path(Sys.getenv("SDY2583_REPO_ROOT", unset = "."), "R", "shared", "bootstrap.R"))

# ------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------

cran_pkgs <- c("dplyr", "readr", "stringr", "tibble", "tidyr")

for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
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
main_fig_dir <- file.path(analysis_dir, "08_main_figures")
flow_fig_dir <- file.path(analysis_dir, "09_flowjo_style_representative_gating")

integrated_matrix_dir <- sd_integrated_dir()
clinical_rdata_dir <- file.path(integrated_matrix_dir, "RData")

out_dir <- file.path(analysis_dir, "10_manuscript_export")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

required_rdata <- c(
  step1 = file.path(rdata_dir, "SDY2583_CP10_STEP1_fcs_inventory_marker_QC.RData"),
  step3a = file.path(rdata_dir, "SDY2583_CP10_STEP3A_metadata_merge_age_QC.RData"),
  step3b = file.path(rdata_dir, "SDY2583_CP10_STEP3B_age_sex_adjusted_statistics.RData"),
  step4 = file.path(rdata_dir, "SDY2583_CP10_STEP4_composite_scores.RData"),
  step4b = file.path(rdata_dir, "SDY2583_CP10_STEP4B_age_sensitivity_caliper_interaction.RData"),
  step5 = file.path(rdata_dir, "SDY2583_CP10_STEP5_threshold_sensitivity_TARGETED.RData"),
  step6 = file.path(clinical_rdata_dir, "SDY2583_CP10_clinical_annotation_STEP6.RData"),
  step7b = file.path(rdata_dir, "SDY2583_CP10_STEP7B_flowjo_style_representative_gating.RData")
)

missing <- required_rdata[!file.exists(required_rdata)]

if (length(missing) > 0) {
  cat("\nMissing RData files:\n")
  print(missing)
  stop("Eksik RData dosyaları var. Önce ilgili adımların tamamlandığından emin olun.")
}

for (fp in required_rdata) {
  load(fp)
}

# Fallback: Step 7A manifest is CSV.
fig7a_manifest_csv <- file.path(main_fig_dir, "SDY2583_CP10_STEP7A_main_figure_manifest.csv")
if (file.exists(fig7a_manifest_csv)) {
  step7a_figure_manifest <- readr::read_csv(fig7a_manifest_csv, show_col_types = FALSE)
} else {
  step7a_figure_manifest <- tibble()
}

# ------------------------------------------------------------
# 3. Safety checks and helpers
# ------------------------------------------------------------

required_objects <- c(
  "cp10_analysis_data",
  "primary_results",
  "sensitivity_summary",
  "primary_composite_results",
  "composite_sensitivity_summary",
  "matched_summary",
  "matched_5y_results",
  "matched_10y_results",
  "interaction_results",
  "matched_direction_summary",
  "threshold_robustness_summary",
  "threshold_extraction_summary",
  "clinical_model_summary",
  "cancer_subgroup_omnibus_results",
  "cancer_subgroup_one_vs_rest_results",
  "therapy_status_binary_results",
  "active_modality_results",
  "therapy_line_binary_results",
  "time_from_start_results",
  "representative_samples",
  "gate_summary",
  "figure_manifest"
)

for (oo in required_objects) {
  if (!exists(oo)) {
    warning("Object not found after RData load: ", oo)
  }
}

fmt_num <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", formatC(x, format = "f", digits = digits))
}

fmt_p <- function(x) {
  ifelse(
    is.na(x), "NA",
    ifelse(x < 0.001, formatC(x, format = "e", digits = 2), signif(x, 3))
  )
}

fmt_pct <- function(x) {
  ifelse(is.na(x), "NA", paste0(formatC(x, format = "f", digits = 1), "%"))
}

safe_pull1 <- function(df, col, default = NA) {
  if (!exists(deparse(substitute(df)))) return(default)
  if (!(col %in% names(df))) return(default)
  if (nrow(df) == 0) return(default)
  df[[col]][1]
}

# ------------------------------------------------------------
# 4. Cohort/QC summary
# ------------------------------------------------------------

cohort_summary <- tibble(
  metric = c(
    "FCS files processed",
    "Unique subjects",
    "Healthy controls",
    "Cancer patients",
    "Valid age values",
    "Model-ready age/sex subjects",
    "Model-ready healthy controls",
    "Model-ready cancer patients",
    "Channel-order mismatch files",
    "Marker mismatch files"
  ),
  value = c(
    nrow(cp10_analysis_data),
    n_distinct(cp10_analysis_data$subject_id),
    sum(cp10_analysis_data$disease_group == "Healthy control", na.rm = TRUE),
    sum(cp10_analysis_data$disease_group == "Cancer patient", na.rm = TRUE),
    sum(!is.na(cp10_analysis_data$age_for_model)),
    sum(cp10_analysis_data$model_ready_age_sex == TRUE, na.rm = TRUE),
    sum(cp10_analysis_data$model_ready_age_sex == TRUE & cp10_analysis_data$disease_group == "Healthy control", na.rm = TRUE),
    sum(cp10_analysis_data$model_ready_age_sex == TRUE & cp10_analysis_data$disease_group == "Cancer patient", na.rm = TRUE),
    sum(cp10_analysis_data$channel_order_mismatch_file == TRUE, na.rm = TRUE),
    sum(cp10_analysis_data$marker_mismatch_file == TRUE, na.rm = TRUE)
  )
)

age_summary <- cp10_analysis_data %>%
  group_by(disease_group) %>%
  summarise(
    n_total = n(),
    n_valid_age = sum(!is.na(age_for_model)),
    mean_age = mean(age_for_model, na.rm = TRUE),
    sd_age = sd(age_for_model, na.rm = TRUE),
    median_age = median(age_for_model, na.rm = TRUE),
    q1_age = as.numeric(quantile(age_for_model, 0.25, na.rm = TRUE)),
    q3_age = as.numeric(quantile(age_for_model, 0.75, na.rm = TRUE)),
    .groups = "drop"
  )

sex_summary <- cp10_analysis_data %>%
  count(disease_group, sex, name = "n") %>%
  arrange(disease_group, sex)

# ------------------------------------------------------------
# 5. Primary results tables
# ------------------------------------------------------------

primary_composite_table <- primary_composite_results %>%
  arrange(FDR_global, p_value) %>%
  transmute(
    composite_score,
    n_model,
    n_healthy,
    n_cancer,
    healthy_mean,
    cancer_mean,
    healthy_median,
    cancer_median,
    beta_cancer_vs_healthy,
    conf_low,
    conf_high,
    p_value,
    FDR_global,
    direction
  )

primary_feature_table <- primary_results %>%
  arrange(FDR_global, p_value) %>%
  transmute(
    feature,
    module,
    interpretation,
    n_model,
    n_healthy,
    n_cancer,
    healthy_mean,
    cancer_mean,
    healthy_median,
    cancer_median,
    beta_cancer_vs_healthy,
    conf_low,
    conf_high,
    p_value,
    FDR_global,
    direction
  )

robust_feature_table <- primary_feature_table %>%
  semi_join(
    sensitivity_summary %>%
      filter(FDR_global_preserved_all_sensitivities == TRUE) %>%
      select(feature),
    by = "feature"
  )

composite_sensitivity_table <- composite_sensitivity_summary %>%
  arrange(primary_FDR_global)

# ------------------------------------------------------------
# 6. Age sensitivity / threshold / clinical tables
# ------------------------------------------------------------

matched_table <- matched_direction_summary %>%
  arrange(primary_FDR_global)

interaction_table <- interaction_results %>%
  arrange(FDR_global, p_value)

threshold_summary_table <- threshold_robustness_summary %>%
  arrange(
    variable_type,
    factor(
      robustness_class,
      levels = c(
        "direction_and_global_FDR_preserved_all_sets",
        "direction_preserved_FDR_not_all_sets",
        "direction_not_preserved"
      )
    ),
    FDR_global__main
  )

clinical_subgroup_omnibus_table <- cancer_subgroup_omnibus_results %>%
  arrange(FDR_global, omnibus_p)

clinical_breast_one_vs_rest_table <- cancer_subgroup_one_vs_rest_results %>%
  filter(cancer_subgroup == "Meme") %>%
  arrange(FDR_global, p_value)

therapy_status_table <- therapy_status_binary_results %>%
  arrange(FDR_global, p_value)

active_modality_table <- active_modality_results %>%
  arrange(FDR_global, p_value)

therapy_line_table <- therapy_line_binary_results %>%
  arrange(FDR_global, p_value)

time_from_start_table <- time_from_start_results %>%
  arrange(FDR_global, p_value)

# ------------------------------------------------------------
# 7. Write CSV tables
# ------------------------------------------------------------

write_csv(cohort_summary, file.path(out_dir, "CP10_Table_00_cohort_QC_summary.csv"))
write_csv(age_summary, file.path(out_dir, "CP10_Table_01_age_summary_by_disease.csv"))
write_csv(sex_summary, file.path(out_dir, "CP10_Table_02_sex_summary_by_disease.csv"))
write_csv(primary_composite_table, file.path(out_dir, "CP10_Table_03_primary_composite_results.csv"))
write_csv(primary_feature_table, file.path(out_dir, "CP10_Table_04_primary_feature_results.csv"))
write_csv(robust_feature_table, file.path(out_dir, "CP10_Table_05_robust_feature_results.csv"))
write_csv(composite_sensitivity_table, file.path(out_dir, "CP10_Table_06_composite_sensitivity_summary.csv"))
write_csv(matched_table, file.path(out_dir, "CP10_Table_07_age_matched_sensitivity_summary.csv"))
write_csv(interaction_table, file.path(out_dir, "CP10_Table_08_disease_by_age_interactions.csv"))
write_csv(threshold_summary_table, file.path(out_dir, "CP10_Table_09_threshold_robustness_summary.csv"))
write_csv(threshold_extraction_summary, file.path(out_dir, "CP10_Table_10_threshold_extraction_summary.csv"))
write_csv(clinical_model_summary, file.path(out_dir, "CP10_Table_11_clinical_model_summary.csv"))
write_csv(clinical_subgroup_omnibus_table, file.path(out_dir, "CP10_Table_12_cancer_subgroup_omnibus.csv"))
write_csv(clinical_breast_one_vs_rest_table, file.path(out_dir, "CP10_Table_13_breast_cancer_one_vs_rest.csv"))
write_csv(therapy_status_table, file.path(out_dir, "CP10_Table_14_therapy_status_results.csv"))
write_csv(active_modality_table, file.path(out_dir, "CP10_Table_15_active_treatment_modality_results.csv"))
write_csv(therapy_line_table, file.path(out_dir, "CP10_Table_16_therapy_line_results.csv"))
write_csv(time_from_start_table, file.path(out_dir, "CP10_Table_17_time_from_start_results.csv"))
write_csv(representative_samples, file.path(out_dir, "CP10_Table_18_flow_representative_samples.csv"))
write_csv(gate_summary, file.path(out_dir, "CP10_Table_19_flow_representative_gate_percentages.csv"))

# ------------------------------------------------------------
# 8. Extract key values for text
# ------------------------------------------------------------

get_comp <- function(score) {
  primary_composite_table %>% filter(composite_score == score) %>% slice(1)
}

get_feature <- function(feature_name) {
  primary_feature_table %>% filter(feature == feature_name) %>% slice(1)
}

integrated_comp <- get_comp("CP10_integrated_myeloid_granulocytic_remodeling_score")
gran_comp <- get_comp("CP10_granulocyte_like_enrichment_score")
balance_comp <- get_comp("CP10_myeloid_granulocytic_to_lymphoid_balance_score")
mono_comp <- get_comp("CP10_monocyte_like_HLA_DR_phenotype_remodeling_score")
apc_comp <- get_comp("CP10_APC_like_myeloid_enrichment_score")
ccr3_comp <- get_comp("CP10_CCR3_eosinophil_like_granulocytic_enrichment_score")

feat_cd13_cd66b <- get_feature("pct_cd13_cd66b_gran_like_within_cd45")
feat_cd66b <- get_feature("pct_cd66b_gran_like_within_cd45")
feat_t <- get_feature("pct_t_like_within_cd45")
feat_lymph <- get_feature("pct_lymphoid_like_any_within_cd45")
feat_cd14 <- get_feature("pct_cd14_mono_like_within_cd45")
feat_cd11c_hladr <- get_feature("pct_cd11c_hladr_apc_like_within_cd45")

matched_summary_5 <- matched_summary %>% filter(caliper_years == 5) %>% slice(1)
matched_summary_10 <- matched_summary %>% filter(caliper_years == 10) %>% slice(1)

threshold_counts <- threshold_summary_table %>%
  count(variable_type, robustness_class, name = "n") %>%
  arrange(variable_type, robustness_class)

n_comp_threshold_robust <- threshold_counts %>%
  filter(variable_type == "composite", robustness_class == "direction_and_global_FDR_preserved_all_sets") %>%
  pull(n)
if (length(n_comp_threshold_robust) == 0) n_comp_threshold_robust <- 0

n_feature_threshold_robust <- threshold_counts %>%
  filter(variable_type == "feature", robustness_class == "direction_and_global_FDR_preserved_all_sets") %>%
  pull(n)
if (length(n_feature_threshold_robust) == 0) n_feature_threshold_robust <- 0

breast_top <- clinical_breast_one_vs_rest_table %>%
  filter(score %in% c(
    "CP10_integrated_myeloid_granulocytic_remodeling_score",
    "CP10_monocyte_like_HLA_DR_phenotype_remodeling_score",
    "CP10_granulocyte_like_enrichment_score",
    "CP10_myeloid_granulocytic_to_lymphoid_balance_score"
  )) %>%
  arrange(FDR_global)

# ------------------------------------------------------------
# 9. Manuscript Methods text
# ------------------------------------------------------------

methods_text <- c(
"CP10 myeloid/granulocytic flow cytometry analysis",
"",
"CP10 was analyzed as a myeloid/granulocytic and leukocyte-lineage composition panel. Raw FCS files were read in R using marker-name rather than fixed channel-position mapping, because two CP10 files had shifted channel order with Time occupying the first parameter position. Compensation matrices were extracted from FCS keywords when available and applied before logicle transformation. Fluorescence channels were transformed using a fixed logicle transformation, and all downstream gates were applied on transformed compensated values.",
"",
"Primary gating first defined viable CD45+ leukocytes using a viability-low gate and CD45 positivity. Within viable CD45+ events, threshold-defined phenotype-like fractions were quantified, including CD3+ T-cell-like, CD19+ B-cell-like, CD3−CD56+ NK-like, CD14+ monocyte-like, CD14+HLA-DR-low, CD14+HLA-DR+, CD11c+HLA-DR+ APC-like myeloid, CD13+CD66b+ granulocyte-like, CCR3/CD193+CD66b+ eosinophil-like granulocytic, and CD123+HLA-DR+ DC-like features. Because CP10 does not include a full lineage-resolution panel, these subsets were interpreted as phenotype-like compartments rather than definitive neutrophils, eosinophils, dendritic cells, or myeloid-derived suppressor cells.",
"",
"Feature-level cancer-versus-healthy comparisons were performed using linear models adjusted for age and sex. Age values outside 18–100 years were excluded from modeling. Primary analyses included subjects with valid disease group, age, and sex. Sensitivity analyses included binary-sex restriction, CD45 event-count quality control, and exclusion of the two files with channel-order/marker mismatch. Multiple-testing correction was performed using the Benjamini–Hochberg false discovery rate globally and, for feature-level models, within functional modules.",
"",
"Composite scores were constructed from directionally oriented z-scored features. Higher composite values represented the cancer-associated remodeling direction observed in the primary models. CP10 composites included granulocyte-like enrichment, integrated myeloid/granulocytic remodeling, myeloid/granulocytic-to-lymphoid balance, monocyte-like HLA-DR phenotype remodeling, APC-like myeloid enrichment, and CCR3/eosinophil-like granulocytic enrichment. Composite scores were tested using the same age- and sex-adjusted primary models and sensitivity models.",
"",
"Age sensitivity was assessed using same-sex nearest-age matching with 5-year and 10-year calipers and by testing disease-by-age interaction terms. Targeted threshold sensitivity was performed by re-extracting CP10 features under three gating configurations: main thresholds, permissive thresholds, and stringent thresholds. For positive-marker gates, permissive thresholds were lowered by 0.2 transformed units and stringent thresholds were raised by 0.2 transformed units. For the viability-low gate, the direction was reversed. A result was considered threshold-robust when the cancer-versus-healthy direction and global FDR significance were preserved across all three threshold configurations.",
"",
"Clinical annotation analyses were performed in cancer patients only after integration of CP10 composite scores into the multi-panel clinical immune score matrix. Models examined cancer subgroup, ongoing active treatment status, treatment modality, therapy line, and time from treatment start, adjusting for age, sex, and cancer subgroup where appropriate. These clinical analyses were treated as secondary and hypothesis-generating unless associations survived global FDR correction."
)

# ------------------------------------------------------------
# 10. Manuscript Results text
# ------------------------------------------------------------

results_text <- c(
"CP10 identifies robust cancer-associated myeloid/granulocytic remodeling",
"",
paste0(
"CP10 data were available for ", cohort_summary$value[cohort_summary$metric == "FCS files processed"],
" FCS files from ", cohort_summary$value[cohort_summary$metric == "Unique subjects"],
" unique subjects, including ", cohort_summary$value[cohort_summary$metric == "Healthy controls"],
" healthy controls and ", cohort_summary$value[cohort_summary$metric == "Cancer patients"],
" cancer patients. After age quality control, ", cohort_summary$value[cohort_summary$metric == "Model-ready age/sex subjects"],
" subjects were available for age- and sex-adjusted modeling (",
cohort_summary$value[cohort_summary$metric == "Model-ready healthy controls"], " healthy controls and ",
cohort_summary$value[cohort_summary$metric == "Model-ready cancer patients"], " cancer patients)."
),
"",
paste0(
"At the composite-score level, CP10 showed a strong cancer-associated myeloid/granulocytic remodeling phenotype. The strongest composite was the granulocyte-like enrichment score (beta=",
fmt_num(gran_comp$beta_cancer_vs_healthy), ", 95% CI ",
fmt_num(gran_comp$conf_low), " to ", fmt_num(gran_comp$conf_high),
", FDR=", fmt_p(gran_comp$FDR_global), "), followed by the integrated myeloid/granulocytic remodeling score (beta=",
fmt_num(integrated_comp$beta_cancer_vs_healthy), ", 95% CI ",
fmt_num(integrated_comp$conf_low), " to ", fmt_num(integrated_comp$conf_high),
", FDR=", fmt_p(integrated_comp$FDR_global), ") and the myeloid/granulocytic-to-lymphoid balance score (beta=",
fmt_num(balance_comp$beta_cancer_vs_healthy), ", 95% CI ",
fmt_num(balance_comp$conf_low), " to ", fmt_num(balance_comp$conf_high),
", FDR=", fmt_p(balance_comp$FDR_global), "). Monocyte-like HLA-DR phenotype remodeling and APC-like myeloid enrichment scores were also higher in cancer patients (FDR=",
fmt_p(mono_comp$FDR_global), " and FDR=", fmt_p(apc_comp$FDR_global), ", respectively). In contrast, the CCR3/eosinophil-like granulocytic composite was not significant in the primary composite analysis (FDR=",
fmt_p(ccr3_comp$FDR_global), ")."
),
"",
paste0(
"Feature-level analyses indicated that the CP10 signal was driven primarily by CD13/CD66b granulocyte-like enrichment and a shift in the CD45+ leukocyte compartment from lymphoid-like toward myeloid/granulocytic-like phenotypes. CD13+CD66b+ granulocyte-like cells were higher in cancer patients (healthy mean=",
fmt_num(feat_cd13_cd66b$healthy_mean), ", cancer mean=", fmt_num(feat_cd13_cd66b$cancer_mean),
", beta=", fmt_num(feat_cd13_cd66b$beta_cancer_vs_healthy),
", FDR=", fmt_p(feat_cd13_cd66b$FDR_global), "), as were CD66b+ granulocyte-like cells (beta=",
fmt_num(feat_cd66b$beta_cancer_vs_healthy), ", FDR=", fmt_p(feat_cd66b$FDR_global), "). Conversely, CD3+ T-cell-like and broader CD3/CD19/CD56-defined lymphoid-like fractions were lower in cancer patients (T-cell-like beta=",
fmt_num(feat_t$beta_cancer_vs_healthy), ", FDR=", fmt_p(feat_t$FDR_global),
"; lymphoid-like beta=", fmt_num(feat_lymph$beta_cancer_vs_healthy), ", FDR=", fmt_p(feat_lymph$FDR_global), ")."
),
"",
paste0(
"CD14+ monocyte-like features were also enriched in cancer (beta=",
fmt_num(feat_cd14$beta_cancer_vs_healthy), ", FDR=", fmt_p(feat_cd14$FDR_global),
"). CD14+HLA-DR-low and CD14+HLA-DR+ phenotypes increased within the CD45+ denominator, whereas within-CD14 HLA-DR redistribution was not sufficiently robust to support a definitive per-monocyte HLA-DR suppression claim. Therefore, this axis was interpreted as expansion of CD14/HLA-DR-defined monocyte-like phenotypes rather than direct evidence of monocyte HLA-DR attenuation. CD11c+HLA-DR+ APC-like myeloid features were also higher in cancer (beta=",
fmt_num(feat_cd11c_hladr$beta_cancer_vs_healthy), ", FDR=", fmt_p(feat_cd11c_hladr$FDR_global), ")."
),
"",
paste0(
"Age sensitivity analyses supported the age independence of the CP10 signal. Same-sex age-caliper matching yielded ",
matched_summary_5$n_pairs, " matched pairs with a 5-year caliper (mean absolute age difference ",
fmt_num(matched_summary_5$mean_abs_age_diff), " years) and ",
matched_summary_10$n_pairs, " matched pairs with a 10-year caliper (mean absolute age difference ",
fmt_num(matched_summary_10$mean_abs_age_diff), " years). Five of six composite scores retained global FDR significance in both matched analyses, and all six preserved the cancer-versus-healthy direction. Disease-by-age interaction testing did not identify any globally FDR-significant interaction effects."
),
"",
paste0(
"Targeted threshold sensitivity further supported the robustness of the main CP10 findings. Despite substantial shifts in absolute gated percentages under permissive and stringent thresholds, ",
n_comp_threshold_robust, " of 6 composite scores and ",
n_feature_threshold_robust, " targeted features preserved both direction and global FDR significance across all threshold configurations. The threshold-robust composite axes included granulocyte-like enrichment, integrated myeloid/granulocytic remodeling, myeloid/granulocytic-to-lymphoid balance, monocyte-like HLA-DR phenotype remodeling, and APC-like myeloid enrichment. The CCR3/eosinophil-like granulocytic composite preserved direction but did not preserve global FDR significance across all threshold configurations."
),
"",
"Clinical annotation of CP10 composite scores",
"",
paste0(
"CP10 composite scores were integrated into the multi-panel clinical immune matrix, expanding the matrix to eight panels and 55 immune composite scores. In cancer-only clinical annotation analyses, CP10 scores were not robustly associated with ongoing active treatment status, active treatment modality, therapy line, or time from treatment start after global FDR correction."
),
"",
paste0(
"Cancer-subgroup omnibus testing showed nominal but not globally FDR-significant variation in CP10 myeloid/granulocytic scores. However, one-vs-rest analyses identified a lower CP10 remodeling magnitude in breast cancer compared with the remaining cancer groups. Specifically, breast cancer showed lower integrated myeloid/granulocytic remodeling, monocyte-like HLA-DR phenotype remodeling, granulocyte-like enrichment, and myeloid/granulocytic-to-lymphoid balance scores after age and sex adjustment, with global FDR-significant one-vs-rest contrasts. This finding should be interpreted as tumor-context-dependent variation in the magnitude of CP10 myeloid/granulocytic remodeling rather than as evidence of treatment resistance or a general claim that breast cancer is less inflammatory."
),
"",
"Overall CP10 interpretation",
"",
"Taken together, CP10 provides a robust myeloid/granulocytic axis for the SDY2583 reanalysis. The dominant cancer-associated pattern is CD13/CD66b granulocyte-like enrichment together with a myeloid/granulocytic-to-lymphoid composition shift and expansion of CD14/HLA-DR/CD11c-defined monocyte-like and APC-like myeloid phenotypes. The signal is preserved after age adjustment, age matching, event-count and technical sensitivity analyses, and targeted threshold sensitivity. Clinical annotation suggests that this axis is not primarily explained by treatment status or modality, although its magnitude varies across tumor contexts."
)

# ------------------------------------------------------------
# 11. Figure legends
# ------------------------------------------------------------

figure_legends <- c(
"CP10 Figure Legends",
"",
"Figure 1. CP10 composite scores in healthy controls and cancer patients. Boxplots show directionally oriented CP10 composite scores in model-ready subjects. Higher values indicate the cancer-associated remodeling direction. Composite scores shown include granulocyte-like enrichment, integrated myeloid/granulocytic remodeling, myeloid/granulocytic-to-lymphoid balance, monocyte-like HLA-DR phenotype remodeling, and APC-like myeloid enrichment. Labels report age- and sex-adjusted cancer-versus-healthy beta estimates and global FDR values.",
"",
"Figure 2. Key CP10 feature-level axes. Boxplots show representative CP10 phenotype-like features within viable CD45+ leukocytes unless otherwise specified. The main cancer-associated pattern includes increased CD13+CD66b+ granulocyte-like and CD66b+ granulocyte-like fractions, increased CD14+ monocyte-like and CD11c+HLA-DR+ APC-like myeloid phenotypes, reduced lymphoid-like fractions, and increased myeloid/granulocytic-to-lymphoid balance.",
"",
"Figure 3. Age- and sex-adjusted CP10 feature-level effect estimates. Forest plot shows cancer-versus-healthy beta estimates with 95% confidence intervals for selected CP10 features. Positive beta values indicate higher values in cancer patients, whereas negative beta values indicate lower values in cancer patients. Filled points denote features passing global FDR correction.",
"",
"Figure 4. CP10 composite-score threshold sensitivity. Cancer-versus-healthy beta estimates are shown for main, permissive, and stringent threshold configurations. The main CP10 composite axes preserved direction and global FDR significance across threshold configurations, supporting threshold robustness of the myeloid/granulocytic remodeling phenotype.",
"",
"Figure 5. Same-sex age-matched sensitivity of CP10 composite effects. Cancer-versus-healthy beta estimates are shown for the primary age- and sex-adjusted model and for same-sex nearest-age matching using 5-year and 10-year calipers. The principal CP10 composite effects remained higher in cancer after age matching.",
"",
"Figure 6. Cancer-subgroup annotation of CP10 remodeling magnitude in breast cancer. Boxplots show selected CP10 composite scores in breast cancer compared with other cancer groups in cancer-only analyses. Breast cancer showed lower integrated myeloid/granulocytic remodeling, granulocyte-like enrichment, monocyte-like HLA-DR phenotype remodeling, and myeloid/granulocytic-to-lymphoid balance scores compared with the remaining cancer groups. This analysis is secondary and should be interpreted as tumor-context-dependent variation in remodeling magnitude.",
"",
"Supplementary/Flow Figure. FlowJo-style representative CP10 gating plots. Representative healthy-control and cancer-patient samples were selected as the subjects closest to the group-specific median integrated CP10 myeloid/granulocytic remodeling score. Each panel displays real FCS events after compensation and logicle transformation. Gates show viability-low CD45+ leukocytes, CD3/CD19 lineage features, CD3/CD56 NK-like features, CD13/CD66b granulocyte-like features, CD14/HLA-DR monocyte-like features, CD11c/HLA-DR APC-like myeloid features, CCR3/CD66b eosinophil-like granulocytic features, and CD123/HLA-DR DC-like features. These gates represent phenotype-like compartments and should not be interpreted as definitive lineage assignments."
)

# ------------------------------------------------------------
# 12. Supplementary guide
# ------------------------------------------------------------

supplementary_guide <- c(
"CP10 Supplementary Output Guide",
"",
paste0("Main output folder: ", analysis_dir),
paste0("Manuscript export folder: ", out_dir),
"",
"Key analysis folders:",
"01/Step 1: FCS inventory and marker/channel QC",
"03_metadata_merge_age_QC: metadata merge and age quality control",
"04_age_sex_adjusted_models: feature-level primary and sensitivity models",
"05_composite_scores: CP10 composite score construction and models",
"06_age_sensitivity_caliper_interaction: age matching and disease-by-age interaction",
"07_threshold_sensitivity: targeted threshold sensitivity",
"08_main_figures: manuscript-style quantitative figures",
"09_flowjo_style_representative_gating: real-FCS representative gating figures",
"10_manuscript_export: final CP10 manuscript-ready tables and text exports",
"",
"Core manuscript tables exported by Step 8:",
"CP10_Table_00_cohort_QC_summary.csv",
"CP10_Table_03_primary_composite_results.csv",
"CP10_Table_04_primary_feature_results.csv",
"CP10_Table_07_age_matched_sensitivity_summary.csv",
"CP10_Table_09_threshold_robustness_summary.csv",
"CP10_Table_13_breast_cancer_one_vs_rest.csv",
"CP10_Table_19_flow_representative_gate_percentages.csv",
"",
"Interpretation safeguards:",
"1. Use phenotype-like terminology for CP10 populations.",
"2. Do not describe CD13+CD66b+ cells as definitive neutrophils.",
"3. Do not describe CCR3+CD66b+ cells as definitive eosinophils.",
"4. Do not describe CD123+HLA-DR+ cells as definitive dendritic cells.",
"5. Do not describe CD14+HLA-DR-low cells as definitive MDSCs.",
"6. Do not claim per-monocyte HLA-DR suppression; CP10 supports CD14/HLA-DR phenotype remodeling within the CD45+ compartment.",
"7. Breast cancer subgroup findings should be described as lower remodeling magnitude relative to other cancer groups, not as treatment resistance or a broad inflammatory claim."
)

# ------------------------------------------------------------
# 13. Write text files
# ------------------------------------------------------------

writeLines(methods_text, file.path(out_dir, "CP10_Methods_text.txt"), useBytes = TRUE)
writeLines(results_text, file.path(out_dir, "CP10_Results_text.txt"), useBytes = TRUE)
writeLines(figure_legends, file.path(out_dir, "CP10_Figure_legends.txt"), useBytes = TRUE)
writeLines(supplementary_guide, file.path(out_dir, "CP10_Supplementary_output_guide.txt"), useBytes = TRUE)

# Combined text
combined_text <- c(
"============================================================",
"CP10 METHODS",
"============================================================",
methods_text,
"",
"============================================================",
"CP10 RESULTS",
"============================================================",
results_text,
"",
"============================================================",
"CP10 FIGURE LEGENDS",
"============================================================",
figure_legends,
"",
"============================================================",
"CP10 SUPPLEMENTARY GUIDE",
"============================================================",
supplementary_guide
)

writeLines(combined_text, file.path(out_dir, "CP10_Manuscript_export_combined.txt"), useBytes = TRUE)

# ------------------------------------------------------------
# 14. Manifest
# ------------------------------------------------------------

export_manifest <- tibble(
  output_type = c(
    "table", "table", "table", "table", "table",
    "table", "table", "table", "text", "text",
    "text", "text", "text"
  ),
  file_name = c(
    "CP10_Table_00_cohort_QC_summary.csv",
    "CP10_Table_03_primary_composite_results.csv",
    "CP10_Table_04_primary_feature_results.csv",
    "CP10_Table_07_age_matched_sensitivity_summary.csv",
    "CP10_Table_09_threshold_robustness_summary.csv",
    "CP10_Table_13_breast_cancer_one_vs_rest.csv",
    "CP10_Table_18_flow_representative_samples.csv",
    "CP10_Table_19_flow_representative_gate_percentages.csv",
    "CP10_Methods_text.txt",
    "CP10_Results_text.txt",
    "CP10_Figure_legends.txt",
    "CP10_Supplementary_output_guide.txt",
    "CP10_Manuscript_export_combined.txt"
  ),
  file_path = file.path(out_dir, file_name)
)

write_csv(export_manifest, file.path(out_dir, "CP10_STEP8_export_manifest.csv"))

save(
  cohort_summary,
  age_summary,
  sex_summary,
  primary_composite_table,
  primary_feature_table,
  robust_feature_table,
  composite_sensitivity_table,
  matched_table,
  interaction_table,
  threshold_summary_table,
  clinical_subgroup_omnibus_table,
  clinical_breast_one_vs_rest_table,
  therapy_status_table,
  active_modality_table,
  therapy_line_table,
  time_from_start_table,
  methods_text,
  results_text,
  figure_legends,
  supplementary_guide,
  export_manifest,
  file = file.path(rdata_dir, "SDY2583_CP10_STEP8_manuscript_export.RData")
)

# ------------------------------------------------------------
# 15. Console output
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SDY2583 CP10 STEP 8 COMPLETE: MANUSCRIPT EXPORT\n")
cat("============================================================\n")

cat("\nExport folder:\n")
print(out_dir)

cat("\nCohort/QC summary:\n")
print(as.data.frame(cohort_summary), row.names = FALSE)

cat("\nPrimary composite table:\n")
print(as.data.frame(primary_composite_table), row.names = FALSE)

cat("\nAge-matched summary:\n")
print(as.data.frame(matched_table), row.names = FALSE)

cat("\nThreshold robustness counts:\n")
print(as.data.frame(threshold_counts), row.names = FALSE)

cat("\nBreast cancer one-vs-rest CP10 table:\n")
print(as.data.frame(clinical_breast_one_vs_rest_table), row.names = FALSE)

cat("\nRepresentative flow gate percentages:\n")
print(as.data.frame(gate_summary), row.names = FALSE)

cat("\nExport manifest:\n")
print(as.data.frame(export_manifest), row.names = FALSE)

cat("\nText files written:\n")
cat(" - CP10_Methods_text.txt\n")
cat(" - CP10_Results_text.txt\n")
cat(" - CP10_Figure_legends.txt\n")
cat(" - CP10_Supplementary_output_guide.txt\n")
cat(" - CP10_Manuscript_export_combined.txt\n")

cat("============================================================\n")
