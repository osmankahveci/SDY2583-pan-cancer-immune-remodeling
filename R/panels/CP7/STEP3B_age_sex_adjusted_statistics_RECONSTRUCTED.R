# ============================================================
# SDY2583 CP7
# STEP 3B RECONSTRUCTED: age/sex-adjusted statistics
#
# Reconstruction basis:
#   - archived CP7 feature-module map;
#   - archived model schema and output columns;
#   - archived main, binary-sex, and event-count sensitivity tables.
#
# This is reconstructed source, not the original archived script.
# Validate all coefficients, confidence intervals, and FDR values against the
# archived CP7 statistics before treating this step as verified.
# ============================================================

rm(list = ls())
source(file.path(Sys.getenv("SDY2583_REPO_ROOT", unset = "."), "R", "shared", "bootstrap.R"))

cran_pkgs <- c("dplyr", "readr", "tibble", "purrr")
for (p in cran_pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(purrr)
})

analysis_dir <- sd_analysis_dir("CP7")
rdata_dir <- file.path(analysis_dir, "11_RData")
out_dir <- file.path(analysis_dir, "06_statistics")
sensitivity_dir <- file.path(analysis_dir, "07_sensitivity")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(sensitivity_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rdata_dir, recursive = TRUE, showWarnings = FALSE)

step3a_file <- file.path(rdata_dir, "SDY2583_CP7_STEP3A_metadata_merge_age_QC_RECONSTRUCTED.RData")
if (!file.exists(step3a_file)) stop("Run reconstructed CP7 Step 3A first: ", step3a_file)
load(step3a_file)
if (!exists("analysis_data")) stop("Step 3A RData does not contain analysis_data.")

composition <- c(
  "pct_cd3_pos_total", "pct_cd8_pos_total", "pct_cd3_cd8_pos_total", "pct_cd8_within_cd3"
)

differentiation <- c(
  "median_CD62L_in_CD3CD8", "median_CD27_in_CD3CD8", "median_CD45RA_in_CD3CD8",
  "pct_naive_like", "pct_tcm_like", "pct_tem_like", "pct_temra_like",
  "pct_temra_cd27neg", "pct_temra_cd27pos", "pct_cd62lneg_cd27neg",
  "pct_cd45ra_pos_cd62lneg_cd27neg", "pct_naive_cd27pos", "pct_naive_cd27neg"
)

checkpoint_medians <- c(
  "median_PD1_in_CD3CD8", "median_TIM3_in_CD3CD8", "median_LAG3_in_CD3CD8",
  "median_TIGIT_in_CD3CD8", "median_ICOS_in_CD3CD8", "median_CD39_in_CD3CD8"
)

checkpoint_single_positive <- c(
  "pct_pd1_pos", "pct_pd1_high", "pct_tim3_pos", "pct_lag3_pos",
  "pct_tigit_pos", "pct_icos_pos", "pct_cd39_pos"
)

checkpoint_coexpression <- c(
  "pct_pd1_tim3_pos", "pct_pd1_lag3_pos", "pct_pd1_tigit_pos",
  "pct_pd1_cd39_pos", "pct_tim3_lag3_pos", "pct_tigit_cd39_pos",
  "pct_icos_cd39_pos", "pct_pd1_tigit_cd39_pos", "pct_pd1_tim3_lag3_pos",
  "pct_pd1_tigit_icos_cd39_pos", "pct_pd1_tigit_cd39_temra_cd27neg"
)

terminal_checkpoint <- c(
  "pct_pd1_pos_within_temra_like", "pct_tim3_pos_within_temra_like",
  "pct_lag3_pos_within_temra_like", "pct_tigit_pos_within_temra_like",
  "pct_cd39_pos_within_temra_like", "pct_pd1_pos_within_temra_cd27neg",
  "pct_tigit_pos_within_temra_cd27neg", "pct_cd39_pos_within_temra_cd27neg",
  "pct_pd1_temra", "pct_tigit_temra", "pct_cd39_temra",
  "pct_pd1_temra_cd27neg", "pct_tigit_temra_cd27neg", "pct_cd39_temra_cd27neg"
)

feature_module_map <- bind_rows(
  tibble(feature = composition, module = "composition"),
  tibble(feature = differentiation, module = "differentiation"),
  tibble(feature = checkpoint_medians, module = "checkpoint_medians"),
  tibble(feature = checkpoint_single_positive, module = "checkpoint_single_positive"),
  tibble(feature = checkpoint_coexpression, module = "checkpoint_coexpression"),
  tibble(feature = terminal_checkpoint, module = "terminal_checkpoint")
)

missing_features <- setdiff(feature_module_map$feature, names(analysis_data))
if (length(missing_features) > 0L) {
  stop("Required CP7 features are missing: ", paste(missing_features, collapse = ", "))
}

prepare_model_data <- function(data, feature, binary_sex_only = FALSE, minimum_events = 0L) {
  out <- data %>%
    transmute(
      y = suppressWarnings(as.numeric(.data[[feature]])),
      disease_group = factor(as.character(disease_group), levels = c("Healthy control", "Cancer patient")),
      age_for_model = suppressWarnings(as.numeric(age_for_model)),
      sex = factor(as.character(sex)),
      n_cd3_cd8_pos = suppressWarnings(as.numeric(n_cd3_cd8_pos))
    ) %>%
    filter(
      !is.na(y), !is.na(disease_group), !is.na(age_for_model), !is.na(sex),
      !is.na(n_cd3_cd8_pos), n_cd3_cd8_pos >= minimum_events
    )

  if (binary_sex_only) {
    out <- out %>%
      filter(as.character(sex) %in% c("Female", "Male")) %>%
      mutate(sex = factor(as.character(sex), levels = c("Female", "Male")))
  }
  out
}

fit_one_feature <- function(data, feature, module, subset_label,
                            binary_sex_only = FALSE, minimum_events = 0L) {
  dat <- prepare_model_data(data, feature, binary_sex_only, minimum_events)
  if (nrow(dat) < 10L || n_distinct(dat$disease_group) < 2L) {
    return(tibble(
      subset_label = subset_label, module = module, feature = feature,
      n_model = nrow(dat), n_healthy = sum(dat$disease_group == "Healthy control"),
      n_cancer = sum(dat$disease_group == "Cancer patient"),
      healthy_mean = NA_real_, cancer_mean = NA_real_, healthy_median = NA_real_,
      cancer_median = NA_real_, crude_mean_difference_cancer_minus_healthy = NA_real_,
      beta_cancer_vs_healthy = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
      t_value = NA_real_, p_value = NA_real_,
      model_formula = "y ~ disease_group + age_for_model + sex", direction = NA_character_
    ))
  }

  fit <- stats::lm(y ~ disease_group + age_for_model + sex, data = dat)
  coef_table <- summary(fit)$coefficients
  coef_name <- grep("^disease_group", rownames(coef_table), value = TRUE)[1]
  if (is.na(coef_name)) stop("Cancer coefficient was not found for feature: ", feature)
  ci <- stats::confint(fit, parm = coef_name, level = 0.95)

  healthy_values <- dat$y[dat$disease_group == "Healthy control"]
  cancer_values <- dat$y[dat$disease_group == "Cancer patient"]
  healthy_mean <- mean(healthy_values, na.rm = TRUE)
  cancer_mean <- mean(cancer_values, na.rm = TRUE)
  beta <- unname(coef_table[coef_name, "Estimate"])

  tibble(
    subset_label = subset_label,
    module = module,
    feature = feature,
    n_model = stats::nobs(fit),
    n_healthy = length(healthy_values),
    n_cancer = length(cancer_values),
    healthy_mean = healthy_mean,
    cancer_mean = cancer_mean,
    healthy_median = median(healthy_values, na.rm = TRUE),
    cancer_median = median(cancer_values, na.rm = TRUE),
    crude_mean_difference_cancer_minus_healthy = cancer_mean - healthy_mean,
    beta_cancer_vs_healthy = beta,
    ci_low = unname(ci[1]),
    ci_high = unname(ci[2]),
    t_value = unname(coef_table[coef_name, "t value"]),
    p_value = unname(coef_table[coef_name, "Pr(>|t|)"]),
    model_formula = "y ~ disease_group + age_for_model + sex",
    direction = ifelse(beta > 0, "higher_in_cancer", ifelse(beta < 0, "lower_in_cancer", "no_difference"))
  )
}

fit_feature_set <- function(data, subset_label, binary_sex_only = FALSE, minimum_events = 0L) {
  results <- purrr::map2_dfr(
    feature_module_map$feature,
    feature_module_map$module,
    ~fit_one_feature(
      data = data, feature = .x, module = .y, subset_label = subset_label,
      binary_sex_only = binary_sex_only, minimum_events = minimum_events
    )
  )

  results %>%
    mutate(fdr_all = p.adjust(p_value, method = "BH")) %>%
    group_by(module) %>%
    mutate(fdr_within_module = p.adjust(p_value, method = "BH")) %>%
    ungroup() %>%
    arrange(p_value, module, feature)
}

main_statistics <- fit_feature_set(
  analysis_data,
  subset_label = "main_valid_age_all_sex_categories",
  binary_sex_only = FALSE,
  minimum_events = 0L
)

binary_sex_statistics <- fit_feature_set(
  analysis_data,
  subset_label = "binary_sex_only",
  binary_sex_only = TRUE,
  minimum_events = 0L
)

event_thresholds <- c(0L, 500L, 1000L, 2000L)
event_threshold_statistics <- purrr::map_dfr(event_thresholds, function(threshold) {
  fit_feature_set(
    analysis_data,
    subset_label = paste0("min_CD3CD8_events_", threshold),
    binary_sex_only = FALSE,
    minimum_events = threshold
  )
})

significant_global <- main_statistics %>% filter(fdr_all < 0.05)
significant_module <- main_statistics %>% filter(fdr_within_module < 0.05)
top10_per_module <- main_statistics %>% group_by(module) %>% slice_min(p_value, n = 10, with_ties = FALSE) %>% ungroup()
compact_significant <- main_statistics %>%
  filter(fdr_all < 0.05 | fdr_within_module < 0.05) %>%
  select(module, feature, n_model, healthy_mean, cancer_mean, beta_cancer_vs_healthy,
         ci_low, ci_high, p_value, fdr_all, fdr_within_module, direction)

readr::write_csv(feature_module_map, file.path(out_dir, "SDY2583_CP7_feature_module_map_RECONSTRUCTED.csv"))
readr::write_csv(main_statistics, file.path(out_dir, "SDY2583_CP7_FULL850_age_sex_adjusted_statistics_MAIN_RECONSTRUCTED.csv"))
readr::write_csv(significant_global, file.path(out_dir, "SDY2583_CP7_FULL850_significant_global_FDR_0p05_MAIN_RECONSTRUCTED.csv"))
readr::write_csv(significant_module, file.path(out_dir, "SDY2583_CP7_FULL850_significant_module_FDR_0p05_MAIN_RECONSTRUCTED.csv"))
readr::write_csv(top10_per_module, file.path(out_dir, "SDY2583_CP7_top10_results_per_module_MAIN_RECONSTRUCTED.csv"))
readr::write_csv(compact_significant, file.path(out_dir, "SDY2583_CP7_compact_significant_results_MAIN_RECONSTRUCTED.csv"))
readr::write_csv(binary_sex_statistics, file.path(sensitivity_dir, "SDY2583_CP7_binary_sex_only_sensitivity_statistics_RECONSTRUCTED.csv"))
readr::write_csv(event_threshold_statistics, file.path(sensitivity_dir, "SDY2583_CP7_event_threshold_sensitivity_statistics_RECONSTRUCTED.csv"))

save(
  feature_module_map, main_statistics, binary_sex_statistics,
  event_threshold_statistics, significant_global, significant_module,
  file = file.path(rdata_dir, "SDY2583_CP7_STEP3B_age_sex_adjusted_statistics_RECONSTRUCTED.RData")
)

cat("CP7 reconstructed main models:", nrow(main_statistics), "features\n")
cat("Expected archived main-model N: 832 (398 healthy, 434 cancer).\n")
