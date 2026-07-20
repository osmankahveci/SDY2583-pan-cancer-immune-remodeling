# ============================================================
# SDY2583 CP23
# STEP 4B SAFE: Age sensitivity / same-sex age-caliper matching
#               / disease-by-age interaction
#
# Input:
#   outputs/CP23/11_RData/
#     SDY2583_CP23_STEP4_composite_scores.RData
#
# Output:
#   outputs/CP23/
#     06_age_sensitivity_caliper_interaction
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
out_dir <- file.path(analysis_dir, "06_age_sensitivity_caliper_interaction")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

step4_rdata <- file.path(rdata_dir, "SDY2583_CP23_STEP4_composite_scores.RData")

if (!file.exists(step4_rdata)) {
  stop("Step 4 RData bulunamadı: ", step4_rdata)
}

load(step4_rdata)

# Reset paths after load.
analysis_dir <- sd_analysis_dir("CP23")
rdata_dir <- file.path(analysis_dir, "11_RData")
out_dir <- file.path(analysis_dir, "06_age_sensitivity_caliper_interaction")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!exists("cp23_score_data")) stop("cp23_score_data bulunamadı.")
if (!exists("score_names")) stop("score_names bulunamadı.")

if (!exists("composite_score_dictionary")) {
  composite_score_dictionary <- tibble(
    composite_score = score_names,
    score_label = score_names
  )
}

# ------------------------------------------------------------
# Model helper
# ------------------------------------------------------------

run_score_model <- function(df, score_name, model_label, sex_var = "sex", caliper_years = NA_real_) {

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
    ),
    caliper_years = caliper_years
  )
}

add_score_fdr <- function(res_df) {
  if (is.null(res_df) || nrow(res_df) == 0) return(tibble())

  res_df %>%
    left_join(composite_score_dictionary, by = "composite_score") %>%
    mutate(FDR_global = p.adjust(p_value, method = "BH")) %>%
    arrange(FDR_global, p_value)
}

# ------------------------------------------------------------
# Same-sex nearest-age caliper matching without replacement
# ------------------------------------------------------------

make_same_sex_age_matched <- function(df, caliper_years = 5) {

  base <- df %>%
    filter(
      model_ready_age_sex == TRUE,
      disease_group %in% c("Healthy control", "Cancer patient"),
      !is.na(age_for_model),
      !is.na(sex)
    ) %>%
    mutate(
      sex_match = as.character(sex)
    )

  # Prefer biologically interpretable same-sex matching.
  # If sex_binary exists, use only Female/Male-coded rows.
  if ("sex_binary" %in% names(base)) {
    base <- base %>%
      filter(!is.na(sex_binary)) %>%
      mutate(sex_match = as.character(sex_binary))
  } else {
    base <- base %>%
      filter(sex_match %in% c("Female", "Male"))
  }

  healthy <- base %>%
    filter(disease_group == "Healthy control") %>%
    arrange(sex_match, age_for_model)

  cancer <- base %>%
    filter(disease_group == "Cancer patient") %>%
    arrange(sex_match, age_for_model)

  used_healthy <- character(0)
  pair_list <- list()
  pair_id <- 0L

  # Match cancer subjects to nearest available healthy subject of same sex.
  for (i in seq_len(nrow(cancer))) {
    c_row <- cancer[i, , drop = FALSE]

    candidates <- healthy %>%
      filter(
        sex_match == c_row$sex_match[1],
        !(subject_id %in% used_healthy)
      ) %>%
      mutate(abs_age_diff = abs(age_for_model - c_row$age_for_model[1])) %>%
      filter(abs_age_diff <= caliper_years) %>%
      arrange(abs_age_diff, age_for_model)

    if (nrow(candidates) == 0) next

    h_row <- candidates[1, , drop = FALSE]
    used_healthy <- c(used_healthy, h_row$subject_id[1])
    pair_id <- pair_id + 1L

    h_row$matched_pair_id <- paste0("pair_", caliper_years, "y_", pair_id)
    c_row$matched_pair_id <- paste0("pair_", caliper_years, "y_", pair_id)
    h_row$abs_age_diff_pair <- abs(h_row$age_for_model[1] - c_row$age_for_model[1])
    c_row$abs_age_diff_pair <- abs(h_row$age_for_model[1] - c_row$age_for_model[1])
    h_row$caliper_years <- caliper_years
    c_row$caliper_years <- caliper_years

    pair_list[[length(pair_list) + 1L]] <- h_row
    pair_list[[length(pair_list) + 1L]] <- c_row
  }

  if (length(pair_list) == 0) return(tibble())

  bind_rows(pair_list) %>%
    arrange(matched_pair_id, disease_group)
}

matched_5y <- make_same_sex_age_matched(cp23_score_data, caliper_years = 5)
matched_10y <- make_same_sex_age_matched(cp23_score_data, caliper_years = 10)

matched_dataset_summary <- bind_rows(
  matched_5y %>%
    summarise(
      caliper_years = 5,
      n_rows = n(),
      n_pairs = n_distinct(matched_pair_id),
      mean_abs_age_diff = mean(abs_age_diff_pair, na.rm = TRUE),
      median_abs_age_diff = median(abs_age_diff_pair, na.rm = TRUE),
      max_abs_age_diff = max(abs_age_diff_pair, na.rm = TRUE),
      n_healthy = sum(disease_group == "Healthy control"),
      n_cancer = sum(disease_group == "Cancer patient")
    ),
  matched_10y %>%
    summarise(
      caliper_years = 10,
      n_rows = n(),
      n_pairs = n_distinct(matched_pair_id),
      mean_abs_age_diff = mean(abs_age_diff_pair, na.rm = TRUE),
      median_abs_age_diff = median(abs_age_diff_pair, na.rm = TRUE),
      max_abs_age_diff = max(abs_age_diff_pair, na.rm = TRUE),
      n_healthy = sum(disease_group == "Healthy control"),
      n_cancer = sum(disease_group == "Cancer patient")
    )
)

matched_5y_results <- bind_rows(lapply(score_names, function(ss) {
  run_score_model(matched_5y, ss, "same_sex_age_matched_5y", sex_var = "sex_match", caliper_years = 5)
})) %>% add_score_fdr()

matched_10y_results <- bind_rows(lapply(score_names, function(ss) {
  run_score_model(matched_10y, ss, "same_sex_age_matched_10y", sex_var = "sex_match", caliper_years = 10)
})) %>% add_score_fdr()

# ------------------------------------------------------------
# Disease-by-age interaction
# ------------------------------------------------------------

run_interaction_model <- function(df, score_name) {

  if (!(score_name %in% names(df))) return(NULL)

  work <- df %>%
    filter(model_ready_age_sex == TRUE) %>%
    transmute(
      value = suppressWarnings(as.numeric(.data[[score_name]])),
      disease_group = disease_group,
      age_for_model = age_for_model,
      sex_model = sex
    ) %>%
    filter(
      !is.na(value),
      !is.na(disease_group),
      !is.na(age_for_model),
      !is.na(sex_model),
      disease_group %in% c("Healthy control", "Cancer patient")
    ) %>%
    mutate(
      disease_group = factor(
        as.character(disease_group),
        levels = c("Healthy control", "Cancer patient")
      ),
      sex_model = factor(sex_model),
      age_z = as.numeric(scale(age_for_model))
    )

  if (nrow(work) < 30) return(NULL)
  if (length(unique(work$disease_group)) < 2) return(NULL)

  fit <- tryCatch(
    lm(value ~ disease_group * age_z + sex_model, data = work),
    error = function(e) NULL
  )

  if (is.null(fit)) return(NULL)

  tt <- tryCatch(broom::tidy(fit, conf.int = TRUE), error = function(e) NULL)
  if (is.null(tt)) return(NULL)

  interaction_term <- "disease_groupCancer patient:age_z"
  if (!(interaction_term %in% tt$term)) return(NULL)

  term_row <- tt %>% filter(term == interaction_term)

  tibble(
    composite_score = score_name,
    n_model = nrow(work),
    n_healthy = sum(work$disease_group == "Healthy control"),
    n_cancer = sum(work$disease_group == "Cancer patient"),
    beta_interaction_cancer_by_age_z = term_row$estimate[1],
    conf_low = term_row$conf.low[1],
    conf_high = term_row$conf.high[1],
    p_value = term_row$p.value[1],
    direction = case_when(
      term_row$estimate[1] > 0 ~ "stronger_cancer_effect_with_older_age",
      term_row$estimate[1] < 0 ~ "weaker_cancer_effect_with_older_age",
      TRUE ~ "no_interaction_direction"
    )
  )
}

interaction_results <- bind_rows(lapply(score_names, function(ss) {
  run_interaction_model(cp23_score_data, ss)
})) %>%
  left_join(composite_score_dictionary, by = "composite_score") %>%
  mutate(FDR_global = p.adjust(p_value, method = "BH")) %>%
  arrange(FDR_global, p_value)

# ------------------------------------------------------------
# Preservation summary
# ------------------------------------------------------------

primary_join <- primary_composite_results %>%
  select(
    composite_score,
    score_label,
    primary_beta = beta_cancer_vs_healthy,
    primary_direction = direction,
    primary_FDR_global = FDR_global
  )

matched_5y_join <- matched_5y_results %>%
  select(
    composite_score,
    matched_5y_beta = beta_cancer_vs_healthy,
    matched_5y_direction = direction,
    matched_5y_FDR_global = FDR_global
  )

matched_10y_join <- matched_10y_results %>%
  select(
    composite_score,
    matched_10y_beta = beta_cancer_vs_healthy,
    matched_10y_direction = direction,
    matched_10y_FDR_global = FDR_global
  )

matched_direction_sensitivity_summary <- primary_join %>%
  left_join(matched_5y_join, by = "composite_score") %>%
  left_join(matched_10y_join, by = "composite_score") %>%
  mutate(
    direction_preserved_5y = primary_direction == matched_5y_direction,
    direction_preserved_10y = primary_direction == matched_10y_direction,
    direction_preserved_both_matched =
      direction_preserved_5y == TRUE & direction_preserved_10y == TRUE,
    FDR_global_preserved_5y =
      primary_FDR_global < 0.05 & matched_5y_FDR_global < 0.05,
    FDR_global_preserved_10y =
      primary_FDR_global < 0.05 & matched_10y_FDR_global < 0.05,
    FDR_global_preserved_both_matched =
      FDR_global_preserved_5y == TRUE & FDR_global_preserved_10y == TRUE
  ) %>%
  arrange(primary_FDR_global)

overall_age_sensitivity_summary <- tibble(
  n_primary_model_subjects = nrow(cp23_score_data %>% filter(model_ready_age_sex == TRUE)),
  n_primary_healthy = sum(cp23_score_data$model_ready_age_sex == TRUE & cp23_score_data$disease_group == "Healthy control", na.rm = TRUE),
  n_primary_cancer = sum(cp23_score_data$model_ready_age_sex == TRUE & cp23_score_data$disease_group == "Cancer patient", na.rm = TRUE),
  n_5y_pairs = ifelse(nrow(matched_5y) > 0, n_distinct(matched_5y$matched_pair_id), 0),
  n_10y_pairs = ifelse(nrow(matched_10y) > 0, n_distinct(matched_10y$matched_pair_id), 0),
  n_scores_tested = length(score_names),
  n_5y_FDR_lt_0p05 = sum(matched_5y_results$FDR_global < 0.05, na.rm = TRUE),
  n_10y_FDR_lt_0p05 = sum(matched_10y_results$FDR_global < 0.05, na.rm = TRUE),
  n_direction_preserved_both_matched =
    sum(matched_direction_sensitivity_summary$direction_preserved_both_matched == TRUE, na.rm = TRUE),
  n_FDR_preserved_both_matched =
    sum(matched_direction_sensitivity_summary$FDR_global_preserved_both_matched == TRUE, na.rm = TRUE),
  n_interaction_FDR_lt_0p05 =
    sum(interaction_results$FDR_global < 0.05, na.rm = TRUE)
)

# ------------------------------------------------------------
# Save outputs
# ------------------------------------------------------------

write_csv(matched_5y, file.path(out_dir, "SDY2583_CP23_same_sex_age_matched_5y_dataset_STEP4B.csv"))
write_csv(matched_10y, file.path(out_dir, "SDY2583_CP23_same_sex_age_matched_10y_dataset_STEP4B.csv"))
write_csv(matched_dataset_summary, file.path(out_dir, "SDY2583_CP23_matched_dataset_summary_STEP4B.csv"))
write_csv(matched_5y_results, file.path(out_dir, "SDY2583_CP23_same_sex_age_matched_5y_results_STEP4B.csv"))
write_csv(matched_10y_results, file.path(out_dir, "SDY2583_CP23_same_sex_age_matched_10y_results_STEP4B.csv"))
write_csv(interaction_results, file.path(out_dir, "SDY2583_CP23_disease_by_age_interaction_results_STEP4B.csv"))
write_csv(matched_direction_sensitivity_summary, file.path(out_dir, "SDY2583_CP23_matched_direction_sensitivity_summary_STEP4B.csv"))
write_csv(overall_age_sensitivity_summary, file.path(out_dir, "SDY2583_CP23_overall_age_sensitivity_summary_STEP4B.csv"))

save(
  matched_5y,
  matched_10y,
  matched_dataset_summary,
  matched_5y_results,
  matched_10y_results,
  interaction_results,
  matched_direction_sensitivity_summary,
  overall_age_sensitivity_summary,
  file = file.path(rdata_dir, "SDY2583_CP23_STEP4B_age_sensitivity_caliper_interaction.RData")
)

cat("\n============================================================\n")
cat("SDY2583 CP23 STEP 4B COMPLETE: AGE SENSITIVITY / MATCHING / INTERACTION\n")
cat("============================================================\n")

cat("\nMatched dataset summary:\n")
print(as.data.frame(matched_dataset_summary), row.names = FALSE)

cat("\nOverall age sensitivity summary:\n")
print(as.data.frame(overall_age_sensitivity_summary), row.names = FALSE)

cat("\nSame-sex age-matched 5-year results:\n")
print(as.data.frame(matched_5y_results), row.names = FALSE)

cat("\nSame-sex age-matched 10-year results:\n")
print(as.data.frame(matched_10y_results), row.names = FALSE)

cat("\nDisease-by-age interaction results:\n")
print(as.data.frame(interaction_results), row.names = FALSE)

cat("\nMatched direction sensitivity summary:\n")
print(as.data.frame(matched_direction_sensitivity_summary), row.names = FALSE)

cat("\nFiles saved in:\n")
print(out_dir)

cat("============================================================\n")
