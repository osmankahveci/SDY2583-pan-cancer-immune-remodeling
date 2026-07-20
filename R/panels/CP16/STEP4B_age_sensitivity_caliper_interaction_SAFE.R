# ============================================================
# SDY2583 CP16
# STEP 4B SAFE: Age sensitivity / same-sex age matching /
# disease-by-age interaction for CP16 composite scores
#
# Input:
#   outputs/CP16/11_RData/
#     SDY2583_CP16_STEP4_composite_scores.RData
#
# Analyses:
#   1) Same-sex nearest-age caliper matching, 5-year caliper
#   2) Same-sex nearest-age caliper matching, 10-year caliper
#   3) Disease-by-age interaction:
#        composite ~ disease_group * age_z + sex
#
# Output:
#   outputs/CP16/06_age_sensitivity_caliper_interaction
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

analysis_dir <- sd_analysis_dir("CP16")
rdata_dir <- file.path(analysis_dir, "11_RData")
out_dir <- file.path(analysis_dir, "06_age_sensitivity_caliper_interaction")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rdata_dir, recursive = TRUE, showWarnings = FALSE)

step4_rdata <- file.path(rdata_dir, "SDY2583_CP16_STEP4_composite_scores.RData")

if (!file.exists(step4_rdata)) {
  stop("Step 4 RData bulunamadı: ", step4_rdata)
}

load(step4_rdata)

# Reset paths after RData load.
analysis_dir <- sd_analysis_dir("CP16")
rdata_dir <- file.path(analysis_dir, "11_RData")
out_dir <- file.path(analysis_dir, "06_age_sensitivity_caliper_interaction")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!exists("cp16_scores_data")) stop("cp16_scores_data bulunamadı.")
if (!exists("composite_score_cols")) stop("composite_score_cols bulunamadı.")
if (!exists("primary_composite_results")) stop("primary_composite_results bulunamadı.")

# ------------------------------------------------------------
# 3. Matching helper
# ------------------------------------------------------------

make_same_sex_age_matched <- function(df, caliper_years = 5) {

  base_df <- df %>%
    filter(
      model_ready_age_sex == TRUE,
      !is.na(sex_binary),
      !is.na(age_for_model),
      disease_group %in% c("Healthy control", "Cancer patient")
    ) %>%
    mutate(
      disease_group = as.character(disease_group),
      sex_binary = as.character(sex_binary),
      age_for_model = as.numeric(age_for_model)
    ) %>%
    arrange(sex_binary, age_for_model, subject_id)

  healthy_pool <- base_df %>%
    filter(disease_group == "Healthy control") %>%
    select(subject_id, sex_binary, age_for_model) %>%
    mutate(used = FALSE)

  cancer_pool <- base_df %>%
    filter(disease_group == "Cancer patient") %>%
    select(subject_id, sex_binary, age_for_model) %>%
    arrange(sex_binary, age_for_model, subject_id)

  pair_rows <- list()
  pair_id <- 0L

  for (i in seq_len(nrow(cancer_pool))) {

    cc <- cancer_pool[i, ]

    candidates <- healthy_pool %>%
      filter(
        used == FALSE,
        sex_binary == cc$sex_binary,
        abs(age_for_model - cc$age_for_model) <= caliper_years
      ) %>%
      mutate(abs_age_diff = abs(age_for_model - cc$age_for_model)) %>%
      arrange(abs_age_diff, age_for_model, subject_id)

    if (nrow(candidates) == 0) next

    hh <- candidates[1, ]
    pair_id <- pair_id + 1L

    pair_rows[[pair_id]] <- tibble(
      pair_id = pair_id,
      caliper_years = caliper_years,
      cancer_subject_id = cc$subject_id,
      healthy_subject_id = hh$subject_id,
      sex_binary = cc$sex_binary,
      cancer_age = cc$age_for_model,
      healthy_age = hh$age_for_model,
      abs_age_diff = abs(cc$age_for_model - hh$age_for_model)
    )

    healthy_pool$used[healthy_pool$subject_id == hh$subject_id] <- TRUE
  }

  pairs <- bind_rows(pair_rows)

  if (nrow(pairs) == 0) {
    return(list(pairs = pairs, matched_data = tibble()))
  }

  matched_ids <- bind_rows(
    pairs %>%
      transmute(pair_id, subject_id = healthy_subject_id, matched_role = "matched_healthy"),
    pairs %>%
      transmute(pair_id, subject_id = cancer_subject_id, matched_role = "matched_cancer")
  )

  matched_data <- df %>%
    inner_join(matched_ids, by = "subject_id") %>%
    mutate(
      disease_group = factor(as.character(disease_group), levels = c("Healthy control", "Cancer patient")),
      sex = factor(as.character(sex))
    )

  list(pairs = pairs, matched_data = matched_data)
}

# ------------------------------------------------------------
# 4. Composite model helper
# ------------------------------------------------------------

run_composite_model <- function(df, score_name, model_label, caliper_years = NA_real_) {

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

  if (nrow(model_df) < 50 ||
      sum(model_df$disease_group == "Healthy control") < 20 ||
      sum(model_df$disease_group == "Cancer patient") < 20) {
    return(NULL)
  }

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
    n_model = nrow(model_df),
    n_healthy = sum(model_df$disease_group == "Healthy control"),
    n_cancer = sum(model_df$disease_group == "Cancer patient"),
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
    ),
    caliper_years = caliper_years
  )
}

run_matched_set <- function(matched_df, model_label, caliper_years) {
  bind_rows(lapply(composite_score_cols, function(sc) {
    run_composite_model(matched_df, sc, model_label, caliper_years)
  })) %>%
    mutate(FDR_global = p.adjust(p_value, method = "BH")) %>%
    arrange(FDR_global, p_value)
}

# ------------------------------------------------------------
# 5. Create matched datasets
# ------------------------------------------------------------

matched_5y <- make_same_sex_age_matched(cp16_scores_data, caliper_years = 5)
matched_10y <- make_same_sex_age_matched(cp16_scores_data, caliper_years = 10)

pairs_5y <- matched_5y$pairs
pairs_10y <- matched_10y$pairs

matched_5y_data <- matched_5y$matched_data
matched_10y_data <- matched_10y$matched_data

matched_dataset_summary <- bind_rows(
  tibble(
    caliper_years = 5,
    n_rows = nrow(matched_5y_data),
    n_pairs = nrow(pairs_5y),
    mean_abs_age_diff = mean(pairs_5y$abs_age_diff, na.rm = TRUE),
    median_abs_age_diff = median(pairs_5y$abs_age_diff, na.rm = TRUE),
    max_abs_age_diff = max(pairs_5y$abs_age_diff, na.rm = TRUE),
    n_healthy = sum(matched_5y_data$disease_group == "Healthy control"),
    n_cancer = sum(matched_5y_data$disease_group == "Cancer patient")
  ),
  tibble(
    caliper_years = 10,
    n_rows = nrow(matched_10y_data),
    n_pairs = nrow(pairs_10y),
    mean_abs_age_diff = mean(pairs_10y$abs_age_diff, na.rm = TRUE),
    median_abs_age_diff = median(pairs_10y$abs_age_diff, na.rm = TRUE),
    max_abs_age_diff = max(pairs_10y$abs_age_diff, na.rm = TRUE),
    n_healthy = sum(matched_10y_data$disease_group == "Healthy control"),
    n_cancer = sum(matched_10y_data$disease_group == "Cancer patient")
  )
)

# ------------------------------------------------------------
# 6. Run matched composite models
# ------------------------------------------------------------

matched_5y_results <- run_matched_set(
  matched_5y_data,
  "same_sex_age_matched_5y",
  5
)

matched_10y_results <- run_matched_set(
  matched_10y_data,
  "same_sex_age_matched_10y",
  10
)

# ------------------------------------------------------------
# 7. Disease-by-age interaction
# ------------------------------------------------------------

run_interaction_model <- function(df, score_name) {

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
      sex = droplevels(factor(as.character(sex))),
      age_z = as.numeric(scale(age_for_model))
    )

  if (nrow(model_df) < 50 ||
      sum(model_df$disease_group == "Healthy control") < 20 ||
      sum(model_df$disease_group == "Cancer patient") < 20) {
    return(NULL)
  }

  fit <- tryCatch(
    lm(value ~ disease_group * age_z + sex, data = model_df),
    error = function(e) NULL
  )

  if (is.null(fit)) return(NULL)

  tt <- tryCatch(broom::tidy(fit, conf.int = TRUE), error = function(e) NULL)
  if (is.null(tt)) return(NULL)

  term <- "disease_groupCancer patient:age_z"
  if (!(term %in% tt$term)) return(NULL)

  out <- tt %>% filter(term == !!term)

  tibble(
    composite_score = score_name,
    n_model = nrow(model_df),
    n_healthy = sum(model_df$disease_group == "Healthy control"),
    n_cancer = sum(model_df$disease_group == "Cancer patient"),
    beta_interaction_cancer_by_age_z = out$estimate[1],
    conf_low = out$conf.low[1],
    conf_high = out$conf.high[1],
    p_value = out$p.value[1],
    direction = case_when(
      out$estimate[1] > 0 ~ "stronger_cancer_effect_with_older_age",
      out$estimate[1] < 0 ~ "weaker_cancer_effect_with_older_age",
      TRUE ~ "no_interaction_direction"
    )
  )
}

interaction_results <- bind_rows(lapply(composite_score_cols, function(sc) {
  run_interaction_model(
    cp16_scores_data %>% filter(model_ready_age_sex == TRUE),
    sc
  )
})) %>%
  mutate(FDR_global = p.adjust(p_value, method = "BH")) %>%
  arrange(FDR_global, p_value)

# ------------------------------------------------------------
# 8. Matched direction summary
# ------------------------------------------------------------

get_matched_value <- function(df, score_name, col_name) {
  x <- df %>%
    filter(composite_score == score_name) %>%
    pull(!!rlang::sym(col_name))
  if (length(x) == 0) return(NA)
  x[1]
}

matched_direction_summary <- bind_rows(lapply(composite_score_cols, function(sc) {

  primary_beta <- primary_composite_results %>%
    filter(composite_score == sc) %>%
    pull(beta_cancer_vs_healthy)

  primary_direction <- primary_composite_results %>%
    filter(composite_score == sc) %>%
    pull(direction)

  primary_FDR <- primary_composite_results %>%
    filter(composite_score == sc) %>%
    pull(FDR_global)

  if (length(primary_beta) == 0) primary_beta <- NA_real_
  if (length(primary_direction) == 0) primary_direction <- NA_character_
  if (length(primary_FDR) == 0) primary_FDR <- NA_real_

  matched_5y_beta <- as.numeric(get_matched_value(matched_5y_results, sc, "beta_cancer_vs_healthy"))
  matched_5y_direction <- as.character(get_matched_value(matched_5y_results, sc, "direction"))
  matched_5y_FDR <- as.numeric(get_matched_value(matched_5y_results, sc, "FDR_global"))

  matched_10y_beta <- as.numeric(get_matched_value(matched_10y_results, sc, "beta_cancer_vs_healthy"))
  matched_10y_direction <- as.character(get_matched_value(matched_10y_results, sc, "direction"))
  matched_10y_FDR <- as.numeric(get_matched_value(matched_10y_results, sc, "FDR_global"))

  tibble(
    composite_score = sc,
    primary_beta = primary_beta[1],
    primary_direction = primary_direction[1],
    primary_FDR_global = primary_FDR[1],
    matched_5y_beta = matched_5y_beta,
    matched_5y_direction = matched_5y_direction,
    matched_5y_FDR_global = matched_5y_FDR,
    matched_10y_beta = matched_10y_beta,
    matched_10y_direction = matched_10y_direction,
    matched_10y_FDR_global = matched_10y_FDR,
    direction_preserved_5y = !is.na(primary_direction[1]) && primary_direction[1] == matched_5y_direction,
    direction_preserved_10y = !is.na(primary_direction[1]) && primary_direction[1] == matched_10y_direction,
    direction_preserved_both_matched =
      !is.na(primary_direction[1]) &&
      primary_direction[1] == matched_5y_direction &&
      primary_direction[1] == matched_10y_direction,
    FDR_global_preserved_5y =
      !is.na(primary_FDR[1]) &&
      primary_FDR[1] < 0.05 &&
      !is.na(matched_5y_FDR) &&
      matched_5y_FDR < 0.05,
    FDR_global_preserved_10y =
      !is.na(primary_FDR[1]) &&
      primary_FDR[1] < 0.05 &&
      !is.na(matched_10y_FDR) &&
      matched_10y_FDR < 0.05,
    FDR_global_preserved_both_matched =
      !is.na(primary_FDR[1]) &&
      primary_FDR[1] < 0.05 &&
      !is.na(matched_5y_FDR) &&
      matched_5y_FDR < 0.05 &&
      !is.na(matched_10y_FDR) &&
      matched_10y_FDR < 0.05
  )
})) %>%
  arrange(primary_FDR_global)

overall_age_sensitivity_summary <- tibble(
  n_primary_model_subjects = primary_composite_results$n_model[1],
  n_primary_healthy = primary_composite_results$n_healthy[1],
  n_primary_cancer = primary_composite_results$n_cancer[1],
  n_5y_pairs = nrow(pairs_5y),
  n_10y_pairs = nrow(pairs_10y),
  n_scores_tested = length(composite_score_cols),
  n_5y_FDR_lt_0p05 = sum(matched_5y_results$FDR_global < 0.05, na.rm = TRUE),
  n_10y_FDR_lt_0p05 = sum(matched_10y_results$FDR_global < 0.05, na.rm = TRUE),
  n_direction_preserved_both_matched = sum(matched_direction_summary$direction_preserved_both_matched == TRUE, na.rm = TRUE),
  n_FDR_preserved_both_matched = sum(matched_direction_summary$FDR_global_preserved_both_matched == TRUE, na.rm = TRUE),
  n_interaction_FDR_lt_0p05 = sum(interaction_results$FDR_global < 0.05, na.rm = TRUE)
)

# ------------------------------------------------------------
# 9. Save outputs
# ------------------------------------------------------------

write_csv(pairs_5y, file.path(out_dir, "SDY2583_CP16_same_sex_age_matched_pairs_5y_STEP4B.csv"))
write_csv(pairs_10y, file.path(out_dir, "SDY2583_CP16_same_sex_age_matched_pairs_10y_STEP4B.csv"))
write_csv(matched_dataset_summary, file.path(out_dir, "SDY2583_CP16_matched_dataset_summary_STEP4B.csv"))
write_csv(matched_5y_results, file.path(out_dir, "SDY2583_CP16_same_sex_age_matched_5y_results_STEP4B.csv"))
write_csv(matched_10y_results, file.path(out_dir, "SDY2583_CP16_same_sex_age_matched_10y_results_STEP4B.csv"))
write_csv(interaction_results, file.path(out_dir, "SDY2583_CP16_disease_by_age_interaction_results_STEP4B.csv"))
write_csv(matched_direction_summary, file.path(out_dir, "SDY2583_CP16_matched_direction_sensitivity_summary_STEP4B.csv"))
write_csv(overall_age_sensitivity_summary, file.path(out_dir, "SDY2583_CP16_overall_age_sensitivity_summary_STEP4B.csv"))

save(
  pairs_5y,
  pairs_10y,
  matched_5y_data,
  matched_10y_data,
  matched_dataset_summary,
  matched_5y_results,
  matched_10y_results,
  interaction_results,
  matched_direction_summary,
  overall_age_sensitivity_summary,
  file = file.path(rdata_dir, "SDY2583_CP16_STEP4B_age_sensitivity_caliper_interaction.RData")
)

# ------------------------------------------------------------
# 10. Console output
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SDY2583 CP16 STEP 4B COMPLETE: AGE SENSITIVITY / MATCHING / INTERACTION\n")
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
print(as.data.frame(matched_direction_summary), row.names = FALSE)

cat("\nFiles saved in:\n")
print(out_dir)

cat("============================================================\n")
