# ============================================================
# SDY2583 CP22
# STEP 4B SAFE: Age sensitivity, same-sex age-caliper matching,
# and disease-by-age interaction for CP22 composite scores
#
# Uses:
# outputs/CP22/11_RData/
#   SDY2583_CP22_STEP4_composite_scores.RData
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

step4_rdata <- file.path(
  analysis_dir,
  "11_RData",
  "SDY2583_CP22_STEP4_composite_scores.RData"
)

if (!file.exists(step4_rdata)) {
  stop("CP22 Step 4 RData bulunamadı: ", step4_rdata)
}

out_dir <- file.path(analysis_dir, "06_age_sensitivity_caliper_interaction")
rdata_dir <- file.path(analysis_dir, "11_RData")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rdata_dir, recursive = TRUE, showWarnings = FALSE)

load(step4_rdata)

if (!exists("analysis_df")) {
  stop("Step 4 RData içinde analysis_df bulunamadı.")
}

if (!exists("primary_composite_results")) {
  stop("Step 4 RData içinde primary_composite_results bulunamadı.")
}

dat <- analysis_df

# ------------------------------------------------------------
# 3. Target composite scores
# ------------------------------------------------------------

target_scores <- primary_composite_results$composite_score
target_scores <- target_scores[target_scores %in% names(dat)]

if (length(target_scores) == 0) {
  stop("Composite score columns bulunamadı.")
}

# ------------------------------------------------------------
# 4. Modeling helpers
# ------------------------------------------------------------

run_lm_score <- function(df, score_name, model_label = "model") {

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
    lm(value ~ disease_group + age_for_model + sex, data = model_df),
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

  out <- tt %>% filter(term == disease_term)

  desc <- model_df %>%
    group_by(disease_group) %>%
    summarise(
      n = n(),
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

  tibble(
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
    direction = case_when(
      out$estimate[1] > 0 ~ "higher_in_cancer",
      out$estimate[1] < 0 ~ "lower_in_cancer",
      TRUE ~ "no_direction"
    )
  )
}

run_interaction_score <- function(df, score_name) {

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
      disease_group = factor(
        as.character(disease_group),
        levels = c("Healthy control", "Cancer patient")
      ),
      sex = droplevels(factor(as.character(sex))),
      age_z = as.numeric(scale(age_for_model))
    )

  n_model <- nrow(model_df)
  n_healthy <- sum(model_df$disease_group == "Healthy control")
  n_cancer <- sum(model_df$disease_group == "Cancer patient")

  if (n_model < 50 || n_healthy < 20 || n_cancer < 20) {
    return(NULL)
  }

  fit <- tryCatch(
    lm(value ~ disease_group * age_z + sex, data = model_df),
    error = function(e) NULL
  )

  if (is.null(fit)) return(NULL)

  tt <- tryCatch(
    broom::tidy(fit, conf.int = TRUE),
    error = function(e) NULL
  )

  if (is.null(tt)) return(NULL)

  interaction_term <- "disease_groupCancer patient:age_z"
  if (!(interaction_term %in% tt$term)) return(NULL)

  out <- tt %>% filter(term == interaction_term)

  tibble(
    composite_score = score_name,
    n_model = n_model,
    n_healthy = n_healthy,
    n_cancer = n_cancer,
    beta_interaction_cancer_by_age_z = out$estimate[1],
    conf_low = out$conf.low[1],
    conf_high = out$conf.high[1],
    p_value = out$p.value[1],
    direction = case_when(
      out$estimate[1] > 0 ~ "stronger_cancer_effect_with_older_age",
      out$estimate[1] < 0 ~ "weaker_cancer_effect_with_older_age",
      TRUE ~ "no_direction"
    )
  )
}

# ------------------------------------------------------------
# 5. Same-sex nearest-age caliper matching without replacement
# ------------------------------------------------------------

make_same_sex_age_matched <- function(df, caliper_years = 5) {

  base_df <- df %>%
    filter(
      feature_ok == TRUE,
      !is.na(disease_group),
      !is.na(age_for_model),
      !is.na(sex)
    ) %>%
    mutate(
      disease_group_chr = as.character(disease_group),
      sex_chr = as.character(sex)
    ) %>%
    filter(disease_group_chr %in% c("Healthy control", "Cancer patient")) %>%
    mutate(
      row_id = row_number()
    )

  cancer_df <- base_df %>%
    filter(disease_group_chr == "Cancer patient") %>%
    arrange(age_for_model)

  control_df <- base_df %>%
    filter(disease_group_chr == "Healthy control") %>%
    arrange(age_for_model)

  used_control_rows <- integer(0)
  pair_rows <- list()
  pair_id <- 0

  for (i in seq_len(nrow(cancer_df))) {

    this_cancer <- cancer_df[i, , drop = FALSE]

    eligible_controls <- control_df %>%
      filter(
        !(row_id %in% used_control_rows),
        sex_chr == this_cancer$sex_chr[1]
      ) %>%
      mutate(
        abs_age_diff = abs(age_for_model - this_cancer$age_for_model[1])
      ) %>%
      filter(abs_age_diff <= caliper_years) %>%
      arrange(abs_age_diff, age_for_model)

    if (nrow(eligible_controls) == 0) {
      next
    }

    chosen_control <- eligible_controls[1, , drop = FALSE]

    used_control_rows <- c(used_control_rows, chosen_control$row_id[1])
    pair_id <- pair_id + 1

    this_cancer$matched_pair_id <- pair_id
    this_cancer$matched_abs_age_diff <- chosen_control$abs_age_diff[1]

    chosen_control$matched_pair_id <- pair_id
    chosen_control$matched_abs_age_diff <- chosen_control$abs_age_diff[1]

    pair_rows[[length(pair_rows) + 1]] <- this_cancer
    pair_rows[[length(pair_rows) + 1]] <- chosen_control
  }

  if (length(pair_rows) == 0) {
    return(tibble())
  }

  matched_df <- bind_rows(pair_rows) %>%
    mutate(
      disease_group = factor(
        disease_group_chr,
        levels = c("Healthy control", "Cancer patient")
      ),
      sex = factor(sex_chr)
    )

  matched_df
}

summarise_matched <- function(matched_df, caliper_years) {

  if (nrow(matched_df) == 0) {
    return(tibble(
      caliper_years = caliper_years,
      n_rows = 0,
      n_pairs = 0,
      mean_abs_age_diff = NA_real_,
      median_abs_age_diff = NA_real_,
      max_abs_age_diff = NA_real_,
      n_healthy = 0,
      n_cancer = 0
    ))
  }

  pair_diff <- matched_df %>%
    distinct(matched_pair_id, matched_abs_age_diff)

  tibble(
    caliper_years = caliper_years,
    n_rows = nrow(matched_df),
    n_pairs = n_distinct(matched_df$matched_pair_id),
    mean_abs_age_diff = mean(pair_diff$matched_abs_age_diff, na.rm = TRUE),
    median_abs_age_diff = median(pair_diff$matched_abs_age_diff, na.rm = TRUE),
    max_abs_age_diff = max(pair_diff$matched_abs_age_diff, na.rm = TRUE),
    n_healthy = sum(matched_df$disease_group == "Healthy control"),
    n_cancer = sum(matched_df$disease_group == "Cancer patient")
  )
}

# ------------------------------------------------------------
# 6. Build matched datasets
# ------------------------------------------------------------

matched_5y <- make_same_sex_age_matched(dat, caliper_years = 5)
matched_10y <- make_same_sex_age_matched(dat, caliper_years = 10)

matched_summary <- bind_rows(
  summarise_matched(matched_5y, 5),
  summarise_matched(matched_10y, 10)
)

# ------------------------------------------------------------
# 7. Matched model results
# ------------------------------------------------------------

matched_5y_results <- bind_rows(
  lapply(target_scores, function(cs) {
    run_lm_score(matched_5y, cs, model_label = "same_sex_age_matched_5y")
  })
) %>%
  mutate(
    caliper_years = 5,
    FDR_global = p.adjust(p_value, method = "BH")
  ) %>%
  arrange(FDR_global, p_value)

matched_10y_results <- bind_rows(
  lapply(target_scores, function(cs) {
    run_lm_score(matched_10y, cs, model_label = "same_sex_age_matched_10y")
  })
) %>%
  mutate(
    caliper_years = 10,
    FDR_global = p.adjust(p_value, method = "BH")
  ) %>%
  arrange(FDR_global, p_value)

matched_results_all <- bind_rows(
  matched_5y_results,
  matched_10y_results
)

# ------------------------------------------------------------
# 8. Disease-by-age interaction
# ------------------------------------------------------------

primary_df <- dat %>%
  filter(
    feature_ok == TRUE,
    !is.na(disease_group),
    !is.na(age_for_model),
    !is.na(sex)
  )

interaction_results <- bind_rows(
  lapply(target_scores, function(cs) {
    run_interaction_score(primary_df, cs)
  })
) %>%
  mutate(
    FDR_global = p.adjust(p_value, method = "BH")
  ) %>%
  arrange(FDR_global, p_value)

# ------------------------------------------------------------
# 9. Direction preservation summary
# ------------------------------------------------------------

matched_direction_summary <- primary_composite_results %>%
  select(
    composite_score,
    primary_beta = beta_cancer_vs_healthy,
    primary_direction = direction,
    primary_FDR_global = FDR_global
  ) %>%
  left_join(
    matched_5y_results %>%
      select(
        composite_score,
        matched_5y_beta = beta_cancer_vs_healthy,
        matched_5y_direction = direction,
        matched_5y_FDR_global = FDR_global
      ),
    by = "composite_score"
  ) %>%
  left_join(
    matched_10y_results %>%
      select(
        composite_score,
        matched_10y_beta = beta_cancer_vs_healthy,
        matched_10y_direction = direction,
        matched_10y_FDR_global = FDR_global
      ),
    by = "composite_score"
  ) %>%
  mutate(
    direction_preserved_5y = primary_direction == matched_5y_direction,
    direction_preserved_10y = primary_direction == matched_10y_direction,
    direction_preserved_both_matched = direction_preserved_5y == TRUE & direction_preserved_10y == TRUE,
    FDR_global_preserved_5y = primary_FDR_global < 0.05 & matched_5y_FDR_global < 0.05,
    FDR_global_preserved_10y = primary_FDR_global < 0.05 & matched_10y_FDR_global < 0.05,
    FDR_global_preserved_both_matched = FDR_global_preserved_5y == TRUE & FDR_global_preserved_10y == TRUE
  ) %>%
  arrange(primary_FDR_global)

overall_age_sensitivity_summary <- tibble(
  n_primary_model_subjects = nrow(primary_df),
  n_primary_healthy = sum(primary_df$disease_group == "Healthy control"),
  n_primary_cancer = sum(primary_df$disease_group == "Cancer patient"),
  n_5y_pairs = ifelse(nrow(matched_5y) > 0, n_distinct(matched_5y$matched_pair_id), 0),
  n_10y_pairs = ifelse(nrow(matched_10y) > 0, n_distinct(matched_10y$matched_pair_id), 0),
  n_scores_tested = length(target_scores),
  n_5y_FDR_lt_0p05 = sum(matched_5y_results$FDR_global < 0.05, na.rm = TRUE),
  n_10y_FDR_lt_0p05 = sum(matched_10y_results$FDR_global < 0.05, na.rm = TRUE),
  n_direction_preserved_both_matched = sum(matched_direction_summary$direction_preserved_both_matched == TRUE, na.rm = TRUE),
  n_FDR_preserved_both_matched = sum(matched_direction_summary$FDR_global_preserved_both_matched == TRUE, na.rm = TRUE),
  n_interaction_FDR_lt_0p05 = sum(interaction_results$FDR_global < 0.05, na.rm = TRUE)
)

# ------------------------------------------------------------
# 10. Save outputs
# ------------------------------------------------------------

write_csv(
  matched_summary,
  file.path(out_dir, "SDY2583_CP22_matched_dataset_summary_STEP4B.csv")
)

write_csv(
  matched_5y_results,
  file.path(out_dir, "SDY2583_CP22_same_sex_age_matched_5y_results_STEP4B.csv")
)

write_csv(
  matched_10y_results,
  file.path(out_dir, "SDY2583_CP22_same_sex_age_matched_10y_results_STEP4B.csv")
)

write_csv(
  matched_results_all,
  file.path(out_dir, "SDY2583_CP22_same_sex_age_matched_results_ALL_STEP4B.csv")
)

write_csv(
  interaction_results,
  file.path(out_dir, "SDY2583_CP22_disease_by_age_interaction_results_STEP4B.csv")
)

write_csv(
  matched_direction_summary,
  file.path(out_dir, "SDY2583_CP22_matched_direction_sensitivity_summary_STEP4B.csv")
)

write_csv(
  overall_age_sensitivity_summary,
  file.path(out_dir, "SDY2583_CP22_overall_age_sensitivity_summary_STEP4B.csv")
)

save(
  matched_5y,
  matched_10y,
  matched_summary,
  matched_5y_results,
  matched_10y_results,
  matched_results_all,
  interaction_results,
  matched_direction_summary,
  overall_age_sensitivity_summary,
  target_scores,
  file = file.path(
    rdata_dir,
    "SDY2583_CP22_STEP4B_age_sensitivity_caliper_interaction.RData"
  )
)

# ------------------------------------------------------------
# 11. Console output
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SDY2583 CP22 STEP 4B COMPLETE: AGE SENSITIVITY / MATCHING / INTERACTION\n")
cat("============================================================\n")

cat("\nMatched dataset summary:\n")
print(as.data.frame(matched_summary), row.names = FALSE)

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
