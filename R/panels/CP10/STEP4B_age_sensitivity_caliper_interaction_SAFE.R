# ============================================================
# SDY2583 CP10
# STEP 4B SAFE: Age sensitivity, same-sex age-caliper matching,
# and disease-by-age interaction for CP10 composite scores
# ============================================================

rm(list = ls())
source(file.path(Sys.getenv("SDY2583_REPO_ROOT", unset = "."), "R", "shared", "bootstrap.R"))

cran_pkgs <- c("dplyr", "readr", "tibble", "broom")
for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(broom)
})

filter <- dplyr::filter
select <- dplyr::select
mutate <- dplyr::mutate
arrange <- dplyr::arrange
summarise <- dplyr::summarise
group_by <- dplyr::group_by
ungroup <- dplyr::ungroup
bind_rows <- dplyr::bind_rows
case_when <- dplyr::case_when
n <- dplyr::n

analysis_dir <- sd_analysis_dir("CP10")
rdata_dir <- file.path(analysis_dir, "11_RData")
out_dir <- file.path(analysis_dir, "06_age_sensitivity_caliper_interaction")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

step4_rdata <- file.path(rdata_dir, "SDY2583_CP10_STEP4_composite_scores.RData")
if (!file.exists(step4_rdata)) stop("Step 4 RData bulunamadı: ", step4_rdata)

load(step4_rdata)

if (!exists("cp10_scores_data")) stop("cp10_scores_data bulunamadı.")
if (!exists("primary_composite_results")) stop("primary_composite_results bulunamadı.")
if (!exists("composite_score_cols")) stop("composite_score_cols bulunamadı.")

target_scores <- composite_score_cols

make_same_sex_age_matched <- function(df, caliper_years = 5) {

  base_df <- df %>%
    filter(
      model_ready_age_sex == TRUE,
      !is.na(age_for_model),
      !is.na(sex_binary),
      !is.na(disease_group),
      disease_group %in% c("Healthy control", "Cancer patient")
    ) %>%
    mutate(
      disease_group = as.character(disease_group),
      sex_binary = as.character(sex_binary)
    )

  healthy <- base_df %>%
    filter(disease_group == "Healthy control") %>%
    arrange(sex_binary, age_for_model, subject_id)

  cancer <- base_df %>%
    filter(disease_group == "Cancer patient") %>%
    arrange(sex_binary, age_for_model, subject_id)

  used_healthy <- character()
  matched_pairs <- list()
  pair_id <- 0L

  for (i in seq_len(nrow(cancer))) {
    cc <- cancer[i, ]

    candidate <- healthy %>%
      filter(
        !(subject_id %in% used_healthy),
        sex_binary == cc$sex_binary,
        abs(age_for_model - cc$age_for_model) <= caliper_years
      ) %>%
      mutate(abs_age_diff = abs(age_for_model - cc$age_for_model)) %>%
      arrange(abs_age_diff, age_for_model, subject_id)

    if (nrow(candidate) > 0) {
      hh <- candidate[1, ]
      pair_id <- pair_id + 1L
      used_healthy <- c(used_healthy, hh$subject_id)

      matched_pairs[[pair_id]] <- bind_rows(
        hh %>% mutate(pair_id = pair_id, matched_role = "healthy", matched_age_diff = abs_age_diff),
        cc %>% mutate(pair_id = pair_id, matched_role = "cancer", matched_age_diff = abs(hh$age_for_model - cc$age_for_model))
      )
    }
  }

  if (length(matched_pairs) == 0) return(tibble())

  bind_rows(matched_pairs) %>%
    mutate(
      disease_group = factor(disease_group, levels = c("Healthy control", "Cancer patient")),
      sex = droplevels(factor(as.character(sex))),
      sex_binary = droplevels(factor(as.character(sex_binary)))
    )
}

run_score_model <- function(df, score_name, model_label, caliper_years = NA_real_) {

  model_df <- df %>%
    transmute(
      value = suppressWarnings(as.numeric(.data[[score_name]])),
      disease_group = disease_group,
      age_for_model = suppressWarnings(as.numeric(age_for_model)),
      sex = sex
    ) %>%
    filter(!is.na(value), !is.na(disease_group), !is.na(age_for_model), !is.na(sex)) %>%
    mutate(
      disease_group = factor(as.character(disease_group), levels = c("Healthy control", "Cancer patient")),
      sex = droplevels(factor(as.character(sex)))
    )

  n_model <- nrow(model_df)
  n_healthy <- sum(model_df$disease_group == "Healthy control")
  n_cancer <- sum(model_df$disease_group == "Cancer patient")
  if (n_model < 50 || n_healthy < 20 || n_cancer < 20) return(NULL)

  fit <- tryCatch(lm(value ~ disease_group + age_for_model + sex, data = model_df), error = function(e) NULL)
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
    n_model = n_model,
    n_healthy = n_healthy,
    n_cancer = n_cancer,
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

run_model_set <- function(df, model_label, caliper_years = NA_real_) {
  bind_rows(lapply(target_scores, function(sc) run_score_model(df, sc, model_label, caliper_years))) %>%
    mutate(FDR_global = p.adjust(p_value, method = "BH")) %>%
    arrange(FDR_global, p_value)
}

matched_5y <- make_same_sex_age_matched(cp10_scores_data, 5)
matched_10y <- make_same_sex_age_matched(cp10_scores_data, 10)

pair_summary_one <- function(matched_df, caliper_years) {
  if (nrow(matched_df) == 0) return(tibble())
  matched_df %>%
    group_by(pair_id) %>%
    summarise(
      caliper_years = caliper_years,
      abs_age_diff = abs(diff(age_for_model)),
      .groups = "drop"
    )
}

matched_summary <- bind_rows(
  pair_summary_one(matched_5y, 5),
  pair_summary_one(matched_10y, 10)
) %>%
  group_by(caliper_years) %>%
  summarise(
    n_rows = n() * 2,
    n_pairs = n(),
    mean_abs_age_diff = mean(abs_age_diff, na.rm = TRUE),
    median_abs_age_diff = median(abs_age_diff, na.rm = TRUE),
    max_abs_age_diff = max(abs_age_diff, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(n_healthy = n_pairs, n_cancer = n_pairs)

matched_5y_results <- run_model_set(matched_5y, "same_sex_age_matched_5y", 5)
matched_10y_results <- run_model_set(matched_10y, "same_sex_age_matched_10y", 10)

run_interaction_one <- function(df, score_name) {

  model_df <- df %>%
    transmute(
      value = suppressWarnings(as.numeric(.data[[score_name]])),
      disease_group = disease_group,
      age_for_model = suppressWarnings(as.numeric(age_for_model)),
      sex = sex
    ) %>%
    filter(!is.na(value), !is.na(disease_group), !is.na(age_for_model), !is.na(sex)) %>%
    mutate(
      disease_group = factor(as.character(disease_group), levels = c("Healthy control", "Cancer patient")),
      sex = droplevels(factor(as.character(sex))),
      age_z = as.numeric(scale(age_for_model))
    )

  n_model <- nrow(model_df)
  n_healthy <- sum(model_df$disease_group == "Healthy control")
  n_cancer <- sum(model_df$disease_group == "Cancer patient")
  if (n_model < 50 || n_healthy < 20 || n_cancer < 20) return(NULL)

  fit <- tryCatch(lm(value ~ disease_group * age_z + sex, data = model_df), error = function(e) NULL)
  if (is.null(fit)) return(NULL)

  tt <- tryCatch(broom::tidy(fit, conf.int = TRUE), error = function(e) NULL)
  if (is.null(tt)) return(NULL)

  term <- "disease_groupCancer patient:age_z"
  if (!(term %in% tt$term)) return(NULL)
  out <- tt %>% filter(term == !!term)

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
      TRUE ~ "no_interaction_direction"
    )
  )
}

interaction_results <- bind_rows(lapply(target_scores, function(sc) {
  run_interaction_one(cp10_scores_data %>% filter(model_ready_age_sex == TRUE), sc)
})) %>%
  mutate(FDR_global = p.adjust(p_value, method = "BH")) %>%
  arrange(FDR_global, p_value)

get_val <- function(df, score_name, col_name) {
  x <- df %>%
    filter(composite_score == score_name) %>%
    pull(!!rlang::sym(col_name))
  if (length(x) == 0) return(NA)
  x[1]
}

matched_direction_summary <- bind_rows(lapply(target_scores, function(sc) {

  p_row <- primary_composite_results %>% filter(composite_score == sc)

  m5_beta <- as.numeric(get_val(matched_5y_results, sc, "beta_cancer_vs_healthy"))
  m5_dir <- as.character(get_val(matched_5y_results, sc, "direction"))
  m5_fdr <- as.numeric(get_val(matched_5y_results, sc, "FDR_global"))

  m10_beta <- as.numeric(get_val(matched_10y_results, sc, "beta_cancer_vs_healthy"))
  m10_dir <- as.character(get_val(matched_10y_results, sc, "direction"))
  m10_fdr <- as.numeric(get_val(matched_10y_results, sc, "FDR_global"))

  tibble(
    composite_score = sc,
    primary_beta = p_row$beta_cancer_vs_healthy[1],
    primary_direction = p_row$direction[1],
    primary_FDR_global = p_row$FDR_global[1],
    matched_5y_beta = m5_beta,
    matched_5y_direction = m5_dir,
    matched_5y_FDR_global = m5_fdr,
    matched_10y_beta = m10_beta,
    matched_10y_direction = m10_dir,
    matched_10y_FDR_global = m10_fdr,
    direction_preserved_5y = !is.na(m5_dir) && m5_dir == p_row$direction[1],
    direction_preserved_10y = !is.na(m10_dir) && m10_dir == p_row$direction[1],
    direction_preserved_both_matched = direction_preserved_5y & direction_preserved_10y,
    FDR_global_preserved_5y = !is.na(m5_fdr) && m5_fdr < 0.05,
    FDR_global_preserved_10y = !is.na(m10_fdr) && m10_fdr < 0.05,
    FDR_global_preserved_both_matched = FDR_global_preserved_5y & FDR_global_preserved_10y
  )
})) %>%
  arrange(primary_FDR_global)

overall_age_sensitivity_summary <- tibble(
  n_primary_model_subjects = nrow(cp10_scores_data %>% filter(model_ready_age_sex == TRUE)),
  n_primary_healthy = sum(cp10_scores_data$model_ready_age_sex == TRUE & cp10_scores_data$disease_group == "Healthy control", na.rm = TRUE),
  n_primary_cancer = sum(cp10_scores_data$model_ready_age_sex == TRUE & cp10_scores_data$disease_group == "Cancer patient", na.rm = TRUE),
  n_5y_pairs = ifelse(any(matched_summary$caliper_years == 5), matched_summary$n_pairs[matched_summary$caliper_years == 5], 0),
  n_10y_pairs = ifelse(any(matched_summary$caliper_years == 10), matched_summary$n_pairs[matched_summary$caliper_years == 10], 0),
  n_scores_tested = length(target_scores),
  n_5y_FDR_lt_0p05 = sum(matched_5y_results$FDR_global < 0.05, na.rm = TRUE),
  n_10y_FDR_lt_0p05 = sum(matched_10y_results$FDR_global < 0.05, na.rm = TRUE),
  n_direction_preserved_both_matched = sum(matched_direction_summary$direction_preserved_both_matched == TRUE, na.rm = TRUE),
  n_FDR_preserved_both_matched = sum(matched_direction_summary$FDR_global_preserved_both_matched == TRUE, na.rm = TRUE),
  n_interaction_FDR_lt_0p05 = sum(interaction_results$FDR_global < 0.05, na.rm = TRUE)
)

write_csv(matched_summary, file.path(out_dir, "SDY2583_CP10_matched_dataset_summary_STEP4B.csv"))
write_csv(overall_age_sensitivity_summary, file.path(out_dir, "SDY2583_CP10_overall_age_sensitivity_summary_STEP4B.csv"))
write_csv(matched_5y_results, file.path(out_dir, "SDY2583_CP10_same_sex_age_matched_5y_results_STEP4B.csv"))
write_csv(matched_10y_results, file.path(out_dir, "SDY2583_CP10_same_sex_age_matched_10y_results_STEP4B.csv"))
write_csv(interaction_results, file.path(out_dir, "SDY2583_CP10_disease_by_age_interaction_results_STEP4B.csv"))
write_csv(matched_direction_summary, file.path(out_dir, "SDY2583_CP10_matched_direction_sensitivity_summary_STEP4B.csv"))

if (nrow(matched_5y) > 0) {
  write_csv(
    matched_5y %>% select(subject_id, disease_group, sex_binary, age_for_model, pair_id, matched_role, matched_age_diff),
    file.path(out_dir, "SDY2583_CP10_same_sex_age_matched_5y_subjects_STEP4B.csv")
  )
}

if (nrow(matched_10y) > 0) {
  write_csv(
    matched_10y %>% select(subject_id, disease_group, sex_binary, age_for_model, pair_id, matched_role, matched_age_diff),
    file.path(out_dir, "SDY2583_CP10_same_sex_age_matched_10y_subjects_STEP4B.csv")
  )
}

save(
  matched_5y,
  matched_10y,
  matched_summary,
  matched_5y_results,
  matched_10y_results,
  interaction_results,
  matched_direction_summary,
  overall_age_sensitivity_summary,
  file = file.path(rdata_dir, "SDY2583_CP10_STEP4B_age_sensitivity_caliper_interaction.RData")
)

cat("\n============================================================\n")
cat("SDY2583 CP10 STEP 4B COMPLETE: AGE SENSITIVITY / MATCHING / INTERACTION\n")
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
