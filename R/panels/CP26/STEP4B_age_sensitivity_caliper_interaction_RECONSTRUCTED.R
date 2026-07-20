# CP26 Step 4B reconstructed archive-compatible matching and interaction.

rm(list = ls())
source(file.path(
  Sys.getenv("SDY2583_REPO_ROOT", unset = "."),
  "R", "shared", "reconstructed_panel_framework.R"
))
rp_install_and_load(c("dplyr", "readr", "tibble", "purrr", "broom"))
source(file.path(
  sd_repo_root(), "R", "panels", "CP26", "MANIFEST_RECONSTRUCTED.R"
))

analysis_dir <- sd_analysis_dir("CP26")
rdata_dir <- file.path(analysis_dir, "11_RData")
out_dir <- file.path(analysis_dir, "07_age_sensitivity_matching_interaction")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

load(file.path(
  rdata_dir,
  "SDY2583_CP26_STEP4_composite_scores_RECONSTRUCTED.RData"
))

outcomes <- unique(c(
  CP26_MANIFEST$targeted_outcomes,
  names(CP26_MANIFEST$score_definitions),
  CP26_MANIFEST$integrated_score
))
stopifnot(length(outcomes) == 20L)

prepare_match <- function(data, caliper_years) {
  matched <- rp_greedy_match(data, caliper_years)
  matched_data <- matched$data |>
    dplyr::mutate(
      matched_pair_id = paste0("pair_", pair_id),
      matching_caliper = caliper_years
    ) |>
    dplyr::group_by(matched_pair_id) |>
    dplyr::mutate(
      abs_age_diff = max(age_for_model, na.rm = TRUE) -
        min(age_for_model, na.rm = TRUE)
    ) |>
    dplyr::ungroup() |>
    dplyr::select(-pair_id)

  result_table <- purrr::map_dfr(outcomes, function(variable_name) {
    rp_fit_one(
      matched$data,
      variable_name,
      subset_label = paste0("same_sex_age_matched_", caliper_years, "y"),
      binary_sex = FALSE,
      event_col = CP26_MANIFEST$primary_event_col,
      min_events = 0
    )
  }) |>
    dplyr::mutate(FDR_all = p.adjust(p_value, method = "BH")) |>
    dplyr::transmute(
      variable = feature,
      n_model,
      n_pairs = n_model / 2,
      healthy_mean,
      cancer_mean,
      beta_cancer_vs_healthy,
      ci_low,
      ci_high,
      p_value,
      FDR_all,
      direction
    )

  list(
    data = matched_data,
    results = result_table,
    summary = tibble::tibble(
      caliper_years = caliper_years,
      n_pairs = nrow(matched$pairs),
      n_rows = nrow(matched_data),
      mean_abs_age_difference = mean(
        matched$pairs$abs_age_difference,
        na.rm = TRUE
      ),
      median_abs_age_difference = stats::median(
        matched$pairs$abs_age_difference,
        na.rm = TRUE
      ),
      max_abs_age_difference = max(
        matched$pairs$abs_age_difference,
        na.rm = TRUE
      )
    )
  )
}

matched_5 <- prepare_match(scored_data, 5)
matched_10 <- prepare_match(scored_data, 10)
matching_summary <- dplyr::bind_rows(matched_5$summary, matched_10$summary)

run_interaction <- function(variable_name) {
  work <- scored_data |>
    dplyr::transmute(
      value = suppressWarnings(as.numeric(.data[[variable_name]])),
      disease_group = factor(
        as.character(disease_group),
        levels = c("Healthy control", "Cancer patient")
      ),
      age_for_model = suppressWarnings(as.numeric(age_for_model)),
      sex = factor(as.character(sex))
    ) |>
    dplyr::filter(stats::complete.cases(.)) |>
    dplyr::mutate(
      sex = droplevels(sex),
      age_z = as.numeric(scale(age_for_model))
    )

  fit <- stats::lm(value ~ disease_group * age_z + sex, data = work)
  tidy_fit <- broom::tidy(fit, conf.int = TRUE)
  term <- tidy_fit |>
    dplyr::filter(term == "disease_groupCancer patient:age_z")
  beta <- term$estimate[1]

  tibble::tibble(
    variable = variable_name,
    n_model = stats::nobs(fit),
    beta_interaction = beta,
    ci_low = term$conf.low[1],
    ci_high = term$conf.high[1],
    p_value_interaction = term$p.value[1],
    interaction_direction = dplyr::case_when(
      beta > 0 ~ "stronger_cancer_effect_with_older_age",
      beta < 0 ~ "weaker_cancer_effect_with_older_age",
      TRUE ~ "no_interaction_direction"
    )
  )
}

interaction_results <- purrr::map_dfr(outcomes, run_interaction) |>
  dplyr::mutate(FDR_all = p.adjust(p_value_interaction, method = "BH")) |>
  dplyr::select(
    variable, n_model, beta_interaction, ci_low, ci_high,
    p_value_interaction, FDR_all, interaction_direction
  )

readr::write_csv(
  matched_5$data,
  file.path(
    out_dir,
    "SDY2583_CP26_same_sex_age_matched_5y_STEP4B_RECONSTRUCTED.csv"
  )
)
readr::write_csv(
  matched_10$data,
  file.path(
    out_dir,
    "SDY2583_CP26_same_sex_age_matched_10y_STEP4B_RECONSTRUCTED.csv"
  )
)
readr::write_csv(
  matched_5$results,
  file.path(out_dir, "SDY2583_CP26_matched_5y_results_STEP4B_RECONSTRUCTED.csv")
)
readr::write_csv(
  matched_10$results,
  file.path(out_dir, "SDY2583_CP26_matched_10y_results_STEP4B_RECONSTRUCTED.csv")
)
readr::write_csv(
  dplyr::bind_rows(
    dplyr::mutate(matched_5$results, matching_caliper = 5),
    dplyr::mutate(matched_10$results, matching_caliper = 10)
  ),
  file.path(out_dir, "SDY2583_CP26_matched_results_STEP4B_RECONSTRUCTED.csv")
)
readr::write_csv(
  matching_summary,
  file.path(out_dir, "SDY2583_CP26_matching_summary_STEP4B_RECONSTRUCTED.csv")
)
readr::write_csv(
  interaction_results,
  file.path(
    out_dir,
    "SDY2583_CP26_disease_by_age_interaction_results_STEP4B_RECONSTRUCTED.csv"
  )
)

save(
  outcomes,
  matched_5,
  matched_10,
  matching_summary,
  interaction_results,
  file = file.path(
    rdata_dir,
    "SDY2583_CP26_STEP4B_age_sensitivity_RECONSTRUCTED.RData"
  )
)
print(matching_summary)
