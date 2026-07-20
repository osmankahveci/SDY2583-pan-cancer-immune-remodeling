# CP26 Step 4B reconstructed matching and disease-by-age sensitivity.
# Uses the 20-variable archived target family.

rm(list = ls())
source(file.path(
  Sys.getenv("SDY2583_REPO_ROOT", unset = "."),
  "R", "shared", "reconstructed_panel_framework.R"
))
rp_install_and_load(c("dplyr", "readr", "tibble", "purrr"))
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

result <- rp_age_sensitivity(
  scored_data,
  outcomes,
  event_col = CP26_MANIFEST$primary_event_col,
  min_events = 0,
  calipers = c(5, 10)
)

matched_5_data <- result$matched_data[["5"]] |>
  dplyr::rename(matched_pair_id = pair_id) |>
  dplyr::mutate(
    abs_age_diff = NA_real_,
    matching_caliper = 5
  )
matched_10_data <- result$matched_data[["10"]] |>
  dplyr::rename(matched_pair_id = pair_id) |>
  dplyr::mutate(
    abs_age_diff = NA_real_,
    matching_caliper = 10
  )

matched_5_results <- result$matched |>
  dplyr::filter(grepl("_5y$", subset_label)) |>
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
    FDR_all = fdr,
    direction
  )
matched_10_results <- result$matched |>
  dplyr::filter(grepl("_10y$", subset_label)) |>
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
    FDR_all = fdr,
    direction
  )

interaction_results <- result$interactions |>
  dplyr::transmute(
    variable = feature,
    n_model,
    beta_interaction = interaction_beta,
    ci_low = NA_real_,
    ci_high = NA_real_,
    p_value_interaction = interaction_p,
    FDR_all = interaction_fdr,
    interaction_direction = dplyr::case_when(
      interaction_beta > 0 ~ "stronger_cancer_effect_with_older_age",
      interaction_beta < 0 ~ "weaker_cancer_effect_with_older_age",
      TRUE ~ "no_interaction_direction"
    )
  )

# Confidence intervals are restored using the same outcome-specific models.
add_interaction_ci <- function(table, data) {
  purrr::map_dfr(table$variable, function(variable_name) {
    work <- data |>
      dplyr::transmute(
        value = suppressWarnings(as.numeric(.data[[variable_name]])),
        disease_group = factor(
          as.character(disease_group),
          levels = c("Healthy control", "Cancer patient")
        ),
        age_z = as.numeric(scale(age_for_model)),
        sex = factor(as.character(sex))
      ) |>
      dplyr::filter(stats::complete.cases(.))
    fit <- stats::lm(value ~ disease_group * age_z + sex, data = work)
    term <- grep(
      "disease_group.*:age_z",
      names(stats::coef(fit)),
      value = TRUE
    )[1]
    interval <- stats::confint(fit, parm = term)
    tibble::tibble(
      variable = variable_name,
      ci_low = unname(interval[1]),
      ci_high = unname(interval[2])
    )
  })
}
interaction_ci <- add_interaction_ci(interaction_results, scored_data)
interaction_results <- interaction_results |>
  dplyr::select(-ci_low, -ci_high) |>
  dplyr::left_join(interaction_ci, by = "variable") |>
  dplyr::select(
    variable, n_model, beta_interaction, ci_low, ci_high,
    p_value_interaction, FDR_all, interaction_direction
  )

readr::write_csv(
  result$stratified,
  file.path(
    out_dir,
    "SDY2583_CP26_age_stratified_results_STEP4B_RECONSTRUCTED.csv"
  )
)
readr::write_csv(
  matched_5_data,
  file.path(
    out_dir,
    "SDY2583_CP26_same_sex_age_matched_5y_STEP4B_RECONSTRUCTED.csv"
  )
)
readr::write_csv(
  matched_10_data,
  file.path(
    out_dir,
    "SDY2583_CP26_same_sex_age_matched_10y_STEP4B_RECONSTRUCTED.csv"
  )
)
readr::write_csv(
  matched_5_results,
  file.path(out_dir, "SDY2583_CP26_matched_5y_results_STEP4B_RECONSTRUCTED.csv")
)
readr::write_csv(
  matched_10_results,
  file.path(out_dir, "SDY2583_CP26_matched_10y_results_STEP4B_RECONSTRUCTED.csv")
)
readr::write_csv(
  dplyr::bind_rows(
    dplyr::mutate(matched_5_results, matching_caliper = 5),
    dplyr::mutate(matched_10_results, matching_caliper = 10)
  ),
  file.path(out_dir, "SDY2583_CP26_matched_results_STEP4B_RECONSTRUCTED.csv")
)
readr::write_csv(
  result$match_summary,
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
  result,
  outcomes,
  matched_5_data,
  matched_10_data,
  matched_5_results,
  matched_10_results,
  interaction_results,
  file = file.path(
    rdata_dir,
    "SDY2583_CP26_STEP4B_age_sensitivity_RECONSTRUCTED.RData"
  )
)
print(result$match_summary)
