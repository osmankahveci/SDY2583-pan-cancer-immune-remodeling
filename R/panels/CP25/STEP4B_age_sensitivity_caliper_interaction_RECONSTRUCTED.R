# CP25 Step 4B reconstructed archive-compatible age sensitivity analyses.
# The archived family contains all 77 features plus all 9 composite scores.

rm(list = ls())
source(file.path(
  Sys.getenv("SDY2583_REPO_ROOT", unset = "."),
  "R", "shared", "reconstructed_panel_framework.R"
))
rp_install_and_load(c(
  "dplyr", "readr", "tibble", "purrr", "broom", "tidyr"
))
source(file.path(
  sd_repo_root(), "R", "panels", "CP25", "MANIFEST_RECONSTRUCTED.R"
))

analysis_dir <- sd_analysis_dir("CP25")
rdata_dir <- file.path(analysis_dir, "11_RData")
out_dir <- file.path(analysis_dir, "10_age_sensitivity")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

load(file.path(
  rdata_dir,
  "SDY2583_CP25_STEP4_composite_scores_RECONSTRUCTED.RData"
))

score_names <- c(
  names(CP25_MANIFEST$score_definitions),
  CP25_MANIFEST$integrated_score
)
outcomes <- unique(c(CP25_FEATURES, score_names))
stopifnot(length(outcomes) == 86L)

cp25_fit_outcome <- function(data, feature_name) {
  work <- data |>
    dplyr::transmute(
      value = suppressWarnings(as.numeric(.data[[feature_name]])),
      disease_group = factor(
        as.character(disease_group),
        levels = c("Healthy control", "Cancer patient")
      ),
      age_for_model = suppressWarnings(as.numeric(age_for_model)),
      sex = factor(as.character(sex))
    ) |>
    dplyr::filter(stats::complete.cases(.)) |>
    dplyr::mutate(sex = droplevels(sex))

  fit <- stats::lm(
    value ~ disease_group + age_for_model + sex,
    data = work
  )
  tidy_fit <- broom::tidy(fit, conf.int = TRUE)
  term <- tidy_fit |>
    dplyr::filter(term == "disease_groupCancer patient")

  healthy <- work$value[work$disease_group == "Healthy control"]
  cancer <- work$value[work$disease_group == "Cancer patient"]
  beta <- term$estimate[1]

  tibble::tibble(
    feature = feature_name,
    n_model = stats::nobs(fit),
    n_healthy = length(healthy),
    n_cancer = length(cancer),
    healthy_mean = mean(healthy),
    cancer_mean = mean(cancer),
    beta_cancer_vs_healthy = beta,
    ci_low = term$conf.low[1],
    ci_high = term$conf.high[1],
    p_value = term$p.value[1],
    direction = dplyr::case_when(
      beta > 0 ~ "higher_in_cancer",
      beta < 0 ~ "lower_in_cancer",
      TRUE ~ "no_difference"
    )
  )
}

cp25_match <- function(data, caliper_years) {
  base <- data |>
    dplyr::filter(
      !is.na(age_for_model),
      as.character(sex) %in% c("Female", "Male"),
      disease_group %in% c("Healthy control", "Cancer patient")
    ) |>
    dplyr::mutate(
      disease_group = factor(
        as.character(disease_group),
        levels = c("Healthy control", "Cancer patient")
      ),
      sex = factor(as.character(sex), levels = c("Female", "Male"))
    )

  cancer <- base |>
    dplyr::filter(disease_group == "Cancer patient") |>
    dplyr::arrange(sex, age_for_model, subject_id)
  healthy <- base |>
    dplyr::filter(disease_group == "Healthy control") |>
    dplyr::arrange(sex, age_for_model, subject_id)

  used <- rep(FALSE, nrow(healthy))
  pair_data <- list()
  pair_table <- list()
  pair_number <- 0L

  for (index in seq_len(nrow(cancer))) {
    cancer_row <- cancer[index, , drop = FALSE]
    candidates <- which(
      !used & as.character(healthy$sex) == as.character(cancer_row$sex[1])
    )
    if (length(candidates) == 0L) next

    candidate_data <- healthy[candidates, , drop = FALSE] |>
      dplyr::mutate(
        candidate_row = candidates,
        abs_age_diff = abs(age_for_model - cancer_row$age_for_model[1])
      ) |>
      dplyr::filter(abs_age_diff <= caliper_years) |>
      dplyr::arrange(abs_age_diff, age_for_model, subject_id)

    if (nrow(candidate_data) == 0L) next

    healthy_row <- candidate_data[1, , drop = FALSE]
    healthy_index <- healthy_row$candidate_row[1]
    used[healthy_index] <- TRUE
    pair_number <- pair_number + 1L
    pair_id <- paste0("pair_", pair_number)
    age_difference <- healthy_row$abs_age_diff[1]

    cancer_row$pair_id <- pair_id
    cancer_row$abs_age_diff <- age_difference
    healthy_selected <- healthy[healthy_index, , drop = FALSE]
    healthy_selected$pair_id <- pair_id
    healthy_selected$abs_age_diff <- age_difference

    # Preserve archived row order: cancer then healthy within each pair.
    pair_data[[length(pair_data) + 1L]] <- cancer_row
    pair_data[[length(pair_data) + 1L]] <- healthy_selected

    pair_table[[pair_number]] <- dplyr::bind_rows(
      cancer_row |>
        dplyr::transmute(
          pair_id,
          subject_id,
          disease_group = as.character(disease_group),
          sex = as.character(sex),
          age_for_model,
          abs_age_diff
        ),
      healthy_selected |>
        dplyr::transmute(
          pair_id,
          subject_id,
          disease_group = as.character(disease_group),
          sex = as.character(sex),
          age_for_model,
          abs_age_diff
        )
    )
  }

  list(
    data = dplyr::bind_rows(pair_data),
    pairs = dplyr::bind_rows(pair_table)
  )
}

age_stratified <- purrr::map_dfr(
  levels(scored_data$age_group_for_model),
  function(age_group) {
    group_data <- scored_data |>
      dplyr::filter(age_group_for_model == age_group)
    result <- purrr::map_dfr(
      outcomes,
      function(feature_name) cp25_fit_outcome(group_data, feature_name)
    )
    result$age_group_for_model <- age_group
    result |>
      dplyr::mutate(FDR_within_age_group = p.adjust(p_value, method = "BH"))
  }
)

matched_5 <- cp25_match(scored_data, 5)
matched_10 <- cp25_match(scored_data, 10)

fit_matched <- function(match_result, caliper_years) {
  purrr::map_dfr(
    outcomes,
    function(feature_name) cp25_fit_outcome(match_result$data, feature_name)
  ) |>
    dplyr::mutate(
      caliper_years = caliper_years,
      FDR_within_caliper = p.adjust(p_value, method = "BH")
    )
}

matched_statistics <- dplyr::bind_rows(
  fit_matched(matched_5, 5),
  fit_matched(matched_10, 10)
)

match_summary <- dplyr::bind_rows(
  tibble::tibble(
    caliper_years = 5,
    n_pairs = nrow(matched_5$pairs) / 2,
    n_rows = nrow(matched_5$pairs),
    mean_abs_age_diff = mean(matched_5$pairs$abs_age_diff, na.rm = TRUE),
    median_abs_age_diff = stats::median(
      matched_5$pairs$abs_age_diff,
      na.rm = TRUE
    ),
    max_abs_age_diff = max(matched_5$pairs$abs_age_diff, na.rm = TRUE)
  ),
  tibble::tibble(
    caliper_years = 10,
    n_pairs = nrow(matched_10$pairs) / 2,
    n_rows = nrow(matched_10$pairs),
    mean_abs_age_diff = mean(matched_10$pairs$abs_age_diff, na.rm = TRUE),
    median_abs_age_diff = stats::median(
      matched_10$pairs$abs_age_diff,
      na.rm = TRUE
    ),
    max_abs_age_diff = max(matched_10$pairs$abs_age_diff, na.rm = TRUE)
  )
)

cp25_interaction <- function(data, feature_name) {
  work <- data |>
    dplyr::transmute(
      value = suppressWarnings(as.numeric(.data[[feature_name]])),
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
      # Archive rule: standardize age within each outcome's complete-case set.
      age_z = as.numeric(scale(age_for_model))
    )

  fit <- stats::lm(
    value ~ disease_group * age_z + sex,
    data = work
  )
  tidy_fit <- broom::tidy(fit, conf.int = TRUE)
  term <- tidy_fit |>
    dplyr::filter(term == "disease_groupCancer patient:age_z")
  beta <- term$estimate[1]

  tibble::tibble(
    feature = feature_name,
    n_model = stats::nobs(fit),
    beta_interaction = beta,
    ci_low = term$conf.low[1],
    ci_high = term$conf.high[1],
    p_value = term$p.value[1],
    direction = dplyr::case_when(
      beta > 0 ~ "stronger_cancer_effect_with_older_age",
      beta < 0 ~ "weaker_cancer_effect_with_older_age",
      TRUE ~ "no_interaction_direction"
    )
  )
}

interaction_statistics <- purrr::map_dfr(
  outcomes,
  function(feature_name) cp25_interaction(scored_data, feature_name)
) |>
  dplyr::mutate(FDR = p.adjust(p_value, method = "BH"))

readr::write_csv(
  age_stratified,
  file.path(out_dir, "SDY2583_CP25_age_stratified_statistics_RECONSTRUCTED.csv")
)
readr::write_csv(
  matched_statistics,
  file.path(
    out_dir,
    "SDY2583_CP25_age_matched_sensitivity_statistics_RECONSTRUCTED.csv"
  )
)
readr::write_csv(
  match_summary,
  file.path(
    out_dir,
    "SDY2583_CP25_age_matched_sensitivity_summary_RECONSTRUCTED.csv"
  )
)
readr::write_csv(
  matched_5$pairs,
  file.path(
    out_dir,
    "SDY2583_CP25_same_sex_age_matched_pairs_5yr_RECONSTRUCTED.csv"
  )
)
readr::write_csv(
  matched_10$pairs,
  file.path(
    out_dir,
    "SDY2583_CP25_same_sex_age_matched_pairs_10yr_RECONSTRUCTED.csv"
  )
)
readr::write_csv(
  interaction_statistics,
  file.path(
    out_dir,
    "SDY2583_CP25_disease_by_age_interaction_statistics_RECONSTRUCTED.csv"
  )
)

save(
  age_stratified,
  matched_statistics,
  match_summary,
  matched_5,
  matched_10,
  interaction_statistics,
  outcomes,
  file = file.path(
    rdata_dir,
    "SDY2583_CP25_STEP4B_age_sensitivity_RECONSTRUCTED.RData"
  )
)
print(match_summary)
