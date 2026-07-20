# CP28 Step 4B reconstructed archive-compatible age sensitivity.

rm(list = ls())
source(file.path(
  Sys.getenv("SDY2583_REPO_ROOT", unset = "."),
  "R", "shared", "reconstructed_panel_framework.R"
))
rp_install_and_load(c("dplyr", "readr", "tibble", "purrr", "broom"))
source(file.path(
  sd_repo_root(), "R", "panels", "CP28", "MANIFEST_RECONSTRUCTED.R"
))

analysis_dir <- sd_analysis_dir("CP28")
rdata_dir <- file.path(analysis_dir, "11_RData")
out_dir <- file.path(analysis_dir, "10_age_sensitivity")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

load(file.path(
  rdata_dir,
  "SDY2583_CP28_STEP4_composite_scores_RECONSTRUCTED.RData"
))
outcomes <- CP28_MANIFEST$age_targets
stopifnot(length(outcomes) == 52L, length(unique(outcomes)) == 52L)

fit_outcome <- function(data, feature_name) {
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

  fit <- stats::lm(value ~ disease_group + age_for_model + sex, data = work)
  term <- broom::tidy(fit, conf.int = TRUE) |>
    dplyr::filter(term == "disease_groupCancer patient")
  beta <- term$estimate[1]

  tibble::tibble(
    feature = feature_name,
    n_model = stats::nobs(fit),
    n_healthy = sum(work$disease_group == "Healthy control"),
    n_cancer = sum(work$disease_group == "Cancer patient"),
    healthy_mean = mean(work$value[work$disease_group == "Healthy control"]),
    cancer_mean = mean(work$value[work$disease_group == "Cancer patient"]),
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

age_stratified_statistics <- purrr::map_dfr(
  levels(scored_data$age_group_for_model),
  function(age_group) {
    result <- purrr::map_dfr(
      outcomes,
      function(feature_name) {
        fit_outcome(
          scored_data |>
            dplyr::filter(age_group_for_model == age_group),
          feature_name
        )
      }
    )
    result$age_group_for_model <- age_group
    result |>
      dplyr::mutate(FDR_within_age_group = p.adjust(p_value, method = "BH"))
  }
)

match_data <- function(data, caliper_years) {
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
  pair_rows <- list()
  pair_number <- 0L

  for (index in seq_len(nrow(cancer))) {
    cancer_row <- cancer[index, , drop = FALSE]
    candidate_index <- which(
      !used & as.character(healthy$sex) == as.character(cancer_row$sex[1])
    )
    if (length(candidate_index) == 0L) next

    candidates <- healthy[candidate_index, , drop = FALSE] |>
      dplyr::mutate(
        candidate_index = candidate_index,
        abs_age_diff = abs(age_for_model - cancer_row$age_for_model[1])
      ) |>
      dplyr::filter(abs_age_diff <= caliper_years) |>
      dplyr::arrange(abs_age_diff, age_for_model, subject_id)
    if (nrow(candidates) == 0L) next

    healthy_row <- candidates[1, , drop = FALSE]
    selected_index <- healthy_row$candidate_index[1]
    used[selected_index] <- TRUE
    pair_number <- pair_number + 1L
    pair_id <- paste0("pair_", pair_number)
    age_difference <- healthy_row$abs_age_diff[1]

    cancer_row$pair_id <- pair_id
    cancer_row$abs_age_diff <- age_difference
    healthy_selected <- healthy[selected_index, , drop = FALSE]
    healthy_selected$pair_id <- pair_id
    healthy_selected$abs_age_diff <- age_difference

    # Archive order: cancer then healthy.
    pair_rows[[length(pair_rows) + 1L]] <- cancer_row
    pair_rows[[length(pair_rows) + 1L]] <- healthy_selected
  }

  dplyr::bind_rows(pair_rows)
}

matched_5 <- match_data(scored_data, 5)
matched_10 <- match_data(scored_data, 10)

fit_matched <- function(data, caliper_years) {
  purrr::map_dfr(outcomes, function(feature_name) {
    fit_outcome(data, feature_name)
  }) |>
    dplyr::mutate(
      caliper_years = caliper_years,
      FDR_within_caliper = p.adjust(p_value, method = "BH")
    )
}

age_matched_statistics <- dplyr::bind_rows(
  fit_matched(matched_5, 5),
  fit_matched(matched_10, 10)
)
age_matched_summary <- dplyr::bind_rows(
  tibble::tibble(
    caliper_years = 5,
    n_pairs = dplyr::n_distinct(matched_5$pair_id),
    n_rows = nrow(matched_5),
    mean_abs_age_diff = mean(matched_5$abs_age_diff),
    median_abs_age_diff = stats::median(matched_5$abs_age_diff),
    max_abs_age_diff = max(matched_5$abs_age_diff)
  ),
  tibble::tibble(
    caliper_years = 10,
    n_pairs = dplyr::n_distinct(matched_10$pair_id),
    n_rows = nrow(matched_10),
    mean_abs_age_diff = mean(matched_10$abs_age_diff),
    median_abs_age_diff = stats::median(matched_10$abs_age_diff),
    max_abs_age_diff = max(matched_10$abs_age_diff)
  )
)

run_interaction <- function(feature_name) {
  work <- scored_data |>
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
      age_z = as.numeric(scale(age_for_model))
    )

  fit <- stats::lm(value ~ disease_group * age_z + sex, data = work)
  term <- broom::tidy(fit, conf.int = TRUE) |>
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

disease_by_age_statistics <- purrr::map_dfr(outcomes, run_interaction) |>
  dplyr::mutate(FDR = p.adjust(p_value, method = "BH"))

readr::write_csv(
  age_stratified_statistics,
  file.path(out_dir, "SDY2583_CP28_age_stratified_statistics_RECONSTRUCTED.csv")
)
readr::write_csv(
  matched_5 |>
    dplyr::select(pair_id, subject_id, disease_group, sex, age_for_model, abs_age_diff),
  file.path(out_dir, "SDY2583_CP28_same_sex_age_matched_pairs_5yr_RECONSTRUCTED.csv")
)
readr::write_csv(
  matched_10 |>
    dplyr::select(pair_id, subject_id, disease_group, sex, age_for_model, abs_age_diff),
  file.path(out_dir, "SDY2583_CP28_same_sex_age_matched_pairs_10yr_RECONSTRUCTED.csv")
)
readr::write_csv(
  age_matched_statistics,
  file.path(out_dir, "SDY2583_CP28_age_matched_sensitivity_statistics_RECONSTRUCTED.csv")
)
readr::write_csv(
  age_matched_summary,
  file.path(out_dir, "SDY2583_CP28_age_matched_sensitivity_summary_RECONSTRUCTED.csv")
)
readr::write_csv(
  disease_by_age_statistics,
  file.path(out_dir, "SDY2583_CP28_disease_by_age_interaction_statistics_RECONSTRUCTED.csv")
)

save(
  outcomes,
  age_stratified_statistics,
  matched_5,
  matched_10,
  age_matched_statistics,
  age_matched_summary,
  disease_by_age_statistics,
  file = file.path(
    rdata_dir,
    "SDY2583_CP28_STEP4B_age_sensitivity_RECONSTRUCTED.RData"
  )
)
print(age_matched_summary)
