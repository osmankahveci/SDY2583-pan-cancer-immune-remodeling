# ============================================================
# SDY2583 CP8
# STEP 4B RECONSTRUCTED: age-stratified, same-sex age-matched,
# and disease-by-age interaction analyses.
#
# Archive-specific implementation confirmed by independent recalculation:
#   - exact 31-outcome family;
#   - BH FDR separately within each age stratum and caliper;
#   - cancer and healthy records sorted by sex, ascending age, subject ID;
#   - nearest-control ties resolved by the first control in that order
#     (therefore lower control age is preferred before subject ID);
#   - matched inference remains age/sex adjusted:
#       outcome ~ disease_group + age_for_model + sex.
# ============================================================

rm(list = ls())
source(file.path(
  Sys.getenv("SDY2583_REPO_ROOT", unset = "."),
  "R", "shared", "reconstructed_panel_framework.R"
))
rp_install_and_load(c("dplyr", "readr", "tibble", "purrr"))
source(file.path(sd_repo_root(), "R", "panels", "CP8", "MANIFEST_RECONSTRUCTED.R"))

analysis_dir <- sd_analysis_dir("CP8")
rdata_dir <- file.path(analysis_dir, "11_RData")
out_dir <- file.path(analysis_dir, "10_age_sensitivity")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

step4_file <- file.path(
  rdata_dir,
  "SDY2583_CP8_STEP4_composite_scores_RECONSTRUCTED.RData"
)
if (!file.exists(step4_file)) stop("Run reconstructed CP8 Step 4 first: ", step4_file)
load(step4_file)
if (!exists("scored_data")) stop("CP8 Step 4 RData does not contain scored_data.")

outcomes <- CP8_MANIFEST$age_sensitivity_outcomes
missing_outcomes <- setdiff(outcomes, names(scored_data))
if (length(missing_outcomes) > 0L) {
  stop("CP8 age-sensitivity outcomes are missing: ", paste(missing_outcomes, collapse = ", "))
}

# The reconstructed metadata framework uses portable factor labels with
# underscores; archived CP8 tables use display labels with spaces/symbols.
analysis_source <- scored_data |>
  dplyr::mutate(
    age_group_archive = dplyr::case_when(
      as.character(age_group_for_model) %in% c("Young <40", "Young_<40") ~ "Young <40",
      as.character(age_group_for_model) %in% c("Middle 40-59", "Middle_40_59") ~ "Middle 40-59",
      as.character(age_group_for_model) %in% c("Older 60+", "Older_60plus") ~ "Older 60+",
      TRUE ~ NA_character_
    )
  )

fit_adjusted <- function(data, feature) {
  d <- data |>
    dplyr::transmute(
      y = suppressWarnings(as.numeric(.data[[feature]])),
      disease_group = factor(
        as.character(disease_group),
        levels = c("Healthy control", "Cancer patient")
      ),
      age_for_model = suppressWarnings(as.numeric(age_for_model)),
      sex = factor(as.character(sex))
    ) |>
    dplyr::filter(stats::complete.cases(.))

  if (nrow(d) < 10L || dplyr::n_distinct(d$disease_group) < 2L) {
    return(tibble::tibble(
      feature = feature,
      n_model = nrow(d),
      n_healthy = sum(d$disease_group == "Healthy control"),
      n_cancer = sum(d$disease_group == "Cancer patient"),
      healthy_mean = NA_real_,
      cancer_mean = NA_real_,
      beta_cancer_vs_healthy = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      p_value = NA_real_,
      direction = NA_character_
    ))
  }

  fit <- stats::lm(y ~ disease_group + age_for_model + sex, data = d)
  coefficient_table <- summary(fit)$coefficients
  coefficient_name <- grep("^disease_group", rownames(coefficient_table), value = TRUE)[1]
  ci <- stats::confint(fit, parm = coefficient_name)
  healthy <- d$y[d$disease_group == "Healthy control"]
  cancer <- d$y[d$disease_group == "Cancer patient"]
  beta <- unname(coefficient_table[coefficient_name, "Estimate"])

  tibble::tibble(
    feature = feature,
    n_model = stats::nobs(fit),
    n_healthy = length(healthy),
    n_cancer = length(cancer),
    healthy_mean = mean(healthy),
    cancer_mean = mean(cancer),
    beta_cancer_vs_healthy = beta,
    ci_low = unname(ci[1]),
    ci_high = unname(ci[2]),
    p_value = unname(coefficient_table[coefficient_name, "Pr(>|t|)"]),
    direction = ifelse(
      beta > 0, "higher_in_cancer",
      ifelse(beta < 0, "lower_in_cancer", "no_difference")
    )
  )
}

age_levels <- c("Young <40", "Middle 40-59", "Older 60+")
age_stratified_statistics <- purrr::map_dfr(age_levels, function(age_level) {
  subset <- analysis_source |>
    dplyr::filter(age_group_archive == age_level)

  purrr::map_dfr(outcomes, ~fit_adjusted(subset, .x)) |>
    dplyr::mutate(
      age_group_for_model = age_level,
      FDR_within_age_group = stats::p.adjust(p_value, method = "BH")
    ) |>
    dplyr::select(
      age_group_for_model, feature, n_model, n_healthy, n_cancer,
      healthy_mean, cancer_mean, beta_cancer_vs_healthy,
      ci_low, ci_high, p_value, FDR_within_age_group, direction
    ) |>
    dplyr::arrange(p_value)
})

cp8_greedy_match <- function(data, caliper_years) {
  eligible <- data |>
    dplyr::filter(
      !is.na(age_for_model),
      as.character(sex) %in% c("Female", "Male"),
      !is.na(disease_group)
    )

  cancer <- eligible |>
    dplyr::filter(disease_group == "Cancer patient") |>
    dplyr::arrange(sex, age_for_model, subject_id)
  healthy <- eligible |>
    dplyr::filter(disease_group == "Healthy control") |>
    dplyr::arrange(sex, age_for_model, subject_id)

  used_healthy <- rep(FALSE, nrow(healthy))
  pair_rows <- vector("list", 0L)

  for (i in seq_len(nrow(cancer))) {
    candidates <- which(
      !used_healthy &
        as.character(healthy$sex) == as.character(cancer$sex[i])
    )
    if (length(candidates) == 0L) next

    distances <- abs(
      healthy$age_for_model[candidates] - cancer$age_for_model[i]
    )

    # which.min preserves the first tied candidate. Since healthy controls are
    # sorted by ascending age and then subject ID, this exactly reproduces the
    # archived CP8 lower-age-first tie rule.
    nearest_local <- which.min(distances)
    healthy_index <- candidates[nearest_local]

    if (distances[nearest_local] <= caliper_years) {
      used_healthy[healthy_index] <- TRUE
      pair_rows[[length(pair_rows) + 1L]] <- tibble::tibble(
        pair_number = length(pair_rows) + 1L,
        cancer_index = i,
        healthy_index = healthy_index,
        abs_age_diff = distances[nearest_local]
      )
    }
  }

  pair_table <- dplyr::bind_rows(pair_rows)
  if (nrow(pair_table) == 0L) {
    return(list(data = eligible[0, ], pair_listing = tibble::tibble()))
  }

  matched_data <- purrr::pmap_dfr(
    pair_table,
    function(pair_number, cancer_index, healthy_index, abs_age_diff) {
      pair_id <- paste0("pair_", pair_number)
      dplyr::bind_rows(
        cancer[cancer_index, ] |>
          dplyr::mutate(pair_id = pair_id, abs_age_diff = abs_age_diff),
        healthy[healthy_index, ] |>
          dplyr::mutate(pair_id = pair_id, abs_age_diff = abs_age_diff)
      )
    }
  )

  pair_listing <- matched_data |>
    dplyr::select(
      pair_id, subject_id, disease_group, sex, age_for_model, abs_age_diff
    )

  list(data = matched_data, pair_listing = pair_listing)
}

calipers <- c(5, 10)
matched_objects <- setNames(
  lapply(calipers, function(caliper) cp8_greedy_match(analysis_source, caliper)),
  as.character(calipers)
)

age_matched_sensitivity_statistics <- purrr::imap_dfr(
  matched_objects,
  function(match_object, caliper_label) {
    purrr::map_dfr(outcomes, ~fit_adjusted(match_object$data, .x)) |>
      dplyr::mutate(
        caliper_years = as.numeric(caliper_label),
        FDR_within_caliper = stats::p.adjust(p_value, method = "BH")
      ) |>
      dplyr::select(
        caliper_years, feature, n_model, n_healthy, n_cancer,
        healthy_mean, cancer_mean, beta_cancer_vs_healthy,
        ci_low, ci_high, p_value, FDR_within_caliper, direction
      ) |>
      dplyr::arrange(p_value)
  }
)

age_matched_sensitivity_summary <- purrr::imap_dfr(
  matched_objects,
  function(match_object, caliper_label) {
    pair_listing <- match_object$pair_listing
    one_per_pair <- pair_listing |>
      dplyr::distinct(pair_id, abs_age_diff)

    tibble::tibble(
      caliper_years = as.numeric(caliper_label),
      n_pairs = dplyr::n_distinct(pair_listing$pair_id),
      n_rows = nrow(pair_listing),
      mean_abs_age_diff = mean(one_per_pair$abs_age_diff),
      median_abs_age_diff = stats::median(one_per_pair$abs_age_diff),
      max_abs_age_diff = max(one_per_pair$abs_age_diff)
    )
  }
)

interaction_source <- analysis_source |>
  dplyr::mutate(age_z = as.numeric(scale(age_for_model)))

fit_interaction <- function(feature) {
  d <- interaction_source |>
    dplyr::transmute(
      y = suppressWarnings(as.numeric(.data[[feature]])),
      disease_group = factor(
        as.character(disease_group),
        levels = c("Healthy control", "Cancer patient")
      ),
      age_z = age_z,
      sex = factor(as.character(sex))
    ) |>
    dplyr::filter(stats::complete.cases(.))

  fit <- stats::lm(y ~ disease_group * age_z + sex, data = d)
  coefficient_table <- summary(fit)$coefficients
  coefficient_name <- grep(
    "disease_group.*:age_z",
    rownames(coefficient_table),
    value = TRUE
  )[1]
  ci <- stats::confint(fit, parm = coefficient_name)
  beta <- unname(coefficient_table[coefficient_name, "Estimate"])

  tibble::tibble(
    feature = feature,
    n_model = stats::nobs(fit),
    beta_interaction = beta,
    ci_low = unname(ci[1]),
    ci_high = unname(ci[2]),
    p_value = unname(coefficient_table[coefficient_name, "Pr(>|t|)"]),
    direction = ifelse(
      beta > 0,
      "stronger_cancer_effect_with_older_age",
      "weaker_cancer_effect_with_older_age"
    )
  )
}

disease_by_age_interaction_statistics <- purrr::map_dfr(
  outcomes,
  fit_interaction
) |>
  dplyr::mutate(FDR = stats::p.adjust(p_value, method = "BH")) |>
  dplyr::arrange(p_value)

readr::write_csv(
  age_stratified_statistics,
  file.path(out_dir, "SDY2583_CP8_age_stratified_statistics_RECONSTRUCTED.csv")
)
readr::write_csv(
  age_matched_sensitivity_statistics,
  file.path(out_dir, "SDY2583_CP8_age_matched_sensitivity_statistics_RECONSTRUCTED.csv")
)
readr::write_csv(
  age_matched_sensitivity_summary,
  file.path(out_dir, "SDY2583_CP8_age_matched_sensitivity_summary_RECONSTRUCTED.csv")
)
readr::write_csv(
  disease_by_age_interaction_statistics,
  file.path(out_dir, "SDY2583_CP8_disease_by_age_interaction_statistics_RECONSTRUCTED.csv")
)

for (caliper_label in names(matched_objects)) {
  readr::write_csv(
    matched_objects[[caliper_label]]$pair_listing,
    file.path(
      out_dir,
      paste0(
        "SDY2583_CP8_same_sex_age_matched_pairs_",
        caliper_label,
        "yr_RECONSTRUCTED.csv"
      )
    )
  )
}

save(
  outcomes,
  age_stratified_statistics,
  matched_objects,
  age_matched_sensitivity_statistics,
  age_matched_sensitivity_summary,
  disease_by_age_interaction_statistics,
  file = file.path(
    rdata_dir,
    "SDY2583_CP8_STEP4B_age_sensitivity_RECONSTRUCTED.RData"
  )
)

print(age_matched_sensitivity_summary)
message(
  "Archived CP8 benchmarks: 31 outcomes; 265 five-year pairs; ",
  "273 ten-year pairs."
)
