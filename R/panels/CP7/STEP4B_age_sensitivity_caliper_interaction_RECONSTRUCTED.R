# ============================================================
# SDY2583 CP7
# STEP 4B RECONSTRUCTED: age-stratified, same-sex age-matched,
# and disease-by-age interaction analyses.
#
# Reconstruction basis:
#   - archived 32-outcome age-sensitivity feature map;
#   - archived 5-year and 10-year matched subject datasets;
#   - archived model formulas and numerical results.
#
# CP7-specific matching details required for exact reproduction:
#   1. retain only Female/Male records with valid age;
#   2. assign row IDs in the retained source order;
#   3. process cancer subjects from oldest to youngest;
#   4. greedily select the nearest unused same-sex healthy control;
#   5. fit matched models as outcome ~ disease_group + matched_pair_id.
# ============================================================

rm(list = ls())
source(file.path(
  Sys.getenv("SDY2583_REPO_ROOT", unset = "."),
  "R", "shared", "reconstructed_panel_framework.R"
))
rp_install_and_load(c("dplyr", "readr", "tibble", "purrr", "tidyr"))
source(file.path(sd_repo_root(), "R", "panels", "CP7", "MANIFEST_RECONSTRUCTED.R"))

analysis_dir <- sd_analysis_dir("CP7")
rdata_dir <- file.path(analysis_dir, "11_RData")
out_dir <- file.path(analysis_dir, "10_age_sensitivity")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

step4_file <- file.path(
  rdata_dir,
  "SDY2583_CP7_STEP4_composite_scores_RECONSTRUCTED.RData"
)
if (!file.exists(step4_file)) stop("Run reconstructed CP7 Step 4 first: ", step4_file)
load(step4_file)
if (!exists("scored_data")) stop("Step 4 RData does not contain scored_data.")

feature_map <- tibble::tibble(
  feature = names(CP7_MANIFEST$age_sensitivity_modules),
  module = unname(CP7_MANIFEST$age_sensitivity_modules)
)
outcomes <- feature_map$feature
missing_outcomes <- setdiff(outcomes, names(scored_data))
if (length(missing_outcomes) > 0L) {
  stop("CP7 age-sensitivity outcomes are missing: ", paste(missing_outcomes, collapse = ", "))
}

fit_adjusted <- function(data, feature, subset_label, module) {
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

  fit <- stats::lm(y ~ disease_group + age_for_model + sex, data = d)
  coefficient_table <- summary(fit)$coefficients
  coefficient_name <- grep("^disease_group", rownames(coefficient_table), value = TRUE)[1]
  ci <- stats::confint(fit, parm = coefficient_name)
  healthy <- d$y[d$disease_group == "Healthy control"]
  cancer <- d$y[d$disease_group == "Cancer patient"]
  beta <- unname(coefficient_table[coefficient_name, "Estimate"])

  tibble::tibble(
    subset_label = subset_label,
    module = module,
    feature = feature,
    n_model = stats::nobs(fit),
    n_healthy = length(healthy),
    n_cancer = length(cancer),
    healthy_mean = mean(healthy),
    cancer_mean = mean(cancer),
    healthy_median = stats::median(healthy),
    cancer_median = stats::median(cancer),
    beta_cancer_vs_healthy = beta,
    ci_low = unname(ci[1]),
    ci_high = unname(ci[2]),
    t_value = unname(coefficient_table[coefficient_name, "t value"]),
    p_value = unname(coefficient_table[coefficient_name, "Pr(>|t|)"]),
    model_formula = "y ~ disease_group + age_for_model + sex",
    direction = ifelse(
      beta > 0, "higher_in_cancer",
      ifelse(beta < 0, "lower_in_cancer", "no_difference")
    )
  )
}

age_levels <- levels(droplevels(scored_data$age_group_for_model))
age_stratified_statistics <- purrr::map_dfr(age_levels, function(age_level) {
  subset <- scored_data |>
    dplyr::filter(as.character(age_group_for_model) == age_level)
  purrr::map2_dfr(
    feature_map$feature,
    feature_map$module,
    ~fit_adjusted(
      subset,
      feature = .x,
      subset_label = paste0("age_stratum_", age_level),
      module = .y
    )
  )
}) |>
  dplyr::group_by(subset_label) |>
  dplyr::mutate(fdr_within_stratum = stats::p.adjust(p_value, method = "BH")) |>
  dplyr::ungroup() |>
  dplyr::arrange(subset_label, p_value)

cp7_greedy_match <- function(data, caliper_years) {
  eligible <- data |>
    dplyr::filter(
      !is.na(age_for_model),
      as.character(sex) %in% c("Female", "Male"),
      !is.na(disease_group)
    ) |>
    dplyr::mutate(.row_id_for_matching = dplyr::row_number())

  cancer <- eligible |>
    dplyr::filter(disease_group == "Cancer patient") |>
    dplyr::arrange(dplyr::desc(age_for_model), .row_id_for_matching)
  healthy <- eligible |>
    dplyr::filter(disease_group == "Healthy control")

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
    nearest_local <- order(
      distances,
      healthy$.row_id_for_matching[candidates]
    )[1]
    healthy_index <- candidates[nearest_local]

    if (distances[nearest_local] <= caliper_years) {
      used_healthy[healthy_index] <- TRUE
      pair_rows[[length(pair_rows) + 1L]] <- tibble::tibble(
        pair_number = length(pair_rows) + 1L,
        cancer_index = i,
        healthy_index = healthy_index,
        matched_age_diff_abs = distances[nearest_local]
      )
    }
  }

  pair_table <- dplyr::bind_rows(pair_rows)
  if (nrow(pair_table) == 0L) {
    return(list(data = eligible[0, ], pairs = pair_table))
  }

  matched_data <- purrr::pmap_dfr(pair_table, function(
    pair_number, cancer_index, healthy_index, matched_age_diff_abs
  ) {
    pair_id <- paste0("pair_", pair_number)
    dplyr::bind_rows(
      cancer[cancer_index, ] |>
        dplyr::mutate(
          matched_pair_id = pair_id,
          matched_age_diff_abs = matched_age_diff_abs
        ),
      healthy[healthy_index, ] |>
        dplyr::mutate(
          matched_pair_id = pair_id,
          matched_age_diff_abs = matched_age_diff_abs
        )
    )
  })

  list(data = matched_data, pairs = pair_table)
}

fit_paired <- function(matched_data, feature, subset_label, module) {
  pair_wide <- matched_data |>
    dplyr::transmute(
      matched_pair_id,
      disease_group = as.character(disease_group),
      value = suppressWarnings(as.numeric(.data[[feature]]))
    ) |>
    dplyr::filter(!is.na(value)) |>
    tidyr::pivot_wider(names_from = disease_group, values_from = value) |>
    dplyr::filter(!is.na(`Healthy control`), !is.na(`Cancer patient`))

  differences <- pair_wide$`Cancer patient` - pair_wide$`Healthy control`
  n_pairs <- length(differences)
  beta <- mean(differences)
  standard_error <- stats::sd(differences) / sqrt(n_pairs)
  t_value <- beta / standard_error
  degrees_freedom <- n_pairs - 1L
  p_value <- 2 * stats::pt(abs(t_value), df = degrees_freedom, lower.tail = FALSE)
  critical_value <- stats::qt(0.975, df = degrees_freedom)

  tibble::tibble(
    subset_label = subset_label,
    module = module,
    feature = feature,
    n_model = 2L * n_pairs,
    n_pairs = n_pairs,
    n_healthy = n_pairs,
    n_cancer = n_pairs,
    healthy_mean = mean(pair_wide$`Healthy control`),
    cancer_mean = mean(pair_wide$`Cancer patient`),
    beta_cancer_vs_healthy = beta,
    ci_low = beta - critical_value * standard_error,
    ci_high = beta + critical_value * standard_error,
    t_value = t_value,
    p_value = p_value,
    model_formula = "y ~ disease_group + matched_pair_id",
    direction = ifelse(
      beta > 0, "higher_in_cancer",
      ifelse(beta < 0, "lower_in_cancer", "no_difference")
    )
  )
}

calipers <- c(5, 10)
matched_objects <- setNames(lapply(calipers, function(caliper) {
  cp7_greedy_match(scored_data, caliper)
}), as.character(calipers))

age_matched_sensitivity_statistics <- purrr::imap_dfr(
  matched_objects,
  function(match_object, caliper_label) {
    subset_label <- paste0(
      "same_sex_nearest_age_caliper_", caliper_label, "y"
    )
    purrr::map2_dfr(
      feature_map$feature,
      feature_map$module,
      ~fit_paired(
        match_object$data,
        feature = .x,
        subset_label = subset_label,
        module = .y
      )
    ) |>
      dplyr::mutate(
        fdr_within_subset = stats::p.adjust(p_value, method = "BH")
      ) |>
      dplyr::arrange(p_value)
  }
)

age_matched_sensitivity_summary <- purrr::imap_dfr(
  matched_objects,
  function(match_object, caliper_label) {
    tibble::tibble(
      subset_label = paste0(
        "same_sex_nearest_age_caliper_", caliper_label, "y"
      ),
      caliper_years = as.numeric(caliper_label),
      n_rows = nrow(match_object$data),
      n_pairs = nrow(match_object$pairs),
      n_healthy = sum(
        match_object$data$disease_group == "Healthy control",
        na.rm = TRUE
      ),
      n_cancer = sum(
        match_object$data$disease_group == "Cancer patient",
        na.rm = TRUE
      ),
      mean_abs_age_difference = mean(
        match_object$pairs$matched_age_diff_abs,
        na.rm = TRUE
      ),
      median_abs_age_difference = stats::median(
        match_object$pairs$matched_age_diff_abs,
        na.rm = TRUE
      ),
      max_abs_age_difference = max(
        match_object$pairs$matched_age_diff_abs,
        na.rm = TRUE
      )
    )
  }
)

interaction_source <- scored_data |>
  dplyr::mutate(age_z = as.numeric(scale(age_for_model)))

fit_interaction <- function(feature, module) {
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

  tibble::tibble(
    module = module,
    feature = feature,
    n_model = stats::nobs(fit),
    beta_disease_age_interaction = unname(
      coefficient_table[coefficient_name, "Estimate"]
    ),
    ci_low = unname(ci[1]),
    ci_high = unname(ci[2]),
    t_value = unname(coefficient_table[coefficient_name, "t value"]),
    p_value = unname(coefficient_table[coefficient_name, "Pr(>|t|)"]),
    model_formula = "y ~ disease_group * age_z + sex",
    interpretation = paste(
      "Positive beta means the cancer-vs-healthy difference",
      "increases with older age."
    )
  )
}

disease_by_age_interaction_statistics <- purrr::map2_dfr(
  feature_map$feature,
  feature_map$module,
  fit_interaction
) |>
  dplyr::mutate(fdr = stats::p.adjust(p_value, method = "BH")) |>
  dplyr::arrange(p_value)

readr::write_csv(
  feature_map,
  file.path(out_dir, "SDY2583_CP7_age_sensitivity_feature_map_RECONSTRUCTED.csv")
)
readr::write_csv(
  age_stratified_statistics,
  file.path(out_dir, "SDY2583_CP7_age_stratified_statistics_RECONSTRUCTED.csv")
)
readr::write_csv(
  age_matched_sensitivity_statistics,
  file.path(out_dir, "SDY2583_CP7_age_matched_sensitivity_statistics_RECONSTRUCTED.csv")
)
readr::write_csv(
  age_matched_sensitivity_summary,
  file.path(out_dir, "SDY2583_CP7_age_matched_sensitivity_summary_RECONSTRUCTED.csv")
)
readr::write_csv(
  disease_by_age_interaction_statistics,
  file.path(out_dir, "SDY2583_CP7_disease_by_age_interaction_statistics_RECONSTRUCTED.csv")
)

for (caliper_label in names(matched_objects)) {
  readr::write_csv(
    matched_objects[[caliper_label]]$data,
    file.path(
      out_dir,
      paste0(
        "SDY2583_CP7_age_matched_dataset_caliper_",
        caliper_label,
        "y_RECONSTRUCTED.csv"
      )
    )
  )
}

save(
  feature_map, outcomes, age_stratified_statistics,
  age_matched_sensitivity_statistics, age_matched_sensitivity_summary,
  disease_by_age_interaction_statistics, matched_objects,
  file = file.path(
    rdata_dir,
    "SDY2583_CP7_STEP4B_age_sensitivity_RECONSTRUCTED.RData"
  )
)

print(age_matched_sensitivity_summary)
message(
  "Archived CP7 matching benchmarks: 5y = ",
  CP7_MANIFEST$benchmarks$matched_pairs_5y,
  " pairs; 10y = ", CP7_MANIFEST$benchmarks$matched_pairs_10y,
  " pairs."
)
