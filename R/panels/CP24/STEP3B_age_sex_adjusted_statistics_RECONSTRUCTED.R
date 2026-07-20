# CP24 full-cohort Step 3B reconstructed official clean-age statistics.
#
# Archive-specific convention:
# - descriptive n/summary values use all subjects in each event-count set;
# - n_model reports subjects with a valid clean age and outcome;
# - lm() additionally performs complete-case removal for sex. Thus the official
#   all-850 table reports n_model = 832 while the fitted model has 828 records.

rm(list = ls())
source(file.path(
  Sys.getenv("SDY2583_REPO_ROOT", unset = "."),
  "R", "shared", "reconstructed_panel_framework.R"
))
rp_install_and_load(c(
  "dplyr", "readr", "tibble", "purrr", "broom", "tidyr"
))
source(file.path(
  sd_repo_root(), "R", "panels", "CP24", "MANIFEST_RECONSTRUCTED.R"
))

analysis_dir <- sd_analysis_dir("CP24")
rdata_dir <- file.path(analysis_dir, "11_RData")
out_dir <- file.path(analysis_dir, "09_statistics_RECONSTRUCTED")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

load(file.path(
  rdata_dir,
  "SDY2583_CP24_STEP3A_metadata_merge_age_QC_RECONSTRUCTED.RData"
))

module_map <- CP24_MANIFEST$module_map

cp24_run_outcome <- function(data, module_name, outcome_name, analysis_set) {
  if (!(outcome_name %in% names(data))) {
    stop("Missing CP24 outcome: ", outcome_name)
  }

  descriptive <- data |>
    dplyr::transmute(
      value = suppressWarnings(as.numeric(.data[[outcome_name]])),
      disease_group = factor(
        as.character(disease_group),
        levels = c("Healthy control", "Cancer patient")
      )
    ) |>
    dplyr::filter(!is.na(value), !is.na(disease_group))

  model_count_data <- data |>
    dplyr::transmute(
      value = suppressWarnings(as.numeric(.data[[outcome_name]])),
      disease_group = disease_group,
      age_for_model = suppressWarnings(as.numeric(age_for_model))
    ) |>
    dplyr::filter(
      !is.na(value),
      !is.na(disease_group),
      !is.na(age_for_model)
    )

  model_data <- data |>
    dplyr::transmute(
      value = suppressWarnings(as.numeric(.data[[outcome_name]])),
      disease_group = factor(
        as.character(disease_group),
        levels = c("Healthy control", "Cancer patient")
      ),
      age_for_model = suppressWarnings(as.numeric(age_for_model)),
      sex = factor(as.character(sex))
    ) |>
    dplyr::filter(stats::complete.cases(.)) |>
    dplyr::mutate(sex = droplevels(sex))

  if (nrow(model_data) < 50L || dplyr::n_distinct(model_data$disease_group) < 2L) {
    stop("Insufficient CP24 model data for ", outcome_name, " in ", analysis_set)
  }

  fit <- stats::lm(
    value ~ disease_group + age_for_model + sex,
    data = model_data
  )
  tidy_fit <- broom::tidy(fit, conf.int = TRUE)
  disease_row <- tidy_fit |>
    dplyr::filter(term == "disease_groupCancer patient")
  if (nrow(disease_row) != 1L) {
    stop("Cancer-vs-healthy coefficient missing for ", outcome_name)
  }

  healthy <- descriptive$value[
    descriptive$disease_group == "Healthy control"
  ]
  cancer <- descriptive$value[
    descriptive$disease_group == "Cancer patient"
  ]

  wilcox_p <- tryCatch(
    stats::wilcox.test(
      value ~ disease_group,
      data = descriptive,
      exact = FALSE
    )$p.value,
    error = function(error) NA_real_
  )

  tibble::tibble(
    analysis_set = analysis_set,
    module = module_name,
    outcome = outcome_name,
    n_healthy = length(healthy),
    healthy_mean = mean(healthy, na.rm = TRUE),
    healthy_sd = stats::sd(healthy, na.rm = TRUE),
    healthy_median = stats::median(healthy, na.rm = TRUE),
    healthy_q1 = as.numeric(stats::quantile(healthy, 0.25, na.rm = TRUE)),
    healthy_q3 = as.numeric(stats::quantile(healthy, 0.75, na.rm = TRUE)),
    n_cancer = length(cancer),
    cancer_mean = mean(cancer, na.rm = TRUE),
    cancer_sd = stats::sd(cancer, na.rm = TRUE),
    cancer_median = stats::median(cancer, na.rm = TRUE),
    cancer_q1 = as.numeric(stats::quantile(cancer, 0.25, na.rm = TRUE)),
    cancer_q3 = as.numeric(stats::quantile(cancer, 0.75, na.rm = TRUE)),
    mean_difference_cancer_minus_healthy = mean(cancer, na.rm = TRUE) -
      mean(healthy, na.rm = TRUE),
    median_difference_cancer_minus_healthy = stats::median(cancer, na.rm = TRUE) -
      stats::median(healthy, na.rm = TRUE),
    # This reproduces the archived reporting convention, not stats::nobs(fit).
    n_model = nrow(model_count_data),
    covariates = "age_years_clean + sex",
    beta_cancer = disease_row$estimate[1],
    std_error = disease_row$std.error[1],
    t_value = disease_row$statistic[1],
    p_value_model = disease_row$p.value[1],
    ci_lower = disease_row$conf.low[1],
    ci_upper = disease_row$conf.high[1],
    p_value_wilcox = wilcox_p
  )
}

cp24_run_model_family <- function(data, analysis_set) {
  result <- purrr::pmap_dfr(
    module_map,
    function(module, feature) {
      cp24_run_outcome(data, module, feature, analysis_set)
    }
  )

  result |>
    dplyr::group_by(module) |>
    dplyr::mutate(p_FDR_module = p.adjust(p_value_model, method = "BH")) |>
    dplyr::ungroup() |>
    dplyr::mutate(p_FDR_all = p.adjust(p_value_model, method = "BH"))
}

main_statistics <- cp24_run_model_family(
  analysis_data,
  "all_850_clean_age"
)

binary_sex_data <- analysis_data |>
  dplyr::filter(as.character(sex) %in% c("Female", "Male")) |>
  dplyr::mutate(sex = droplevels(factor(as.character(sex))))
binary_sex_statistics <- cp24_run_model_family(
  binary_sex_data,
  "binary_sex_only_clean_age"
)

event_sets <- list(
  all_850_clean_age = analysis_data,
  cd3cd8_ge500_clean_age = analysis_data |>
    dplyr::filter(!is.na(n_cd3_cd8_pos), n_cd3_cd8_pos >= 500),
  cd3cd8_ge1000_clean_age = analysis_data |>
    dplyr::filter(!is.na(n_cd3_cd8_pos), n_cd3_cd8_pos >= 1000)
)

event_statistics <- purrr::imap_dfr(
  event_sets,
  function(data, label) cp24_run_model_family(data, label)
)

direction_summary <- event_statistics |>
  dplyr::mutate(direction = dplyr::case_when(
    beta_cancer > 0 ~ "higher_in_cancer",
    beta_cancer < 0 ~ "lower_in_cancer",
    TRUE ~ "no_direction"
  )) |>
  dplyr::group_by(module, outcome) |>
  dplyr::summarise(
    direction_preserved = dplyr::n_distinct(direction, na.rm = TRUE) == 1L,
    module_fdr_all_sets = all(p_FDR_module < 0.05, na.rm = TRUE),
    global_fdr_all_sets = all(p_FDR_all < 0.05, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(
  module_map,
  file.path(out_dir, "SDY2583_CP24_module_feature_map_RECONSTRUCTED.csv")
)
readr::write_csv(
  main_statistics,
  file.path(out_dir, "SDY2583_CP24_CLEANAGE_model_stats_RECONSTRUCTED.csv")
)
readr::write_csv(
  binary_sex_statistics,
  file.path(out_dir, "SDY2583_CP24_binary_sex_sensitivity_RECONSTRUCTED.csv")
)
readr::write_csv(
  event_statistics,
  file.path(
    out_dir,
    "SDY2583_CP24_event_count_sensitivity_all_outcomes_RECONSTRUCTED.csv"
  )
)
readr::write_csv(
  direction_summary,
  file.path(
    out_dir,
    "SDY2583_CP24_event_count_sensitivity_direction_summary_RECONSTRUCTED.csv"
  )
)

save(
  main_statistics,
  binary_sex_statistics,
  event_statistics,
  direction_summary,
  module_map,
  cp24_run_outcome,
  cp24_run_model_family,
  file = file.path(
    rdata_dir,
    "SDY2583_CP24_STEP3B_statistics_RECONSTRUCTED.RData"
  )
)

cat("CP24 official modeled rows:", nrow(main_statistics), "\n")
cat("CP24 unique modeled outcomes:", dplyr::n_distinct(main_statistics$outcome), "\n")
