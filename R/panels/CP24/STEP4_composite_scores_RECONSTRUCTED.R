# CP24 full-cohort Step 4 reconstructed official hierarchical composite scores.
#
# Official hierarchy:
# - five domain scores are means of oriented raw-feature z scores;
# - core differentiation is the arithmetic mean of the first four domain scores;
# - integrated remodeling is the arithmetic mean of core differentiation and
#   the PD-1-associated terminal score.

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
out_dir <- file.path(analysis_dir, "08_composite_scores_RECONSTRUCTED")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

load(file.path(
  rdata_dir,
  "SDY2583_CP24_STEP3A_metadata_merge_age_QC_RECONSTRUCTED.RData"
))
load(file.path(
  rdata_dir,
  "SDY2583_CP24_STEP3B_statistics_RECONSTRUCTED.RData"
))

if (!exists("cp24_run_outcome")) {
  stop("CP24 Step 3B model helper was not found.")
}

scored_data <- cp24_build_official_scores(analysis_data)
score_names <- CP24_MANIFEST$score_names

cp24_run_score_family <- function(data, analysis_set) {
  result <- purrr::map_dfr(score_names, function(score_name) {
    cp24_run_outcome(
      data,
      module_name = "composite_scores",
      outcome_name = score_name,
      analysis_set = analysis_set
    )
  }) |>
    dplyr::rename(score = outcome) |>
    dplyr::select(-module) |>
    dplyr::mutate(p_FDR_scores = p.adjust(p_value_model, method = "BH"))

  result
}

score_statistics_with_set <- cp24_run_score_family(
  scored_data,
  "all_850_clean_age"
)
score_statistics <- score_statistics_with_set |>
  dplyr::select(-analysis_set)

binary_score_statistics <- cp24_run_score_family(
  scored_data |>
    dplyr::filter(as.character(sex) %in% c("Female", "Male")) |>
    dplyr::mutate(sex = droplevels(factor(as.character(sex)))),
  "binary_sex_only_clean_age"
)

event_score_sets <- list(
  all_850_clean_age = scored_data,
  cd3cd8_ge500_clean_age = scored_data |>
    dplyr::filter(!is.na(n_cd3_cd8_pos), n_cd3_cd8_pos >= 500),
  cd3cd8_ge1000_clean_age = scored_data |>
    dplyr::filter(!is.na(n_cd3_cd8_pos), n_cd3_cd8_pos >= 1000)
)

event_score_statistics <- purrr::imap_dfr(
  event_score_sets,
  function(data, label) cp24_run_score_family(data, label)
)

feature_sets <- purrr::imap_dfr(
  CP24_RAW_SCORE_DEFINITIONS,
  function(definition, score_name) {
    tibble::tibble(
      score = score_name,
      positive_features = paste(
        names(definition)[definition > 0],
        collapse = "; "
      ),
      negative_features = paste(
        names(definition)[definition < 0],
        collapse = "; "
      ),
      score_level = "oriented_raw_feature_mean"
    )
  }
) |>
  dplyr::bind_rows(
    tibble::tibble(
      score = "CP24_core_differentiation_score",
      positive_features = paste(CP24_CORE_DOMAIN_SCORES, collapse = "; "),
      negative_features = "",
      score_level = "domain_score_mean"
    ),
    tibble::tibble(
      score = "CP24_integrated_remodeling_score",
      positive_features = paste(
        c(
          "CP24_core_differentiation_score",
          "CP24_PD1_associated_terminal_score"
        ),
        collapse = "; "
      ),
      negative_features = "",
      score_level = "top_level_score_mean"
    )
  )

readr::write_csv(
  feature_sets,
  file.path(out_dir, "SDY2583_CP24_composite_score_feature_sets_RECONSTRUCTED.csv")
)
readr::write_csv(
  scored_data,
  file.path(out_dir, "SDY2583_CP24_score_dataset_RECONSTRUCTED.csv")
)
readr::write_csv(
  score_statistics,
  file.path(
    out_dir,
    "SDY2583_CP24_composite_score_results_CLEANAGE_RECONSTRUCTED.csv"
  )
)
readr::write_csv(
  binary_score_statistics,
  file.path(
    out_dir,
    "SDY2583_CP24_composite_binary_sex_sensitivity_RECONSTRUCTED.csv"
  )
)
readr::write_csv(
  event_score_statistics,
  file.path(
    out_dir,
    "SDY2583_CP24_composite_event_sensitivity_RECONSTRUCTED.csv"
  )
)

save(
  scored_data,
  score_statistics,
  score_statistics_with_set,
  binary_score_statistics,
  event_score_statistics,
  feature_sets,
  score_names,
  cp24_run_score_family,
  file = file.path(
    rdata_dir,
    "SDY2583_CP24_STEP4_composite_scores_RECONSTRUCTED.RData"
  )
)

print(score_statistics)
