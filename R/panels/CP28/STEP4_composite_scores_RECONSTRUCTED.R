# CP28 Step 4 reconstructed archive-compatible composite scores.

rm(list = ls())
source(file.path(
  Sys.getenv("SDY2583_REPO_ROOT", unset = "."),
  "R", "shared", "reconstructed_panel_framework.R"
))
rp_install_and_load(c("dplyr", "readr", "tibble", "purrr"))
source(file.path(
  sd_repo_root(), "R", "panels", "CP28", "MANIFEST_RECONSTRUCTED.R"
))

analysis_dir <- sd_analysis_dir("CP28")
rdata_dir <- file.path(analysis_dir, "11_RData")
out_dir <- file.path(analysis_dir, "08_composite_scores")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

load(file.path(
  rdata_dir,
  "SDY2583_CP28_STEP3A_metadata_merge_age_QC_RECONSTRUCTED.RData"
))

scored_data <- rp_build_scores(
  analysis_data,
  CP28_MANIFEST$score_definitions,
  CP28_MANIFEST$integrated_score
)
score_names <- CP28_SCORE_NAMES

format_score_statistics <- function(result) {
  result |>
    dplyr::transmute(
      feature = score,
      n_model, n_healthy, n_cancer,
      healthy_mean, cancer_mean,
      beta_cancer_vs_healthy,
      ci_low, ci_high, p_value, direction,
      FDR = fdr
    )
}

score_statistics <- rp_fit_scores(scored_data, score_names) |>
  format_score_statistics()
binary_score_statistics <- rp_fit_scores(
  scored_data,
  score_names,
  binary_sex = TRUE
) |>
  format_score_statistics()

run_event_score <- function(score_name) {
  if (score_name == CP28_MANIFEST$integrated_score) {
    filtered <- scored_data |>
      dplyr::filter(
        n_cd3_pos_dump_low >= CP28_MANIFEST$event_qc,
        n_cd3_cd8_pos >= CP28_MANIFEST$event_qc,
        n_cd3neg_cd56pos_nk_like >= CP28_MANIFEST$event_qc
      )
    event_column <- NULL
  } else {
    filtered <- scored_data
    event_column <- unname(CP28_MANIFEST$score_event_qc[[score_name]])
    if (is.na(event_column) || !nzchar(event_column)) event_column <- NULL
  }

  rp_fit_one(
    filtered,
    score_name,
    subset_label = "event_QC_ge300",
    binary_sex = FALSE,
    event_col = event_column,
    min_events = ifelse(is.null(event_column), 0, CP28_MANIFEST$event_qc)
  ) |>
    dplyr::rename(feature = feature)
}

event_score_statistics <- purrr::map_dfr(score_names, run_event_score) |>
  dplyr::mutate(FDR = p.adjust(p_value, method = "BH")) |>
  dplyr::select(
    feature, n_model, n_healthy, n_cancer,
    healthy_mean, cancer_mean, beta_cancer_vs_healthy,
    ci_low, ci_high, p_value, direction, FDR
  )

feature_sets <- purrr::imap_dfr(
  CP28_MANIFEST$score_definitions,
  function(definition, score_name) {
    tibble::tibble(
      composite_score = score_name,
      feature = c(definition$positive, definition$negative),
      orientation = c(
        rep("positive", length(definition$positive)),
        rep("negative", length(definition$negative))
      )
    )
  }
)

score_direction_sensitivity <- score_statistics |>
  dplyr::select(
    feature,
    beta_main = beta_cancer_vs_healthy,
    direction_main = direction,
    FDR_main = FDR
  ) |>
  dplyr::left_join(
    binary_score_statistics |>
      dplyr::select(
        feature,
        beta_binary = beta_cancer_vs_healthy,
        direction_binary = direction,
        FDR_binary = FDR
      ),
    by = "feature"
  ) |>
  dplyr::left_join(
    event_score_statistics |>
      dplyr::select(
        feature,
        beta_event_QC = beta_cancer_vs_healthy,
        direction_event_QC = direction,
        FDR_event_QC = FDR
      ),
    by = "feature"
  ) |>
  dplyr::mutate(
    direction_preserved_binary = direction_main == direction_binary,
    direction_preserved_event_QC = direction_main == direction_event_QC,
    direction_preserved_all = direction_preserved_binary &
      direction_preserved_event_QC
  )

readr::write_csv(
  feature_sets,
  file.path(out_dir, "SDY2583_CP28_composite_score_feature_sets_RECONSTRUCTED.csv")
)
readr::write_csv(
  scored_data,
  file.path(out_dir, "SDY2583_CP28_analysis_data_with_composite_scores_RECONSTRUCTED.csv")
)
readr::write_csv(
  score_statistics,
  file.path(out_dir, "SDY2583_CP28_composite_score_statistics_RECONSTRUCTED.csv")
)
readr::write_csv(
  binary_score_statistics,
  file.path(
    out_dir,
    "SDY2583_CP28_composite_score_binary_sex_sensitivity_statistics_RECONSTRUCTED.csv"
  )
)
readr::write_csv(
  event_score_statistics,
  file.path(
    out_dir,
    "SDY2583_CP28_composite_score_event_QC_sensitivity_statistics_RECONSTRUCTED.csv"
  )
)
readr::write_csv(
  score_direction_sensitivity,
  file.path(out_dir, "SDY2583_CP28_composite_score_direction_sensitivity_RECONSTRUCTED.csv")
)

save(
  scored_data,
  score_statistics,
  binary_score_statistics,
  event_score_statistics,
  feature_sets,
  score_direction_sensitivity,
  score_names,
  file = file.path(
    rdata_dir,
    "SDY2583_CP28_STEP4_composite_scores_RECONSTRUCTED.RData"
  )
)
print(score_statistics)
