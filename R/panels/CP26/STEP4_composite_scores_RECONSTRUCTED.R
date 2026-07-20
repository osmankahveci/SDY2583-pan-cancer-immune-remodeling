# CP26 Step 4 reconstructed archive-compatible composite scores.
# Main scores are standardized within the 832 model-ready subjects.

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
out_dir <- file.path(analysis_dir, "06_composite_scores")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

load(file.path(
  rdata_dir,
  "SDY2583_CP26_STEP3A_metadata_merge_age_QC_RECONSTRUCTED.RData"
))

model_ready_data <- analysis_data |>
  dplyr::filter(
    !is.na(disease_group),
    !is.na(age_for_model),
    !is.na(sex)
  )
scored_data <- cp26_build_official_scores(model_ready_data)
score_names <- c(
  CP26_MANIFEST$integrated_score,
  names(CP26_MANIFEST$score_definitions)
)

score_statistics <- rp_fit_scores(scored_data, score_names)
binary_score_statistics <- rp_fit_scores(
  scored_data,
  score_names,
  binary_sex = TRUE
)
event_score_statistics <- rp_fit_scores(
  scored_data,
  score_names,
  event_col = CP26_MANIFEST$primary_event_col,
  min_events = CP26_MANIFEST$benchmarks$event_qc
)

feature_sets <- purrr::imap_dfr(
  CP26_MANIFEST$score_definitions,
  function(definition, score_name) {
    tibble::tibble(
      score = score_name,
      feature = c(definition$positive, definition$negative),
      orientation = c(
        rep(1, length(definition$positive)),
        rep(-1, length(definition$negative))
      ),
      feature_available = c(
        definition$positive,
        definition$negative
      ) %in% names(scored_data)
    )
  }
) |>
  dplyr::bind_rows(
    tibble::tibble(
      score = CP26_MANIFEST$integrated_score,
      feature = c(
        CP26_MANIFEST$integrated_definition$positive,
        CP26_MANIFEST$integrated_definition$negative
      ),
      orientation = c(
        rep(1, length(CP26_MANIFEST$integrated_definition$positive)),
        rep(-1, length(CP26_MANIFEST$integrated_definition$negative))
      ),
      feature_available = c(
        CP26_MANIFEST$integrated_definition$positive,
        CP26_MANIFEST$integrated_definition$negative
      ) %in% names(scored_data)
    )
  )

readr::write_csv(
  feature_sets,
  file.path(out_dir, "SDY2583_CP26_score_feature_map_STEP4_RECONSTRUCTED.csv")
)
readr::write_csv(
  scored_data,
  file.path(
    out_dir,
    "SDY2583_CP26_analysis_data_with_composite_scores_STEP4_RECONSTRUCTED.csv"
  )
)
readr::write_csv(
  score_statistics,
  file.path(out_dir, "SDY2583_CP26_composite_score_statistics_STEP4_RECONSTRUCTED.csv")
)
readr::write_csv(
  binary_score_statistics,
  file.path(
    out_dir,
    "SDY2583_CP26_composite_score_binary_sex_sensitivity_STEP4_RECONSTRUCTED.csv"
  )
)
readr::write_csv(
  event_score_statistics,
  file.path(
    out_dir,
    "SDY2583_CP26_composite_score_event_QC_sensitivity_STEP4_RECONSTRUCTED.csv"
  )
)

save(
  scored_data,
  score_statistics,
  binary_score_statistics,
  event_score_statistics,
  feature_sets,
  score_names,
  file = file.path(
    rdata_dir,
    "SDY2583_CP26_STEP4_composite_scores_RECONSTRUCTED.RData"
  )
)
print(score_statistics)
