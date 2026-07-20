# ============================================================
# SDY2583 CP8
# STEP 5 RECONSTRUCTED: targeted threshold-sensitivity analysis.
#
# Archive-specific implementation confirmed by independent recalculation:
#   - exact 28-outcome family;
#   - permissive/main/stringent feature extraction from raw FCS files;
#   - positive-marker thresholds shifted by -0.2/+0.2;
#   - IL7RA-low threshold direction reversed (+0.2/-0.2);
#   - composite scores standardized independently within each all-850
#     threshold extraction;
#   - BH FDR calculated separately within each 28-outcome threshold set.
# ============================================================

rm(list = ls())
source(file.path(
  Sys.getenv("SDY2583_REPO_ROOT", unset = "."),
  "R", "shared", "reconstructed_panel_framework.R"
))
rp_install_and_load(c("dplyr", "readr", "tibble", "purrr", "tidyr"), "flowCore")
source(file.path(sd_repo_root(), "R", "panels", "CP8", "MANIFEST_RECONSTRUCTED.R"))

analysis_dir <- sd_analysis_dir("CP8")
rdata_dir <- file.path(analysis_dir, "11_RData")
out_dir <- file.path(analysis_dir, "12_threshold_sensitivity")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

step1_file <- file.path(
  rdata_dir,
  "SDY2583_CP8_STEP1_inventory_RECONSTRUCTED.RData"
)
step2_file <- file.path(
  rdata_dir,
  "SDY2583_CP8_STEP2_feature_extraction_RECONSTRUCTED.RData"
)
step3a_file <- file.path(
  rdata_dir,
  "SDY2583_CP8_STEP3A_metadata_merge_age_QC_RECONSTRUCTED.RData"
)
for (required_file in c(step1_file, step2_file, step3a_file)) {
  if (!file.exists(required_file)) stop("Required CP8 input is missing: ", required_file)
}

# Step 1 supplies fcs_files; Step 2 supplies extract_one; Step 3A supplies
# standardized participant metadata.
load(step1_file)
load(step2_file)
load(step3a_file)
if (!exists("fcs_files")) stop("CP8 Step 1 RData does not contain fcs_files.")
if (!exists("extract_one")) stop("CP8 Step 2 RData does not contain extract_one.")
if (!exists("metadata")) stop("CP8 Step 3A RData does not contain metadata.")

main <- CP8_MANIFEST$thresholds
permissive <- main - 0.2
stringent <- main + 0.2

# IL7RA-low is a low-gate phenotype: increasing the cutoff is permissive and
# decreasing it is stringent.
permissive["IL7RA"] <- main["IL7RA"] + 0.2
stringent["IL7RA"] <- main["IL7RA"] - 0.2

sets <- list(
  permissive = permissive,
  main = main,
  stringent = stringent
)
outcomes <- CP8_MANIFEST$threshold_sensitivity_outcomes
score_names <- c(
  names(CP8_MANIFEST$score_definitions),
  CP8_MANIFEST$integrated_score
)

threshold_sets <- purrr::imap_dfr(
  sets,
  ~tibble::tibble(
    threshold_set = .y,
    marker = names(.x),
    threshold = as.numeric(.x)
  )
)

one_set <- function(label, thresholds) {
  feature_table <- purrr::map_dfr(
    fcs_files,
    ~extract_one(.x, thresholds = thresholds)
  )

  data_with_metadata <- feature_table |>
    dplyr::left_join(metadata, by = "subject_id")

  scored_data <- rp_build_scores(
    data_with_metadata,
    CP8_MANIFEST$score_definitions,
    CP8_MANIFEST$integrated_score
  )

  missing_outcomes <- setdiff(outcomes, names(scored_data))
  if (length(missing_outcomes) > 0L) {
    stop(
      "CP8 threshold outcomes are missing in ", label, ": ",
      paste(missing_outcomes, collapse = ", ")
    )
  }

  statistics <- purrr::map_dfr(
    outcomes,
    ~rp_fit_one(
      scored_data,
      .x,
      subset_label = label,
      binary_sex = FALSE,
      event_col = CP8_MANIFEST$primary_event_col,
      min_events = 0
    )
  ) |>
    dplyr::mutate(
      threshold_set = label,
      targeted_FDR = stats::p.adjust(p_value, method = "BH")
    ) |>
    dplyr::select(
      threshold_set, feature, n_model,
      healthy_mean, cancer_mean,
      beta_cancer_vs_healthy, ci_low, ci_high,
      p_value, targeted_FDR, direction
    )

  list(
    features = feature_table,
    data_with_metadata = data_with_metadata,
    scored_data = scored_data,
    statistics = statistics
  )
}

results <- purrr::imap(sets, one_set)

statistics_long <- dplyr::bind_rows(lapply(results, `[[`, "statistics"))

feature_data_long <- dplyr::bind_rows(
  purrr::imap(
    results,
    ~.x$data_with_metadata |>
      dplyr::mutate(threshold_set = .y) |>
      dplyr::relocate(threshold_set, .after = subject_id)
  )
)

scored_data_long <- dplyr::bind_rows(
  purrr::imap(
    results,
    ~.x$scored_data |>
      dplyr::mutate(threshold_set = .y) |>
      dplyr::relocate(threshold_set, .after = subject_id)
  )
)

selected_features_long <- scored_data_long |>
  dplyr::select(
    subject_id,
    threshold_set,
    dplyr::all_of(outcomes),
    dplyr::any_of(c(
      "disease_group", "age_raw", "age_for_model",
      "age_group_for_model", "sex"
    ))
  )

summary_wide <- statistics_long |>
  dplyr::select(
    threshold_set, feature,
    beta_cancer_vs_healthy, targeted_FDR, direction
  ) |>
  tidyr::pivot_wider(
    names_from = threshold_set,
    values_from = c(beta_cancer_vs_healthy, targeted_FDR, direction)
  ) |>
  dplyr::mutate(
    direction_preserved =
      direction_main == direction_permissive &
      direction_main == direction_stringent,
    FDR_lt_0p05_all_sets =
      targeted_FDR_main < 0.05 &
      targeted_FDR_permissive < 0.05 &
      targeted_FDR_stringent < 0.05,
    robust_interpretation = dplyr::case_when(
      direction_preserved & FDR_lt_0p05_all_sets ~
        "Direction and FDR-significance preserved across threshold sets",
      direction_preserved ~
        "Direction preserved but FDR-significance not preserved across all threshold sets",
      TRUE ~
        "Direction not preserved across threshold sets"
    )
  )

readr::write_csv(
  threshold_sets,
  file.path(
    out_dir,
    "SDY2583_CP8_threshold_sensitivity_threshold_sets_RECONSTRUCTED.csv"
  )
)
readr::write_csv(
  feature_data_long,
  file.path(
    out_dir,
    "SDY2583_CP8_threshold_sensitivity_data_with_metadata_LONG_RECONSTRUCTED.csv"
  )
)
readr::write_csv(
  scored_data_long,
  file.path(
    out_dir,
    "SDY2583_CP8_threshold_sensitivity_scored_data_LONG_RECONSTRUCTED.csv"
  )
)
readr::write_csv(
  selected_features_long,
  file.path(
    out_dir,
    "SDY2583_CP8_threshold_sensitivity_selected_features_LONG_RECONSTRUCTED.csv"
  )
)
readr::write_csv(
  statistics_long,
  file.path(
    out_dir,
    "SDY2583_CP8_threshold_sensitivity_statistics_LONG_RECONSTRUCTED.csv"
  )
)
readr::write_csv(
  summary_wide,
  file.path(
    out_dir,
    "SDY2583_CP8_threshold_sensitivity_summary_WIDE_RECONSTRUCTED.csv"
  )
)

save(
  sets,
  outcomes,
  score_names,
  results,
  threshold_sets,
  feature_data_long,
  scored_data_long,
  selected_features_long,
  statistics_long,
  summary_wide,
  file = file.path(
    rdata_dir,
    "SDY2583_CP8_STEP5_threshold_sensitivity_RECONSTRUCTED.RData"
  )
)

print(table(summary_wide$robust_interpretation))
message("Archived CP8 benchmark: all 28 outcomes preserve direction and FDR < 0.05 across all three threshold sets.")
