# ============================================================
# SDY2583 CP7
# STEP 5B RECONSTRUCTED: archived targeted threshold family
#
# Step 5 performs the three FCS extractions and retains an intentionally broad
# set of candidate outcomes. The archived official targeted-FDR tables used an
# exact 11-outcome family. This step selects that family, recomputes BH FDR
# within each threshold set, and writes archive-compatible official outputs.
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
out_dir <- file.path(analysis_dir, "12_threshold_sensitivity")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

step5_file <- file.path(
  rdata_dir,
  "SDY2583_CP7_STEP5_threshold_sensitivity_RECONSTRUCTED.RData"
)
if (!file.exists(step5_file)) stop("Run reconstructed CP7 Step 5 first: ", step5_file)
load(step5_file)
if (!exists("results") || !exists("statistics_long")) {
  stop("Step 5 RData must contain results and statistics_long.")
}

threshold_targets <- CP7_MANIFEST$threshold_targeted_outcomes
if (length(threshold_targets) != CP7_MANIFEST$benchmarks$threshold_target_count) {
  stop("CP7 threshold target family must contain exactly 11 outcomes.")
}
missing_targets <- setdiff(threshold_targets, unique(statistics_long$feature))
if (length(missing_targets) > 0L) {
  stop("Threshold statistics are missing: ", paste(missing_targets, collapse = ", "))
}

threshold_order <- c("permissive", "main", "stringent")
statistics_official <- statistics_long |>
  dplyr::filter(feature %in% threshold_targets) |>
  dplyr::mutate(
    threshold_set = factor(threshold_set, levels = threshold_order),
    feature = factor(feature, levels = threshold_targets)
  ) |>
  dplyr::group_by(threshold_set) |>
  dplyr::mutate(targeted_FDR = stats::p.adjust(p_value, method = "BH")) |>
  dplyr::ungroup() |>
  dplyr::arrange(threshold_set, feature) |>
  dplyr::transmute(
    threshold_set = as.character(threshold_set),
    feature = as.character(feature),
    n_model,
    healthy_mean,
    cancer_mean,
    beta_cancer_vs_healthy,
    ci_low,
    ci_high,
    p_value,
    targeted_FDR,
    direction
  )

summary_wide_official <- statistics_official |>
  dplyr::select(
    threshold_set, feature, beta_cancer_vs_healthy,
    targeted_FDR, direction
  ) |>
  tidyr::pivot_wider(
    names_from = threshold_set,
    values_from = c(beta_cancer_vs_healthy, targeted_FDR, direction)
  ) |>
  dplyr::mutate(
    feature = factor(feature, levels = threshold_targets),
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
      TRUE ~ "Direction not preserved across threshold sets"
    )
  ) |>
  dplyr::arrange(feature) |>
  dplyr::mutate(feature = as.character(feature)) |>
  dplyr::select(
    feature,
    beta_cancer_vs_healthy_main,
    targeted_FDR_main,
    direction_main,
    beta_cancer_vs_healthy_permissive,
    targeted_FDR_permissive,
    direction_permissive,
    beta_cancer_vs_healthy_stringent,
    targeted_FDR_stringent,
    direction_stringent,
    direction_preserved,
    FDR_lt_0p05_all_sets,
    robust_interpretation
  )

scored_data_long <- purrr::imap_dfr(results, function(x, threshold_label) {
  x$data |>
    dplyr::mutate(threshold_set = threshold_label, .after = subject_id)
}) |>
  dplyr::mutate(
    threshold_set = factor(threshold_set, levels = threshold_order)
  ) |>
  dplyr::arrange(threshold_set, subject_id) |>
  dplyr::mutate(threshold_set = as.character(threshold_set))

readr::write_csv(
  statistics_official,
  file.path(
    out_dir,
    "SDY2583_CP7_threshold_sensitivity_statistics_LONG_RECONSTRUCTED.csv"
  )
)
readr::write_csv(
  summary_wide_official,
  file.path(
    out_dir,
    "SDY2583_CP7_threshold_sensitivity_summary_WIDE_RECONSTRUCTED.csv"
  )
)
readr::write_csv(
  scored_data_long,
  file.path(
    out_dir,
    "SDY2583_CP7_threshold_sensitivity_scored_data_LONG_RECONSTRUCTED.csv"
  )
)

save(
  threshold_targets, statistics_official, summary_wide_official,
  scored_data_long,
  file = file.path(
    rdata_dir,
    "SDY2583_CP7_STEP5B_archive_target_family_RECONSTRUCTED.RData"
  )
)

print(summary_wide_official)
message(
  "CP7 official threshold family reproduced with ",
  nrow(summary_wide_official), " outcomes."
)
