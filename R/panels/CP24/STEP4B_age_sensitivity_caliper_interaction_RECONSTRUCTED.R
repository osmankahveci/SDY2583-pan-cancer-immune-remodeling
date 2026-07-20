# CP24 reconstructed age-stratified, same-sex matching, and interaction checks.
# These are transparent reconstructed extensions; the official CP24 archive did
# not contain direct matched/interaction comparator tables.

rm(list = ls())
source(file.path(
  Sys.getenv("SDY2583_REPO_ROOT", unset = "."),
  "R", "shared", "reconstructed_panel_framework.R"
))
rp_install_and_load(c("dplyr", "readr", "tibble", "purrr"))
source(file.path(
  sd_repo_root(), "R", "panels", "CP24", "MANIFEST_RECONSTRUCTED.R"
))

analysis_dir <- sd_analysis_dir("CP24")
rdata_dir <- file.path(analysis_dir, "11_RData")
out_dir <- file.path(analysis_dir, "10_age_sensitivity_RECONSTRUCTED")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

load(file.path(
  rdata_dir,
  "SDY2583_CP24_STEP4_composite_scores_RECONSTRUCTED.RData"
))

outcomes <- unique(c(
  CP24_MANIFEST$score_names,
  "pct_cd3_pos_total",
  "pct_cd3_cd8_pos_total",
  "pct_cd8_within_cd3",
  "median_CD62L_in_CD3CD8",
  "median_CD27_in_CD3CD8",
  "pct_temra_like",
  "pct_temra_cd27neg",
  "pct_cd62lneg_cd27neg",
  "pct_cd57_cx3cr1_cd95_pos",
  "pct_pd1_pos"
))

result <- rp_age_sensitivity(
  scored_data,
  outcomes,
  event_col = CP24_MANIFEST$primary_event_col,
  min_events = 0,
  calipers = c(5, 10)
)

readr::write_csv(
  result$stratified,
  file.path(out_dir, "SDY2583_CP24_age_stratified_statistics_RECONSTRUCTED.csv")
)
readr::write_csv(
  result$matched,
  file.path(out_dir, "SDY2583_CP24_age_matched_statistics_RECONSTRUCTED.csv")
)
readr::write_csv(
  result$match_summary,
  file.path(out_dir, "SDY2583_CP24_age_matching_summary_RECONSTRUCTED.csv")
)
readr::write_csv(
  result$interactions,
  file.path(out_dir, "SDY2583_CP24_disease_by_age_interactions_RECONSTRUCTED.csv")
)

save(
  result,
  outcomes,
  file = file.path(
    rdata_dir,
    "SDY2583_CP24_STEP4B_age_sensitivity_RECONSTRUCTED.RData"
  )
)
print(result$match_summary)
