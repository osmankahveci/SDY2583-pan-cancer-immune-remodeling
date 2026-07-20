# CP24 full-cohort reconstructed quantitative figures.

rm(list = ls())
source(file.path(
  Sys.getenv("SDY2583_REPO_ROOT", unset = "."),
  "R", "shared", "figure_framework.R"
))
rp_install_and_load(c(
  "dplyr", "readr", "tidyr", "ggplot2", "patchwork", "stringr"
))
source(file.path(
  sd_repo_root(), "R", "panels", "CP24", "MANIFEST_RECONSTRUCTED.R"
))

analysis_dir <- sd_analysis_dir("CP24")
rdata_dir <- file.path(analysis_dir, "11_RData")
out_dir <- file.path(analysis_dir, "10_figures_RECONSTRUCTED")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

load(file.path(
  rdata_dir,
  "SDY2583_CP24_STEP4_composite_scores_RECONSTRUCTED.RData"
))
load(file.path(
  rdata_dir,
  "SDY2583_CP24_STEP3B_statistics_RECONSTRUCTED.RData"
))

main_figure_statistics <- main_statistics |>
  dplyr::rename(feature = outcome)
score_figure_statistics <- score_statistics |>
  dplyr::rename(feature = score)

rp_quantitative_figure(
  scored_data,
  c("pct_cd3_pos_total", "pct_cd3_cd8_pos_total", "pct_cd8_within_cd3"),
  main_figure_statistics,
  file.path(out_dir, "Figure1_CP24_composition_RECONSTRUCTED.pdf")
)
rp_quantitative_figure(
  scored_data,
  c(
    "median_CD62L_in_CD3CD8",
    "median_CD27_in_CD3CD8",
    "pct_naive_like",
    "pct_temra_like",
    "pct_temra_cd27neg",
    "pct_cd62lneg_cd27neg"
  ),
  main_figure_statistics,
  file.path(out_dir, "Figure2_CP24_terminal_differentiation_RECONSTRUCTED.pdf")
)
rp_quantitative_figure(
  scored_data,
  c(
    "pct_cd57_pos",
    "pct_cx3cr1_pos",
    "pct_cd95_pos",
    "pct_cd57_cx3cr1_cd95_pos",
    "pct_pd1_pos"
  ),
  main_figure_statistics,
  file.path(out_dir, "Figure3_CP24_extended_PD1_RECONSTRUCTED.pdf")
)
rp_quantitative_figure(
  scored_data,
  CP24_MANIFEST$score_names,
  score_figure_statistics,
  file.path(out_dir, "Figure4_CP24_composite_scores_RECONSTRUCTED.pdf"),
  ncol = 2
)
