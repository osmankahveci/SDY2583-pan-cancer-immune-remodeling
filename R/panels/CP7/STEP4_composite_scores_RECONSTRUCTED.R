# ============================================================
# SDY2583 CP7
# STEP 4 RECONSTRUCTED: composite scores
#
# Reconstruction basis:
#   - archived CP7 composite-score feature sets;
#   - archived subject-level score table (832 model-ready subjects);
#   - archived adjusted composite statistics.
#
# IMPORTANT: the primary CP7 composite scores were standardized in the
# valid-age, non-missing disease/sex analysis set (n = 832), not in all 850
# feature-extracted subjects. Threshold sensitivity intentionally uses its own
# all-850 standardization within each threshold set.
# ============================================================

rm(list = ls())
source(file.path(
  Sys.getenv("SDY2583_REPO_ROOT", unset = "."),
  "R", "shared", "reconstructed_panel_framework.R"
))
rp_install_and_load(c("dplyr", "readr", "tibble", "purrr"))
source(file.path(sd_repo_root(), "R", "panels", "CP7", "MANIFEST_RECONSTRUCTED.R"))

analysis_dir <- sd_analysis_dir("CP7")
rdata_dir <- file.path(analysis_dir, "11_RData")
out_dir <- file.path(analysis_dir, "08_composite_scores")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rdata_dir, recursive = TRUE, showWarnings = FALSE)

step3a_file <- file.path(
  rdata_dir,
  "SDY2583_CP7_STEP3A_metadata_merge_age_QC_RECONSTRUCTED.RData"
)
if (!file.exists(step3a_file)) stop("Run reconstructed CP7 Step 3A first: ", step3a_file)
load(step3a_file)
if (!exists("analysis_data")) stop("Step 3A RData does not contain analysis_data.")

score_definitions <- CP7_MANIFEST$score_definitions
component_scores <- names(score_definitions)
integrated_score <- CP7_MANIFEST$integrated_score
required_features <- unique(unlist(lapply(
  score_definitions,
  function(x) c(x$positive, x$negative)
)))
missing_features <- setdiff(required_features, names(analysis_data))
if (length(missing_features) > 0L) {
  stop("Required CP7 composite features are missing: ", paste(missing_features, collapse = ", "))
}

# This filter is required to reproduce the archived subject-level scores and
# the primary integrated beta of 0.351119363746345.
score_source_data <- analysis_data |>
  dplyr::filter(
    !is.na(age_for_model),
    !is.na(disease_group),
    !is.na(sex)
  )
if (nrow(score_source_data) != CP7_MANIFEST$benchmarks$n_model) {
  warning(
    "CP7 score-source row count is ", nrow(score_source_data),
    "; archived benchmark is ", CP7_MANIFEST$benchmarks$n_model, "."
  )
}

scored_data <- rp_build_scores(
  score_source_data,
  score_definitions,
  integrated_name = integrated_score
)

score_interpretations <- c(
  CP7_CD8_differentiation_remodeling_score =
    "Higher score indicates stronger CD8 differentiation / terminal-like remodeling.",
  CP7_CD39_enrichment_score =
    "Higher score indicates stronger CD39-associated regulatory/checkpoint enrichment.",
  CP7_PD1_TIGIT_attenuation_score =
    "Higher score indicates lower PD1/TIGIT expression or positivity.",
  CP7_checkpoint_coexpression_score =
    "Higher score indicates stronger multi-checkpoint co-expression enrichment.",
  CP7_terminal_checkpoint_remodeling_score =
    "Higher score indicates terminal-like checkpoint remodeling with CD39/LAG3 enrichment and PD1/TIGIT attenuation.",
  CP7_integrated_checkpoint_remodeling_score =
    "Higher score indicates integrated CP7 checkpoint-remodeling phenotype."
)

score_feature_sets <- purrr::imap_dfr(score_definitions, function(definition, score_name) {
  tibble::tibble(
    score = score_name,
    positive_features = if (length(definition$positive)) {
      paste(definition$positive, collapse = " | ")
    } else NA_character_,
    negative_features = if (length(definition$negative)) {
      paste(definition$negative, collapse = " | ")
    } else NA_character_,
    interpretation = unname(score_interpretations[[score_name]])
  )
}) |>
  dplyr::bind_rows(tibble::tibble(
    score = integrated_score,
    positive_features = "z-score mean of component scores",
    negative_features = NA_character_,
    interpretation = unname(score_interpretations[[integrated_score]])
  ))

score_names <- c(integrated_score, component_scores)
fit_score <- function(score_name) {
  x <- rp_fit_one(scored_data, score_name)
  x |>
    dplyr::transmute(
      score = score_name,
      n_model, n_healthy, n_cancer,
      healthy_mean, cancer_mean,
      healthy_median, cancer_median,
      beta_cancer_vs_healthy,
      ci_low, ci_high, t_value, p_value,
      model_formula = "score ~ disease_group + age_for_model + sex",
      direction
    )
}

composite_statistics <- purrr::map_dfr(score_names, fit_score) |>
  dplyr::mutate(fdr = stats::p.adjust(p_value, method = "BH")) |>
  dplyr::arrange(p_value)

readr::write_csv(
  score_feature_sets,
  file.path(out_dir, "SDY2583_CP7_composite_score_feature_sets_RECONSTRUCTED.csv")
)
readr::write_csv(
  scored_data,
  file.path(out_dir, "SDY2583_CP7_analysis_data_with_composite_scores_RECONSTRUCTED.csv")
)
readr::write_csv(
  composite_statistics,
  file.path(out_dir, "SDY2583_CP7_composite_score_statistics_RECONSTRUCTED.csv")
)

save(
  scored_data, score_source_data, score_definitions, score_feature_sets,
  composite_statistics,
  file = file.path(
    rdata_dir,
    "SDY2583_CP7_STEP4_composite_scores_RECONSTRUCTED.RData"
  )
)

print(composite_statistics)
message(
  "CP7 archived primary composite benchmark: n = 832; integrated beta = ",
  CP7_MANIFEST$benchmarks$integrated_beta
)
