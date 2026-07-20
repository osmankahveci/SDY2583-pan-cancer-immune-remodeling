# ============================================================
# SDY2583 CP7
# STEP 4 RECONSTRUCTED: composite scores
#
# Reconstruction basis:
#   - archived CP7 composite-score feature-set table;
#   - archived score orientation and interpretation records;
#   - archived composite-score model output schema.
#
# This is reconstructed source, not the original archived script.
# Validate subject-level scores and model coefficients against the archived
# CP7 composite outputs before treating this step as verified.
# ============================================================

rm(list = ls())
source(file.path(Sys.getenv("SDY2583_REPO_ROOT", unset = "."), "R", "shared", "bootstrap.R"))

cran_pkgs <- c("dplyr", "readr", "tibble", "purrr")
for (p in cran_pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(purrr)
})

analysis_dir <- sd_analysis_dir("CP7")
rdata_dir <- file.path(analysis_dir, "11_RData")
out_dir <- file.path(analysis_dir, "08_composite_scores")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rdata_dir, recursive = TRUE, showWarnings = FALSE)

step3a_file <- file.path(rdata_dir, "SDY2583_CP7_STEP3A_metadata_merge_age_QC_RECONSTRUCTED.RData")
if (!file.exists(step3a_file)) stop("Run reconstructed CP7 Step 3A first: ", step3a_file)
load(step3a_file)
if (!exists("analysis_data")) stop("Step 3A RData does not contain analysis_data.")

score_definitions <- list(
  CP7_CD8_differentiation_remodeling_score = list(
    positive = c(
      "pct_temra_like", "pct_temra_cd27pos", "pct_temra_cd27neg",
      "pct_cd62lneg_cd27neg", "pct_cd45ra_pos_cd62lneg_cd27neg", "pct_cd8_within_cd3"
    ),
    negative = c(
      "median_CD62L_in_CD3CD8", "pct_naive_like", "pct_tcm_like",
      "pct_cd3_pos_total", "pct_cd3_cd8_pos_total"
    ),
    interpretation = "Higher score indicates stronger CD8 differentiation / terminal-like remodeling."
  ),
  CP7_CD39_enrichment_score = list(
    positive = c(
      "median_CD39_in_CD3CD8", "pct_cd39_pos", "pct_cd39_temra",
      "pct_cd39_temra_cd27neg", "pct_cd39_pos_within_temra_like",
      "pct_cd39_pos_within_temra_cd27neg", "pct_pd1_cd39_pos",
      "pct_tigit_cd39_pos", "pct_icos_cd39_pos"
    ),
    negative = character(),
    interpretation = "Higher score indicates stronger CD39-associated regulatory/checkpoint enrichment."
  ),
  CP7_PD1_TIGIT_attenuation_score = list(
    positive = character(),
    negative = c(
      "pct_pd1_pos", "pct_pd1_high", "median_PD1_in_CD3CD8", "pct_tigit_pos",
      "median_TIGIT_in_CD3CD8", "pct_pd1_tigit_pos",
      "pct_tigit_pos_within_temra_like", "pct_tigit_pos_within_temra_cd27neg",
      "pct_pd1_pos_within_temra_like", "pct_pd1_pos_within_temra_cd27neg"
    ),
    interpretation = "Higher score indicates lower PD1/TIGIT expression or positivity."
  ),
  CP7_checkpoint_coexpression_score = list(
    positive = c(
      "pct_tim3_lag3_pos", "pct_pd1_tim3_lag3_pos", "pct_icos_cd39_pos",
      "pct_pd1_tim3_pos", "pct_pd1_cd39_pos", "pct_pd1_tigit_cd39_pos",
      "pct_pd1_tigit_icos_cd39_pos"
    ),
    negative = character(),
    interpretation = "Higher score indicates stronger multi-checkpoint co-expression enrichment."
  ),
  CP7_terminal_checkpoint_remodeling_score = list(
    positive = c(
      "pct_cd39_temra", "pct_cd39_temra_cd27neg", "pct_cd39_pos_within_temra_like",
      "pct_cd39_pos_within_temra_cd27neg", "pct_lag3_pos_within_temra_like"
    ),
    negative = c(
      "pct_tigit_pos_within_temra_like", "pct_tigit_pos_within_temra_cd27neg",
      "pct_pd1_pos_within_temra_like", "pct_pd1_pos_within_temra_cd27neg"
    ),
    interpretation = paste(
      "Higher score indicates terminal-like checkpoint remodeling with CD39/LAG3",
      "enrichment and PD1/TIGIT attenuation."
    )
  )
)

all_required_features <- unique(unlist(lapply(score_definitions, function(x) c(x$positive, x$negative))))
missing_features <- setdiff(all_required_features, names(analysis_data))
if (length(missing_features) > 0L) {
  stop("Required CP7 composite features are missing: ", paste(missing_features, collapse = ", "))
}

safe_z <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (all(is.na(x)) || stats::sd(x, na.rm = TRUE) == 0) return(rep(NA_real_, length(x)))
  as.numeric(scale(x))
}

build_score <- function(data, positive, negative) {
  oriented <- list()
  if (length(positive) > 0L) {
    for (feature in positive) oriented[[paste0("plus__", feature)]] <- safe_z(data[[feature]])
  }
  if (length(negative) > 0L) {
    for (feature in negative) oriented[[paste0("minus__", feature)]] <- -safe_z(data[[feature]])
  }
  matrix_data <- as.data.frame(oriented, check.names = FALSE)
  if (ncol(matrix_data) == 0L) return(rep(NA_real_, nrow(data)))
  score <- rowMeans(matrix_data, na.rm = TRUE)
  score[rowSums(!is.na(matrix_data)) == 0L] <- NA_real_
  score
}

scored_data <- analysis_data
for (score_name in names(score_definitions)) {
  definition <- score_definitions[[score_name]]
  scored_data[[score_name]] <- build_score(
    scored_data,
    positive = definition$positive,
    negative = definition$negative
  )
}

component_scores <- names(score_definitions)
integrated_components <- as.data.frame(
  lapply(scored_data[component_scores], safe_z),
  check.names = FALSE
)
scored_data$CP7_integrated_checkpoint_remodeling_score <- rowMeans(integrated_components, na.rm = TRUE)
scored_data$CP7_integrated_checkpoint_remodeling_score[
  rowSums(!is.na(integrated_components)) == 0L
] <- NA_real_

score_feature_sets <- purrr::imap_dfr(score_definitions, function(definition, score_name) {
  tibble(
    score = score_name,
    positive_features = paste(definition$positive, collapse = " | "),
    negative_features = paste(definition$negative, collapse = " | "),
    interpretation = definition$interpretation
  )
}) %>%
  bind_rows(tibble(
    score = "CP7_integrated_checkpoint_remodeling_score",
    positive_features = "z-score mean of component scores",
    negative_features = "",
    interpretation = "Higher score indicates integrated CP7 checkpoint-remodeling phenotype."
  ))

score_names <- c("CP7_integrated_checkpoint_remodeling_score", component_scores)

fit_score <- function(score_name) {
  dat <- scored_data %>%
    transmute(
      score = suppressWarnings(as.numeric(.data[[score_name]])),
      disease_group = factor(as.character(disease_group), levels = c("Healthy control", "Cancer patient")),
      age_for_model = suppressWarnings(as.numeric(age_for_model)),
      sex = factor(as.character(sex))
    ) %>%
    filter(!is.na(score), !is.na(disease_group), !is.na(age_for_model), !is.na(sex))

  fit <- stats::lm(score ~ disease_group + age_for_model + sex, data = dat)
  coefficient_table <- summary(fit)$coefficients
  coefficient_name <- grep("^disease_group", rownames(coefficient_table), value = TRUE)[1]
  ci <- stats::confint(fit, parm = coefficient_name, level = 0.95)
  healthy <- dat$score[dat$disease_group == "Healthy control"]
  cancer <- dat$score[dat$disease_group == "Cancer patient"]
  beta <- unname(coefficient_table[coefficient_name, "Estimate"])

  tibble(
    score = score_name,
    n_model = stats::nobs(fit),
    n_healthy = length(healthy),
    n_cancer = length(cancer),
    healthy_mean = mean(healthy, na.rm = TRUE),
    cancer_mean = mean(cancer, na.rm = TRUE),
    healthy_median = median(healthy, na.rm = TRUE),
    cancer_median = median(cancer, na.rm = TRUE),
    beta_cancer_vs_healthy = beta,
    ci_low = unname(ci[1]),
    ci_high = unname(ci[2]),
    t_value = unname(coefficient_table[coefficient_name, "t value"]),
    p_value = unname(coefficient_table[coefficient_name, "Pr(>|t|)"]),
    model_formula = "score ~ disease_group + age_for_model + sex",
    direction = ifelse(beta > 0, "higher_in_cancer", ifelse(beta < 0, "lower_in_cancer", "no_difference"))
  )
}

composite_statistics <- purrr::map_dfr(score_names, fit_score) %>%
  mutate(fdr = p.adjust(p_value, method = "BH")) %>%
  arrange(p_value)

readr::write_csv(score_feature_sets, file.path(out_dir, "SDY2583_CP7_composite_score_feature_sets_RECONSTRUCTED.csv"))
readr::write_csv(scored_data, file.path(out_dir, "SDY2583_CP7_analysis_data_with_composite_scores_RECONSTRUCTED.csv"))
readr::write_csv(composite_statistics, file.path(out_dir, "SDY2583_CP7_composite_score_statistics_RECONSTRUCTED.csv"))

save(
  scored_data, score_definitions, score_feature_sets, composite_statistics,
  file = file.path(rdata_dir, "SDY2583_CP7_STEP4_composite_scores_RECONSTRUCTED.RData")
)

print(composite_statistics)
message("Archived benchmark: all six CP7 scores were higher in cancer; integrated beta approximately 0.3511.")
