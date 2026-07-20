# CP24 full-cohort reconstructed manifest.
#
# The recovered codebook remains a separate pilot/QC artifact. This manifest
# reproduces the archived 850-subject official clean-age analysis schema.

CP24_MODULES <- list(
  composition = c(
    "pct_cd3_pos_total",
    "pct_cd8_pos_total",
    "pct_cd3_cd8_pos_total",
    "pct_cd8_within_cd3"
  ),
  marker_medians = c(
    "median_CD62L_in_CD3CD8",
    "median_CD27_in_CD3CD8",
    "median_PD1_in_CD3CD8",
    "median_CD57_in_CD3CD8",
    "median_CX3CR1_in_CD3CD8",
    "median_CD95_in_CD3CD8",
    "median_CD45RA_in_CD3CD8",
    "median_CXCR3_in_CD3CD8",
    "median_CXCR5_in_CD3CD8"
  ),
  cd45ra_cd62l_quadrants = c(
    "pct_naive_like",
    "pct_tcm_like",
    "pct_tem_like",
    "pct_temra_like",
    "pct_cd62l_pos_within_cd3cd8",
    "pct_cd62l_neg_within_cd3cd8",
    "pct_cd45ra_pos_within_cd3cd8",
    "pct_cd45ra_neg_within_cd3cd8"
  ),
  cd27_refined_temra = c(
    "pct_cd27_pos_within_cd3cd8",
    "pct_cd27_neg_within_cd3cd8",
    "pct_temra_cd27neg",
    "pct_temra_cd27pos",
    "pct_cd62lneg_cd27neg",
    "pct_cd45ra_pos_cd62lneg_cd27neg",
    "pct_naive_cd27pos",
    "pct_naive_cd27neg"
  ),
  extended_differentiation = c(
    "pct_cd57_pos",
    "pct_cx3cr1_pos",
    "pct_cd95_pos",
    "pct_cd57_cx3cr1_pos",
    "pct_cd57_cd95_pos",
    "pct_cx3cr1_cd95_pos",
    "pct_cd57_cx3cr1_cd95_pos"
  ),
  pd1_module = c(
    "median_PD1_in_CD3CD8",
    "pct_pd1_pos",
    "pct_pd1_high",
    "pct_pd1_pos_within_naive_like",
    "pct_pd1_pos_within_tcm_like",
    "pct_pd1_pos_within_tem_like",
    "pct_pd1_pos_within_temra_like",
    "pct_pd1_pos_within_temra_cd27neg",
    "pct_pd1_temra",
    "pct_pd1_temra_cd27neg",
    "pct_pd1high_temra",
    "pct_pd1high_temra_cd27neg",
    "pct_pd1_cd57_pos",
    "pct_pd1_cx3cr1_pos",
    "pct_pd1_cd95_pos",
    "pct_pd1_cd57_cx3cr1_pos",
    "pct_pd1_cd57_cd95_pos",
    "pct_pd1_cx3cr1_cd95_pos",
    "pct_pd1_cd57_cx3cr1_cd95_pos"
  )
)

# This is deliberately a row-level map: median PD-1 belongs to both the marker
# median and PD-1 modules in the archived 55-test family.
CP24_MODULE_MAP <- dplyr::bind_rows(lapply(names(CP24_MODULES), function(module_name) {
  tibble::tibble(module = module_name, feature = CP24_MODULES[[module_name]])
}))

CP24_RAW_SCORE_DEFINITIONS <- list(
  CP24_marker_loss_score = c(
    median_CD62L_in_CD3CD8 = -1,
    median_CD27_in_CD3CD8 = -1
  ),
  CP24_quadrant_shift_score = c(
    pct_naive_like = -1,
    pct_temra_like = 1,
    pct_cd62l_neg_within_cd3cd8 = 1
  ),
  CP24_CD27_refined_terminal_score = c(
    pct_cd27_neg_within_cd3cd8 = 1,
    pct_temra_cd27neg = 1,
    pct_cd62lneg_cd27neg = 1,
    pct_cd45ra_pos_cd62lneg_cd27neg = 1
  ),
  CP24_extended_effector_score = c(
    pct_cd57_pos = 1,
    pct_cx3cr1_pos = 1,
    pct_cd95_pos = 1,
    pct_cd57_cx3cr1_cd95_pos = 1
  ),
  CP24_PD1_associated_terminal_score = c(
    pct_pd1_temra_cd27neg = 1,
    pct_pd1high_temra_cd27neg = 1,
    pct_pd1_cd57_cx3cr1_cd95_pos = 1
  )
)

CP24_CORE_DOMAIN_SCORES <- c(
  "CP24_marker_loss_score",
  "CP24_quadrant_shift_score",
  "CP24_CD27_refined_terminal_score",
  "CP24_extended_effector_score"
)

CP24_SCORE_NAMES <- c(
  names(CP24_RAW_SCORE_DEFINITIONS),
  "CP24_core_differentiation_score",
  "CP24_integrated_remodeling_score"
)

cp24_z_safe <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  s <- stats::sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - m) / s
}

cp24_row_mean_safe <- function(df, columns) {
  columns <- columns[columns %in% names(df)]
  if (length(columns) == 0L) return(rep(NA_real_, nrow(df)))
  mat <- as.data.frame(df[, columns, drop = FALSE], check.names = FALSE)
  value <- rowMeans(mat, na.rm = TRUE)
  value[rowSums(!is.na(mat)) == 0L] <- NA_real_
  value
}

cp24_build_official_scores <- function(data) {
  scored <- data

  for (score_name in names(CP24_RAW_SCORE_DEFINITIONS)) {
    definition <- CP24_RAW_SCORE_DEFINITIONS[[score_name]]
    z_columns <- character()

    for (feature_name in names(definition)) {
      if (!(feature_name %in% names(scored))) next
      direction <- as.numeric(definition[[feature_name]])
      z_name <- paste0("z_oriented__", feature_name)
      scored[[z_name]] <- cp24_z_safe(scored[[feature_name]]) * direction
      z_columns <- c(z_columns, z_name)
    }

    scored[[score_name]] <- cp24_row_mean_safe(scored, z_columns)
  }

  scored$CP24_core_differentiation_score <- cp24_row_mean_safe(
    scored,
    CP24_CORE_DOMAIN_SCORES
  )
  scored$CP24_integrated_remodeling_score <- cp24_row_mean_safe(
    scored,
    c(
      "CP24_core_differentiation_score",
      "CP24_PD1_associated_terminal_score"
    )
  )

  scored
}

CP24_MANIFEST <- list(
  panel = "CP24",
  primary_event_col = "n_cd3_cd8_pos",
  thresholds = c(
    CD3 = 2,
    CD8 = 2.2,
    CD45RA = 2,
    CD62L = 2.2,
    CD27 = 2,
    CD57 = 1.5,
    CX3CR1 = 1.5,
    CD95 = 1.5,
    PD1 = 1.5,
    PD1_HIGH = 2,
    CXCR3 = 1.5,
    CXCR5 = 1.5
  ),
  expected_markers = c(
    "CD3", "CD8", "CD45RA", "CD62L", "CD27", "CD57",
    "CX3CR1", "CD95", "PD-1", "CXCR3", "CXCR5"
  ),
  module_map = CP24_MODULE_MAP,
  # Retained for compatibility with downstream script discovery. The exact
  # archived hierarchy is implemented by cp24_build_official_scores().
  score_definitions = CP24_RAW_SCORE_DEFINITIONS,
  score_names = CP24_SCORE_NAMES,
  integrated_score = "CP24_integrated_remodeling_score",
  targeted_outcomes = unique(unlist(CP24_MODULES)),
  benchmarks = list(
    n_files = 850,
    n_clean_age = 832,
    n_lm_complete = 828,
    n_healthy_clean_age = 398,
    n_cancer_clean_age = 434,
    median_total_events = 302822,
    median_primary_events = 11836,
    event_n_all = 850,
    event_n_500 = 849,
    event_n_1000 = 844,
    event_model_n_all = 832,
    event_model_n_500 = 831,
    event_model_n_1000 = 826,
    integrated_beta = 0.113850267291264,
    core_beta = 0.167915840013875
  )
)
