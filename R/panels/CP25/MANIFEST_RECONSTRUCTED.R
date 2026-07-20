# CP25 reconstructed manifest from archived thresholds, feature schema,
# composite definitions, and final Methods/Results record.

CP25_FEATURES <- c(
  "pct_cd3_pos_total", "pct_cd4_pos_total", "pct_dump_low_total",
  "pct_cd3_cd4_raw_total", "pct_cd3_cd4_primary_total",
  "pct_cd4_within_cd3_raw", "pct_cd4_within_cd3_primary",
  "pct_dump_low_within_cd3_cd4_raw",
  "median_CD45RA_in_CD3CD4", "median_CD27_in_CD3CD4",
  "median_CD25_in_CD3CD4", "median_IL7RA_in_CD3CD4",
  "median_CTLA4_in_CD3CD4", "median_ICOS_in_CD3CD4",
  "median_CD39_in_CD3CD4", "median_PD1_in_CD3CD4",
  "median_LAG3_in_CD3CD4",
  "pct_naive_like", "pct_memory_like", "pct_memory_cd27pos",
  "pct_memory_cd27neg", "pct_cd45ra_cd27neg", "pct_cd27_pos",
  "pct_cd27_neg",
  "pct_cd25_pos", "pct_cd25_high", "pct_il7ra_pos", "pct_il7ra_low",
  "pct_cd25pos_il7ralow_treg_like",
  "pct_cd25high_il7ralow_treg_enriched",
  "pct_cd25pos_il7ralow_within_memory",
  "pct_ctla4_pos", "pct_icos_pos", "pct_cd39_pos", "pct_pd1_pos",
  "pct_pd1_high", "pct_lag3_pos", "pct_any_checkpoint_reg_marker",
  "pct_any_inhibitory_checkpoint",
  "pct_pd1_lag3_pos", "pct_pd1_cd39_pos", "pct_lag3_cd39_pos",
  "pct_ctla4_cd39_pos", "pct_icos_cd39_pos", "pct_ctla4_icos_pos",
  "pct_pd1_icos_pos",
  "pct_ctla4_treg_like", "pct_cd39_treg_like", "pct_icos_treg_like",
  "pct_pd1_treg_like", "pct_lag3_treg_like",
  "pct_ctla4_cd39_treg_like", "pct_icos_cd39_treg_like",
  "pct_pd1_lag3_treg_like",
  "pct_ctla4pos_within_treg_like", "pct_cd39pos_within_treg_like",
  "pct_icospos_within_treg_like", "pct_pd1pos_within_treg_like",
  "pct_lag3pos_within_treg_like",
  "pct_ctla4_cd39pos_within_treg_like",
  "pct_icos_cd39pos_within_treg_like",
  "pct_pd1_lag3pos_within_treg_like",
  "pct_ctla4pos_within_treg_enriched",
  "pct_cd39pos_within_treg_enriched",
  "pct_icospos_within_treg_enriched",
  "pct_pd1pos_within_treg_enriched",
  "pct_lag3pos_within_treg_enriched",
  "pct_ctla4pos_within_memory", "pct_cd39pos_within_memory",
  "pct_icospos_within_memory", "pct_pd1pos_within_memory",
  "pct_lag3pos_within_memory",
  "pct_ctla4pos_within_cd39pos", "pct_icospos_within_cd39pos",
  "pct_pd1pos_within_cd39pos", "pct_lag3pos_within_cd39pos",
  "pct_lag3pos_within_pd1pos", "pct_cd39pos_within_icospos"
)

module_for <- function(feature) {
  dplyr::case_when(
    grepl(
      "^pct_(cd3|cd4|dump_low|cd3_cd4)|pct_cd4_within|pct_dump_low_within",
      feature
    ) ~ "composition",
    grepl("^median_", feature) ~ "marker_medians",
    feature %in% c(
      "pct_naive_like", "pct_memory_like", "pct_memory_cd27pos",
      "pct_memory_cd27neg", "pct_cd45ra_cd27neg", "pct_cd27_pos",
      "pct_cd27_neg"
    ) ~ "differentiation",
    feature %in% c(
      "pct_cd25_pos", "pct_cd25_high", "pct_il7ra_pos", "pct_il7ra_low",
      "pct_cd25pos_il7ralow_treg_like",
      "pct_cd25high_il7ralow_treg_enriched",
      "pct_cd25pos_il7ralow_within_memory"
    ) ~ "regulatory_like",
    grepl("_treg_like$", feature) ~ "checkpoint_positive_regulatory_like",
    grepl("within_treg_like", feature) ~ "within_regulatory_like",
    grepl("within_treg_enriched", feature) ~ "within_treg_enriched",
    grepl("within_memory", feature) ~ "within_memory",
    grepl(
      "within_cd39pos|within_pd1pos|within_icospos",
      feature
    ) ~ "conditional_checkpoint",
    grepl(
      "pd1_lag3|pd1_cd39|lag3_cd39|ctla4_cd39|icos_cd39|ctla4_icos|pd1_icos",
      feature
    ) ~ "checkpoint_coexpression",
    TRUE ~ "single_checkpoint_regulatory"
  )
}

# Exact archived threshold-sensitivity family: 64 features plus 9 scores.
CP25_THRESHOLD_TARGET_FEATURES <- c(
  "pct_cd3_cd4_raw_total", "pct_cd3_cd4_primary_total",
  "pct_cd3_pos_total", "pct_cd4_pos_total", "pct_cd25_high",
  "pct_il7ra_pos", "pct_il7ra_low", "pct_cd39_treg_like",
  "pct_cd25pos_il7ralow_treg_like", "pct_lag3_treg_like",
  "pct_cd25high_il7ralow_treg_enriched",
  "pct_icospos_within_treg_enriched", "pct_cd4_within_cd3_primary",
  "pct_icos_treg_like", "pct_icos_cd39_treg_like",
  "pct_pd1_lag3_treg_like", "pct_cd4_within_cd3_raw",
  "median_IL7RA_in_CD3CD4", "pct_cd25pos_il7ralow_within_memory",
  "pct_pd1_treg_like", "pct_lag3pos_within_pd1pos",
  "pct_icospos_within_treg_like", "pct_cd25_pos", "pct_cd27_pos",
  "pct_cd27_neg", "pct_memory_cd27neg",
  "pct_lag3pos_within_treg_enriched", "pct_naive_like",
  "pct_ctla4_cd39_treg_like", "pct_pd1_lag3_pos", "pct_memory_like",
  "pct_lag3pos_within_treg_like", "pct_lag3_cd39_pos", "pct_lag3_pos",
  "pct_any_inhibitory_checkpoint", "pct_cd45ra_cd27neg",
  "median_LAG3_in_CD3CD4", "median_CD45RA_in_CD3CD4",
  "pct_ctla4_cd39_pos", "pct_dump_low_within_cd3_cd4_raw",
  "pct_pd1_cd39_pos", "pct_cd39pos_within_icospos",
  "pct_ctla4_treg_like", "pct_cd39_pos", "pct_pd1_high",
  "pct_icos_cd39_pos", "median_PD1_in_CD3CD4",
  "pct_ctla4pos_within_treg_enriched", "pct_pd1_icos_pos",
  "pct_ctla4pos_within_treg_like", "pct_pd1_pos",
  "median_CD25_in_CD3CD4", "median_CD39_in_CD3CD4",
  "pct_any_checkpoint_reg_marker", "pct_ctla4_icos_pos",
  "pct_ctla4_pos", "pct_cd39pos_within_treg_like",
  "pct_cd39pos_within_treg_enriched", "median_CD27_in_CD3CD4",
  "pct_icos_pos", "median_CTLA4_in_CD3CD4", "median_ICOS_in_CD3CD4",
  "pct_pd1pos_within_treg_like", "pct_pd1pos_within_treg_enriched"
)

CP25_MANIFEST <- list(
  panel = "CP25",
  primary_event_col = "n_cd3_cd4_primary",
  thresholds = c(
    CD3 = 2, CD4 = 2, DUMP_LOW = 1.5, CD45RA = 2, CD27 = 2,
    CD25 = 1.5, CD25_HIGH = 2, IL7RA = 1.5, CTLA4 = 1.5,
    ICOS = 1.5, CD39 = 1.5, PD1 = 1.5, PD1_HIGH = 2, LAG3 = 1.5
  ),
  expected_markers = c(
    "CD3", "CD4", "Viability_CD8_CD13_CD19_TCRgd", "CD45RA", "CD27",
    "CD25", "IL7RA", "CTLA4", "ICOS", "CD39", "PD-1", "LAG-3"
  ),
  module_map = tibble::tibble(
    feature = CP25_FEATURES,
    module = vapply(CP25_FEATURES, module_for, character(1))
  ),
  score_definitions = list(
    CP25_CD4_composition_depletion_score = list(
      positive = character(),
      negative = c(
        "pct_cd3_pos_total", "pct_cd4_pos_total", "pct_cd3_cd4_raw_total",
        "pct_cd3_cd4_primary_total", "pct_cd4_within_cd3_raw",
        "pct_cd4_within_cd3_primary", "pct_dump_low_within_cd3_cd4_raw"
      )
    ),
    CP25_CD45RA_CD27_differentiation_score = list(
      positive = c(
        "pct_memory_like", "pct_memory_cd27neg", "pct_cd45ra_cd27neg",
        "pct_cd27_neg"
      ),
      negative = c(
        "pct_naive_like", "pct_cd27_pos", "median_CD45RA_in_CD3CD4",
        "median_CD27_in_CD3CD4"
      )
    ),
    CP25_regulatory_like_enrichment_score = list(
      positive = c(
        "pct_cd25_pos", "pct_cd25_high", "pct_il7ra_low",
        "pct_cd25pos_il7ralow_treg_like",
        "pct_cd25high_il7ralow_treg_enriched",
        "pct_cd25pos_il7ralow_within_memory"
      ),
      negative = c("pct_il7ra_pos", "median_IL7RA_in_CD3CD4")
    ),
    CP25_total_checkpoint_regulatory_marker_score = list(
      positive = c(
        "pct_ctla4_pos", "pct_icos_pos", "pct_cd39_pos", "pct_pd1_pos",
        "pct_pd1_high", "pct_lag3_pos", "pct_any_checkpoint_reg_marker",
        "pct_any_inhibitory_checkpoint", "median_CTLA4_in_CD3CD4",
        "median_ICOS_in_CD3CD4", "median_CD39_in_CD3CD4",
        "median_PD1_in_CD3CD4", "median_LAG3_in_CD3CD4"
      ),
      negative = character()
    ),
    CP25_checkpoint_coexpression_score = list(
      positive = c(
        "pct_pd1_lag3_pos", "pct_pd1_cd39_pos", "pct_lag3_cd39_pos",
        "pct_ctla4_cd39_pos", "pct_icos_cd39_pos", "pct_ctla4_icos_pos",
        "pct_pd1_icos_pos"
      ),
      negative = character()
    ),
    CP25_regulatory_like_checkpoint_enrichment_score = list(
      positive = c(
        "pct_ctla4_treg_like", "pct_cd39_treg_like", "pct_icos_treg_like",
        "pct_pd1_treg_like", "pct_lag3_treg_like",
        "pct_ctla4_cd39_treg_like", "pct_icos_cd39_treg_like",
        "pct_pd1_lag3_treg_like"
      ),
      negative = character()
    ),
    CP25_CD39_ICOS_regulatory_axis_score = list(
      positive = c(
        "pct_cd39_pos", "pct_icos_pos", "pct_icos_cd39_pos",
        "pct_cd39_treg_like", "pct_icos_treg_like",
        "pct_icos_cd39_treg_like", "median_CD39_in_CD3CD4",
        "median_ICOS_in_CD3CD4"
      ),
      negative = character()
    ),
    CP25_PD1_LAG3_regulatory_axis_score = list(
      positive = c(
        "pct_pd1_pos", "pct_pd1_high", "pct_lag3_pos", "pct_pd1_lag3_pos",
        "pct_pd1_treg_like", "pct_lag3_treg_like",
        "pct_pd1_lag3_treg_like", "pct_lag3pos_within_pd1pos",
        "median_PD1_in_CD3CD4", "median_LAG3_in_CD3CD4"
      ),
      negative = character()
    )
  ),
  integrated_score = "CP25_integrated_CD4_regulatory_checkpoint_remodeling_score",
  targeted_outcomes = CP25_FEATURES,
  threshold_target_features = CP25_THRESHOLD_TARGET_FEATURES,
  benchmarks = list(
    n_files = 850,
    n_model = 832,
    median_total_events = 302748,
    median_raw_events = 26488,
    median_primary_events = 16027,
    event_qc = 300,
    matched_pairs_5y = 265,
    matched_pairs_10y = 273,
    integrated_beta = 0.4392101623,
    robust_both = 48,
    threshold_direction_only = 23,
    threshold_direction_fail = 2,
    targeted_n = 73
  )
)
