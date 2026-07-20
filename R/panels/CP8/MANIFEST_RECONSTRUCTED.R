# CP8 reconstructed manifest. Values are transcribed from archived CP8
# threshold, module-map, composite-definition, age-sensitivity, and
# threshold-sensitivity records.

CP8_MANIFEST <- list(
  panel = "CP8",
  primary_event_col = "n_cd3_cd4_pos",
  thresholds = c(
    CD3 = 2.0, CD4 = 2.0, CD45RA = 2.0, CD62L = 2.2, CD27 = 2.0,
    CCR6 = 1.5, IL7RA = 1.5, CCR4 = 1.5, CD25 = 1.5,
    CD25_HIGH = 2.0, CXCR5 = 1.5
  ),
  expected_markers = c(
    "CD3", "CD4", "CD45RA", "CD62L", "CD27", "CCR6",
    "IL7RA", "CCR4", "CD25", "CXCR5"
  ),
  module_map = tibble::tribble(
    ~module, ~feature,
    "composition", "pct_cd3_pos_total",
    "composition", "pct_cd4_pos_total",
    "composition", "pct_cd3_cd4_pos_total",
    "composition", "pct_cd4_within_cd3",
    "marker_medians", "median_CD45RA_in_CD3CD4",
    "marker_medians", "median_CD62L_in_CD3CD4",
    "marker_medians", "median_CD27_in_CD3CD4",
    "marker_medians", "median_CCR6_in_CD3CD4",
    "marker_medians", "median_IL7RA_in_CD3CD4",
    "marker_medians", "median_CCR4_in_CD3CD4",
    "marker_medians", "median_CD25_in_CD3CD4",
    "marker_medians", "median_CXCR5_in_CD3CD4",
    "cd4_differentiation", "pct_naive_like",
    "cd4_differentiation", "pct_tcm_like",
    "cd4_differentiation", "pct_tem_like",
    "cd4_differentiation", "pct_temra_like",
    "cd4_differentiation", "pct_memory_like",
    "cd27_refined", "pct_cd27_pos",
    "cd27_refined", "pct_cd27_neg",
    "cd27_refined", "pct_naive_cd27pos",
    "cd27_refined", "pct_memory_cd27pos",
    "cd27_refined", "pct_memory_cd27neg",
    "cd27_refined", "pct_cd62lneg_cd27neg",
    "single_helper_markers", "pct_ccr6_pos",
    "single_helper_markers", "pct_il7ra_pos",
    "single_helper_markers", "pct_il7ra_low",
    "single_helper_markers", "pct_ccr4_pos",
    "single_helper_markers", "pct_cd25_pos",
    "single_helper_markers", "pct_cd25_high",
    "single_helper_markers", "pct_cxcr5_pos",
    "treg_like", "pct_cd25pos_il7ralow_treg_like",
    "treg_like", "pct_cd25high_il7ralow_treg_enriched",
    "treg_like", "pct_ccr4_cd25pos_il7ralow",
    "treg_like", "pct_ccr4_cd25high_il7ralow",
    "tfh_th17_helper", "pct_cxcr5pos_tfh_like",
    "tfh_th17_helper", "pct_ccr6pos_th17_like",
    "tfh_th17_helper", "pct_ccr6pos_cxcr5pos",
    "tfh_th17_helper", "pct_ccr4pos_ccr6pos",
    "tfh_th17_helper", "pct_ccr4pos_cxcr5pos",
    "within_subset", "pct_cd25pos_within_cxcr5pos",
    "within_subset", "pct_ccr6pos_within_cxcr5pos",
    "within_subset", "pct_cxcr5pos_within_ccr6pos",
    "within_subset", "pct_ccr4pos_within_treg_like",
    "within_subset", "pct_cd25pos_il7ralow_within_memory"
  ),
  score_definitions = list(
    CP8_CD4_composition_depletion_score = list(
      positive = character(),
      negative = c(
        "pct_cd3_pos_total", "pct_cd4_pos_total",
        "pct_cd3_cd4_pos_total", "pct_cd4_within_cd3"
      )
    ),
    CP8_CD4_memory_differentiation_score = list(
      positive = c(
        "pct_memory_like", "pct_tem_like", "pct_temra_like", "pct_cd27_neg",
        "pct_memory_cd27neg", "pct_cd62lneg_cd27neg"
      ),
      negative = c(
        "pct_naive_like", "pct_naive_cd27pos",
        "median_CD45RA_in_CD3CD4", "median_CD62L_in_CD3CD4"
      )
    ),
    CP8_Treg_like_enrichment_score = list(
      positive = c(
        "pct_cd25pos_il7ralow_treg_like",
        "pct_cd25high_il7ralow_treg_enriched",
        "pct_ccr4_cd25pos_il7ralow",
        "pct_ccr4_cd25high_il7ralow",
        "pct_cd25pos_il7ralow_within_memory",
        "pct_ccr4pos_within_treg_like",
        "pct_il7ra_low", "pct_cd25_high", "median_CD25_in_CD3CD4"
      ),
      negative = c("pct_il7ra_pos", "median_IL7RA_in_CD3CD4")
    ),
    CP8_CCR6_CCR4_helper_remodeling_score = list(
      positive = c(
        "pct_ccr6_pos", "median_CCR6_in_CD3CD4",
        "pct_ccr4_pos", "median_CCR4_in_CD3CD4",
        "pct_ccr6pos_th17_like", "pct_ccr4pos_ccr6pos",
        "pct_ccr6pos_within_cxcr5pos"
      ),
      negative = character()
    ),
    CP8_Tfh_attenuation_score = list(
      positive = character(),
      negative = c(
        "pct_cxcr5_pos", "pct_cxcr5pos_tfh_like",
        "median_CXCR5_in_CD3CD4", "pct_ccr4pos_cxcr5pos",
        "pct_cxcr5pos_within_ccr6pos"
      )
    )
  ),
  integrated_score = "CP8_integrated_CD4_helper_regulatory_remodeling_score",

  # Compact subset used for selected figures and summaries.
  targeted_outcomes = c(
    "pct_cd3_cd4_pos_total", "pct_cd4_within_cd3", "pct_naive_like",
    "pct_memory_like", "pct_tem_like", "median_CD62L_in_CD3CD4",
    "pct_il7ra_low", "pct_cd25_high",
    "pct_cd25pos_il7ralow_treg_like",
    "pct_cd25high_il7ralow_treg_enriched",
    "pct_ccr4_cd25pos_il7ralow", "pct_ccr6_pos",
    "pct_ccr4_pos", "pct_cxcr5_pos"
  ),

  # Exact archived 31-outcome family used for age strata, matching, and
  # disease-by-age interaction analyses.
  age_sensitivity_outcomes = c(
    "CP8_integrated_CD4_helper_regulatory_remodeling_score",
    "CP8_CD4_composition_depletion_score",
    "CP8_CD4_memory_differentiation_score",
    "CP8_Treg_like_enrichment_score",
    "CP8_CCR6_CCR4_helper_remodeling_score",
    "CP8_Tfh_attenuation_score",
    "pct_cd3_pos_total", "pct_cd4_pos_total", "pct_cd3_cd4_pos_total",
    "pct_cd4_within_cd3", "pct_naive_like", "pct_memory_like",
    "pct_tem_like", "pct_cd62lneg_cd27neg", "pct_memory_cd27neg",
    "median_CD62L_in_CD3CD4",
    "pct_cd25pos_il7ralow_treg_like",
    "pct_cd25high_il7ralow_treg_enriched",
    "pct_ccr4_cd25pos_il7ralow", "pct_ccr4_cd25high_il7ralow",
    "pct_cd25pos_il7ralow_within_memory", "pct_il7ra_low",
    "pct_il7ra_pos", "pct_ccr6_pos", "pct_ccr4_pos", "pct_cxcr5_pos",
    "pct_ccr6pos_th17_like", "pct_cxcr5pos_tfh_like",
    "pct_ccr4pos_ccr6pos", "pct_ccr6pos_within_cxcr5pos",
    "pct_cxcr5pos_within_ccr6pos"
  ),

  # Exact archived 28-outcome threshold-sensitivity family. BH FDR is
  # calculated separately within permissive, main, and stringent sets.
  threshold_sensitivity_outcomes = c(
    "pct_cd3_pos_total", "pct_cd4_pos_total", "pct_cd3_cd4_pos_total",
    "pct_cd4_within_cd3", "pct_naive_like", "pct_memory_like",
    "pct_tem_like", "pct_cd62lneg_cd27neg", "median_CD62L_in_CD3CD4",
    "pct_cd25pos_il7ralow_treg_like",
    "pct_cd25high_il7ralow_treg_enriched",
    "pct_ccr4_cd25pos_il7ralow", "pct_ccr4_cd25high_il7ralow",
    "pct_cd25pos_il7ralow_within_memory", "pct_il7ra_low",
    "pct_il7ra_pos", "pct_ccr6_pos", "pct_ccr4_pos", "pct_cxcr5_pos",
    "pct_ccr6pos_th17_like", "pct_cxcr5pos_tfh_like",
    "pct_ccr4pos_ccr6pos",
    "CP8_integrated_CD4_helper_regulatory_remodeling_score",
    "CP8_CD4_composition_depletion_score",
    "CP8_CD4_memory_differentiation_score",
    "CP8_Treg_like_enrichment_score",
    "CP8_CCR6_CCR4_helper_remodeling_score",
    "CP8_Tfh_attenuation_score"
  ),

  benchmarks = list(
    n_files = 850,
    n_model = 832,
    n_healthy_model = 398,
    n_cancer_model = 434,
    main_feature_count = 44,
    age_outcome_count = 31,
    threshold_outcome_count = 28,
    median_total_events = 299655,
    median_primary_events = 26156,
    integrated_beta = 0.527791485673167,
    matched_pairs_5y = 265,
    matched_pairs_10y = 273
  )
)
