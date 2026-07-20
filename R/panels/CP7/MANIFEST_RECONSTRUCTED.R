# CP7 reconstructed manifest for downstream sensitivity, annotation, and figures.
# All feature lists below are tied to the archived CP7 output tables.
CP7_MANIFEST <- list(
  panel = "CP7",
  primary_event_col = "n_cd3_cd8_pos",
  thresholds = c(
    CD3 = 2.0, CD8 = 2.2, CD45RA = 2.0, CD62L = 2.2, CD27 = 2.0,
    PD1 = 1.5, PD1_HIGH = 2.0, TIM3 = 1.5, LAG3 = 1.5,
    TIGIT = 1.5, ICOS = 1.5, CD39 = 1.5
  ),
  score_definitions = list(
    CP7_CD8_differentiation_remodeling_score = list(
      positive = c(
        "pct_temra_like", "pct_temra_cd27pos", "pct_temra_cd27neg",
        "pct_cd62lneg_cd27neg", "pct_cd45ra_pos_cd62lneg_cd27neg",
        "pct_cd8_within_cd3"
      ),
      negative = c(
        "median_CD62L_in_CD3CD8", "pct_naive_like", "pct_tcm_like",
        "pct_cd3_pos_total", "pct_cd3_cd8_pos_total"
      )
    ),
    CP7_CD39_enrichment_score = list(
      positive = c(
        "median_CD39_in_CD3CD8", "pct_cd39_pos", "pct_cd39_temra",
        "pct_cd39_temra_cd27neg", "pct_cd39_pos_within_temra_like",
        "pct_cd39_pos_within_temra_cd27neg", "pct_pd1_cd39_pos",
        "pct_tigit_cd39_pos", "pct_icos_cd39_pos"
      ),
      negative = character()
    ),
    CP7_PD1_TIGIT_attenuation_score = list(
      positive = character(),
      negative = c(
        "pct_pd1_pos", "pct_pd1_high", "median_PD1_in_CD3CD8",
        "pct_tigit_pos", "median_TIGIT_in_CD3CD8", "pct_pd1_tigit_pos",
        "pct_tigit_pos_within_temra_like",
        "pct_tigit_pos_within_temra_cd27neg",
        "pct_pd1_pos_within_temra_like",
        "pct_pd1_pos_within_temra_cd27neg"
      )
    ),
    CP7_checkpoint_coexpression_score = list(
      positive = c(
        "pct_tim3_lag3_pos", "pct_pd1_tim3_lag3_pos",
        "pct_icos_cd39_pos", "pct_pd1_tim3_pos", "pct_pd1_cd39_pos",
        "pct_pd1_tigit_cd39_pos", "pct_pd1_tigit_icos_cd39_pos"
      ),
      negative = character()
    ),
    CP7_terminal_checkpoint_remodeling_score = list(
      positive = c(
        "pct_cd39_temra", "pct_cd39_temra_cd27neg",
        "pct_cd39_pos_within_temra_like",
        "pct_cd39_pos_within_temra_cd27neg",
        "pct_lag3_pos_within_temra_like"
      ),
      negative = c(
        "pct_tigit_pos_within_temra_like",
        "pct_tigit_pos_within_temra_cd27neg",
        "pct_pd1_pos_within_temra_like",
        "pct_pd1_pos_within_temra_cd27neg"
      )
    )
  ),
  integrated_score = "CP7_integrated_checkpoint_remodeling_score",

  # Exact 32-outcome family used in the archived age-stratified, matched,
  # and disease-by-age interaction analyses.
  age_sensitivity_modules = c(
    pct_cd3_pos_total = "composition",
    pct_cd8_pos_total = "composition",
    pct_cd3_cd8_pos_total = "composition",
    pct_cd8_within_cd3 = "composition",
    median_CD62L_in_CD3CD8 = "differentiation",
    pct_temra_like = "differentiation",
    pct_temra_cd27pos = "differentiation",
    median_TIGIT_in_CD3CD8 = "checkpoint_medians",
    median_CD39_in_CD3CD8 = "checkpoint_medians",
    pct_pd1_pos = "checkpoint_single_positive",
    pct_pd1_high = "checkpoint_single_positive",
    pct_tigit_pos = "checkpoint_single_positive",
    pct_cd39_pos = "checkpoint_single_positive",
    pct_pd1_tigit_pos = "checkpoint_coexpression",
    pct_tim3_lag3_pos = "checkpoint_coexpression",
    pct_icos_cd39_pos = "checkpoint_coexpression",
    pct_pd1_tim3_lag3_pos = "checkpoint_coexpression",
    pct_pd1_tigit_cd39_temra_cd27neg = "checkpoint_coexpression",
    pct_pd1_pos_within_temra_like = "terminal_checkpoint",
    pct_lag3_pos_within_temra_like = "terminal_checkpoint",
    pct_tigit_pos_within_temra_like = "terminal_checkpoint",
    pct_cd39_pos_within_temra_like = "terminal_checkpoint",
    pct_tigit_pos_within_temra_cd27neg = "terminal_checkpoint",
    pct_cd39_pos_within_temra_cd27neg = "terminal_checkpoint",
    pct_cd39_temra = "terminal_checkpoint",
    pct_cd39_temra_cd27neg = "terminal_checkpoint",
    CP7_integrated_checkpoint_remodeling_score = "composite_score",
    CP7_terminal_checkpoint_remodeling_score = "composite_score",
    CP7_CD8_differentiation_remodeling_score = "composite_score",
    CP7_CD39_enrichment_score = "composite_score",
    CP7_PD1_TIGIT_attenuation_score = "composite_score",
    CP7_checkpoint_coexpression_score = "composite_score"
  ),

  # Exact 11-outcome family used for the targeted threshold-sensitivity FDR.
  threshold_targeted_outcomes = c(
    "pct_cd3_pos_total", "pct_cd3_cd8_pos_total", "pct_cd8_within_cd3",
    "pct_temra_like", "pct_cd39_pos", "pct_cd39_temra",
    "pct_tim3_lag3_pos", "CP7_PD1_TIGIT_attenuation_score",
    "CP7_CD39_enrichment_score",
    "CP7_terminal_checkpoint_remodeling_score",
    "CP7_integrated_checkpoint_remodeling_score"
  ),

  # Retained for quantitative figures and clinical annotation; this is not the
  # threshold-testing family and must not determine targeted FDR denominators.
  targeted_outcomes = c(
    "pct_cd3_pos_total", "pct_cd8_pos_total", "pct_cd3_cd8_pos_total",
    "pct_cd8_within_cd3", "median_CD62L_in_CD3CD8", "pct_naive_like",
    "pct_tcm_like", "pct_temra_like", "pct_temra_cd27pos",
    "pct_temra_cd27neg", "pct_cd39_pos", "median_CD39_in_CD3CD8",
    "pct_tim3_lag3_pos", "pct_pd1_pos", "pct_tigit_pos",
    "pct_pd1_tigit_pos", "pct_cd39_temra", "pct_cd39_temra_cd27neg",
    "pct_cd39_pos_within_temra_like",
    "pct_cd39_pos_within_temra_cd27neg",
    "pct_tigit_pos_within_temra_like", "pct_pd1_pos_within_temra_like"
  ),
  benchmarks = list(
    n_files = 850,
    n_model = 832,
    matched_pairs_5y = 302,
    matched_pairs_10y = 332,
    integrated_beta = 0.351119363746345,
    threshold_target_count = 11
  )
)
