# CP26 reconstructed NK-panel manifest.
CP26_MANIFEST <- list(
 panel="CP26", primary_event_col="n_nk_like",
 thresholds=c(DUMP_LOW=1.5,CD45=2,CD56=1.5,CD56_HIGH=2.4,CD16=1.5,CD16_HIGH=2.2,NKG2A=1.5,NKG2C=1.5,NKG2D=1.5,CD158=1.5,CD57=1.5,NKp44=1.5,CD161=1.5,CD107a=1.5),
 expected_markers=c("Viability_CD3_CD13_CD19_CD123","CD45","CD56","CD16","NKG2A","NKG2C","NKG2D","CD158","CD57","NKp44","CD161","CD107a"),
 fallback_channels=list("Viability_CD3_CD13_CD19_CD123"="BV510-A"),
 module_map=tibble::tribble(
 ~feature,~module,
 "pct_nk_like_total","composition_NK",
 "pct_nk_like_within_cd45_dump_low","composition_NK",
 "pct_cd56pos_within_nk_like","cd56_cd16_balance",
 "pct_cd16pos_within_nk_like","cd56_cd16_balance",
 "pct_cd56_cd16_double_pos_within_nk_like","cd56_cd16_balance",
 "pct_cd56bright_like_within_nk_like","cd56_cd16_balance",
 "pct_cd56pos_cd16high_like_within_nk_like","cd56_cd16_balance",
 "pct_nkg2a_pos_within_nk_like","receptor_phenotype_NK",
 "pct_nkg2c_pos_within_nk_like","receptor_phenotype_NK",
 "pct_nkg2d_pos_within_nk_like","receptor_phenotype_NK",
 "pct_cd158_pos_within_nk_like","receptor_phenotype_NK",
 "pct_nkp44_pos_within_nk_like","receptor_phenotype_NK",
 "pct_cd161_pos_within_nk_like","receptor_phenotype_NK",
 "pct_cd57_pos_within_nk_like","activation_maturation_NK",
 "pct_cd107a_pos_within_nk_like","activation_maturation_NK"),
 score_definitions=list(
  CP26_NK_like_composition_depletion_score=list(positive=character(),negative=c("pct_nk_like_total","pct_nk_like_within_cd45_dump_low")),
  CP26_CD56_CD16_mature_NK_depletion_score=list(positive=character(),negative=c("pct_cd56pos_within_nk_like","pct_cd56_cd16_double_pos_within_nk_like","pct_cd56pos_cd16high_like_within_nk_like")),
  CP26_NK_receptor_maturation_attenuation_score=list(positive=character(),negative=c("pct_nkg2d_pos_within_nk_like","pct_cd161_pos_within_nk_like","pct_cd57_pos_within_nk_like","pct_cd158_pos_within_nk_like")),
  CP26_NKG2C_adaptive_like_shift_score=list(positive=c("pct_nkg2c_pos_within_nk_like"),negative=character())
 ),
 integrated_score="CP26_integrated_NK_remodeling_score",
 targeted_outcomes=c("pct_nk_like_total","pct_nk_like_within_cd45_dump_low","pct_cd56pos_within_nk_like","pct_cd16pos_within_nk_like","pct_cd56_cd16_double_pos_within_nk_like","pct_cd56bright_like_within_nk_like","pct_cd56pos_cd16high_like_within_nk_like","pct_nkg2a_pos_within_nk_like","pct_nkg2c_pos_within_nk_like","pct_nkg2d_pos_within_nk_like","pct_cd158_pos_within_nk_like","pct_nkp44_pos_within_nk_like","pct_cd161_pos_within_nk_like","pct_cd57_pos_within_nk_like","pct_cd107a_pos_within_nk_like"),
 benchmarks=list(n_files=850,n_model=832,n_healthy=398,n_cancer=434,median_total_events=299279.5,median_cd45_dump_low=3839,median_nk_like=3505,event_qc=300,matched_pairs_5y=263,matched_pairs_10y=272,integrated_beta=.364755491769086,threshold_target_n=20,robust_both=10,direction_only=7,direction_fail=3)
)
