# CP7 reconstructed manifest for downstream sensitivity, annotation, and figures.
CP7_MANIFEST <- list(
 panel="CP7", primary_event_col="n_cd3_cd8_pos",
 thresholds=c(CD3=2,CD8=2.2,CD45RA=2,CD62L=2.2,CD27=2,PD1=1.5,PD1_HIGH=2,TIM3=1.5,LAG3=1.5,TIGIT=1.5,ICOS=1.5,CD39=1.5),
 score_definitions=list(
  CP7_CD8_differentiation_remodeling_score=list(positive=c("pct_temra_like","pct_temra_cd27pos","pct_temra_cd27neg","pct_cd62lneg_cd27neg","pct_cd45ra_pos_cd62lneg_cd27neg","pct_cd8_within_cd3"),negative=c("median_CD62L_in_CD3CD8","pct_naive_like","pct_tcm_like","pct_cd3_pos_total","pct_cd3_cd8_pos_total")),
  CP7_CD39_enrichment_score=list(positive=c("median_CD39_in_CD3CD8","pct_cd39_pos","pct_cd39_temra","pct_cd39_temra_cd27neg","pct_cd39_pos_within_temra_like","pct_cd39_pos_within_temra_cd27neg","pct_pd1_cd39_pos","pct_tigit_cd39_pos","pct_icos_cd39_pos"),negative=character()),
  CP7_PD1_TIGIT_attenuation_score=list(positive=character(),negative=c("pct_pd1_pos","pct_pd1_high","median_PD1_in_CD3CD8","pct_tigit_pos","median_TIGIT_in_CD3CD8","pct_pd1_tigit_pos","pct_tigit_pos_within_temra_like","pct_tigit_pos_within_temra_cd27neg","pct_pd1_pos_within_temra_like","pct_pd1_pos_within_temra_cd27neg")),
  CP7_checkpoint_coexpression_score=list(positive=c("pct_tim3_lag3_pos","pct_pd1_tim3_lag3_pos","pct_icos_cd39_pos","pct_pd1_tim3_pos","pct_pd1_cd39_pos","pct_pd1_tigit_cd39_pos","pct_pd1_tigit_icos_cd39_pos"),negative=character()),
  CP7_terminal_checkpoint_remodeling_score=list(positive=c("pct_cd39_temra","pct_cd39_temra_cd27neg","pct_cd39_pos_within_temra_like","pct_cd39_pos_within_temra_cd27neg","pct_lag3_pos_within_temra_like"),negative=c("pct_tigit_pos_within_temra_like","pct_tigit_pos_within_temra_cd27neg","pct_pd1_pos_within_temra_like","pct_pd1_pos_within_temra_cd27neg"))
 ),
 integrated_score="CP7_integrated_checkpoint_remodeling_score",
 targeted_outcomes=c("pct_cd3_pos_total","pct_cd8_pos_total","pct_cd3_cd8_pos_total","pct_cd8_within_cd3","median_CD62L_in_CD3CD8","pct_naive_like","pct_tcm_like","pct_temra_like","pct_temra_cd27pos","pct_temra_cd27neg","pct_cd39_pos","median_CD39_in_CD3CD8","pct_tim3_lag3_pos","pct_pd1_pos","pct_tigit_pos","pct_pd1_tigit_pos","pct_cd39_temra","pct_cd39_temra_cd27neg","pct_cd39_pos_within_temra_like","pct_cd39_pos_within_temra_cd27neg","pct_tigit_pos_within_temra_like","pct_pd1_pos_within_temra_like"),
 benchmarks=list(n_files=850,n_model=832,matched_pairs_5y=302,matched_pairs_10y=332,integrated_beta=.351119363746345)
)
