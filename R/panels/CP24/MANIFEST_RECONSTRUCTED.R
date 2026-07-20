# CP24 full-cohort reconstructed manifest. The recovered codebook remains a
# separate pilot/QC artifact; this manifest describes the archived 850-subject
# final analysis reconstructed from Methods and final output schemas.
CP24_MODULES<-list(
 composition=c("pct_cd3_pos_total","pct_cd8_pos_total","pct_cd3_cd8_pos_total","pct_cd8_within_cd3"),
 marker_medians=c("median_CD62L_in_CD3CD8","median_CD27_in_CD3CD8","median_PD1_in_CD3CD8","median_CD57_in_CD3CD8","median_CX3CR1_in_CD3CD8","median_CD95_in_CD3CD8","median_CD45RA_in_CD3CD8","median_CXCR3_in_CD3CD8","median_CXCR5_in_CD3CD8"),
 cd45ra_cd62l_quadrants=c("pct_naive_like","pct_tcm_like","pct_tem_like","pct_temra_like","pct_cd62l_pos_within_cd3cd8","pct_cd62l_neg_within_cd3cd8","pct_cd45ra_pos_within_cd3cd8","pct_cd45ra_neg_within_cd3cd8"),
 cd27_refined_temra=c("pct_cd27_pos_within_cd3cd8","pct_cd27_neg_within_cd3cd8","pct_temra_cd27neg","pct_temra_cd27pos","pct_cd62lneg_cd27neg","pct_cd45ra_pos_cd62lneg_cd27neg","pct_naive_cd27pos","pct_naive_cd27neg"),
 extended_differentiation=c("pct_cd57_pos","pct_cx3cr1_pos","pct_cd95_pos","pct_cd57_cx3cr1_pos","pct_cd57_cd95_pos","pct_cx3cr1_cd95_pos","pct_cd57_cx3cr1_cd95_pos"),
 pd1_module=c("pct_pd1_pos","pct_pd1_high","pct_pd1_temra","pct_pd1_temra_cd27neg","pct_pd1_cd57_pos","pct_pd1_cx3cr1_pos","pct_pd1_cd95_pos","pct_pd1_cd57_cx3cr1_cd95_pos","pct_pd1_pos_within_temra","pct_pd1_pos_within_temra_cd27neg")
)
CP24_MODULE_MAP<-dplyr::bind_rows(lapply(names(CP24_MODULES),function(m)tibble::tibble(module=m,feature=CP24_MODULES[[m]])))
CP24_MANIFEST<-list(
 panel="CP24",primary_event_col="n_cd3_cd8_pos",
 thresholds=c(CD3=2,CD8=2.2,CD45RA=2,CD62L=2.2,CD27=2,CD57=1.5,CX3CR1=1.5,CD95=1.5,PD1=1.5,PD1_HIGH=2,CXCR3=1.5,CXCR5=1.5),
 expected_markers=c("CD3","CD8","CD45RA","CD62L","CD27","CD57","CX3CR1","CD95","PD-1","CXCR3","CXCR5"),module_map=CP24_MODULE_MAP,
 score_definitions=list(
  CP24_marker_loss_score=list(positive=character(),negative=c("median_CD62L_in_CD3CD8","median_CD27_in_CD3CD8")),
  CP24_quadrant_shift_score=list(positive=c("pct_temra_like","pct_cd62l_neg_within_cd3cd8"),negative=c("pct_naive_like","pct_cd62l_pos_within_cd3cd8")),
  CP24_CD27_refined_terminal_score=list(positive=c("pct_cd27_neg_within_cd3cd8","pct_temra_cd27neg","pct_cd62lneg_cd27neg","pct_cd45ra_pos_cd62lneg_cd27neg"),negative=c("pct_cd27_pos_within_cd3cd8","pct_naive_cd27pos")),
  CP24_extended_effector_score=list(positive=c("pct_cd57_pos","pct_cx3cr1_pos","pct_cd95_pos","pct_cd57_cx3cr1_pos","pct_cd57_cd95_pos","pct_cx3cr1_cd95_pos","pct_cd57_cx3cr1_cd95_pos"),negative=character()),
  CP24_PD1_associated_terminal_score=list(positive=c("pct_pd1_pos","pct_pd1_high","pct_pd1_temra","pct_pd1_temra_cd27neg","pct_pd1_pos_within_temra","pct_pd1_pos_within_temra_cd27neg","pct_pd1_cd57_cx3cr1_cd95_pos"),negative=character()),
  CP24_core_differentiation_score=list(positive=c("pct_temra_like","pct_cd27_neg_within_cd3cd8","pct_temra_cd27neg","pct_cd62lneg_cd27neg","pct_cd57_cx3cr1_cd95_pos"),negative=c("median_CD62L_in_CD3CD8","median_CD27_in_CD3CD8","pct_naive_like"))
 ),integrated_score="CP24_integrated_remodeling_score",
 targeted_outcomes=unique(unlist(CP24_MODULES)),
 benchmarks=list(n_files=850,n_model=832,n_healthy=398,n_cancer=434,median_total_events=302822,median_primary_events=11836,event_n_all=850,event_n_500=849,event_n_1000=844,integrated_beta=.11,core_beta=.17)
)
