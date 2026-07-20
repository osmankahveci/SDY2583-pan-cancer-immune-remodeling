# CP28 reconstructed T/NK-interface manifest.
CP28_MODULE_FEATURES <- list(
 composition_T_NK_interface=c("pct_dump_low_total","pct_cd3_pos_total","pct_cd3_pos_within_dump_low","pct_cd3_neg_within_dump_low","pct_cd8_pos_total","pct_cd8_pos_within_dump_low","pct_cd3_cd8_pos_total","pct_cd3_cd8_pos_within_dump_low","pct_cd8_within_cd3","pct_cd3_within_cd8","pct_cd3neg_cd56pos_nk_like_total","pct_cd3neg_cd56pos_nk_like_within_dump_low","pct_cd3neg_cd56pos_nk_like_within_cd3neg"),
 cd8_differentiation_CD45RA_CD27=c("median_CD45RA_in_CD3CD8","median_CD27_in_CD3CD8","pct_cd8_naive_like","pct_cd8_memory_like","pct_cd8_temra_like","pct_cd8_cd45ra_pos","pct_cd8_cd27_pos","pct_cd8_cd27_neg"),
 cd8_effector_NK_like_axis=c("median_CD57_in_CD3CD8","median_CD56_in_CD3CD8","median_CD161_in_CD3CD8","pct_cd57_pos_within_cd8","pct_cd56_pos_within_cd8","pct_cd161_pos_within_cd8","pct_cd57_cd27neg_within_cd8","pct_cd57_temra_within_cd8","pct_cd56_cd57_within_cd8","pct_cd56_temra_within_cd8","pct_cd57pos_within_cd8_temra","pct_cd56pos_within_cd8_temra"),
 cd3_innate_like_T_cell_axis=c("median_CD8_in_CD3","median_CD56_in_CD3","median_CD57_in_CD3","median_CD161_in_CD3","median_Vdelta2_in_CD3","median_Va24Ja11_in_CD3","median_Valpha7_in_CD3","pct_cd8_pos_within_cd3","pct_cd56_pos_within_cd3","pct_cd57_pos_within_cd3","pct_cd161_pos_within_cd3","pct_vdelta2_pos_within_cd3","pct_va24ja11_pos_within_cd3","pct_valpha7_pos_within_cd3","pct_valpha7_cd161_pos_within_cd3","pct_cd56_cd57_pos_within_cd3","pct_cd8_cd57_pos_within_cd3","pct_cd8_cd56_pos_within_cd3","pct_cd8_cd161_pos_within_cd3","pct_vdelta2_cd8_pos_within_cd3","pct_va24ja11_cd8_pos_within_cd3","pct_valpha7_cd161_cd8_pos_within_cd3","pct_cd8pos_within_vdelta2_cd3","pct_cd8pos_within_va24ja11_cd3","pct_cd161pos_within_valpha7_cd3","pct_valpha7pos_within_cd161_cd3"),
 cd3neg_CD56pos_NK_like_axis=c("pct_cd3neg_cd56pos_nk_like_total","pct_cd3neg_cd56pos_nk_like_within_dump_low","pct_cd3neg_cd56pos_nk_like_within_cd3neg","pct_cd57pos_within_cd3neg_cd56pos_nk_like","pct_cd8pos_within_cd3neg_cd56pos_nk_like","pct_cd57_cd8_pos_within_cd3neg_cd56pos_nk_like")
)
CP28_MODULE_MAP <- dplyr::bind_rows(lapply(names(CP28_MODULE_FEATURES),function(m)tibble::tibble(module=m,feature=CP28_MODULE_FEATURES[[m]])))
CP28_MANIFEST <- list(
 panel="CP28", primary_event_col="n_cd3_cd8_pos",
 thresholds=c(DUMP_LOW=1.5,CD3=2,CD8=2,CD45RA=2,CD27=2,CD57=1.5,CD56=1.5,CD161=1.5,VDELTA2=1.5,VA24_JA11_TCR=1.5,VALPHA7=1.5),
 expected_markers=c("DUMP","CD3","CD8","CD45RA","CD27","CD57","CD56","CD161","Vdelta2","Va24-Ja11 TCR","Valpha7"),
 fallback_channels=list(DUMP="BV510-A"), module_map=CP28_MODULE_MAP,
 score_definitions=list(
  CP28_T_cell_composition_depletion_score=list(positive=c("pct_cd3_neg_within_dump_low"),negative=c("pct_cd3_pos_total","pct_cd3_pos_within_dump_low","pct_cd8_pos_total","pct_cd3_cd8_pos_total","pct_cd3_cd8_pos_within_dump_low")),
  CP28_CD8_enrichment_within_T_score=list(positive=c("median_CD8_in_CD3","pct_cd8_within_cd3","pct_cd8_pos_within_cd3"),negative=character()),
  CP28_CD8_CD45RA_TEMRA_differentiation_score=list(positive=c("median_CD45RA_in_CD3CD8","pct_cd8_cd45ra_pos","pct_cd8_temra_like","pct_cd8_cd27_neg"),negative=c("median_CD27_in_CD3CD8","pct_cd8_cd27_pos","pct_cd8_memory_like")),
  CP28_CD8_CD57_effector_like_score=list(positive=c("median_CD57_in_CD3","median_CD57_in_CD3CD8","pct_cd57_pos_within_cd3","pct_cd8_cd57_pos_within_cd3","pct_cd57_pos_within_cd8","pct_cd57_cd27neg_within_cd8","pct_cd57_temra_within_cd8"),negative=character()),
  CP28_CD8_CD161_attenuation_score=list(positive=character(),negative=c("median_CD161_in_CD3","median_CD161_in_CD3CD8","pct_cd161_pos_within_cd3","pct_cd161_pos_within_cd8","pct_cd161pos_within_valpha7_cd3","pct_valpha7pos_within_cd161_cd3","pct_cd161pos_within_valpha7_cd8","pct_valpha7pos_within_cd161_cd8")),
  CP28_Va24_Valpha7_CD8_innate_like_score=list(positive=c("median_Valpha7_in_CD3","pct_valpha7_pos_within_cd3","pct_cd8pos_within_va24ja11_cd3","pct_va24ja11_cd8_pos_within_cd3","pct_valpha7_cd161_cd8_pos_within_cd3","pct_valpha7_cd161_pos_within_cd8"),negative=character()),
  CP28_CD56_NK_like_T_score=list(positive=c("median_CD56_in_CD3","median_CD56_in_CD3CD8","pct_cd56_pos_within_cd3","pct_cd8_cd56_pos_within_cd3","pct_cd56_pos_within_cd8","pct_cd56_cd57_pos_within_cd3","pct_cd56_cd57_within_cd8","pct_cd56_temra_within_cd8"),negative=character()),
  CP28_CD3neg_CD56pos_NK_like_depletion_score=list(positive=character(),negative=c("pct_cd3neg_cd56pos_nk_like_total","pct_cd3neg_cd56pos_nk_like_within_dump_low","pct_cd3neg_cd56pos_nk_like_within_cd3neg"))
 ),
 integrated_score="CP28_integrated_TNK_interface_remodeling_score",
 targeted_outcomes=unique(unlist(CP28_MODULE_FEATURES)),
 event_qc_col="n_cd3_cd8_pos",event_qc=1000,
 benchmarks=list(n_files=850,n_model=832,median_total_events=309700.5,median_dump_low=77505,median_cd3_dump_low=25096,median_cd3cd8=7440,median_nk_like=3066.5,matched_pairs_5y=265,matched_pairs_10y=273,integrated_beta=.257)
)
