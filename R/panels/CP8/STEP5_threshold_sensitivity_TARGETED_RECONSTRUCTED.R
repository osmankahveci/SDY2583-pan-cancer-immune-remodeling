# CP8 Step 5 reconstructed targeted threshold sensitivity.
rm(list=ls())
source(file.path(Sys.getenv("SDY2583_REPO_ROOT", unset="."),"R","shared","reconstructed_panel_framework.R"))
rp_install_and_load(c("dplyr","readr","tibble","purrr","tidyr"),"flowCore")
source(file.path(sd_repo_root(),"R","panels","CP8","MANIFEST_RECONSTRUCTED.R"))
analysis_dir<-sd_analysis_dir("CP8"); rdata_dir<-file.path(analysis_dir,"11_RData"); out_dir<-file.path(analysis_dir,"12_threshold_sensitivity"); dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
load(file.path(rdata_dir,"SDY2583_CP8_STEP2_feature_extraction_RECONSTRUCTED.RData"))
load(file.path(rdata_dir,"SDY2583_CP8_STEP3A_metadata_merge_age_QC_RECONSTRUCTED.RData"))
metadata_for_merge<-metadata
main<-CP8_MANIFEST$thresholds
permissive<-main-0.2; stringent<-main+0.2
# IL7RA-low is the primary low-gate phenotype, so its permissive/stringent
# direction is reversed relative to positive gates, matching the archive.
permissive["IL7RA"]<-main["IL7RA"]+0.2; stringent["IL7RA"]<-main["IL7RA"]-0.2
sets<-list(permissive=permissive,main=main,stringent=stringent)
threshold_sets<-purrr::imap_dfr(sets,~tibble::tibble(threshold_set=.y,marker=names(.x),threshold=as.numeric(.x)))

one_set<-function(label,thresholds){
  ft<-purrr::map_dfr(fcs_files,~extract_one(.x,thresholds))
  dat<-ft |> dplyr::left_join(metadata_for_merge,by="subject_id")
  scored<-rp_build_scores(dat,CP8_MANIFEST$score_definitions,CP8_MANIFEST$integrated_score)
  outcomes<-unique(c(CP8_MANIFEST$targeted_outcomes,CP8_MANIFEST$integrated_score,names(CP8_MANIFEST$score_definitions)))
  stats<-purrr::map_dfr(outcomes,~rp_fit_one(scored,.x,label,FALSE,CP8_MANIFEST$primary_event_col,0)) |> dplyr::mutate(targeted_fdr=p.adjust(p_value,"BH"),threshold_set=label)
  list(features=ft,data=scored,statistics=stats)
}
results<-purrr::imap(sets,one_set)
statistics_long<-dplyr::bind_rows(lapply(results,`[[`,"statistics"))
selected_features_long<-dplyr::bind_rows(purrr::imap(results,~.x$data |> dplyr::select(subject_id,dplyr::all_of(unique(c(CP8_MANIFEST$targeted_outcomes,CP8_MANIFEST$integrated_score,names(CP8_MANIFEST$score_definitions))))) |> dplyr::mutate(threshold_set=.y)))
summary_wide<-statistics_long |> dplyr::select(threshold_set,feature,beta_cancer_vs_healthy,p_value,targeted_fdr,direction) |> tidyr::pivot_wider(names_from=threshold_set,values_from=c(beta_cancer_vs_healthy,p_value,targeted_fdr,direction)) |> dplyr::mutate(direction_preserved=direction_main==direction_permissive & direction_main==direction_stringent,fdr_preserved=targeted_fdr_main<.05 & targeted_fdr_permissive<.05 & targeted_fdr_stringent<.05,robust_class=dplyr::case_when(direction_preserved & fdr_preserved~"direction_and_FDR_preserved",direction_preserved~"direction_preserved_FDR_not_preserved",TRUE~"direction_not_preserved"))
readr::write_csv(threshold_sets,file.path(out_dir,"SDY2583_CP8_threshold_sensitivity_threshold_sets_RECONSTRUCTED.csv")); readr::write_csv(selected_features_long,file.path(out_dir,"SDY2583_CP8_threshold_sensitivity_selected_features_LONG_RECONSTRUCTED.csv")); readr::write_csv(statistics_long,file.path(out_dir,"SDY2583_CP8_threshold_sensitivity_statistics_LONG_RECONSTRUCTED.csv")); readr::write_csv(summary_wide,file.path(out_dir,"SDY2583_CP8_threshold_sensitivity_summary_WIDE_RECONSTRUCTED.csv"))
save(results,statistics_long,summary_wide,threshold_sets,file=file.path(rdata_dir,"SDY2583_CP8_STEP5_threshold_sensitivity_RECONSTRUCTED.RData")); print(table(summary_wide$robust_class))
