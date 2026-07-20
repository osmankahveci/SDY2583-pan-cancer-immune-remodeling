# Validate reconstructed CP8 outputs against fixed archive benchmarks and,
# when SDY2583_CP8_REFERENCE_DIR is set, archived CSV tables.
source(file.path(Sys.getenv("SDY2583_REPO_ROOT",unset="."),"R","shared","reconstructed_panel_framework.R"))
rp_install_and_load(c("dplyr","readr","tibble","purrr"))
source(file.path(sd_repo_root(),"R","panels","CP8","MANIFEST_RECONSTRUCTED.R"))
a<-sd_analysis_dir("CP8"); out<-file.path(a,"validation"); dir.create(out,recursive=TRUE,showWarnings=FALSE)
step2<-readr::read_csv(file.path(a,"02_feature_extraction","SDY2583_CP8_feature_extraction_summary_STEP2_RECONSTRUCTED.csv"),show_col_types=FALSE)
meta<-readr::read_csv(file.path(a,"04_metadata_merge","SDY2583_CP8_metadata_merge_summary_RECONSTRUCTED.csv"),show_col_types=FALSE)
stats<-readr::read_csv(file.path(a,"06_statistics","SDY2583_CP8_FULL850_age_sex_adjusted_statistics_MAIN_RECONSTRUCTED.csv"),show_col_types=FALSE)
scores<-readr::read_csv(file.path(a,"08_composite_scores","SDY2583_CP8_composite_score_statistics_RECONSTRUCTED.csv"),show_col_types=FALSE)
match<-readr::read_csv(file.path(a,"10_age_sensitivity","SDY2583_CP8_age_matched_sensitivity_summary_RECONSTRUCTED.csv"),show_col_types=FALSE)
near<-function(x,y,tol=.02) is.finite(x)&&abs(x-y)<=tol*max(1,abs(y))
checks<-tibble::tribble(~check,~observed,~expected,~pass,
"files",step2$n_files[1],850,step2$n_files[1]==850,
"successful extractions",step2$n_success[1],850,step2$n_success[1]==850,
"median total events",step2$median_total_events[1],299655,near(step2$median_total_events[1],299655,.001),
"median primary events",step2$median_cd3_cd4_events[1],26156,near(step2$median_cd3_cd4_events[1],26156,.001),
"valid ages",meta$n_valid_age[1],832,meta$n_valid_age[1]==832,
"feature count",nrow(stats),45,nrow(stats)==45,
"5y matched pairs",match$n_pairs[match$caliper_years==5][1],265,match$n_pairs[match$caliper_years==5][1]==265,
"10y matched pairs",match$n_pairs[match$caliper_years==10][1],273,match$n_pairs[match$caliper_years==10][1]==273,
"integrated beta",scores$beta_cancer_vs_healthy[scores$score==CP8_MANIFEST$integrated_score][1],.528,near(scores$beta_cancer_vs_healthy[scores$score==CP8_MANIFEST$integrated_score][1],.528,.01))
ref<-path.expand(Sys.getenv("SDY2583_CP8_REFERENCE_DIR",unset=""))
if(nzchar(ref)&&dir.exists(ref)){
  pairs<-list(
    c(file.path(a,"02_feature_extraction","SDY2583_CP8_FULL_850_feature_table_STEP2_RECONSTRUCTED.csv"),"SDY2583_CP8_FULL_850_feature_table_STEP2.csv","subject_id"),
    c(file.path(a,"06_statistics","SDY2583_CP8_FULL850_age_sex_adjusted_statistics_MAIN_RECONSTRUCTED.csv"),"SDY2583_CP8_FULL850_age_sex_adjusted_statistics_MAIN.csv","module,feature"),
    c(file.path(a,"08_composite_scores","SDY2583_CP8_composite_score_statistics_RECONSTRUCTED.csv"),"SDY2583_CP8_composite_score_statistics.csv","score")
  )
  for(p in pairs){hit<-list.files(ref,p[2],full.names=TRUE,recursive=TRUE,fixed=TRUE); cmp<-if(length(hit)) rp_compare_csv(p[1],hit[1],strsplit(p[3],",",fixed=TRUE)[[1]]) else tibble::tibble(pass=FALSE,detail="reference missing"); checks<-dplyr::bind_rows(checks,tibble::tibble(check=p[2],observed=cmp$detail[1],expected="table match",pass=cmp$pass[1]))}
}
readr::write_csv(checks,file.path(out,"SDY2583_CP8_reconstruction_validation_report.csv")); print(checks,n=nrow(checks)); if(any(!checks$pass)) stop("CP8 validation failed.")
