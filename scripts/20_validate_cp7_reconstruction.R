# Validate the complete reconstructed CP7 pipeline against fixed archive
# benchmarks and optional local archived tables.
source(file.path(Sys.getenv("SDY2583_REPO_ROOT",unset="."),"R","shared","reconstructed_panel_framework.R")); rp_install_and_load(c("dplyr","readr","tibble","purrr","stringr")); source(file.path(sd_repo_root(),"R","panels","CP7","MANIFEST_RECONSTRUCTED.R"))
a<-sd_analysis_dir("CP7"); out<-file.path(a,"validation"); dir.create(out,recursive=TRUE,showWarnings=FALSE); strict<-tolower(Sys.getenv("SDY2583_VALIDATION_STRICT",unset="true")) %in% c("true","1","yes"); tol<-as.numeric(Sys.getenv("SDY2583_VALIDATION_TOLERANCE",unset="1e-6")); if(!is.finite(tol))tol<-1e-6
read_req<-function(p){if(!file.exists(p))stop("Missing reconstructed output: ",p);readr::read_csv(p,show_col_types=FALSE)}; near<-function(x,y,t=tol)is.finite(x)&&abs(x-y)<=t*max(1,abs(y))
step1<-read_req(file.path(a,"01_channel_marker_QC","SDY2583_CP7_STEP1_QC_summary.csv")); step2<-read_req(file.path(a,"02_feature_extraction","SDY2583_CP7_FULL_850_STEP2_QC_summary_RECONSTRUCTED.csv")); meta<-read_req(file.path(a,"04_metadata_merge","SDY2583_CP7_metadata_merge_summary_RECONSTRUCTED.csv")); st<-read_req(file.path(a,"06_statistics","SDY2583_CP7_FULL850_age_sex_adjusted_statistics_MAIN_RECONSTRUCTED.csv")); sc<-read_req(file.path(a,"08_composite_scores","SDY2583_CP7_composite_score_statistics_RECONSTRUCTED.csv")); mt<-read_req(file.path(a,"10_age_sensitivity","SDY2583_CP7_age_matched_sensitivity_summary_RECONSTRUCTED.csv")); ts<-read_req(file.path(a,"12_threshold_sensitivity","SDY2583_CP7_threshold_sensitivity_summary_WIDE_RECONSTRUCTED.csv"))
get_beta<-function(feature)st$beta_cancer_vs_healthy[st$feature==feature][1]
checks<-tibble::tribble(~check,~observed,~expected,~pass,
"FCS files",step1$n_fcs_files[1],850,step1$n_fcs_files[1]==850,
"unique subjects",step1$n_unique_subjects[1],850,step1$n_unique_subjects[1]==850,
"successful extractions",step2$feature_ok[1],850,step2$feature_ok[1]==850,
"failed extractions",step2$feature_failed[1],0,step2$feature_failed[1]==0,
"median total events",step2$median_total_events[1],301989.5,near(step2$median_total_events[1],301989.5),
"median CD3CD8 events",step2$median_cd3cd8_events[1],12411.5,near(step2$median_cd3cd8_events[1],12411.5),
"minimum CD3CD8 events",step2$min_cd3cd8_events[1],407,step2$min_cd3cd8_events[1]==407,
"maximum CD3CD8 events",step2$max_cd3cd8_events[1],59837,step2$max_cd3cd8_events[1]==59837,
"valid ages",meta$n_valid_age[1],832,meta$n_valid_age[1]==832,
"main feature count",nrow(st),55,nrow(st)==55,
"CD3 beta",get_beta("pct_cd3_pos_total"),-4.71800773276706,near(get_beta("pct_cd3_pos_total"),-4.71800773276706),
"CD8 within CD3 beta",get_beta("pct_cd8_within_cd3"),4.28759233054575,near(get_beta("pct_cd8_within_cd3"),4.28759233054575),
"TEMRA beta",get_beta("pct_temra_like"),4.36723790655398,near(get_beta("pct_temra_like"),4.36723790655398),
"CD39 beta",get_beta("pct_cd39_pos"),3.23434446149965,near(get_beta("pct_cd39_pos"),3.23434446149965),
"integrated beta",sc$beta_cancer_vs_healthy[sc$score==CP7_MANIFEST$integrated_score][1],CP7_MANIFEST$benchmarks$integrated_beta,near(sc$beta_cancer_vs_healthy[sc$score==CP7_MANIFEST$integrated_score][1],CP7_MANIFEST$benchmarks$integrated_beta),
"5y matched pairs",mt$n_pairs[mt$caliper_years==5][1],302,mt$n_pairs[mt$caliper_years==5][1]==302,
"10y matched pairs",mt$n_pairs[mt$caliper_years==10][1],332,mt$n_pairs[mt$caliper_years==10][1]==332,
"threshold configurations",sum(c("direction_main","direction_permissive","direction_stringent")%in%names(ts)),3,all(c("direction_main","direction_permissive","direction_stringent")%in%names(ts)),
"quantitative figures",length(list.files(file.path(a,"09_main_figures_bracketed"),"\\.pdf$")),5,length(list.files(file.path(a,"09_main_figures_bracketed"),"\\.pdf$"))>=5,
"representative atlas",file.exists(file.path(a,"05_flow_figures_flowjo_style_TALL_LAYOUT","FigureS_CP7_representative_checkpoint_gating_RECONSTRUCTED.pdf")),TRUE,file.exists(file.path(a,"05_flow_figures_flowjo_style_TALL_LAYOUT","FigureS_CP7_representative_checkpoint_gating_RECONSTRUCTED.pdf")))
ref<-path.expand(Sys.getenv("SDY2583_CP7_REFERENCE_DIR",unset="")); if(nzchar(ref)&&dir.exists(ref)){pairs<-list(c(file.path(a,"02_feature_extraction","SDY2583_CP7_FULL_850_feature_table_STEP2_RECONSTRUCTED.csv"),"SDY2583_CP7_FULL_850_feature_table_STEP2.csv","subject_id"),c(file.path(a,"06_statistics","SDY2583_CP7_FULL850_age_sex_adjusted_statistics_MAIN_RECONSTRUCTED.csv"),"SDY2583_CP7_FULL850_age_sex_adjusted_statistics_MAIN.csv","module,feature"),c(file.path(a,"08_composite_scores","SDY2583_CP7_composite_score_statistics_RECONSTRUCTED.csv"),"SDY2583_CP7_composite_score_statistics.csv","score"),c(file.path(a,"10_age_sensitivity","SDY2583_CP7_age_matched_sensitivity_statistics_RECONSTRUCTED.csv"),"SDY2583_CP7_age_matched_sensitivity_statistics.csv","subset_label,feature"),c(file.path(a,"12_threshold_sensitivity","SDY2583_CP7_threshold_sensitivity_summary_WIDE_RECONSTRUCTED.csv"),"SDY2583_CP7_threshold_sensitivity_summary_WIDE.csv","feature"));for(p in pairs){hit<-list.files(ref,p[2],full.names=TRUE,recursive=TRUE,fixed=TRUE);cmp<-if(length(hit))rp_compare_csv(p[1],hit[1],strsplit(p[3],",",fixed=TRUE)[[1]],tol)else tibble::tibble(pass=FALSE,detail="reference missing");checks<-dplyr::bind_rows(checks,tibble::tibble(check=p[2],observed=cmp$detail[1],expected="table match",pass=cmp$pass[1]))}}
readr::write_csv(checks,file.path(out,"SDY2583_CP7_reconstruction_validation_report.csv"));print(checks,n=nrow(checks));if(any(!checks$pass)&&strict)stop("CP7 validation failed.")
