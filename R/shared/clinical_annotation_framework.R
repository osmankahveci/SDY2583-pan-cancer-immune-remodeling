# Generic cancer-subgroup and treatment-annotation framework for reconstructed
# SDY2583 panels. Participant-level outputs remain local and ignored by Git.
source(file.path(Sys.getenv("SDY2583_REPO_ROOT", unset="."),"R","shared","reconstructed_panel_framework.R"))

rp_as01 <- function(x) {
  if (is.logical(x)) return(as.integer(x))
  if (is.numeric(x)) return(ifelse(is.na(x),NA_integer_,as.integer(x!=0)))
  y <- tolower(trimws(as.character(x)))
  dplyr::case_when(y %in% c("1","yes","y","true","present","positive","evet","var")~1L,
                   y %in% c("0","no","n","false","absent","negative","hayir","hayır","yok")~0L,
                   TRUE~NA_integer_)
}

rp_attach_clinical_metadata <- function(scored_data) {
  source_file <- rp_find_metadata_file()
  raw <- readr::read_csv(source_file,show_col_types=FALSE,progress=FALSE)
  subject_col <- rp_first_col(raw,c("subject_id","dbg_id","subject"))
  subgroup_col <- rp_first_col(raw,c("cancer_subgroup_model","cancer_subgroup","cancer_type_model","cancer_type","tumor_type"))
  therapy_col <- rp_first_col(raw,c("therapy_status_model","therapy_status_4level","therapy_status","treatment_status"))
  line_col <- rp_first_col(raw,c("therapy_line_number_model","therapy_line_number","treatment_line_number","line_of_therapy"))
  time_col <- rp_first_col(raw,c("time_from_start_days_model","time_from_start_days","time_from_treatment_start_days","days_from_treatment_start"))
  clinical <- tibble::tibble(subject_id=as.character(raw[[subject_col]]))
  clinical$cancer_subgroup <- if(is.na(subgroup_col)) NA_character_ else as.character(raw[[subgroup_col]])
  clinical$therapy_status_4level <- if(is.na(therapy_col)) NA_character_ else as.character(raw[[therapy_col]])
  clinical$therapy_line_number <- if(is.na(line_col)) NA_real_ else suppressWarnings(as.numeric(raw[[line_col]]))
  clinical$time_from_start_days <- if(is.na(time_col)) NA_real_ else suppressWarnings(as.numeric(raw[[time_col]]))
  exposures <- c("chemotherapy","targeted_therapy","any_immunotherapy","ici_immunotherapy","endocrine_hormonal","adc","experimental","radiotherapy")
  for(x in exposures) clinical[[x]] <- if(x %in% names(raw)) rp_as01(raw[[x]]) else NA_integer_
  clinical <- clinical |> dplyr::filter(!is.na(subject_id),subject_id!="") |> dplyr::group_by(subject_id) |> dplyr::summarise(dplyr::across(dplyr::everything(),~dplyr::first(stats::na.omit(.x),default=NA)),.groups="drop")
  scored_data |> dplyr::left_join(clinical,by="subject_id")
}

rp_fit_cancer_one_vs_rest <- function(data,outcomes,min_group=15) {
  cancer <- data |> dplyr::filter(disease_group=="Cancer patient",!is.na(cancer_subgroup),cancer_subgroup!="")
  groups <- names(which(table(cancer$cancer_subgroup)>=min_group))
  purrr::map_dfr(groups,function(group) {
    d <- cancer |> dplyr::mutate(group_indicator=factor(ifelse(cancer_subgroup==group,"Target","Other"),levels=c("Other","Target")))
    purrr::map_dfr(outcomes,function(outcome) {
      m <- d |> dplyr::transmute(y=as.numeric(.data[[outcome]]),group_indicator,age_for_model=as.numeric(age_for_model),sex=factor(as.character(sex))) |> dplyr::filter(stats::complete.cases(.))
      if(nrow(m)<20||dplyr::n_distinct(m$group_indicator)<2) return(tibble::tibble(cancer_subgroup=group,feature=outcome,n_model=nrow(m),beta=NA_real_,ci_low=NA_real_,ci_high=NA_real_,p_value=NA_real_))
      fit<-stats::lm(y~group_indicator+age_for_model+sex,data=m); ct<-summary(fit)$coefficients; cn<-grep("^group_indicator",rownames(ct),value=TRUE)[1]; ci<-stats::confint(fit,cn)
      tibble::tibble(cancer_subgroup=group,feature=outcome,n_model=stats::nobs(fit),beta=unname(ct[cn,"Estimate"]),ci_low=unname(ci[1]),ci_high=unname(ci[2]),p_value=unname(ct[cn,"Pr(>|t|)"]))
    })
  }) |> dplyr::group_by(cancer_subgroup) |> dplyr::mutate(fdr=p.adjust(p_value,"BH"),direction=ifelse(beta>0,"higher_in_target",ifelse(beta<0,"lower_in_target","no_difference"))) |> dplyr::ungroup()
}

rp_fit_exposure_models <- function(data,outcomes) {
  exposures<-c("chemotherapy","targeted_therapy","any_immunotherapy","ici_immunotherapy","endocrine_hormonal","adc","experimental","radiotherapy")
  purrr::map_dfr(exposures,function(exposure) {
    if(!exposure %in% names(data)||sum(!is.na(data[[exposure]]))<20||dplyr::n_distinct(stats::na.omit(data[[exposure]]))<2) return(NULL)
    purrr::map_dfr(outcomes,function(outcome) {
      d<-data |> dplyr::filter(disease_group=="Cancer patient") |> dplyr::transmute(y=as.numeric(.data[[outcome]]),exposure=factor(.data[[exposure]],levels=c(0,1)),age_for_model=as.numeric(age_for_model),sex=factor(as.character(sex))) |> dplyr::filter(stats::complete.cases(.))
      if(nrow(d)<20||dplyr::n_distinct(d$exposure)<2) return(tibble::tibble(exposure=exposure,feature=outcome,n_model=nrow(d),beta=NA_real_,ci_low=NA_real_,ci_high=NA_real_,p_value=NA_real_))
      fit<-stats::lm(y~exposure+age_for_model+sex,data=d); ct<-summary(fit)$coefficients; cn<-grep("^exposure",rownames(ct),value=TRUE)[1]; ci<-stats::confint(fit,cn)
      tibble::tibble(exposure=exposure,feature=outcome,n_model=stats::nobs(fit),beta=unname(ct[cn,"Estimate"]),ci_low=unname(ci[1]),ci_high=unname(ci[2]),p_value=unname(ct[cn,"Pr(>|t|)"]))
    })
  }) |> dplyr::group_by(exposure) |> dplyr::mutate(fdr=p.adjust(p_value,"BH"),direction=ifelse(beta>0,"higher_with_exposure",ifelse(beta<0,"lower_with_exposure","no_difference"))) |> dplyr::ungroup()
}

rp_run_clinical_annotation <- function(panel,scored_data,outcomes,min_group=15) {
  data<-rp_attach_clinical_metadata(scored_data)
  subgroup<-rp_fit_cancer_one_vs_rest(data,outcomes,min_group)
  exposure<-rp_fit_exposure_models(data,outcomes)
  out_dir<-file.path(sd_analysis_dir(panel),"clinical_annotation"); dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
  readr::write_csv(data,file.path(out_dir,paste0("SDY2583_",panel,"_analysis_data_with_clinical_annotation_RECONSTRUCTED.csv")))
  readr::write_csv(subgroup,file.path(out_dir,paste0("SDY2583_",panel,"_cancer_subgroup_one_vs_rest_RECONSTRUCTED.csv")))
  readr::write_csv(exposure,file.path(out_dir,paste0("SDY2583_",panel,"_therapy_exposure_models_RECONSTRUCTED.csv")))
  list(data=data,subgroup=subgroup,exposure=exposure)
}
