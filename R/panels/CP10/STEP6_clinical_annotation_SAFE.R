# ============================================================
# SDY2583 CP10 STEP 6 SAFE: Clinical annotation + ALL8 matrix
# ============================================================

rm(list = ls())
source(file.path(Sys.getenv("SDY2583_REPO_ROOT", unset = "."), "R", "shared", "bootstrap.R"))

pkgs <- c("dplyr","readr","stringr","tibble","broom","tidyr","purrr")
for (p in pkgs) if (!requireNamespace(p, quietly=TRUE)) install.packages(p)

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(stringr); library(tibble)
  library(broom); library(tidyr); library(purrr)
})

filter <- dplyr::filter; select <- dplyr::select; mutate <- dplyr::mutate
arrange <- dplyr::arrange; summarise <- dplyr::summarise; group_by <- dplyr::group_by
ungroup <- dplyr::ungroup; count <- dplyr::count; distinct <- dplyr::distinct
case_when <- dplyr::case_when; bind_rows <- dplyr::bind_rows; left_join <- dplyr::left_join
n_distinct <- dplyr::n_distinct; n <- dplyr::n

# -----------------------------
# Paths
# -----------------------------
cp10_dir <- sd_analysis_dir("CP10")
cp10_rdata_dir <- file.path(cp10_dir, "11_RData")
matrix_dir <- sd_integrated_dir()
out_dir <- file.path(matrix_dir, "06_CP10_clinical_annotation")
rdata_out_dir <- file.path(matrix_dir, "RData")
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(rdata_out_dir, recursive=TRUE, showWarnings=FALSE)

load(file.path(cp10_rdata_dir, "SDY2583_CP10_STEP4_composite_scores_for_integration.RData"))
load(file.path(cp10_rdata_dir, "SDY2583_CP10_STEP4_composite_scores.RData"))

if (!exists("cp10_scores_for_integration")) stop("cp10_scores_for_integration yok.")
cp10_score_cols <- grep("^CP10_.*_score$", names(cp10_scores_for_integration), value=TRUE)
if (length(cp10_score_cols) == 0) stop("CP10 score kolonu yok.")

# -----------------------------
# Helpers
# -----------------------------
first_existing_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0) NA_character_ else hit[1]
}

as_01 <- function(x) {
  if (is.logical(x)) return(as.integer(x))
  if (is.numeric(x)) return(ifelse(is.na(x), NA_integer_, ifelse(x != 0, 1L, 0L)))
  y <- tolower(trimws(as.character(x)))
  case_when(
    y %in% c("1","yes","y","true","t","present","positive","var","evet") ~ 1L,
    y %in% c("0","no","n","false","f","absent","negative","yok","hayir","hayır") ~ 0L,
    TRUE ~ NA_integer_
  )
}

clean_term <- function(x) stringr::str_replace_all(x, "[^A-Za-z0-9_]+", "_")

# -----------------------------
# 1) Select previous integrated matrix
# -----------------------------
csvs <- list.files(matrix_dir, pattern="\\.csv$", full.names=TRUE, recursive=TRUE)
csvs <- csvs[grepl("integrated|matrix|ALL6|ALL7|ALL8", basename(csvs), ignore.case=TRUE)]

inspect_csv <- function(fp) {
  x <- tryCatch(readr::read_csv(fp, show_col_types=FALSE), error=function(e) NULL)
  if (is.null(x)) return(tibble(file_path=fp, file_name=basename(fp), readable=FALSE, n_rows=NA_integer_, n_cols=NA_integer_, has_subject_id=FALSE, has_disease=FALSE, has_age=FALSE, has_sex=FALSE, has_clinical=FALSE, n_score_cols=NA_integer_, priority=-Inf))
  nm <- names(x)
  score_cols <- grep("^CP[0-9]+_.*_score$", nm, value=TRUE)
  has_disease <- any(c("disease_group_model","disease_group_clinical","disease_group","group","condition") %in% nm)
  has_age <- any(c("age_for_model","age_clinical","age_years","age") %in% nm)
  has_sex <- any(c("sex_for_clinical_model","sex_clinical","sex","gender") %in% nm)
  has_clinical <- any(c("cancer_subgroup","cancer_subgroup_model","therapy_status_4level","therapy_status_model","chemotherapy","targeted_therapy") %in% nm)
  priority <- 0
  priority <- priority + ifelse("subject_id" %in% nm, 1000, 0)
  priority <- priority + ifelse(nrow(x) >= 800 & nrow(x) <= 900, 300, 0)
  priority <- priority + ifelse(has_disease, 100, 0) + ifelse(has_age,100,0) + ifelse(has_sex,100,0) + ifelse(has_clinical,100,0)
  priority <- priority + length(score_cols)
  priority <- priority + ifelse(grepl("ALL7", basename(fp), ignore.case=TRUE), 150, 0)
  priority <- priority + ifelse(grepl("with_CP22", basename(fp), ignore.case=TRUE), 150, 0)
  priority <- priority - ifelse(grepl("summary|dictionary|availability|counts|results|variable_sources|exposure|therapy|subgroup|CP10_clinical", basename(fp), ignore.case=TRUE), 600, 0)
  tibble(file_path=fp, file_name=basename(fp), readable=TRUE, n_rows=nrow(x), n_cols=ncol(x),
         has_subject_id="subject_id" %in% nm, has_disease=has_disease, has_age=has_age,
         has_sex=has_sex, has_clinical=has_clinical, n_score_cols=length(score_cols), priority=priority)
}

candidate_inventory <- bind_rows(lapply(csvs, inspect_csv)) %>% arrange(desc(priority), desc(n_score_cols), desc(n_cols))
write_csv(candidate_inventory, file.path(out_dir, "SDY2583_CP10_integrated_matrix_candidate_inventory_STEP6.csv"))

valid <- candidate_inventory %>%
  filter(readable, has_subject_id, has_disease, has_age, has_sex, n_rows >= 800, n_rows <= 900)

if (nrow(valid) == 0) {
  print(as.data.frame(candidate_inventory), row.names=FALSE)
  stop("Uygun integrated matrix bulunamadı.")
}

selected_matrix_file <- valid$file_path[1]
cat("\nSelected previous integrated matrix:\n"); print(selected_matrix_file)

mat_prev <- readr::read_csv(selected_matrix_file, show_col_types=FALSE)

# Remove old CP10 score columns if script is rerun
mat_no_cp10 <- mat_prev %>% select(-any_of(grep("^CP10_.*_score$", names(mat_prev), value=TRUE)))
mat_all8 <- mat_no_cp10 %>% left_join(cp10_scores_for_integration, by="subject_id")

all_score_cols <- grep("^CP[0-9]+_.*_score$", names(mat_all8), value=TRUE)
updated_matrix_file <- file.path(matrix_dir, "SDY2583_integrated_clinical_immune_score_matrix_ALL8_with_CP10.csv")
write_csv(mat_all8, updated_matrix_file)

score_availability_by_panel <- bind_rows(lapply(unique(str_extract(all_score_cols, "^CP[0-9]+")), function(pp) {
  cols <- all_score_cols[str_detect(all_score_cols, paste0("^", pp, "_"))]
  tibble(panel=pp, n_score_cols=length(cols),
         n_subjects_with_any_score=sum(rowSums(!is.na(mat_all8[, cols, drop=FALSE])) > 0))
})) %>% arrange(panel)

all8_summary <- tibble(
  n_subjects=nrow(mat_all8),
  n_unique_subjects=n_distinct(mat_all8$subject_id),
  n_with_any_CP10_score=sum(rowSums(!is.na(mat_all8[, cp10_score_cols, drop=FALSE])) > 0),
  n_with_all_CP10_scores=sum(rowSums(!is.na(mat_all8[, cp10_score_cols, drop=FALSE])) == length(cp10_score_cols)),
  n_total_score_cols_ALL8=length(all_score_cols),
  median_available_score_cols=median(rowSums(!is.na(mat_all8[, all_score_cols, drop=FALSE])), na.rm=TRUE)
)

write_csv(all8_summary, file.path(out_dir, "SDY2583_CP10_ALL8_integrated_matrix_summary_STEP6.csv"))
write_csv(score_availability_by_panel, file.path(out_dir, "SDY2583_CP10_score_availability_by_panel_after_addition_STEP6.csv"))

# -----------------------------
# 2) Standardize clinical variables
# -----------------------------
df <- mat_all8

disease_col <- first_existing_col(df, c("disease_group_model","disease_group_clinical","disease_group","group","condition"))
age_col <- first_existing_col(df, c("age_for_model","age_clinical","age_years","age"))
sex_col <- first_existing_col(df, c("sex_for_clinical_model","sex_clinical","sex","gender"))
subgroup_col <- first_existing_col(df, c("cancer_subgroup_model","cancer_subgroup","cancer_type_model","cancer_type","tumor_type"))
therapy_col <- first_existing_col(df, c("therapy_status_model","therapy_status_4level","therapy_status","treatment_status"))
line_col <- first_existing_col(df, c("therapy_line_number_model","therapy_line_number","treatment_line_number","line_of_therapy"))
time_col <- first_existing_col(df, c("time_from_start_days_model","time_from_start_days","time_from_treatment_start_days","days_from_treatment_start"))

if (any(is.na(c(disease_col, age_col, sex_col, subgroup_col, therapy_col)))) {
  stop("Zorunlu klinik kolonlardan biri bulunamadı: disease/age/sex/subgroup/therapy_status")
}

clinical_df <- df %>%
  mutate(
    disease_group_model = as.character(.data[[disease_col]]),
    disease_group_model = case_when(
      disease_group_model %in% c("Healthy control","Healthy","healthy","control","Control") ~ "Healthy control",
      disease_group_model %in% c("Cancer patient","Cancer","cancer","patient","Patient") ~ "Cancer patient",
      str_detect(disease_group_model, regex("healthy|control", ignore_case=TRUE)) ~ "Healthy control",
      str_detect(disease_group_model, regex("cancer|patient", ignore_case=TRUE)) ~ "Cancer patient",
      TRUE ~ disease_group_model
    ),
    disease_group_model = factor(disease_group_model, levels=c("Healthy control","Cancer patient")),
    age_clinical_model = suppressWarnings(as.numeric(.data[[age_col]])),
    age_clinical_model = ifelse(age_clinical_model < 18 | age_clinical_model > 100, NA_real_, age_clinical_model),
    sex_clinical_model = case_when(
      as.character(.data[[sex_col]]) %in% c("Female","F","female","f") ~ "Female",
      as.character(.data[[sex_col]]) %in% c("Male","M","male","m") ~ "Male",
      TRUE ~ as.character(.data[[sex_col]])
    ),
    sex_clinical_model = factor(sex_clinical_model),
    cancer_subgroup_model = factor(as.character(.data[[subgroup_col]])),
    therapy_status_model = as.character(.data[[therapy_col]]),
    therapy_status_model = case_when(
      str_detect(therapy_status_model, regex("ongoing_active|active", ignore_case=TRUE)) ~ "ongoing_active_treatment",
      str_detect(therapy_status_model, regex("no_ongoing", ignore_case=TRUE)) ~ "no_ongoing_therapy",
      str_detect(therapy_status_model, regex("no_treatment_data|missing|unknown|NA", ignore_case=TRUE)) ~ "no_treatment_data",
      TRUE ~ therapy_status_model
    )
  )

clinical_df$therapy_line_number_model <- if (!is.na(line_col)) suppressWarnings(as.numeric(clinical_df[[line_col]])) else NA_real_
clinical_df$therapy_line_group <- factor(case_when(
  is.na(clinical_df$therapy_line_number_model) ~ NA_character_,
  clinical_df$therapy_line_number_model == 1 ~ "first_line",
  clinical_df$therapy_line_number_model >= 2 ~ "later_line",
  TRUE ~ NA_character_
), levels=c("first_line","later_line"))

clinical_df$time_from_start_days_model <- if (!is.na(time_col)) suppressWarnings(as.numeric(clinical_df[[time_col]])) else NA_real_

exposure_vars <- c("chemotherapy","targeted_therapy","any_immunotherapy","ici_immunotherapy","endocrine_hormonal","adc","experimental","radiotherapy")
for (ee in exposure_vars) clinical_df[[paste0(ee, "_model")]] <- if (ee %in% names(clinical_df)) as_01(clinical_df[[ee]]) else NA_integer_

variable_sources <- tibble(
  variable=c("disease_group","age","sex","cancer_subgroup","therapy_status","therapy_line_number","time_from_start_days"),
  source_column=c(disease_col, age_col, sex_col, subgroup_col, therapy_col, line_col, time_col)
)
write_csv(variable_sources, file.path(out_dir, "SDY2583_CP10_clinical_variable_sources_STEP6.csv"))

# -----------------------------
# 3) Summary tables
# -----------------------------
cancer_df <- clinical_df %>%
  filter(disease_group_model == "Cancer patient") %>%
  mutate(core_covariates_ok = !is.na(age_clinical_model) & !is.na(sex_clinical_model) & !is.na(cancer_subgroup_model))

clinical_model_summary <- tibble(
  n_cancer_total=nrow(cancer_df),
  n_with_valid_age=sum(!is.na(cancer_df$age_clinical_model)),
  n_with_valid_sex=sum(!is.na(cancer_df$sex_clinical_model)),
  n_with_cancer_subgroup=sum(!is.na(cancer_df$cancer_subgroup_model)),
  n_with_therapy_status=sum(!is.na(cancer_df$therapy_status_model)),
  n_with_all_core_covariates=sum(cancer_df$core_covariates_ok, na.rm=TRUE),
  n_ongoing_active_treatment=sum(cancer_df$therapy_status_model == "ongoing_active_treatment", na.rm=TRUE),
  n_no_ongoing_therapy=sum(cancer_df$therapy_status_model == "no_ongoing_therapy", na.rm=TRUE),
  n_no_treatment_data=sum(cancer_df$therapy_status_model == "no_treatment_data", na.rm=TRUE)
)

cancer_subgroup_counts <- cancer_df %>% count(cancer_subgroup_model, name="n") %>% arrange(desc(n))
therapy_status_counts <- cancer_df %>% count(therapy_status_model, name="n") %>% arrange(desc(n))

active_treatment_df <- cancer_df %>% filter(therapy_status_model == "ongoing_active_treatment", core_covariates_ok)
active_exposure_counts <- bind_rows(lapply(exposure_vars, function(ee) {
  vv <- paste0(ee, "_model")
  ne <- sum(active_treatment_df[[vv]] == 1, na.rm=TRUE)
  nu <- sum(active_treatment_df[[vv]] == 0, na.rm=TRUE)
  tibble(exposure=ee, n_exposed=ne, n_unexposed=nu,
         analysis_tier=case_when(
           ne >= 20 & nu >= 20 ~ "primary_or_exploratory_powered",
           ne >= 10 & nu >= 20 ~ "exploratory_underpowered",
           TRUE ~ "descriptive_only"
         ))
}))

therapy_line_counts <- cancer_df %>% count(therapy_line_number_model, therapy_line_group, name="n") %>% arrange(therapy_line_number_model)
time_from_start_summary <- cancer_df %>%
  summarise(n_with_time=sum(!is.na(time_from_start_days_model)),
            min_days=min(time_from_start_days_model, na.rm=TRUE),
            q1_days=as.numeric(quantile(time_from_start_days_model, 0.25, na.rm=TRUE)),
            median_days=median(time_from_start_days_model, na.rm=TRUE),
            mean_days=mean(time_from_start_days_model, na.rm=TRUE),
            q3_days=as.numeric(quantile(time_from_start_days_model, 0.75, na.rm=TRUE)),
            max_days=max(time_from_start_days_model, na.rm=TRUE))

write_csv(clinical_model_summary, file.path(out_dir, "SDY2583_CP10_clinical_model_summary_STEP6.csv"))
write_csv(cancer_subgroup_counts, file.path(out_dir, "SDY2583_CP10_cancer_subgroup_counts_STEP6.csv"))
write_csv(therapy_status_counts, file.path(out_dir, "SDY2583_CP10_therapy_status_counts_STEP6.csv"))
write_csv(active_exposure_counts, file.path(out_dir, "SDY2583_CP10_active_treatment_exposure_counts_STEP6.csv"))
write_csv(therapy_line_counts, file.path(out_dir, "SDY2583_CP10_therapy_line_counts_STEP6.csv"))
write_csv(time_from_start_summary, file.path(out_dir, "SDY2583_CP10_time_from_start_summary_STEP6.csv"))

# -----------------------------
# 4) Model helpers
# -----------------------------
drop_bad_covars <- function(model_df, covars) {
  keep <- covars
  for (cc in covars) {
    if (is.character(model_df[[cc]]) || is.factor(model_df[[cc]])) {
      if (n_distinct(model_df[[cc]]) < 2) keep <- setdiff(keep, cc)
    }
  }
  keep
}

run_coef_model <- function(df, score, predictor, covars=c("age_clinical_model","sex_clinical_model","cancer_subgroup_model")) {
  model_df <- df %>%
    select(all_of(c(score, predictor, covars))) %>%
    mutate(value=suppressWarnings(as.numeric(.data[[score]]))) %>%
    filter(!is.na(value), !is.na(.data[[predictor]]))
  for (cc in covars) model_df <- model_df %>% filter(!is.na(.data[[cc]]))
  if (nrow(model_df) < 30 || n_distinct(model_df[[predictor]]) < 2) return(NULL)
  covars2 <- drop_bad_covars(model_df, covars)
  f <- as.formula(paste("value ~", paste(c(predictor, covars2), collapse=" + ")))
  fit <- tryCatch(lm(f, data=model_df), error=function(e) NULL)
  if (is.null(fit)) return(NULL)
  tt <- tryCatch(broom::tidy(fit, conf.int=TRUE), error=function(e) NULL)
  if (is.null(tt)) return(NULL)
  row <- tt %>% filter(str_detect(term, paste0("^", predictor)))
  if (nrow(row) == 0) return(NULL)
  row <- row[1,]
  tibble(score=score, n_model=nrow(model_df), beta=row$estimate, conf_low=row$conf.low,
         conf_high=row$conf.high, p_value=row$p.value,
         direction=case_when(beta > 0 ~ "higher_in_exposed_or_group", beta < 0 ~ "lower_in_exposed_or_group", TRUE ~ "no_direction"))
}

run_numeric_model <- function(df, score, predictor, covars=c("age_clinical_model","sex_clinical_model","cancer_subgroup_model")) {
  model_df <- df %>%
    select(all_of(c(score, predictor, covars))) %>%
    mutate(value=suppressWarnings(as.numeric(.data[[score]])),
           predictor_value=suppressWarnings(as.numeric(.data[[predictor]]))) %>%
    filter(!is.na(value), !is.na(predictor_value))
  for (cc in covars) model_df <- model_df %>% filter(!is.na(.data[[cc]]))
  if (nrow(model_df) < 30) return(NULL)
  covars2 <- drop_bad_covars(model_df, covars)
  model_df <- model_df %>% mutate(predictor_z=as.numeric(scale(log1p(predictor_value))))
  f <- as.formula(paste("value ~", paste(c("predictor_z", covars2), collapse=" + ")))
  fit <- tryCatch(lm(f, data=model_df), error=function(e) NULL)
  if (is.null(fit)) return(NULL)
  tt <- tryCatch(broom::tidy(fit, conf.int=TRUE), error=function(e) NULL)
  if (is.null(tt)) return(NULL)
  row <- tt %>% filter(term == "predictor_z")
  if (nrow(row) == 0) return(NULL)
  tibble(score=score, n_model=nrow(model_df), median_time_days=median(model_df$predictor_value, na.rm=TRUE),
         beta_per_log_time_z=row$estimate[1], conf_low=row$conf.low[1], conf_high=row$conf.high[1],
         p_value=row$p.value[1],
         direction=case_when(beta_per_log_time_z > 0 ~ "higher_with_longer_time_from_start",
                             beta_per_log_time_z < 0 ~ "lower_with_longer_time_from_start",
                             TRUE ~ "no_direction"))
}

# -----------------------------
# 5) Cancer subgroup omnibus
# -----------------------------
run_omnibus <- function(score) {
  model_df <- cancer_df %>%
    transmute(value=suppressWarnings(as.numeric(.data[[score]])),
              cancer_subgroup_model=cancer_subgroup_model,
              age_clinical_model=age_clinical_model,
              sex_clinical_model=sex_clinical_model) %>%
    filter(!is.na(value), !is.na(cancer_subgroup_model), !is.na(age_clinical_model), !is.na(sex_clinical_model)) %>%
    mutate(cancer_subgroup_model=droplevels(factor(cancer_subgroup_model)),
           sex_clinical_model=droplevels(factor(sex_clinical_model)))
  if (nrow(model_df) < 50 || n_distinct(model_df$cancer_subgroup_model) < 2) return(NULL)
  fit <- tryCatch(lm(value ~ cancer_subgroup_model + age_clinical_model + sex_clinical_model, data=model_df), error=function(e) NULL)
  if (is.null(fit)) return(NULL)
  dd <- tryCatch(drop1(fit, test="F"), error=function(e) NULL)
  if (is.null(dd)) return(NULL)
  dd <- as.data.frame(dd); dd$term <- rownames(dd)
  row <- dd %>% filter(term == "cancer_subgroup_model")
  tibble(score=score, n_model=nrow(model_df), n_subgroups=n_distinct(model_df$cancer_subgroup_model),
         omnibus_F=row$`F value`[1], omnibus_p=row$`Pr(>F)`[1])
}

subgroup_omnibus_results <- bind_rows(lapply(cp10_score_cols, run_omnibus)) %>%
  mutate(FDR_global=p.adjust(omnibus_p, method="BH")) %>% arrange(FDR_global, omnibus_p)

# -----------------------------
# 6) One-vs-rest subgroup
# -----------------------------
subgroups <- cancer_subgroup_counts %>% filter(!is.na(cancer_subgroup_model), n >= 15) %>% pull(cancer_subgroup_model) %>% as.character()

run_one_vs_rest <- function(score, sg) {
  df1 <- cancer_df %>%
    filter(!is.na(cancer_subgroup_model), !is.na(age_clinical_model), !is.na(sex_clinical_model)) %>%
    mutate(subgroup_binary=factor(ifelse(as.character(cancer_subgroup_model)==sg, "target_subgroup", "rest"), levels=c("rest","target_subgroup")))
  n_sg <- sum(df1$subgroup_binary == "target_subgroup")
  n_rest <- sum(df1$subgroup_binary == "rest")
  if (n_sg < 15 || n_rest < 20) return(NULL)
  res <- run_coef_model(df1, score, "subgroup_binary", covars=c("age_clinical_model","sex_clinical_model"))
  if (is.null(res)) return(NULL)
  res %>% mutate(cancer_subgroup=sg, n_subgroup=n_sg, n_rest=n_rest,
                 beta_subgroup_vs_rest=beta,
                 direction=case_when(beta_subgroup_vs_rest > 0 ~ "higher_in_subgroup_vs_rest",
                                     beta_subgroup_vs_rest < 0 ~ "lower_in_subgroup_vs_rest",
                                     TRUE ~ "no_direction")) %>%
    select(score, cancer_subgroup, n_model, n_subgroup, n_rest, beta_subgroup_vs_rest, conf_low, conf_high, p_value, direction)
}

subgroup_one_vs_rest_results <- bind_rows(lapply(cp10_score_cols, function(sc) bind_rows(lapply(subgroups, function(sg) run_one_vs_rest(sc, sg))))) %>%
  group_by(cancer_subgroup) %>% mutate(FDR_within_subgroup=p.adjust(p_value, method="BH")) %>% ungroup() %>%
  group_by(score) %>% mutate(FDR_within_score=p.adjust(p_value, method="BH")) %>% ungroup() %>%
  mutate(FDR_global=p.adjust(p_value, method="BH")) %>% arrange(FDR_global, p_value)

# -----------------------------
# 7) Therapy status
# -----------------------------
therapy_status_df <- cancer_df %>%
  filter(therapy_status_model %in% c("ongoing_active_treatment","no_ongoing_therapy"), core_covariates_ok) %>%
  mutate(therapy_status_binary=factor(ifelse(therapy_status_model=="ongoing_active_treatment","ongoing_active_treatment","reference"),
                                      levels=c("reference","ongoing_active_treatment")))

therapy_status_results <- bind_rows(lapply(cp10_score_cols, function(sc) {
  res <- run_coef_model(therapy_status_df, sc, "therapy_status_binary")
  if (is.null(res)) return(NULL)
  res %>% mutate(n_ongoing=sum(therapy_status_df$therapy_status_binary=="ongoing_active_treatment"),
                 n_no_ongoing=sum(therapy_status_df$therapy_status_binary=="reference"),
                 beta_ongoing_vs_no_ongoing=beta,
                 direction=case_when(beta_ongoing_vs_no_ongoing > 0 ~ "higher_in_ongoing_active_treatment",
                                     beta_ongoing_vs_no_ongoing < 0 ~ "lower_in_ongoing_active_treatment",
                                     TRUE ~ "no_direction")) %>%
    select(score, n_model, n_ongoing, n_no_ongoing, beta_ongoing_vs_no_ongoing, conf_low, conf_high, p_value, direction)
})) %>% mutate(FDR_global=p.adjust(p_value, method="BH")) %>% arrange(FDR_global, p_value)

# -----------------------------
# 8) Active treatment modality
# -----------------------------
modality_results <- bind_rows(lapply(exposure_vars, function(ee) {
  vv <- paste0(ee, "_model")
  edf <- active_treatment_df %>%
    filter(!is.na(.data[[vv]])) %>%
    mutate(exposure_binary=factor(ifelse(.data[[vv]]==1, "exposed", "reference"), levels=c("reference","exposed")))
  n_exp <- sum(edf$exposure_binary=="exposed")
  n_unexp <- sum(edf$exposure_binary=="reference")
  tier <- case_when(n_exp >= 20 & n_unexp >= 20 ~ "primary_or_exploratory_powered",
                    n_exp >= 10 & n_unexp >= 20 ~ "exploratory_underpowered",
                    TRUE ~ "descriptive_only")
  if (n_exp < 5 || n_unexp < 10) return(NULL)
  bind_rows(lapply(cp10_score_cols, function(sc) {
    res <- run_coef_model(edf, sc, "exposure_binary")
    if (is.null(res)) return(NULL)
    res %>% mutate(exposure=ee, n_exposed=n_exp, n_unexposed=n_unexp, analysis_tier=tier,
                   beta_exposed_vs_unexposed=beta,
                   direction=case_when(beta_exposed_vs_unexposed > 0 ~ "higher_in_exposed",
                                       beta_exposed_vs_unexposed < 0 ~ "lower_in_exposed",
                                       TRUE ~ "no_direction")) %>%
      select(score, exposure, n_model, n_exposed, n_unexposed, analysis_tier, beta_exposed_vs_unexposed, conf_low, conf_high, p_value, direction)
  }))
})) %>%
  group_by(exposure) %>% mutate(FDR_within_exposure=p.adjust(p_value, method="BH")) %>% ungroup() %>%
  group_by(score) %>% mutate(FDR_within_score=p.adjust(p_value, method="BH")) %>% ungroup() %>%
  mutate(FDR_global=p.adjust(p_value, method="BH")) %>% arrange(FDR_global, p_value)

# -----------------------------
# 9) Therapy line and time
# -----------------------------
therapy_line_df <- cancer_df %>%
  filter(therapy_status_model=="ongoing_active_treatment", core_covariates_ok, !is.na(therapy_line_group)) %>%
  mutate(therapy_line_binary=factor(as.character(therapy_line_group), levels=c("first_line","later_line")))

therapy_line_results <- bind_rows(lapply(cp10_score_cols, function(sc) {
  if (nrow(therapy_line_df) < 40 ||
      sum(therapy_line_df$therapy_line_binary=="first_line") < 10 ||
      sum(therapy_line_df$therapy_line_binary=="later_line") < 10) return(NULL)
  res <- run_coef_model(therapy_line_df, sc, "therapy_line_binary")
  if (is.null(res)) return(NULL)
  res %>% mutate(n_first_line=sum(therapy_line_df$therapy_line_binary=="first_line"),
                 n_later_line=sum(therapy_line_df$therapy_line_binary=="later_line"),
                 beta_later_vs_first=beta,
                 direction=case_when(beta_later_vs_first > 0 ~ "higher_in_later_line",
                                     beta_later_vs_first < 0 ~ "lower_in_later_line",
                                     TRUE ~ "no_direction")) %>%
    select(score, n_model, n_first_line, n_later_line, beta_later_vs_first, conf_low, conf_high, p_value, direction)
})) %>% mutate(FDR_global=p.adjust(p_value, method="BH")) %>% arrange(FDR_global, p_value)

time_df <- cancer_df %>% filter(therapy_status_model=="ongoing_active_treatment", core_covariates_ok, !is.na(time_from_start_days_model))
time_results <- bind_rows(lapply(cp10_score_cols, function(sc) run_numeric_model(time_df, sc, "time_from_start_days_model"))) %>%
  mutate(FDR_global=p.adjust(p_value, method="BH")) %>% arrange(FDR_global, p_value)

# -----------------------------
# 10) Save model outputs
# -----------------------------
write_csv(subgroup_omnibus_results, file.path(out_dir, "SDY2583_CP10_cancer_subgroup_omnibus_results_STEP6.csv"))
write_csv(subgroup_one_vs_rest_results, file.path(out_dir, "SDY2583_CP10_cancer_subgroup_one_vs_rest_results_STEP6.csv"))
write_csv(therapy_status_results, file.path(out_dir, "SDY2583_CP10_therapy_status_ongoing_vs_no_ongoing_results_STEP6.csv"))
write_csv(modality_results, file.path(out_dir, "SDY2583_CP10_active_treatment_modality_results_STEP6.csv"))
write_csv(therapy_line_results, file.path(out_dir, "SDY2583_CP10_therapy_line_later_vs_first_results_STEP6.csv"))
write_csv(time_results, file.path(out_dir, "SDY2583_CP10_time_from_start_results_STEP6.csv"))

save(
  mat_all8, clinical_df, cancer_df, cp10_score_cols,
  selected_matrix_file, updated_matrix_file,
  all8_summary, score_availability_by_panel, variable_sources,
  clinical_model_summary, cancer_subgroup_counts, therapy_status_counts,
  active_exposure_counts, therapy_line_counts, time_from_start_summary,
  subgroup_omnibus_results, subgroup_one_vs_rest_results,
  therapy_status_results, modality_results, therapy_line_results, time_results,
  file=file.path(rdata_out_dir, "SDY2583_CP10_clinical_annotation_STEP6.RData")
)

# -----------------------------
# 11) Console output
# -----------------------------
cat("\n============================================================\n")
cat("SDY2583 CP10 STEP 6 COMPLETE: CLINICAL ANNOTATION\n")
cat("============================================================\n")

cat("\nSelected previous integrated matrix:\n"); print(selected_matrix_file)
cat("\nUpdated ALL8 matrix saved as:\n"); print(updated_matrix_file)

cat("\nClinical variable sources:\n")
print(as.data.frame(variable_sources), row.names=FALSE)

cat("\nALL8 integrated matrix summary:\n")
print(as.data.frame(all8_summary), row.names=FALSE)

cat("\nScore availability by panel after CP10 addition:\n")
print(as.data.frame(score_availability_by_panel), row.names=FALSE)

cat("\nCP10 clinical model summary:\n")
print(as.data.frame(clinical_model_summary), row.names=FALSE)

cat("\nCP10 cancer subgroup counts:\n")
print(as.data.frame(cancer_subgroup_counts), row.names=FALSE)

cat("\nCP10 therapy status counts:\n")
print(as.data.frame(therapy_status_counts), row.names=FALSE)

cat("\nCP10 active-treatment exposure counts:\n")
print(as.data.frame(active_exposure_counts), row.names=FALSE)

cat("\nCP10 therapy line counts:\n")
print(as.data.frame(therapy_line_counts), row.names=FALSE)

cat("\nCP10 time-from-start summary:\n")
print(as.data.frame(time_from_start_summary), row.names=FALSE)

cat("\nTop CP10 cancer-subgroup omnibus results:\n")
print(as.data.frame(subgroup_omnibus_results %>% slice_head(n=10)), row.names=FALSE)

cat("\nTop CP10 cancer-subgroup one-vs-rest results:\n")
print(as.data.frame(subgroup_one_vs_rest_results %>% slice_head(n=30)), row.names=FALSE)

cat("\nTop CP10 therapy-status ongoing-vs-no-ongoing results:\n")
print(as.data.frame(therapy_status_results %>% slice_head(n=10)), row.names=FALSE)

cat("\nTop CP10 active-treatment modality results:\n")
print(as.data.frame(modality_results %>% slice_head(n=40)), row.names=FALSE)

cat("\nTop CP10 therapy-line later-vs-first results:\n")
print(as.data.frame(therapy_line_results %>% slice_head(n=10)), row.names=FALSE)

cat("\nTop CP10 time-from-start results:\n")
print(as.data.frame(time_results %>% slice_head(n=10)), row.names=FALSE)

cat("\nFiles saved in:\n"); print(out_dir)
cat("============================================================\n")
