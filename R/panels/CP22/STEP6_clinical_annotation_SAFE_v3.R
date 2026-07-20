# ============================================================
# SDY2583 CP22
# STEP 6 SAFE: Clinical annotation integration and clinical models
#
# Adds CP22 composite scores to the integrated clinical immune matrix
# and runs CP22 clinical annotation analyses:
#   1) Cancer subgroup
#   2) Therapy status
#   3) Active treatment modality
#   4) Therapy line and time-from-treatment-start
# ============================================================

rm(list = ls())
source(file.path(Sys.getenv("SDY2583_REPO_ROOT", unset = "."), "R", "shared", "bootstrap.R"))

cran_pkgs <- c("dplyr", "readr", "stringr", "tibble", "broom", "purrr", "tidyr")
for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(broom)
  library(purrr)
  library(tidyr)
})

# Namespace safety:
# In some R sessions, filter/select can be masked by stats/flowCore or other packages.
# These aliases force the script to use dplyr verbs.
filter <- dplyr::filter
select <- dplyr::select
mutate <- dplyr::mutate
arrange <- dplyr::arrange
group_by <- dplyr::group_by
ungroup <- dplyr::ungroup
summarise <- dplyr::summarise
count <- dplyr::count
left_join <- dplyr::left_join
distinct <- dplyr::distinct
transmute <- dplyr::transmute
pull <- dplyr::pull
case_when <- dplyr::case_when
n <- dplyr::n
n_distinct <- dplyr::n_distinct

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------

cp22_analysis_dir <- sd_analysis_dir("CP22")
cp22_rdata_dir <- file.path(cp22_analysis_dir, "11_RData")

cp22_step4_for_integration <- file.path(
  cp22_rdata_dir,
  "SDY2583_CP22_STEP4_composite_scores_for_integration.RData"
)

cp22_step4_full <- file.path(
  cp22_rdata_dir,
  "SDY2583_CP22_STEP4_composite_scores.RData"
)

if (file.exists(cp22_step4_for_integration)) {
  cp22_step4_rdata <- cp22_step4_for_integration
} else if (file.exists(cp22_step4_full)) {
  cp22_step4_rdata <- cp22_step4_full
} else {
  stop("CP22 Step 4 RData bulunamadı.")
}

clinical_matrix_dir <- sd_integrated_dir()
if (!dir.exists(clinical_matrix_dir)) {
  stop("Integrated clinical immune matrix klasörü bulunamadı: ", clinical_matrix_dir)
}

out_dir <- file.path(clinical_matrix_dir, "05_CP22_clinical_annotation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

rdata_out_dir <- file.path(clinical_matrix_dir, "RData")
dir.create(rdata_out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Load CP22 composite scores
# ------------------------------------------------------------

load(cp22_step4_rdata)

if (exists("cp22_scores_data")) {
  cp22_data <- cp22_scores_data
} else if (exists("analysis_df")) {
  cp22_data <- analysis_df
} else {
  stop("CP22 Step 4 RData içinde cp22_scores_data veya analysis_df bulunamadı.")
}

if (!("subject_id" %in% names(cp22_data))) stop("CP22 data içinde subject_id yok.")

cp22_score_cols <- grep("^CP22_.*_score$", names(cp22_data), value = TRUE)
if (length(cp22_score_cols) == 0) stop("CP22 composite score kolonu bulunamadı.")

cp22_scores_for_join <- cp22_data %>%
  select(subject_id, all_of(cp22_score_cols)) %>%
  distinct(subject_id, .keep_all = TRUE)

# ------------------------------------------------------------
# Load previous ALL6/ALL7 integrated matrix
# Robust selector:
# The folder also contains summaries/dictionaries/results CSVs.
# Some of those have "ALL6" or "matrix" in the filename but do not
# contain subject_id. Therefore we inspect CSV contents and select
# the real subject-level matrix.
# ------------------------------------------------------------

csv_candidates <- list.files(
  clinical_matrix_dir,
  pattern = "\\.csv$",
  full.names = TRUE,
  recursive = TRUE
)

csv_candidates <- csv_candidates[
  grepl("integrated|matrix|ALL6|ALL7", basename(csv_candidates), ignore.case = TRUE)
]

if (length(csv_candidates) == 0) {
  stop("Integrated matrix CSV adayı bulunamadı: ", clinical_matrix_dir)
}

inspect_candidate_matrix <- function(fp) {

  dat0 <- tryCatch(
    readr::read_csv(fp, show_col_types = FALSE),
    error = function(e) NULL
  )

  if (is.null(dat0)) {
    return(tibble(
      file_path = fp,
      file_name = basename(fp),
      readable = FALSE,
      n_rows = NA_integer_,
      n_cols = NA_integer_,
      has_subject_id = FALSE,
      n_score_cols = NA_integer_,
      has_disease = FALSE,
      has_clinical = FALSE,
      has_cp22 = FALSE,
      priority = -Inf
    ))
  }

  nm <- names(dat0)
  score_cols_tmp <- grep("^CP[0-9]+_.*_score$", nm, value = TRUE)

  has_disease_tmp <- any(c(
    "disease_group_model", "disease_group", "disease_group_clinical",
    "group", "condition"
  ) %in% nm)

  has_clinical_tmp <- any(c(
    "cancer_subgroup", "therapy_status_4level", "chemotherapy",
    "targeted_therapy", "any_immunotherapy", "ici_immunotherapy"
  ) %in% nm)

  has_cp22_tmp <- any(grepl("^CP22_", score_cols_tmp))

  # Subject-level matrix should have subject_id, ~850 rows, many columns,
  # disease information, clinical annotation, and many CP score columns.
  priority_tmp <- 0
  priority_tmp <- priority_tmp + ifelse("subject_id" %in% nm, 1000, 0)
  priority_tmp <- priority_tmp + ifelse(nrow(dat0) >= 800 & nrow(dat0) <= 900, 300, 0)
  priority_tmp <- priority_tmp + ifelse(ncol(dat0) >= 30, 100, 0)
  priority_tmp <- priority_tmp + ifelse(has_disease_tmp, 100, 0)
  priority_tmp <- priority_tmp + ifelse(has_clinical_tmp, 100, 0)
  priority_tmp <- priority_tmp + length(score_cols_tmp)
  priority_tmp <- priority_tmp + ifelse(grepl("ALL6", basename(fp), ignore.case = TRUE), 50, 0)
  priority_tmp <- priority_tmp + ifelse(grepl("ALL7", basename(fp), ignore.case = TRUE), 20, 0)
  priority_tmp <- priority_tmp - ifelse(grepl("summary|dictionary|availability|counts|results|variable_sources|exposure|therapy|subgroup", basename(fp), ignore.case = TRUE), 500, 0)

  tibble(
    file_path = fp,
    file_name = basename(fp),
    readable = TRUE,
    n_rows = nrow(dat0),
    n_cols = ncol(dat0),
    has_subject_id = "subject_id" %in% nm,
    n_score_cols = length(score_cols_tmp),
    has_disease = has_disease_tmp,
    has_clinical = has_clinical_tmp,
    has_cp22 = has_cp22_tmp,
    priority = priority_tmp
  )
}

candidate_inventory <- dplyr::bind_rows(
  lapply(csv_candidates, inspect_candidate_matrix)
) %>%
  arrange(desc(priority), desc(n_score_cols), desc(n_cols))

readr::write_csv(
  candidate_inventory,
  file.path(out_dir, "SDY2583_CP22_integrated_matrix_candidate_inventory_STEP6.csv")
)

valid_matrix_candidates <- candidate_inventory %>%
  filter(
    readable == TRUE,
    has_subject_id == TRUE,
    n_rows >= 800,
    n_rows <= 900,
    has_disease == TRUE,
    has_clinical == TRUE,
    n_score_cols >= 20
  )

if (nrow(valid_matrix_candidates) == 0) {
  cat("\nCandidate inventory:\n")
  print(as.data.frame(candidate_inventory), row.names = FALSE)
  stop("Gerçek subject-level integrated matrix bulunamadı. Candidate inventory dosyasını kontrol et.")
}

# Prefer ALL6 if available, otherwise the highest-priority valid matrix.
valid_all6 <- valid_matrix_candidates %>%
  filter(grepl("ALL6", file_name, ignore.case = TRUE))

if (nrow(valid_all6) > 0) {
  selected_matrix_file <- valid_all6$file_path[1]
} else {
  selected_matrix_file <- valid_matrix_candidates$file_path[1]
}

cat("\nSelected previous integrated matrix:\n")
print(selected_matrix_file)

cat("\nTop integrated matrix candidates:\n")
print(as.data.frame(head(candidate_inventory, 10)), row.names = FALSE)

integrated_prev <- readr::read_csv(selected_matrix_file, show_col_types = FALSE)

if (!("subject_id" %in% names(integrated_prev))) {
  stop("Selected integrated matrix içinde subject_id yok: ", selected_matrix_file)
}

integrated_base <- integrated_prev %>% select(-any_of(cp22_score_cols))
integrated_all7 <- integrated_base %>% left_join(cp22_scores_for_join, by = "subject_id")

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

first_existing <- function(df, possible_cols) {
  hit <- possible_cols[possible_cols %in% names(df)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

as_01 <- function(x) {
  if (is.logical(x)) return(as.integer(x))
  if (is.numeric(x)) return(ifelse(is.na(x), NA_integer_, ifelse(x != 0, 1L, 0L)))
  y <- tolower(trimws(as.character(x)))
  case_when(
    y %in% c("1", "yes", "y", "true", "t", "present", "positive", "var", "evet") ~ 1L,
    y %in% c("0", "no", "n", "false", "f", "absent", "negative", "yok", "hayir", "hayır") ~ 0L,
    TRUE ~ NA_integer_
  )
}

safe_print <- function(x, n = 30) {
  if (is.null(x)) print(NULL) else print(as.data.frame(head(x, n)), row.names = FALSE)
}

add_fdr <- function(df, p_col = "p_value") {
  if (is.null(df) || nrow(df) == 0) return(df)
  df %>% mutate(FDR_global = p.adjust(.data[[p_col]], method = "BH"))
}

# ------------------------------------------------------------
# Harmonize clinical variables
# ------------------------------------------------------------

disease_col <- first_existing(integrated_all7, c("disease_group_clinical", "disease_group", "group", "condition"))
age_col <- first_existing(integrated_all7, c("age_clinical", "age_for_model", "age_years", "age"))
sex_col <- first_existing(integrated_all7, c("sex_clinical", "sex", "gender"))

if (is.na(disease_col)) stop("Disease group column bulunamadı.")
if (is.na(age_col)) stop("Age column bulunamadı.")
if (is.na(sex_col)) stop("Sex column bulunamadı.")
if (!("cancer_subgroup" %in% names(integrated_all7))) stop("cancer_subgroup kolonu yok.")
if (!("therapy_status_4level" %in% names(integrated_all7))) stop("therapy_status_4level kolonu yok.")

clinical_variable_sources <- tibble(
  variable = c("disease_group", "age", "sex"),
  source_column = c(disease_col, age_col, sex_col)
)

integrated_all7 <- integrated_all7 %>%
  mutate(
    disease_group_model = as.character(.data[[disease_col]]),
    disease_group_model = case_when(
      disease_group_model %in% c("Cancer patient", "Cancer", "cancer", "patient", "Patient") ~ "Cancer patient",
      disease_group_model %in% c("Healthy control", "Healthy", "healthy", "control", "Control") ~ "Healthy control",
      str_detect(disease_group_model, regex("cancer|patient", ignore_case = TRUE)) ~ "Cancer patient",
      str_detect(disease_group_model, regex("healthy|control", ignore_case = TRUE)) ~ "Healthy control",
      TRUE ~ disease_group_model
    ),
    disease_group_model = factor(disease_group_model, levels = c("Healthy control", "Cancer patient")),
    age_for_clinical_model = suppressWarnings(as.numeric(.data[[age_col]])),
    age_for_clinical_model = ifelse(age_for_clinical_model < 18 | age_for_clinical_model > 100, NA_real_, age_for_clinical_model),
    sex_for_clinical_model = as.character(.data[[sex_col]]),
    sex_for_clinical_model = case_when(
      sex_for_clinical_model %in% c("Female", "F", "female", "f") ~ "Female",
      sex_for_clinical_model %in% c("Male", "M", "male", "m") ~ "Male",
      TRUE ~ sex_for_clinical_model
    ),
    sex_for_clinical_model = factor(sex_for_clinical_model),
    cancer_subgroup_model = as.character(cancer_subgroup),
    cancer_subgroup_model = ifelse(is.na(cancer_subgroup_model) | cancer_subgroup_model == "", NA_character_, cancer_subgroup_model),
    cancer_subgroup_model = factor(cancer_subgroup_model),
    therapy_status_model = as.character(therapy_status_4level),
    therapy_status_model = ifelse(is.na(therapy_status_model) | therapy_status_model == "", NA_character_, therapy_status_model),
    therapy_status_model = factor(therapy_status_model)
  )

exposure_cols <- intersect(
  c("chemotherapy", "targeted_therapy", "any_immunotherapy", "ici_immunotherapy",
    "endocrine_hormonal", "adc", "experimental", "radiotherapy"),
  names(integrated_all7)
)

for (ec in exposure_cols) {
  integrated_all7[[paste0(ec, "_01")]] <- as_01(integrated_all7[[ec]])
}

therapy_line_col <- first_existing(integrated_all7, c("therapy_line_number", "therapy_line_model", "therapy_line"))
time_col <- first_existing(integrated_all7, c("time_from_start_days", "time_from_treatment_start_days", "time_from_start"))

if (!is.na(therapy_line_col)) {
  integrated_all7 <- integrated_all7 %>%
    mutate(
      therapy_line_number_model = suppressWarnings(as.numeric(.data[[therapy_line_col]])),
      therapy_line_group = case_when(
        therapy_line_number_model == 1 ~ "first_line",
        therapy_line_number_model > 1 ~ "later_line",
        TRUE ~ NA_character_
      ),
      therapy_line_group = factor(therapy_line_group, levels = c("first_line", "later_line"))
    )
} else {
  integrated_all7$therapy_line_number_model <- NA_real_
  integrated_all7$therapy_line_group <- factor(NA_character_, levels = c("first_line", "later_line"))
}

if (!is.na(time_col)) {
  integrated_all7 <- integrated_all7 %>%
    mutate(
      time_from_start_days_model = suppressWarnings(as.numeric(.data[[time_col]])),
      time_from_start_log1p_z = as.numeric(scale(log1p(time_from_start_days_model)))
    )
} else {
  integrated_all7$time_from_start_days_model <- NA_real_
  integrated_all7$time_from_start_log1p_z <- NA_real_
}

# ------------------------------------------------------------
# Save updated ALL7 matrix
# ------------------------------------------------------------

score_cols_all <- grep("^CP[0-9]+_.*_score$", names(integrated_all7), value = TRUE)
panel_from_score <- str_extract(score_cols_all, "^CP[0-9]+")

score_column_dictionary_ALL7 <- tibble(score_column = score_cols_all, panel = panel_from_score) %>%
  arrange(panel, score_column)

score_availability_by_panel_ALL7 <- score_column_dictionary_ALL7 %>%
  group_by(panel) %>%
  summarise(
    n_score_cols = n(),
    n_subjects_with_any_score = sum(rowSums(!is.na(integrated_all7[, score_column, drop = FALSE])) > 0),
    .groups = "drop"
  )

available_panel_count <- apply(!is.na(integrated_all7[, score_cols_all, drop = FALSE]), 1, function(x) {
  length(unique(panel_from_score[x]))
})

integrated_matrix_summary_ALL7 <- tibble(
  n_subjects = nrow(integrated_all7),
  n_unique_subjects = n_distinct(integrated_all7$subject_id),
  n_healthy = sum(integrated_all7$disease_group_model == "Healthy control", na.rm = TRUE),
  n_cancer = sum(integrated_all7$disease_group_model == "Cancer patient", na.rm = TRUE),
  n_with_any_CP22_score = sum(rowSums(!is.na(integrated_all7[, cp22_score_cols, drop = FALSE])) > 0),
  n_with_all_CP22_scores = sum(rowSums(!is.na(integrated_all7[, cp22_score_cols, drop = FALSE])) == length(cp22_score_cols)),
  n_total_score_cols_ALL7 = length(score_cols_all),
  median_available_score_cols = median(rowSums(!is.na(integrated_all7[, score_cols_all, drop = FALSE])), na.rm = TRUE),
  median_available_panels = median(available_panel_count, na.rm = TRUE)
)

all7_csv <- file.path(clinical_matrix_dir, "SDY2583_integrated_clinical_immune_score_matrix_ALL7_with_CP22.csv")
write_csv(integrated_all7, all7_csv)
write_csv(score_column_dictionary_ALL7, file.path(clinical_matrix_dir, "SDY2583_ALL7_score_column_dictionary_with_CP22.csv"))
write_csv(score_availability_by_panel_ALL7, file.path(clinical_matrix_dir, "SDY2583_ALL7_score_availability_by_panel_with_CP22.csv"))
write_csv(integrated_matrix_summary_ALL7, file.path(clinical_matrix_dir, "SDY2583_ALL7_integrated_matrix_summary_with_CP22.csv"))

save(
  integrated_all7, score_column_dictionary_ALL7, score_availability_by_panel_ALL7,
  integrated_matrix_summary_ALL7, selected_matrix_file, cp22_score_cols,
  file = file.path(rdata_out_dir, "SDY2583_integrated_clinical_immune_score_matrix_ALL7_with_CP22.RData")
)

# ------------------------------------------------------------
# CP22 clinical datasets
# ------------------------------------------------------------

cp22_cancer <- integrated_all7 %>% filter(disease_group_model == "Cancer patient")

cp22_cancer_model_ready <- cp22_cancer %>%
  filter(
    !is.na(age_for_clinical_model),
    !is.na(sex_for_clinical_model),
    !is.na(cancer_subgroup_model)
  )

active_treatment_df <- cp22_cancer_model_ready %>%
  filter(therapy_status_model == "ongoing_active_treatment")

# ------------------------------------------------------------
# Model helpers
# ------------------------------------------------------------

run_cancer_subgroup_omnibus <- function(df, score) {
  model_df <- df %>%
    transmute(
      value = suppressWarnings(as.numeric(.data[[score]])),
      cancer_subgroup_model = droplevels(factor(cancer_subgroup_model)),
      age_for_clinical_model = age_for_clinical_model,
      sex_for_clinical_model = droplevels(factor(sex_for_clinical_model))
    ) %>%
    filter(!is.na(value), !is.na(cancer_subgroup_model), !is.na(age_for_clinical_model), !is.na(sex_for_clinical_model))

  if (nrow(model_df) < 50 || nlevels(model_df$cancer_subgroup_model) < 2) return(NULL)

  fit <- tryCatch(lm(value ~ cancer_subgroup_model + age_for_clinical_model + sex_for_clinical_model, data = model_df), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  aa <- tryCatch(anova(fit), error = function(e) NULL)
  if (is.null(aa) || !("cancer_subgroup_model" %in% rownames(aa))) return(NULL)
  row <- aa["cancer_subgroup_model", , drop = FALSE]

  tibble(
    score = score,
    n_model = nrow(model_df),
    n_subgroups = nlevels(model_df$cancer_subgroup_model),
    omnibus_F = as.numeric(row$`F value`[1]),
    omnibus_p = as.numeric(row$`Pr(>F)`[1])
  )
}

run_one_vs_rest <- function(df, score, subgroup) {
  model_df <- df %>%
    transmute(
      value = suppressWarnings(as.numeric(.data[[score]])),
      subgroup_raw = as.character(cancer_subgroup_model),
      age_for_clinical_model = age_for_clinical_model,
      sex_for_clinical_model = droplevels(factor(sex_for_clinical_model))
    ) %>%
    filter(!is.na(value), !is.na(subgroup_raw), !is.na(age_for_clinical_model), !is.na(sex_for_clinical_model)) %>%
    mutate(is_subgroup = factor(ifelse(subgroup_raw == subgroup, "target", "rest"), levels = c("rest", "target")))

  n_target <- sum(model_df$is_subgroup == "target")
  n_rest <- sum(model_df$is_subgroup == "rest")
  if (nrow(model_df) < 50 || n_target < 10 || n_rest < 20) return(NULL)

  fit <- tryCatch(lm(value ~ is_subgroup + age_for_clinical_model + sex_for_clinical_model, data = model_df), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  tt <- tryCatch(broom::tidy(fit, conf.int = TRUE), error = function(e) NULL)
  if (is.null(tt) || !("is_subgrouptarget" %in% tt$term)) return(NULL)
  out <- tt %>% filter(term == "is_subgrouptarget")

  tibble(
    score = score,
    cancer_subgroup = subgroup,
    n_model = nrow(model_df),
    n_subgroup = n_target,
    n_rest = n_rest,
    beta_subgroup_vs_rest = out$estimate[1],
    conf_low = out$conf.low[1],
    conf_high = out$conf.high[1],
    p_value = out$p.value[1],
    direction = case_when(
      out$estimate[1] > 0 ~ "higher_in_subgroup_vs_rest",
      out$estimate[1] < 0 ~ "lower_in_subgroup_vs_rest",
      TRUE ~ "no_direction"
    )
  )
}

run_therapy_status_binary <- function(df, score) {
  model_df <- df %>%
    filter(therapy_status_model %in% c("no_ongoing_therapy", "ongoing_active_treatment")) %>%
    transmute(
      value = suppressWarnings(as.numeric(.data[[score]])),
      therapy_binary = factor(as.character(therapy_status_model), levels = c("no_ongoing_therapy", "ongoing_active_treatment")),
      age_for_clinical_model = age_for_clinical_model,
      sex_for_clinical_model = droplevels(factor(sex_for_clinical_model)),
      cancer_subgroup_model = droplevels(factor(cancer_subgroup_model))
    ) %>%
    filter(!is.na(value), !is.na(therapy_binary), !is.na(age_for_clinical_model), !is.na(sex_for_clinical_model), !is.na(cancer_subgroup_model))

  n_ongoing <- sum(model_df$therapy_binary == "ongoing_active_treatment")
  n_no <- sum(model_df$therapy_binary == "no_ongoing_therapy")
  if (nrow(model_df) < 50 || n_ongoing < 20 || n_no < 20) return(NULL)

  fit <- tryCatch(lm(value ~ therapy_binary + age_for_clinical_model + sex_for_clinical_model + cancer_subgroup_model, data = model_df), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  tt <- tryCatch(broom::tidy(fit, conf.int = TRUE), error = function(e) NULL)
  term <- "therapy_binaryongoing_active_treatment"
  if (is.null(tt) || !(term %in% tt$term)) return(NULL)
  out <- tt %>% filter(term == !!term)

  tibble(
    score = score,
    n_model = nrow(model_df),
    n_ongoing = n_ongoing,
    n_no_ongoing = n_no,
    beta_ongoing_vs_no_ongoing = out$estimate[1],
    conf_low = out$conf.low[1],
    conf_high = out$conf.high[1],
    p_value = out$p.value[1],
    direction = case_when(
      out$estimate[1] > 0 ~ "higher_in_ongoing_active_treatment",
      out$estimate[1] < 0 ~ "lower_in_ongoing_active_treatment",
      TRUE ~ "no_direction"
    )
  )
}

run_binary_exposure_model <- function(df, score, exposure_col, exposure_label) {
  model_df <- df %>%
    transmute(
      value = suppressWarnings(as.numeric(.data[[score]])),
      exposure_num = suppressWarnings(as.numeric(.data[[exposure_col]])),
      age_for_clinical_model = age_for_clinical_model,
      sex_for_clinical_model = droplevels(factor(sex_for_clinical_model)),
      cancer_subgroup_model = droplevels(factor(cancer_subgroup_model))
    ) %>%
    filter(!is.na(value), !is.na(exposure_num), !is.na(age_for_clinical_model), !is.na(sex_for_clinical_model), !is.na(cancer_subgroup_model))

  n_exposed <- sum(model_df$exposure_num == 1)
  n_unexposed <- sum(model_df$exposure_num == 0)
  if (nrow(model_df) < 50 || n_exposed < 10 || n_unexposed < 20) return(NULL)

  fit <- tryCatch(lm(value ~ exposure_num + age_for_clinical_model + sex_for_clinical_model + cancer_subgroup_model, data = model_df), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  tt <- tryCatch(broom::tidy(fit, conf.int = TRUE), error = function(e) NULL)
  if (is.null(tt) || !("exposure_num" %in% tt$term)) return(NULL)
  out <- tt %>% filter(term == "exposure_num")

  tibble(
    score = score,
    exposure = exposure_label,
    n_model = nrow(model_df),
    n_exposed = n_exposed,
    n_unexposed = n_unexposed,
    beta_exposed_vs_unexposed = out$estimate[1],
    conf_low = out$conf.low[1],
    conf_high = out$conf.high[1],
    p_value = out$p.value[1],
    direction = case_when(
      out$estimate[1] > 0 ~ "higher_in_exposed",
      out$estimate[1] < 0 ~ "lower_in_exposed",
      TRUE ~ "no_direction"
    )
  )
}

run_therapy_line_binary <- function(df, score) {
  model_df <- df %>%
    transmute(
      value = suppressWarnings(as.numeric(.data[[score]])),
      therapy_line_group = therapy_line_group,
      age_for_clinical_model = age_for_clinical_model,
      sex_for_clinical_model = droplevels(factor(sex_for_clinical_model)),
      cancer_subgroup_model = droplevels(factor(cancer_subgroup_model))
    ) %>%
    filter(!is.na(value), !is.na(therapy_line_group), !is.na(age_for_clinical_model), !is.na(sex_for_clinical_model), !is.na(cancer_subgroup_model))

  n_first <- sum(model_df$therapy_line_group == "first_line")
  n_later <- sum(model_df$therapy_line_group == "later_line")
  if (nrow(model_df) < 50 || n_first < 10 || n_later < 10) return(NULL)

  fit <- tryCatch(lm(value ~ therapy_line_group + age_for_clinical_model + sex_for_clinical_model + cancer_subgroup_model, data = model_df), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  tt <- tryCatch(broom::tidy(fit, conf.int = TRUE), error = function(e) NULL)
  term <- "therapy_line_grouplater_line"
  if (is.null(tt) || !(term %in% tt$term)) return(NULL)
  out <- tt %>% filter(term == !!term)

  tibble(
    score = score,
    n_model = nrow(model_df),
    n_first_line = n_first,
    n_later_line = n_later,
    beta_later_vs_first = out$estimate[1],
    conf_low = out$conf.low[1],
    conf_high = out$conf.high[1],
    p_value = out$p.value[1],
    direction = case_when(
      out$estimate[1] > 0 ~ "higher_in_later_line",
      out$estimate[1] < 0 ~ "lower_in_later_line",
      TRUE ~ "no_direction"
    )
  )
}

run_time_from_start <- function(df, score) {
  model_df <- df %>%
    transmute(
      value = suppressWarnings(as.numeric(.data[[score]])),
      time_from_start_log1p_z = time_from_start_log1p_z,
      time_from_start_days_model = time_from_start_days_model,
      age_for_clinical_model = age_for_clinical_model,
      sex_for_clinical_model = droplevels(factor(sex_for_clinical_model)),
      cancer_subgroup_model = droplevels(factor(cancer_subgroup_model))
    ) %>%
    filter(!is.na(value), !is.na(time_from_start_log1p_z), !is.na(age_for_clinical_model), !is.na(sex_for_clinical_model), !is.na(cancer_subgroup_model))

  if (nrow(model_df) < 50) return(NULL)

  fit <- tryCatch(lm(value ~ time_from_start_log1p_z + age_for_clinical_model + sex_for_clinical_model + cancer_subgroup_model, data = model_df), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  tt <- tryCatch(broom::tidy(fit, conf.int = TRUE), error = function(e) NULL)
  term <- "time_from_start_log1p_z"
  if (is.null(tt) || !(term %in% tt$term)) return(NULL)
  out <- tt %>% filter(term == !!term)

  tibble(
    score = score,
    n_model = nrow(model_df),
    median_time_days = median(model_df$time_from_start_days_model, na.rm = TRUE),
    beta_per_log_time_z = out$estimate[1],
    conf_low = out$conf.low[1],
    conf_high = out$conf.high[1],
    p_value = out$p.value[1],
    direction = case_when(
      out$estimate[1] > 0 ~ "higher_with_longer_time_from_start",
      out$estimate[1] < 0 ~ "lower_with_longer_time_from_start",
      TRUE ~ "no_direction"
    )
  )
}

# ------------------------------------------------------------
# Summaries
# ------------------------------------------------------------

cp22_clinical_model_summary <- tibble(
  n_cancer_total = nrow(cp22_cancer),
  n_with_valid_age = sum(!is.na(cp22_cancer$age_for_clinical_model)),
  n_with_valid_sex = sum(!is.na(cp22_cancer$sex_for_clinical_model)),
  n_with_cancer_subgroup = sum(!is.na(cp22_cancer$cancer_subgroup_model)),
  n_with_therapy_status = sum(!is.na(cp22_cancer$therapy_status_model)),
  n_with_all_core_covariates = nrow(cp22_cancer_model_ready),
  n_ongoing_active_treatment = sum(cp22_cancer$therapy_status_model == "ongoing_active_treatment", na.rm = TRUE),
  n_no_ongoing_therapy = sum(cp22_cancer$therapy_status_model == "no_ongoing_therapy", na.rm = TRUE),
  n_no_treatment_data = sum(cp22_cancer$therapy_status_model == "no_treatment_data", na.rm = TRUE)
)

cp22_cancer_subgroup_counts <- cp22_cancer %>%
  count(cancer_subgroup_model, name = "n") %>%
  arrange(desc(n))

cp22_therapy_status_counts <- cp22_cancer %>%
  count(therapy_status_model, name = "n") %>%
  arrange(desc(n))

active_exposure_counts <- bind_rows(lapply(exposure_cols, function(ec) {
  x <- active_treatment_df[[paste0(ec, "_01")]]
  n_exposed <- sum(x == 1, na.rm = TRUE)
  n_unexposed <- sum(x == 0, na.rm = TRUE)
  tibble(
    exposure = ec,
    n_exposed = n_exposed,
    n_unexposed = n_unexposed,
    analysis_tier = case_when(
      n_exposed >= 20 & n_unexposed >= 20 ~ "primary_or_exploratory_powered",
      n_exposed >= 10 & n_unexposed >= 20 ~ "exploratory_underpowered",
      TRUE ~ "descriptive_only"
    )
  )
}))

therapy_line_counts <- active_treatment_df %>%
  count(therapy_line_number_model, therapy_line_group, name = "n") %>%
  arrange(therapy_line_number_model)

time_from_start_summary <- active_treatment_df %>%
  summarise(
    n_with_time = sum(!is.na(time_from_start_days_model)),
    min_days = ifelse(n_with_time > 0, min(time_from_start_days_model, na.rm = TRUE), NA_real_),
    q1_days = ifelse(n_with_time > 0, as.numeric(quantile(time_from_start_days_model, 0.25, na.rm = TRUE)), NA_real_),
    median_days = ifelse(n_with_time > 0, median(time_from_start_days_model, na.rm = TRUE), NA_real_),
    mean_days = ifelse(n_with_time > 0, mean(time_from_start_days_model, na.rm = TRUE), NA_real_),
    q3_days = ifelse(n_with_time > 0, as.numeric(quantile(time_from_start_days_model, 0.75, na.rm = TRUE)), NA_real_),
    max_days = ifelse(n_with_time > 0, max(time_from_start_days_model, na.rm = TRUE), NA_real_)
  )

# ------------------------------------------------------------
# Run clinical models
# ------------------------------------------------------------

cancer_subgroup_omnibus_results <- bind_rows(lapply(cp22_score_cols, function(sc) {
  run_cancer_subgroup_omnibus(cp22_cancer_model_ready, sc)
})) %>%
  add_fdr("omnibus_p") %>%
  arrange(FDR_global, omnibus_p)

subgroups <- cp22_cancer_subgroup_counts %>%
  filter(!is.na(cancer_subgroup_model), n >= 10) %>%
  pull(cancer_subgroup_model) %>%
  as.character()

one_vs_rest_results <- bind_rows(lapply(cp22_score_cols, function(sc) {
  bind_rows(lapply(subgroups, function(sg) run_one_vs_rest(cp22_cancer_model_ready, sc, sg)))
}))

if (nrow(one_vs_rest_results) > 0) {
  one_vs_rest_results <- one_vs_rest_results %>%
    group_by(cancer_subgroup) %>%
    mutate(FDR_within_subgroup = p.adjust(p_value, method = "BH")) %>%
    ungroup() %>%
    group_by(score) %>%
    mutate(FDR_within_score = p.adjust(p_value, method = "BH")) %>%
    ungroup() %>%
    mutate(FDR_global = p.adjust(p_value, method = "BH")) %>%
    arrange(FDR_global, p_value)
}

therapy_status_binary_results <- bind_rows(lapply(cp22_score_cols, function(sc) {
  run_therapy_status_binary(cp22_cancer_model_ready, sc)
})) %>%
  add_fdr("p_value") %>%
  arrange(FDR_global, p_value)

active_modality_results <- bind_rows(lapply(cp22_score_cols, function(sc) {
  bind_rows(lapply(exposure_cols, function(ec) {
    run_binary_exposure_model(active_treatment_df, sc, paste0(ec, "_01"), ec)
  }))
}))

if (nrow(active_modality_results) > 0) {
  active_modality_results <- active_modality_results %>%
    left_join(active_exposure_counts, by = "exposure") %>%
    group_by(exposure) %>%
    mutate(FDR_within_exposure = p.adjust(p_value, method = "BH")) %>%
    ungroup() %>%
    group_by(score) %>%
    mutate(FDR_within_score = p.adjust(p_value, method = "BH")) %>%
    ungroup() %>%
    mutate(FDR_global = p.adjust(p_value, method = "BH")) %>%
    arrange(FDR_global, p_value)
}

therapy_line_binary_results <- bind_rows(lapply(cp22_score_cols, function(sc) {
  run_therapy_line_binary(active_treatment_df, sc)
})) %>%
  add_fdr("p_value") %>%
  arrange(FDR_global, p_value)

time_from_start_results <- bind_rows(lapply(cp22_score_cols, function(sc) {
  run_time_from_start(active_treatment_df, sc)
})) %>%
  add_fdr("p_value") %>%
  arrange(FDR_global, p_value)

# ------------------------------------------------------------
# Save outputs
# ------------------------------------------------------------

write_csv(clinical_variable_sources, file.path(out_dir, "SDY2583_CP22_clinical_variable_sources_STEP6.csv"))
write_csv(cp22_clinical_model_summary, file.path(out_dir, "SDY2583_CP22_clinical_model_summary_STEP6.csv"))
write_csv(cp22_cancer_subgroup_counts, file.path(out_dir, "SDY2583_CP22_cancer_subgroup_counts_STEP6.csv"))
write_csv(cp22_therapy_status_counts, file.path(out_dir, "SDY2583_CP22_therapy_status_counts_STEP6.csv"))
write_csv(active_exposure_counts, file.path(out_dir, "SDY2583_CP22_active_treatment_exposure_counts_STEP6.csv"))
write_csv(therapy_line_counts, file.path(out_dir, "SDY2583_CP22_therapy_line_counts_STEP6.csv"))
write_csv(time_from_start_summary, file.path(out_dir, "SDY2583_CP22_time_from_start_summary_STEP6.csv"))
write_csv(cancer_subgroup_omnibus_results, file.path(out_dir, "SDY2583_CP22_cancer_subgroup_omnibus_results_STEP6.csv"))
write_csv(one_vs_rest_results, file.path(out_dir, "SDY2583_CP22_cancer_subgroup_one_vs_rest_results_STEP6.csv"))
write_csv(therapy_status_binary_results, file.path(out_dir, "SDY2583_CP22_therapy_status_ongoing_vs_no_ongoing_results_STEP6.csv"))
write_csv(active_modality_results, file.path(out_dir, "SDY2583_CP22_active_treatment_modality_results_STEP6.csv"))
write_csv(therapy_line_binary_results, file.path(out_dir, "SDY2583_CP22_therapy_line_later_vs_first_results_STEP6.csv"))
write_csv(time_from_start_results, file.path(out_dir, "SDY2583_CP22_time_from_start_results_STEP6.csv"))

save(
  integrated_all7, cp22_scores_for_join, cp22_score_cols,
  cp22_cancer, cp22_cancer_model_ready, active_treatment_df,
  clinical_variable_sources, cp22_clinical_model_summary,
  cp22_cancer_subgroup_counts, cp22_therapy_status_counts,
  active_exposure_counts, therapy_line_counts, time_from_start_summary,
  cancer_subgroup_omnibus_results, one_vs_rest_results,
  therapy_status_binary_results, active_modality_results,
  therapy_line_binary_results, time_from_start_results,
  file = file.path(rdata_out_dir, "SDY2583_CP22_clinical_annotation_STEP6.RData")
)

# ------------------------------------------------------------
# Console output
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SDY2583 CP22 STEP 6 COMPLETE: CLINICAL ANNOTATION\n")
cat("============================================================\n")

cat("\nSelected previous integrated matrix:\n")
print(selected_matrix_file)

cat("\nCP22 Step 4 source:\n")
print(cp22_step4_rdata)

cat("\nClinical variable sources:\n")
safe_print(clinical_variable_sources, 20)

cat("\nALL7 integrated matrix summary:\n")
safe_print(integrated_matrix_summary_ALL7, 20)

cat("\nScore availability by panel after CP22 addition:\n")
safe_print(score_availability_by_panel_ALL7, 30)

cat("\nCP22 clinical model summary:\n")
safe_print(cp22_clinical_model_summary, 20)

cat("\nCP22 cancer subgroup counts:\n")
safe_print(cp22_cancer_subgroup_counts, 30)

cat("\nCP22 therapy status counts:\n")
safe_print(cp22_therapy_status_counts, 20)

cat("\nCP22 active-treatment exposure counts:\n")
safe_print(active_exposure_counts, 20)

cat("\nCP22 therapy line counts:\n")
safe_print(therapy_line_counts, 30)

cat("\nCP22 time-from-start summary:\n")
safe_print(time_from_start_summary, 20)

cat("\nTop CP22 cancer-subgroup omnibus results:\n")
safe_print(cancer_subgroup_omnibus_results, 20)

cat("\nTop CP22 cancer-subgroup one-vs-rest results:\n")
safe_print(one_vs_rest_results, 30)

cat("\nTop CP22 therapy-status ongoing-vs-no-ongoing results:\n")
safe_print(therapy_status_binary_results, 20)

cat("\nTop CP22 active-treatment modality results:\n")
safe_print(active_modality_results, 30)

cat("\nTop CP22 therapy-line later-vs-first results:\n")
safe_print(therapy_line_binary_results, 20)

cat("\nTop CP22 time-from-start results:\n")
safe_print(time_from_start_results, 20)

cat("\nFiles saved in:\n")
print(out_dir)

cat("\nUpdated ALL7 matrix saved as:\n")
print(all7_csv)

cat("============================================================\n")
