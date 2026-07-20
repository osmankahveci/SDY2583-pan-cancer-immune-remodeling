# ============================================================
# SDY2583 CP16
# STEP 6 SAFE: Clinical annotation and ALL9 matrix integration
#
# Input:
#   1) Previous integrated clinical immune matrix:
#      data/derived/clinical_integration/
#        SDY2583_integrated_clinical_immune_score_matrix_ALL8_with_CP10.csv
#
#   2) CP16 composite scores:
#      outputs/CP16/11_RData/
#        SDY2583_CP16_STEP4_composite_scores.RData
#
# Output:
#   data/derived/clinical_integration/
#      SDY2583_integrated_clinical_immune_score_matrix_ALL9_with_CP16.csv
#
#   data/derived/clinical_integration/07_CP16_clinical_annotation
#
# Clinical models:
#   - Cancer subgroup omnibus
#   - Cancer subgroup one-vs-rest
#   - Ongoing active treatment vs no ongoing therapy
#   - Active-treatment modality exposure models
#   - Later-line vs first-line therapy
#   - Time from treatment start
#
# CP16-specific caution:
#   - These analyses are secondary clinical annotations.
#   - Do not infer treatment resistance, response, progression, PFS, or OS.
#   - Use "treatment context" and "tumor-context variation".
# ============================================================

rm(list = ls())
source(file.path(Sys.getenv("SDY2583_REPO_ROOT", unset = "."), "R", "shared", "bootstrap.R"))

# ------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------

cran_pkgs <- c("dplyr", "readr", "stringr", "tibble", "broom", "tidyr", "purrr")

for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(broom)
  library(tidyr)
  library(purrr)
})

# Namespace safety
filter <- dplyr::filter
select <- dplyr::select
mutate <- dplyr::mutate
arrange <- dplyr::arrange
summarise <- dplyr::summarise
group_by <- dplyr::group_by
ungroup <- dplyr::ungroup
count <- dplyr::count
distinct <- dplyr::distinct
case_when <- dplyr::case_when
bind_rows <- dplyr::bind_rows
left_join <- dplyr::left_join
n_distinct <- dplyr::n_distinct
n <- dplyr::n

# ------------------------------------------------------------
# 2. Paths
# ------------------------------------------------------------

matrix_dir <- sd_integrated_dir()
cp16_dir <- sd_analysis_dir("CP16")
cp16_rdata_dir <- file.path(cp16_dir, "11_RData")

out_dir <- file.path(matrix_dir, "07_CP16_clinical_annotation")
rdata_out_dir <- file.path(matrix_dir, "RData")

dir.create(matrix_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rdata_out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 3. Load previous integrated matrix
# ------------------------------------------------------------

preferred_matrix <- file.path(matrix_dir, "SDY2583_integrated_clinical_immune_score_matrix_ALL8_with_CP10.csv")

if (file.exists(preferred_matrix)) {
  selected_previous_matrix <- preferred_matrix
} else {
  all_candidates <- list.files(
    matrix_dir,
    pattern = "\\.csv$",
    full.names = TRUE,
    recursive = FALSE
  )

  all_candidates <- all_candidates[
    grepl("integrated_clinical_immune_score_matrix|integrated.*matrix|ALL", basename(all_candidates), ignore.case = TRUE)
  ]

  all_candidates <- all_candidates[
    !grepl("summary|counts|annotation|subgroup|therapy|manifest|dictionary|results", basename(all_candidates), ignore.case = TRUE)
  ]

  if (length(all_candidates) == 0) {
    stop("Önceki integrated matrix CSV bulunamadı.")
  }

  # Prefer ALL8, then highest ALL number.
  candidate_rank <- tibble(file = all_candidates) %>%
    mutate(
      base = basename(file),
      all_number = suppressWarnings(as.numeric(str_match(base, "ALL([0-9]+)")[, 2])),
      all_number = ifelse(is.na(all_number), -Inf, all_number),
      has_cp10 = grepl("CP10", base, ignore.case = TRUE)
    ) %>%
    arrange(desc(has_cp10), desc(all_number), desc(file.info(file)$mtime))

  selected_previous_matrix <- candidate_rank$file[1]
}

integrated_prev <- readr::read_csv(selected_previous_matrix, show_col_types = FALSE)

if (!("subject_id" %in% names(integrated_prev))) {
  stop("Integrated matrix içinde subject_id yok.")
}

# ------------------------------------------------------------
# 4. Load CP16 composite scores
# ------------------------------------------------------------

step4_rdata <- file.path(cp16_rdata_dir, "SDY2583_CP16_STEP4_composite_scores.RData")

if (!file.exists(step4_rdata)) {
  stop("CP16 Step 4 RData bulunamadı: ", step4_rdata)
}

load(step4_rdata)

if (!exists("cp16_scores_for_integration")) {
  if (exists("cp16_scores_data") && exists("composite_score_cols")) {
    cp16_scores_for_integration <- cp16_scores_data %>%
      select(subject_id, all_of(composite_score_cols))
  } else {
    stop("CP16 score object bulunamadı.")
  }
}

if (!exists("composite_score_cols")) {
  composite_score_cols <- grep("^CP16_.*_score$", names(cp16_scores_for_integration), value = TRUE)
}

cp16_score_cols <- composite_score_cols
cp16_score_cols <- cp16_score_cols[cp16_score_cols %in% names(cp16_scores_for_integration)]

cp16_scores_for_integration <- cp16_scores_for_integration %>%
  mutate(subject_id = as.character(subject_id)) %>%
  select(subject_id, all_of(cp16_score_cols)) %>%
  distinct(subject_id, .keep_all = TRUE)

# Remove any prior CP16 score columns before re-merging.
integrated_prev_clean <- integrated_prev %>%
  mutate(subject_id = as.character(subject_id)) %>%
  select(-any_of(cp16_score_cols))

integrated_all9 <- integrated_prev_clean %>%
  left_join(cp16_scores_for_integration, by = "subject_id")

all_score_cols <- grep("^CP[0-9]+_.*_score$", names(integrated_all9), value = TRUE)

integrated_all9 <- integrated_all9 %>%
  mutate(
    n_available_score_cols_ALL9 = rowSums(!is.na(across(all_of(all_score_cols)))),
    has_any_CP16_score = rowSums(!is.na(across(all_of(cp16_score_cols)))) > 0,
    has_all_CP16_scores = rowSums(!is.na(across(all_of(cp16_score_cols)))) == length(cp16_score_cols)
  )

updated_matrix_file <- file.path(matrix_dir, "SDY2583_integrated_clinical_immune_score_matrix_ALL9_with_CP16.csv")
write_csv(integrated_all9, updated_matrix_file)

# ------------------------------------------------------------
# 5. Clinical variable standardization
# ------------------------------------------------------------

first_existing_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

standardize_disease <- function(x) {
  y <- as.character(x)
  case_when(
    y %in% c("Healthy control", "Healthy", "healthy", "Control", "control", "HC") ~ "Healthy control",
    y %in% c("Cancer patient", "Cancer", "cancer", "Patient", "patient") ~ "Cancer patient",
    str_detect(y, regex("healthy|control", ignore_case = TRUE)) ~ "Healthy control",
    str_detect(y, regex("cancer|patient|tumou?r", ignore_case = TRUE)) ~ "Cancer patient",
    TRUE ~ y
  )
}

standardize_sex <- function(x) {
  y <- as.character(x)
  case_when(
    y %in% c("Female", "F", "female", "f") ~ "Female",
    y %in% c("Male", "M", "male", "m") ~ "Male",
    is.na(y) | y == "" | y == "NA" ~ NA_character_,
    TRUE ~ y
  )
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

disease_col <- first_existing_col(integrated_all9, c("disease_group_model", "disease_group_clinical", "disease_group", "group", "condition"))
age_col <- first_existing_col(integrated_all9, c("age_for_model", "age_clinical", "age_years", "age"))
sex_col <- first_existing_col(integrated_all9, c("sex_for_clinical_model", "sex_clinical", "sex", "gender"))
cancer_subgroup_col <- first_existing_col(integrated_all9, c("cancer_subgroup_model", "cancer_subgroup", "cancer_type_model", "cancer_type", "tumor_type"))
therapy_status_col <- first_existing_col(integrated_all9, c("therapy_status_model", "therapy_status_4level", "therapy_status", "treatment_status"))
therapy_line_col <- first_existing_col(integrated_all9, c("therapy_line_number_model", "therapy_line_number", "treatment_line_number", "line_of_therapy"))
time_col <- first_existing_col(integrated_all9, c("time_from_start_days_model", "time_from_start_days", "time_from_treatment_start_days", "days_from_treatment_start"))

if (is.na(disease_col)) stop("Disease column bulunamadı.")
if (is.na(age_col)) stop("Age column bulunamadı.")
if (is.na(sex_col)) stop("Sex column bulunamadı.")
if (is.na(cancer_subgroup_col)) warning("Cancer subgroup column bulunamadı.")
if (is.na(therapy_status_col)) warning("Therapy status column bulunamadı.")

clinical_df <- integrated_all9 %>%
  mutate(
    disease_group_model = standardize_disease(.data[[disease_col]]),
    disease_group_model = factor(disease_group_model, levels = c("Healthy control", "Cancer patient")),
    age_clinical_model = suppressWarnings(as.numeric(.data[[age_col]])),
    age_clinical_model = ifelse(age_clinical_model >= 18 & age_clinical_model <= 100, age_clinical_model, NA_real_),
    sex_for_clinical_model = standardize_sex(.data[[sex_col]]),
    sex_for_clinical_model = factor(sex_for_clinical_model)
  )

if (!is.na(cancer_subgroup_col)) {
  clinical_df <- clinical_df %>%
    mutate(cancer_subgroup_model = as.character(.data[[cancer_subgroup_col]]))
} else {
  clinical_df$cancer_subgroup_model <- NA_character_
}

if (!is.na(therapy_status_col)) {
  clinical_df <- clinical_df %>%
    mutate(therapy_status_model = as.character(.data[[therapy_status_col]]))
} else {
  clinical_df$therapy_status_model <- NA_character_
}

if (!is.na(therapy_line_col)) {
  clinical_df <- clinical_df %>%
    mutate(therapy_line_number_model = suppressWarnings(as.numeric(.data[[therapy_line_col]])))
} else {
  clinical_df$therapy_line_number_model <- NA_real_
}

if (!is.na(time_col)) {
  clinical_df <- clinical_df %>%
    mutate(time_from_start_days_model = suppressWarnings(as.numeric(.data[[time_col]])))
} else {
  clinical_df$time_from_start_days_model <- NA_real_
}

exposure_vars <- c(
  "chemotherapy",
  "targeted_therapy",
  "any_immunotherapy",
  "ici_immunotherapy",
  "endocrine_hormonal",
  "adc",
  "experimental",
  "radiotherapy"
)

for (ee in exposure_vars) {
  if (ee %in% names(clinical_df)) {
    clinical_df[[ee]] <- as_01(clinical_df[[ee]])
  } else {
    clinical_df[[ee]] <- NA_integer_
  }
}

clinical_df <- clinical_df %>%
  mutate(
    therapy_status_binary = case_when(
      therapy_status_model == "ongoing_active_treatment" ~ "ongoing_active_treatment",
      therapy_status_model == "no_ongoing_therapy" ~ "no_ongoing_therapy",
      TRUE ~ NA_character_
    ),
    therapy_status_binary = factor(therapy_status_binary, levels = c("no_ongoing_therapy", "ongoing_active_treatment")),
    therapy_line_group = case_when(
      therapy_line_number_model == 1 ~ "first_line",
      therapy_line_number_model >= 2 ~ "later_line",
      TRUE ~ NA_character_
    ),
    therapy_line_group = factor(therapy_line_group, levels = c("first_line", "later_line")),
    log_time_from_start_z = ifelse(
      !is.na(time_from_start_days_model) & time_from_start_days_model >= 0,
      log1p(time_from_start_days_model),
      NA_real_
    ),
    log_time_from_start_z = as.numeric(scale(log_time_from_start_z))
  )

# ------------------------------------------------------------
# 6. Matrix summaries
# ------------------------------------------------------------

score_availability_by_panel <- tibble(
  score_col = all_score_cols,
  panel = str_extract(score_col, "^CP[0-9]+")
) %>%
  group_by(panel) %>%
  summarise(
    n_score_cols = n(),
    n_subjects_with_any_score = sum(rowSums(!is.na(integrated_all9[, score_col, drop = FALSE])) > 0),
    .groups = "drop"
  ) %>%
  arrange(panel)

all9_matrix_summary <- tibble(
  n_subjects = nrow(integrated_all9),
  n_unique_subjects = n_distinct(integrated_all9$subject_id),
  n_with_any_CP16_score = sum(integrated_all9$has_any_CP16_score, na.rm = TRUE),
  n_with_all_CP16_scores = sum(integrated_all9$has_all_CP16_scores, na.rm = TRUE),
  n_total_score_cols_ALL9 = length(all_score_cols),
  median_available_score_cols = median(integrated_all9$n_available_score_cols_ALL9, na.rm = TRUE)
)

clinical_variable_sources <- tibble(
  variable = c("disease_group", "age", "sex", "cancer_subgroup", "therapy_status", "therapy_line_number", "time_from_start_days"),
  source_column = c(disease_col, age_col, sex_col, cancer_subgroup_col, therapy_status_col, therapy_line_col, time_col)
)

# ------------------------------------------------------------
# 7. Clinical model helper
# ------------------------------------------------------------

fmt_score_label <- function(x) x

run_lm_safe <- function(model_df, formula_obj, term_name, score_name, model_type, extra = list()) {

  fit <- tryCatch(lm(formula_obj, data = model_df), error = function(e) NULL)
  if (is.null(fit)) return(NULL)

  tt <- tryCatch(broom::tidy(fit, conf.int = TRUE), error = function(e) NULL)
  if (is.null(tt)) return(NULL)
  if (!(term_name %in% tt$term)) return(NULL)

  out <- tt %>% filter(term == term_name)

  bind_cols(
    tibble(
      score = score_name,
      model_type = model_type,
      n_model = nrow(model_df),
      beta = out$estimate[1],
      conf_low = out$conf.low[1],
      conf_high = out$conf.high[1],
      p_value = out$p.value[1],
      direction = case_when(
        out$estimate[1] > 0 ~ "higher",
        out$estimate[1] < 0 ~ "lower",
        TRUE ~ "no_direction"
      )
    ),
    as_tibble(extra)
  )
}

# ------------------------------------------------------------
# 8. Cancer-only data
# ------------------------------------------------------------

cancer_df <- clinical_df %>%
  filter(disease_group_model == "Cancer patient") %>%
  filter(!is.na(age_clinical_model), !is.na(sex_for_clinical_model)) %>%
  mutate(
    cancer_subgroup_model = ifelse(is.na(cancer_subgroup_model) | cancer_subgroup_model == "", NA_character_, cancer_subgroup_model),
    cancer_subgroup_model = factor(cancer_subgroup_model)
  )

cp16_clinical_model_summary <- tibble(
  n_cancer_total = sum(clinical_df$disease_group_model == "Cancer patient", na.rm = TRUE),
  n_with_valid_age = sum(clinical_df$disease_group_model == "Cancer patient" & !is.na(clinical_df$age_clinical_model), na.rm = TRUE),
  n_with_valid_sex = sum(clinical_df$disease_group_model == "Cancer patient" & !is.na(clinical_df$sex_for_clinical_model), na.rm = TRUE),
  n_with_cancer_subgroup = sum(clinical_df$disease_group_model == "Cancer patient" & !is.na(clinical_df$cancer_subgroup_model), na.rm = TRUE),
  n_with_therapy_status = sum(clinical_df$disease_group_model == "Cancer patient" & !is.na(clinical_df$therapy_status_model), na.rm = TRUE),
  n_with_all_core_covariates = nrow(cancer_df),
  n_ongoing_active_treatment = sum(cancer_df$therapy_status_binary == "ongoing_active_treatment", na.rm = TRUE),
  n_no_ongoing_therapy = sum(cancer_df$therapy_status_binary == "no_ongoing_therapy", na.rm = TRUE),
  n_no_treatment_data = sum(cancer_df$therapy_status_model == "no_treatment_data", na.rm = TRUE)
)

cp16_cancer_subgroup_counts <- cancer_df %>%
  count(cancer_subgroup_model, name = "n") %>%
  arrange(desc(n))

cp16_therapy_status_counts <- cancer_df %>%
  count(therapy_status_model, name = "n") %>%
  arrange(desc(n))

active_treatment_df <- cancer_df %>%
  filter(therapy_status_model == "ongoing_active_treatment")

cp16_active_treatment_exposure_counts <- bind_rows(lapply(exposure_vars, function(ee) {
  tibble(
    exposure = ee,
    n_exposed = sum(active_treatment_df[[ee]] == 1, na.rm = TRUE),
    n_unexposed = sum(active_treatment_df[[ee]] == 0, na.rm = TRUE),
    analysis_tier = case_when(
      sum(active_treatment_df[[ee]] == 1, na.rm = TRUE) >= 20 &
        sum(active_treatment_df[[ee]] == 0, na.rm = TRUE) >= 20 ~ "primary_or_exploratory_powered",
      sum(active_treatment_df[[ee]] == 1, na.rm = TRUE) >= 10 &
        sum(active_treatment_df[[ee]] == 0, na.rm = TRUE) >= 20 ~ "exploratory_underpowered",
      TRUE ~ "descriptive_only"
    )
  )
}))

cp16_therapy_line_counts <- cancer_df %>%
  count(therapy_line_number_model, therapy_line_group, name = "n") %>%
  arrange(therapy_line_number_model)

cp16_time_from_start_summary <- cancer_df %>%
  summarise(
    n_with_time = sum(!is.na(time_from_start_days_model)),
    min_days = min(time_from_start_days_model, na.rm = TRUE),
    q1_days = as.numeric(quantile(time_from_start_days_model, 0.25, na.rm = TRUE)),
    median_days = median(time_from_start_days_model, na.rm = TRUE),
    mean_days = mean(time_from_start_days_model, na.rm = TRUE),
    q3_days = as.numeric(quantile(time_from_start_days_model, 0.75, na.rm = TRUE)),
    max_days = max(time_from_start_days_model, na.rm = TRUE)
  )

# ------------------------------------------------------------
# 9. Cancer subgroup omnibus models
# ------------------------------------------------------------

run_subgroup_omnibus <- function(score_name) {

  df <- cancer_df %>%
    transmute(
      value = suppressWarnings(as.numeric(.data[[score_name]])),
      age_clinical_model,
      sex_for_clinical_model,
      cancer_subgroup_model
    ) %>%
    filter(!is.na(value), !is.na(age_clinical_model), !is.na(sex_for_clinical_model), !is.na(cancer_subgroup_model)) %>%
    mutate(
      cancer_subgroup_model = droplevels(factor(cancer_subgroup_model)),
      sex_for_clinical_model = droplevels(factor(sex_for_clinical_model))
    )

  if (nrow(df) < 50 || nlevels(df$cancer_subgroup_model) < 2) return(NULL)

  fit_full <- tryCatch(lm(value ~ cancer_subgroup_model + age_clinical_model + sex_for_clinical_model, data = df), error = function(e) NULL)
  fit_reduced <- tryCatch(lm(value ~ age_clinical_model + sex_for_clinical_model, data = df), error = function(e) NULL)
  if (is.null(fit_full) || is.null(fit_reduced)) return(NULL)

  aa <- tryCatch(anova(fit_reduced, fit_full), error = function(e) NULL)
  if (is.null(aa) || nrow(aa) < 2) return(NULL)

  tibble(
    score = score_name,
    n_model = nrow(df),
    n_subgroups = nlevels(df$cancer_subgroup_model),
    omnibus_F = aa$F[2],
    omnibus_p = aa$`Pr(>F)`[2]
  )
}

cancer_subgroup_omnibus_results <- bind_rows(lapply(cp16_score_cols, run_subgroup_omnibus)) %>%
  mutate(FDR_global = p.adjust(omnibus_p, method = "BH")) %>%
  arrange(FDR_global, omnibus_p)

# ------------------------------------------------------------
# 10. Cancer subgroup one-vs-rest
# ------------------------------------------------------------

run_subgroup_one_vs_rest <- function(score_name, subgroup_name) {

  df <- cancer_df %>%
    transmute(
      value = suppressWarnings(as.numeric(.data[[score_name]])),
      age_clinical_model,
      sex_for_clinical_model,
      cancer_subgroup_model = as.character(cancer_subgroup_model)
    ) %>%
    filter(!is.na(value), !is.na(age_clinical_model), !is.na(sex_for_clinical_model), !is.na(cancer_subgroup_model)) %>%
    mutate(
      subgroup_vs_rest = ifelse(cancer_subgroup_model == subgroup_name, "subgroup", "rest"),
      subgroup_vs_rest = factor(subgroup_vs_rest, levels = c("rest", "subgroup")),
      sex_for_clinical_model = droplevels(factor(sex_for_clinical_model))
    )

  n_subgroup <- sum(df$subgroup_vs_rest == "subgroup")
  n_rest <- sum(df$subgroup_vs_rest == "rest")

  if (nrow(df) < 50 || n_subgroup < 10 || n_rest < 20) return(NULL)

  run_lm_safe(
    df,
    value ~ subgroup_vs_rest + age_clinical_model + sex_for_clinical_model,
    "subgroup_vs_restsubgroup",
    score_name,
    "cancer_subgroup_one_vs_rest",
    extra = list(
      cancer_subgroup = subgroup_name,
      n_subgroup = n_subgroup,
      n_rest = n_rest
    )
  ) %>%
    mutate(
      direction = case_when(
        beta > 0 ~ "higher_in_subgroup_vs_rest",
        beta < 0 ~ "lower_in_subgroup_vs_rest",
        TRUE ~ "no_direction"
      )
    )
}

subgroups_to_test <- cancer_df %>%
  filter(!is.na(cancer_subgroup_model)) %>%
  count(cancer_subgroup_model, name = "n") %>%
  filter(n >= 10) %>%
  pull(cancer_subgroup_model) %>%
  as.character()

cancer_subgroup_one_vs_rest_results <- bind_rows(lapply(cp16_score_cols, function(sc) {
  bind_rows(lapply(subgroups_to_test, function(sg) {
    run_subgroup_one_vs_rest(sc, sg)
  }))
})) %>%
  group_by(score) %>%
  mutate(FDR_within_score = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  group_by(cancer_subgroup) %>%
  mutate(FDR_within_subgroup = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  mutate(FDR_global = p.adjust(p_value, method = "BH")) %>%
  arrange(FDR_global, p_value)

# ------------------------------------------------------------
# 11. Therapy status: ongoing vs no ongoing
# ------------------------------------------------------------

therapy_status_binary_results <- bind_rows(lapply(cp16_score_cols, function(sc) {

  df <- cancer_df %>%
    transmute(
      value = suppressWarnings(as.numeric(.data[[sc]])),
      age_clinical_model,
      sex_for_clinical_model,
      cancer_subgroup_model,
      therapy_status_binary
    ) %>%
    filter(
      !is.na(value),
      !is.na(age_clinical_model),
      !is.na(sex_for_clinical_model),
      !is.na(cancer_subgroup_model),
      !is.na(therapy_status_binary)
    ) %>%
    mutate(
      therapy_status_binary = factor(therapy_status_binary, levels = c("no_ongoing_therapy", "ongoing_active_treatment")),
      cancer_subgroup_model = droplevels(factor(cancer_subgroup_model)),
      sex_for_clinical_model = droplevels(factor(sex_for_clinical_model))
    )

  if (nrow(df) < 50 ||
      sum(df$therapy_status_binary == "ongoing_active_treatment") < 20 ||
      sum(df$therapy_status_binary == "no_ongoing_therapy") < 20) return(NULL)

  run_lm_safe(
    df,
    value ~ therapy_status_binary + age_clinical_model + sex_for_clinical_model + cancer_subgroup_model,
    "therapy_status_binaryongoing_active_treatment",
    sc,
    "therapy_status_ongoing_vs_no_ongoing",
    extra = list(
      n_ongoing = sum(df$therapy_status_binary == "ongoing_active_treatment"),
      n_no_ongoing = sum(df$therapy_status_binary == "no_ongoing_therapy")
    )
  ) %>%
    mutate(
      direction = case_when(
        beta > 0 ~ "higher_in_ongoing_active_treatment",
        beta < 0 ~ "lower_in_ongoing_active_treatment",
        TRUE ~ "no_direction"
      )
    )
})) %>%
  mutate(FDR_global = p.adjust(p_value, method = "BH")) %>%
  arrange(FDR_global, p_value)

# ------------------------------------------------------------
# 12. Active-treatment modality models
# ------------------------------------------------------------

run_modality_model <- function(score_name, exposure_var) {

  df <- active_treatment_df %>%
    transmute(
      value = suppressWarnings(as.numeric(.data[[score_name]])),
      exposure = .data[[exposure_var]],
      age_clinical_model,
      sex_for_clinical_model,
      cancer_subgroup_model
    ) %>%
    filter(
      !is.na(value),
      !is.na(exposure),
      !is.na(age_clinical_model),
      !is.na(sex_for_clinical_model),
      !is.na(cancer_subgroup_model)
    ) %>%
    mutate(
      exposure = as.integer(exposure),
      cancer_subgroup_model = droplevels(factor(cancer_subgroup_model)),
      sex_for_clinical_model = droplevels(factor(sex_for_clinical_model))
    )

  n_exposed <- sum(df$exposure == 1)
  n_unexposed <- sum(df$exposure == 0)

  tier <- case_when(
    n_exposed >= 20 & n_unexposed >= 20 ~ "primary_or_exploratory_powered",
    n_exposed >= 10 & n_unexposed >= 20 ~ "exploratory_underpowered",
    TRUE ~ "descriptive_only"
  )

  if (nrow(df) < 30 || n_exposed < 5 || n_unexposed < 10) return(NULL)

  res <- run_lm_safe(
    df,
    value ~ exposure + age_clinical_model + sex_for_clinical_model + cancer_subgroup_model,
    "exposure",
    score_name,
    "active_treatment_modality",
    extra = list(
      exposure = exposure_var,
      n_exposed = n_exposed,
      n_unexposed = n_unexposed,
      analysis_tier = tier
    )
  )

  if (is.null(res)) return(NULL)

  res %>%
    mutate(
      direction = case_when(
        beta > 0 ~ "higher_in_exposed",
        beta < 0 ~ "lower_in_exposed",
        TRUE ~ "no_direction"
      )
    )
}

active_modality_results <- bind_rows(lapply(cp16_score_cols, function(sc) {
  bind_rows(lapply(exposure_vars, function(ee) {
    run_modality_model(sc, ee)
  }))
})) %>%
  group_by(exposure) %>%
  mutate(FDR_within_exposure = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  group_by(score) %>%
  mutate(FDR_within_score = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  mutate(FDR_global = p.adjust(p_value, method = "BH")) %>%
  arrange(FDR_global, p_value)

# ------------------------------------------------------------
# 13. Therapy line: later vs first
# ------------------------------------------------------------

therapy_line_binary_results <- bind_rows(lapply(cp16_score_cols, function(sc) {

  df <- cancer_df %>%
    transmute(
      value = suppressWarnings(as.numeric(.data[[sc]])),
      age_clinical_model,
      sex_for_clinical_model,
      cancer_subgroup_model,
      therapy_line_group
    ) %>%
    filter(
      !is.na(value),
      !is.na(age_clinical_model),
      !is.na(sex_for_clinical_model),
      !is.na(cancer_subgroup_model),
      !is.na(therapy_line_group)
    ) %>%
    mutate(
      therapy_line_group = factor(therapy_line_group, levels = c("first_line", "later_line")),
      cancer_subgroup_model = droplevels(factor(cancer_subgroup_model)),
      sex_for_clinical_model = droplevels(factor(sex_for_clinical_model))
    )

  if (nrow(df) < 40 ||
      sum(df$therapy_line_group == "first_line") < 15 ||
      sum(df$therapy_line_group == "later_line") < 15) return(NULL)

  run_lm_safe(
    df,
    value ~ therapy_line_group + age_clinical_model + sex_for_clinical_model + cancer_subgroup_model,
    "therapy_line_grouplater_line",
    sc,
    "therapy_line_later_vs_first",
    extra = list(
      n_first_line = sum(df$therapy_line_group == "first_line"),
      n_later_line = sum(df$therapy_line_group == "later_line")
    )
  ) %>%
    mutate(
      direction = case_when(
        beta > 0 ~ "higher_in_later_line",
        beta < 0 ~ "lower_in_later_line",
        TRUE ~ "no_direction"
      )
    )
})) %>%
  mutate(FDR_global = p.adjust(p_value, method = "BH")) %>%
  arrange(FDR_global, p_value)

# ------------------------------------------------------------
# 14. Time from treatment start
# ------------------------------------------------------------

time_from_start_results <- bind_rows(lapply(cp16_score_cols, function(sc) {

  df <- active_treatment_df %>%
    transmute(
      value = suppressWarnings(as.numeric(.data[[sc]])),
      log_time_from_start_z,
      age_clinical_model,
      sex_for_clinical_model,
      cancer_subgroup_model
    ) %>%
    filter(
      !is.na(value),
      !is.na(log_time_from_start_z),
      !is.na(age_clinical_model),
      !is.na(sex_for_clinical_model),
      !is.na(cancer_subgroup_model)
    ) %>%
    mutate(
      cancer_subgroup_model = droplevels(factor(cancer_subgroup_model)),
      sex_for_clinical_model = droplevels(factor(sex_for_clinical_model))
    )

  if (nrow(df) < 40) return(NULL)

  run_lm_safe(
    df,
    value ~ log_time_from_start_z + age_clinical_model + sex_for_clinical_model + cancer_subgroup_model,
    "log_time_from_start_z",
    sc,
    "time_from_treatment_start",
    extra = list(
      n_with_time = nrow(df),
      median_time_days = median(active_treatment_df$time_from_start_days_model, na.rm = TRUE)
    )
  ) %>%
    mutate(
      direction = case_when(
        beta > 0 ~ "higher_with_longer_time_from_start",
        beta < 0 ~ "lower_with_longer_time_from_start",
        TRUE ~ "no_direction"
      )
    )
})) %>%
  mutate(FDR_global = p.adjust(p_value, method = "BH")) %>%
  arrange(FDR_global, p_value)

# ------------------------------------------------------------
# 15. Save outputs
# ------------------------------------------------------------

write_csv(clinical_variable_sources, file.path(out_dir, "SDY2583_CP16_clinical_variable_sources_STEP6.csv"))
write_csv(all9_matrix_summary, file.path(out_dir, "SDY2583_CP16_ALL9_integrated_matrix_summary_STEP6.csv"))
write_csv(score_availability_by_panel, file.path(out_dir, "SDY2583_CP16_score_availability_by_panel_after_ALL9_STEP6.csv"))
write_csv(cp16_clinical_model_summary, file.path(out_dir, "SDY2583_CP16_clinical_model_summary_STEP6.csv"))
write_csv(cp16_cancer_subgroup_counts, file.path(out_dir, "SDY2583_CP16_cancer_subgroup_counts_STEP6.csv"))
write_csv(cp16_therapy_status_counts, file.path(out_dir, "SDY2583_CP16_therapy_status_counts_STEP6.csv"))
write_csv(cp16_active_treatment_exposure_counts, file.path(out_dir, "SDY2583_CP16_active_treatment_exposure_counts_STEP6.csv"))
write_csv(cp16_therapy_line_counts, file.path(out_dir, "SDY2583_CP16_therapy_line_counts_STEP6.csv"))
write_csv(cp16_time_from_start_summary, file.path(out_dir, "SDY2583_CP16_time_from_start_summary_STEP6.csv"))

write_csv(cancer_subgroup_omnibus_results, file.path(out_dir, "SDY2583_CP16_cancer_subgroup_omnibus_results_STEP6.csv"))
write_csv(cancer_subgroup_one_vs_rest_results, file.path(out_dir, "SDY2583_CP16_cancer_subgroup_one_vs_rest_results_STEP6.csv"))
write_csv(therapy_status_binary_results, file.path(out_dir, "SDY2583_CP16_therapy_status_ongoing_vs_no_ongoing_results_STEP6.csv"))
write_csv(active_modality_results, file.path(out_dir, "SDY2583_CP16_active_treatment_modality_results_STEP6.csv"))
write_csv(therapy_line_binary_results, file.path(out_dir, "SDY2583_CP16_therapy_line_later_vs_first_results_STEP6.csv"))
write_csv(time_from_start_results, file.path(out_dir, "SDY2583_CP16_time_from_treatment_start_results_STEP6.csv"))

save(
  selected_previous_matrix,
  updated_matrix_file,
  integrated_all9,
  clinical_df,
  cp16_score_cols,
  clinical_variable_sources,
  all9_matrix_summary,
  score_availability_by_panel,
  cp16_clinical_model_summary,
  cp16_cancer_subgroup_counts,
  cp16_therapy_status_counts,
  cp16_active_treatment_exposure_counts,
  cp16_therapy_line_counts,
  cp16_time_from_start_summary,
  cancer_subgroup_omnibus_results,
  cancer_subgroup_one_vs_rest_results,
  therapy_status_binary_results,
  active_modality_results,
  therapy_line_binary_results,
  time_from_start_results,
  file = file.path(rdata_out_dir, "SDY2583_CP16_clinical_annotation_STEP6.RData")
)

# ------------------------------------------------------------
# 16. Console output
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SDY2583 CP16 STEP 6 COMPLETE: CLINICAL ANNOTATION\n")
cat("============================================================\n")

cat("\nSelected previous integrated matrix:\n")
print(selected_previous_matrix)

cat("\nUpdated ALL9 matrix saved as:\n")
print(updated_matrix_file)

cat("\nClinical variable sources:\n")
print(as.data.frame(clinical_variable_sources), row.names = FALSE)

cat("\nALL9 integrated matrix summary:\n")
print(as.data.frame(all9_matrix_summary), row.names = FALSE)

cat("\nScore availability by panel after CP16 addition:\n")
print(as.data.frame(score_availability_by_panel), row.names = FALSE)

cat("\nCP16 clinical model summary:\n")
print(as.data.frame(cp16_clinical_model_summary), row.names = FALSE)

cat("\nCP16 cancer subgroup counts:\n")
print(as.data.frame(cp16_cancer_subgroup_counts), row.names = FALSE)

cat("\nCP16 therapy status counts:\n")
print(as.data.frame(cp16_therapy_status_counts), row.names = FALSE)

cat("\nCP16 active-treatment exposure counts:\n")
print(as.data.frame(cp16_active_treatment_exposure_counts), row.names = FALSE)

cat("\nCP16 therapy line counts:\n")
print(as.data.frame(cp16_therapy_line_counts), row.names = FALSE)

cat("\nCP16 time-from-start summary:\n")
print(as.data.frame(cp16_time_from_start_summary), row.names = FALSE)

cat("\nTop CP16 cancer-subgroup omnibus results:\n")
print(as.data.frame(cancer_subgroup_omnibus_results), row.names = FALSE)

cat("\nTop CP16 cancer-subgroup one-vs-rest results:\n")
print(as.data.frame(cancer_subgroup_one_vs_rest_results %>% slice_head(n = 40)), row.names = FALSE)

cat("\nTop CP16 therapy-status ongoing-vs-no-ongoing results:\n")
print(as.data.frame(therapy_status_binary_results), row.names = FALSE)

cat("\nTop CP16 active-treatment modality results:\n")
print(as.data.frame(active_modality_results %>% slice_head(n = 50)), row.names = FALSE)

cat("\nTop CP16 therapy-line later-vs-first results:\n")
print(as.data.frame(therapy_line_binary_results), row.names = FALSE)

cat("\nTop CP16 time-from-start results:\n")
print(as.data.frame(time_from_start_results), row.names = FALSE)

cat("\nFiles saved in:\n")
print(out_dir)

cat("============================================================\n")
