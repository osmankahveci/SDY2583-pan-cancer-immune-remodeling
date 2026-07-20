source(file.path(Sys.getenv("SDY2583_REPO_ROOT", unset = "."), "R", "shared", "bootstrap.R"))

# ============================================================
# SDY2583 CP23
# STEP 6 SAFE: Clinical annotation and integrated matrix update
#
# Panel: CP23 monocyte/macrophage-like myeloid remodeling panel
# Input:
#   - CP23 Step 4 composite score RData
#   - Previous integrated clinical immune matrix, preferably ALL9_with_CP16
# Output:
#   - ALL10 integrated matrix with CP23 scores
#   - CP23 clinical annotation tables
# ============================================================

cran_pkgs <- c("dplyr", "readr", "tibble", "tidyr", "purrr", "stringr", "broom")
for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(broom)
})

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

analysis_dir <- sd_analysis_dir("CP23")
rdata_dir <- file.path(analysis_dir, "11_RData")

integrated_dir <- sd_integrated_dir()
dir.create(integrated_dir, recursive = TRUE, showWarnings = FALSE)

out_dir <- file.path(integrated_dir, "08_CP23_clinical_annotation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

pick_col <- function(df, candidates, required = FALSE, label = "column") {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) > 0) return(hit[1])
  if (required) {
    stop("Required ", label, " not found. Tried: ", paste(candidates, collapse = ", "))
  }
  return(NA_character_)
}

to_num <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

clean_chr <- function(x) {
  y <- as.character(x)
  y[y %in% c("", "NA", "NaN", "NULL", "null", "None")] <- NA_character_
  y
}

get_all_number <- function(path) {
  m <- stringr::str_match(basename(path), "ALL([0-9]+)")
  suppressWarnings(as.integer(m[, 2]))
}

safe_lm_effect <- function(df, score_col, predictor_col, model_type, extra_cols = list()) {
  dat <- df %>%
    transmute(
      y = as.numeric(.data[[score_col]]),
      x = as.numeric(.data[[predictor_col]]),
      age = as.numeric(age_clinical_model),
      sex = factor(sex_clinical_model)
    ) %>%
    filter(is.finite(y), is.finite(x), is.finite(age), !is.na(sex))

  if (nrow(dat) < 20 || length(unique(dat$x)) < 2) {
    return(tibble(
      score = score_col,
      model_type = model_type,
      n_model = nrow(dat),
      beta = NA_real_,
      conf_low = NA_real_,
      conf_high = NA_real_,
      p_value = NA_real_
    ) %>% bind_cols(as_tibble(extra_cols)))
  }

  fit <- tryCatch(lm(y ~ x + age + sex, data = dat), error = function(e) NULL)
  if (is.null(fit)) {
    return(tibble(
      score = score_col,
      model_type = model_type,
      n_model = nrow(dat),
      beta = NA_real_,
      conf_low = NA_real_,
      conf_high = NA_real_,
      p_value = NA_real_
    ) %>% bind_cols(as_tibble(extra_cols)))
  }

  td <- broom::tidy(fit, conf.int = TRUE)
  row <- td %>% filter(term == "x")

  tibble(
    score = score_col,
    model_type = model_type,
    n_model = nrow(dat),
    beta = row$estimate[1],
    conf_low = row$conf.low[1],
    conf_high = row$conf.high[1],
    p_value = row$p.value[1]
  ) %>% bind_cols(as_tibble(extra_cols))
}

safe_lm_time <- function(df, score_col, model_type, extra_cols = list()) {
  dat <- df %>%
    transmute(
      y = as.numeric(.data[[score_col]]),
      time_days = as.numeric(time_from_start_days_model),
      age = as.numeric(age_clinical_model),
      sex = factor(sex_clinical_model)
    ) %>%
    filter(is.finite(y), is.finite(time_days), time_days >= 0, is.finite(age), !is.na(sex)) %>%
    mutate(log_time_z = as.numeric(scale(log1p(time_days))))

  if (nrow(dat) < 20 || length(unique(dat$log_time_z)) < 2) {
    return(tibble(
      score = score_col,
      model_type = model_type,
      n_model = nrow(dat),
      median_time_days = ifelse(nrow(dat) > 0, median(dat$time_days, na.rm = TRUE), NA_real_),
      beta = NA_real_,
      conf_low = NA_real_,
      conf_high = NA_real_,
      p_value = NA_real_
    ) %>% bind_cols(as_tibble(extra_cols)))
  }

  fit <- tryCatch(lm(y ~ log_time_z + age + sex, data = dat), error = function(e) NULL)
  if (is.null(fit)) {
    return(tibble(
      score = score_col,
      model_type = model_type,
      n_model = nrow(dat),
      median_time_days = median(dat$time_days, na.rm = TRUE),
      beta = NA_real_,
      conf_low = NA_real_,
      conf_high = NA_real_,
      p_value = NA_real_
    ) %>% bind_cols(as_tibble(extra_cols)))
  }

  td <- broom::tidy(fit, conf.int = TRUE)
  row <- td %>% filter(term == "log_time_z")

  tibble(
    score = score_col,
    model_type = model_type,
    n_model = nrow(dat),
    median_time_days = median(dat$time_days, na.rm = TRUE),
    beta = row$estimate[1],
    conf_low = row$conf.low[1],
    conf_high = row$conf.high[1],
    p_value = row$p.value[1]
  ) %>% bind_cols(as_tibble(extra_cols))
}

# ------------------------------------------------------------
# Load CP23 Step 4 composite score objects
# ------------------------------------------------------------

step4_candidates <- list.files(
  rdata_dir,
  pattern = "CP23.*STEP4.*composite.*\\.RData$|CP23.*STEP4.*\\.RData$",
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(step4_candidates) == 0) {
  step4_candidates <- list.files(rdata_dir, pattern = "CP23.*\\.RData$", full.names = TRUE, ignore.case = TRUE)
}

if (length(step4_candidates) == 0) {
  stop("No CP23 RData file found in: ", rdata_dir)
}

# Prefer Step4 composite RData if multiple files exist.
step4_file <- step4_candidates[grepl("STEP4|composite", basename(step4_candidates), ignore.case = TRUE)][1]
if (is.na(step4_file)) step4_file <- step4_candidates[1]

load(step4_file)

env_names <- ls()
df_names <- env_names[vapply(env_names, function(obj) {
  x <- get(obj)
  is.data.frame(x) && ("subject_id" %in% names(x))
}, logical(1))]

if (length(df_names) == 0) {
  stop("No subject-level data frame with subject_id found after loading Step 4 RData.")
}

score_df_candidates <- df_names[vapply(df_names, function(obj) {
  x <- get(obj)
  length(grep("^CP23_.*_score$", names(x), value = TRUE)) >= 5
}, logical(1))]

if (length(score_df_candidates) == 0) {
  stop("No data frame with at least 5 CP23 composite score columns found.")
}

score_counts <- vapply(score_df_candidates, function(obj) {
  x <- get(obj)
  length(grep("^CP23_.*_score$", names(x), value = TRUE))
}, integer(1))

score_df_name <- score_df_candidates[which.max(score_counts)]
cp23_score_data <- get(score_df_name)

cp23_score_cols <- grep("^CP23_.*_score$", names(cp23_score_data), value = TRUE)

cp23_scores_to_add <- cp23_score_data %>%
  select(subject_id, all_of(cp23_score_cols)) %>%
  distinct(subject_id, .keep_all = TRUE)

# Score labels
if (exists("composite_score_dictionary")) {
  cp23_score_dictionary <- composite_score_dictionary %>%
    select(any_of(c("composite_score", "score_label"))) %>%
    distinct()
} else {
  cp23_score_dictionary <- tibble(
    composite_score = cp23_score_cols,
    score_label = cp23_score_cols
  )
}

if (!("score_label" %in% names(cp23_score_dictionary))) {
  cp23_score_dictionary$score_label <- cp23_score_dictionary$composite_score
}

score_label_lookup <- cp23_score_dictionary %>%
  filter(composite_score %in% cp23_score_cols) %>%
  distinct(composite_score, .keep_all = TRUE)

# ------------------------------------------------------------
# Select previous integrated matrix
# ------------------------------------------------------------

preferred_previous <- file.path(
  integrated_dir,
  "SDY2583_integrated_clinical_immune_score_matrix_ALL9_with_CP16.csv"
)

if (file.exists(preferred_previous)) {
  previous_matrix_file <- preferred_previous
} else {
  integrated_files <- list.files(
    integrated_dir,
    pattern = "SDY2583_integrated_clinical_immune_score_matrix_ALL[0-9]+_with_.*\\.csv$",
    full.names = TRUE
  )

  if (length(integrated_files) == 0) {
    stop("No previous integrated matrix found in: ", integrated_dir)
  }

  # Avoid using an already CP23-updated output as the previous matrix.
  integrated_files_no_cp23 <- integrated_files[!grepl("CP23", basename(integrated_files), ignore.case = TRUE)]
  if (length(integrated_files_no_cp23) > 0) integrated_files <- integrated_files_no_cp23

  all_nums <- get_all_number(integrated_files)
  previous_matrix_file <- integrated_files[which.max(all_nums)]
}

previous_matrix <- read_csv(previous_matrix_file, show_col_types = FALSE)

previous_all <- get_all_number(previous_matrix_file)
new_all <- ifelse(is.finite(previous_all), previous_all + 1L, 10L)

updated_matrix_file <- file.path(
  integrated_dir,
  paste0("SDY2583_integrated_clinical_immune_score_matrix_ALL", new_all, "_with_CP23.csv")
)

# Remove old CP23 score columns if the script is rerun.
previous_matrix_no_cp23 <- previous_matrix %>%
  select(-any_of(grep("^CP23_.*_score$", names(previous_matrix), value = TRUE)))

updated_matrix <- previous_matrix_no_cp23 %>%
  left_join(cp23_scores_to_add, by = "subject_id")

write_csv(updated_matrix, updated_matrix_file)

# ------------------------------------------------------------
# Clinical variable standardization
# ------------------------------------------------------------

disease_col <- pick_col(updated_matrix, c("disease_group_model", "disease_group", "disease", "group"), TRUE, "disease group")
age_col <- pick_col(updated_matrix, c("age_clinical", "age_for_model", "age", "age_raw"), TRUE, "age")
sex_col <- pick_col(updated_matrix, c("sex_for_clinical_model", "sex_for_model", "sex"), TRUE, "sex")
subgroup_col <- pick_col(updated_matrix, c("cancer_subgroup_model", "cancer_subgroup", "cancer_type", "tumor_type"), FALSE, "cancer subgroup")
therapy_status_col <- pick_col(updated_matrix, c("therapy_status_model", "therapy_status_4level", "therapy_status"), FALSE, "therapy status")
therapy_line_col <- pick_col(updated_matrix, c("therapy_line_number_model", "therapy_line_number", "line_number"), FALSE, "therapy line")
time_col <- pick_col(updated_matrix, c("time_from_start_days_model", "time_from_start_days", "time_from_treatment_start_days"), FALSE, "time from start")

exposure_candidates <- c(
  "chemotherapy", "targeted_therapy", "any_immunotherapy", "ici_immunotherapy",
  "endocrine_hormonal", "adc", "experimental", "radiotherapy"
)
exposure_cols <- exposure_candidates[exposure_candidates %in% names(updated_matrix)]

clinical_df <- updated_matrix %>%
  mutate(
    disease_group_for_clinical = clean_chr(.data[[disease_col]]),
    age_clinical_model = to_num(.data[[age_col]]),
    sex_clinical_model = clean_chr(.data[[sex_col]]),
    cancer_subgroup_model = if (!is.na(subgroup_col)) clean_chr(.data[[subgroup_col]]) else NA_character_,
    therapy_status_model = if (!is.na(therapy_status_col)) clean_chr(.data[[therapy_status_col]]) else NA_character_,
    therapy_line_number_model = if (!is.na(therapy_line_col)) to_num(.data[[therapy_line_col]]) else NA_real_,
    time_from_start_days_model = if (!is.na(time_col)) to_num(.data[[time_col]]) else NA_real_
  ) %>%
  mutate(
    disease_group_for_clinical = case_when(
      str_detect(tolower(disease_group_for_clinical), "cancer|patient") ~ "Cancer patient",
      str_detect(tolower(disease_group_for_clinical), "healthy|control") ~ "Healthy control",
      TRUE ~ disease_group_for_clinical
    ),
    sex_clinical_model = ifelse(is.na(sex_clinical_model), "Not Specified", sex_clinical_model),
    sex_clinical_model = factor(sex_clinical_model)
  )

cancer_df <- clinical_df %>%
  filter(disease_group_for_clinical == "Cancer patient") %>%
  filter(is.finite(age_clinical_model), !is.na(sex_clinical_model))

# ------------------------------------------------------------
# Availability summaries
# ------------------------------------------------------------

clinical_variable_sources <- tibble(
  variable = c(
    "disease_group", "age", "sex", "cancer_subgroup",
    "therapy_status", "therapy_line_number", "time_from_start_days"
  ),
  source_column = c(
    disease_col, age_col, sex_col, subgroup_col,
    therapy_status_col, therapy_line_col, time_col
  )
)

integrated_summary <- tibble(
  n_subjects = nrow(updated_matrix),
  n_unique_subjects = n_distinct(updated_matrix$subject_id),
  n_with_any_CP23_score = sum(rowSums(!is.na(updated_matrix[, cp23_score_cols, drop = FALSE])) > 0),
  n_with_all_CP23_scores = sum(rowSums(!is.na(updated_matrix[, cp23_score_cols, drop = FALSE])) == length(cp23_score_cols)),
  n_total_score_cols_ALL = length(grep("_score$", names(updated_matrix), value = TRUE)),
  median_available_score_cols = median(rowSums(!is.na(updated_matrix[, grep("_score$", names(updated_matrix), value = TRUE), drop = FALSE])))
)

score_availability_by_panel <- tibble(
  panel = stringr::str_extract(grep("_score$", names(updated_matrix), value = TRUE), "^CP[0-9]+")
) %>%
  filter(!is.na(panel)) %>%
  count(panel, name = "n_score_cols") %>%
  arrange(panel) %>%
  rowwise() %>%
  mutate(
    n_subjects_with_any_score = {
      panel_cols <- grep(paste0("^", panel, "_.*_score$"), names(updated_matrix), value = TRUE)
      sum(rowSums(!is.na(updated_matrix[, panel_cols, drop = FALSE])) > 0)
    }
  ) %>%
  ungroup()

cp23_clinical_model_summary <- tibble(
  n_cancer_total = sum(clinical_df$disease_group_for_clinical == "Cancer patient", na.rm = TRUE),
  n_with_valid_age = sum(clinical_df$disease_group_for_clinical == "Cancer patient" & is.finite(clinical_df$age_clinical_model), na.rm = TRUE),
  n_with_valid_sex = sum(clinical_df$disease_group_for_clinical == "Cancer patient" & !is.na(clinical_df$sex_clinical_model), na.rm = TRUE),
  n_with_cancer_subgroup = sum(clinical_df$disease_group_for_clinical == "Cancer patient" & !is.na(clinical_df$cancer_subgroup_model), na.rm = TRUE),
  n_with_therapy_status = sum(clinical_df$disease_group_for_clinical == "Cancer patient" & !is.na(clinical_df$therapy_status_model), na.rm = TRUE),
  n_with_all_core_covariates = nrow(cancer_df),
  n_ongoing_active_treatment = sum(cancer_df$therapy_status_model == "ongoing_active_treatment", na.rm = TRUE),
  n_no_ongoing_therapy = sum(cancer_df$therapy_status_model == "no_ongoing_therapy", na.rm = TRUE),
  n_no_treatment_data = sum(cancer_df$therapy_status_model == "no_treatment_data", na.rm = TRUE)
)

cancer_subgroup_counts <- cancer_df %>%
  filter(!is.na(cancer_subgroup_model)) %>%
  count(cancer_subgroup_model, sort = TRUE)

therapy_status_counts <- cancer_df %>%
  filter(!is.na(therapy_status_model)) %>%
  count(therapy_status_model, sort = TRUE)

# Exposure summary in active-treatment patients
active_df <- cancer_df %>%
  filter(therapy_status_model == "ongoing_active_treatment")

active_treatment_exposure_counts <- bind_rows(lapply(exposure_cols, function(exposure) {
  x <- active_df[[exposure]]
  exposed <- sum(tolower(as.character(x)) %in% c("true", "1", "yes", "y", "exposed"), na.rm = TRUE)
  unexposed <- sum(!(tolower(as.character(x)) %in% c("true", "1", "yes", "y", "exposed")) & !is.na(x), na.rm = TRUE)

  tier <- case_when(
    exposed >= 20 & unexposed >= 20 ~ "primary_or_exploratory_powered",
    exposed >= 10 & unexposed >= 20 ~ "exploratory_underpowered",
    TRUE ~ "descriptive_only"
  )

  tibble(
    exposure = exposure,
    n_exposed = exposed,
    n_unexposed = unexposed,
    analysis_tier = tier
  )
}))

therapy_line_counts <- cancer_df %>%
  mutate(
    therapy_line_group = case_when(
      therapy_line_number_model == 1 ~ "first_line",
      therapy_line_number_model >= 2 ~ "later_line",
      TRUE ~ NA_character_
    )
  ) %>%
  count(therapy_line_number_model, therapy_line_group, sort = FALSE)

time_from_start_summary <- cancer_df %>%
  filter(is.finite(time_from_start_days_model)) %>%
  summarise(
    n_with_time = n(),
    min_days = min(time_from_start_days_model, na.rm = TRUE),
    q1_days = quantile(time_from_start_days_model, 0.25, na.rm = TRUE),
    median_days = median(time_from_start_days_model, na.rm = TRUE),
    mean_days = mean(time_from_start_days_model, na.rm = TRUE),
    q3_days = quantile(time_from_start_days_model, 0.75, na.rm = TRUE),
    max_days = max(time_from_start_days_model, na.rm = TRUE)
  )

# ------------------------------------------------------------
# Clinical models
# ------------------------------------------------------------

score_cols <- cp23_score_cols

# Cancer subgroup omnibus
cancer_subgroup_omnibus_results <- bind_rows(lapply(score_cols, function(sc) {
  dat <- cancer_df %>%
    transmute(
      y = as.numeric(.data[[sc]]),
      cancer_subgroup_model = factor(cancer_subgroup_model),
      age = age_clinical_model,
      sex = factor(sex_clinical_model)
    ) %>%
    filter(is.finite(y), !is.na(cancer_subgroup_model), is.finite(age), !is.na(sex)) %>%
    droplevels()

  if (nrow(dat) < 30 || nlevels(dat$cancer_subgroup_model) < 2) {
    return(tibble(score = sc, n_model = nrow(dat), n_subgroups = nlevels(dat$cancer_subgroup_model),
                  omnibus_F = NA_real_, omnibus_p = NA_real_))
  }

  fit <- tryCatch(lm(y ~ cancer_subgroup_model + age + sex, data = dat), error = function(e) NULL)
  if (is.null(fit)) {
    return(tibble(score = sc, n_model = nrow(dat), n_subgroups = nlevels(dat$cancer_subgroup_model),
                  omnibus_F = NA_real_, omnibus_p = NA_real_))
  }

  av <- anova(fit)
  tibble(
    score = sc,
    n_model = nrow(dat),
    n_subgroups = nlevels(dat$cancer_subgroup_model),
    omnibus_F = av["cancer_subgroup_model", "F value"],
    omnibus_p = av["cancer_subgroup_model", "Pr(>F)"]
  )
})) %>%
  left_join(score_label_lookup, by = c("score" = "composite_score")) %>%
  mutate(FDR_global = p.adjust(omnibus_p, method = "BH")) %>%
  arrange(FDR_global, omnibus_p)

# Cancer subgroup one-vs-rest
subgroups <- cancer_df %>%
  filter(!is.na(cancer_subgroup_model)) %>%
  count(cancer_subgroup_model, sort = TRUE) %>%
  pull(cancer_subgroup_model)

cancer_subgroup_one_vs_rest_results <- bind_rows(lapply(score_cols, function(sc) {
  bind_rows(lapply(subgroups, function(sg) {
    dat <- cancer_df %>%
      filter(!is.na(cancer_subgroup_model)) %>%
      mutate(in_subgroup = as.integer(cancer_subgroup_model == sg))

    n_subgroup <- sum(dat$in_subgroup == 1, na.rm = TRUE)
    n_rest <- sum(dat$in_subgroup == 0, na.rm = TRUE)

    safe_lm_effect(
      dat,
      score_col = sc,
      predictor_col = "in_subgroup",
      model_type = "cancer_subgroup_one_vs_rest",
      extra_cols = list(cancer_subgroup = sg, n_subgroup = n_subgroup, n_rest = n_rest)
    )
  }))
})) %>%
  left_join(score_label_lookup, by = c("score" = "composite_score")) %>%
  group_by(score) %>%
  mutate(FDR_within_score = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  mutate(
    FDR_global = p.adjust(p_value, method = "BH"),
    direction = ifelse(beta >= 0, "higher_in_subgroup_vs_rest", "lower_in_subgroup_vs_rest")
  ) %>%
  arrange(FDR_global, p_value)

# Therapy status: ongoing active vs no ongoing therapy
therapy_status_model_df <- cancer_df %>%
  filter(therapy_status_model %in% c("ongoing_active_treatment", "no_ongoing_therapy")) %>%
  mutate(ongoing_active = as.integer(therapy_status_model == "ongoing_active_treatment"))

therapy_status_ongoing_vs_no_results <- bind_rows(lapply(score_cols, function(sc) {
  safe_lm_effect(
    therapy_status_model_df,
    score_col = sc,
    predictor_col = "ongoing_active",
    model_type = "therapy_status_ongoing_vs_no_ongoing",
    extra_cols = list(
      n_ongoing = sum(therapy_status_model_df$ongoing_active == 1, na.rm = TRUE),
      n_no_ongoing = sum(therapy_status_model_df$ongoing_active == 0, na.rm = TRUE)
    )
  )
})) %>%
  left_join(score_label_lookup, by = c("score" = "composite_score")) %>%
  mutate(
    FDR_global = p.adjust(p_value, method = "BH"),
    direction = ifelse(beta >= 0, "higher_in_ongoing_active_treatment", "lower_in_ongoing_active_treatment")
  ) %>%
  arrange(FDR_global, p_value)

# Active-treatment modality
active_treatment_modality_results <- bind_rows(lapply(score_cols, function(sc) {
  bind_rows(lapply(exposure_cols, function(exposure) {
    dat <- active_df
    x <- dat[[exposure]]
    dat$exposure_bin <- as.integer(tolower(as.character(x)) %in% c("true", "1", "yes", "y", "exposed"))

    n_exposed <- sum(dat$exposure_bin == 1, na.rm = TRUE)
    n_unexposed <- sum(dat$exposure_bin == 0, na.rm = TRUE)

    tier <- case_when(
      n_exposed >= 20 & n_unexposed >= 20 ~ "primary_or_exploratory_powered",
      n_exposed >= 10 & n_unexposed >= 20 ~ "exploratory_underpowered",
      TRUE ~ "descriptive_only"
    )

    safe_lm_effect(
      dat,
      score_col = sc,
      predictor_col = "exposure_bin",
      model_type = "active_treatment_modality",
      extra_cols = list(
        exposure = exposure,
        n_exposed = n_exposed,
        n_unexposed = n_unexposed,
        analysis_tier = tier
      )
    )
  }))
})) %>%
  left_join(score_label_lookup, by = c("score" = "composite_score")) %>%
  group_by(exposure) %>%
  mutate(FDR_within_exposure = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  group_by(score) %>%
  mutate(FDR_within_score = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  mutate(
    FDR_global = p.adjust(p_value, method = "BH"),
    direction = ifelse(beta >= 0, "higher_in_exposed", "lower_in_exposed")
  ) %>%
  arrange(FDR_global, p_value)

# Therapy line: later vs first
therapy_line_model_df <- cancer_df %>%
  filter(is.finite(therapy_line_number_model), therapy_line_number_model >= 1) %>%
  mutate(
    line_group = case_when(
      therapy_line_number_model == 1 ~ "first_line",
      therapy_line_number_model >= 2 ~ "later_line",
      TRUE ~ NA_character_
    ),
    later_line = as.integer(line_group == "later_line")
  ) %>%
  filter(!is.na(line_group))

therapy_line_later_vs_first_results <- bind_rows(lapply(score_cols, function(sc) {
  safe_lm_effect(
    therapy_line_model_df,
    score_col = sc,
    predictor_col = "later_line",
    model_type = "therapy_line_later_vs_first",
    extra_cols = list(
      n_first_line = sum(therapy_line_model_df$later_line == 0, na.rm = TRUE),
      n_later_line = sum(therapy_line_model_df$later_line == 1, na.rm = TRUE)
    )
  )
})) %>%
  left_join(score_label_lookup, by = c("score" = "composite_score")) %>%
  mutate(
    FDR_global = p.adjust(p_value, method = "BH"),
    direction = ifelse(beta >= 0, "higher_in_later_line", "lower_in_later_line")
  ) %>%
  arrange(FDR_global, p_value)

# Time from treatment start
time_from_start_results <- bind_rows(lapply(score_cols, function(sc) {
  safe_lm_time(
    cancer_df,
    score_col = sc,
    model_type = "time_from_treatment_start",
    extra_cols = list()
  )
})) %>%
  left_join(score_label_lookup, by = c("score" = "composite_score")) %>%
  mutate(
    FDR_global = p.adjust(p_value, method = "BH"),
    direction = ifelse(beta >= 0, "higher_with_longer_time_from_start", "lower_with_longer_time_from_start")
  ) %>%
  arrange(FDR_global, p_value)

# ------------------------------------------------------------
# Save outputs
# ------------------------------------------------------------

write_csv(clinical_variable_sources, file.path(out_dir, "SDY2583_CP23_clinical_variable_sources_STEP6.csv"))
write_csv(integrated_summary, file.path(out_dir, "SDY2583_CP23_ALL10_integrated_matrix_summary_STEP6.csv"))
write_csv(score_availability_by_panel, file.path(out_dir, "SDY2583_score_availability_by_panel_after_CP23_STEP6.csv"))
write_csv(cp23_clinical_model_summary, file.path(out_dir, "SDY2583_CP23_clinical_model_summary_STEP6.csv"))
write_csv(cancer_subgroup_counts, file.path(out_dir, "SDY2583_CP23_cancer_subgroup_counts_STEP6.csv"))
write_csv(therapy_status_counts, file.path(out_dir, "SDY2583_CP23_therapy_status_counts_STEP6.csv"))
write_csv(active_treatment_exposure_counts, file.path(out_dir, "SDY2583_CP23_active_treatment_exposure_counts_STEP6.csv"))
write_csv(therapy_line_counts, file.path(out_dir, "SDY2583_CP23_therapy_line_counts_STEP6.csv"))
write_csv(time_from_start_summary, file.path(out_dir, "SDY2583_CP23_time_from_start_summary_STEP6.csv"))

write_csv(cancer_subgroup_omnibus_results, file.path(out_dir, "SDY2583_CP23_cancer_subgroup_omnibus_results_STEP6.csv"))
write_csv(cancer_subgroup_one_vs_rest_results, file.path(out_dir, "SDY2583_CP23_cancer_subgroup_one_vs_rest_results_STEP6.csv"))
write_csv(therapy_status_ongoing_vs_no_results, file.path(out_dir, "SDY2583_CP23_therapy_status_ongoing_vs_no_results_STEP6.csv"))
write_csv(active_treatment_modality_results, file.path(out_dir, "SDY2583_CP23_active_treatment_modality_results_STEP6.csv"))
write_csv(therapy_line_later_vs_first_results, file.path(out_dir, "SDY2583_CP23_therapy_line_later_vs_first_results_STEP6.csv"))
write_csv(time_from_start_results, file.path(out_dir, "SDY2583_CP23_time_from_start_results_STEP6.csv"))

save(
  previous_matrix_file,
  updated_matrix_file,
  clinical_variable_sources,
  integrated_summary,
  score_availability_by_panel,
  cp23_clinical_model_summary,
  cancer_subgroup_counts,
  therapy_status_counts,
  active_treatment_exposure_counts,
  therapy_line_counts,
  time_from_start_summary,
  cancer_subgroup_omnibus_results,
  cancer_subgroup_one_vs_rest_results,
  therapy_status_ongoing_vs_no_results,
  active_treatment_modality_results,
  therapy_line_later_vs_first_results,
  time_from_start_results,
  cp23_score_cols,
  cp23_score_dictionary,
  file = file.path(out_dir, "SDY2583_CP23_STEP6_clinical_annotation.RData")
)

# ------------------------------------------------------------
# Console report
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SDY2583 CP23 STEP 6 COMPLETE: CLINICAL ANNOTATION\n")
cat("============================================================\n")

cat("\nSelected previous integrated matrix:\n")
print(previous_matrix_file)

cat("\nUpdated ALL", new_all, " matrix saved as:\n", sep = "")
print(updated_matrix_file)

cat("\nClinical variable sources:\n")
print(as.data.frame(clinical_variable_sources), row.names = FALSE)

cat("\nALL", new_all, " integrated matrix summary:\n", sep = "")
print(as.data.frame(integrated_summary), row.names = FALSE)

cat("\nScore availability by panel after CP23 addition:\n")
print(as.data.frame(score_availability_by_panel), row.names = FALSE)

cat("\nCP23 clinical model summary:\n")
print(as.data.frame(cp23_clinical_model_summary), row.names = FALSE)

cat("\nCP23 cancer subgroup counts:\n")
print(as.data.frame(cancer_subgroup_counts), row.names = FALSE)

cat("\nCP23 therapy status counts:\n")
print(as.data.frame(therapy_status_counts), row.names = FALSE)

cat("\nCP23 active-treatment exposure counts:\n")
print(as.data.frame(active_treatment_exposure_counts), row.names = FALSE)

cat("\nCP23 therapy line counts:\n")
print(as.data.frame(therapy_line_counts), row.names = FALSE)

cat("\nCP23 time-from-start summary:\n")
print(as.data.frame(time_from_start_summary), row.names = FALSE)

cat("\nTop CP23 cancer-subgroup omnibus results:\n")
print(as.data.frame(cancer_subgroup_omnibus_results %>% arrange(FDR_global, omnibus_p)), row.names = FALSE)

cat("\nTop CP23 cancer-subgroup one-vs-rest results:\n")
print(as.data.frame(cancer_subgroup_one_vs_rest_results %>% arrange(FDR_global, p_value) %>% head(30)), row.names = FALSE)

cat("\nTop CP23 therapy-status ongoing-vs-no-ongoing results:\n")
print(as.data.frame(therapy_status_ongoing_vs_no_results %>% arrange(FDR_global, p_value)), row.names = FALSE)

cat("\nTop CP23 active-treatment modality results:\n")
print(as.data.frame(active_treatment_modality_results %>% arrange(FDR_global, p_value) %>% head(50)), row.names = FALSE)

cat("\nTop CP23 therapy-line later-vs-first results:\n")
print(as.data.frame(therapy_line_later_vs_first_results %>% arrange(FDR_global, p_value)), row.names = FALSE)

cat("\nTop CP23 time-from-start results:\n")
print(as.data.frame(time_from_start_results %>% arrange(FDR_global, p_value)), row.names = FALSE)

cat("\nFiles saved in:\n")
print(out_dir)
cat("============================================================\n")
