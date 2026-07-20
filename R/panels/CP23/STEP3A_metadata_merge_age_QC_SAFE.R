# ============================================================
# SDY2583 CP23
# STEP 3A SAFE: Metadata merge / age QC
#
# Purpose:
#   Merge CP23 feature table with the integrated clinical immune
#   matrix, standardize disease/age/sex fields, perform age QC,
#   and summarize event-count/technical QC by disease group.
#
# Input:
#   outputs/CP23/11_RData/
#     SDY2583_CP23_STEP2_feature_extraction_myeloid_macrophage_like.RData
#
# Preferred metadata source:
#   data/derived/clinical_integration/
#     SDY2583_integrated_clinical_immune_score_matrix_ALL9_with_CP16.csv
#
# Output:
#   outputs/CP23/03_metadata_merge_age_QC
#   outputs/CP23/11_RData
# ============================================================

rm(list = ls())
source(file.path(Sys.getenv("SDY2583_REPO_ROOT", unset = "."), "R", "shared", "bootstrap.R"))

# ------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------

cran_pkgs <- c("dplyr", "readr", "stringr", "tibble", "tidyr", "purrr")

for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
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

analysis_dir <- sd_analysis_dir("CP23")
rdata_dir <- file.path(analysis_dir, "11_RData")
out_dir <- file.path(analysis_dir, "03_metadata_merge_age_QC")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rdata_dir, recursive = TRUE, showWarnings = FALSE)

step2_rdata <- file.path(rdata_dir, "SDY2583_CP23_STEP2_feature_extraction_myeloid_macrophage_like.RData")

if (!file.exists(step2_rdata)) {
  stop("Step 2 RData bulunamadı: ", step2_rdata)
}

load(step2_rdata)

# Reset paths after load.
analysis_dir <- sd_analysis_dir("CP23")
rdata_dir <- file.path(analysis_dir, "11_RData")
out_dir <- file.path(analysis_dir, "03_metadata_merge_age_QC")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!exists("cp23_features")) stop("cp23_features bulunamadı.")
if (!exists("feature_dictionary_full")) stop("feature_dictionary_full bulunamadı.")

# ------------------------------------------------------------
# 3. Metadata source discovery
# ------------------------------------------------------------

matrix_dir <- sd_integrated_dir()

candidate_metadata_files <- c(
  file.path(matrix_dir, "SDY2583_integrated_clinical_immune_score_matrix_ALL9_with_CP16.csv"),
  file.path(matrix_dir, "SDY2583_integrated_clinical_immune_score_matrix_ALL8_with_CP10.csv"),
  file.path(matrix_dir, "SDY2583_integrated_clinical_immune_score_matrix_ALL7_with_CP22.csv")
)

candidate_metadata_files <- candidate_metadata_files[file.exists(candidate_metadata_files)]

if (length(candidate_metadata_files) == 0 && dir.exists(matrix_dir)) {
  all_csv <- list.files(matrix_dir, pattern = "\\.csv$", full.names = TRUE)
  all_csv <- all_csv[
    grepl("integrated_clinical_immune_score_matrix|integrated.*matrix|ALL", basename(all_csv), ignore.case = TRUE) &
      !grepl("summary|counts|annotation|subgroup|therapy|manifest|dictionary|results", basename(all_csv), ignore.case = TRUE)
  ]

  if (length(all_csv) > 0) {
    ranked <- tibble(file = all_csv) %>%
      mutate(
        base = basename(file),
        all_number = suppressWarnings(as.numeric(str_match(base, "ALL([0-9]+)")[, 2])),
        all_number = ifelse(is.na(all_number), -Inf, all_number),
        mtime = file.info(file)$mtime
      ) %>%
      arrange(desc(all_number), desc(mtime))

    candidate_metadata_files <- ranked$file
  }
}

if (length(candidate_metadata_files) == 0) {
  stop("Integrated clinical immune metadata matrix bulunamadı.")
}

metadata_source_file <- candidate_metadata_files[1]
metadata_source_label <- "integrated_clinical_immune_matrix_csv"

metadata_raw <- readr::read_csv(metadata_source_file, show_col_types = FALSE)

if (!("subject_id" %in% names(metadata_raw))) {
  stop("Metadata içinde subject_id yok: ", metadata_source_file)
}

# ------------------------------------------------------------
# 4. Standardization helpers
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
    y %in% c("Not Specified", "Not specified", "not specified", "Unknown", "unknown") ~ "Not Specified",
    is.na(y) | y == "" | y == "NA" ~ NA_character_,
    TRUE ~ y
  )
}

# ------------------------------------------------------------
# 5. Standardize metadata fields
# ------------------------------------------------------------

disease_col <- first_existing_col(metadata_raw, c("disease_group_model", "disease_group_clinical", "disease_group", "group", "condition"))
age_col <- first_existing_col(metadata_raw, c("age_for_model", "age_clinical", "age_years", "age", "age_raw"))
sex_col <- first_existing_col(metadata_raw, c("sex_for_clinical_model", "sex_clinical", "sex", "gender"))
cancer_subgroup_col <- first_existing_col(metadata_raw, c("cancer_subgroup_model", "cancer_subgroup", "cancer_type_model", "cancer_type", "tumor_type"))
therapy_status_col <- first_existing_col(metadata_raw, c("therapy_status_model", "therapy_status_4level", "therapy_status", "treatment_status"))
therapy_line_col <- first_existing_col(metadata_raw, c("therapy_line_number_model", "therapy_line_number", "treatment_line_number", "line_of_therapy"))
time_col <- first_existing_col(metadata_raw, c("time_from_start_days_model", "time_from_start_days", "time_from_treatment_start_days", "days_from_treatment_start"))

if (is.na(disease_col)) stop("Disease column bulunamadı.")
if (is.na(age_col)) stop("Age column bulunamadı.")
if (is.na(sex_col)) stop("Sex column bulunamadı.")

metadata_standardized <- metadata_raw %>%
  mutate(
    subject_id = as.character(subject_id),
    disease_group = standardize_disease(.data[[disease_col]]),
    disease_group = factor(disease_group, levels = c("Healthy control", "Cancer patient")),
    age_raw = suppressWarnings(as.numeric(.data[[age_col]])),
    age_for_model = ifelse(age_raw >= 18 & age_raw <= 100, age_raw, NA_real_),
    invalid_or_missing_age_for_model = is.na(age_for_model),
    sex = standardize_sex(.data[[sex_col]]),
    sex = factor(sex),
    sex_binary = ifelse(as.character(sex) %in% c("Female", "Male"), as.character(sex), NA_character_),
    sex_binary = factor(sex_binary, levels = c("Female", "Male"))
  )

if (!is.na(cancer_subgroup_col)) {
  metadata_standardized <- metadata_standardized %>%
    mutate(cancer_subgroup = as.character(.data[[cancer_subgroup_col]]))
} else {
  metadata_standardized$cancer_subgroup <- NA_character_
}

if (!is.na(therapy_status_col)) {
  metadata_standardized <- metadata_standardized %>%
    mutate(therapy_status_4level = as.character(.data[[therapy_status_col]]))
} else {
  metadata_standardized$therapy_status_4level <- NA_character_
}

if (!is.na(therapy_line_col)) {
  metadata_standardized <- metadata_standardized %>%
    mutate(therapy_line_number = suppressWarnings(as.numeric(.data[[therapy_line_col]])))
} else {
  metadata_standardized$therapy_line_number <- NA_real_
}

if (!is.na(time_col)) {
  metadata_standardized <- metadata_standardized %>%
    mutate(time_from_start_days = suppressWarnings(as.numeric(.data[[time_col]])))
} else {
  metadata_standardized$time_from_start_days <- NA_real_
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
  if (!(ee %in% names(metadata_standardized))) {
    metadata_standardized[[ee]] <- NA
  }
}

metadata_keep <- metadata_standardized %>%
  select(
    subject_id,
    disease_group,
    age_raw,
    age_for_model,
    invalid_or_missing_age_for_model,
    sex,
    sex_binary,
    cancer_subgroup,
    therapy_status_4level,
    all_of(exposure_vars),
    therapy_line_number,
    time_from_start_days
  ) %>%
  distinct(subject_id, .keep_all = TRUE)

# ------------------------------------------------------------
# 6. Merge with CP23 features
# ------------------------------------------------------------

cp23_analysis_data <- cp23_features %>%
  mutate(subject_id = as.character(subject_id)) %>%
  left_join(metadata_keep, by = "subject_id") %>%
  mutate(
    model_ready_age_sex =
      feature_ok == TRUE &
      !is.na(disease_group) &
      !is.na(age_for_model) &
      !is.na(sex),
    model_ready_binary_sex =
      feature_ok == TRUE &
      !is.na(disease_group) &
      !is.na(age_for_model) &
      !is.na(sex_binary)
  )

# ------------------------------------------------------------
# 7. Summaries
# ------------------------------------------------------------

metadata_source_summary <- tibble(
  metadata_source_label = metadata_source_label,
  metadata_source_file = metadata_source_file,
  n_metadata_rows = nrow(metadata_raw),
  n_metadata_subjects = n_distinct(metadata_raw$subject_id),
  n_standardized_metadata_subjects = n_distinct(metadata_keep$subject_id)
)

clinical_variable_availability <- tibble(
  variable = c(
    "cancer_subgroup",
    "therapy_status_4level",
    exposure_vars,
    "therapy_line_number",
    "time_from_start_days"
  ),
  present = c(
    "cancer_subgroup" %in% names(cp23_analysis_data),
    "therapy_status_4level" %in% names(cp23_analysis_data),
    exposure_vars %in% names(cp23_analysis_data),
    "therapy_line_number" %in% names(cp23_analysis_data),
    "time_from_start_days" %in% names(cp23_analysis_data)
  ),
  n_nonmissing = c(
    sum(!is.na(cp23_analysis_data$cancer_subgroup)),
    sum(!is.na(cp23_analysis_data$therapy_status_4level)),
    vapply(exposure_vars, function(ee) sum(!is.na(cp23_analysis_data[[ee]])), numeric(1)),
    sum(!is.na(cp23_analysis_data$therapy_line_number)),
    sum(!is.na(cp23_analysis_data$time_from_start_days))
  )
)

merge_summary <- cp23_analysis_data %>%
  summarise(
    n_cp23_rows = n(),
    n_unique_subjects = n_distinct(subject_id),
    n_feature_ok = sum(feature_ok == TRUE, na.rm = TRUE),
    n_with_disease_group = sum(!is.na(disease_group)),
    n_with_valid_age_for_model = sum(!is.na(age_for_model)),
    n_invalid_or_missing_age_for_model = sum(is.na(age_for_model)),
    n_with_sex = sum(!is.na(sex)),
    n_with_binary_sex = sum(!is.na(sex_binary)),
    n_model_ready_age_sex = sum(model_ready_age_sex == TRUE, na.rm = TRUE),
    n_channel_order_mismatch_files = sum(channel_order_mismatch_file == TRUE, na.rm = TRUE),
    n_marker_order_mismatch_files = sum(marker_order_mismatch_file == TRUE, na.rm = TRUE),
    n_marker_set_mismatch_files = sum(marker_set_mismatch_file == TRUE, na.rm = TRUE)
  )

disease_counts <- cp23_analysis_data %>%
  count(disease_group, name = "n") %>%
  arrange(disease_group)

model_ready_by_disease <- cp23_analysis_data %>%
  count(disease_group, model_ready_age_sex, name = "n") %>%
  arrange(disease_group, model_ready_age_sex)

age_summary_by_disease <- cp23_analysis_data %>%
  group_by(disease_group) %>%
  summarise(
    n_total = n(),
    n_valid_age = sum(!is.na(age_for_model)),
    mean_age = mean(age_for_model, na.rm = TRUE),
    sd_age = sd(age_for_model, na.rm = TRUE),
    median_age = median(age_for_model, na.rm = TRUE),
    q1_age = as.numeric(quantile(age_for_model, 0.25, na.rm = TRUE)),
    q3_age = as.numeric(quantile(age_for_model, 0.75, na.rm = TRUE)),
    min_age = min(age_for_model, na.rm = TRUE),
    max_age = max(age_for_model, na.rm = TRUE),
    .groups = "drop"
  )

sex_counts_by_disease <- cp23_analysis_data %>%
  count(disease_group, sex, name = "n") %>%
  arrange(disease_group, sex)

cd45_event_qc_by_disease <- cp23_analysis_data %>%
  count(disease_group, cd45_dump_low_event_qc_bin, name = "n") %>%
  arrange(disease_group, factor(cd45_dump_low_event_qc_bin, levels = c("<50", "50-99", "100-299", "300-999", ">=1000")))

hladr_myeloid_event_qc_by_disease <- cp23_analysis_data %>%
  count(disease_group, hladr_myeloid_event_qc_bin, name = "n") %>%
  arrange(disease_group, factor(hladr_myeloid_event_qc_bin, levels = c("<50", "50-99", "100-299", "300-999", ">=1000")))

technical_mismatch_by_disease <- cp23_analysis_data %>%
  count(
    disease_group,
    channel_order_mismatch_file,
    marker_order_mismatch_file,
    marker_set_mismatch_file,
    name = "n"
  ) %>%
  arrange(disease_group, desc(channel_order_mismatch_file))

selected_features_for_descriptive <- c(
  "pct_cd45_dump_low_within_total",
  "pct_cd14_mono_like_within_cd45_dump_low",
  "pct_cd14_hladr_low_within_cd45_dump_low",
  "pct_hladr_low_within_cd14_mono_like",
  "pct_cd33_hladr_myeloid_like_within_cd45_dump_low",
  "pct_cd15_gran_like_within_cd45_dump_low",
  "pct_cd14_cd206_pos_within_cd45_dump_low",
  "pct_cd14_cd169_pos_within_cd45_dump_low",
  "pct_cd14_cd9_pos_within_cd45_dump_low",
  "pct_cd14_cd84_pos_within_cd45_dump_low",
  "pct_fceri_hladr_apc_like_within_cd45_dump_low"
)

selected_features_for_descriptive <- selected_features_for_descriptive[
  selected_features_for_descriptive %in% names(cp23_analysis_data)
]

selected_feature_medians_by_disease <- cp23_analysis_data %>%
  filter(feature_ok == TRUE) %>%
  group_by(disease_group) %>%
  summarise(
    n = n(),
    across(
      all_of(selected_features_for_descriptive),
      ~ median(.x, na.rm = TRUE),
      .names = "median_{.col}"
    ),
    .groups = "drop"
  )

metadata_subjects_missing_cp23_fcs <- metadata_keep %>%
  anti_join(cp23_analysis_data %>% distinct(subject_id), by = "subject_id") %>%
  select(subject_id, disease_group, age_raw, age_for_model, sex, cancer_subgroup)

cp23_subjects_missing_metadata <- cp23_analysis_data %>%
  filter(is.na(disease_group) | is.na(age_raw) | is.na(sex)) %>%
  select(subject_id, file_name, disease_group, age_raw, age_for_model, sex) %>%
  arrange(subject_id)

# ------------------------------------------------------------
# 8. Save
# ------------------------------------------------------------

write_csv(cp23_analysis_data, file.path(out_dir, "SDY2583_CP23_analysis_data_STEP3A.csv"))
write_csv(metadata_source_summary, file.path(out_dir, "SDY2583_CP23_metadata_source_summary_STEP3A.csv"))
write_csv(merge_summary, file.path(out_dir, "SDY2583_CP23_merge_summary_STEP3A.csv"))
write_csv(disease_counts, file.path(out_dir, "SDY2583_CP23_disease_counts_STEP3A.csv"))
write_csv(model_ready_by_disease, file.path(out_dir, "SDY2583_CP23_model_ready_by_disease_STEP3A.csv"))
write_csv(age_summary_by_disease, file.path(out_dir, "SDY2583_CP23_age_summary_by_disease_STEP3A.csv"))
write_csv(sex_counts_by_disease, file.path(out_dir, "SDY2583_CP23_sex_counts_by_disease_STEP3A.csv"))
write_csv(cd45_event_qc_by_disease, file.path(out_dir, "SDY2583_CP23_cd45_event_QC_by_disease_STEP3A.csv"))
write_csv(hladr_myeloid_event_qc_by_disease, file.path(out_dir, "SDY2583_CP23_hladr_myeloid_event_QC_by_disease_STEP3A.csv"))
write_csv(technical_mismatch_by_disease, file.path(out_dir, "SDY2583_CP23_technical_mismatch_by_disease_STEP3A.csv"))
write_csv(selected_feature_medians_by_disease, file.path(out_dir, "SDY2583_CP23_selected_feature_medians_by_disease_STEP3A.csv"))
write_csv(clinical_variable_availability, file.path(out_dir, "SDY2583_CP23_clinical_variable_availability_STEP3A.csv"))
write_csv(metadata_subjects_missing_cp23_fcs, file.path(out_dir, "SDY2583_CP23_metadata_subjects_missing_CP23_FCS_STEP3A.csv"))
write_csv(cp23_subjects_missing_metadata, file.path(out_dir, "SDY2583_CP23_subjects_missing_metadata_STEP3A.csv"))

save(
  cp23_analysis_data,
  metadata_raw,
  metadata_keep,
  metadata_source_file,
  metadata_source_summary,
  clinical_variable_availability,
  merge_summary,
  disease_counts,
  model_ready_by_disease,
  age_summary_by_disease,
  sex_counts_by_disease,
  cd45_event_qc_by_disease,
  hladr_myeloid_event_qc_by_disease,
  technical_mismatch_by_disease,
  selected_feature_medians_by_disease,
  metadata_subjects_missing_cp23_fcs,
  cp23_subjects_missing_metadata,
  feature_dictionary_full,
  file = file.path(rdata_dir, "SDY2583_CP23_STEP3A_metadata_merge_age_QC.RData")
)

# ------------------------------------------------------------
# 9. Console output
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SDY2583 CP23 STEP 3A COMPLETE: METADATA MERGE / AGE QC\n")
cat("============================================================\n")

cat("\nMetadata source summary:\n")
print(as.data.frame(metadata_source_summary), row.names = FALSE)

cat("\nMerge summary:\n")
print(as.data.frame(merge_summary), row.names = FALSE)

cat("\nDisease counts:\n")
print(as.data.frame(disease_counts), row.names = FALSE)

cat("\nModel-ready by disease:\n")
print(as.data.frame(model_ready_by_disease), row.names = FALSE)

cat("\nAge summary by disease:\n")
print(as.data.frame(age_summary_by_disease), row.names = FALSE)

cat("\nSex counts by disease:\n")
print(as.data.frame(sex_counts_by_disease), row.names = FALSE)

cat("\nCD45 dump-low event QC by disease:\n")
print(as.data.frame(cd45_event_qc_by_disease), row.names = FALSE)

cat("\nCD33+HLA-DR+ myeloid/APC-like event QC by disease:\n")
print(as.data.frame(hladr_myeloid_event_qc_by_disease), row.names = FALSE)

cat("\nTechnical mismatch by disease:\n")
print(as.data.frame(technical_mismatch_by_disease), row.names = FALSE)

cat("\nSelected feature medians by disease, unadjusted descriptive only:\n")
print(as.data.frame(selected_feature_medians_by_disease), row.names = FALSE)

cat("\nClinical variable availability:\n")
print(as.data.frame(clinical_variable_availability), row.names = FALSE)

cat("\nMetadata subjects missing CP23 FCS:\n")
print(as.data.frame(metadata_subjects_missing_cp23_fcs), row.names = FALSE)

cat("\nCP23 subjects missing metadata / key fields:\n")
print(as.data.frame(cp23_subjects_missing_metadata), row.names = FALSE)

cat("\nOutputs saved in:\n")
print(out_dir)

cat("============================================================\n")
