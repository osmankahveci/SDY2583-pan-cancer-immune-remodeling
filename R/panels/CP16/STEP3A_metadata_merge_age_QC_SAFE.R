# ============================================================
# SDY2583 CP16
# STEP 3A SAFE: Metadata merge / age quality control
#
# Input:
#   outputs/CP16/11_RData/
#     SDY2583_CP16_STEP2_feature_extraction_APC_DC_myeloid.RData
#
# Metadata source priority:
#   1) Integrated clinical immune matrix ALL8 with CP10
#   2) Integrated clinical immune matrix ALL7 with CP22
#   3) Any integrated clinical immune matrix CSV
#   4) Existing SDY2583 metadata RData if available
#
# Output:
#   outputs/CP16/03_metadata_merge_age_QC
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

analysis_dir <- sd_analysis_dir("CP16")
rdata_dir <- file.path(analysis_dir, "11_RData")
out_dir <- file.path(analysis_dir, "03_metadata_merge_age_QC")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rdata_dir, recursive = TRUE, showWarnings = FALSE)

step2_rdata <- file.path(rdata_dir, "SDY2583_CP16_STEP2_feature_extraction_APC_DC_myeloid.RData")

if (!file.exists(step2_rdata)) {
  stop("Step 2 RData bulunamadı: ", step2_rdata)
}

load(step2_rdata)

# IMPORTANT: loaded RData may contain old path variables. Reset paths after load.
analysis_dir <- sd_analysis_dir("CP16")
rdata_dir <- file.path(analysis_dir, "11_RData")
out_dir <- file.path(analysis_dir, "03_metadata_merge_age_QC")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!exists("cp16_features")) stop("cp16_features bulunamadı.")
if (!exists("feature_dictionary")) stop("feature_dictionary bulunamadı.")

# ------------------------------------------------------------
# 3. Helper functions
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
    y %in% c("Cancer patient", "Cancer", "cancer", "Patient", "patient", "Cancer patient ") ~ "Cancer patient",
    str_detect(y, regex("healthy|control", ignore_case = TRUE)) ~ "Healthy control",
    str_detect(y, regex("cancer|patient|tumou?r", ignore_case = TRUE)) ~ "Cancer patient",
    TRUE ~ y
  )
}

standardize_sex <- function(x) {
  y <- as.character(x)
  case_when(
    y %in% c("Female", "F", "female", "f", "Woman", "woman") ~ "Female",
    y %in% c("Male", "M", "male", "m", "Man", "man") ~ "Male",
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

# ------------------------------------------------------------
# 4. Find metadata source
# ------------------------------------------------------------

matrix_dir <- sd_integrated_dir()

candidate_files <- c(
  file.path(matrix_dir, "SDY2583_integrated_clinical_immune_score_matrix_ALL8_with_CP10.csv"),
  file.path(matrix_dir, "SDY2583_integrated_clinical_immune_score_matrix_ALL7_with_CP22.csv")
)

if (dir.exists(matrix_dir)) {
  all_matrix_csv <- list.files(
    matrix_dir,
    pattern = "\\.csv$",
    full.names = TRUE,
    recursive = TRUE
  )

  all_matrix_csv <- all_matrix_csv[
    grepl("integrated|matrix|ALL", basename(all_matrix_csv), ignore.case = TRUE)
  ]

  candidate_files <- unique(c(candidate_files, all_matrix_csv))
}

candidate_files <- candidate_files[file.exists(candidate_files)]

inspect_metadata_file <- function(fp) {
  dat <- tryCatch(readr::read_csv(fp, show_col_types = FALSE), error = function(e) NULL)

  if (is.null(dat)) {
    return(tibble(
      file_path = fp,
      readable = FALSE,
      n_rows = NA_integer_,
      n_cols = NA_integer_,
      has_subject_id = FALSE,
      has_disease = FALSE,
      has_age = FALSE,
      has_sex = FALSE,
      has_clinical = FALSE,
      priority = -Inf
    ))
  }

  nm <- names(dat)

  has_disease <- any(c("disease_group_model", "disease_group_clinical", "disease_group", "group", "condition") %in% nm)
  has_age <- any(c("age_for_model", "age_clinical", "age_years", "age") %in% nm)
  has_sex <- any(c("sex_for_clinical_model", "sex_clinical", "sex", "gender") %in% nm)
  has_clinical <- any(c("cancer_subgroup", "cancer_subgroup_model", "therapy_status_4level", "therapy_status_model", "chemotherapy") %in% nm)
  score_cols <- grep("^CP[0-9]+_.*_score$", nm, value = TRUE)

  priority <- 0
  priority <- priority + ifelse("subject_id" %in% nm, 1000, 0)
  priority <- priority + ifelse(nrow(dat) >= 800 & nrow(dat) <= 900, 300, 0)
  priority <- priority + ifelse(has_disease, 100, 0)
  priority <- priority + ifelse(has_age, 100, 0)
  priority <- priority + ifelse(has_sex, 100, 0)
  priority <- priority + ifelse(has_clinical, 100, 0)
  priority <- priority + length(score_cols)
  priority <- priority + ifelse(grepl("ALL8_with_CP10", basename(fp), ignore.case = TRUE), 300, 0)
  priority <- priority + ifelse(grepl("ALL7_with_CP22", basename(fp), ignore.case = TRUE), 200, 0)
  priority <- priority - ifelse(grepl("summary|counts|results|dictionary|manifest|annotation|subgroup|therapy|exposure", basename(fp), ignore.case = TRUE), 500, 0)

  tibble(
    file_path = fp,
    readable = TRUE,
    n_rows = nrow(dat),
    n_cols = ncol(dat),
    has_subject_id = "subject_id" %in% nm,
    has_disease = has_disease,
    has_age = has_age,
    has_sex = has_sex,
    has_clinical = has_clinical,
    n_score_cols = length(score_cols),
    priority = priority
  )
}

candidate_inventory <- bind_rows(lapply(candidate_files, inspect_metadata_file)) %>%
  arrange(desc(priority), desc(n_score_cols), desc(n_cols))

write_csv(
  candidate_inventory,
  file.path(out_dir, "SDY2583_CP16_metadata_candidate_inventory_STEP3A.csv")
)

valid_candidates <- candidate_inventory %>%
  filter(
    readable == TRUE,
    has_subject_id == TRUE,
    has_disease == TRUE,
    has_age == TRUE,
    has_sex == TRUE,
    n_rows >= 800,
    n_rows <= 900
  )

if (nrow(valid_candidates) == 0) {
  cat("\nMetadata candidate inventory:\n")
  print(as.data.frame(candidate_inventory), row.names = FALSE)
  stop("Uygun metadata/integrated matrix CSV bulunamadı.")
}

metadata_source_file <- valid_candidates$file_path[1]
metadata_source_label <- "integrated_clinical_immune_matrix_csv"

metadata_raw <- readr::read_csv(metadata_source_file, show_col_types = FALSE)

# ------------------------------------------------------------
# 5. Standardize metadata
# ------------------------------------------------------------

disease_col <- first_existing_col(metadata_raw, c("disease_group_model", "disease_group_clinical", "disease_group", "group", "condition"))
age_col <- first_existing_col(metadata_raw, c("age_for_model", "age_clinical", "age_years", "age"))
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
  if (ee %in% names(metadata_standardized)) {
    metadata_standardized[[ee]] <- as_01(metadata_standardized[[ee]])
  } else {
    metadata_standardized[[ee]] <- NA_integer_
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
    chemotherapy,
    targeted_therapy,
    any_immunotherapy,
    ici_immunotherapy,
    endocrine_hormonal,
    adc,
    experimental,
    radiotherapy,
    therapy_line_number,
    time_from_start_days
  ) %>%
  distinct(subject_id, .keep_all = TRUE)

# ------------------------------------------------------------
# 6. Merge CP16 features with metadata
# ------------------------------------------------------------

cp16_analysis_data <- cp16_features %>%
  mutate(subject_id = as.character(subject_id)) %>%
  left_join(metadata_keep, by = "subject_id") %>%
  mutate(
    model_ready_age_sex = feature_ok == TRUE &
      !is.na(disease_group) &
      !is.na(age_for_model) &
      !is.na(sex),
    model_ready_binary_sex = model_ready_age_sex == TRUE &
      !is.na(sex_binary),
    event_QC_cd45_ge_300 = feature_ok == TRUE &
      !is.na(n_cd45_dump_low) &
      n_cd45_dump_low >= 300,
    event_QC_hladr_ge_300 = feature_ok == TRUE &
      !is.na(n_hladr_apc_core) &
      n_hladr_apc_core >= 300,
    event_QC_cd45_ge_1000 = feature_ok == TRUE &
      !is.na(n_cd45_dump_low) &
      n_cd45_dump_low >= 1000,
    event_QC_hladr_ge_1000 = feature_ok == TRUE &
      !is.na(n_hladr_apc_core) &
      n_hladr_apc_core >= 1000
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

merge_summary <- tibble(
  n_cp16_rows = nrow(cp16_analysis_data),
  n_unique_subjects = n_distinct(cp16_analysis_data$subject_id),
  n_feature_ok = sum(cp16_analysis_data$feature_ok == TRUE, na.rm = TRUE),
  n_with_disease_group = sum(!is.na(cp16_analysis_data$disease_group)),
  n_with_valid_age_for_model = sum(!is.na(cp16_analysis_data$age_for_model)),
  n_invalid_or_missing_age_for_model = sum(cp16_analysis_data$invalid_or_missing_age_for_model == TRUE, na.rm = TRUE),
  n_with_sex = sum(!is.na(cp16_analysis_data$sex)),
  n_with_binary_sex = sum(!is.na(cp16_analysis_data$sex_binary)),
  n_model_ready_age_sex = sum(cp16_analysis_data$model_ready_age_sex == TRUE, na.rm = TRUE),
  n_channel_order_mismatch_files = sum(cp16_analysis_data$channel_order_mismatch_file == TRUE, na.rm = TRUE),
  n_marker_order_mismatch_files = sum(cp16_analysis_data$marker_order_mismatch_file == TRUE, na.rm = TRUE),
  n_marker_set_mismatch_files = sum(cp16_analysis_data$marker_set_mismatch_file == TRUE, na.rm = TRUE)
)

missing_cp16_subjects_from_metadata <- cp16_analysis_data %>%
  filter(is.na(disease_group) | is.na(age_raw) | is.na(sex)) %>%
  select(subject_id, file_name, disease_group, age_raw, sex) %>%
  arrange(subject_id)

metadata_subjects_missing_cp16 <- metadata_keep %>%
  anti_join(cp16_analysis_data %>% select(subject_id), by = "subject_id") %>%
  select(subject_id, disease_group, age_raw, age_for_model, sex, cancer_subgroup) %>%
  arrange(subject_id)

disease_counts <- cp16_analysis_data %>%
  count(disease_group, name = "n") %>%
  arrange(disease_group)

model_ready_by_disease <- cp16_analysis_data %>%
  count(disease_group, model_ready_age_sex, name = "n") %>%
  arrange(disease_group, model_ready_age_sex)

age_summary_by_disease <- cp16_analysis_data %>%
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

sex_counts_by_disease <- cp16_analysis_data %>%
  count(disease_group, sex, name = "n") %>%
  arrange(disease_group, sex)

cd45_event_qc_by_disease <- cp16_analysis_data %>%
  count(disease_group, cd45_dump_low_event_qc_bin, name = "n") %>%
  arrange(disease_group, cd45_dump_low_event_qc_bin)

hladr_event_qc_by_disease <- cp16_analysis_data %>%
  count(disease_group, hladr_apc_core_event_qc_bin, name = "n") %>%
  arrange(disease_group, hladr_apc_core_event_qc_bin)

technical_mismatch_by_disease <- cp16_analysis_data %>%
  count(
    disease_group,
    channel_order_mismatch_file,
    marker_order_mismatch_file,
    marker_set_mismatch_file,
    name = "n"
  ) %>%
  arrange(disease_group, desc(channel_order_mismatch_file))

selected_feature_medians_by_disease <- cp16_analysis_data %>%
  group_by(disease_group) %>%
  summarise(
    n = n(),
    median_pct_cd45_dump_low_within_total = median(pct_cd45_dump_low_within_total, na.rm = TRUE),
    median_pct_hladr_apc_core_within_cd45_dump_low = median(pct_hladr_apc_core_within_cd45_dump_low, na.rm = TRUE),
    median_pct_cd11c_hladr_apc_like_within_cd45_dump_low = median(pct_cd11c_hladr_apc_like_within_cd45_dump_low, na.rm = TRUE),
    median_pct_cd14_mono_like_within_cd45_dump_low = median(pct_cd14_mono_like_within_cd45_dump_low, na.rm = TRUE),
    median_pct_cd14_hladr_low_within_cd45_dump_low = median(pct_cd14_hladr_low_within_cd45_dump_low, na.rm = TRUE),
    median_pct_cd1c_cd14low_cdc2_like_within_cd45_dump_low = median(pct_cd1c_cd14low_cdc2_like_within_cd45_dump_low, na.rm = TRUE),
    median_pct_cd141_clec9a_cdc1_like_within_cd45_dump_low = median(pct_cd141_clec9a_cdc1_like_within_cd45_dump_low, na.rm = TRUE),
    median_pct_cd123_pdc_like_within_cd45_dump_low = median(pct_cd123_pdc_like_within_cd45_dump_low, na.rm = TRUE),
    median_pct_cd13_cd11c_hladr_apc_like_within_cd45_dump_low = median(pct_cd13_cd11c_hladr_apc_like_within_cd45_dump_low, na.rm = TRUE),
    .groups = "drop"
  )

clinical_variable_availability <- tibble(
  variable = c(
    "cancer_subgroup",
    "therapy_status_4level",
    "chemotherapy",
    "targeted_therapy",
    "any_immunotherapy",
    "ici_immunotherapy",
    "therapy_line_number",
    "time_from_start_days"
  ),
  present = c(
    "cancer_subgroup" %in% names(cp16_analysis_data),
    "therapy_status_4level" %in% names(cp16_analysis_data),
    "chemotherapy" %in% names(cp16_analysis_data),
    "targeted_therapy" %in% names(cp16_analysis_data),
    "any_immunotherapy" %in% names(cp16_analysis_data),
    "ici_immunotherapy" %in% names(cp16_analysis_data),
    "therapy_line_number" %in% names(cp16_analysis_data),
    "time_from_start_days" %in% names(cp16_analysis_data)
  ),
  n_nonmissing = c(
    sum(!is.na(cp16_analysis_data$cancer_subgroup)),
    sum(!is.na(cp16_analysis_data$therapy_status_4level)),
    sum(!is.na(cp16_analysis_data$chemotherapy)),
    sum(!is.na(cp16_analysis_data$targeted_therapy)),
    sum(!is.na(cp16_analysis_data$any_immunotherapy)),
    sum(!is.na(cp16_analysis_data$ici_immunotherapy)),
    sum(!is.na(cp16_analysis_data$therapy_line_number)),
    sum(!is.na(cp16_analysis_data$time_from_start_days))
  )
)

# ------------------------------------------------------------
# 8. Save outputs
# ------------------------------------------------------------

write_csv(cp16_analysis_data, file.path(out_dir, "SDY2583_CP16_analysis_data_STEP3A.csv"))
write_csv(metadata_source_summary, file.path(out_dir, "SDY2583_CP16_metadata_source_summary_STEP3A.csv"))
write_csv(merge_summary, file.path(out_dir, "SDY2583_CP16_merge_summary_STEP3A.csv"))
write_csv(missing_cp16_subjects_from_metadata, file.path(out_dir, "SDY2583_CP16_missing_or_incomplete_metadata_subjects_STEP3A.csv"))
write_csv(metadata_subjects_missing_cp16, file.path(out_dir, "SDY2583_CP16_metadata_subjects_missing_CP16_FCS_STEP3A.csv"))
write_csv(disease_counts, file.path(out_dir, "SDY2583_CP16_disease_counts_STEP3A.csv"))
write_csv(model_ready_by_disease, file.path(out_dir, "SDY2583_CP16_model_ready_by_disease_STEP3A.csv"))
write_csv(age_summary_by_disease, file.path(out_dir, "SDY2583_CP16_age_summary_by_disease_STEP3A.csv"))
write_csv(sex_counts_by_disease, file.path(out_dir, "SDY2583_CP16_sex_counts_by_disease_STEP3A.csv"))
write_csv(cd45_event_qc_by_disease, file.path(out_dir, "SDY2583_CP16_CD45_event_QC_by_disease_STEP3A.csv"))
write_csv(hladr_event_qc_by_disease, file.path(out_dir, "SDY2583_CP16_HLADR_event_QC_by_disease_STEP3A.csv"))
write_csv(technical_mismatch_by_disease, file.path(out_dir, "SDY2583_CP16_technical_mismatch_by_disease_STEP3A.csv"))
write_csv(selected_feature_medians_by_disease, file.path(out_dir, "SDY2583_CP16_selected_feature_medians_by_disease_STEP3A.csv"))
write_csv(clinical_variable_availability, file.path(out_dir, "SDY2583_CP16_clinical_variable_availability_STEP3A.csv"))

save(
  cp16_analysis_data,
  metadata_source_summary,
  merge_summary,
  missing_cp16_subjects_from_metadata,
  metadata_subjects_missing_cp16,
  disease_counts,
  model_ready_by_disease,
  age_summary_by_disease,
  sex_counts_by_disease,
  cd45_event_qc_by_disease,
  hladr_event_qc_by_disease,
  technical_mismatch_by_disease,
  selected_feature_medians_by_disease,
  clinical_variable_availability,
  feature_dictionary,
  thresholds,
  file = file.path(rdata_dir, "SDY2583_CP16_STEP3A_metadata_merge_age_QC.RData")
)

# ------------------------------------------------------------
# 9. Console output
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SDY2583 CP16 STEP 3A COMPLETE: METADATA MERGE / AGE QC\n")
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

cat("\nHLA-DR APC-like core event QC by disease:\n")
print(as.data.frame(hladr_event_qc_by_disease), row.names = FALSE)

cat("\nTechnical mismatch by disease:\n")
print(as.data.frame(technical_mismatch_by_disease), row.names = FALSE)

cat("\nSelected feature medians by disease, unadjusted descriptive only:\n")
print(as.data.frame(selected_feature_medians_by_disease), row.names = FALSE)

cat("\nClinical variable availability:\n")
print(as.data.frame(clinical_variable_availability), row.names = FALSE)

cat("\nMetadata subjects missing CP16 FCS:\n")
print(as.data.frame(metadata_subjects_missing_cp16), row.names = FALSE)

cat("\nOutputs saved in:\n")
print(out_dir)

cat("============================================================\n")
