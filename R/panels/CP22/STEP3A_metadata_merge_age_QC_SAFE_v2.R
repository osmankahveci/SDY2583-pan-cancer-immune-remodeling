# ============================================================
# SDY2583 CP22
# STEP 3A SAFE v2: Metadata merge + age quality control
# Fix: robust age extraction including age_for_model
# ============================================================

rm(list = ls())
source(file.path(Sys.getenv("SDY2583_REPO_ROOT", unset = "."), "R", "shared", "bootstrap.R"))

cran_pkgs <- c("dplyr", "readr", "stringr", "tibble", "purrr")
for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(purrr)
})

analysis_dir <- sd_analysis_dir("CP22")
step2_rdata <- file.path(analysis_dir, "11_RData", "SDY2583_CP22_STEP2_feature_extraction_Bcell_humoral.RData")
if (!file.exists(step2_rdata)) {
  step2_rdata <- file.path(analysis_dir, "11_RData", "SDY2583_CP22_STEP2_feature_extraction.RData")
}
if (!file.exists(step2_rdata)) stop("CP22 Step 2 RData bulunamadı.")

out_dir <- file.path(analysis_dir, "03_metadata_merge_age_QC")
rdata_dir <- file.path(analysis_dir, "11_RData")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rdata_dir, recursive = TRUE, showWarnings = FALSE)

load(step2_rdata)
if (!exists("cp22_feature_table")) stop("Step 2 RData içinde cp22_feature_table bulunamadı.")
cp22_features <- cp22_feature_table

metadata_candidates <- path.expand(c(
  file.path(sd_metadata_dir(), "SDY2583_CP7_STEP3A_metadata_merge_age_QC.RData"),
  file.path(sd_analysis_dir("CP7"), "11_RData", "SDY2583_CP7_STEP3A_metadata_merge_age_QC.RData"),
  file.path(sd_integration_rdata_dir(), "SDY2583_CP7_STEP3A_metadata_merge_age_QC.RData"),
  file.path(sd_integration_rdata_dir(), "SDY2583_CP7_STEP4_composite_scores.RData"),
  file.path(sd_integrated_dir(), "SDY2583_integrated_clinical_immune_score_matrix_ALL6_v4_CP24_ID_REPAIR.csv")
))
metadata_candidates <- metadata_candidates[file.exists(metadata_candidates)]
if (length(metadata_candidates) == 0) {
  stop("Metadata kaynağı bulunamadı. CP7 Step3A RData veya integrated clinical immune matrix dosyasını beklenen klasöre koy.")
}
selected_metadata_file <- metadata_candidates[1]

cat("\nSelected metadata source:\n")
print(selected_metadata_file)

load_rdata_objects <- function(file_path) {
  e <- new.env(parent = emptyenv())
  nm <- load(file_path, envir = e)
  list(env = e, objects = nm)
}

is_good_metadata_candidate <- function(obj) {
  if (!is.data.frame(obj)) return(FALSE)
  cn <- names(obj)
  has_id <- any(c("subject_id", "participant_id", "file_name") %in% cn)
  has_group <- any(grepl("disease|group|condition", cn, ignore.case = TRUE))
  has_age <- any(grepl("age_for_model|age_year|age_clinical|^age$|age", cn, ignore.case = TRUE))
  has_sex <- any(grepl("^sex$|sex_clinical|gender", cn, ignore.case = TRUE))
  has_id && (has_group || has_age || has_sex)
}

extract_metadata_from_rdata <- function(file_path) {
  loaded <- load_rdata_objects(file_path)
  e <- loaded$env
  obj_names <- loaded$objects
  candidates <- list()
  for (nm in obj_names) {
    obj <- get(nm, envir = e)
    if (is_good_metadata_candidate(obj)) candidates[[nm]] <- obj
  }
  if (length(candidates) == 0) return(NULL)
  preferred <- c("cp7_analysis_data", "analysis_df", "analysis_data", "cp7_scores_data", "metadata", "clinical_master")
  for (p in preferred) {
    if (p %in% names(candidates)) return(list(data = candidates[[p]], object_name = p))
  }
  sizes <- vapply(candidates, nrow, numeric(1))
  pick <- names(which.max(sizes))
  list(data = candidates[[pick]], object_name = pick)
}

if (grepl("\\.RData$", selected_metadata_file, ignore.case = TRUE)) {
  md <- extract_metadata_from_rdata(selected_metadata_file)
  if (is.null(md)) stop("RData içinde uygun metadata data.frame bulunamadı.")
  metadata_raw <- md$data
  selected_metadata_object <- md$object_name
} else if (grepl("\\.csv$", selected_metadata_file, ignore.case = TRUE)) {
  metadata_raw <- readr::read_csv(selected_metadata_file, show_col_types = FALSE)
  selected_metadata_object <- "csv_file"
} else {
  stop("Metadata file type not supported.")
}

cat("\nSelected metadata object:\n")
print(selected_metadata_object)

get_first_existing_col <- function(df, possible_cols) {
  hit <- possible_cols[possible_cols %in% names(df)]
  if (length(hit) == 0) return(list(values = rep(NA, nrow(df)), source_col = NA_character_))
  list(values = df[[hit[1]]], source_col = hit[1])
}

subject_obj <- get_first_existing_col(metadata_raw, c("subject_id", "participant_id", "Subject_ID", "subject", "id"))
disease_obj <- get_first_existing_col(metadata_raw, c("disease_group", "disease_group_clinical", "group", "condition", "diagnosis_group"))
sex_obj <- get_first_existing_col(metadata_raw, c("sex", "sex_clinical", "gender", "Sex"))
age_obj <- get_first_existing_col(metadata_raw, c("age_for_model", "age_years", "age_clinical", "age", "Age", "age_at_sampling", "age_at_draw", "age_years_raw"))

metadata_variable_sources <- tibble(
  variable = c("subject_id", "disease_group", "sex", "age"),
  source_column = c(subject_obj$source_col, disease_obj$source_col, sex_obj$source_col, age_obj$source_col)
)

cat("\nMetadata variable sources:\n")
print(as.data.frame(metadata_variable_sources), row.names = FALSE)

metadata_core <- tibble(
  subject_id = as.character(subject_obj$values),
  disease_group_raw = as.character(disease_obj$values),
  sex_raw = as.character(sex_obj$values),
  age_source_value = suppressWarnings(as.numeric(age_obj$values))
) %>%
  filter(!is.na(subject_id), subject_id != "") %>%
  distinct(subject_id, .keep_all = TRUE) %>%
  mutate(
    disease_group = case_when(
      disease_group_raw %in% c("Cancer patient", "Cancer", "cancer", "patient", "Patient") ~ "Cancer patient",
      disease_group_raw %in% c("Healthy control", "Healthy", "healthy", "control", "Control") ~ "Healthy control",
      str_detect(disease_group_raw, regex("cancer|patient", ignore_case = TRUE)) ~ "Cancer patient",
      str_detect(disease_group_raw, regex("healthy|control", ignore_case = TRUE)) ~ "Healthy control",
      TRUE ~ disease_group_raw
    ),
    sex = case_when(
      sex_raw %in% c("Female", "F", "female", "f") ~ "Female",
      sex_raw %in% c("Male", "M", "male", "m") ~ "Male",
      TRUE ~ sex_raw
    ),
    age_years = age_source_value,
    age_for_model = ifelse(age_years < 18 | age_years > 100, NA_real_, age_years)
  )

metadata_age_diagnostic <- metadata_core %>%
  summarise(
    n_metadata_subjects = n(),
    n_age_source_nonmissing = sum(!is.na(age_source_value)),
    n_age_for_model_nonmissing = sum(!is.na(age_for_model)),
    min_age_for_model = ifelse(sum(!is.na(age_for_model)) > 0, min(age_for_model, na.rm = TRUE), NA_real_),
    max_age_for_model = ifelse(sum(!is.na(age_for_model)) > 0, max(age_for_model, na.rm = TRUE), NA_real_)
  )

clinical_file <- file.path(
  sd_integrated_dir(),
  "SDY2583_integrated_clinical_immune_score_matrix_ALL6_v4_CP24_ID_REPAIR.csv"
)
clinical_annotation <- NULL

if (file.exists(clinical_file)) {
  clinical_raw <- readr::read_csv(clinical_file, show_col_types = FALSE)
  clinical_cols_keep <- intersect(c(
    "subject_id", "disease_group_clinical", "age_clinical", "sex_clinical",
    "cancer_subgroup", "therapy_status_4level",
    "chemotherapy", "targeted_therapy", "any_immunotherapy", "ici_immunotherapy",
    "endocrine_hormonal", "adc", "experimental", "radiotherapy",
    "therapy_line_number", "therapy_line_raw", "time_from_start_days", "time_from_start_raw"
  ), names(clinical_raw))
  clinical_annotation <- clinical_raw %>%
    select(all_of(clinical_cols_keep)) %>%
    distinct(subject_id, .keep_all = TRUE)
}

cp22_analysis_data <- cp22_features %>%
  left_join(metadata_core, by = "subject_id")

if (!is.null(clinical_annotation)) {
  clinical_annotation_no_dupes <- clinical_annotation %>%
    select(-any_of(c("disease_group_clinical", "age_clinical", "sex_clinical")))
  cp22_analysis_data <- cp22_analysis_data %>%
    left_join(clinical_annotation_no_dupes, by = "subject_id")
}

cp22_analysis_data <- cp22_analysis_data %>%
  mutate(
    disease_group = factor(disease_group, levels = c("Healthy control", "Cancer patient")),
    sex = factor(sex),
    model_ready_age_sex = feature_ok == TRUE & !is.na(disease_group) & !is.na(age_for_model) & !is.na(sex)
  )

safe_min <- function(x) if (sum(!is.na(x)) == 0) NA_real_ else min(x, na.rm = TRUE)
safe_max <- function(x) if (sum(!is.na(x)) == 0) NA_real_ else max(x, na.rm = TRUE)
safe_quantile <- function(x, p) if (sum(!is.na(x)) == 0) NA_real_ else as.numeric(quantile(x, p, na.rm = TRUE))

merge_summary <- tibble(
  n_cp22_rows = nrow(cp22_analysis_data),
  n_unique_subjects = n_distinct(cp22_analysis_data$subject_id),
  n_feature_ok = sum(cp22_analysis_data$feature_ok == TRUE, na.rm = TRUE),
  n_with_disease_group = sum(!is.na(cp22_analysis_data$disease_group)),
  n_with_valid_age_for_model = sum(!is.na(cp22_analysis_data$age_for_model)),
  n_invalid_age_for_model = sum(is.na(cp22_analysis_data$age_for_model)),
  n_with_sex = sum(!is.na(cp22_analysis_data$sex)),
  n_model_ready_age_sex = sum(cp22_analysis_data$model_ready_age_sex, na.rm = TRUE),
  n_with_cancer_subgroup = if ("cancer_subgroup" %in% names(cp22_analysis_data)) sum(!is.na(cp22_analysis_data$cancer_subgroup)) else NA_integer_,
  n_with_therapy_status = if ("therapy_status_4level" %in% names(cp22_analysis_data)) sum(!is.na(cp22_analysis_data$therapy_status_4level)) else NA_integer_
)

disease_counts <- cp22_analysis_data %>% count(disease_group, sort = FALSE)
model_ready_by_disease <- cp22_analysis_data %>% count(disease_group, model_ready_age_sex, sort = FALSE)

age_summary_by_disease <- cp22_analysis_data %>%
  group_by(disease_group) %>%
  summarise(
    n_total = n(),
    n_valid_age = sum(!is.na(age_for_model)),
    mean_age = ifelse(n_valid_age > 0, mean(age_for_model, na.rm = TRUE), NA_real_),
    sd_age = ifelse(n_valid_age > 1, sd(age_for_model, na.rm = TRUE), NA_real_),
    median_age = ifelse(n_valid_age > 0, median(age_for_model, na.rm = TRUE), NA_real_),
    q1_age = safe_quantile(age_for_model, 0.25),
    q3_age = safe_quantile(age_for_model, 0.75),
    min_age = safe_min(age_for_model),
    max_age = safe_max(age_for_model),
    .groups = "drop"
  )

sex_counts_by_disease <- cp22_analysis_data %>% count(disease_group, sex, sort = FALSE)

b_event_qc_by_disease <- cp22_analysis_data %>%
  mutate(
    b_event_qc_bin = case_when(
      is.na(n_cd19_b) ~ "missing",
      n_cd19_b < 50 ~ "<50",
      n_cd19_b < 100 ~ "50-99",
      n_cd19_b < 300 ~ "100-299",
      n_cd19_b < 1000 ~ "300-999",
      TRUE ~ ">=1000"
    )
  ) %>%
  count(disease_group, b_event_qc_bin, sort = FALSE)

dump_variant_by_disease <- cp22_analysis_data %>% count(disease_group, dump_marker_variant, sort = FALSE)

write_csv(cp22_analysis_data, file.path(out_dir, "SDY2583_CP22_analysis_data_STEP3A_metadata_age_QC.csv"))
write_csv(metadata_variable_sources, file.path(out_dir, "SDY2583_CP22_metadata_variable_sources_STEP3A.csv"))
write_csv(metadata_age_diagnostic, file.path(out_dir, "SDY2583_CP22_metadata_age_diagnostic_STEP3A.csv"))
write_csv(merge_summary, file.path(out_dir, "SDY2583_CP22_metadata_merge_summary_STEP3A.csv"))
write_csv(disease_counts, file.path(out_dir, "SDY2583_CP22_disease_counts_STEP3A.csv"))
write_csv(model_ready_by_disease, file.path(out_dir, "SDY2583_CP22_model_ready_by_disease_STEP3A.csv"))
write_csv(age_summary_by_disease, file.path(out_dir, "SDY2583_CP22_age_summary_by_disease_STEP3A.csv"))
write_csv(sex_counts_by_disease, file.path(out_dir, "SDY2583_CP22_sex_counts_by_disease_STEP3A.csv"))
write_csv(b_event_qc_by_disease, file.path(out_dir, "SDY2583_CP22_Bcell_event_QC_by_disease_STEP3A.csv"))
write_csv(dump_variant_by_disease, file.path(out_dir, "SDY2583_CP22_dump_variant_by_disease_STEP3A.csv"))

save(
  cp22_analysis_data, cp22_features, metadata_core, metadata_variable_sources,
  metadata_age_diagnostic, selected_metadata_file, selected_metadata_object,
  merge_summary, disease_counts, model_ready_by_disease, age_summary_by_disease,
  sex_counts_by_disease, b_event_qc_by_disease, dump_variant_by_disease,
  file = file.path(rdata_dir, "SDY2583_CP22_STEP3A_metadata_merge_age_QC.RData")
)

cat("\n============================================================\n")
cat("SDY2583 CP22 STEP 3A v2 COMPLETE: METADATA MERGE + AGE QC\n")
cat("============================================================\n")

cat("\nSelected metadata source:\n")
print(selected_metadata_file)

cat("\nSelected metadata object:\n")
print(selected_metadata_object)

cat("\nMetadata variable sources:\n")
print(as.data.frame(metadata_variable_sources), row.names = FALSE)

cat("\nMetadata age diagnostic:\n")
print(as.data.frame(metadata_age_diagnostic), row.names = FALSE)

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

cat("\nB-cell event QC by disease:\n")
print(as.data.frame(b_event_qc_by_disease), row.names = FALSE)

cat("\nDump marker variant by disease:\n")
print(as.data.frame(dump_variant_by_disease), row.names = FALSE)

cat("\nFiles saved in:\n")
print(out_dir)

cat("============================================================\n")
