# ============================================================
# SDY2583 CP10
# STEP 3A SAFE: Metadata merge + age quality control
#
# Input:
#   outputs/CP10/11_RData/
#     SDY2583_CP10_STEP2_feature_extraction_myeloid_granulocyte.RData
#
# Metadata source priority:
#   1) Integrated clinical immune matrix ALL7 with CP22
#   2) Integrated clinical immune matrix ALL6/ALL7 subject-level CSV
#   3) SDY2583_MASTER_METADATA cp7_analysis_data RData
#
# Output:
#   cp10_analysis_data
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

analysis_dir <- sd_analysis_dir("CP10")
rdata_dir <- file.path(analysis_dir, "11_RData")
out_dir <- file.path(analysis_dir, "03_metadata_merge_age_QC")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rdata_dir, recursive = TRUE, showWarnings = FALSE)

step2_rdata <- file.path(
  rdata_dir,
  "SDY2583_CP10_STEP2_feature_extraction_myeloid_granulocyte.RData"
)

if (!file.exists(step2_rdata)) {
  stop("CP10 Step 2 RData bulunamadı: ", step2_rdata)
}

load(step2_rdata)

if (!exists("cp10_features")) {
  stop("Step 2 RData içinde cp10_features bulunamadı.")
}

# ------------------------------------------------------------
# 3. Helper functions
# ------------------------------------------------------------

first_existing_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
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

standardize_metadata <- function(meta_raw) {

  if (!("subject_id" %in% names(meta_raw))) {
    stop("Metadata içinde subject_id yok.")
  }

  disease_col <- first_existing_col(
    meta_raw,
    c("disease_group_model", "disease_group_clinical", "disease_group", "group", "condition")
  )

  age_col <- first_existing_col(
    meta_raw,
    c("age_for_model", "age_clinical", "age_years", "age")
  )

  sex_col <- first_existing_col(
    meta_raw,
    c("sex_for_clinical_model", "sex_clinical", "sex", "gender")
  )

  if (is.na(disease_col)) stop("Metadata içinde disease/group kolonu bulunamadı.")
  if (is.na(age_col)) stop("Metadata içinde age kolonu bulunamadı.")
  if (is.na(sex_col)) stop("Metadata içinde sex kolonu bulunamadı.")

  out <- meta_raw %>%
    mutate(
      disease_group_raw = as.character(.data[[disease_col]]),
      disease_group = case_when(
        disease_group_raw %in% c("Healthy control", "Healthy", "healthy", "control", "Control") ~ "Healthy control",
        disease_group_raw %in% c("Cancer patient", "Cancer", "cancer", "patient", "Patient") ~ "Cancer patient",
        str_detect(disease_group_raw, regex("healthy|control", ignore_case = TRUE)) ~ "Healthy control",
        str_detect(disease_group_raw, regex("cancer|patient", ignore_case = TRUE)) ~ "Cancer patient",
        TRUE ~ disease_group_raw
      ),
      disease_group = factor(disease_group, levels = c("Healthy control", "Cancer patient")),
      age_source_value = suppressWarnings(as.numeric(.data[[age_col]])),
      age_for_model = ifelse(
        is.na(age_source_value) | age_source_value < 18 | age_source_value > 100,
        NA_real_,
        age_source_value
      ),
      age_invalid_reason = case_when(
        is.na(age_source_value) ~ "missing",
        age_source_value < 18 ~ "age_below_18",
        age_source_value > 100 ~ "age_above_100",
        TRUE ~ "valid"
      ),
      sex_raw = as.character(.data[[sex_col]]),
      sex = case_when(
        sex_raw %in% c("Female", "F", "female", "f") ~ "Female",
        sex_raw %in% c("Male", "M", "male", "m") ~ "Male",
        TRUE ~ sex_raw
      ),
      sex = factor(sex)
    )

  # Add clinical variables if available.
  optional_cols <- intersect(
    c(
      "cancer_subgroup",
      "therapy_status_4level",
      "chemotherapy",
      "targeted_therapy",
      "any_immunotherapy",
      "ici_immunotherapy",
      "endocrine_hormonal",
      "adc",
      "experimental",
      "radiotherapy",
      "therapy_line_number",
      "time_from_start_days",
      "time_from_treatment_start_days"
    ),
    names(out)
  )

  out <- out %>%
    select(
      subject_id,
      disease_group,
      age_source_value,
      age_for_model,
      age_invalid_reason,
      sex,
      any_of(optional_cols)
    ) %>%
    distinct(subject_id, .keep_all = TRUE)

  out
}

# ------------------------------------------------------------
# 4. Find metadata source
# ------------------------------------------------------------

clinical_matrix_dir <- sd_integrated_dir()
master_metadata_dir <- sd_metadata_dir()

candidate_metadata_inventory <- tibble()

metadata_source_label <- NA_character_
metadata_source_file <- NA_character_
metadata_raw <- NULL

# 4A. Prefer integrated clinical immune matrix CSV.
if (dir.exists(clinical_matrix_dir)) {

  csv_candidates <- list.files(
    clinical_matrix_dir,
    pattern = "\\.csv$",
    full.names = TRUE,
    recursive = TRUE
  )

  csv_candidates <- csv_candidates[
    grepl("integrated|matrix|ALL6|ALL7", basename(csv_candidates), ignore.case = TRUE)
  ]

  inspect_csv <- function(fp) {
    dat0 <- tryCatch(readr::read_csv(fp, show_col_types = FALSE), error = function(e) NULL)

    if (is.null(dat0)) {
      return(tibble(
        file_path = fp,
        file_name = basename(fp),
        readable = FALSE,
        n_rows = NA_integer_,
        n_cols = NA_integer_,
        has_subject_id = FALSE,
        has_disease = FALSE,
        has_age = FALSE,
        has_sex = FALSE,
        has_clinical = FALSE,
        n_score_cols = NA_integer_,
        priority = -Inf
      ))
    }

    nm <- names(dat0)

    has_disease <- any(c("disease_group_model", "disease_group_clinical", "disease_group", "group", "condition") %in% nm)
    has_age <- any(c("age_for_model", "age_clinical", "age_years", "age") %in% nm)
    has_sex <- any(c("sex_for_clinical_model", "sex_clinical", "sex", "gender") %in% nm)
    has_clinical <- any(c("cancer_subgroup", "therapy_status_4level", "chemotherapy", "targeted_therapy") %in% nm)
    score_cols <- grep("^CP[0-9]+_.*_score$", nm, value = TRUE)

    priority <- 0
    priority <- priority + ifelse("subject_id" %in% nm, 1000, 0)
    priority <- priority + ifelse(nrow(dat0) >= 800 & nrow(dat0) <= 900, 300, 0)
    priority <- priority + ifelse(has_disease, 100, 0)
    priority <- priority + ifelse(has_age, 100, 0)
    priority <- priority + ifelse(has_sex, 100, 0)
    priority <- priority + ifelse(has_clinical, 100, 0)
    priority <- priority + length(score_cols)
    priority <- priority + ifelse(grepl("ALL7", basename(fp), ignore.case = TRUE), 80, 0)
    priority <- priority + ifelse(grepl("ALL6", basename(fp), ignore.case = TRUE), 50, 0)
    priority <- priority - ifelse(grepl("summary|dictionary|availability|counts|results|variable_sources|exposure|therapy|subgroup", basename(fp), ignore.case = TRUE), 500, 0)

    tibble(
      file_path = fp,
      file_name = basename(fp),
      readable = TRUE,
      n_rows = nrow(dat0),
      n_cols = ncol(dat0),
      has_subject_id = "subject_id" %in% nm,
      has_disease = has_disease,
      has_age = has_age,
      has_sex = has_sex,
      has_clinical = has_clinical,
      n_score_cols = length(score_cols),
      priority = priority
    )
  }

  if (length(csv_candidates) > 0) {
    candidate_metadata_inventory <- bind_rows(lapply(csv_candidates, inspect_csv)) %>%
      arrange(desc(priority), desc(n_score_cols), desc(n_cols))

    valid_csv <- candidate_metadata_inventory %>%
      filter(
        readable == TRUE,
        has_subject_id == TRUE,
        has_disease == TRUE,
        has_age == TRUE,
        has_sex == TRUE,
        n_rows >= 800,
        n_rows <= 900
      )

    if (nrow(valid_csv) > 0) {
      metadata_source_file <- valid_csv$file_path[1]
      metadata_source_label <- "integrated_clinical_immune_matrix_csv"
      metadata_raw <- readr::read_csv(metadata_source_file, show_col_types = FALSE)
    }
  }
}

# 4B. Fallback: CP7 master metadata RData
if (is.null(metadata_raw) && dir.exists(master_metadata_dir)) {

  rdata_candidates <- list.files(
    master_metadata_dir,
    pattern = "\\.RData$|\\.rda$",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )

  # Prefer files with metadata/merge/age_QC/cp7 in their name.
  rdata_candidates <- rdata_candidates[
    grepl("metadata|merge|age|CP7|cp7", basename(rdata_candidates), ignore.case = TRUE)
  ]

  for (fp in rdata_candidates) {
    env <- new.env()
    try_load <- tryCatch({
      load(fp, envir = env)
      TRUE
    }, error = function(e) FALSE)

    if (!try_load) next

    obj_names <- ls(env)

    for (oo in obj_names) {
      x <- get(oo, envir = env)
      if (is.data.frame(x) && "subject_id" %in% names(x)) {
        nm <- names(x)
        has_disease <- any(c("disease_group", "disease_group_clinical", "group", "condition") %in% nm)
        has_age <- any(c("age_for_model", "age", "age_years") %in% nm)
        has_sex <- any(c("sex", "sex_clinical", "gender") %in% nm)
        if (has_disease && has_age && has_sex && nrow(x) >= 800) {
          metadata_raw <- x
          metadata_source_file <- paste0(fp, "::", oo)
          metadata_source_label <- "master_metadata_RData"
          break
        }
      }
    }
    if (!is.null(metadata_raw)) break
  }
}

if (is.null(metadata_raw)) {
  if (nrow(candidate_metadata_inventory) > 0) {
    readr::write_csv(
      candidate_metadata_inventory,
      file.path(out_dir, "SDY2583_CP10_metadata_candidate_inventory_STEP3A.csv")
    )
    cat("\nMetadata candidate inventory:\n")
    print(as.data.frame(candidate_metadata_inventory), row.names = FALSE)
  }
  stop("Metadata source bulunamadı. Integrated matrix veya master metadata klasörünü kontrol et.")
}

metadata_std <- standardize_metadata(metadata_raw)

# ------------------------------------------------------------
# 5. Merge CP10 features with metadata
# ------------------------------------------------------------

cp10_analysis_data <- cp10_features %>%
  left_join(metadata_std, by = "subject_id") %>%
  mutate(
    model_ready_age_sex = feature_ok == TRUE &
      !is.na(disease_group) &
      !is.na(age_for_model) &
      !is.na(sex),
    sex_binary = ifelse(as.character(sex) %in% c("Female", "Male"), as.character(sex), NA_character_),
    sex_binary = factor(sex_binary, levels = c("Female", "Male")),
    cd45_event_qc_bin = case_when(
      is.na(n_cd45_viable) ~ "missing",
      n_cd45_viable < 50 ~ "<50",
      n_cd45_viable < 100 ~ "50-99",
      n_cd45_viable < 300 ~ "100-299",
      n_cd45_viable < 1000 ~ "300-999",
      TRUE ~ ">=1000"
    )
  )

# ------------------------------------------------------------
# 6. Summaries
# ------------------------------------------------------------

metadata_source_summary <- tibble(
  metadata_source_label = metadata_source_label,
  metadata_source_file = metadata_source_file,
  n_metadata_rows = nrow(metadata_raw),
  n_metadata_subjects = n_distinct(metadata_raw$subject_id),
  n_standardized_metadata_subjects = n_distinct(metadata_std$subject_id)
)

merge_summary <- tibble(
  n_cp10_rows = nrow(cp10_analysis_data),
  n_unique_subjects = n_distinct(cp10_analysis_data$subject_id),
  n_feature_ok = sum(cp10_analysis_data$feature_ok == TRUE, na.rm = TRUE),
  n_with_disease_group = sum(!is.na(cp10_analysis_data$disease_group)),
  n_with_valid_age_for_model = sum(!is.na(cp10_analysis_data$age_for_model)),
  n_invalid_or_missing_age_for_model = sum(is.na(cp10_analysis_data$age_for_model)),
  n_with_sex = sum(!is.na(cp10_analysis_data$sex)),
  n_with_binary_sex = sum(!is.na(cp10_analysis_data$sex_binary)),
  n_model_ready_age_sex = sum(cp10_analysis_data$model_ready_age_sex, na.rm = TRUE),
  n_channel_order_mismatch_files = sum(cp10_analysis_data$channel_order_mismatch_file == TRUE, na.rm = TRUE),
  n_marker_mismatch_files = sum(cp10_analysis_data$marker_mismatch_file == TRUE, na.rm = TRUE)
)

disease_counts <- cp10_analysis_data %>%
  count(disease_group, name = "n") %>%
  arrange(disease_group)

model_ready_by_disease <- cp10_analysis_data %>%
  count(disease_group, model_ready_age_sex, name = "n") %>%
  arrange(disease_group, model_ready_age_sex)

age_summary_by_disease <- cp10_analysis_data %>%
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

sex_counts_by_disease <- cp10_analysis_data %>%
  count(disease_group, sex, name = "n") %>%
  arrange(disease_group, sex)

cd45_event_qc_by_disease <- cp10_analysis_data %>%
  count(disease_group, cd45_event_qc_bin, name = "n") %>%
  arrange(disease_group, cd45_event_qc_bin)

technical_mismatch_by_disease <- cp10_analysis_data %>%
  count(disease_group, channel_order_mismatch_file, marker_mismatch_file, name = "n") %>%
  arrange(desc(channel_order_mismatch_file), desc(marker_mismatch_file), disease_group)

selected_feature_medians_by_disease <- cp10_analysis_data %>%
  group_by(disease_group) %>%
  summarise(
    n = n(),
    median_pct_t_like_within_cd45 = median(pct_t_like_within_cd45, na.rm = TRUE),
    median_pct_b_like_within_cd45 = median(pct_b_like_within_cd45, na.rm = TRUE),
    median_pct_nk_like_within_cd45 = median(pct_nk_like_within_cd45, na.rm = TRUE),
    median_pct_cd14_mono_like_within_cd45 = median(pct_cd14_mono_like_within_cd45, na.rm = TRUE),
    median_pct_hladr_low_within_cd14_mono_like = median(pct_hladr_low_within_cd14_mono_like, na.rm = TRUE),
    median_pct_cd66b_gran_like_within_cd45 = median(pct_cd66b_gran_like_within_cd45, na.rm = TRUE),
    median_pct_ccr3_cd66b_eosinophil_like_within_cd45 = median(pct_ccr3_cd66b_eosinophil_like_within_cd45, na.rm = TRUE),
    median_pct_cd123_hladr_dc_like_within_cd45 = median(pct_cd123_hladr_dc_like_within_cd45, na.rm = TRUE),
    median_ratio_myeloid_granulocytic_to_lymphoid_like = median(ratio_myeloid_granulocytic_to_lymphoid_like, na.rm = TRUE),
    .groups = "drop"
  )

# Optional clinical variable availability
clinical_availability <- tibble(
  variable = c(
    "cancer_subgroup",
    "therapy_status_4level",
    "chemotherapy",
    "targeted_therapy",
    "any_immunotherapy",
    "ici_immunotherapy",
    "therapy_line_number",
    "time_from_start_days"
  )
) %>%
  mutate(
    present = variable %in% names(cp10_analysis_data),
    n_nonmissing = purrr::map_int(variable, function(v) {
      if (v %in% names(cp10_analysis_data)) sum(!is.na(cp10_analysis_data[[v]])) else 0L
    })
  )

# ------------------------------------------------------------
# 7. Save outputs
# ------------------------------------------------------------

readr::write_csv(metadata_source_summary, file.path(out_dir, "SDY2583_CP10_metadata_source_summary_STEP3A.csv"))
readr::write_csv(merge_summary, file.path(out_dir, "SDY2583_CP10_metadata_merge_summary_STEP3A.csv"))
readr::write_csv(cp10_analysis_data, file.path(out_dir, "SDY2583_CP10_analysis_data_STEP3A.csv"))
readr::write_csv(disease_counts, file.path(out_dir, "SDY2583_CP10_disease_counts_STEP3A.csv"))
readr::write_csv(model_ready_by_disease, file.path(out_dir, "SDY2583_CP10_model_ready_by_disease_STEP3A.csv"))
readr::write_csv(age_summary_by_disease, file.path(out_dir, "SDY2583_CP10_age_summary_by_disease_STEP3A.csv"))
readr::write_csv(sex_counts_by_disease, file.path(out_dir, "SDY2583_CP10_sex_counts_by_disease_STEP3A.csv"))
readr::write_csv(cd45_event_qc_by_disease, file.path(out_dir, "SDY2583_CP10_CD45_event_QC_by_disease_STEP3A.csv"))
readr::write_csv(technical_mismatch_by_disease, file.path(out_dir, "SDY2583_CP10_technical_mismatch_by_disease_STEP3A.csv"))
readr::write_csv(selected_feature_medians_by_disease, file.path(out_dir, "SDY2583_CP10_selected_feature_medians_by_disease_STEP3A.csv"))
readr::write_csv(clinical_availability, file.path(out_dir, "SDY2583_CP10_clinical_variable_availability_STEP3A.csv"))

if (nrow(candidate_metadata_inventory) > 0) {
  readr::write_csv(candidate_metadata_inventory, file.path(out_dir, "SDY2583_CP10_metadata_candidate_inventory_STEP3A.csv"))
}

save(
  cp10_analysis_data,
  cp10_features,
  metadata_raw,
  metadata_std,
  metadata_source_summary,
  merge_summary,
  disease_counts,
  model_ready_by_disease,
  age_summary_by_disease,
  sex_counts_by_disease,
  cd45_event_qc_by_disease,
  technical_mismatch_by_disease,
  selected_feature_medians_by_disease,
  clinical_availability,
  feature_dictionary,
  thresholds,
  file = file.path(rdata_dir, "SDY2583_CP10_STEP3A_metadata_merge_age_QC.RData")
)

# ------------------------------------------------------------
# 8. Console output
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SDY2583 CP10 STEP 3A COMPLETE: METADATA MERGE / AGE QC\n")
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

cat("\nCD45 event QC by disease:\n")
print(as.data.frame(cd45_event_qc_by_disease), row.names = FALSE)

cat("\nTechnical mismatch by disease:\n")
print(as.data.frame(technical_mismatch_by_disease), row.names = FALSE)

cat("\nSelected feature medians by disease, unadjusted descriptive only:\n")
print(as.data.frame(selected_feature_medians_by_disease), row.names = FALSE)

cat("\nClinical variable availability:\n")
print(as.data.frame(clinical_availability), row.names = FALSE)

cat("\nOutputs saved in:\n")
print(out_dir)

cat("============================================================\n")
