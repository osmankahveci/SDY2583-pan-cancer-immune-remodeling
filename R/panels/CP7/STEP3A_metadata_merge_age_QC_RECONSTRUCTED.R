# ============================================================
# SDY2583 CP7
# STEP 3A RECONSTRUCTED: metadata merge and age quality control
#
# Reconstruction basis:
#   - archived CP7 merged analysis-table schema;
#   - archived CP7 Methods/Results record;
#   - metadata conventions in the recovered panel pipelines.
#
# This is reconstructed source, not the original archived script.
# It expects a local subject-level metadata/integrated matrix containing a
# DBG-style subject_id plus disease, age, and sex fields. Participant-level
# metadata are never written to the repository.
# ============================================================

rm(list = ls())
source(file.path(Sys.getenv("SDY2583_REPO_ROOT", unset = "."), "R", "shared", "bootstrap.R"))

cran_pkgs <- c("dplyr", "readr", "stringr", "tibble", "purrr")
for (p in cran_pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(purrr)
})

analysis_dir <- sd_analysis_dir("CP7")
rdata_dir <- file.path(analysis_dir, "11_RData")
out_dir <- file.path(analysis_dir, "04_metadata_merge")
age_qc_dir <- file.path(analysis_dir, "05_age_QC")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(age_qc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rdata_dir, recursive = TRUE, showWarnings = FALSE)

step2_file <- file.path(rdata_dir, "SDY2583_CP7_STEP2_feature_extraction_RECONSTRUCTED.RData")
if (!file.exists(step2_file)) stop("Run reconstructed CP7 Step 2 first: ", step2_file)
load(step2_file)
if (!exists("feature_table")) stop("Step 2 RData does not contain feature_table.")

first_existing_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0L) NA_character_ else hit[[1L]]
}

standardize_disease <- function(x) {
  y <- trimws(as.character(x))
  dplyr::case_when(
    stringr::str_detect(y, stringr::regex("healthy|control", ignore_case = TRUE)) ~ "Healthy control",
    stringr::str_detect(y, stringr::regex("cancer|patient|tumou?r", ignore_case = TRUE)) ~ "Cancer patient",
    TRUE ~ y
  )
}

standardize_sex <- function(x) {
  y <- trimws(as.character(x))
  dplyr::case_when(
    tolower(y) %in% c("female", "f", "woman") ~ "Female",
    tolower(y) %in% c("male", "m", "man") ~ "Male",
    is.na(y) | y == "" | tolower(y) %in% c("na", "unknown") ~ NA_character_,
    TRUE ~ y
  )
}

# An exact metadata matrix can be supplied explicitly. Otherwise, scan the
# configured integrated and metadata directories and rank compatible CSVs.
explicit_metadata_file <- path.expand(Sys.getenv("SDY2583_METADATA_MATRIX_FILE", unset = ""))
search_dirs <- unique(c(sd_integrated_dir(), sd_metadata_dir()))
search_dirs <- search_dirs[dir.exists(search_dirs)]

candidate_files <- character()
if (nzchar(explicit_metadata_file) && file.exists(explicit_metadata_file)) {
  candidate_files <- explicit_metadata_file
} else if (length(search_dirs) > 0L) {
  candidate_files <- unique(unlist(lapply(search_dirs, function(dd) {
    list.files(dd, pattern = "\\.csv$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  })))
}

if (length(candidate_files) == 0L) {
  stop(
    "No metadata CSV was found. Set SDY2583_METADATA_MATRIX_FILE to a local ",
    "subject-level matrix with subject_id, disease, age, and sex columns."
  )
}

inspect_candidate <- function(path) {
  header <- tryCatch(
    readr::read_csv(path, n_max = 5, show_col_types = FALSE, progress = FALSE),
    error = function(e) NULL
  )
  if (is.null(header)) {
    return(tibble(file_path = path, readable = FALSE, score = -Inf,
                  subject_col = NA_character_, disease_col = NA_character_,
                  age_col = NA_character_, sex_col = NA_character_))
  }

  subject_col <- first_existing_col(header, c("subject_id", "dbg_id", "subject"))
  disease_col <- first_existing_col(header, c(
    "disease_group", "disease_group_model", "disease_group_clinical",
    "group", "condition", "arm_name"
  ))
  age_col <- first_existing_col(header, c(
    "age_for_model", "age_raw", "age_clinical", "age_years", "age"
  ))
  sex_col <- first_existing_col(header, c(
    "sex", "sex_for_clinical_model", "sex_clinical", "gender"
  ))

  score <- 0
  score <- score + ifelse(!is.na(subject_col), 1000, 0)
  score <- score + ifelse(!is.na(disease_col), 300, 0)
  score <- score + ifelse(!is.na(age_col), 300, 0)
  score <- score + ifelse(!is.na(sex_col), 300, 0)
  score <- score + ifelse(grepl("integrated|clinical|metadata|matrix", basename(path), ignore.case = TRUE), 100, 0)
  score <- score - ifelse(grepl("summary|statistics|result|figure|dictionary|manifest", basename(path), ignore.case = TRUE), 500, 0)

  tibble(
    file_path = path,
    readable = TRUE,
    score = score,
    subject_col = subject_col,
    disease_col = disease_col,
    age_col = age_col,
    sex_col = sex_col
  )
}

candidate_inventory <- purrr::map_dfr(candidate_files, inspect_candidate) %>%
  arrange(desc(score), file_path)
readr::write_csv(candidate_inventory, file.path(out_dir, "SDY2583_CP7_metadata_candidate_inventory_RECONSTRUCTED.csv"))

valid_candidates <- candidate_inventory %>%
  filter(readable, !is.na(subject_col), !is.na(disease_col), !is.na(age_col), !is.na(sex_col))
if (nrow(valid_candidates) == 0L) {
  stop("No candidate metadata file contained subject_id, disease, age, and sex fields.")
}

metadata_source_file <- valid_candidates$file_path[[1L]]
metadata_raw <- readr::read_csv(metadata_source_file, show_col_types = FALSE, progress = FALSE)
subject_col <- valid_candidates$subject_col[[1L]]
disease_col <- valid_candidates$disease_col[[1L]]
age_col <- valid_candidates$age_col[[1L]]
sex_col <- valid_candidates$sex_col[[1L]]

metadata_standardized <- metadata_raw %>%
  transmute(
    subject_id = as.character(.data[[subject_col]]),
    disease_group = standardize_disease(.data[[disease_col]]),
    age_raw = suppressWarnings(as.numeric(.data[[age_col]])),
    sex = standardize_sex(.data[[sex_col]])
  ) %>%
  filter(!is.na(subject_id), subject_id != "") %>%
  group_by(subject_id) %>%
  summarise(
    disease_group = dplyr::first(stats::na.omit(disease_group), default = NA_character_),
    age_raw = dplyr::first(stats::na.omit(age_raw), default = NA_real_),
    sex = dplyr::first(stats::na.omit(sex), default = NA_character_),
    .groups = "drop"
  ) %>%
  mutate(
    disease_group = factor(disease_group, levels = c("Healthy control", "Cancer patient")),
    age_for_model = ifelse(age_raw >= 18 & age_raw <= 100, age_raw, NA_real_),
    age_group_for_model = case_when(
      is.na(age_for_model) ~ NA_character_,
      age_for_model < 40 ~ "Young_<40",
      age_for_model < 60 ~ "Middle_40_59",
      TRUE ~ "Older_60plus"
    ),
    age_group_for_model = factor(
      age_group_for_model,
      levels = c("Young_<40", "Middle_40_59", "Older_60plus")
    ),
    sex = factor(sex),
    sex_binary = factor(
      ifelse(as.character(sex) %in% c("Female", "Male"), as.character(sex), NA_character_),
      levels = c("Female", "Male")
    )
  )

analysis_data <- feature_table %>%
  mutate(subject_id = as.character(subject_id)) %>%
  left_join(metadata_standardized, by = "subject_id")

feature_columns <- names(analysis_data)[
  grepl("^(pct_|median_)", names(analysis_data)) & vapply(analysis_data, is.numeric, logical(1))
]
feature_dictionary <- tibble(feature = feature_columns)

metadata_merge_summary <- tibble(
  metadata_source_file = normalizePath(metadata_source_file, mustWork = FALSE),
  n_feature_rows = nrow(feature_table),
  n_unique_feature_subjects = n_distinct(feature_table$subject_id),
  n_metadata_subjects = nrow(metadata_standardized),
  n_disease_mapped = sum(!is.na(analysis_data$disease_group)),
  n_age_raw_mapped = sum(!is.na(analysis_data$age_raw)),
  n_valid_age = sum(!is.na(analysis_data$age_for_model)),
  n_invalid_or_missing_age = sum(is.na(analysis_data$age_for_model)),
  n_sex_mapped = sum(!is.na(analysis_data$sex)),
  n_healthy = sum(analysis_data$disease_group == "Healthy control", na.rm = TRUE),
  n_cancer = sum(analysis_data$disease_group == "Cancer patient", na.rm = TRUE)
)

age_qc_summary <- tibble(
  metric = c(
    "n_total", "n_valid_age", "n_invalid_or_missing_age", "minimum_valid_age",
    "q1_valid_age", "median_valid_age", "mean_valid_age", "q3_valid_age", "maximum_valid_age"
  ),
  value = c(
    nrow(analysis_data),
    sum(!is.na(analysis_data$age_for_model)),
    sum(is.na(analysis_data$age_for_model)),
    min(analysis_data$age_for_model, na.rm = TRUE),
    as.numeric(quantile(analysis_data$age_for_model, 0.25, na.rm = TRUE)),
    median(analysis_data$age_for_model, na.rm = TRUE),
    mean(analysis_data$age_for_model, na.rm = TRUE),
    as.numeric(quantile(analysis_data$age_for_model, 0.75, na.rm = TRUE)),
    max(analysis_data$age_for_model, na.rm = TRUE)
  )
)

readr::write_csv(metadata_standardized, file.path(out_dir, "SDY2583_CP7_subject_metadata_final_RECONSTRUCTED.csv"))
readr::write_csv(analysis_data, file.path(out_dir, "SDY2583_CP7_FULL_850_analysis_data_with_metadata_RECONSTRUCTED.csv"))
readr::write_csv(metadata_merge_summary, file.path(out_dir, "SDY2583_CP7_metadata_merge_summary_RECONSTRUCTED.csv"))
readr::write_csv(age_qc_summary, file.path(age_qc_dir, "SDY2583_CP7_age_QC_summary_RECONSTRUCTED.csv"))
readr::write_csv(feature_dictionary, file.path(out_dir, "SDY2583_CP7_feature_dictionary_RECONSTRUCTED.csv"))

save(
  analysis_data, metadata_standardized, metadata_merge_summary,
  age_qc_summary, feature_dictionary, metadata_source_file,
  file = file.path(rdata_dir, "SDY2583_CP7_STEP3A_metadata_merge_age_QC_RECONSTRUCTED.RData")
)

print(metadata_merge_summary)
message("Archived CP7 benchmark: 850 mapped subjects; 832 valid ages; 408 healthy and 442 cancer overall.")
