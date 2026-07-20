# Reconstruct a portable subject/FCS metadata matrix directly from the local
# ImmPort SDY2583 tabular download. Participant-level outputs are written only
# under data/derived and outputs/, both excluded from Git.
rm(list = ls())
source(file.path(Sys.getenv("SDY2583_REPO_ROOT", unset = "."), "R", "shared", "bootstrap.R"))

for (p in c("dplyr", "readr", "stringr", "tibble", "purrr", "tidyr")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(stringr)
  library(tibble); library(purrr); library(tidyr)
})

root <- sd_immport_download_dir()
out_dir <- sd_metadata_dir()
if (!dir.exists(root)) {
  stop("ImmPort download directory does not exist: ", root,
       "\nSet SDY2583_IMMPORT_DOWNLOAD_DIR to the unpacked SDY2583 tabular package.")
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
all_files <- list.files(root, full.names = TRUE, recursive = TRUE, ignore.case = TRUE)

find_one <- function(patterns, required = TRUE) {
  hits <- all_files[vapply(
    all_files,
    function(x) any(grepl(patterns, basename(x), ignore.case = TRUE)),
    logical(1)
  )]
  hits <- hits[!is.na(file.info(hits)$isdir) & !file.info(hits)$isdir]
  if (length(hits) == 0L) {
    if (required) {
      stop("Required ImmPort tabular file not found for patterns: ",
           paste(patterns, collapse = ", "))
    }
    return(NA_character_)
  }
  hits[order(nchar(hits), hits)][1]
}

read_tab <- function(path) {
  if (is.na(path)) return(NULL)
  tryCatch(
    readr::read_tsv(path, show_col_types = FALSE, progress = FALSE, guess_max = 100000),
    error = function(e) read.delim(
      path, sep = "\t", header = TRUE, check.names = FALSE,
      stringsAsFactors = FALSE
    )
  )
}

first_col <- function(df, candidates, regex = NULL) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) > 0L) return(hit[1])
  if (!is.null(regex)) {
    hit <- names(df)[grepl(regex, names(df), ignore.case = TRUE)]
    if (length(hit) > 0L) return(hit[1])
  }
  NA_character_
}

first_non_missing <- function(x) {
  y <- x[!is.na(x)]
  if (is.character(y)) y <- y[nzchar(trimws(y))]
  if (length(y) == 0L) return(x[NA_integer_][1])
  y[1]
}

rename_key <- function(df, old_name, new_name) {
  if (is.na(old_name) || !old_name %in% names(df)) return(df)
  names(df)[names(df) == old_name] <- new_name
  df
}

arm_file <- find_one(c("^arm_2_subject\\.txt$", "arm.*2.*subject"))
subject_file <- find_one(c("^subject\\.txt$", "^subjects\\.txt$"), required = FALSE)
flow_file <- find_one(c(
  "Subject_2_Flow_cytometry_result",
  "subject.*flow.*cytometry.*result",
  "flow_cytometry_result"
))

arm <- read_tab(arm_file)
subject <- read_tab(subject_file)
flow <- read_tab(flow_file)

arm_subject_col <- first_col(
  arm, c("SUBJECT_ACCESSION", "Subject Accession", "subject_accession"),
  "subject.*accession"
)
arm_name_col <- first_col(arm, c("ARM_NAME", "ARM Name", "arm_name"), "arm.*name")
arm_age_col <- first_col(arm, c("AGE", "Age", "Subject Age", "age"), "age")
arm_sex_col <- first_col(
  arm, c("GENDER", "Gender", "SEX", "Sex", "gender", "sex"),
  "gender|sex"
)
if (is.na(arm_subject_col) || is.na(arm_name_col)) {
  stop("arm_2_subject lacks subject accession or arm name.")
}

subject_core <- arm %>%
  transmute(
    subject_accession = as.character(.data[[arm_subject_col]]),
    arm_name = as.character(.data[[arm_name_col]]),
    age_from_arm = if (is.na(arm_age_col)) NA_real_ else
      suppressWarnings(as.numeric(.data[[arm_age_col]])),
    sex_from_arm = if (is.na(arm_sex_col)) NA_character_ else
      as.character(.data[[arm_sex_col]])
  )

if (!is.null(subject)) {
  s_subject_col <- first_col(
    subject, c("SUBJECT_ACCESSION", "Subject Accession", "subject_accession"),
    "subject.*accession"
  )
  s_age_col <- first_col(subject, c("AGE", "Age", "age"), "age")
  s_sex_col <- first_col(
    subject, c("GENDER", "Gender", "SEX", "Sex", "gender", "sex"),
    "gender|sex"
  )
  if (!is.na(s_subject_col)) {
    subject_from_subject_file <- subject %>%
      transmute(
        subject_accession = as.character(.data[[s_subject_col]]),
        age_from_subject = if (is.na(s_age_col)) NA_real_ else
          suppressWarnings(as.numeric(.data[[s_age_col]])),
        sex_from_subject = if (is.na(s_sex_col)) NA_character_ else
          as.character(.data[[s_sex_col]])
      )
    subject_core <- subject_core %>%
      left_join(subject_from_subject_file, by = "subject_accession")
  }
}

if (!"age_from_subject" %in% names(subject_core)) subject_core$age_from_subject <- NA_real_
if (!"sex_from_subject" %in% names(subject_core)) subject_core$sex_from_subject <- NA_character_

standardize_sex <- function(x) {
  y <- trimws(tolower(as.character(x)))
  case_when(
    y %in% c("f", "female", "woman") ~ "Female",
    y %in% c("m", "male", "man") ~ "Male",
    is.na(y) | y == "" ~ NA_character_,
    TRUE ~ stringr::str_to_title(y)
  )
}

subject_core <- subject_core %>%
  mutate(
    age_raw = coalesce(age_from_subject, age_from_arm),
    age_for_model = ifelse(age_raw >= 18 & age_raw <= 100, age_raw, NA_real_),
    sex = standardize_sex(coalesce(sex_from_subject, sex_from_arm)),
    disease_group = case_when(
      str_detect(arm_name, regex("healthy|control", ignore_case = TRUE)) ~ "Healthy control",
      str_detect(arm_name, regex("cancer|patient|tumou?r", ignore_case = TRUE)) ~ "Cancer patient",
      TRUE ~ arm_name
    ),
    age_group_for_model = factor(
      case_when(
        is.na(age_for_model) ~ NA_character_,
        age_for_model < 40 ~ "Young_<40",
        age_for_model < 60 ~ "Middle_40_59",
        TRUE ~ "Older_60plus"
      ),
      levels = c("Young_<40", "Middle_40_59", "Older_60plus")
    )
  ) %>%
  group_by(subject_accession) %>%
  summarise(across(everything(), first_non_missing), .groups = "drop")

f_subject_col <- first_col(
  flow, c("SUBJECT_ACCESSION", "Subject Accession", "subject_accession"),
  "subject.*accession"
)
f_file_col <- first_col(
  flow,
  c("FILE_NAME", "File Name", "file_name", "RESULT_FILE_NAME", "Result File Name"),
  "file.*name"
)
f_panel_col <- first_col(flow, c("PANEL", "Panel", "panel"), "panel")
if (is.na(f_subject_col) || is.na(f_file_col)) {
  stop("Flow result mapping lacks subject accession or file name.")
}

flow_map <- flow %>%
  transmute(
    subject_accession = as.character(.data[[f_subject_col]]),
    result_file_name = basename(as.character(.data[[f_file_col]])),
    panel_from_table = if (is.na(f_panel_col)) NA_character_ else
      as.character(.data[[f_panel_col]])
  ) %>%
  mutate(
    subject_id = str_extract(result_file_name, "DBG[0-9]+"),
    panel = coalesce(str_extract(result_file_name, "CP[0-9]+"), panel_from_table)
  ) %>%
  filter(!is.na(subject_id), !is.na(panel)) %>%
  distinct(subject_id, panel, .keep_all = TRUE) %>%
  select(-panel_from_table)

metadata_matrix <- flow_map %>%
  left_join(subject_core, by = "subject_accession") %>%
  mutate(
    subject_id = as.character(subject_id),
    disease_group = factor(
      disease_group, levels = c("Healthy control", "Cancer patient")
    ),
    sex = factor(sex),
    sex_binary = factor(
      ifelse(as.character(sex) %in% c("Female", "Male"), as.character(sex), NA_character_),
      levels = c("Female", "Male")
    )
  )

# Optional local cancer-subgroup/therapy annotation matrix.
clinical_file <- path.expand(Sys.getenv("SDY2583_CLINICAL_ANNOTATION_FILE", unset = ""))
if (nzchar(clinical_file) && file.exists(clinical_file)) {
  clinical <- readr::read_csv(clinical_file, show_col_types = FALSE, progress = FALSE)
  key <- first_col(
    clinical,
    c("subject_id", "subject_accession", "SUBJECT_ACCESSION"),
    "subject"
  )
  if (!is.na(key)) {
    if (grepl("accession", key, ignore.case = TRUE)) {
      clinical <- rename_key(clinical, key, "subject_accession")
      clinical$subject_accession <- as.character(clinical$subject_accession)
      metadata_matrix <- metadata_matrix %>%
        left_join(clinical, by = "subject_accession", suffix = c("", "_clinical"))
    } else {
      clinical <- rename_key(clinical, key, "subject_id")
      clinical$subject_id <- as.character(clinical$subject_id)
      metadata_matrix <- metadata_matrix %>%
        left_join(clinical, by = "subject_id", suffix = c("", "_clinical"))
    }
  }
}

summary <- tibble(
  n_rows = nrow(metadata_matrix),
  n_subjects = n_distinct(metadata_matrix$subject_id),
  n_panels = n_distinct(metadata_matrix$panel),
  n_disease_mapped = sum(!is.na(metadata_matrix$disease_group)),
  n_valid_age = sum(!is.na(metadata_matrix$age_for_model)),
  n_sex_mapped = sum(!is.na(metadata_matrix$sex)),
  arm_file = arm_file,
  subject_file = subject_file,
  flow_file = flow_file
)

readr::write_csv(
  metadata_matrix,
  file.path(out_dir, "SDY2583_subject_FCS_metadata_matrix_RECONSTRUCTED.csv")
)
readr::write_csv(
  subject_core,
  file.path(out_dir, "SDY2583_subject_core_metadata_RECONSTRUCTED.csv")
)
readr::write_csv(
  flow_map,
  file.path(out_dir, "SDY2583_subject_to_FCS_panel_map_RECONSTRUCTED.csv")
)
readr::write_csv(
  summary,
  file.path(out_dir, "SDY2583_metadata_integration_summary_RECONSTRUCTED.csv")
)
cat(
  "Set SDY2583_METADATA_MATRIX_FILE=",
  file.path(out_dir, "SDY2583_subject_FCS_metadata_matrix_RECONSTRUCTED.csv"),
  "\n", sep = ""
)
print(summary)
