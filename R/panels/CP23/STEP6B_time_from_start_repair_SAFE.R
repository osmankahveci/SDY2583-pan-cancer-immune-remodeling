source(file.path(Sys.getenv("SDY2583_REPO_ROOT", unset = "."), "R", "shared", "bootstrap.R"))

# ============================================================
# SDY2583 CP23
# STEP 6B SAFE: Repair time-from-treatment-start clinical model
#
# Why this exists:
#   CP23 Step 6 completed correctly, but the time-from-start table
#   printed as 0 rows because an empty extra_cols bind created an
#   empty tibble in the helper function. This script recomputes only
#   the CP23 time-from-start models and overwrites the standard CSV.
#
# Input:
#   data/derived/clinical_integration/
#     SDY2583_integrated_clinical_immune_score_matrix_ALL10_with_CP23.csv
#
# Output:
#   08_CP23_clinical_annotation/
#     SDY2583_CP23_time_from_start_results_STEP6.csv
#     SDY2583_CP23_time_from_start_results_STEP6_REPAIRED.csv
# ============================================================

cran_pkgs <- c("dplyr", "readr", "tibble", "stringr", "broom")
for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(stringr)
  library(broom)
})

filter <- dplyr::filter
select <- dplyr::select
mutate <- dplyr::mutate
arrange <- dplyr::arrange
case_when <- dplyr::case_when
left_join <- dplyr::left_join
distinct <- dplyr::distinct

integrated_dir <- sd_integrated_dir()
out_dir <- file.path(integrated_dir, "08_CP23_clinical_annotation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

all10_file <- file.path(
  integrated_dir,
  "SDY2583_integrated_clinical_immune_score_matrix_ALL10_with_CP23.csv"
)

if (!file.exists(all10_file)) {
  stop("ALL10 CP23 integrated matrix not found: ", all10_file)
}

dat0 <- readr::read_csv(all10_file, show_col_types = FALSE)

pick_col <- function(df, candidates, required = FALSE, label = "column") {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) > 0) return(hit[1])
  if (required) {
    stop("Required ", label, " not found. Tried: ", paste(candidates, collapse = ", "))
  }
  return(NA_character_)
}

to_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

clean_chr <- function(x) {
  y <- as.character(x)
  y[y %in% c("", "NA", "NaN", "NULL", "null", "None")] <- NA_character_
  y
}

disease_col <- pick_col(dat0, c("disease_group_model", "disease_group", "disease", "group"), TRUE, "disease group")
age_col <- pick_col(dat0, c("age_clinical", "age_for_model", "age", "age_raw"), TRUE, "age")
sex_col <- pick_col(dat0, c("sex_for_clinical_model", "sex_for_model", "sex"), TRUE, "sex")
time_col <- pick_col(dat0, c("time_from_start_days_model", "time_from_start_days", "time_from_treatment_start_days"), TRUE, "time from treatment start")

score_cols <- grep("^CP23_.*_score$", names(dat0), value = TRUE)
if (length(score_cols) == 0) stop("No CP23 score columns found in ALL10 matrix.")

score_labels <- tibble(
  score = score_cols,
  score_label = dplyr::case_when(
    score == "CP23_integrated_monocyte_macrophage_like_myeloid_remodeling_score" ~
      "Integrated CP23 monocyte/macrophage-like myeloid remodeling",
    score == "CP23_CD33_HLA_DR_myeloid_repatterning_score" ~
      "CD33/HLA-DR myeloid repatterning",
    score == "CP23_CD9_CD84_activation_attenuation_score" ~
      "CD9/CD84 activation attenuation",
    score == "CP23_CD45_dump_low_myeloid_enrichment_score" ~
      "CD45+ dump-low myeloid enrichment",
    score == "CP23_FcERI_myeloid_attenuation_score" ~
      "FcERI-associated myeloid attenuation",
    TRUE ~ score
  )
)

clinical_df <- dat0 %>%
  mutate(
    disease_group_for_clinical = clean_chr(.data[[disease_col]]),
    age_clinical_model = to_num(.data[[age_col]]),
    sex_clinical_model = clean_chr(.data[[sex_col]]),
    time_from_start_days_model = to_num(.data[[time_col]])
  ) %>%
  mutate(
    disease_group_for_clinical = case_when(
      stringr::str_detect(tolower(disease_group_for_clinical), "cancer|patient") ~ "Cancer patient",
      stringr::str_detect(tolower(disease_group_for_clinical), "healthy|control") ~ "Healthy control",
      TRUE ~ disease_group_for_clinical
    ),
    sex_clinical_model = ifelse(is.na(sex_clinical_model), "Not Specified", sex_clinical_model),
    sex_clinical_model = factor(sex_clinical_model)
  )

cancer_df <- clinical_df %>%
  filter(disease_group_for_clinical == "Cancer patient") %>%
  filter(is.finite(age_clinical_model), !is.na(sex_clinical_model))

safe_time_model <- function(score_col) {
  dat <- cancer_df %>%
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
      model_type = "time_from_treatment_start",
      n_model = nrow(dat),
      median_time_days = ifelse(nrow(dat) > 0, median(dat$time_days, na.rm = TRUE), NA_real_),
      beta = NA_real_,
      conf_low = NA_real_,
      conf_high = NA_real_,
      p_value = NA_real_
    ))
  }

  fit <- tryCatch(lm(y ~ log_time_z + age + sex, data = dat), error = function(e) NULL)

  if (is.null(fit)) {
    return(tibble(
      score = score_col,
      model_type = "time_from_treatment_start",
      n_model = nrow(dat),
      median_time_days = median(dat$time_days, na.rm = TRUE),
      beta = NA_real_,
      conf_low = NA_real_,
      conf_high = NA_real_,
      p_value = NA_real_
    ))
  }

  td <- broom::tidy(fit, conf.int = TRUE)
  row <- td %>% filter(term == "log_time_z")

  tibble(
    score = score_col,
    model_type = "time_from_treatment_start",
    n_model = nrow(dat),
    median_time_days = median(dat$time_days, na.rm = TRUE),
    beta = row$estimate[1],
    conf_low = row$conf.low[1],
    conf_high = row$conf.high[1],
    p_value = row$p.value[1]
  )
}

time_from_start_results <- dplyr::bind_rows(lapply(score_cols, safe_time_model)) %>%
  left_join(score_labels, by = "score") %>%
  mutate(
    FDR_global = p.adjust(p_value, method = "BH"),
    direction = ifelse(beta >= 0, "higher_with_longer_time_from_start", "lower_with_longer_time_from_start")
  ) %>%
  arrange(FDR_global, p_value)

standard_csv <- file.path(out_dir, "SDY2583_CP23_time_from_start_results_STEP6.csv")
repair_csv <- file.path(out_dir, "SDY2583_CP23_time_from_start_results_STEP6_REPAIRED.csv")

readr::write_csv(time_from_start_results, standard_csv)
readr::write_csv(time_from_start_results, repair_csv)

# Also update Step 6 RData if it exists, preserving the previous objects.
step6_rdata <- file.path(out_dir, "SDY2583_CP23_STEP6_clinical_annotation.RData")
if (file.exists(step6_rdata)) {
  e <- new.env(parent = emptyenv())
  load(step6_rdata, envir = e)
  e$time_from_start_results <- time_from_start_results
  save(list = ls(e), envir = e, file = step6_rdata)
}

cat("\n============================================================\n")
cat("SDY2583 CP23 STEP 6B COMPLETE: TIME-FROM-START REPAIR\n")
cat("============================================================\n")

cat("\nInput ALL10 matrix:\n")
print(all10_file)

cat("\nCP23 score columns used:\n")
print(score_cols)

cat("\nTime-from-start model results:\n")
print(as.data.frame(time_from_start_results), row.names = FALSE)

cat("\nFiles saved:\n")
print(standard_csv)
print(repair_csv)

cat("\nRData updated if present:\n")
print(step6_rdata)

cat("============================================================\n")
