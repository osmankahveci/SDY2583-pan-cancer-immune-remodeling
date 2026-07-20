# ============================================================
# SDY2583 CP22
# STEP 1 SAFE SCRIPT v2: FCS inventory and marker-channel QC
# Fixes filter/read_ok namespace issue by avoiding unqualified filter()
# ============================================================

rm(list = ls())
source(file.path(Sys.getenv("SDY2583_REPO_ROOT", unset = "."), "R", "shared", "bootstrap.R"))

# ------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------

cran_pkgs <- c("dplyr", "readr", "stringr", "tibble")

for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("flowCore", quietly = TRUE)) {
  BiocManager::install("flowCore", ask = FALSE, update = FALSE)
}

suppressPackageStartupMessages({
  library(flowCore)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
})

# ------------------------------------------------------------
# 2. Locate CP22 FCS folder
# ------------------------------------------------------------

fcs_dir <- sd_fcs_dir("CP22")
if (!dir.exists(fcs_dir)) stop("Configured CP22 FCS directory does not exist: ", fcs_dir)

cat("\nSelected CP22 FCS folder:\n")
print(fcs_dir)

# ------------------------------------------------------------
# 3. Output folders
# ------------------------------------------------------------

analysis_dir <- sd_analysis_dir("CP22")
out_dir <- file.path(analysis_dir, "01_fcs_inventory")
rdata_dir <- file.path(analysis_dir, "11_RData")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rdata_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 4. FCS file list
# ------------------------------------------------------------

fcs_files <- list.files(
  fcs_dir,
  pattern = "\\.fcs$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(fcs_files) == 0) {
  stop("FCS file bulunamadı.")
}

extract_subject_id <- function(x) {
  b <- basename(x)
  id <- stringr::str_extract(b, "DBG[0-9]+")
  ifelse(is.na(id), stringr::str_remove(b, "\\.fcs$"), id)
}

fcs_inventory <- tibble::tibble(
  file_path = fcs_files,
  file_name = basename(fcs_files),
  subject_id = extract_subject_id(fcs_files),
  file_size_mb = as.numeric(file.info(fcs_files)$size) / 1024^2
)

# ------------------------------------------------------------
# 5. First file marker-channel map
# ------------------------------------------------------------

first_file <- fcs_files[1]

cat("\nReading first FCS file:\n")
print(first_file)

ff_first <- flowCore::read.FCS(
  first_file,
  transformation = FALSE,
  truncate_max_range = FALSE
)

expr_first <- flowCore::exprs(ff_first)
param_first <- Biobase::pData(flowCore::parameters(ff_first))
keywords_first <- flowCore::keyword(ff_first)

first_channels <- as.character(param_first$name)
first_markers <- as.character(param_first$desc)

idx_empty_marker <- is.na(first_markers) | first_markers == ""
first_markers[idx_empty_marker] <- first_channels[idx_empty_marker]

marker_map <- tibble::tibble(
  parameter_index = seq_along(first_channels),
  channel_name = first_channels,
  marker_desc = first_markers,
  is_scatter_or_time = grepl("FSC|SSC|Time", first_channels, ignore.case = TRUE) |
    grepl("FSC|SSC|Time", first_markers, ignore.case = TRUE)
)

spill_keys <- c("$SPILLOVER", "SPILLOVER", "SPILL")
available_spill_keywords <- spill_keys[spill_keys %in% names(keywords_first)]
spillover_present_first <- length(available_spill_keywords) > 0

# ------------------------------------------------------------
# 6. Per-file QC
# ------------------------------------------------------------

inspect_one <- function(fp) {
  tryCatch({
    ff <- flowCore::read.FCS(
      fp,
      transformation = FALSE,
      truncate_max_range = FALSE
    )

    ex <- flowCore::exprs(ff)
    pp <- Biobase::pData(flowCore::parameters(ff))
    kk <- flowCore::keyword(ff)

    ch <- as.character(pp$name)
    mk <- as.character(pp$desc)
    idx <- is.na(mk) | mk == ""
    mk[idx] <- ch[idx]

    tibble::tibble(
      file_path = fp,
      file_name = basename(fp),
      subject_id = extract_subject_id(fp),
      read_ok = TRUE,
      n_events = as.integer(nrow(ex)),
      n_channels = as.integer(ncol(ex)),
      has_spillover_keyword = any(spill_keys %in% names(kk)),
      channel_names_match_first_file = identical(ch, first_channels),
      marker_desc_match_first_file = identical(mk, first_markers),
      error_message = NA_character_
    )
  }, error = function(e) {
    tibble::tibble(
      file_path = fp,
      file_name = basename(fp),
      subject_id = extract_subject_id(fp),
      read_ok = FALSE,
      n_events = NA_integer_,
      n_channels = NA_integer_,
      has_spillover_keyword = NA,
      channel_names_match_first_file = NA,
      marker_desc_match_first_file = NA,
      error_message = as.character(e$message)
    )
  })
}

cat("\nInspecting CP22 FCS files...\n")

qc_list <- vector("list", length(fcs_files))

for (i in seq_along(fcs_files)) {
  if (i %% 50 == 0) {
    cat("  inspected", i, "of", length(fcs_files), "files\n")
  }
  qc_list[[i]] <- inspect_one(fcs_files[i])
}

fcs_qc <- dplyr::bind_rows(qc_list)

# Hard diagnostic before summaries
cat("\nFCS QC columns:\n")
print(names(fcs_qc))
cat("\nFCS QC dimensions:\n")
print(dim(fcs_qc))

if (!("read_ok" %in% names(fcs_qc))) {
  stop("Internal error: fcs_qc does not contain read_ok. Inspect qc_list object.")
}

# ------------------------------------------------------------
# 7. Summaries
# ------------------------------------------------------------

inventory_summary <- tibble::tibble(
  n_fcs_files = length(fcs_files),
  n_unique_subjects = dplyr::n_distinct(fcs_inventory$subject_id),
  n_read_ok = sum(fcs_qc$read_ok == TRUE, na.rm = TRUE),
  n_read_failed = sum(fcs_qc$read_ok == FALSE, na.rm = TRUE),
  n_channel_mismatch = sum(fcs_qc$channel_names_match_first_file == FALSE, na.rm = TRUE),
  n_marker_mismatch = sum(fcs_qc$marker_desc_match_first_file == FALSE, na.rm = TRUE),
  n_with_spillover_keyword = sum(fcs_qc$has_spillover_keyword == TRUE, na.rm = TRUE),
  median_file_size_mb = median(fcs_inventory$file_size_mb, na.rm = TRUE),
  median_events = median(fcs_qc$n_events, na.rm = TRUE),
  min_events = min(fcs_qc$n_events, na.rm = TRUE),
  max_events = max(fcs_qc$n_events, na.rm = TRUE)
)

fcs_qc_ok <- fcs_qc[fcs_qc$read_ok == TRUE & !is.na(fcs_qc$read_ok), , drop = FALSE]

event_distribution <- tibble::tibble(
  n_files = nrow(fcs_qc_ok),
  mean_events = mean(fcs_qc_ok$n_events, na.rm = TRUE),
  sd_events = stats::sd(fcs_qc_ok$n_events, na.rm = TRUE),
  median_events = median(fcs_qc_ok$n_events, na.rm = TRUE),
  q1_events = as.numeric(stats::quantile(fcs_qc_ok$n_events, 0.25, na.rm = TRUE)),
  q3_events = as.numeric(stats::quantile(fcs_qc_ok$n_events, 0.75, na.rm = TRUE)),
  min_events = min(fcs_qc_ok$n_events, na.rm = TRUE),
  max_events = max(fcs_qc_ok$n_events, na.rm = TRUE)
)

problem_files <- fcs_qc[
  (fcs_qc$read_ok == FALSE) |
    (fcs_qc$channel_names_match_first_file == FALSE) |
    (fcs_qc$marker_desc_match_first_file == FALSE),
  ,
  drop = FALSE
]

problem_files <- problem_files[!is.na(problem_files$read_ok) |
                                 !is.na(problem_files$channel_names_match_first_file) |
                                 !is.na(problem_files$marker_desc_match_first_file), , drop = FALSE]

# ------------------------------------------------------------
# 8. Save outputs
# ------------------------------------------------------------

readr::write_csv(fcs_inventory, file.path(out_dir, "SDY2583_CP22_FCS_inventory_STEP1.csv"))
readr::write_csv(marker_map, file.path(out_dir, "SDY2583_CP22_marker_channel_map_STEP1.csv"))
readr::write_csv(fcs_qc, file.path(out_dir, "SDY2583_CP22_FCS_QC_per_file_STEP1.csv"))
readr::write_csv(inventory_summary, file.path(out_dir, "SDY2583_CP22_FCS_inventory_summary_STEP1.csv"))
readr::write_csv(event_distribution, file.path(out_dir, "SDY2583_CP22_event_distribution_STEP1.csv"))
readr::write_csv(problem_files, file.path(out_dir, "SDY2583_CP22_problem_files_STEP1.csv"))

save(
  fcs_dir,
  fcs_files,
  fcs_inventory,
  marker_map,
  fcs_qc,
  inventory_summary,
  event_distribution,
  problem_files,
  first_file,
  spillover_present_first,
  available_spill_keywords,
  file = file.path(rdata_dir, "SDY2583_CP22_STEP1_fcs_inventory_marker_QC.RData")
)

# ------------------------------------------------------------
# 9. Console output
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SDY2583 CP22 STEP 1 COMPLETE\n")
cat("============================================================\n")

cat("\nInventory summary:\n")
print(as.data.frame(inventory_summary), row.names = FALSE)

cat("\nMarker-channel map:\n")
print(as.data.frame(marker_map), row.names = FALSE)

cat("\nEvent distribution:\n")
print(as.data.frame(event_distribution), row.names = FALSE)

cat("\nProblem files:\n")
if (nrow(problem_files) == 0) {
  cat("No problem files detected.\n")
} else {
  print(as.data.frame(problem_files), row.names = FALSE)
}

cat("\nSpillover keyword present in first file:\n")
print(spillover_present_first)

cat("\nAvailable spillover keyword names in first file:\n")
print(available_spill_keywords)

cat("\nFiles saved in:\n")
print(out_dir)

cat("============================================================\n")
