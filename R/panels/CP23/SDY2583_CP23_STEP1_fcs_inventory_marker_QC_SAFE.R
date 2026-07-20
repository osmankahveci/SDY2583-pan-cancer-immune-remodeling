# ============================================================
# SDY2583 CP23
# STEP 1 SAFE: FCS inventory / marker-channel QC
#
# Purpose:
#   Inventory all CP23 FCS files, read marker/channel structure,
#   verify spillover keywords, detect channel-order and marker-set
#   mismatches, and save a standardized QC package.
#
# Output:
#   outputs/CP23/01_fcs_inventory_marker_QC
#   outputs/CP23/11_RData
# ============================================================

rm(list = ls())
source(file.path(Sys.getenv("SDY2583_REPO_ROOT", unset = "."), "R", "shared", "bootstrap.R"))

# ------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------

cran_pkgs <- c("dplyr", "readr", "stringr", "tibble", "purrr", "tidyr")

for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
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
  library(purrr)
  library(tidyr)
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
n_distinct <- dplyr::n_distinct
n <- dplyr::n

# ------------------------------------------------------------
# 2. Paths
# ------------------------------------------------------------

# Configure with SDY2583_CP23_FCS_DIR; otherwise data/raw/CP23 is used.
candidate_fcs_dirs <- sd_fcs_dir("CP23")

fcs_dir_hits <- candidate_fcs_dirs[
  dir.exists(candidate_fcs_dirs) &
    vapply(candidate_fcs_dirs, function(dd) {
      length(list.files(dd, pattern = "\\.fcs$", recursive = TRUE, ignore.case = TRUE)) > 0
    }, logical(1))
]

if (length(fcs_dir_hits) == 0) {
  stop(
    "CP23 FCS files were not found in: ", sd_fcs_dir("CP23"),
    "\nSet SDY2583_CP23_FCS_DIR to the folder containing CP23 FCS files."
  )
}

selected_fcs_dir <- fcs_dir_hits[1]

analysis_dir <- sd_analysis_dir("CP23")
out_dir <- file.path(analysis_dir, "01_fcs_inventory_marker_QC")
rdata_dir <- file.path(analysis_dir, "11_RData")

dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rdata_dir, recursive = TRUE, showWarnings = FALSE)

fcs_files <- list.files(
  selected_fcs_dir,
  pattern = "\\.fcs$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(fcs_files) == 0) {
  stop("Seçilen CP23 klasöründe FCS dosyası yok: ", selected_fcs_dir)
}

fcs_files <- sort(fcs_files)

# ------------------------------------------------------------
# 3. Helper functions
# ------------------------------------------------------------

extract_subject_id <- function(x) {
  b <- basename(x)
  id <- stringr::str_extract(b, "DBG[0-9]+")
  ifelse(is.na(id), stringr::str_remove(b, "\\.fcs$"), id)
}

get_marker_map_safe <- function(ff) {
  pp <- Biobase::pData(flowCore::parameters(ff))
  channels <- as.character(pp$name)
  markers <- as.character(pp$desc)
  idx <- is.na(markers) | markers == ""
  markers[idx] <- channels[idx]

  tibble(
    parameter_index = seq_along(channels),
    channel_name = channels,
    marker_desc = markers
  )
}

has_spillover_keyword <- function(ff) {
  kk <- flowCore::keyword(ff)
  any(c("SPILL", "$SPILLOVER", "SPILLOVER") %in% names(kk))
}

spillover_keyword_name <- function(ff) {
  kk <- flowCore::keyword(ff)
  nm <- c("SPILL", "$SPILLOVER", "SPILLOVER")
  hit <- nm[nm %in% names(kk)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

read_one_inventory <- function(fp) {
  tryCatch({
    ff <- flowCore::read.FCS(
      fp,
      transformation = FALSE,
      truncate_max_range = FALSE,
      emptyValue = FALSE
    )

    ex <- flowCore::exprs(ff)
    marker_map <- get_marker_map_safe(ff)

    tibble(
      file_name = basename(fp),
      subject_id = extract_subject_id(fp),
      file_path = fp,
      read_ok = TRUE,
      read_error = NA_character_,
      n_events = nrow(ex),
      n_channels = ncol(ex),
      channel_signature = paste(marker_map$channel_name, collapse = "||"),
      marker_signature = paste(marker_map$marker_desc, collapse = "||"),
      marker_set_signature = paste(sort(unique(marker_map$marker_desc)), collapse = "||"),
      has_spillover_keyword = has_spillover_keyword(ff),
      spillover_keyword_name = spillover_keyword_name(ff)
    )
  }, error = function(e) {
    tibble(
      file_name = basename(fp),
      subject_id = extract_subject_id(fp),
      file_path = fp,
      read_ok = FALSE,
      read_error = as.character(e$message),
      n_events = NA_integer_,
      n_channels = NA_integer_,
      channel_signature = NA_character_,
      marker_signature = NA_character_,
      marker_set_signature = NA_character_,
      has_spillover_keyword = NA,
      spillover_keyword_name = NA_character_
    )
  })
}

read_marker_map_one <- function(fp) {
  tryCatch({
    ff <- flowCore::read.FCS(
      fp,
      transformation = FALSE,
      truncate_max_range = FALSE,
      emptyValue = FALSE
    )

    get_marker_map_safe(ff) %>%
      mutate(
        file_name = basename(fp),
        subject_id = extract_subject_id(fp),
        file_path = fp
      ) %>%
      select(file_name, subject_id, file_path, parameter_index, channel_name, marker_desc)
  }, error = function(e) {
    tibble(
      file_name = basename(fp),
      subject_id = extract_subject_id(fp),
      file_path = fp,
      parameter_index = NA_integer_,
      channel_name = NA_character_,
      marker_desc = NA_character_
    )
  })
}

# ------------------------------------------------------------
# 4. Inventory all FCS files
# ------------------------------------------------------------

cat("\nReading CP23 FCS inventory...\n")

inventory_rows <- vector("list", length(fcs_files))

for (i in seq_along(fcs_files)) {
  if (i %% 50 == 0) cat("  processed ", i, " / ", length(fcs_files), "\n", sep = "")
  inventory_rows[[i]] <- read_one_inventory(fcs_files[i])
}

fcs_inventory <- bind_rows(inventory_rows)

failed_files <- fcs_inventory %>%
  filter(read_ok == FALSE) %>%
  select(file_name, subject_id, file_path, read_error)

read_ok_inventory <- fcs_inventory %>%
  filter(read_ok == TRUE)

if (nrow(read_ok_inventory) == 0) {
  stop("Hiçbir CP23 FCS dosyası okunamadı.")
}

# Reference file:
# Prefer the file with the modal channel count and highest event count.
modal_n_channels <- read_ok_inventory %>%
  count(n_channels, name = "n_files") %>%
  arrange(desc(n_files), n_channels) %>%
  slice(1) %>%
  pull(n_channels)

reference_file <- read_ok_inventory %>%
  filter(n_channels == modal_n_channels) %>%
  arrange(desc(n_events), file_name) %>%
  slice(1) %>%
  pull(file_path)

reference_ff <- flowCore::read.FCS(
  reference_file,
  transformation = FALSE,
  truncate_max_range = FALSE,
  emptyValue = FALSE
)

reference_marker_map <- get_marker_map_safe(reference_ff)

reference_channel_signature <- paste(reference_marker_map$channel_name, collapse = "||")
reference_marker_signature <- paste(reference_marker_map$marker_desc, collapse = "||")
reference_marker_set_signature <- paste(sort(unique(reference_marker_map$marker_desc)), collapse = "||")

fcs_inventory <- fcs_inventory %>%
  mutate(
    channel_order_mismatch_file =
      read_ok == TRUE &
      !is.na(channel_signature) &
      channel_signature != reference_channel_signature,
    marker_order_mismatch_file =
      read_ok == TRUE &
      !is.na(marker_signature) &
      marker_signature != reference_marker_signature,
    marker_set_mismatch_file =
      read_ok == TRUE &
      !is.na(marker_set_signature) &
      marker_set_signature != reference_marker_set_signature
  )

# ------------------------------------------------------------
# 5. Marker maps and summaries
# ------------------------------------------------------------

cat("\nReading CP23 marker maps...\n")

marker_map_rows <- vector("list", length(fcs_files))

for (i in seq_along(fcs_files)) {
  if (i %% 50 == 0) cat("  marker maps processed ", i, " / ", length(fcs_files), "\n", sep = "")
  marker_map_rows[[i]] <- read_marker_map_one(fcs_files[i])
}

all_marker_maps <- bind_rows(marker_map_rows)

marker_presence_summary <- all_marker_maps %>%
  filter(!is.na(marker_desc)) %>%
  group_by(marker_desc) %>%
  summarise(
    n_files = n_distinct(file_name),
    example_channel = first(channel_name),
    channels_seen = paste(sort(unique(channel_name)), collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(marker_desc)

channel_count_summary <- fcs_inventory %>%
  filter(read_ok == TRUE) %>%
  count(n_channels, name = "n_files") %>%
  arrange(n_channels)

spillover_summary <- fcs_inventory %>%
  filter(read_ok == TRUE) %>%
  count(has_spillover_keyword, spillover_keyword_name, name = "n_files") %>%
  arrange(desc(n_files))

mismatch_files <- fcs_inventory %>%
  filter(
    read_ok == TRUE,
    channel_order_mismatch_file == TRUE |
      marker_order_mismatch_file == TRUE |
      marker_set_mismatch_file == TRUE
  ) %>%
  select(
    file_name,
    subject_id,
    n_events,
    n_channels,
    channel_order_mismatch_file,
    marker_order_mismatch_file,
    marker_set_mismatch_file,
    has_spillover_keyword,
    spillover_keyword_name
  ) %>%
  arrange(file_name)

source_folder_summary <- tibble(
  selected_fcs_dir = selected_fcs_dir,
  n_fcs_files_found = length(fcs_files),
  analysis_dir = analysis_dir,
  out_dir = out_dir,
  rdata_dir = rdata_dir
)

fcs_inventory_summary <- fcs_inventory %>%
  summarise(
    n_fcs_files = n(),
    n_unique_subjects = n_distinct(subject_id),
    n_read_ok = sum(read_ok == TRUE, na.rm = TRUE),
    n_read_failed = sum(read_ok == FALSE, na.rm = TRUE),
    n_with_spillover_keyword = sum(has_spillover_keyword == TRUE, na.rm = TRUE),
    n_without_spillover_keyword = sum(has_spillover_keyword == FALSE, na.rm = TRUE),
    n_distinct_channel_counts = n_distinct(n_channels[read_ok == TRUE]),
    modal_n_channels = modal_n_channels,
    median_events = median(n_events, na.rm = TRUE),
    min_events = min(n_events, na.rm = TRUE),
    q1_events = as.numeric(quantile(n_events, 0.25, na.rm = TRUE)),
    mean_events = mean(n_events, na.rm = TRUE),
    q3_events = as.numeric(quantile(n_events, 0.75, na.rm = TRUE)),
    max_events = max(n_events, na.rm = TRUE),
    n_channel_order_mismatch_files = sum(channel_order_mismatch_file == TRUE, na.rm = TRUE),
    n_marker_order_mismatch_files = sum(marker_order_mismatch_file == TRUE, na.rm = TRUE),
    n_marker_set_mismatch_files = sum(marker_set_mismatch_file == TRUE, na.rm = TRUE)
  )

# ------------------------------------------------------------
# 6. Save outputs
# ------------------------------------------------------------

write_csv(source_folder_summary, file.path(out_dir, "SDY2583_CP23_source_folder_summary_STEP1.csv"))
write_csv(fcs_inventory, file.path(out_dir, "SDY2583_CP23_fcs_inventory_STEP1.csv"))
write_csv(fcs_inventory_summary, file.path(out_dir, "SDY2583_CP23_fcs_inventory_summary_STEP1.csv"))
write_csv(channel_count_summary, file.path(out_dir, "SDY2583_CP23_channel_count_summary_STEP1.csv"))
write_csv(spillover_summary, file.path(out_dir, "SDY2583_CP23_spillover_summary_STEP1.csv"))
write_csv(reference_marker_map, file.path(out_dir, "SDY2583_CP23_reference_marker_map_STEP1.csv"))
write_csv(all_marker_maps, file.path(out_dir, "SDY2583_CP23_all_marker_maps_STEP1.csv"))
write_csv(marker_presence_summary, file.path(out_dir, "SDY2583_CP23_marker_presence_summary_STEP1.csv"))
write_csv(mismatch_files, file.path(out_dir, "SDY2583_CP23_mismatch_files_STEP1.csv"))
write_csv(failed_files, file.path(out_dir, "SDY2583_CP23_failed_files_STEP1.csv"))

save(
  selected_fcs_dir,
  fcs_files,
  analysis_dir,
  out_dir,
  rdata_dir,
  fcs_inventory,
  fcs_inventory_summary,
  source_folder_summary,
  channel_count_summary,
  spillover_summary,
  reference_file,
  reference_marker_map,
  reference_channel_signature,
  reference_marker_signature,
  reference_marker_set_signature,
  all_marker_maps,
  marker_presence_summary,
  mismatch_files,
  failed_files,
  file = file.path(rdata_dir, "SDY2583_CP23_STEP1_fcs_inventory_marker_QC.RData")
)

# ------------------------------------------------------------
# 7. Console output
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SDY2583 CP23 STEP 1 COMPLETE: FCS INVENTORY / MARKER QC\n")
cat("============================================================\n")

cat("\nSource folder summary:\n")
print(as.data.frame(source_folder_summary), row.names = FALSE)

cat("\nFCS inventory summary:\n")
print(as.data.frame(fcs_inventory_summary), row.names = FALSE)

cat("\nChannel count summary:\n")
print(as.data.frame(channel_count_summary), row.names = FALSE)

cat("\nSpillover summary:\n")
print(as.data.frame(spillover_summary), row.names = FALSE)

cat("\nReference file:\n")
print(basename(reference_file))

cat("\nReference marker map:\n")
print(as.data.frame(reference_marker_map), row.names = FALSE)

cat("\nMarker presence summary:\n")
print(as.data.frame(marker_presence_summary), row.names = FALSE)

cat("\nMismatch files:\n")
print(as.data.frame(mismatch_files), row.names = FALSE)

cat("\nFailed files:\n")
print(as.data.frame(failed_files), row.names = FALSE)

cat("\nOutputs saved in:\n")
print(out_dir)

cat("============================================================\n")
