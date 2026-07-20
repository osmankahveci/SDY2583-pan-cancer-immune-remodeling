# ============================================================
# SDY2583 CP16
# STEP 1 SAFE: FCS inventory / marker-channel QC
#
# Purpose:
#   - Locate CP16 FCS files
#   - Read all FCS headers/events safely
#   - Build marker-channel map
#   - Check channel count, marker presence, spillover keyword
#   - Detect channel-order and marker-order mismatches
#
# Output:
#   outputs/CP16/01_fcs_inventory_marker_QC
#   outputs/CP16/11_RData
#
# After this step, send console output:
#   Source folder summary
#   FCS inventory summary
#   Reference marker map
#   Marker presence summary
#   Mismatch files
#   Failed files
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

analysis_dir <- sd_analysis_dir("CP16")
out_dir <- file.path(analysis_dir, "01_fcs_inventory_marker_QC")
rdata_dir <- file.path(analysis_dir, "11_RData")

dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rdata_dir, recursive = TRUE, showWarnings = FALSE)

# Configure with SDY2583_CP16_FCS_DIR; otherwise data/raw/CP16 is used.
candidate_dirs <- sd_fcs_dir("CP16")
existing_dirs <- candidate_dirs[dir.exists(candidate_dirs)]

find_cp16_files <- function(search_dir) {
  all_fcs <- list.files(
    search_dir,
    pattern = "\\.fcs$",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )
  all_fcs[grepl("CP16", basename(all_fcs), ignore.case = TRUE)]
}

fcs_files <- character()

for (dd in existing_dirs) {
  tmp <- find_cp16_files(dd)
  if (length(tmp) > length(fcs_files)) {
    fcs_files <- tmp
    selected_fcs_dir <- dd
  }
}

fcs_files <- sort(unique(fcs_files))

if (length(fcs_files) == 0) {
  stop(
    "CP16 FCS files were not found in: ", sd_fcs_dir("CP16"),
    "\nSet SDY2583_CP16_FCS_DIR to the folder containing CP16 FCS files."
  )
}

# ------------------------------------------------------------
# 3. Helper functions
# ------------------------------------------------------------

extract_subject_id <- function(x) {
  b <- basename(x)
  id <- stringr::str_extract(b, "DBG[0-9]+")
  ifelse(is.na(id), stringr::str_remove(b, "\\.fcs$"), id)
}

extract_panel_from_filename <- function(x) {
  b <- basename(x)
  panel <- stringr::str_extract(b, "CP[0-9]+")
  ifelse(is.na(panel), NA_character_, panel)
}

get_marker_map <- function(ff) {
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

has_spill_keyword <- function(ff) {
  kk <- flowCore::keyword(ff)
  any(c("SPILL", "$SPILLOVER", "SPILLOVER") %in% names(kk))
}

get_spill_keyword_name <- function(ff) {
  kk <- flowCore::keyword(ff)
  hit <- c("SPILL", "$SPILLOVER", "SPILLOVER")[c("SPILL", "$SPILLOVER", "SPILLOVER") %in% names(kk)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

safe_read_header <- function(fp) {
  tryCatch({
    ff <- flowCore::read.FCS(
      fp,
      transformation = FALSE,
      truncate_max_range = FALSE
    )

    mm <- get_marker_map(ff)
    ex <- flowCore::exprs(ff)

    tibble(
      file_name = basename(fp),
      file_path = fp,
      subject_id = extract_subject_id(fp),
      panel_from_file = extract_panel_from_filename(fp),
      read_ok = TRUE,
      read_error = NA_character_,
      n_events = nrow(ex),
      n_channels = ncol(ex),
      has_spillover_keyword = has_spill_keyword(ff),
      spillover_keyword_name = get_spill_keyword_name(ff),
      channel_signature = paste(mm$channel_name, collapse = "|"),
      marker_signature = paste(mm$marker_desc, collapse = "|"),
      marker_set_signature = paste(sort(unique(mm$marker_desc)), collapse = "|")
    )
  }, error = function(e) {
    tibble(
      file_name = basename(fp),
      file_path = fp,
      subject_id = extract_subject_id(fp),
      panel_from_file = extract_panel_from_filename(fp),
      read_ok = FALSE,
      read_error = as.character(e$message),
      n_events = NA_integer_,
      n_channels = NA_integer_,
      has_spillover_keyword = NA,
      spillover_keyword_name = NA_character_,
      channel_signature = NA_character_,
      marker_signature = NA_character_,
      marker_set_signature = NA_character_
    )
  })
}

safe_get_marker_map_file <- function(fp) {
  tryCatch({
    ff <- flowCore::read.FCS(
      fp,
      transformation = FALSE,
      truncate_max_range = FALSE
    )
    get_marker_map(ff) %>%
      mutate(
        file_name = basename(fp),
        file_path = fp,
        subject_id = extract_subject_id(fp)
      )
  }, error = function(e) {
    tibble(
      parameter_index = NA_integer_,
      channel_name = NA_character_,
      marker_desc = NA_character_,
      file_name = basename(fp),
      file_path = fp,
      subject_id = extract_subject_id(fp),
      marker_map_error = as.character(e$message)
    )
  })
}

# ------------------------------------------------------------
# 4. Read inventory
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SDY2583 CP16 STEP 1 STARTED: FCS INVENTORY / MARKER QC\n")
cat("============================================================\n")

cat("\nSelected FCS search folder:\n")
print(selected_fcs_dir)

cat("\nNumber of CP16 FCS files found:\n")
print(length(fcs_files))

inventory_rows <- vector("list", length(fcs_files))

for (i in seq_along(fcs_files)) {
  if (i %% 50 == 0) cat("Read", i, "of", length(fcs_files), "FCS files\n")
  inventory_rows[[i]] <- safe_read_header(fcs_files[i])
}

fcs_inventory <- bind_rows(inventory_rows)

# ------------------------------------------------------------
# 5. Reference marker map and mismatch checks
# ------------------------------------------------------------

ok_files <- fcs_inventory %>%
  filter(read_ok == TRUE)

if (nrow(ok_files) == 0) {
  stop("Hiçbir CP16 FCS dosyası okunamadı.")
}

# Choose the first successfully read file with the modal channel count as reference.
modal_n_channels <- ok_files %>%
  count(n_channels, name = "n") %>%
  arrange(desc(n), n_channels) %>%
  slice(1) %>%
  pull(n_channels)

reference_file <- ok_files %>%
  filter(n_channels == modal_n_channels) %>%
  arrange(file_name) %>%
  slice(1) %>%
  pull(file_path)

reference_marker_map <- safe_get_marker_map_file(reference_file) %>%
  select(parameter_index, channel_name, marker_desc)

reference_channel_signature <- paste(reference_marker_map$channel_name, collapse = "|")
reference_marker_signature <- paste(reference_marker_map$marker_desc, collapse = "|")
reference_marker_set_signature <- paste(sort(unique(reference_marker_map$marker_desc)), collapse = "|")

fcs_inventory <- fcs_inventory %>%
  mutate(
    channel_order_mismatch_file = ifelse(
      read_ok == TRUE,
      channel_signature != reference_channel_signature,
      NA
    ),
    marker_order_mismatch_file = ifelse(
      read_ok == TRUE,
      marker_signature != reference_marker_signature,
      NA
    ),
    marker_set_mismatch_file = ifelse(
      read_ok == TRUE,
      marker_set_signature != reference_marker_set_signature,
      NA
    )
  )

# Full marker maps for all files. This is useful for CP16 Step1B if needed.
marker_map_all_files <- bind_rows(lapply(fcs_files, safe_get_marker_map_file))

marker_presence_summary <- marker_map_all_files %>%
  filter(!is.na(marker_desc)) %>%
  group_by(marker_desc) %>%
  summarise(
    n_files = n_distinct(file_name),
    example_channel = dplyr::first(channel_name),
    channels_seen = paste(sort(unique(channel_name)), collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(marker_desc)

channel_presence_summary <- marker_map_all_files %>%
  filter(!is.na(channel_name)) %>%
  group_by(channel_name) %>%
  summarise(
    n_files = n_distinct(file_name),
    markers_seen = paste(sort(unique(marker_desc)), collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(channel_name)

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

failed_files <- fcs_inventory %>%
  filter(read_ok == FALSE) %>%
  select(file_name, subject_id, file_path, read_error)

# ------------------------------------------------------------
# 6. Summary tables
# ------------------------------------------------------------

fcs_inventory_summary <- tibble(
  n_fcs_files = nrow(fcs_inventory),
  n_unique_subjects = n_distinct(fcs_inventory$subject_id),
  n_read_ok = sum(fcs_inventory$read_ok == TRUE, na.rm = TRUE),
  n_read_failed = sum(fcs_inventory$read_ok == FALSE, na.rm = TRUE),
  n_with_spillover_keyword = sum(fcs_inventory$has_spillover_keyword == TRUE, na.rm = TRUE),
  n_without_spillover_keyword = sum(fcs_inventory$has_spillover_keyword == FALSE, na.rm = TRUE),
  n_distinct_channel_counts = n_distinct(fcs_inventory$n_channels[!is.na(fcs_inventory$n_channels)]),
  modal_n_channels = modal_n_channels,
  median_events = median(fcs_inventory$n_events, na.rm = TRUE),
  min_events = min(fcs_inventory$n_events, na.rm = TRUE),
  q1_events = as.numeric(quantile(fcs_inventory$n_events, 0.25, na.rm = TRUE)),
  mean_events = mean(fcs_inventory$n_events, na.rm = TRUE),
  q3_events = as.numeric(quantile(fcs_inventory$n_events, 0.75, na.rm = TRUE)),
  max_events = max(fcs_inventory$n_events, na.rm = TRUE),
  n_channel_order_mismatch_files = sum(fcs_inventory$channel_order_mismatch_file == TRUE, na.rm = TRUE),
  n_marker_order_mismatch_files = sum(fcs_inventory$marker_order_mismatch_file == TRUE, na.rm = TRUE),
  n_marker_set_mismatch_files = sum(fcs_inventory$marker_set_mismatch_file == TRUE, na.rm = TRUE)
)

channel_count_summary <- fcs_inventory %>%
  count(n_channels, name = "n_files") %>%
  arrange(n_channels)

spillover_summary <- fcs_inventory %>%
  count(has_spillover_keyword, spillover_keyword_name, name = "n_files") %>%
  arrange(desc(n_files))

event_summary_by_channel_count <- fcs_inventory %>%
  filter(read_ok == TRUE) %>%
  group_by(n_channels) %>%
  summarise(
    n_files = n(),
    median_events = median(n_events, na.rm = TRUE),
    min_events = min(n_events, na.rm = TRUE),
    max_events = max(n_events, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(n_channels)

source_folder_summary <- tibble(
  selected_fcs_dir = selected_fcs_dir,
  n_fcs_files_found = length(fcs_files),
  analysis_dir = analysis_dir,
  out_dir = out_dir,
  rdata_dir = rdata_dir
)

# ------------------------------------------------------------
# 7. Save outputs
# ------------------------------------------------------------

write_csv(source_folder_summary, file.path(out_dir, "SDY2583_CP16_source_folder_summary_STEP1.csv"))
write_csv(fcs_inventory, file.path(out_dir, "SDY2583_CP16_fcs_inventory_STEP1.csv"))
write_csv(fcs_inventory_summary, file.path(out_dir, "SDY2583_CP16_fcs_inventory_summary_STEP1.csv"))
write_csv(reference_marker_map, file.path(out_dir, "SDY2583_CP16_reference_marker_map_STEP1.csv"))
write_csv(marker_map_all_files, file.path(out_dir, "SDY2583_CP16_marker_map_all_files_STEP1.csv"))
write_csv(marker_presence_summary, file.path(out_dir, "SDY2583_CP16_marker_presence_summary_STEP1.csv"))
write_csv(channel_presence_summary, file.path(out_dir, "SDY2583_CP16_channel_presence_summary_STEP1.csv"))
write_csv(channel_count_summary, file.path(out_dir, "SDY2583_CP16_channel_count_summary_STEP1.csv"))
write_csv(spillover_summary, file.path(out_dir, "SDY2583_CP16_spillover_summary_STEP1.csv"))
write_csv(event_summary_by_channel_count, file.path(out_dir, "SDY2583_CP16_event_summary_by_channel_count_STEP1.csv"))
write_csv(mismatch_files, file.path(out_dir, "SDY2583_CP16_mismatch_files_STEP1.csv"))
write_csv(failed_files, file.path(out_dir, "SDY2583_CP16_failed_files_STEP1.csv"))

save(
  fcs_files,
  selected_fcs_dir,
  analysis_dir,
  out_dir,
  rdata_dir,
  fcs_inventory,
  fcs_inventory_summary,
  reference_file,
  reference_marker_map,
  marker_map_all_files,
  marker_presence_summary,
  channel_presence_summary,
  channel_count_summary,
  spillover_summary,
  event_summary_by_channel_count,
  mismatch_files,
  failed_files,
  file = file.path(rdata_dir, "SDY2583_CP16_STEP1_fcs_inventory_marker_QC.RData")
)

# ------------------------------------------------------------
# 8. Console output
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SDY2583 CP16 STEP 1 COMPLETE: FCS INVENTORY / MARKER QC\n")
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
