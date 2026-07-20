# ============================================================
# SDY2583 CP7
# STEP 1 RECONSTRUCTED: FCS inventory / marker-channel QC
#
# Reconstruction basis:
#   - archived 850-file CP7 inventory;
#   - archived marker-channel and channel-pattern tables;
#   - CP7 Methods/Results record;
#   - shared QC conventions used by the recovered panel scripts.
#
# This is reconstructed source, not the original archived script.
# Validate generated tables against the archived CP7 outputs before release.
# ============================================================

rm(list = ls())
source(file.path(Sys.getenv("SDY2583_REPO_ROOT", unset = "."), "R", "shared", "bootstrap.R"))

cran_pkgs <- c("dplyr", "readr", "stringr", "tibble", "purrr", "tidyr")
for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
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

analysis_dir <- sd_analysis_dir("CP7")
out_inventory_dir <- file.path(analysis_dir, "00_file_inventory")
out_qc_dir <- file.path(analysis_dir, "01_channel_marker_QC")
rdata_dir <- file.path(analysis_dir, "11_RData")
dir.create(out_inventory_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_qc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rdata_dir, recursive = TRUE, showWarnings = FALSE)

fcs_dir <- sd_fcs_dir("CP7")
fcs_files <- list.files(
  fcs_dir,
  pattern = "\\.fcs$",
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
)
fcs_files <- sort(unique(fcs_files[grepl("CP7", basename(fcs_files), ignore.case = TRUE)]))
if (length(fcs_files) == 0L) {
  stop(
    "No CP7 FCS files were found in: ", fcs_dir,
    "\nSet SDY2583_CP7_FCS_DIR to the folder containing the CP7 FCS files."
  )
}

extract_subject_id <- function(path) {
  id <- stringr::str_extract(basename(path), "DBG[0-9]+")
  ifelse(is.na(id), stringr::str_remove(basename(path), "\\.fcs$"), id)
}

get_marker_map <- function(ff) {
  pp <- Biobase::pData(flowCore::parameters(ff))
  marker <- as.character(pp$desc)
  channel <- as.character(pp$name)
  marker[is.na(marker) | marker == ""] <- channel[is.na(marker) | marker == ""]
  tibble(
    parameter_index = seq_along(channel),
    parameter = channel,
    marker_desc = marker,
    range = suppressWarnings(as.numeric(pp$range)),
    minRange = suppressWarnings(as.numeric(pp$minRange)),
    maxRange = suppressWarnings(as.numeric(pp$maxRange))
  )
}

spill_keyword_name <- function(ff) {
  keys <- names(flowCore::keyword(ff))
  hit <- c("SPILL", "$SPILLOVER", "SPILLOVER")
  hit <- hit[hit %in% keys]
  if (length(hit) == 0L) NA_character_ else hit[[1L]]
}

safe_header <- function(path) {
  tryCatch({
    ff <- flowCore::read.FCS(path, transformation = FALSE, truncate_max_range = FALSE)
    mm <- get_marker_map(ff)
    ex <- flowCore::exprs(ff)
    tibble(
      file_name = basename(path),
      file_path = normalizePath(path, mustWork = FALSE),
      subject_id = extract_subject_id(path),
      panel = "CP7",
      file_size_MB = round(file.info(path)$size / 1024^2, 3),
      read_ok = TRUE,
      read_error = NA_character_,
      n_events = nrow(ex),
      n_channels = ncol(ex),
      spillover_keyword = spill_keyword_name(ff),
      channel_signature = paste(mm$parameter, collapse = "|"),
      marker_signature = paste(mm$marker_desc, collapse = "|"),
      marker_set_signature = paste(sort(unique(mm$marker_desc)), collapse = "|")
    )
  }, error = function(e) {
    tibble(
      file_name = basename(path),
      file_path = normalizePath(path, mustWork = FALSE),
      subject_id = extract_subject_id(path),
      panel = "CP7",
      file_size_MB = round(file.info(path)$size / 1024^2, 3),
      read_ok = FALSE,
      read_error = conditionMessage(e),
      n_events = NA_integer_,
      n_channels = NA_integer_,
      spillover_keyword = NA_character_,
      channel_signature = NA_character_,
      marker_signature = NA_character_,
      marker_set_signature = NA_character_
    )
  })
}

safe_marker_map <- function(path) {
  tryCatch({
    ff <- flowCore::read.FCS(path, transformation = FALSE, truncate_max_range = FALSE)
    get_marker_map(ff) %>%
      mutate(
        file_name = basename(path),
        file_path = normalizePath(path, mustWork = FALSE),
        subject_id = extract_subject_id(path)
      )
  }, error = function(e) {
    tibble(
      parameter_index = NA_integer_, parameter = NA_character_, marker_desc = NA_character_,
      range = NA_real_, minRange = NA_real_, maxRange = NA_real_,
      file_name = basename(path), file_path = normalizePath(path, mustWork = FALSE),
      subject_id = extract_subject_id(path), marker_map_error = conditionMessage(e)
    )
  })
}

message("Reading ", length(fcs_files), " CP7 FCS files...")
fcs_inventory_qc <- purrr::map_dfr(fcs_files, safe_header)
if (!any(fcs_inventory_qc$read_ok)) stop("No CP7 FCS file could be read.")

ok <- fcs_inventory_qc %>% filter(read_ok)
modal_n_channels <- ok %>% count(n_channels, name = "n") %>% arrange(desc(n), n_channels) %>% slice(1) %>% pull(n_channels)
reference_path <- ok %>% filter(n_channels == modal_n_channels) %>% arrange(file_name) %>% slice(1) %>% pull(file_path)
reference_map <- safe_marker_map(reference_path) %>% select(parameter_index, parameter, marker_desc, range, minRange, maxRange)
reference_channel_signature <- paste(reference_map$parameter, collapse = "|")
reference_marker_signature <- paste(reference_map$marker_desc, collapse = "|")
reference_marker_set_signature <- paste(sort(unique(reference_map$marker_desc)), collapse = "|")

fcs_inventory_qc <- fcs_inventory_qc %>%
  mutate(
    channel_order_mismatch = read_ok & channel_signature != reference_channel_signature,
    marker_order_mismatch = read_ok & marker_signature != reference_marker_signature,
    marker_set_mismatch = read_ok & marker_set_signature != reference_marker_set_signature
  )

marker_map_all_files <- purrr::map_dfr(fcs_files, safe_marker_map)
channel_patterns <- fcs_inventory_qc %>%
  filter(read_ok) %>%
  distinct(channel_signature, marker_signature) %>%
  arrange(channel_signature, marker_signature) %>%
  mutate(cp7_channel_pattern_id = row_number())

fcs_inventory_qc <- fcs_inventory_qc %>%
  left_join(channel_patterns, by = c("channel_signature", "marker_signature"))

marker_channel_long_table <- marker_map_all_files %>%
  left_join(
    fcs_inventory_qc %>% select(file_name, cp7_channel_pattern_id),
    by = "file_name"
  ) %>%
  select(cp7_channel_pattern_id, file_name, parameter, marker_desc, range, minRange, maxRange)

channel_pattern_table <- fcs_inventory_qc %>%
  filter(read_ok) %>%
  count(cp7_channel_pattern_id, channel_signature, marker_signature, name = "n_files") %>%
  arrange(cp7_channel_pattern_id)

file_inventory <- fcs_inventory_qc %>%
  select(file_name, file_path, subject_id, panel, file_size_MB)

file_qc_channel_summary <- fcs_inventory_qc %>%
  select(
    file_name, file_path, subject_id, panel, read_ok, read_error,
    n_events, n_channels, spillover_keyword, cp7_channel_pattern_id,
    channel_order_mismatch, marker_order_mismatch, marker_set_mismatch
  )

summary_table <- tibble(
  n_fcs_files = nrow(fcs_inventory_qc),
  n_unique_subjects = n_distinct(fcs_inventory_qc$subject_id),
  n_read_ok = sum(fcs_inventory_qc$read_ok, na.rm = TRUE),
  n_read_failed = sum(!fcs_inventory_qc$read_ok, na.rm = TRUE),
  modal_n_channels = modal_n_channels,
  n_channel_patterns = n_distinct(fcs_inventory_qc$cp7_channel_pattern_id, na.rm = TRUE),
  median_events = median(fcs_inventory_qc$n_events, na.rm = TRUE),
  min_events = min(fcs_inventory_qc$n_events, na.rm = TRUE),
  max_events = max(fcs_inventory_qc$n_events, na.rm = TRUE),
  n_spillover_keywords = sum(!is.na(fcs_inventory_qc$spillover_keyword)),
  n_channel_order_mismatch = sum(fcs_inventory_qc$channel_order_mismatch, na.rm = TRUE),
  n_marker_order_mismatch = sum(fcs_inventory_qc$marker_order_mismatch, na.rm = TRUE),
  n_marker_set_mismatch = sum(fcs_inventory_qc$marker_set_mismatch, na.rm = TRUE)
)

readr::write_csv(file_inventory, file.path(out_inventory_dir, "SDY2583_CP7_file_inventory.csv"))
readr::write_csv(marker_channel_long_table, file.path(out_qc_dir, "SDY2583_CP7_marker_channel_long_table.csv"))
readr::write_csv(channel_pattern_table, file.path(out_qc_dir, "SDY2583_CP7_channel_pattern_table.csv"))
readr::write_csv(file_qc_channel_summary, file.path(out_qc_dir, "SDY2583_CP7_file_QC_channel_summary.csv"))
readr::write_csv(summary_table, file.path(out_qc_dir, "SDY2583_CP7_STEP1_QC_summary.csv"))

save(
  fcs_files, file_inventory, fcs_inventory_qc, marker_map_all_files,
  marker_channel_long_table, channel_pattern_table, reference_map, summary_table,
  file = file.path(rdata_dir, "SDY2583_CP7_STEP1_fcs_inventory_marker_QC_RECONSTRUCTED.RData")
)

print(summary_table)
message("CP7 Step 1 reconstructed inventory/QC completed.")
