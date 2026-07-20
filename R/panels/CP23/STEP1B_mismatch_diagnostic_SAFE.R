# ============================================================
# SDY2583 CP23
# STEP 1B SAFE: Mismatch diagnostic
#
# Purpose:
#   Diagnose CP23 channel-order / marker-order mismatch files.
#   This determines whether marker sets are identical and whether
#   the mismatch is a simple channel-position shift.
#
# Input:
#   outputs/CP23/11_RData/
#     SDY2583_CP23_STEP1_fcs_inventory_marker_QC.RData
#
# Output:
#   outputs/CP23/01B_mismatch_diagnostic
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
  library(tidyr)
  library(purrr)
})

filter <- dplyr::filter
select <- dplyr::select
mutate <- dplyr::mutate
arrange <- dplyr::arrange
summarise <- dplyr::summarise
group_by <- dplyr::group_by
ungroup <- dplyr::ungroup
count <- dplyr::count
bind_rows <- dplyr::bind_rows
n_distinct <- dplyr::n_distinct
n <- dplyr::n

# ------------------------------------------------------------
# 2. Paths and load Step 1
# ------------------------------------------------------------

analysis_dir <- sd_analysis_dir("CP23")
rdata_dir <- file.path(analysis_dir, "11_RData")
out_dir <- file.path(analysis_dir, "01B_mismatch_diagnostic")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

step1_rdata <- file.path(rdata_dir, "SDY2583_CP23_STEP1_fcs_inventory_marker_QC.RData")

if (!file.exists(step1_rdata)) {
  stop("Step 1 RData bulunamadı: ", step1_rdata)
}

load(step1_rdata)

# Reset paths after RData load.
analysis_dir <- sd_analysis_dir("CP23")
rdata_dir <- file.path(analysis_dir, "11_RData")
out_dir <- file.path(analysis_dir, "01B_mismatch_diagnostic")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!exists("fcs_inventory")) stop("fcs_inventory bulunamadı.")
if (!exists("reference_file")) stop("reference_file bulunamadı.")
if (!exists("reference_marker_map")) stop("reference_marker_map bulunamadı.")

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

read_marker_map_one <- function(fp) {
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
}

# ------------------------------------------------------------
# 4. Mismatch files
# ------------------------------------------------------------

mismatch_inventory <- fcs_inventory %>%
  filter(
    read_ok == TRUE,
    channel_order_mismatch_file == TRUE |
      marker_order_mismatch_file == TRUE |
      marker_set_mismatch_file == TRUE
  ) %>%
  arrange(file_name)

if (nrow(mismatch_inventory) == 0) {
  mismatch_diagnostic_summary <- tibble(
    note = "No mismatch files detected in Step 1."
  )

  write_csv(mismatch_diagnostic_summary, file.path(out_dir, "SDY2583_CP23_mismatch_diagnostic_summary_STEP1B.csv"))

  cat("\n============================================================\n")
  cat("SDY2583 CP23 STEP 1B COMPLETE: NO MISMATCH FILES\n")
  cat("============================================================\n")
  print(as.data.frame(mismatch_diagnostic_summary), row.names = FALSE)

  save(
    mismatch_diagnostic_summary,
    file = file.path(rdata_dir, "SDY2583_CP23_STEP1B_mismatch_diagnostic.RData")
  )

  quit(save = "no")
}

reference_map <- read_marker_map_one(reference_file) %>%
  select(parameter_index, reference_channel = channel_name, reference_marker = marker_desc)

mismatch_maps <- bind_rows(lapply(mismatch_inventory$file_path, read_marker_map_one))

side_by_side_marker_maps <- mismatch_maps %>%
  select(file_name, subject_id, file_path, parameter_index, mismatch_channel = channel_name, mismatch_marker = marker_desc) %>%
  left_join(reference_map, by = "parameter_index") %>%
  mutate(
    same_channel_at_position = mismatch_channel == reference_channel,
    same_marker_at_position = mismatch_marker == reference_marker
  ) %>%
  select(
    file_name,
    subject_id,
    parameter_index,
    reference_channel,
    reference_marker,
    mismatch_channel,
    mismatch_marker,
    same_channel_at_position,
    same_marker_at_position
  ) %>%
  arrange(file_name, parameter_index)

# File-level diagnostics
mismatch_diagnostic_summary <- side_by_side_marker_maps %>%
  group_by(file_name, subject_id) %>%
  summarise(
    n_parameters = n(),
    marker_set_identical_to_reference =
      setequal(mismatch_marker, reference_marker),
    channel_set_identical_to_reference =
      setequal(mismatch_channel, reference_channel),
    n_same_channel_positions = sum(same_channel_at_position, na.rm = TRUE),
    n_same_marker_positions = sum(same_marker_at_position, na.rm = TRUE),
    first_mismatch_position = suppressWarnings(min(parameter_index[!same_channel_at_position | !same_marker_at_position], na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    first_mismatch_position = ifelse(is.infinite(first_mismatch_position), NA_real_, first_mismatch_position)
  ) %>%
  left_join(
    mismatch_inventory %>%
      select(
        file_name,
        n_events,
        n_channels,
        channel_order_mismatch_file,
        marker_order_mismatch_file,
        marker_set_mismatch_file,
        has_spillover_keyword,
        spillover_keyword_name
      ),
    by = "file_name"
  ) %>%
  arrange(file_name)

# Marker-to-channel map ignoring parameter order:
reference_marker_to_channel <- reference_map %>%
  select(reference_marker, reference_channel) %>%
  distinct()

mismatch_marker_to_channel <- mismatch_maps %>%
  select(file_name, subject_id, marker_desc, channel_name) %>%
  distinct()

marker_channel_identity_by_marker <- mismatch_marker_to_channel %>%
  left_join(reference_marker_to_channel, by = c("marker_desc" = "reference_marker")) %>%
  mutate(
    same_channel_for_marker_as_reference = channel_name == reference_channel
  ) %>%
  arrange(file_name, marker_desc)

marker_channel_identity_summary <- marker_channel_identity_by_marker %>%
  group_by(file_name, subject_id) %>%
  summarise(
    n_markers = n(),
    n_markers_same_channel_as_reference = sum(same_channel_for_marker_as_reference, na.rm = TRUE),
    all_markers_same_channel_as_reference = all(same_channel_for_marker_as_reference, na.rm = TRUE),
    markers_changed_channel = paste(marker_desc[!same_channel_for_marker_as_reference], collapse = "; "),
    .groups = "drop"
  )

# ------------------------------------------------------------
# 5. Save
# ------------------------------------------------------------

write_csv(mismatch_diagnostic_summary, file.path(out_dir, "SDY2583_CP23_mismatch_diagnostic_summary_STEP1B.csv"))
write_csv(side_by_side_marker_maps, file.path(out_dir, "SDY2583_CP23_side_by_side_marker_maps_STEP1B.csv"))
write_csv(marker_channel_identity_by_marker, file.path(out_dir, "SDY2583_CP23_marker_channel_identity_by_marker_STEP1B.csv"))
write_csv(marker_channel_identity_summary, file.path(out_dir, "SDY2583_CP23_marker_channel_identity_summary_STEP1B.csv"))

save(
  mismatch_inventory,
  mismatch_diagnostic_summary,
  side_by_side_marker_maps,
  marker_channel_identity_by_marker,
  marker_channel_identity_summary,
  file = file.path(rdata_dir, "SDY2583_CP23_STEP1B_mismatch_diagnostic.RData")
)

# ------------------------------------------------------------
# 6. Console output
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SDY2583 CP23 STEP 1B COMPLETE: MISMATCH DIAGNOSTIC\n")
cat("============================================================\n")

cat("\nReference file:\n")
print(basename(reference_file))

cat("\nMismatch diagnostic summary:\n")
print(as.data.frame(mismatch_diagnostic_summary), row.names = FALSE)

cat("\nSide-by-side marker maps:\n")
print(as.data.frame(side_by_side_marker_maps), row.names = FALSE)

cat("\nMarker-channel identity summary:\n")
print(as.data.frame(marker_channel_identity_summary), row.names = FALSE)

cat("\nOutputs saved in:\n")
print(out_dir)

cat("============================================================\n")
