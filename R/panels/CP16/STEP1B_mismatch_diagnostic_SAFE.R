# ============================================================
# SDY2583 CP16
# STEP 1B SAFE: Mismatch diagnostic
#
# Purpose:
#   - Inspect CP16 files with channel/marker-order mismatch
#   - Determine whether the marker set is identical but shifted
#   - Print side-by-side reference vs mismatch marker maps
#
# Output:
#   outputs/CP16/01B_mismatch_diagnostic
# ============================================================

rm(list = ls())
source(file.path(Sys.getenv("SDY2583_REPO_ROOT", unset = "."), "R", "shared", "bootstrap.R"))

# ------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------

cran_pkgs <- c("dplyr", "readr", "stringr", "tibble")

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
})

# Namespace safety
filter <- dplyr::filter
select <- dplyr::select
mutate <- dplyr::mutate
arrange <- dplyr::arrange
bind_rows <- dplyr::bind_rows
left_join <- dplyr::left_join
n_distinct <- dplyr::n_distinct

# ------------------------------------------------------------
# 2. Paths and load Step 1
# ------------------------------------------------------------

analysis_dir <- sd_analysis_dir("CP16")
rdata_dir <- file.path(analysis_dir, "11_RData")
out_dir <- file.path(analysis_dir, "01B_mismatch_diagnostic")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

step1_rdata <- file.path(rdata_dir, "SDY2583_CP16_STEP1_fcs_inventory_marker_QC.RData")

if (!file.exists(step1_rdata)) {
  stop("Step 1 RData bulunamadı: ", step1_rdata)
}

load(step1_rdata)

if (!exists("fcs_inventory")) stop("fcs_inventory bulunamadı.")
if (!exists("reference_marker_map")) stop("reference_marker_map bulunamadı.")
if (!exists("reference_file")) stop("reference_file bulunamadı.")

# ------------------------------------------------------------
# 3. Helpers
# ------------------------------------------------------------

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

read_marker_map_file <- function(fp) {
  ff <- flowCore::read.FCS(
    fp,
    transformation = FALSE,
    truncate_max_range = FALSE
  )
  get_marker_map(ff)
}

# ------------------------------------------------------------
# 4. Identify mismatch files
# ------------------------------------------------------------

mismatch_df <- fcs_inventory %>%
  filter(
    read_ok == TRUE,
    channel_order_mismatch_file == TRUE |
      marker_order_mismatch_file == TRUE |
      marker_set_mismatch_file == TRUE
  ) %>%
  arrange(file_name)

if (nrow(mismatch_df) == 0) {
  cat("\nNo mismatch files found. Step 1B diagnostic not required.\n")
} else {

  all_diagnostics <- list()
  all_maps <- list()

  ref_map <- reference_marker_map %>%
    mutate(reference_pair = paste(channel_name, marker_desc, sep = " :: "))

  for (i in seq_len(nrow(mismatch_df))) {

    fp <- mismatch_df$file_path[i]
    fn <- mismatch_df$file_name[i]

    mm <- read_marker_map_file(fp) %>%
      mutate(mismatch_pair = paste(channel_name, marker_desc, sep = " :: "))

    side_by_side <- ref_map %>%
      select(parameter_index, reference_channel = channel_name, reference_marker = marker_desc, reference_pair) %>%
      full_join(
        mm %>%
          select(parameter_index, mismatch_channel = channel_name, mismatch_marker = marker_desc, mismatch_pair),
        by = "parameter_index"
      ) %>%
      mutate(
        file_name = fn,
        same_channel_at_position = reference_channel == mismatch_channel,
        same_marker_at_position = reference_marker == mismatch_marker
      ) %>%
      select(
        file_name,
        parameter_index,
        reference_channel,
        reference_marker,
        mismatch_channel,
        mismatch_marker,
        same_channel_at_position,
        same_marker_at_position
      )

    diagnostic <- tibble(
      file_name = fn,
      subject_id = mismatch_df$subject_id[i],
      n_events = mismatch_df$n_events[i],
      n_channels = mismatch_df$n_channels[i],
      marker_set_identical_to_reference =
        identical(sort(unique(ref_map$marker_desc)), sort(unique(mm$marker_desc))),
      channel_set_identical_to_reference =
        identical(sort(unique(ref_map$channel_name)), sort(unique(mm$channel_name))),
      n_same_channel_positions = sum(side_by_side$same_channel_at_position, na.rm = TRUE),
      n_same_marker_positions = sum(side_by_side$same_marker_at_position, na.rm = TRUE),
      n_parameters = nrow(side_by_side),
      first_mismatch_position = min(side_by_side$parameter_index[
        side_by_side$same_channel_at_position == FALSE |
          side_by_side$same_marker_at_position == FALSE
      ], na.rm = TRUE)
    )

    all_diagnostics[[i]] <- diagnostic
    all_maps[[i]] <- side_by_side
  }

  mismatch_diagnostic_summary <- bind_rows(all_diagnostics)
  mismatch_side_by_side_maps <- bind_rows(all_maps)

  write_csv(
    mismatch_diagnostic_summary,
    file.path(out_dir, "SDY2583_CP16_mismatch_diagnostic_summary_STEP1B.csv")
  )

  write_csv(
    mismatch_side_by_side_maps,
    file.path(out_dir, "SDY2583_CP16_mismatch_side_by_side_marker_maps_STEP1B.csv")
  )

  save(
    mismatch_diagnostic_summary,
    mismatch_side_by_side_maps,
    mismatch_df,
    reference_marker_map,
    reference_file,
    file = file.path(rdata_dir, "SDY2583_CP16_STEP1B_mismatch_diagnostic.RData")
  )

  cat("\n============================================================\n")
  cat("SDY2583 CP16 STEP 1B COMPLETE: MISMATCH DIAGNOSTIC\n")
  cat("============================================================\n")

  cat("\nReference file:\n")
  print(basename(reference_file))

  cat("\nMismatch diagnostic summary:\n")
  print(as.data.frame(mismatch_diagnostic_summary), row.names = FALSE)

  cat("\nSide-by-side marker maps:\n")
  print(as.data.frame(mismatch_side_by_side_maps), row.names = FALSE)

  cat("\nOutputs saved in:\n")
  print(out_dir)

  cat("============================================================\n")
}
