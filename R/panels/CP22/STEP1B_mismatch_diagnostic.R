# ============================================================
# SDY2583 CP22
# STEP 1B: Marker/channel mismatch diagnostic for problem files
# Run after STEP 1 v2 completed successfully
# ============================================================

rm(list = ls())
source(file.path(Sys.getenv("SDY2583_REPO_ROOT", unset = "."), "R", "shared", "bootstrap.R"))

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
  library(readr)
  library(stringr)
  library(tibble)
})

analysis_dir <- sd_analysis_dir("CP22")
rdata_file <- file.path(
  analysis_dir,
  "11_RData",
  "SDY2583_CP22_STEP1_fcs_inventory_marker_QC.RData"
)

out_dir <- file.path(analysis_dir, "01_fcs_inventory")

if (!file.exists(rdata_file)) {
  stop("Step 1 RData not found: ", rdata_file)
}

load(rdata_file)

if (!exists("problem_files")) {
  stop("problem_files object not found in Step 1 RData.")
}

if (!exists("first_file")) {
  stop("first_file object not found in Step 1 RData.")
}

get_marker_map <- function(fp, label) {

  ff <- flowCore::read.FCS(
    fp,
    transformation = FALSE,
    truncate_max_range = FALSE
  )

  pp <- Biobase::pData(flowCore::parameters(ff))

  channels <- as.character(pp$name)
  markers <- as.character(pp$desc)

  idx_empty <- is.na(markers) | markers == ""
  markers[idx_empty] <- channels[idx_empty]

  tibble::tibble(
    file_label = label,
    file_path = fp,
    file_name = basename(fp),
    parameter_index = seq_along(channels),
    channel_name = channels,
    marker_desc = markers
  )
}

reference_map <- get_marker_map(first_file, "REFERENCE_FIRST_FILE")

if (nrow(problem_files) == 0) {

  cat("\nNo problem files detected in Step 1. Nothing to diagnose.\n")

} else {

  problem_paths <- unique(problem_files$file_path)

  all_problem_maps <- vector("list", length(problem_paths))

  for (i in seq_along(problem_paths)) {
    fp <- problem_paths[i]
    cat("\nReading problem file", i, "of", length(problem_paths), ":\n")
    print(fp)

    all_problem_maps[[i]] <- tryCatch(
      get_marker_map(fp, paste0("PROBLEM_", i)),
      error = function(e) {
        tibble::tibble(
          file_label = paste0("PROBLEM_", i),
          file_path = fp,
          file_name = basename(fp),
          parameter_index = NA_integer_,
          channel_name = NA_character_,
          marker_desc = NA_character_,
          error_message = as.character(e$message)
        )
      }
    )
  }

  problem_marker_maps <- dplyr::bind_rows(all_problem_maps)

  comparison_rows <- list()

  for (fp in problem_paths) {

    this_map <- problem_marker_maps[problem_marker_maps$file_path == fp, , drop = FALSE]

    if (!("parameter_index" %in% names(this_map)) || all(is.na(this_map$parameter_index))) {
      comparison_rows[[length(comparison_rows) + 1]] <- tibble::tibble(
        file_path = fp,
        file_name = basename(fp),
        parameter_index = NA_integer_,
        reference_channel = NA_character_,
        file_channel = NA_character_,
        reference_marker = NA_character_,
        file_marker = NA_character_,
        channel_same = NA,
        marker_same = NA
      )
      next
    }

    for (j in seq_len(nrow(reference_map))) {

      ref_channel <- reference_map$channel_name[j]
      ref_marker <- reference_map$marker_desc[j]

      file_channel <- NA_character_
      file_marker <- NA_character_

      if (j <= nrow(this_map)) {
        file_channel <- this_map$channel_name[j]
        file_marker <- this_map$marker_desc[j]
      }

      comparison_rows[[length(comparison_rows) + 1]] <- tibble::tibble(
        file_path = fp,
        file_name = basename(fp),
        parameter_index = j,
        reference_channel = ref_channel,
        file_channel = file_channel,
        reference_marker = ref_marker,
        file_marker = file_marker,
        channel_same = identical(ref_channel, file_channel),
        marker_same = identical(ref_marker, file_marker)
      )
    }
  }

  marker_channel_comparison <- dplyr::bind_rows(comparison_rows)

  mismatch_only <- marker_channel_comparison[
    marker_channel_comparison$channel_same == FALSE |
      marker_channel_comparison$marker_same == FALSE |
      is.na(marker_channel_comparison$channel_same) |
      is.na(marker_channel_comparison$marker_same),
    ,
    drop = FALSE
  ]

  mismatch_summary_by_file <- mismatch_only |>
    dplyr::group_by(file_name) |>
    dplyr::summarise(
      n_channel_differences = sum(channel_same == FALSE, na.rm = TRUE),
      n_marker_differences = sum(marker_same == FALSE, na.rm = TRUE),
      differing_parameters = paste(parameter_index, collapse = "; "),
      .groups = "drop"
    )

  readr::write_csv(
    reference_map,
    file.path(out_dir, "SDY2583_CP22_REFERENCE_FIRST_FILE_marker_map_STEP1B.csv")
  )

  readr::write_csv(
    problem_marker_maps,
    file.path(out_dir, "SDY2583_CP22_problem_file_marker_maps_STEP1B.csv")
  )

  readr::write_csv(
    marker_channel_comparison,
    file.path(out_dir, "SDY2583_CP22_problem_file_marker_channel_comparison_STEP1B.csv")
  )

  readr::write_csv(
    mismatch_only,
    file.path(out_dir, "SDY2583_CP22_problem_file_mismatches_ONLY_STEP1B.csv")
  )

  readr::write_csv(
    mismatch_summary_by_file,
    file.path(out_dir, "SDY2583_CP22_problem_file_mismatch_summary_STEP1B.csv")
  )

  save(
    reference_map,
    problem_marker_maps,
    marker_channel_comparison,
    mismatch_only,
    mismatch_summary_by_file,
    file = file.path(
      analysis_dir,
      "11_RData",
      "SDY2583_CP22_STEP1B_problem_file_mismatch_diagnostic.RData"
    )
  )

  cat("\n============================================================\n")
  cat("SDY2583 CP22 STEP 1B MISMATCH DIAGNOSTIC COMPLETE\n")
  cat("============================================================\n")

  cat("\nProblem file mismatch summary:\n")
  print(as.data.frame(mismatch_summary_by_file), row.names = FALSE)

  cat("\nMismatch rows only:\n")
  print(as.data.frame(mismatch_only), row.names = FALSE)

  cat("\nFiles saved in:\n")
  print(out_dir)

  cat("============================================================\n")
}
