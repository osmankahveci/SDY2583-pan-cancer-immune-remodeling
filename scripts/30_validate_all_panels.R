# Aggregate panel-specific validation reports into one release gate.
root <- normalizePath(
  path.expand(Sys.getenv("SDY2583_REPO_ROOT", unset = getwd())),
  mustWork = FALSE
)

for (p in c("dplyr", "readr", "tibble", "purrr", "stringr")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tibble)
  library(purrr); library(stringr)
})

panels <- c("CP7", "CP8", "CP10", "CP16", "CP22", "CP23", "CP24", "CP25", "CP26", "CP28")

find_report <- function(panel) {
  directory <- file.path(root, "outputs", panel, "validation")
  if (!dir.exists(directory)) return(NA_character_)
  files <- list.files(
    directory,
    pattern = "(reconstruction_validation_report|recovered_pipeline_validation_report)\\.csv$",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )
  if (length(files) == 0L) return(NA_character_)
  files[order(file.info(files)$mtime, decreasing = TRUE)][1]
}

reports <- setNames(vapply(panels, find_report, character(1)), panels)
summary_table <- purrr::imap_dfr(reports, function(path, panel) {
  if (is.na(path) || !file.exists(path)) {
    return(tibble(
      panel = panel, report_file = NA_character_, n_checks = 0L,
      n_pass = 0L, n_fail = 1L, all_pass = FALSE
    ))
  }
  x <- readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
  pass_candidates <- c("pass", "PASS")
  pass_col <- pass_candidates[pass_candidates %in% names(x)]
  if (length(pass_col) == 0L) {
    return(tibble(
      panel = panel, report_file = path, n_checks = nrow(x),
      n_pass = 0L, n_fail = max(1L, nrow(x)), all_pass = FALSE
    ))
  }
  passed <- as.logical(x[[pass_col[1]]])
  tibble(
    panel = panel,
    report_file = path,
    n_checks = length(passed),
    n_pass = sum(passed %in% TRUE),
    n_fail = sum(!passed | is.na(passed)),
    all_pass = length(passed) > 0L && all(passed %in% TRUE)
  )
})

metadata_file <- file.path(
  root, "data", "derived", "metadata",
  "SDY2583_subject_FCS_metadata_matrix_RECONSTRUCTED.csv"
)
integration_file <- file.path(
  root, "data", "derived", "clinical_integration",
  "SDY2583_integrated_clinical_immune_score_matrix_ALL_PANELS_RECONSTRUCTED.csv"
)
session_file <- file.path(root, "outputs", "session-info.txt")

global_paths <- c(metadata_file, integration_file, session_file)
global_checks <- tibble(
  panel = c("GLOBAL_METADATA", "GLOBAL_SCORE_MATRIX", "SESSION_INFO"),
  report_file = global_paths,
  n_checks = 1L,
  n_pass = as.integer(file.exists(global_paths)),
  n_fail = 1L - as.integer(file.exists(global_paths)),
  all_pass = file.exists(global_paths)
)

summary_table <- bind_rows(summary_table, global_checks)
out_dir <- file.path(root, "outputs", "validation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
readr::write_csv(
  summary_table,
  file.path(out_dir, "SDY2583_all_panel_validation_summary.csv")
)
print(summary_table, n = nrow(summary_table))
if (any(!summary_table$all_pass)) {
  stop(
    "All-panel validation gate failed. Review ",
    file.path(out_dir, "SDY2583_all_panel_validation_summary.csv")
  )
}
