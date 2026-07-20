# Build a cross-panel subject-level score matrix from locally generated panel
# outputs. Participant-level integrated data remain excluded from Git.
rm(list = ls())
source(file.path(Sys.getenv("SDY2583_REPO_ROOT", unset = "."), "R", "shared", "bootstrap.R"))

for (p in c("dplyr", "readr", "stringr", "tibble", "purrr", "tidyr")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(stringr)
  library(tibble); library(purrr); library(tidyr)
})

panels <- c("CP7", "CP8", "CP10", "CP16", "CP22", "CP23", "CP24", "CP25", "CP26", "CP28")
score_pattern <- "^CP[0-9]+_.*_score$|^CP[0-9]+_integrated_"

first_non_missing <- function(x) {
  y <- x[!is.na(x)]
  if (is.character(y)) y <- y[nzchar(trimws(y))]
  if (length(y) == 0L) return(x[NA_integer_][1])
  y[1]
}

find_score_file <- function(panel) {
  files <- list.files(
    sd_analysis_dir(panel), "\\.csv$", full.names = TRUE,
    recursive = TRUE, ignore.case = TRUE
  )
  if (length(files) == 0L) return(NA_character_)

  base <- basename(files)
  rank <- integer(length(files))
  rank <- rank + 5000L * grepl("analysis_data_with_composite_scores", base, ignore.case = TRUE)
  rank <- rank + 4500L * grepl("analysis_data_with_composite", base, ignore.case = TRUE)
  rank <- rank + 4000L * grepl("score_dataset", base, ignore.case = TRUE)
  rank <- rank + 3000L * grepl("with_composite_scores", base, ignore.case = TRUE)
  rank <- rank - 5000L * grepl(
    "statistics|results?|summary|feature_sets?|feature_map|manifest|validation",
    base, ignore.case = TRUE
  )

  ordered <- files[order(rank, file.info(files)$mtime, decreasing = TRUE)]
  for (path in ordered) {
    header <- tryCatch(
      readr::read_csv(path, n_max = 2, show_col_types = FALSE, progress = FALSE),
      error = function(e) NULL
    )
    if (is.null(header)) next
    has_subject <- any(c("subject_id", "subject_accession") %in% names(header))
    has_score <- any(grepl(score_pattern, names(header)))
    if (has_subject && has_score) return(path)
  }
  NA_character_
}

read_scores <- function(panel) {
  path <- find_score_file(panel)
  if (is.na(path) || !file.exists(path)) return(NULL)
  x <- tryCatch(
    readr::read_csv(path, show_col_types = FALSE, progress = FALSE),
    error = function(e) NULL
  )
  if (is.null(x)) return(NULL)

  subject_candidates <- c("subject_id", "subject_accession")
  subject_col <- subject_candidates[subject_candidates %in% names(x)]
  if (length(subject_col) == 0L) return(NULL)
  subject_col <- subject_col[1]
  score_cols <- grep(score_pattern, names(x), value = TRUE)
  if (length(score_cols) == 0L) return(NULL)

  x %>%
    transmute(
      subject_id = as.character(.data[[subject_col]]),
      across(all_of(score_cols), ~suppressWarnings(as.numeric(.x)))
    ) %>%
    filter(!is.na(subject_id), subject_id != "") %>%
    group_by(subject_id) %>%
    summarise(across(everything(), first_non_missing), .groups = "drop")
}

score_tables <- setNames(lapply(panels, read_scores), panels)
available <- names(score_tables)[!vapply(score_tables, is.null, logical(1))]
if (length(available) == 0L) {
  stop("No panel composite-score datasets were found under the configured outputs directories.")
}

wide <- Reduce(
  function(x, y) full_join(x, y, by = "subject_id"),
  score_tables[available]
)

metadata_file <- path.expand(Sys.getenv(
  "SDY2583_METADATA_MATRIX_FILE",
  unset = file.path(sd_metadata_dir(), "SDY2583_subject_FCS_metadata_matrix_RECONSTRUCTED.csv")
))
if (file.exists(metadata_file)) {
  metadata <- readr::read_csv(metadata_file, show_col_types = FALSE, progress = FALSE)
  metadata <- metadata %>%
    mutate(subject_id = as.character(subject_id)) %>%
    filter(!is.na(subject_id), subject_id != "") %>%
    group_by(subject_id) %>%
    summarise(across(everything(), first_non_missing), .groups = "drop")
  wide <- metadata %>% right_join(wide, by = "subject_id")
}

out_dir <- sd_integrated_dir()
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
manifest <- tibble(
  panel = panels,
  score_file = vapply(panels, find_score_file, character(1)),
  included = panels %in% available,
  n_score_columns = vapply(
    score_tables,
    function(x) if (is.null(x)) 0L else ncol(x) - 1L,
    integer(1)
  )
)

readr::write_csv(
  wide,
  file.path(out_dir, "SDY2583_integrated_clinical_immune_score_matrix_ALL_PANELS_RECONSTRUCTED.csv")
)
readr::write_csv(
  manifest,
  file.path(out_dir, "SDY2583_cross_panel_score_integration_manifest_RECONSTRUCTED.csv")
)
cat("Integrated ", length(available), " panels and ", nrow(wide), " subjects.\n", sep = "")
