# Parse-only and repository-structure audit that does not require SDY2583 data
# or analysis packages. Intended for local use and GitHub Actions.
root <- normalizePath(
  path.expand(Sys.getenv("SDY2583_REPO_ROOT", unset = getwd())),
  mustWork = TRUE
)
Sys.setenv(SDY2583_REPO_ROOT = root)

r_files <- list.files(
  root,
  pattern = "\\.[Rr]$",
  full.names = TRUE,
  recursive = TRUE
)
r_files <- r_files[!grepl("/(data|outputs|renv/library|packrat/lib)/", r_files)]
if (length(r_files) == 0L) stop("No R source files were found.")

parse_results <- lapply(r_files, function(path) {
  error <- tryCatch({
    parse(file = path, keep.source = TRUE)
    NA_character_
  }, error = function(e) conditionMessage(e))
  data.frame(
    file = substring(path, nchar(root) + 2L),
    parse_ok = is.na(error),
    parse_error = error,
    stringsAsFactors = FALSE
  )
})
parse_results <- do.call(rbind, parse_results)

read_text <- function(path) {
  paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
}
strip_comment_lines <- function(text) {
  lines <- strsplit(text, "\n", fixed = TRUE)[[1]]
  lines <- lines[!grepl("^\\s*#", lines)]
  paste(lines, collapse = "\n")
}
source_text <- setNames(vapply(r_files, read_text, character(1)), r_files)
# Do not scan the audit's own rule declarations, and ignore comment-only legacy
# path examples. All files are still parsed above.
scan_text <- source_text[basename(names(source_text)) != "02_static_source_audit.R"]
scan_text <- vapply(scan_text, strip_comment_lines, character(1))

forbidden_patterns <- c(
  hard_coded_macos_user = "/Users/[A-Za-z0-9._-]+/",
  hard_coded_home_desktop = "~/Desktop|~/Downloads",
  interactive_file_choose = "file\\.choose\\s*\\(|choose\\.files\\s*\\("
)
forbidden_parts <- lapply(names(forbidden_patterns), function(name) {
  pattern <- forbidden_patterns[[name]]
  hit <- vapply(scan_text, grepl, logical(1), pattern = pattern, perl = TRUE)
  if (!any(hit)) {
    return(data.frame(
      rule = character(0), file = character(0), pass = logical(0),
      stringsAsFactors = FALSE
    ))
  }
  data.frame(
    rule = rep(name, sum(hit)),
    file = substring(names(scan_text)[hit], nchar(root) + 2L),
    pass = rep(FALSE, sum(hit)),
    stringsAsFactors = FALSE
  )
})
forbidden_results <- do.call(rbind, forbidden_parts)

runner_files <- list.files(
  file.path(root, "scripts"),
  pattern = "^([0-9]+_run_|01_run_all_panels).*\\.[Rr]$",
  full.names = TRUE
)
runner_reference_results <- do.call(rbind, lapply(runner_files, function(path) {
  text <- read_text(path)
  matches <- regmatches(
    text,
    gregexpr(
      "(?:R/panels|R/integration|scripts)/[A-Za-z0-9_./-]+\\.[Rr]",
      text,
      perl = TRUE
    )
  )[[1]]
  matches <- unique(matches[matches != ""])
  if (length(matches) == 0L) {
    return(data.frame(
      runner = substring(path, nchar(root) + 2L),
      reference = NA_character_,
      exists = TRUE,
      stringsAsFactors = FALSE
    ))
  }
  data.frame(
    runner = rep(substring(path, nchar(root) + 2L), length(matches)),
    reference = matches,
    exists = file.exists(file.path(root, matches)),
    stringsAsFactors = FALSE
  )
}))

required_files <- c(
  "README.md", "PROVENANCE.md", "DATA_ACCESS.md", "CITATION.cff", "LICENSE",
  "R/shared/bootstrap.R", "R/shared/reconstructed_panel_framework.R",
  "R/integration/STEP0_build_subject_metadata_matrix_RECONSTRUCTED.R",
  "R/integration/STEP1_build_cross_panel_score_matrix_RECONSTRUCTED.R",
  "scripts/01_run_all_panels.R", "scripts/30_validate_all_panels.R",
  "scripts/99_session_info.R"
)
required_results <- data.frame(
  file = required_files,
  exists = file.exists(file.path(root, required_files)),
  stringsAsFactors = FALSE
)

panels <- c("CP7", "CP8", "CP10", "CP16", "CP22", "CP23", "CP24", "CP25", "CP26", "CP28")
panel_results <- do.call(rbind, lapply(panels, function(panel) {
  directory <- file.path(root, "R", "panels", panel)
  files <- if (dir.exists(directory)) {
    list.files(directory, "\\.[Rr]$", full.names = FALSE)
  } else {
    character()
  }
  data.frame(
    panel = panel,
    directory_exists = dir.exists(directory),
    n_r_files = length(files),
    has_source = length(files) > 0L,
    stringsAsFactors = FALSE
  )
}))

# Binary/raw data must never be committed outside source/config paths.
binary_forbidden_extensions <- "\\.(fcs|rdata|rds|zip|xlsx|xls)$"
tracked_like_files <- list.files(root, full.names = TRUE, recursive = TRUE, all.files = TRUE)
tracked_like_files <- tracked_like_files[!grepl("/\\.git/", tracked_like_files)]
forbidden_binary_data <- tracked_like_files[
  grepl(binary_forbidden_extensions, tracked_like_files, ignore.case = TRUE) &
    !grepl("/(config|scripts|R)/", tracked_like_files)
]

# Aggregate CSV/TSV result tables are allowed only in the curated final-results
# directory or for the three legacy aggregate bootstrap summaries. Their header
# must not expose participant-level identifiers or matched subject IDs.
tabular_files <- tracked_like_files[
  grepl("\\.(csv|tsv)$", tracked_like_files, ignore.case = TRUE) &
    !grepl("/(config|scripts|R)/", tracked_like_files)
]
relative_tabular <- substring(tabular_files, nchar(root) + 2L)
legacy_aggregate <- grepl(
  "^results/bootstrap_(component_stability_summary|panel_component_summary|principal_score_effect_stability)\\.csv$",
  relative_tabular
)
curated_aggregate <- grepl("^results/final/[^/]+\\.(csv|tsv)$", relative_tabular)
allowed_tabular <- curated_aggregate | legacy_aggregate
forbidden_tabular_location <- tabular_files[!allowed_tabular]

prohibited_header_pattern <- paste(
  c(
    "(^|,)(subject_id|subject_accession|participant_id|master_id)(,|$)",
    "(^|,)(cancer_id|healthy_id|matched_pair_id|pair_id)(,|$)",
    "(^|,)(individual_pc1|individual_pc2|pca_coordinate)(,|$)"
  ),
  collapse = "|"
)
aggregate_header_results <- do.call(rbind, lapply(tabular_files[allowed_tabular], function(path) {
  first_line <- readLines(path, n = 1L, warn = FALSE, encoding = "UTF-8")
  if (length(first_line) == 0L) first_line <- ""
  data.frame(
    file = substring(path, nchar(root) + 2L),
    allowed_location = TRUE,
    prohibited_identifier_header = grepl(
      prohibited_header_pattern,
      tolower(first_line),
      perl = TRUE
    ),
    stringsAsFactors = FALSE
  )
}))
if (is.null(aggregate_header_results)) {
  aggregate_header_results <- data.frame(
    file = character(), allowed_location = logical(),
    prohibited_identifier_header = logical(), stringsAsFactors = FALSE
  )
}

out_dir <- file.path(root, "outputs", "static-audit")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(parse_results, file.path(out_dir, "parse-results.csv"), row.names = FALSE)
write.csv(forbidden_results, file.path(out_dir, "forbidden-pattern-results.csv"), row.names = FALSE)
write.csv(
  runner_reference_results,
  file.path(out_dir, "runner-reference-results.csv"),
  row.names = FALSE
)
write.csv(
  required_results,
  file.path(out_dir, "required-file-results.csv"),
  row.names = FALSE
)
write.csv(
  panel_results,
  file.path(out_dir, "panel-source-results.csv"),
  row.names = FALSE
)
write.csv(
  aggregate_header_results,
  file.path(out_dir, "aggregate-result-header-audit.csv"),
  row.names = FALSE
)

failures <- character()
if (any(!parse_results$parse_ok)) {
  failures <- c(
    failures,
    paste0(
      "R parse failures: ",
      paste(parse_results$file[!parse_results$parse_ok], collapse = ", ")
    )
  )
}
if (nrow(forbidden_results) > 0L) {
  failures <- c(
    failures,
    paste0(
      "Forbidden path/interactive patterns: ",
      paste(forbidden_results$file, collapse = ", ")
    )
  )
}
if (any(!runner_reference_results$exists)) {
  failures <- c(
    failures,
    paste0(
      "Missing runner references: ",
      paste(
        runner_reference_results$reference[!runner_reference_results$exists],
        collapse = ", "
      )
    )
  )
}
if (any(!required_results$exists)) {
  failures <- c(
    failures,
    paste0(
      "Missing required files: ",
      paste(required_results$file[!required_results$exists], collapse = ", ")
    )
  )
}
if (any(!panel_results$has_source)) {
  failures <- c(
    failures,
    paste0(
      "Panels without R source: ",
      paste(panel_results$panel[!panel_results$has_source], collapse = ", ")
    )
  )
}
if (length(forbidden_binary_data) > 0L) {
  failures <- c(
    failures,
    paste0(
      "Potential raw/binary data committed: ",
      paste(substring(forbidden_binary_data, nchar(root) + 2L), collapse = ", ")
    )
  )
}
if (length(forbidden_tabular_location) > 0L) {
  failures <- c(
    failures,
    paste0(
      "Tabular data committed outside the aggregate allowlist: ",
      paste(substring(forbidden_tabular_location, nchar(root) + 2L), collapse = ", ")
    )
  )
}
if (nrow(aggregate_header_results) > 0L &&
    any(aggregate_header_results$prohibited_identifier_header)) {
  failures <- c(
    failures,
    paste0(
      "Aggregate result tables contain prohibited participant identifiers: ",
      paste(
        aggregate_header_results$file[
          aggregate_header_results$prohibited_identifier_header
        ],
        collapse = ", "
      )
    )
  )
}

cat("Parsed ", nrow(parse_results), " R files.\n", sep = "")
cat("Checked ", nrow(runner_reference_results), " runner references.\n", sep = "")
cat("Audited ", nrow(aggregate_header_results), " aggregate result tables.\n", sep = "")
cat("All ten panel source directories contain R code.\n")
if (length(failures) > 0L) stop(paste(failures, collapse = "\n"))
cat("Static source audit passed.\n")
