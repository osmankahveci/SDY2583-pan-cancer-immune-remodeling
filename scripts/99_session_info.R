output_dir <- file.path(
  normalizePath(Sys.getenv("SDY2583_REPO_ROOT", unset = getwd()), mustWork = FALSE),
  "outputs"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

session_file <- file.path(output_dir, "session-info.txt")
writeLines(capture.output(sessionInfo()), session_file)
message("Session information written to: ", session_file)
