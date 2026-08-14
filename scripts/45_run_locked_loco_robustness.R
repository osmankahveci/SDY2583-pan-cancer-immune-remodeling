#!/usr/bin/env Rscript

# Runs the training-frozen locked-score construction and then performs the
# leave-one-cancer-type-out robustness analysis in the same R session. Only
# aggregate LOCO outputs are written; participant-level matrices stay in memory.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else "scripts/45_run_locked_loco_robustness.R"
repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)

registry_dir <- file.path(repo_root, "config", "locked_validation")
component_files <- sort(list.files(
  registry_dir,
  pattern = "^score_component_registry_CP[0-9]+\\.csv$",
  full.names = TRUE
))
if (length(component_files) != 10L) {
  stop("Expected 10 panel-specific component registries; found ", length(component_files))
}

registry_lines <- lapply(component_files, readLines, warn = FALSE)
combined_lines <- c(
  registry_lines[[1]],
  unlist(lapply(registry_lines[-1], function(x) x[-1]), use.names = FALSE)
)

tmp_root <- tempfile("sdy2583_locked_loco_registry_")
tmp_registry_dir <- file.path(tmp_root, "config", "locked_validation")
dir.create(tmp_registry_dir, recursive = TRUE, showWarnings = FALSE)
writeLines(combined_lines, file.path(tmp_registry_dir, "score_component_registry.csv"), useBytes = TRUE)

summary_registry <- file.path(registry_dir, "score_summary_registry.csv")
if (!file.exists(summary_registry)) stop("Missing score summary registry: ", summary_registry)
file.copy(summary_registry, tmp_registry_dir, overwrite = TRUE)

if (!nzchar(Sys.getenv("SDY2583_LOCKED_OUTPUT_DIR", ""))) {
  Sys.setenv(SDY2583_LOCKED_OUTPUT_DIR = file.path(repo_root, "results", "local", "locked_internal_validation"))
}
if (!nzchar(Sys.getenv("SDY2583_LOCKED_LOCO_OUTPUT_DIR", ""))) {
  Sys.setenv(SDY2583_LOCKED_LOCO_OUTPUT_DIR = file.path(repo_root, "results", "local", "locked_loco_robustness"))
}

# The locked engine uses SDY2583_REPO_ROOT only to locate its ephemeral combined
# registry. Preserve the actual public repository root separately for LOCO.
Sys.setenv(
  SDY2583_REPO_ROOT = tmp_root,
  SDY2583_PUBLIC_REPO_ROOT = repo_root
)

source(file.path(repo_root, "R", "integration", "locked_internal_validation.R"), chdir = FALSE)
source(file.path(repo_root, "R", "integration", "locked_loco_robustness.R"), chdir = FALSE)

unlink(tmp_root, recursive = TRUE, force = TRUE)
