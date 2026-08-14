#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else "scripts/44_run_locked_internal_validation.R"
repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
source(file.path(repo_root, "R", "integration", "locked_internal_validation.R"), chdir = FALSE)
