# Shared path helpers for the SDY2583 reanalysis scripts.
#
# Run scripts from the repository root, or set SDY2583_REPO_ROOT explicitly.
# Any path can be overridden with the environment variables documented in
# config/paths.example.R. Raw and participant-level data are never stored in
# this repository.

sd_repo_root <- function() {
  normalizePath(
    path.expand(Sys.getenv("SDY2583_REPO_ROOT", unset = getwd())),
    mustWork = FALSE
  )
}

sd_env_path <- function(variable, ...) {
  configured <- Sys.getenv(variable, unset = "")
  path <- if (nzchar(configured)) {
    path.expand(configured)
  } else {
    file.path(sd_repo_root(), ...)
  }
  normalizePath(path, mustWork = FALSE)
}

sd_analysis_dir <- function(panel) {
  sd_env_path(
    paste0("SDY2583_", toupper(panel), "_ANALYSIS_DIR"),
    "outputs",
    toupper(panel)
  )
}

sd_fcs_dir <- function(panel) {
  sd_env_path(
    paste0("SDY2583_", toupper(panel), "_FCS_DIR"),
    "data",
    "raw",
    toupper(panel)
  )
}

sd_integrated_dir <- function() {
  sd_env_path(
    "SDY2583_INTEGRATED_DIR",
    "data",
    "derived",
    "clinical_integration"
  )
}

sd_metadata_dir <- function() {
  sd_env_path(
    "SDY2583_METADATA_DIR",
    "data",
    "derived",
    "metadata"
  )
}

sd_integration_rdata_dir <- function() {
  sd_env_path(
    "SDY2583_INTEGRATION_RDATA_DIR",
    "data",
    "derived",
    "integration_rdata"
  )
}

sd_immport_download_dir <- function() {
  sd_env_path(
    "SDY2583_IMMPORT_DOWNLOAD_DIR",
    "data",
    "raw",
    "immport_download"
  )
}

sd_codebook_dir <- function() {
  sd_env_path("SDY2583_CP24_CODEBOOK_DIR", "outputs", "CP24", "codebook")
}
