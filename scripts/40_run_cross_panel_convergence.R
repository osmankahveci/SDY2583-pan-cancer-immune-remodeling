# Run the SDY2583 cross-panel convergence analysis from the repository root.
repo_root <- normalizePath(
  Sys.getenv("SDY2583_REPO_ROOT", unset = getwd()),
  mustWork = FALSE
)
Sys.setenv(SDY2583_REPO_ROOT = repo_root)

local_config <- file.path(repo_root, "config", "paths.R")
if (file.exists(local_config)) source(local_config)

source(file.path(
  repo_root,
  "R", "integration", "STEP2_cross_panel_convergence_analysis.R"
))
