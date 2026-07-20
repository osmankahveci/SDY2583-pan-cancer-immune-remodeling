# Run the currently reconstructed CP7 pipeline in dependency order.
# Execute from the repository root with:
#   Rscript scripts/10_run_cp7_reconstructed.R
#
# Required local configuration:
#   - SDY2583_CP7_FCS_DIR or data/raw/CP7
#   - SDY2583_METADATA_MATRIX_FILE, or a compatible CSV under the configured
#     integrated/metadata directories
#
# Optional validation archive:
#   - SDY2583_CP7_REFERENCE_DIR

repo_root <- normalizePath(
  path.expand(Sys.getenv("SDY2583_REPO_ROOT", unset = getwd())),
  mustWork = FALSE
)
Sys.setenv(SDY2583_REPO_ROOT = repo_root)

local_config <- file.path(repo_root, "config", "paths.R")
if (file.exists(local_config)) source(local_config)

steps <- c(
  "R/panels/CP7/STEP1_fcs_inventory_marker_QC_RECONSTRUCTED.R",
  "R/panels/CP7/STEP2_feature_extraction_CD8_checkpoint_RECONSTRUCTED.R",
  "R/panels/CP7/STEP3A_metadata_merge_age_QC_RECONSTRUCTED.R",
  "R/panels/CP7/STEP3B_age_sex_adjusted_statistics_RECONSTRUCTED.R",
  "R/panels/CP7/STEP4_composite_scores_RECONSTRUCTED.R"
)

for (relative_path in steps) {
  script_path <- file.path(repo_root, relative_path)
  if (!file.exists(script_path)) stop("Pipeline step is missing: ", script_path)
  cat("\n============================================================\n")
  cat("Running:", relative_path, "\n")
  cat("============================================================\n")
  source(script_path, local = new.env(parent = globalenv()))
}

validation_script <- file.path(repo_root, "scripts", "20_validate_cp7_reconstruction.R")
if (file.exists(validation_script)) {
  cat("\n============================================================\n")
  cat("Running CP7 validation\n")
  cat("============================================================\n")
  source(validation_script, local = new.env(parent = globalenv()))
}

cat("\nCP7 reconstructed pipeline completed.\n")
