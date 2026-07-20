# CP8 Step 1 reconstructed inventory and marker/channel QC.
rm(list = ls())
source(file.path(Sys.getenv("SDY2583_REPO_ROOT", unset = "."), "R", "shared", "reconstructed_panel_framework.R"))
rp_install_and_load(c("dplyr", "readr", "stringr", "tibble", "purrr", "tidyr"), "flowCore")
source(file.path(sd_repo_root(), "R", "panels", "CP8", "MANIFEST_RECONSTRUCTED.R"))
result <- rp_inventory("CP8", CP8_MANIFEST$expected_markers)
print(result$summary)
message("Archived benchmark: 850/850 files; three channel signatures (841, 2, and 7 files).")
