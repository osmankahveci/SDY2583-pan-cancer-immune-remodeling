# CP26 Step 1 reconstructed inventory and marker QC, including the archived
# BV510-A dump-channel fallback required for seven annotation-variant files.
rm(list=ls())
source(file.path(Sys.getenv("SDY2583_REPO_ROOT",unset="."),"R","shared","reconstructed_panel_framework.R"))
rp_install_and_load(c("dplyr","readr","stringr","tibble","purrr","tidyr"),"flowCore")
source(file.path(sd_repo_root(),"R","panels","CP26","MANIFEST_RECONSTRUCTED.R"))
result<-rp_inventory("CP26",CP26_MANIFEST$expected_markers,CP26_MANIFEST$fallback_channels)
print(result$summary)
message("CP26 dump-channel resolution permits marker-name matching with BV510-A fallback.")
