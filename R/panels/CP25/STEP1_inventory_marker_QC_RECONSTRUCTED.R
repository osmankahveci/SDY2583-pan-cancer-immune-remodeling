# CP25 Step 1 reconstructed inventory and marker QC.
rm(list=ls())
source(file.path(Sys.getenv("SDY2583_REPO_ROOT",unset="."),"R","shared","reconstructed_panel_framework.R"))
rp_install_and_load(c("dplyr","readr","stringr","tibble","purrr","tidyr"),"flowCore")
source(file.path(sd_repo_root(),"R","panels","CP25","MANIFEST_RECONSTRUCTED.R"))
result<-rp_inventory("CP25",CP25_MANIFEST$expected_markers,list("Viability_CD8_CD13_CD19_TCRgd"="BV510-A")); print(result$summary)
