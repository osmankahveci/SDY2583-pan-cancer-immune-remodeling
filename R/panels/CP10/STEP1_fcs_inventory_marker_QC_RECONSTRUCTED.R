# CP10 Step 1 reconstructed inventory and marker/channel QC.
# Produces the exact RData filename expected by the recovered CP10 Step 2.
rm(list=ls())
source(file.path(Sys.getenv("SDY2583_REPO_ROOT",unset="."),"R","shared","reconstructed_panel_framework.R"))
rp_install_and_load(c("dplyr","readr","stringr","tibble","purrr","tidyr"),"flowCore")
expected<-c("Viability","CD45","CD3","CD19","CD56","CD14","HLA-DR","CD11c","CD13","CD66b","CCR3","CD123")
result<-rp_inventory("CP10",expected)
fcs_files<-result$fcs_files; fcs_inventory<-result$inventory; reference_marker_map<-result$reference_map; fcs_inventory_summary<-result$summary
analysis_dir<-sd_analysis_dir("CP10"); rdata_dir<-file.path(analysis_dir,"11_RData"); dir.create(rdata_dir,recursive=TRUE,showWarnings=FALSE)
save(fcs_files,fcs_inventory,reference_marker_map,fcs_inventory_summary,file=file.path(rdata_dir,"SDY2583_CP10_STEP1_fcs_inventory_marker_QC.RData"))
print(fcs_inventory_summary)
