# CP8 Step 7B reconstructed representative real-FCS gating atlas.
rm(list=ls())
source(file.path(Sys.getenv("SDY2583_REPO_ROOT",unset="."),"R","shared","figure_framework.R"))
rp_install_and_load(c("dplyr","readr","tibble","ggplot2","patchwork","stringr"),"flowCore")
source(file.path(sd_repo_root(),"R","panels","CP8","MANIFEST_RECONSTRUCTED.R"))
analysis_dir<-sd_analysis_dir("CP8"); rdata_dir<-file.path(analysis_dir,"11_RData"); out_dir<-file.path(analysis_dir,"05_flow_figures_flowjo_style_TALL_LAYOUT_WHITE_BG"); dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
load(file.path(rdata_dir,"SDY2583_CP8_STEP1_inventory_RECONSTRUCTED.RData")); load(file.path(rdata_dir,"SDY2583_CP8_STEP4_composite_scores_RECONSTRUCTED.RData"))
representatives<-rp_select_representative_subjects(scored_data,CP8_MANIFEST$integrated_score)
readr::write_csv(representatives,file.path(out_dir,"SDY2583_CP8_representative_samples_RECONSTRUCTED.csv"))
plots<-list()
for(i in seq_len(nrow(representatives))){
  sid<-representatives$subject_id[i]; group<-as.character(representatives$disease_group[i]); path<-fcs_files[grepl(paste0("^",sid,"_CP8"),basename(fcs_files))][1]
  obj<-rp_read_transform(path); map<-obj$marker_map; ex<-flowCore::exprs(obj$ff)
  ch<-c(CD3=rp_find_channel(map,c("^CD3$")),CD4=rp_find_channel(map,c("^CD4$"),"BB515-A"),CD45RA=rp_find_channel(map,c("CD45RA")),CD62L=rp_find_channel(map,c("CD62L")),IL7RA=rp_find_channel(map,c("IL7RA","CD127")),CD25=rp_find_channel(map,c("^CD25$")),CCR4=rp_find_channel(map,c("CCR4")),CCR6=rp_find_channel(map,c("CCR6")))
  primary<-rp_vec(ex,ch["CD3"])>2 & rp_vec(ex,ch["CD4"])>2
  plots[[length(plots)+1]]<-rp_flow_density_plot(ex,ch["CD3"],ch["CD4"],2,2,paste(group,sid,"CD3/CD4"))
  plots[[length(plots)+1]]<-rp_flow_density_plot(ex[primary,,drop=FALSE],ch["CD45RA"],ch["CD62L"],2,2.2,paste(group,"CD45RA/CD62L"))
  plots[[length(plots)+1]]<-rp_flow_density_plot(ex[primary,,drop=FALSE],ch["IL7RA"],ch["CD25"],1.5,1.5,paste(group,"IL7RA/CD25"))
  plots[[length(plots)+1]]<-rp_flow_density_plot(ex[primary,,drop=FALSE],ch["CCR4"],ch["CCR6"],1.5,1.5,paste(group,"CCR4/CCR6"))
}
rp_save_flow_atlas(plots,file.path(out_dir,"FigureS_CP8_representative_gating_RECONSTRUCTED.pdf"),ncol=2,width=9,height=16)
