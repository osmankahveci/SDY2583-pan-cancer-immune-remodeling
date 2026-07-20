# Shared figure helpers for reconstructed SDY2583 panels.
source(file.path(Sys.getenv("SDY2583_REPO_ROOT", unset="."),"R","shared","reconstructed_panel_framework.R"))

rp_quantitative_figure <- function(data, features, statistics, output_file, labels=NULL, ncol=2) {
  rp_install_and_load(c("dplyr","tidyr","ggplot2","readr","stringr","patchwork"))
  labels <- labels %||% stats::setNames(features,features)
  long <- data |> dplyr::select(subject_id,disease_group,dplyr::all_of(features)) |>
    tidyr::pivot_longer(dplyr::all_of(features),names_to="feature",values_to="value") |>
    dplyr::mutate(feature=factor(feature,levels=features,labels=unname(labels[features])))
  ann <- statistics |> dplyr::filter(feature %in% features) |>
    dplyr::mutate(label=dplyr::case_when(fdr_all<.001~"FDR < 0.001",fdr_all<.01~paste0("FDR = ",formatC(fdr_all,format="g",digits=2)),fdr_all<.05~paste0("FDR = ",formatC(fdr_all,format="g",digits=2)),TRUE~"ns")) |>
    dplyr::select(feature,label)
  ymax <- long |> dplyr::group_by(feature) |> dplyr::summarise(y=max(value,na.rm=TRUE),.groups="drop")
  ann <- ann |> dplyr::left_join(ymax,by="feature")
  p <- ggplot2::ggplot(long,ggplot2::aes(disease_group,value))+
    ggplot2::geom_boxplot(outlier.shape=NA,width=.6)+ggplot2::geom_jitter(width=.14,alpha=.22,size=.55)+
    ggplot2::facet_wrap(~feature,scales="free_y",ncol=ncol)+ggplot2::theme_bw(base_size=10)+
    ggplot2::theme(axis.title.x=ggplot2::element_blank(),axis.text.x=ggplot2::element_text(angle=20,hjust=1),strip.text=ggplot2::element_text(face="bold"))+
    ggplot2::labs(y="Feature value")
  if(nrow(ann)>0) p<-p+ggplot2::geom_text(data=ann,ggplot2::aes(x=1.5,y=y,label=label),inherit.aes=FALSE,vjust=-.3,size=3)
  ggplot2::ggsave(output_file,p,width=8.3,height=max(5,3.1*ceiling(length(features)/ncol)),units="in")
  invisible(p)
}

rp_select_representative_subjects <- function(data, score) {
  d <- data |> dplyr::filter(!is.na(.data[[score]]),!is.na(disease_group)) |>
    dplyr::group_by(disease_group) |> dplyr::mutate(group_median=stats::median(.data[[score]],na.rm=TRUE),distance=abs(.data[[score]]-group_median)) |>
    dplyr::slice_min(distance,n=1,with_ties=FALSE) |> dplyr::ungroup()
  d |> dplyr::select(subject_id,result_file_name,disease_group,dplyr::all_of(score),group_median,distance)
}

rp_flow_density_plot <- function(exprs, x, y, x_threshold=NULL, y_threshold=NULL, title=NULL, sample_n=60000) {
  rp_install_and_load(c("ggplot2","tibble","dplyr"))
  n <- nrow(exprs); idx <- if(n>sample_n) sample.int(n,sample_n) else seq_len(n)
  dat <- tibble::tibble(x=as.numeric(exprs[idx,x]),y=as.numeric(exprs[idx,y]))
  p <- ggplot2::ggplot(dat,ggplot2::aes(x,y))+ggplot2::stat_bin2d(bins=90)+ggplot2::stat_density_2d(bins=8,linewidth=.25)+
    ggplot2::theme_bw(base_size=10)+ggplot2::labs(x=x,y=y,title=title)
  if(!is.null(x_threshold)) p<-p+ggplot2::geom_vline(xintercept=x_threshold,linetype="dashed")
  if(!is.null(y_threshold)) p<-p+ggplot2::geom_hline(yintercept=y_threshold,linetype="dashed")
  p
}

rp_save_flow_atlas <- function(plots, output_file, ncol=2, width=9, height=NULL) {
  rp_install_and_load(c("patchwork","ggplot2"))
  if(is.null(height)) height<-4*ceiling(length(plots)/ncol)
  atlas<-patchwork::wrap_plots(plots,ncol=ncol)
  ggplot2::ggsave(output_file,atlas,width=width,height=height,units="in")
  invisible(atlas)
}
