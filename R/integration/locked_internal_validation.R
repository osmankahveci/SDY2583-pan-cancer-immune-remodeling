# SDY2583 locked internal held-out validation (v1.1.0)
# Aggregate-only public implementation. Score membership is fixed from the
# project registry; direction/centering/scaling and PCA are learned in TRAINING
# only and frozen before VALIDATION. This is internal, not external, validation.

req <- c("dplyr","tidyr","purrr","readr","stringr","tibble","lmtest","sandwich")
miss <- req[!vapply(req, requireNamespace, logical(1), quietly=TRUE)]
if(length(miss)) stop("Missing packages: ", paste(miss, collapse=", "))
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(readr)
  library(stringr); library(tibble); library(lmtest); library(sandwich)
})

args <- commandArgs(trailingOnly=FALSE)
fa <- grep("^--file=", args, value=TRUE)
sp <- if(length(fa)) sub("^--file=", "", fa[1]) else NA_character_
root0 <- if(!is.na(sp)) normalizePath(file.path(dirname(sp),"../.."),mustWork=FALSE) else getwd()
root <- Sys.getenv("SDY2583_REPO_ROOT", root0)
inroot <- Sys.getenv("SDY2583_LOCKED_INPUT_ROOT", "")
outdir <- Sys.getenv("SDY2583_LOCKED_OUTPUT_DIR", file.path(root,"results","local","locked_internal_validation"))
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

comp <- read_csv(file.path(root,"config","locked_validation","score_component_registry.csv"),show_col_types=FALSE)
reg <- read_csv(file.path(root,"config","locked_validation","score_summary_registry.csv"),show_col_types=FALSE)

resolve <- function(env,names){
  x <- Sys.getenv(env, "")
  if(nzchar(x)){if(!file.exists(x)) stop(env," not found"); return(normalizePath(x))}
  if(!nzchar(inroot) || !dir.exists(inroot)) stop("Set ",env," or SDY2583_LOCKED_INPUT_ROOT")
  h <- list.files(inroot,recursive=TRUE,full.names=TRUE); h <- h[basename(h)%in%names]
  if(length(h)!=1) stop("Could not uniquely resolve ",env,": ",paste(h,collapse=" | "))
  normalizePath(h[1])
}

mfile <- resolve("SDY2583_LOCKED_MATRIX_FILE",c(
  "SDY2583_integrated_clinical_immune_score_matrix_ALL10_with_CP23.csv",
  "SDY2583_integrated_clinical_immune_score_matrix_ALL10_with_CP23(1).csv"))
fn <- list(
 CP7="SDY2583_CP7_analysis_data_with_composite_scores.csv",
 CP8="SDY2583_CP8_analysis_data_with_composite_scores.csv",
 CP10="SDY2583_CP10_composite_scores_data_STEP4.csv",
 CP16="SDY2583_CP16_scores_data_STEP4.csv",
 CP22="SDY2583_CP22_analysis_data_with_composite_scores_STEP4.csv",
 CP23="SDY2583_CP23_score_data_STEP4.csv",
 CP24="SDY2583_CP24_FULL_850_STEP5_score_dataset.csv",
 CP25="SDY2583_CP25_analysis_data_with_composite_scores.csv",
 CP26="SDY2583_CP26_analysis_data_with_composite_scores_STEP4.csv",
 CP28="SDY2583_CP28_analysis_data_with_composite_scores.csv")
paths <- imap_chr(fn,~resolve(paste0("SDY2583_LOCKED_",.y,"_FILE"),.x))

integ <- read_csv(mfile,show_col_types=FALSE)
meta_names <- c("subject_id","cohort","disease_group_model","age_for_clinical_model","sex_for_clinical_model")
if(length(setdiff(meta_names,names(integ)))) stop("Integrated matrix missing required metadata")
meta <- integ %>% select(all_of(meta_names)) %>% distinct(subject_id,.keep_all=TRUE)
if(!all(c("TRAINING","VALIDATION")%in%meta$cohort)) stop("TRAINING/VALIDATION labels missing")

pdata <- imap(paths,function(p,panel){
  x <- read_csv(p,show_col_types=FALSE); if(!"subject_id"%in%names(x)) stop(panel,": subject_id missing")
  if(panel=="CP24"){
    mp <- integ %>% select(canonical_subject_id=subject_id,panel_subject_id=CP24_subject_id_original_panel,
      cohort,disease_group_model,age_for_clinical_model,sex_for_clinical_model) %>%
      distinct(panel_subject_id,.keep_all=TRUE)
    x <- x %>% rename(subject_id_original_panel=subject_id) %>%
      left_join(mp,by=c("subject_id_original_panel"="panel_subject_id")) %>%
      mutate(subject_id=canonical_subject_id) %>% select(-canonical_subject_id)
  } else {
    x <- x %>% select(-any_of(c("cohort","disease_group_model","age_for_clinical_model","sex_for_clinical_model"))) %>%
      left_join(meta,by="subject_id")
  }
  if(any(is.na(x$cohort))) stop(panel,": cohort mapping failed"); x
})

fit <- function(y,m){
  d <- tibble(y=as.numeric(y),disease=as.character(m$disease_group_model),age=as.numeric(m$age_for_clinical_model),
    sex=factor(as.character(m$sex_for_clinical_model))) %>%
    mutate(cancer=ifelse(disease=="Cancer patient",1,0)) %>%
    filter(!is.na(y),!is.na(age),!is.na(sex),disease%in%c("Cancer patient","Healthy control"))
  if(nrow(d)<20 || length(unique(d$cancer))<2) return(NULL)
  f <- lm(y~cancer+age+sex,data=d); z <- coeftest(f,vcov.=vcovHC(f,type="HC3"))
  if(!"cancer"%in%rownames(z)) return(NULL)
  b<-unname(z["cancer",1]); se<-unname(z["cancer",2]); p<-unname(z["cancer",4]); q<-qnorm(.975)
  list(n=nrow(d),beta=b,ci_low=b-q*se,ci_high=b+q*se,p=p)
}

cache <- new.env(parent=emptyenv()); pars <- list(); k <- 0L
minreq <- function(s,n){r<-reg%>%filter(`Composite score`==s)%>%pull(`Coverage rule`); if(length(r)&&str_detect(r[1],"At least half")) ceiling(n/2) else 1L}
build <- function(panel,s){
  key<-paste(panel,s,sep="::"); if(exists(key,cache,inherits=FALSE)) return(get(key,cache))
  x<-pdata[[panel]]; dd<-comp%>%filter(Panel==panel,`Composite score`==s); vv<-vector("list",nrow(dd))
  for(j in seq_len(nrow(dd))){
    c0<-dd$`Component feature or component score`[j]; rd<-as.numeric(dd$Direction[j])
    nested<-any(comp$Panel==panel & comp$`Composite score`==c0)
    raw<-if(nested) build(panel,c0) else if(c0%in%names(x)) as.numeric(x[[c0]]) else stop(panel,"/",s,": missing ",c0)
    tr<-x$cohort=="TRAINING"; mu<-mean(raw[tr],na.rm=TRUE); sd0<-sd(raw[tr],na.rm=TRUE); a<-fit(raw[tr],x[tr,,drop=FALSE])
    tb<-if(is.null(a)) NA_real_ else a$beta; ld<-if(is.na(tb)) ifelse(is.na(rd)||rd>=0,1,-1) else ifelse(tb>=0,1,-1)
    vv[[j]]<-if(is.na(sd0)||sd0==0) rep(NA_real_,length(raw)) else ((raw-mu)/sd0)*ld
    k<<-k+1L; pars[[k]]<<-tibble(panel=panel,score=s,component=c0,registry_direction=rd,
      training_beta_for_orientation=tb,locked_direction=ld,training_mean=mu,training_sd=sd0,
      direction_changed=!is.na(rd)&&rd!=ld)
  }
  M<-do.call(cbind,vv); n0<-rowSums(!is.na(M)); y<-rowMeans(M,na.rm=TRUE); y[n0<minreq(s,ncol(M))|is.nan(y)]<-NA_real_
  assign(key,y,cache); y
}
walk2(reg$Panel,reg$`Composite score`,build)
parameters <- bind_rows(pars)

models <- map2_dfr(reg$Panel,reg$`Composite score`,function(panel,s){
  x<-pdata[[panel]]; y<-get(paste(panel,s,sep="::"),cache)
  map_dfr(c("TRAINING","VALIDATION"),function(cc){ii<-x$cohort==cc;a<-fit(y[ii],x[ii,,drop=FALSE]);
    if(is.null(a)) tibble(panel=panel,score=s,cohort=cc,n=NA,beta=NA,ci_low=NA,ci_high=NA,p=NA)
    else tibble(panel=panel,score=s,cohort=cc,n=a$n,beta=a$beta,ci_low=a$ci_low,ci_high=a$ci_high,p=a$p)})
}) %>% left_join(reg%>%transmute(panel=Panel,score=`Composite score`,principal=`Principal score`),by=c("panel","score"))
models <- models %>% group_by(cohort) %>% mutate(q_bh_66=p.adjust(p,"BH")) %>% ungroup()
prin <- reg%>%filter(`Principal score`=="Yes")%>%transmute(panel=Panel,score=`Composite score`)
p10 <- models%>%semi_join(prin,by=c("panel","score"))%>%group_by(cohort)%>%mutate(q_bh_10=p.adjust(p,"BH"))%>%ungroup()
pw <- p10%>%select(panel,score,cohort,n,beta,ci_low,ci_high,p,q_bh_10)%>%
  pivot_wider(names_from=cohort,values_from=c(n,beta,ci_low,ci_high,p,q_bh_10),names_sep="_")%>%
  mutate(same_direction=sign(beta_TRAINING)==sign(beta_VALIDATION),validation_ci_excludes_zero=ci_low_VALIDATION>0|ci_high_VALIDATION<0)
a66 <- models%>%select(panel,score,principal,cohort,n,beta,ci_low,ci_high,p,q_bh_66)%>%
  pivot_wider(names_from=cohort,values_from=c(n,beta,ci_low,ci_high,p,q_bh_66),names_sep="_")%>%
  mutate(same_direction=sign(beta_TRAINING)==sign(beta_VALIDATION),validation_FDR_sig=q_bh_66_VALIDATION<.05)

S <- meta
for(i in seq_len(nrow(prin))){p<-prin$panel[i];s<-prin$score[i];x<-pdata[[p]];add<-tibble(subject_id=x$subject_id,!!s:=get(paste(p,s,sep="::"),cache))%>%distinct(subject_id,.keep_all=TRUE);S<-left_join(S,add,by="subject_id")}
pc<-prin$score; tr<-S%>%filter(cohort=="TRAINING"); va<-S%>%filter(cohort=="VALIDATION")
cen<-sapply(tr[pc],mean,na.rm=TRUE); scl<-sapply(tr[pc],sd,na.rm=TRUE)
scfun<-function(d){o<-as.data.frame(d[pc]);for(nm in pc)o[[nm]]<-(o[[nm]]-cen[[nm]])/scl[[nm]];o}
Xt<-scfun(tr);Xv<-scfun(va);ct<-complete.cases(Xt);cv<-complete.cases(Xv);pf<-prcomp(Xt[ct,,drop=FALSE],center=FALSE,scale.=FALSE)
pt<-as.numeric(predict(pf,Xt[ct,,drop=FALSE])[,1]);pv<-as.numeric(predict(pf,Xv[cv,,drop=FALSE])[,1]);sg<-ifelse(fit(pt,tr[ct,,drop=FALSE])$beta>=0,1,-1);pt<-pt*sg;pv<-pv*sg
mt<-fit(pt,tr[ct,,drop=FALSE]);mv<-fit(pv,va[cv,,drop=FALSE])
pload<-tibble(score=rownames(pf$rotation),loading_PC1=pf$rotation[,1]*sg,panel=str_extract(rownames(pf$rotation),"^CP[0-9]+"))
psum<-tibble(training_complete_n=sum(ct),validation_complete_n=sum(cv),PC1_training_explained_variance_fraction=(pf$sdev[1]^2)/sum(pf$sdev^2),
 PC1_training_cancer_beta=mt$beta,PC1_training_ci_low=mt$ci_low,PC1_training_ci_high=mt$ci_high,PC1_training_p=mt$p,
 PC1_validation_cancer_beta=mv$beta,PC1_validation_ci_low=mv$ci_low,PC1_validation_ci_high=mv$ci_high,PC1_validation_p=mv$p,
 PC1_orientation="PC1 sign set in TRAINING; VALIDATION projected with frozen TRAINING centering/scaling/loadings")

resid1<-function(s){
 T<-tr%>%transmute(y=.data[[s]],age=as.numeric(age_for_clinical_model),sex=factor(as.character(sex_for_clinical_model)),cancer=ifelse(disease_group_model=="Cancer patient",1,0))
 V<-va%>%transmute(y=.data[[s]],age=as.numeric(age_for_clinical_model),sex=factor(as.character(sex_for_clinical_model),levels=levels(T$sex)),cancer=ifelse(disease_group_model=="Cancer patient",1,0))
 f<-lm(y~age+sex+cancer,data=T[complete.cases(T),,drop=FALSE]);list(train=T$y-predict(f,newdata=T),validation=V$y-predict(f,newdata=V))}
RL<-setNames(lapply(pc,resid1),pc);RT<-as.data.frame(lapply(RL,`[[`,"train"),check.names=FALSE);RV<-as.data.frame(lapply(RL,`[[`,"validation"),check.names=FALSE)
cr<-map_dfr(combn(pc,2,simplify=FALSE),function(z){a<-complete.cases(RT[,z]);b<-complete.cases(RV[,z]);t1<-cor.test(RT[[z[1]]][a],RT[[z[2]]][a],method="spearman",exact=FALSE);t2<-cor.test(RV[[z[1]]][b],RV[[z[2]]][b],method="spearman",exact=FALSE);tibble(score1=z[1],score2=z[2],panel1=str_extract(z[1],"^CP[0-9]+"),panel2=str_extract(z[2],"^CP[0-9]+"),n_train=sum(a),rho_train=unname(t1$estimate),p_train=t1$p.value,n_val=sum(b),rho_val=unname(t2$estimate),p_val=t2$p.value)})%>%
 mutate(q_train=p.adjust(p_train,"BH"),q_val=p.adjust(p_val,"BH"),same_direction=sign(rho_train)==sign(rho_val))
cp<-cor.test(cr$rho_train,cr$rho_val);cs<-cor.test(cr$rho_train,cr$rho_val,method="spearman",exact=FALSE)
crsum<-tibble(n_pairs=nrow(cr),same_direction_pairs=sum(cr$same_direction),same_direction_fraction=mean(cr$same_direction),validation_pairs_FDR_lt_0_05=sum(cr$q_val<.05),training_validation_rho_concordance_Pearson_r=unname(cp$estimate),training_validation_rho_concordance_Pearson_p=cp$p.value,training_validation_rho_concordance_Spearman_rho=unname(cs$estimate),training_validation_rho_concordance_Spearman_p=cs$p.value)
orient<-parameters%>%group_by(panel)%>%summarise(n_component_assignments=n(),n_direction_changed=sum(direction_changed),.groups="drop")
over<-tibble(metric=c("Original training subjects","Original validation subjects","Principal scores same direction in validation","Principal scores validation BH-FDR < 0.05","All 66 scores same direction in validation","All 66 scores validation BH-FDR < 0.05","PC1 training explained variance fraction","Cross-panel correlations same direction","Cross-panel validation correlations BH-FDR < 0.05"),value=c(sum(meta$cohort=="TRAINING"),sum(meta$cohort=="VALIDATION"),sum(pw$same_direction),sum(pw$q_bh_10_VALIDATION<.05),sum(a66$same_direction),sum(a66$validation_FDR_sig),psum$PC1_training_explained_variance_fraction,crsum$same_direction_pairs,crsum$validation_pairs_FDR_lt_0_05))

write_csv(pw,file.path(outdir,"locked_validation_principal_scores.csv"));write_csv(a66,file.path(outdir,"locked_validation_all66_scores.csv"));write_csv(parameters,file.path(outdir,"locked_component_parameters.csv"));write_csv(pload,file.path(outdir,"locked_validation_pc1_loadings.csv"));write_csv(psum,file.path(outdir,"locked_validation_pca_summary.csv"));write_csv(cr,file.path(outdir,"locked_validation_cross_panel_correlations.csv"));write_csv(crsum,file.path(outdir,"locked_validation_cross_panel_summary.csv"));write_csv(orient,file.path(outdir,"locked_validation_orientation_audit.csv"));write_csv(over,file.path(outdir,"locked_validation_overall_summary.csv"))
cat("Locked internal validation complete. Aggregate outputs: ",outdir,"\n",sep="")
