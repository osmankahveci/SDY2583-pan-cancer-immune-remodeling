# RECONSTRUCTED representative supplementary flow-atlas workflow for SDY2583.
# Raw FCS files and participant identifiers remain local. Run from repo root.

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
if (!requireNamespace("flowCore", quietly = TRUE))
  BiocManager::install("flowCore", ask = FALSE, update = FALSE)
for (p in c("ggplot2", "patchwork", "ggrastr", "MASS"))
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = "https://cloud.r-project.org")
suppressPackageStartupMessages({library(flowCore); library(ggplot2); library(patchwork)})
set.seed(2583)

Q <- function(parent,x,y,xt,yt,title,labels,xlab=x,ylab=y,parent_text=parent)
  list(type="q",parent=parent,x=x,y=y,xt=xt,yt=yt,title=title,labels=labels,
       xlab=xlab,ylab=ylab,parent_text=parent_text)
H <- function(parent,x,t,title,neg,pos,xlab=x,parent_text=parent)
  list(type="h",parent=parent,x=x,t=t,title=title,neg=neg,pos=pos,
       xlab=xlab,parent_text=parent_text)
C <- function(...) c(...)

cfg <- list(
CP7=list(S="S1",env="SDY2583_CP7_REPRESENTATIVE_FCS",comp=TRUE,ncol=2,w=9.5,h=11.5,
 ch=C(TIM3="BV421-A",DUMP="BV510-A",CD3="BV605-A",CD62L="BV650-A",CD27="BV711-A",LAG3="BV786-A",ICOS="BB515-A",CD8="PerCP-Cy5-5-A",TIGIT="PE-A",PD1="PE-CF594-A",CD45RA="PE-Cy5-A",CD39="PE-Cy7-A"),
 t=C(CD3=2,CD8=2.2,CD45RA=2,CD62L=2.2,CD27=2,PD1=1.5,TIM3=1.5,LAG3=1.5,TIGIT=1.5,ICOS=1.5,CD39=1.5),
 min=C(cd3cd8=300,temra=50),title="Representative CP7 CD8 differentiation and checkpoint atlas",
 cap="Percentages use all events in the displayed parent. The combined viability/CD4/CD13/CD19/TCRγδ channel is not an independent CD4 marker. Gates document the automated analytical architecture, not exact manual FlowJo polygons.",
 plots=list(
  Q("all","CD3","CD8",2,2.2,"CD3 × CD8",c("CD3−CD8+","CD3+CD8+","CD3−CD8−","CD3+CD8−"),parent_text="Total events"),
  Q("cd3cd8","CD45RA","CD62L",2,2.2,"CD45RA × CD62L",c("CD45RA−CD62L+","CD45RA+CD62L+","CD45RA−CD62L−","CD45RA+CD62L−"),parent_text="Within CD3+CD8+ events"),
  Q("cd3cd8","PD1","CD39",1.5,1.5,"PD-1 × CD39",c("PD-1−CD39+","PD-1+CD39+","PD-1−CD39−","PD-1+CD39−"),xlab="PD-1",parent_text="Within CD3+CD8+ events"),
  Q("cd3cd8","TIM3","LAG3",1.5,1.5,"TIM-3 × LAG-3",c("TIM-3−LAG-3+","TIM-3+LAG-3+","TIM-3−LAG-3−","TIM-3+LAG-3−"),xlab="TIM-3",ylab="LAG-3",parent_text="Within CD3+CD8+ events"),
  Q("temra","TIGIT","CD39",1.5,1.5,"TIGIT × CD39",c("TIGIT−CD39+","TIGIT+CD39+","TIGIT−CD39−","TIGIT+CD39−"),parent_text="Within TEMRA-like events"),
  Q("temra","ICOS","CD39",1.5,1.5,"ICOS × CD39",c("ICOS−CD39+","ICOS+CD39+","ICOS−CD39−","ICOS+CD39−"),parent_text="Within TEMRA-like events"))),
CP8=list(S="S2",env="SDY2583_CP8_REPRESENTATIVE_FCS",comp=TRUE,ncol=2,w=11.69,h=8.27,
 ch=C(CD3="BV605-A",CD4="BB515-A",CD45RA="BV786-A",CD62L="BV650-A",CD27="BV711-A",CCR6="PerCP-Cy5-5-A",IL7RA="PE-A",CCR4="PE-CF594-A",CD25="PE-Cy5-A",CXCR5="PE-Cy7-A"),
 t=C(CD3=2,CD4=2,CD45RA=2,CD62L=2.2,CD27=2,CCR6=1.5,IL7RA=1.5,CCR4=1.5,CD25=1.5,CXCR5=1.5),min=C(cd3cd4=300),
 title="Representative CP8 CD4 helper/regulatory-like atlas",cap="Percentages use all events in the displayed parent. Marker-defined states are phenotype-like analytical compartments and do not establish functional helper or regulatory identity.",
 plots=list(
  Q("all","CD3","CD4",2,2,"CD3 × CD4",c("CD3−CD4+","CD3+CD4+","CD3−CD4−","CD3+CD4−"),parent_text="Total events"),
  Q("cd3cd4","CD25","IL7RA",1.5,1.5,"CD25 × IL7RA",c("CD25−IL7RA+","CD25+IL7RA+","CD25−IL7RA-low","CD25+IL7RA-low"),parent_text="Within CD3+CD4+ events"),
  Q("cd3cd4","CCR4","CCR6",1.5,1.5,"CCR4 × CCR6",c("CCR4−CCR6+","CCR4+CCR6+","CCR4−CCR6−","CCR4+CCR6−"),parent_text="Within CD3+CD4+ events"),
  Q("cd3cd4","CCR6","CXCR5",1.5,1.5,"CCR6 × CXCR5",c("CCR6−CXCR5+","CCR6+CXCR5+","CCR6−CXCR5−","CCR6+CXCR5−"),parent_text="Within CD3+CD4+ events"))),
CP10=list(S="S3",env="SDY2583_CP10_REPRESENTATIVE_FCS",comp=TRUE,ncol=2,w=11.69,h=8.27,
 ch=C(Viability="BV510-A",CD45="BV786-A",CD14="BB515-A",HLA_DR="BV650-A",CD11c="PE-Cy7-A",CD13="PerCP-Cy5-5-A",CD66b="PE-A"),
 t=C(Viability=1.5,CD45=2,CD14=1.5,HLA_DR=1.5,CD11c=1.5,CD13=1.5,CD66b=1.5),min=C(primary=300),
 title="Representative CP10 myeloid/granulocytic atlas",cap="Percentages use all events in the displayed parent. Marker-defined compartments are descriptive phenotypes rather than definitive lineage or functional assignments.",
 plots=list(
  Q("all","Viability","CD45",1.5,2,"Viability × CD45",c("Viability-low CD45+","Viability-high CD45+","Viability-low CD45−","Viability-high CD45−"),xlab="Viability / dump",parent_text="Total events"),
  Q("primary","CD13","CD66b",1.5,1.5,"CD13 × CD66b",c("CD13−CD66b+","CD13+CD66b+","CD13−CD66b−","CD13+CD66b−"),parent_text="Within viable CD45+ events"),
  Q("primary","CD14","HLA_DR",1.5,1.5,"CD14 × HLA-DR",c("CD14−HLA-DR+","CD14+HLA-DR+","CD14−HLA-DR−","CD14+HLA-DR−"),ylab="HLA-DR",parent_text="Within viable CD45+ events"),
  Q("primary","CD11c","HLA_DR",1.5,1.5,"CD11c × HLA-DR",c("CD11c−HLA-DR+","CD11c+HLA-DR+","CD11c−HLA-DR−","CD11c+HLA-DR−"),ylab="HLA-DR",parent_text="Within viable CD45+ events"))),
CP16=list(S="S4",env="SDY2583_CP16_REPRESENTATIVE_FCS",comp=TRUE,ncol=2,w=9.5,h=11.5,
 ch=C(DUMP="BV510-A",CD45="BB515-A",HLA_DR="BV650-A",CD11c="PE-Cy7-A",CD14="PE-CF594-A",CD1c="BV711-A",FceRI="PE-A",CD123="PerCP-Cy5-5-A",CD141="BV786-A",CLEC9A="BV421-A"),
 t=C(DUMP=1.5,CD45=2,HLA_DR=1.5,CD11c=1.5,CD14=1.5,CD1c=1.5,FceRI=1.5,CD123=1.5,CD141=1.5,CLEC9A=1.5),min=C(primary=50,apc=50),
 title="Representative CP16 APC/DC-like atlas",cap="Percentages use all events in the displayed parent. CP16 populations are APC/DC-like or monocyte-like phenotype compartments; low event counts and limited lineage markers preclude definitive dendritic-cell identity. CD141/CLEC9A is supportive and threshold-sensitive.",
 plots=list(
  Q("all","DUMP","CD45",1.5,2,"Dump-low × CD45",c("Dump-low CD45+","Dump-high CD45+","Dump-low CD45−","Dump-high CD45−"),xlab="Viability / dump",parent_text="Total events"),
  Q("primary","CD11c","HLA_DR",1.5,1.5,"CD11c × HLA-DR",c("CD11c−HLA-DR+","CD11c+HLA-DR+","CD11c−HLA-DR−","CD11c+HLA-DR−"),ylab="HLA-DR",parent_text="Within CD45+ dump-low events"),
  Q("primary","CD14","HLA_DR",1.5,1.5,"CD14 × HLA-DR",c("CD14−HLA-DR+","CD14+HLA-DR+","CD14−HLA-DR−","CD14+HLA-DR−"),ylab="HLA-DR",parent_text="Within CD45+ dump-low events"),
  Q("apc","CD1c","FceRI",1.5,1.5,"CD1c × FcεRI",c("CD1c−FcεRI+","CD1c+FcεRI+","CD1c−FcεRI−","CD1c+FcεRI−"),ylab="FcεRI",parent_text="Within HLA-DR+ APC-like events"),
  Q("apc","CD11c","CD123",1.5,1.5,"CD11c × CD123",c("CD11c−CD123+","CD11c+CD123+","CD11c−CD123−","CD11c+CD123−"),parent_text="Within HLA-DR+ APC-like events"),
  Q("apc","CD141","CLEC9A",1.5,1.5,"CD141 × CLEC9A",c("CD141−CLEC9A+","CD141+CLEC9A+","CD141−CLEC9A−","CD141+CLEC9A−"),parent_text="Within HLA-DR+ APC-like events"))),
CP22=list(S="S5",env="SDY2583_CP22_REPRESENTATIVE_FCS",comp=FALSE,ncol=2,w=9.5,h=11.5,
 ch=C(CD138="BV421-A",DUMP="BV510-A",IgG="BV605-A",CD39="BV650-A",CD24="BV711-A",CD10="BV786-A",CD19="BB515-A",IgD="PerCP-Cy5-5-A",IgA="PE-A",IgM="PE-CF594-A",CD38="PE-Cy5-A",CD27="PE-Cy7-A"),
 t=C(DUMP_LOW=1.5,CD19=2,IgD=1.5,IgA=1.5,IgG=1.5,CD27=1.5,CD38=1.5,CD38_HIGH=2,CD10=1.5,CD24=1.5,CD39=1.5),min=C(bcells=300),
 title="Representative CP22 B-cell and humoral architecture atlas",cap="Percentages use all events in the displayed parent. IgA/IgG panels illustrate B-cell isotype architecture, not total circulating IgA or IgG abundance. Plasmablast-, transitional-, and regulatory-like labels are phenotypic and not functional validation.",
 plots=list(
  Q("all","DUMP","CD19",1.5,2,"Dump-low × CD19",c("Dump-low CD19+","Dump-high CD19+","Dump-low CD19−","Dump-high CD19−"),xlab="Viability / dump",parent_text="Total events"),
  Q("bcells","IgD","CD27",1.5,1.5,"IgD × CD27",c("IgD−CD27+","IgD+CD27+","IgD−CD27−","IgD+CD27−"),parent_text="Within dump-low CD19+ B cells"),
  Q("bcells","IgA","IgG",1.5,1.5,"IgA × IgG",c("IgA−IgG+","IgA+IgG+","IgA−IgG−","IgA+IgG−"),parent_text="Within dump-low CD19+ B cells"),
  Q("bcells","CD27","CD38",1.5,2,"CD27 × CD38",c("CD27−CD38-high","CD27+CD38-high","CD27−CD38-low","CD27+CD38-low"),parent_text="Within dump-low CD19+ B cells"),
  Q("bcells","CD10","CD24",1.5,1.5,"CD10 × CD24",c("CD10−CD24+","CD10+CD24+","CD10−CD24−","CD10+CD24−"),parent_text="Within dump-low CD19+ B cells"),
  Q("bcells","CD24","CD39",1.5,1.5,"CD24 × CD39",c("CD24−CD39+","CD24+CD39+","CD24−CD39−","CD24+CD39−"),parent_text="Within dump-low CD19+ B cells"))),
CP23=list(S="S6",env="SDY2583_CP23_REPRESENTATIVE_FCS",comp=TRUE,ncol=2,w=9.5,h=11.5,
 ch=C(CD16="BV421-A",DUMP="BV510-A",FceRI="BV605-A",HLA_DR="BV650-A",CD33="BV711-A",CD45="BV786-A",CD14="BB515-A",CD9="PerCP-Cy5-5-A",CD84="PE-A",CD15="PE-CF594-A",CD206="PE-Cy5-A",CD169="PE-Cy7-A"),
 t=C(DUMP=1.5,CD45=2,HLA_DR=1.5,CD14=1.5,CD16=1.5,CD33=1.5,FceRI=1.5,CD9=1.5,CD84=1.5),min=C(primary=50,cd14=50),
 title="Representative CP23 monocyte/myeloid phenotype atlas",cap="Percentages use all events in the displayed parent. Populations are threshold-defined phenotype-like compartments and are not definitive macrophage, MDSC, basophil, or M1/M2 assignments.",
 plots=list(
  Q("all","DUMP","CD45",1.5,2,"Dump-low × CD45",c("Dump-low CD45+","Dump-high CD45+","Dump-low CD45−","Dump-high CD45−"),xlab="Viability / dump",parent_text="Total events"),
  Q("primary","CD14","HLA_DR",1.5,1.5,"CD14 × HLA-DR",c("CD14−HLA-DR+","CD14+HLA-DR+","CD14−HLA-DR−","CD14+HLA-DR−"),ylab="HLA-DR",parent_text="Within CD45+ dump-low events"),
  Q("primary","CD14","CD16",1.5,1.5,"CD14 × CD16",c("CD14−CD16+","CD14+CD16+","CD14−CD16−","CD14+CD16−"),parent_text="Within CD45+ dump-low events"),
  Q("primary","CD33","HLA_DR",1.5,1.5,"CD33 × HLA-DR",c("CD33−HLA-DR+","CD33+HLA-DR+","CD33−HLA-DR−","CD33+HLA-DR−"),ylab="HLA-DR",parent_text="Within CD45+ dump-low events"),
  Q("cd14","CD9","CD84",1.5,1.5,"CD9 × CD84",c("CD9−CD84+","CD9+CD84+","CD9−CD84−","CD9+CD84−"),parent_text="Within CD14+ monocyte-like events"),
  Q("primary","FceRI","HLA_DR",1.5,1.5,"FcεRI × HLA-DR",c("FcεRI−HLA-DR+","FcεRI+HLA-DR+","FcεRI−HLA-DR−","FcεRI+HLA-DR−"),xlab="FcεRI",ylab="HLA-DR",parent_text="Within CD45+ dump-low events"))),
CP24=list(S="S7",env="SDY2583_CP24_REPRESENTATIVE_FCS",comp=TRUE,ncol=4,w=16.54,h=8.6,
 ch=C(CXCR3="BV421-A",DUMP="BV510-A",CD3="BV605-A",CD62L="BV650-A",CD95="BV711-A",CD57="BV786-A",CD27="BB515-A",CD8="PerCP-Cy5-5-A",CX3CR1="PE-A",PD1="PE-CF594-A",CD45RA="PE-Cy5-A",CXCR5="PE-Cy7-A"),
 t=C(CD3=2,CD8=2.2,CD45RA=2,CD62L=2.2,CD57=1.5,CX3CR1=1.8,CD95=2,CD27=2,PD1=1.5),min=C(cd3cd8=300),
 title="Representative CP24 CD8 differentiation and effector atlas",cap="Percentages use all events in the displayed parent. The combined viability/CD4/CD13/CD19/TCRγδ channel is not an independent CD4 marker. Quantitative inference derives from the full cohort, not this representative sample.",
 plots=list(
  Q("all","CD3","CD8",2,2.2,"CD3 × CD8",c("CD3−CD8+","CD3+CD8+","CD3−CD8−","CD3+CD8−"),parent_text="Total events"),
  Q("cd3cd8","CD45RA","CD62L",2,2.2,"CD45RA × CD62L",c("TCM-like","Naive-like","TEM-like","TEMRA-like"),parent_text="Within CD3+CD8+ events"),
  Q("cd3cd8","CD62L","CD27",2.2,2,"CD62L × CD27",c("CD62L−CD27+","CD62L+CD27+","CD62L−CD27−","CD62L+CD27−"),parent_text="Within CD3+CD8+ events"),
  Q("cd3cd8","CD57","CX3CR1",1.5,1.8,"CD57 × CX3CR1",c("CD57−CX3CR1+","CD57+CX3CR1+","CD57−CX3CR1−","CD57+CX3CR1−"),parent_text="Within CD3+CD8+ events"),
  Q("cd3cd8","CD57","CD95",1.5,2,"CD57 × CD95",c("CD57−CD95+","CD57+CD95+","CD57−CD95−","CD57+CD95−"),parent_text="Within CD3+CD8+ events"),
  Q("cd3cd8","CX3CR1","CD95",1.8,2,"CX3CR1 × CD95",c("CX3CR1−CD95+","CX3CR1+CD95+","CX3CR1−CD95−","CX3CR1+CD95−"),parent_text="Within CD3+CD8+ events"),
  Q("cd3cd8","PD1","CD27",1.5,2,"PD-1 × CD27",c("PD-1−CD27+","PD-1+CD27+","PD-1−CD27−","PD-1+CD27−"),xlab="PD-1",parent_text="Within CD3+CD8+ events"))),
CP25=list(S="S8",env="SDY2583_CP25_REPRESENTATIVE_FCS",comp=TRUE,ncol=2,w=9.5,h=11.5,
 ch=C(CD25="BV421-A",DUMP="BV510-A",CD3="BV605-A",ICOS="BV650-A",CD27="BV711-A",LAG3="BV786-A",CTLA4="BB515-A",CD4="PerCP-Cy5-5-A",IL7RA="PE-A",PD1="PE-CF594-A",CD45RA="PE-Cy5-A",CD39="PE-Cy7-A"),
 t=C(CD3=2,CD4=2,DUMP=1.5,CD25=1.5,IL7RA=1.5,CTLA4=1.5,ICOS=1.5,CD39=1.5,PD1=1.5,LAG3=1.5),min=C(cd3cd4=300,regulatory=50),
 title="Representative CP25 CD4 regulatory-checkpoint atlas",cap="Percentages use all events in the displayed parent. FOXP3 is absent; CD25+IL7RA-low cells are regulatory-like or Treg-enriched rather than definitive regulatory T cells.",
 plots=list(
  Q("dump_low","CD3","CD4",2,2,"CD3 × CD4",c("CD3−CD4+","CD3+CD4+","CD3−CD4−","CD3+CD4−"),parent_text="Within dump-low events"),
  Q("cd3cd4","CD25","IL7RA",1.5,1.5,"CD25 × IL7RA",c("CD25−IL7RA+","CD25+IL7RA+","CD25−IL7RA-low","CD25+IL7RA-low"),parent_text="Within CD3+CD4+ dump-low events"),
  Q("regulatory","PD1","LAG3",1.5,1.5,"PD-1 × LAG-3",c("PD-1−LAG-3+","PD-1+LAG-3+","PD-1−LAG-3−","PD-1+LAG-3−"),xlab="PD-1",ylab="LAG-3",parent_text="Within CD25+IL7RA-low regulatory-like CD4 events"),
  Q("regulatory","ICOS","CD39",1.5,1.5,"ICOS × CD39",c("ICOS−CD39+","ICOS+CD39+","ICOS−CD39−","ICOS+CD39−"),parent_text="Within CD25+IL7RA-low regulatory-like CD4 events"),
  H("regulatory","CTLA4",1.5,"CTLA-4 expression","CTLA-4−","CTLA-4+",xlab="CTLA-4",parent_text="Within CD25+IL7RA-low regulatory-like CD4 events"))),
CP26=list(S="S9",env="SDY2583_CP26_REPRESENTATIVE_FCS",comp=TRUE,ncol=2,w=9.5,h=11.5,
 ch=C(CD16="BV421-A",DUMP="BV510-A",NKG2A="BV605-A",CD158="BV650-A",NKG2C="BV711-A",CD57="BV786-A",CD45="BB515-A",NKp44="PerCP-Cy5-5-A",CD161="PE-A",CD56="PE-CF594-A",CD107a="PE-Cy5-A",NKG2D="PE-Cy7-A"),
 t=C(DUMP=1.5,CD45=2,CD56=1.5,CD16=1.5,NKG2C=1.5,NKG2D=1.5,CD57=1.5,CD161=1.5),min=C(primary=50,nk_like=300),
 title="Representative CP26 NK-like receptor and maturation atlas",cap="Percentages use all events in the displayed parent. NKG2D and CD161 are principal receptor axes; NKG2C and CD57 are supportive maturation/adaptive-like readouts and require cautious interpretation.",
 plots=list(
  Q("all","DUMP","CD45",1.5,2,"Dump-low × CD45",c("Dump-low CD45+","Dump-high CD45+","Dump-low CD45−","Dump-high CD45−"),xlab="Viability / dump",parent_text="Total events"),
  Q("primary","CD16","CD56",1.5,1.5,"CD16 × CD56",c("CD16−CD56+","CD16+CD56+","CD16−CD56−","CD16+CD56−"),parent_text="Within CD45+ dump-low events"),
  H("nk_like","NKG2D",1.5,"NKG2D expression","NKG2D−","NKG2D+",parent_text="Within NK-like total events"),
  H("nk_like","CD161",1.5,"CD161 expression","CD161−","CD161+",parent_text="Within NK-like total events"),
  H("nk_like","NKG2C",1.5,"NKG2C expression","NKG2C−","NKG2C+",parent_text="Within NK-like total events"),
  H("nk_like","CD57",1.5,"CD57 expression","CD57−","CD57+",parent_text="Within NK-like total events"))),
CP28=list(S="S10",env="SDY2583_CP28_REPRESENTATIVE_FCS",comp=TRUE,ncol=3,w=13.2,h=5.5,
 ch=C(DUMP="BV510-A",CD3="BV605-A",VA24_JA11_TCR="BV650-A",VDELTA2="BV711-A",CD57="BV786-A",CD27="BB515-A",CD8="PerCP-Cy5-5-A",CD161="PE-A",CD56="PE-CF594-A",CD45RA="PE-Cy5-A",VALPHA7="PE-Cy7-A"),
 t=C(DUMP=1.5,CD3=2,CD8=2,CD56=1.5,VA24_JA11_TCR=1.5),min=C(dump_low=300,cd3=300),
 title="Representative CP28 T/NK-interface atlas",cap="Percentages use all events in the displayed parent. CD3−CD56+ events are NK-like because CP28 lacks CD16. Vα24-Jα11/CD8 features are innate-like T-cell readouts and do not establish definitive iNKT identity.",
 plots=list(
  Q("dump_low","CD3","CD56",2,1.5,"CD3 × CD56",c("CD3−CD56+ NK-like","CD3+CD56+ NK-like T","CD3−CD56−","CD3+CD56−"),parent_text="Within dump-low events"),
  Q("dump_low","CD3","CD8",2,2,"CD3 × CD8",c("CD3−CD8+","CD3+CD8+","CD3−CD8−","CD3+CD8−"),parent_text="Within dump-low events"),
  Q("cd3","VA24_JA11_TCR","CD8",1.5,2,"Vα24-Jα11 TCR × CD8",c("Vα24-Jα11−CD8+","Vα24-Jα11+CD8+ innate-like","Vα24-Jα11−CD8−","Vα24-Jα11+CD8−"),xlab="Vα24-Jα11 TCR",parent_text="Within CD3+ dump-low events")))
)

panel <- toupper(trimws(Sys.getenv("SDY2583_FLOW_PANEL", "")))
if (!panel %in% names(cfg)) stop("Set SDY2583_FLOW_PANEL to: ", paste(names(cfg), collapse=", "))
z <- cfg[[panel]]
fcs <- path.expand(Sys.getenv(z$env, ""))
if (!nzchar(fcs) || !file.exists(fcs)) stop("Set ", z$env, " to an existing local FCS file.")
out <- path.expand(Sys.getenv("SDY2583_FLOW_OUTPUT_DIR", file.path(getwd(),"outputs","supplementary_flow")))
dir.create(out, recursive=TRUE, showWarnings=FALSE)

parse_spill <- function(x){
 if(length(x)!=1L || !is.character(x)) return(NULL); a<-strsplit(x,",",fixed=TRUE)[[1]]; n<-suppressWarnings(as.integer(a[1]));
 if(is.na(n)||length(a)<1+n+n*n) return(NULL); v<-suppressWarnings(as.numeric(a[seq.int(2+n,1+n+n*n)])); if(anyNA(v)) return(NULL);
 matrix(v,n,n,byrow=TRUE,dimnames=list(a[2:(1+n)],a[2:(1+n)]))}
get_spill <- function(ff){
 s<-tryCatch(flowCore::spillover(ff),error=function(e)NULL); if(is.matrix(s)||inherits(s,"compensation")) return(s)
 if(is.list(s)){k<-vapply(s,function(x)is.matrix(x)||inherits(x,"compensation"),logical(1));if(any(k))return(s[[which(k)[1]]])}
 kw<-flowCore::keyword(ff);for(n in c("$SPILLOVER","SPILLOVER","$SPILL","SPILL")){o<-kw[[n]];if(is.null(o))next;if(is.matrix(o)||inherits(o,"compensation"))return(o);if(is.list(o)&&length(o)==1)o<-o[[1]];q<-parse_spill(o);if(!is.null(q))return(q)};NULL}

ff<-flowCore::read.FCS(fcs,transformation=FALSE,truncate_max_range=FALSE)
pd<-Biobase::pData(flowCore::parameters(ff)); desc<-if("desc"%in%names(pd))as.character(pd$desc) else rep(NA_character_,ncol(ff)); desc[is.na(desc)|desc==""]<-colnames(ff)[is.na(desc)|desc==""]
audit<-data.frame(channel=colnames(ff),marker_description=desc,stringsAsFactors=FALSE)
norm<-function(x)toupper(gsub("[^A-Z0-9]","",x))
alias<-list(DUMP=c("VIABILITY","DUMP"),HLA_DR=c("HLADR"),FceRI=c("FCERI","FCER1"),PD1=c("PD1"),TIM3=c("TIM3"),LAG3=c("LAG3"),CTLA4=c("CTLA4"),VA24_JA11_TCR=c("VA24JA11TCR","VALPHA24JALPHA11TCR"),VDELTA2=c("VDELTA2","VD2"),VALPHA7=c("VALPHA7","VA7"))
resolve<-function(m){
 keys<-unique(c(norm(m),alias[[m]])); hit<-which(vapply(norm(desc),function(d)any(vapply(keys,function(k)grepl(k,d,fixed=TRUE),logical(1))),logical(1)))
 if(length(hit)==1)return(colnames(ff)[hit]); if(length(hit)>1 && z$ch[[m]]%in%colnames(ff)[hit])return(z$ch[[m]]); if(z$ch[[m]]%in%colnames(ff))return(z$ch[[m]]); stop("Missing/ambiguous marker: ",m)}
map<-vapply(names(z$ch),resolve,character(1));audit$resolved_marker<-NA_character_;for(m in names(map))audit$resolved_marker[audit$channel==map[[m]]]<-m
write.csv(audit,file.path(out,paste0(panel,"_channel_audit.csv")),row.names=FALSE)
if(z$comp){sp<-get_spill(ff);if(is.null(sp))stop("No usable spillover matrix for ",panel);ff<-flowCore::compensate(ff,sp)}
tr<-flowCore::logicleTransform("fixed_logicle",w=.5,t=262144,m=4.5,a=0)
ff<-flowCore::transform(ff,flowCore::transformList(unique(map),tr)); e<-flowCore::exprs(ff);dat<-as.data.frame(lapply(map,function(ch)as.numeric(e[,ch])),check.names=FALSE);names(dat)<-names(map)

parents<-switch(panel,
 CP7={cd<-dat[dat$CD3>2&dat$CD8>2.2,,drop=FALSE];list(all=dat,cd3cd8=cd,temra=cd[cd$CD45RA>2&cd$CD62L<=2.2,,drop=FALSE])},
 CP8={cd<-dat[dat$CD3>2&dat$CD4>2,,drop=FALSE];list(all=dat,cd3cd4=cd)},
 CP10={p<-dat[dat$Viability<=1.5&dat$CD45>2,,drop=FALSE];list(all=dat,primary=p)},
 CP16={p<-dat[dat$DUMP<=1.5&dat$CD45>2,,drop=FALSE];list(all=dat,primary=p,apc=p[p$HLA_DR>1.5,,drop=FALSE])},
 CP22={b<-dat[dat$DUMP<=1.5&dat$CD19>2,,drop=FALSE];list(all=dat,bcells=b)},
 CP23={p<-dat[dat$DUMP<=1.5&dat$CD45>2,,drop=FALSE];list(all=dat,primary=p,cd14=p[p$CD14>1.5,,drop=FALSE])},
 CP24={cd<-dat[dat$CD3>2&dat$CD8>2.2,,drop=FALSE];list(all=dat,cd3cd8=cd)},
 CP25={d<-dat[dat$DUMP<=1.5,,drop=FALSE];cd<-d[d$CD3>2&d$CD4>2,,drop=FALSE];list(all=dat,dump_low=d,cd3cd4=cd,regulatory=cd[cd$CD25>1.5&cd$IL7RA<=1.5,,drop=FALSE])},
 CP26={p<-dat[dat$DUMP<=1.5&dat$CD45>2,,drop=FALSE];list(all=dat,primary=p,nk_like=p[p$CD56>1.5|p$CD16>1.5,,drop=FALSE])},
 CP28={d<-dat[dat$DUMP<=1.5,,drop=FALSE];c3<-d[d$CD3>2,,drop=FALSE];list(all=dat,dump_low=d,cd3=c3)} )
for(n in names(z$min))if(nrow(parents[[n]])<z$min[[n]])stop("Insufficient events in parent ",n)
write.csv(data.frame(population=names(parents),events=vapply(parents,nrow,integer(1))),file.path(out,paste0(panel,"_parent_event_counts.csv")),row.names=FALSE)

pal<-c("#E8F3FA","#A9CEE6","#5A99C7","#1F79A8","#28B4BE","#38C96B","#E8DF3F","#F29A24","#D82420")
lims<-function(x,t){q<-quantile(x[is.finite(x)],c(.001,.999),na.rm=TRUE);r<-range(c(q,t));s<-diff(r);if(!is.finite(s)||s<=0)s<-1;r+c(-.07,.07)*s}
points2d<-function(d,x,y,xl,yl){d<-d[is.finite(d[[x]])&is.finite(d[[y]])&d[[x]]>=xl[1]&d[[x]]<=xl[2]&d[[y]]>=yl[1]&d[[y]]<=yl[2],c(x,y),drop=FALSE];if(nrow(d)>80000)d<-d[sample.int(nrow(d),80000),,drop=FALSE];k<-tryCatch(MASS::kde2d(d[[x]],d[[y]],n=160,lims=c(xl,yl)),error=function(e)NULL);if(is.null(k))d$.den<-1 else{ix<-pmax(1,pmin(length(k$x),findInterval(d[[x]],k$x)));iy<-pmax(1,pmin(length(k$y),findInterval(d[[y]],k$y)));d$.den<-k$z[cbind(ix,iy)]};d[order(d$.den),,drop=FALSE]}
qplot1<-function(s,id){d<-parents[[s$parent]][,c(s$x,s$y),drop=FALSE];d<-d[complete.cases(d),,drop=FALSE];x<-d[[s$x]];y<-d[[s$y]];n<-length(x);pc<-c(sum(x<=s$xt&y>s$yt),sum(x>s$xt&y>s$yt),sum(x<=s$xt&y<=s$yt),sum(x>s$xt&y<=s$yt))*100/n;xl<-lims(x,s$xt);yl<-lims(y,s$yt);pp<-points2d(d,s$x,s$y,xl,yl);lab<-data.frame(x=c(xl[1]+.18*diff(xl),xl[1]+.82*diff(xl),xl[1]+.18*diff(xl),xl[1]+.82*diff(xl)),y=c(yl[1]+.83*diff(yl),yl[1]+.83*diff(yl),yl[1]+.17*diff(yl),yl[1]+.17*diff(yl)),label=sprintf("%s\n%.1f%%",s$labels,pc));p<-ggplot(pp,aes(x=.data[[s$x]],y=.data[[s$y]]))+ggrastr::geom_point_rast(aes(color=.den),size=.17,alpha=.7,raster.dpi=600)+scale_color_gradientn(colours=pal,guide="none")+geom_vline(xintercept=s$xt,linetype="dashed")+geom_hline(yintercept=s$yt,linetype="dashed")+geom_label(data=lab,aes(x=x,y=y,label=label),inherit.aes=FALSE,size=2.2,fontface="bold",fill="white",linewidth=.2)+coord_cartesian(xlim=xl,ylim=yl,expand=FALSE)+labs(title=paste0(id,". ",s$title),subtitle=paste0(s$parent_text,"\n",format(n,big.mark=",")," events"),x=s$xlab,y=s$ylab)+theme_classic(base_size=9)+theme(plot.title=element_text(face="bold"),plot.subtitle=element_text(size=7),axis.title=element_text(face="bold",size=8));list(plot=p,stats=data.frame(panel=id,title=s$title,type="quadrant",parent=s$parent_text,n=n,pct1=pc[1],pct2=pc[2],pct3=pc[3],pct4=pc[4]))}
hplot1<-function(s,id){x<-parents[[s$parent]][[s$x]];x<-x[is.finite(x)];n<-length(x);pos<-100*mean(x>s$t);p<-ggplot(data.frame(x=x),aes(x=x))+geom_histogram(aes(y=after_stat(density)),bins=80,fill="#DCECF7",color="grey35")+geom_density(color="#1F79A8")+geom_vline(xintercept=s$t,linetype="dashed")+labs(title=paste0(id,". ",s$title),subtitle=paste0(s$parent_text,"\n",format(n,big.mark=",")," events"),x=s$xlab,y="Density")+theme_classic(base_size=9)+theme(plot.title=element_text(face="bold"),plot.subtitle=element_text(size=7),axis.title=element_text(face="bold",size=8));list(plot=p,stats=data.frame(panel=id,title=s$title,type="histogram",parent=s$parent_text,n=n,pct1=100-pos,pct2=pos,pct3=NA_real_,pct4=NA_real_))}
res<-lapply(seq_along(z$plots),function(i)if(z$plots[[i]]$type=="q")qplot1(z$plots[[i]],LETTERS[i])else hplot1(z$plots[[i]],LETTERS[i]));write.csv(do.call(rbind,lapply(res,`[[`,"stats")),file.path(out,paste0(panel,"_plot_percentages.csv")),row.names=FALSE)
fig<-patchwork::wrap_plots(lapply(res,`[[`,"plot"),ncol=z$ncol)+patchwork::plot_annotation(title=z$title,subtitle="Representative local FCS sample; fixed logicle transformation; dashed lines denote analytical thresholds",caption=z$cap)
prefix<-paste0("Supplementary_Figure_",z$S,"_",panel,"_flow_atlas")
ggsave(file.path(out,paste0(prefix,".pdf")),fig,width=z$w,height=z$h,device=if(capabilities("cairo"))cairo_pdf else "pdf")
ggsave(file.path(out,paste0(prefix,"_600dpi.png")),fig,width=z$w,height=z$h,dpi=600,bg="white")
message(panel," supplementary flow atlas completed: ",out)
