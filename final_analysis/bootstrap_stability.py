from __future__ import annotations
import argparse,json,os,re
from pathlib import Path
import numpy as np,pandas as pd
from scipy import stats
from .common import PANEL_ORDER,PRINCIPAL,IMMUNOTYPES,bh

def ols(x,y,k):
    inv=np.linalg.pinv(x.T@x); b=inv@x.T@y; e=y-x@b; h=np.sum((x@inv)*x,axis=1); z=x@inv[:,k]; v=np.sum((z[:,None]**2)*(e**2)/np.maximum(1-h[:,None],1e-8)**2,axis=0); se=np.sqrt(v); p=2*stats.t.sf(abs(b[k]/se),x.shape[0]-x.shape[1]); return b[k],se,p
def rpca(x,y):
    e=y-x@(np.linalg.pinv(x)@y); z=(e-e.mean(0))/e.std(0,ddof=1); _,s,vt=np.linalg.svd(z,full_matrices=False); ev=(s**2)/(len(z)-1); return vt.T,ev/ev.sum(),z@vt.T
def pct(x,q): return float(np.quantile(x,q))
def run(matrix,out,resamples=2000,seed=2583):
    path=Path(matrix or os.getenv("SDY2583_ALL10_MATRIX","")).expanduser().resolve(); out=Path(out or os.getenv("SDY2583_FINAL_RESULTS_DIR","results/local-final")).expanduser().resolve(); out.mkdir(parents=True,exist_ok=True)
    if not path.exists(): raise SystemExit(f"Matrix not found: {path}")
    d=pd.read_csv(path); scores=[c for c in d if re.match(r"^CP\d+_.*_score$",c)]
    if len(scores)!=66: raise SystemExit(f"Expected 66 score columns; found {len(scores)}")
    pcs=[PRINCIPAL[p] for p in PANEL_ORDER]; need=["subject_id","age_for_clinical_model","sex_for_clinical_model","disease_group_model","original_cluster",*scores]; c=d[need].dropna().copy()
    if len(c)!=829: raise SystemExit(f"Expected complete n=829; found {len(c)}")
    c["cancer"]=c.disease_group_model.eq("Cancer patient").astype(float); sd=pd.get_dummies(c.sex_for_clinical_model,drop_first=True,dtype=float); x=np.column_stack([np.ones(len(c)),c.age_for_clinical_model,c.cancer,sd]); y66=c[scores].to_numpy(float); idx={v:i for i,v in enumerate(scores)}; ids=[idx[v] for v in pcs]; y10=y66[:,ids]; cl=c.original_cluster.astype(str).to_numpy(); cancer=c.cancer.astype(int).to_numpy(); hi=np.where(cancer==0)[0]; ci=np.where(cancer==1)[0]
    rl,refev,rs=rpca(x,y10)
    if rl[:,0].mean()<0: rl[:,0]*=-1; rs[:,0]*=-1
    my=[PANEL_ORDER.index(p) for p in ["CP10","CP16","CP23","CP26"]]; tn=[PANEL_ORDER.index(p) for p in ["CP7","CP24","CP25","CP28"]]
    if rl[my,1].mean()-rl[tn,1].mean()<0: rl[:,1]*=-1; rs[:,1]*=-1
    rc=np.array([[rs[cl==g,j].mean() for j in range(2)] for g in IMMUNOTYPES]); b10,se10,p10=ols(x,y10,2); q10=bh(p10); b66,se66,p66=ols(x,y66,2); q66=bh(p66); B=resamples; rng=np.random.default_rng(seed)
    bl=np.empty((B,10,2)); be=np.empty((B,2)); cos=np.empty((B,2)); sub=np.empty((B,2)); pos=np.empty(B,bool); pole=np.empty(B,bool); cc=np.empty((B,2)); bb10=np.empty((B,10)); bp10=np.empty((B,10)); bq10=np.empty((B,10)); bb66=np.empty((B,66)); bp66=np.empty((B,66)); bq66=np.empty((B,66))
    for z in range(B):
        ix=np.r_[rng.choice(hi,len(hi),True),rng.choice(ci,len(ci),True)]; rng.shuffle(ix); xb=x[ix]; l,e,s=rpca(xb,y10[ix]); ident=abs(rl[:,0]@l[:,0])+abs(rl[:,1]@l[:,1]); swap=abs(rl[:,0]@l[:,1])+abs(rl[:,1]@l[:,0])
        if swap>ident: l[:,[0,1]]=l[:,[1,0]]; s[:,[0,1]]=s[:,[1,0]]; e[[0,1]]=e[[1,0]]
        for j in range(2):
            if rl[:,j]@l[:,j]<0: l[:,j]*=-1; s[:,j]*=-1
        bl[z]=l[:,:2]; be[z]=e[:2]; cos[z]=[rl[:,j]@l[:,j] for j in range(2)]; sv=np.linalg.svd(rl[:,:2].T@l[:,:2],compute_uv=False); sub[z]=[sv.min(),sv.mean()]; pos[z]=np.all(l[:,0]>0); pole[z]=l[my,1].mean()-l[tn,1].mean()>0; scl=cl[ix]; cen=np.array([[s[scl==g,j].mean() for j in range(2)] for g in IMMUNOTYPES]); cc[z]=[np.corrcoef(rc[:,j],cen[:,j])[0,1] for j in range(2)]; a,_,p=ols(xb,y10[ix],2); bb10[z]=a; bp10[z]=p; bq10[z]=bh(p); a,_,p=ols(xb,y66[ix],2); bb66[z]=a; bp66[z]=p; bq66[z]=bh(p)
    glob=[dict(metric="bootstrap_resamples",reference_value=B,bootstrap_median=np.nan,bootstrap_95pct_low=np.nan,bootstrap_95pct_high=np.nan),dict(metric="common_complete_case_n",reference_value=len(c),bootstrap_median=np.nan,bootstrap_95pct_low=np.nan,bootstrap_95pct_high=np.nan)]
    for name,ref,v in [("pc1_explained_variance",refev[0],be[:,0]),("pc2_explained_variance",refev[1],be[:,1]),("pc1_loading_cosine_congruence",1,cos[:,0]),("pc2_loading_cosine_congruence",1,cos[:,1]),("two_component_minimum_canonical_correlation",1,sub[:,0]),("two_component_mean_canonical_correlation",1,sub[:,1]),("pc1_immunotype_centroid_correlation",1,cc[:,0]),("pc2_immunotype_centroid_correlation",1,cc[:,1])]: glob.append(dict(metric=name,reference_value=ref,bootstrap_median=pct(v,.5),bootstrap_95pct_low=pct(v,.025),bootstrap_95pct_high=pct(v,.975)))
    glob += [dict(metric="proportion_all_positive_pc1_loadings",reference_value=pos.mean(),bootstrap_median=np.nan,bootstrap_95pct_low=np.nan,bootstrap_95pct_high=np.nan),dict(metric="proportion_retaining_pc2_biological_orientation",reference_value=pole.mean(),bootstrap_median=np.nan,bootstrap_95pct_low=np.nan,bootstrap_95pct_high=np.nan)]; pd.DataFrame(glob).to_csv(out/"bootstrap_global_summary.csv",index=False)
    lr=[]
    for j in range(2):
        for i,p in enumerate(PANEL_ORDER):
            v=bl[:,i,j]; lo,up=pct(v,.025),pct(v,.975); lr.append(dict(panel=p,score=PRINCIPAL[p],component=f"PC{j+1}",reference_loading=rl[i,j],bootstrap_median_loading=pct(v,.5),bootstrap_95pct_low=lo,bootstrap_95pct_high=up,same_sign_probability=np.mean(np.sign(v)==np.sign(rl[i,j])),bootstrap_interval_excludes_zero=(lo>0 or up<0)))
    pd.DataFrame(lr).to_csv(out/"bootstrap_pca_loading_stability.csv",index=False)
    pr=[]
    for i,p in enumerate(PANEL_ORDER):
        v=bb10[:,i]; lo,up=pct(v,.025),pct(v,.975); pr.append(dict(panel=p,score=PRINCIPAL[p],reference_beta=b10[i],reference_hc3_se=se10[i],reference_p=p10[i],reference_q_across_10=q10[i],bootstrap_median_beta=pct(v,.5),bootstrap_95pct_low=lo,bootstrap_95pct_high=up,direction_stability_probability=np.mean(np.sign(v)==np.sign(b10[i])),nominal_p_below_0_05_probability=np.mean(bp10[:,i]<.05),fdr_q_below_0_05_probability_across_10=np.mean(bq10[:,i]<.05),bootstrap_interval_excludes_zero=(lo>0 or up<0)))
    pt=pd.DataFrame(pr); pt.to_csv(out/"bootstrap_principal_score_effect_stability.csv",index=False)
    ar=[]
    for i,score in enumerate(scores):
        p=re.match(r"^(CP\d+)_",score).group(1); v=bb66[:,i]; lo,up=pct(v,.025),pct(v,.975); ar.append(dict(panel=p,score=score,is_principal_integrated_score=score in pcs,reference_beta=b66[i],reference_hc3_se=se66[i],reference_p=p66[i],reference_q_across_66=q66[i],bootstrap_median_beta=pct(v,.5),bootstrap_95pct_low=lo,bootstrap_95pct_high=up,direction_stability_probability=np.mean(np.sign(v)==np.sign(b66[i])),fdr_q_below_0_05_probability_across_66=np.mean(bq66[:,i]<.05),bootstrap_interval_excludes_zero=(lo>0 or up<0)))
    at=pd.DataFrame(ar); at.to_csv(out/"bootstrap_all66_component_effect_stability.csv",index=False)
    panel=[]
    for p in PANEL_ORDER:
        z=at[at.panel.eq(p)]; q=z[z.is_principal_integrated_score].iloc[0]; panel.append(dict(panel=p,number_of_scores=len(z),scores_direction_stability_at_least_0_95=(z.direction_stability_probability>=.95).sum(),scores_bootstrap_interval_excluding_zero=z.bootstrap_interval_excludes_zero.sum(),scores_fdr_significant_in_at_least_0_80_resamples=(z.fdr_q_below_0_05_probability_across_66>=.8).sum(),minimum_direction_stability_probability=z.direction_stability_probability.min(),minimum_fdr_significance_probability=z.fdr_q_below_0_05_probability_across_66.min(),principal_score_direction_stability_probability=q.direction_stability_probability,principal_score_fdr_significance_probability=q.fdr_q_below_0_05_probability_across_66))
    pd.DataFrame(panel).to_csv(out/"bootstrap_panel_component_summary.csv",index=False)
    s={"resamples":B,"seed":seed,"n_complete":len(c),"pc1_loading_cosine_median":pct(cos[:,0],.5),"pc2_loading_cosine_median":pct(cos[:,1],.5),"minimum_subspace_correlation_median":pct(sub[:,0],.5),"all66_direction_stability_at_least_0_95":int((at.direction_stability_probability>=.95).sum()),"all66_bootstrap_intervals_excluding_zero":int(at.bootstrap_interval_excludes_zero.sum()),"all66_fdr_significant_at_least_0_80_resamples":int((at.fdr_q_below_0_05_probability_across_66>=.8).sum()),"principal_intervals_excluding_zero":int(pt.bootstrap_interval_excludes_zero.sum())}; (out/"bootstrap_component_stability_summary.json").write_text(json.dumps(s,indent=2)); print(json.dumps(s,indent=2))
def main():
    p=argparse.ArgumentParser(); p.add_argument("--matrix",default=os.getenv("SDY2583_ALL10_MATRIX","")); p.add_argument("--output-dir",default=os.getenv("SDY2583_FINAL_RESULTS_DIR","results/local-final")); p.add_argument("--resamples",type=int,default=int(os.getenv("SDY2583_BOOTSTRAP_RESAMPLES","2000"))); p.add_argument("--seed",type=int,default=int(os.getenv("SDY2583_BOOTSTRAP_SEED","2583"))); a=p.parse_args(); run(a.matrix,a.output_dir,a.resamples,a.seed)
if __name__=="__main__": main()
