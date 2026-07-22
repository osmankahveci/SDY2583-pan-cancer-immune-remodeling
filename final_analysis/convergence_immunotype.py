from __future__ import annotations
import argparse,json
from itertools import combinations
import numpy as np,pandas as pd,patsy
import statsmodels.formula.api as smf
from scipy import stats
from .common import PANEL_ORDER,PRINCIPAL,IMMUNOTYPES,bh,master,paths,residual

def cors(r,cohort,method):
    rows=[]
    for a,b in combinations(PANEL_ORDER,2):
        x,y=r[a],r[b]; ok=x.notna()&y.notna(); n=int(ok.sum())
        if n<20: e=p=np.nan
        elif method=="spearman": e,p=stats.spearmanr(x[ok],y[ok])
        else: e,p=stats.pearsonr(x[ok],y[ok])
        rows.append(dict(cohort=cohort,method=method,panel_1=a,panel_2=b,score_1=PRINCIPAL[a],score_2=PRINCIPAL[b],n_pair=n,estimate=e,p_value=p))
    t=pd.DataFrame(rows); t["q_fdr"]=bh(t.p_value); return t

def wald_p(fit):
    names=list(fit.params.index); ids=[i for i,n in enumerate(names) if n.startswith("C(immunotype)")]
    R=np.zeros((len(ids),len(names)))
    for j,i in enumerate(ids): R[j,i]=1
    return float(fit.wald_test(R,scalar=True).pvalue)

def marginal(fit,d,cluster):
    nd=d.copy(); nd["immunotype"]=cluster
    X=np.asarray(patsy.build_design_matrices([fit.model.data.design_info],nd,return_type="dataframe")[0]); x=X.mean(0)
    est=float(x@fit.params.to_numpy()); se=float(np.sqrt(x@np.asarray(fit.cov_params())@x)); return est,est-1.96*se,est+1.96*se

def models(d,name,cancer=False,treatment=False):
    oms=[]; means=[]
    for panel in PANEL_ORDER:
        score=PRINCIPAL[panel]; cols=[score,"immunotype","age","sex"]; terms=["C(immunotype)","age","C(sex)"]
        if cancer: cols+= ["cancer_subgroup_adjusted"]; terms+=["C(cancer_subgroup_adjusted)"]
        else: cols+=["disease"]; terms+=["C(disease)"]
        if treatment: cols+=["therapy_status"]; terms+=["C(therapy_status)"]
        w=d[cols].dropna(); fit=smf.ols(f"Q('{score}') ~ "+" + ".join(terms),w).fit(cov_type="HC3")
        oms.append(dict(subset=name,panel=panel,score=score,n=len(w),robust_wald_p=wald_p(fit),r_squared=fit.rsquared))
        for g in IMMUNOTYPES:
            if g in set(w.immunotype):
                est,lo,hi=marginal(fit,w,g); means.append(dict(subset=name,panel=panel,score=score,immunotype=g,adjusted_mean=est,ci95_low=lo,ci95_high=hi,n_model=len(w)))
    o=pd.DataFrame(oms); o["q_fdr_across_10"]=bh(o.robust_wald_p); return o,pd.DataFrame(means)

def run(matrix=None,out=None):
    m,o=paths(matrix,out); d=master(m); cancer=d.disease.eq("Cancer patient"); healthy=d.disease.eq("Healthy control")
    residuals={}
    for label,sub,covs in [("all",d,["age","sex","disease"]),("cancer",d[cancer],["age","sex"]),("healthy",d[healthy],["age","sex"])]:
        r=pd.DataFrame(index=d.index)
        for p in PANEL_ORDER: r[p]=residual(sub,PRINCIPAL[p],covs).reindex(d.index)
        residuals[label]=r
    tables=[]
    for label,r in residuals.items(): tables += [cors(r,label,"spearman"),cors(r,label,"pearson")]
    cor=pd.concat(tables,ignore_index=True); cor.to_csv(o/"cross_panel_principal_correlations.csv",index=False)
    fs=cor.query("cohort=='all' and method=='spearman'").copy()
    cs=cor.query("cohort=='cancer' and method=='spearman'")[["panel_1","panel_2","estimate"]].rename(columns={"estimate":"estimate_cancer"})
    hs=cor.query("cohort=='healthy' and method=='spearman'")[["panel_1","panel_2","estimate"]].rename(columns={"estimate":"estimate_healthy"})
    rob=fs.merge(cs,on=["panel_1","panel_2"]).merge(hs,on=["panel_1","panel_2"]); rob["same_direction_cancer"]=np.sign(rob.estimate)==np.sign(rob.estimate_cancer); rob["same_direction_healthy"]=np.sign(rob.estimate)==np.sign(rob.estimate_healthy); rob.to_csv(o/"cross_panel_direction_robustness.csv",index=False)
    common=d.dropna(subset=["age","sex","disease","immunotype","cohort",*PRINCIPAL.values()])
    om=[]; mm=[]
    for name,sub in [("full",common),("training",common[common.cohort.eq("TRAINING")]),("validation",common[common.cohort.eq("VALIDATION")])]:
        a,b=models(sub,name); om.append(a); mm.append(b)
    cc=common[common.disease.eq("Cancer patient")]; a,b=models(cc,"cancer_only",True); om.append(a); mm.append(b)
    ct=cc[cc.therapy_status.notna()&cc.therapy_status.ne("no_treatment_data")]; a,b=models(ct,"cancer_treatment_adjusted",True,True); om.append(a); mm.append(b)
    pd.concat(om).to_csv(o/"immunotype_principal_score_omnibus.csv",index=False); pd.concat(mm).to_csv(o/"immunotype_adjusted_marginal_means.csv",index=False)
    summary={"n_master":len(d),"n_healthy":int(healthy.sum()),"n_cancer":int(cancer.sum()),"spearman_fdr":int((fs.q_fdr<.05).sum()),"pearson_fdr":int((cor.query("cohort=='all' and method=='pearson'").q_fdr<.05).sum()),"significant_spearman_same_direction_both_strata":int(((rob.q_fdr<.05)&rob.same_direction_cancer&rob.same_direction_healthy).sum())}
    (o/"convergence_immunotype_summary.json").write_text(json.dumps(summary,indent=2)); print(json.dumps(summary,indent=2))

def main():
    p=argparse.ArgumentParser(); p.add_argument("--matrix"); p.add_argument("--output-dir"); a=p.parse_args(); run(a.matrix,a.output_dir)
if __name__=="__main__": main()
