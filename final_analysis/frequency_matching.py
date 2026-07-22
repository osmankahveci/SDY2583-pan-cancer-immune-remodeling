from __future__ import annotations
import argparse,json,os
from pathlib import Path
import numpy as np,pandas as pd,statsmodels.api as sm,statsmodels.formula.api as smf
from scipy import stats
from statsmodels.othermod.betareg import BetaModel
from .common import bh
SPECS=[("CD8 differentiation/checkpoint remodeling","CP24","pct_temra_like","TEMRA-like CD8 cells within CD3+CD8+ (%)"),("CD4 helper/regulatory-like remodeling","CP8","pct_cd25pos_il7ralow_treg_like","CD25+IL7RA-low regulatory-like cells (%)"),("B-cell/humoral repatterning","CP22","pct_switched_memory_like_within_b","Switched-memory-like cells within B cells (%)"),("NK-cell attenuation","CP26","pct_nk_like_within_cd45_dump_low","NK-like cells within CD45+ dump-low events (%)"),("T/NK-interface remodeling","CP28","pct_cd8_within_cd3","CD8+ cells within CD3+ cells (%)"),("Myeloid/granulocytic remodeling","CP10","pct_cd13_cd66b_gran_like_within_cd45","CD13+CD66b+ granulocyte-like cells within CD45+ (%)"),("APC/DC-like depletion","CP16","pct_hladr_apc_core_within_cd45_dump_low","HLA-DR+ APC-core within CD45+ dump-low events (%)"),("Monocyte/macrophage-like remodeling","CP23","pct_cd45_dump_low_within_total","CD45+ dump-low events within total events (%)")]
PANELS=["CP24","CP8","CP10","CP16","CP22","CP23","CP26","CP28"]
def req(x,label):
    if not x: raise SystemExit(f"Missing {label}")
    p=Path(x).expanduser().resolve()
    if not p.exists(): raise SystemExit(f"{label} not found: {p}")
    return p
def panel(path,p,master):
    d=pd.read_csv(path)
    if p=="CP24": d=d.merge(master[["subject_id","CP24_subject_id_original_panel"]].rename(columns={"subject_id":"master_id","CP24_subject_id_original_panel":"subject_id"}),on="subject_id",how="left")
    else: d=d.rename(columns={"subject_id":"master_id"})
    return d
def fit(d,y):
    x=d[[y,"cancer","age","sex"]].dropna().rename(columns={y:"pct"}); x=x[x.pct.between(0,100)].copy(); n=len(x)
    lin=smf.ols("pct~cancer+age+C(sex)",x).fit(cov_type="HC3"); p=x.pct/100; x["ya"]=np.clip((p*(n-1)+.5)/n,1e-8,1-1e-8); x["elogit"]=np.log(x.ya/(1-x.ya)); el=smf.ols("elogit~cancer+age+C(sex)",x).fit(cov_type="HC3")
    status="beta_regression_logit_HC3"
    try: be=BetaModel.from_formula("ya~cancer+age+C(sex)",x).fit(disp=False,maxiter=2000,cov_type="HC3")
    except Exception: status="fractional_logit_HC3_fallback"; be=smf.glm("ya~cancer+age+C(sex)",x,family=sm.families.Binomial()).fit(cov_type="HC3")
    return dict(n=n,healthy_mean_pct=x.loc[x.cancer.eq(0),"pct"].mean(),cancer_mean_pct=x.loc[x.cancer.eq(1),"pct"].mean(),linear_beta_percentage_points=lin.params.cancer,linear_hc3_se=lin.bse.cancer,linear_p=lin.pvalues.cancer,empirical_logit_beta=el.params.cancer,empirical_logit_hc3_se=el.bse.cancer,empirical_logit_p=el.pvalues.cancer,beta_regression_beta=be.params.cancer,beta_regression_hc3_se=be.bse.cancer,beta_regression_p=be.pvalues.cancer,beta_model_status=status)
def smd(a,b): return float((a.mean()-b.mean())/np.sqrt((a.var(ddof=1)+b.var(ddof=1))/2))
def bsmd(a,b): return float((a-b)/np.sqrt((a*(1-a)+b*(1-b))/2))
def match(d,cal):
    ca=d[d.cancer.eq(1)].sort_values(["age","master_id"]); he=d[d.cancer.eq(0)]; avail=set(he.master_id); rows=[]
    for _,c in ca.iterrows():
        z=he[he.sex.eq(c.sex)&he.master_id.isin(avail)].copy(); z["gap"]=(z.age-c.age).abs(); z=z[z.gap<=cal].sort_values(["gap","age","master_id"])
        if z.empty: continue
        h=z.iloc[0]; avail.remove(h.master_id); rows.append(dict(pair_id=len(rows)+1,cancer_id=c.master_id,healthy_id=h.master_id,sex=c.sex,cancer_age=c.age,healthy_age=h.age,signed_age_difference=c.age-h.age,absolute_age_difference=abs(c.age-h.age)))
    return pd.DataFrame(rows)
def diag(base,pairs,cal):
    ca=base[base.cancer.eq(1)]; he=base[base.cancer.eq(0)]; f=(pairs.sex=="Female").mean(); g=pairs.absolute_age_difference
    return dict(analysis=f"same_sex_nearest_age_{cal}y",caliper_years=cal,matched_pairs=len(pairs),matched_subjects=2*len(pairs),unmatched_cancer=len(ca)-len(pairs),unmatched_healthy=len(he)-len(pairs),cancer_mean_age=pairs.cancer_age.mean(),healthy_mean_age=pairs.healthy_age.mean(),age_smd=smd(pairs.cancer_age,pairs.healthy_age),female_proportion_cancer=f,female_proportion_healthy=f,sex_smd=0,mean_signed_age_difference=pairs.signed_age_difference.mean(),mean_absolute_age_difference=g.mean(),median_absolute_age_difference=g.median(),q25_absolute_age_difference=g.quantile(.25),q75_absolute_age_difference=g.quantile(.75),p90_absolute_age_difference=g.quantile(.9),p95_absolute_age_difference=g.quantile(.95),maximum_absolute_age_difference=g.max())
def paired(pairs,tables,label):
    rows=[]
    for axis,p,y,name in SPECS:
        v=tables[p][["master_id",y]]; a=pairs[["pair_id","cancer_id"]].merge(v,left_on="cancer_id",right_on="master_id")[["pair_id",y]].rename(columns={y:"a"}); b=pairs[["pair_id","healthy_id"]].merge(v,left_on="healthy_id",right_on="master_id")[["pair_id",y]].rename(columns={y:"b"}); dif=a.merge(b,on="pair_id").dropna(); x=dif.a-dif.b; n=len(x); mean=x.mean(); se=x.std(ddof=1)/np.sqrt(n); pv=2*stats.t.sf(abs(mean/se),n-1); tc=stats.t.ppf(.975,n-1); rows.append(dict(caliper=label,axis=axis,panel=p,outcome=y,label=name,n_pairs=n,mean_paired_difference_percentage_points=mean,paired_se=se,ci95_low=mean-tc*se,ci95_high=mean+tc*se,paired_p=pv))
    r=pd.DataFrame(rows); r["paired_q_fdr"]=bh(r.paired_p); return r
def npct(n,d): return f"{int(n)} ({100*n/d:.1f}%)"
def table1(m):
    h=m[m.disease_group_model.eq("Healthy control")]; c=m[m.disease_group_model.eq("Cancer patient")]; rows=[]
    def add(s,k,a,b,z): rows.append(dict(section=s,characteristic=k,overall_n850=a,healthy_controls_n408=b,cancer_patients_n442=z))
    def age(x,mean):
        v=x.age_for_clinical_model.dropna(); return f"{v.mean():.1f} ({v.std(ddof=1):.1f})" if mean else f"{v.median():.1f} [{v.quantile(.25):.1f}, {v.quantile(.75):.1f}]"
    add("Cohort","Participants, n","850","408","442")
    for label,fn in [("Valid age for adjusted models, n (%)",lambda x:x.age_for_clinical_model.notna().sum()),("Invalid or missing age, n (%)",lambda x:x.age_for_clinical_model.isna().sum())]: add("Demographics",label,npct(fn(m),len(m)),npct(fn(h),len(h)),npct(fn(c),len(c)))
    add("Demographics","Age, years, mean (SD)",age(m,1),age(h,1),age(c,1)); add("Demographics","Age, years, median [IQR]",age(m,0),age(h,0),age(c,0))
    for v in ["Female","Male","Not Specified"]: add("Demographics",f"Sex: {v}, n (%)",npct(m.sex_for_clinical_model.eq(v).sum(),len(m)),npct(h.sex_for_clinical_model.eq(v).sum(),len(h)),npct(c.sex_for_clinical_model.eq(v).sum(),len(c)))
    for v,n in c.cancer_subgroup_model.value_counts().items(): add("Cancer subgroups",f"{v}, n (% of cancer cohort)","—","—",npct(n,len(c)))
    known=c.therapy_status_model.ne("no_treatment_data"); nk=known.sum(); add("Treatment at sampling","Treatment metadata available, n (%)","—","—",npct(nk,len(c))); add("Treatment at sampling","Ongoing active treatment, n (% of treatment-known)","—","—",npct(c.therapy_status_model.eq("ongoing_active_treatment").sum(),nk)); add("Treatment at sampling","No ongoing therapy, n (% of treatment-known)","—","—",npct(c.therapy_status_model.eq("no_ongoing_therapy").sum(),nk)); add("Treatment at sampling","Treatment metadata unavailable, n (%)","—","—",npct(c.therapy_status_model.eq("no_treatment_data").sum(),len(c)))
    for lab,col in [("Chemotherapy","chemotherapy_01"),("Targeted therapy","targeted_therapy_01"),("Any immunotherapy","any_immunotherapy_01"),("Immune-checkpoint inhibitor","ici_immunotherapy_01"),("Endocrine/hormonal therapy","endocrine_hormonal_01"),("Antibody-drug conjugate","adc_01"),("Experimental therapy","experimental_01"),("Radiotherapy","radiotherapy_01")]: add("Treatment modalities",f"{lab}, n (% of treatment-known)","—","—",npct(c[col].eq(1).sum(),nk))
    nl=c.therapy_line_group.notna().sum(); add("Treatment line","Treatment line available, n (%)","—","—",npct(nl,len(c))); add("Treatment line","First line, n (% of line-known)","—","—",npct(c.therapy_line_group.eq("first_line").sum(),nl)); add("Treatment line","Later line, n (% of line-known)","—","—",npct(c.therapy_line_group.eq("later_line").sum(),nl)); add("Treatment line","Treatment line unavailable, n (%)","—","—",npct(c.therapy_line_group.isna().sum(),len(c)))
    return pd.DataFrame(rows)
def run(args):
    mp=req(args.matrix,"ALL10 matrix"); m=pd.read_csv(mp); out=Path(args.output_dir).expanduser().resolve(); out.mkdir(parents=True,exist_ok=True); tables={p:panel(req(getattr(args,p.lower()),f"{p} features"),p,m) for p in PANELS}; meta=m[["subject_id","disease_group_model","age_for_clinical_model","sex_for_clinical_model"]].rename(columns={"subject_id":"master_id","age_for_clinical_model":"age","sex_for_clinical_model":"sex"}); meta["cancer"]=meta.disease_group_model.eq("Cancer patient").astype(int)
    for p,d in tables.items(): tables[p]=d.drop(columns=[x for x in ["age","sex","cancer","disease_group_model"] if x in d]).merge(meta,on="master_id",how="left")
    rows=[]
    for axis,p,y,label in SPECS: z=fit(tables[p],y); z.update(axis=axis,panel=p,outcome=y,label=label); rows.append(z)
    sen=pd.DataFrame(rows)
    for x in ["linear","empirical_logit","beta_regression"]: sen[x+"_q_fdr"]=bh(sen[x+"_p"])
    sen["direction_concordant"]=(np.sign(sen.linear_beta_percentage_points)==np.sign(sen.empirical_logit_beta))&(np.sign(sen.linear_beta_percentage_points)==np.sign(sen.beta_regression_beta)); sen["fdr_significant_all_models"]=(sen.linear_q_fdr<.05)&(sen.empirical_logit_q_fdr<.05)&(sen.beta_regression_q_fdr<.05); sen.to_csv(out/"selected_frequency_model_sensitivity.csv",index=False)
    base=meta[meta.age.notna()&meta.sex.isin(["Male","Female"])]; ca=base[base.cancer.eq(1)]; he=base[base.cancer.eq(0)]; fc=(ca.sex=="Female").mean(); fh=(he.sex=="Female").mean(); pre=dict(analysis="before_matching",caliper_years=np.nan,matched_pairs=np.nan,matched_subjects=len(base),unmatched_cancer=np.nan,unmatched_healthy=np.nan,cancer_mean_age=ca.age.mean(),healthy_mean_age=he.age.mean(),age_smd=smd(ca.age,he.age),female_proportion_cancer=fc,female_proportion_healthy=fh,sex_smd=bsmd(fc,fh),mean_signed_age_difference=np.nan,mean_absolute_age_difference=np.nan,median_absolute_age_difference=np.nan,q25_absolute_age_difference=np.nan,q75_absolute_age_difference=np.nan,p90_absolute_age_difference=np.nan,p95_absolute_age_difference=np.nan,maximum_absolute_age_difference=np.nan)
    p5,p10=match(base,5),match(base,10); pd.DataFrame([pre,diag(base,p5,5),diag(base,p10,10)]).to_csv(out/"age_matching_diagnostics.csv",index=False); pr=pd.concat([paired(p5,tables,"5_year"),paired(p10,tables,"10_year")]); pr.to_csv(out/"age_matched_paired_outcome_results.csv",index=False); table1(m).to_csv(out/"table1_cohort_characteristics.csv",index=False); pd.DataFrame([dict(panel=p,fcs_files_analyzed=847 if p=="CP16" else 850,unique_subjects=847 if p=="CP16" else 850) for p in ["CP7","CP8","CP10","CP16","CP22","CP23","CP24","CP25","CP26","CP28"]]).to_csv(out/"panel_fcs_availability.csv",index=False)
    s={"eight_outcomes_direction_concordant":int(sen.direction_concordant.sum()),"eight_outcomes_fdr_significant_all_models":int(sen.fdr_significant_all_models.sum()),"five_year_matched_pairs":len(p5),"ten_year_matched_pairs":len(p10),"five_year_outcomes_fdr_significant":int((pr.query("caliper=='5_year'").paired_q_fdr<.05).sum()),"ten_year_outcomes_fdr_significant":int((pr.query("caliper=='10_year'").paired_q_fdr<.05).sum())}; (out/"frequency_matching_table1_summary.json").write_text(json.dumps(s,indent=2)); print(json.dumps(s,indent=2))
def main():
    p=argparse.ArgumentParser(); p.add_argument("--matrix",default=os.getenv("SDY2583_ALL10_MATRIX","")); [p.add_argument("--"+x.lower(),default=os.getenv(f"SDY2583_{x}_FEATURES","")) for x in PANELS]; p.add_argument("--output-dir",default=os.getenv("SDY2583_FINAL_RESULTS_DIR","results/local-final")); run(p.parse_args())
if __name__=="__main__": main()
