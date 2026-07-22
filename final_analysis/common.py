from __future__ import annotations
import os
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

SEED=2583
PANEL_ORDER=["CP7","CP8","CP10","CP16","CP22","CP23","CP24","CP25","CP26","CP28"]
PRINCIPAL={
"CP7":"CP7_integrated_checkpoint_remodeling_score",
"CP8":"CP8_integrated_CD4_helper_regulatory_remodeling_score",
"CP10":"CP10_integrated_myeloid_granulocytic_remodeling_score",
"CP16":"CP16_integrated_APC_DC_myeloid_remodeling_score",
"CP22":"CP22_integrated_humoral_B_cell_remodeling_score",
"CP23":"CP23_integrated_monocyte_macrophage_like_myeloid_remodeling_score",
"CP24":"CP24_integrated_remodeling_score",
"CP25":"CP25_integrated_CD4_regulatory_checkpoint_remodeling_score",
"CP26":"CP26_integrated_NK_remodeling_score",
"CP28":"CP28_integrated_TNK_interface_remodeling_score"}
IMMUNOTYPES=["G1_Naive","G2_Primed","G3_Progressive","G4_Chronic","G5_Suppressive"]

def bh(p):
    p=np.asarray(p,float); out=np.full(len(p),np.nan); ok=np.isfinite(p)
    if not ok.any(): return out
    order=np.argsort(p[ok]); v=p[ok][order]; n=len(v)
    q=np.minimum.accumulate((v*n/np.arange(1,n+1))[::-1])[::-1]
    out[np.where(ok)[0][order]]=np.minimum(q,1); return out

def paths(matrix=None,out=None):
    matrix=matrix or os.getenv("SDY2583_ALL10_MATRIX","")
    if not matrix: raise SystemExit("Set SDY2583_ALL10_MATRIX or pass --matrix")
    m=Path(matrix).expanduser().resolve(); o=Path(out or os.getenv("SDY2583_FINAL_RESULTS_DIR","results/local-final")).expanduser().resolve()
    if not m.exists(): raise SystemExit(f"Matrix not found: {m}")
    o.mkdir(parents=True,exist_ok=True); return m,o

def master(path):
    d=pd.read_csv(path).rename(columns={"age_for_clinical_model":"age","sex_for_clinical_model":"sex","disease_group_model":"disease","original_cluster":"immunotype","therapy_status_model":"therapy_status","cancer_subgroup_model":"cancer_subgroup_adjusted"})
    missing=[c for c in ["subject_id","age","sex","disease","immunotype","cohort",*PRINCIPAL.values()] if c not in d]
    if missing: raise SystemExit(f"Missing columns: {missing}")
    d["immunotype"]=pd.Categorical(d["immunotype"],categories=IMMUNOTYPES,ordered=True); return d

def residual(d,y,covars):
    w=d[[y,*covars]].dropna(); r=pd.Series(np.nan,index=d.index,dtype=float)
    terms=[c if pd.api.types.is_numeric_dtype(w[c]) else f"C({c})" for c in covars if w[c].nunique()>1]
    if len(w)>=20 and w[y].nunique()>1: r.loc[w.index]=smf.ols(f"Q('{y}') ~ "+(" + ".join(terms) or "1"),w).fit().resid
    return r

def complete_z(d):
    cols=["subject_id","age","sex","disease","immunotype","cohort",*PRINCIPAL.values()]
    c=d[cols].dropna().copy(); r=pd.DataFrame(index=c.index)
    for p in PANEL_ORDER: r[p]=residual(c,PRINCIPAL[p],["age","sex","disease"])
    c=c.loc[r.dropna().index]; r=r.loc[c.index]
    return c,(r-r.mean())/r.std(ddof=1)
