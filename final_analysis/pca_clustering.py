from __future__ import annotations
import argparse,json
import numpy as np,pandas as pd
from sklearn.cluster import AgglomerativeClustering,KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_rand_score,calinski_harabasz_score,davies_bouldin_score,normalized_mutual_info_score,silhouette_score
from sklearn.mixture import GaussianMixture
from sklearn.neighbors import NearestNeighbors
from .common import PANEL_ORDER,PRINCIPAL,IMMUNOTYPES,SEED,master,paths,complete_z

def orient(load,scores):
    if load[:,0].mean()<0: load[:,0]*=-1; scores[:,0]*=-1
    my=[PANEL_ORDER.index(x) for x in ["CP10","CP16","CP23","CP26"]]; tn=[PANEL_ORDER.index(x) for x in ["CP7","CP24","CP25","CP28"]]
    if load[my,1].mean()-load[tn,1].mean()<0: load[:,1]*=-1; scores[:,1]*=-1
    return load,scores

def pca_fit(d):
    c,z=complete_z(d); p=PCA(n_components=len(PANEL_ORDER),svd_solver="full"); s=p.fit_transform(z); l,s=orient(p.components_.T.copy(),s); return c,z,p,l,s

def hopkins(x,f=.1):
    rng=np.random.default_rng(SEED); n,k=x.shape; m=max(10,min(n-1,round(n*f))); idx=rng.choice(n,m,False); lo=x.min(0); hi=x.max(0); u=rng.uniform(lo,hi,(m,k)); nn=NearestNeighbors(n_neighbors=2).fit(x); du=nn.kneighbors(u,1,return_distance=True)[0][:,0]; dw=nn.kneighbors(x[idx],2,return_distance=True)[0][:,1]; return float(du.sum()/(du.sum()+dw.sum()))

def run(matrix=None,out=None):
    m,o=paths(matrix,out); d=master(m); c,z,p,l,s=pca_fit(d)
    load=[]
    for j in range(2):
        for i,panel in enumerate(PANEL_ORDER): load.append(dict(component=f"PC{j+1}",panel=panel,score=PRINCIPAL[panel],loading=l[i,j],explained_variance_ratio=p.explained_variance_ratio_[j]))
    pd.DataFrame(load).to_csv(o/"pca_principal_score_loadings.csv",index=False)
    sf=pd.DataFrame(s[:,:2],columns=["PC1","PC2"],index=c.index); cen=[]
    for g in IMMUNOTYPES:
        for pc in ["PC1","PC2"]:
            v=sf.loc[c.immunotype.eq(g),pc]; se=v.std(ddof=1)/np.sqrt(len(v)); cen.append(dict(immunotype=g,component=pc,n=len(v),centroid=v.mean(),ci95_low=v.mean()-1.96*se,ci95_high=v.mean()+1.96*se))
    ct=pd.DataFrame(cen); ct.to_csv(o/"pca_immunotype_centroids.csv",index=False)
    fullc=ct.pivot(index="immunotype",columns="component",values="centroid").loc[IMMUNOTYPES]; rep=[]
    for subset in ["TRAINING","VALIDATION"]:
        sc,sz,sp,sl,ss=pca_fit(d[d.cohort.eq(subset)])
        identity=abs(l[:,0]@sl[:,0])+abs(l[:,1]@sl[:,1]); swapped=abs(l[:,0]@sl[:,1])+abs(l[:,1]@sl[:,0])
        if swapped>identity: sl[:,[0,1]]=sl[:,[1,0]]; ss[:,[0,1]]=ss[:,[1,0]]
        for j in range(2):
            if l[:,j]@sl[:,j]<0: sl[:,j]*=-1; ss[:,j]*=-1
        subcen=np.array([[ss[sc.immunotype.eq(g).to_numpy(),j].mean() for j in range(2)] for g in IMMUNOTYPES])
        for j in range(2): rep.append(dict(subset=subset,component=f"PC{j+1}",n=len(sc),loading_correlation_with_full=np.corrcoef(l[:,j],sl[:,j])[0,1],centroid_correlation_with_full=np.corrcoef(fullc.iloc[:,j],subcen[:,j])[0,1],explained_variance_ratio=sp.explained_variance_ratio_[j]))
    pd.DataFrame(rep).to_csv(o/"pca_training_validation_replication.csv",index=False)
    x=z.to_numpy(); original=pd.Categorical(c.immunotype,categories=IMMUNOTYPES).codes; rows=[]
    for k in range(2,9):
        for alg,labels in [("kmeans",KMeans(k,n_init=100,random_state=SEED).fit_predict(x)),("ward",AgglomerativeClustering(k,linkage="ward").fit_predict(x))]:
            rows.append(dict(algorithm=alg,k=k,n=len(x),silhouette=silhouette_score(x,labels),calinski_harabasz=calinski_harabasz_score(x,labels),davies_bouldin=davies_bouldin_score(x,labels),adjusted_rand_vs_original_immunotype=adjusted_rand_score(original,labels),normalized_mutual_information_vs_original_immunotype=normalized_mutual_info_score(original,labels),bic=np.nan,aic=np.nan))
    for k in range(1,7):
        gm=GaussianMixture(k,covariance_type="full",random_state=SEED,n_init=20).fit(x); lab=gm.predict(x); valid=k>1 and len(np.unique(lab))>1
        rows.append(dict(algorithm="gaussian_mixture",k=k,n=len(x),silhouette=silhouette_score(x,lab) if valid else np.nan,calinski_harabasz=calinski_harabasz_score(x,lab) if valid else np.nan,davies_bouldin=davies_bouldin_score(x,lab) if valid else np.nan,adjusted_rand_vs_original_immunotype=adjusted_rand_score(original,lab),normalized_mutual_information_vs_original_immunotype=normalized_mutual_info_score(original,lab),bic=gm.bic(x),aic=gm.aic(x)))
    cl=pd.DataFrame(rows); cl.to_csv(o/"clustering_stability_metrics.csv",index=False)
    summary={"n_complete":len(c),"pc1_explained_variance":p.explained_variance_ratio_[0],"pc2_explained_variance":p.explained_variance_ratio_[1],"hopkins_statistic":hopkins(x),"maximum_silhouette_non_gmm":cl.query("algorithm!='gaussian_mixture'").silhouette.max()}
    (o/"pca_clustering_summary.json").write_text(json.dumps(summary,indent=2)); print(json.dumps(summary,indent=2))

def main():
    p=argparse.ArgumentParser(); p.add_argument("--matrix"); p.add_argument("--output-dir"); a=p.parse_args(); run(a.matrix,a.output_dir)
if __name__=="__main__": main()
