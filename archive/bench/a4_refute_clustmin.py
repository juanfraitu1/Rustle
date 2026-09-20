#!/usr/bin/env python
"""ADVERSARIAL REFUTATION of the claim:
   clust_min (min endpoint clustering coeff) is a NEW, label-free, alignment-free
   de-bridge differentiator because it AGREES with co-threading frac_ref at
   Spearman +0.682 (+0.600 deg-controlled), grouped-perm p=0.0005.

Three required checks:
  (1) CIRCULAR / re-expresses frac_ref? The agreement target frac_ref IS co-threading.
      Does clust_min separate the INDEPENDENT label (compara vs antisense|overlap)?
      Does the 'agreement' just mean it reconstructs the lever we already have?
  (2) CONFOUNDED with seq_id / gene LENGTH / node DEGREE? Partial Spearman of
      clust_min vs frac_ref controlling those; AUC-added-beyond-degree.
  (3) Does the 12-positive Compara orientation survive a GROUPED permutation null
      and component-grouped resampling?
"""
import numpy as np, pandas as pd, networkx as nx
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score
from sklearn.linear_model import LinearRegression, LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import GroupKFold

RNG = np.random.default_rng(0)
F='/mnt/c/Users/jfris/Desktop/Rustle/bench/family_def_features.tsv'
df=pd.read_csv(F,sep='\t')
gs=df.groupby('group').size(); df['gsize']=df['group'].map(gs)

# Rebuild clustering/k-core natively (confirm reproduction)
G=nx.Graph()
for a,b in zip(df.a,df.b): G.add_edge(a,b)
clust=nx.clustering(G); core=nx.core_number(G)
df['clust_min']=[min(clust[a],clust[b]) for a,b in zip(df.a,df.b)]
df['clust_max']=[max(clust[a],clust[b]) for a,b in zip(df.a,df.b)]
df['core_min']=[min(core[a],core[b]) for a,b in zip(df.a,df.b)]

co=df.dropna(subset=['frac_ref']).copy()
big=co[co.gsize>=100].copy()

print("#"*78)
print("# CHECK 0: reproduce the claimed agreement number")
print("#"*78)
r=spearmanr(big.clust_min, big.frac_ref)[0]
print(f"  spearman(clust_min, frac_ref) gsize>=100 = {r:+.3f}   (claim: +0.682)")

print("\n"+"#"*78)
print("# CHECK 1: CIRCULAR? agreement target frac_ref IS a co-threading feature.")
print("#  -> 'agreement' = reconstructing the EXISTING lever, not a NEW differentiator.")
print("#  -> the only non-circular test is the INDEPENDENT label.")
print("#"*78)
# Independent negatives orthogonal to co-threading+topology
neg=(df.antisense==1)|(df.overlap==1)
mask=(df.compara_pos==1)|neg
y=df.compara_pos[mask].astype(int).values
print("  INDEPENDENT label = compara_pos(+) vs antisense|overlap(-):")
for f in ['clust_min','clust_max','core_min','jaccard_nbr','frac_ref','co_junc']:
    s=df[f][mask].values; m=np.isfinite(s)
    print(f"    AUC {f:11s} = {roc_auc_score(y[m],s[m]):.3f}")
print("  -> If clust_min ~0.5 it does NOT distinguish real-paralog from artifact;")
print("     it only re-expresses frac_ref (which is high-AUC) within big comps.")

# Direct circularity: is clust_min essentially frac_ref re-expressed?
print("\n  Direct: how much of clust_min is explained by co-threading features?")
ct=['frac_ref','frac_mem','aln_cov_small','contig','co_junc']
sub=co.dropna(subset=ct).copy()
Xc=StandardScaler().fit_transform(sub[ct].values)
from sklearn.linear_model import LinearRegression as LR
r2=LR().fit(Xc, sub.clust_min.values).score(Xc, sub.clust_min.values)
print(f"    R^2(clust_min ~ co-threading feats) = {r2:.3f}  (high => circular re-expression)")

print("\n"+"#"*78)
print("# CHECK 2: CONFOUNDED with seq_id / LENGTH / DEGREE?")
print("#"*78)
def partial_spearman(a,b,controls,d):
    s=d[[a,b]+controls].dropna(); R=s.rank()
    def res(c):
        X=R[controls].values; yv=R[c].values
        return yv-LinearRegression().fit(X,yv).predict(X)
    return spearmanr(res(a),res(b))[0]
for ctrl,name in [(['deg_min','deg_max'],'degree'),
                  (['len_min','len_max'],'cDNA length'),
                  (['seq_id'],'seq_id'),
                  (['deg_min','deg_max','common_nbr'],'degree+common_nbr')]:
    r=partial_spearman('clust_min','frac_ref',ctrl,big)
    print(f"  partial spearman(clust_min,frac_ref | {name:18s}) = {r:+.3f}")
print("  (claim says +0.600 after deg control; but note frac_ref target = co-threading)")

# Does clust_min separate the INDEPENDENT label BEYOND degree? grouped CV
print("\n  Independent-label grouped CV AUC (compara vs antisense|overlap):")
sub=df[mask].copy(); sub['y']=sub.compara_pos.astype(int)
def cv_auc(cols):
    s=sub.dropna(subset=cols)
    X=StandardScaler().fit_transform(s[cols].values); yv=s.y.values; g=s.group.values
    if yv.sum()==0: return np.nan
    aucs=[]; gkf=GroupKFold(min(5,s.group.nunique()))
    for tr,te in gkf.split(X,yv,g):
        if yv[tr].sum() in (0,len(tr)) or len(np.unique(yv[te]))<2: continue
        m=LogisticRegression(max_iter=500).fit(X[tr],yv[tr])
        aucs.append(roc_auc_score(yv[te],m.predict_proba(X[te])[:,1]))
    return (np.mean(aucs),len(aucs))
for cols in [['deg_min','deg_max'],['clust_min'],['deg_min','deg_max','clust_min']]:
    res=cv_auc(cols); print(f"    {str(cols):40s} -> {res}")

print("\n"+"#"*78)
print("# CHECK 3: does 12-positive Compara orientation survive GROUPED perm null + resampling?")
print("#"*78)
# (3a) AUC of clust_min for compara vs antisense|overlap, grouped-label permutation null.
s=df[mask].dropna(subset=['clust_min']).copy()
yv=s.compara_pos.astype(int).values; sc=s.clust_min.values; grp=s.group.values
obs_auc=roc_auc_score(yv,sc)
# grouped permutation: shuffle the POSITIVE label among whole groups (a group is all-pos or not).
# Build group-level label then reassign.
gpos=s.groupby('group').compara_pos.max()  # which groups carry a positive
gids=gpos.index.values; glab=gpos.values
npos_groups=int(glab.sum())
null=[]
for _ in range(5000):
    perm=RNG.permutation(len(gids))
    pg=set(gids[perm[:npos_groups]])
    yperm=np.array([1 if g in pg else 0 for g in grp])
    if yperm.sum()==0 or yperm.sum()==len(yperm): continue
    null.append(roc_auc_score(yperm,sc))
null=np.array(null)
p=(np.sum(null>=obs_auc)+1)/(len(null)+1)
print(f"  clust_min AUC (compara vs antisense|overlap) obs = {obs_auc:.3f}")
print(f"    grouped-label perm null mean {null.mean():.3f} sd {null.std():.3f}  p(one-sided)={p:.4f}")
print(f"    n positive groups = {npos_groups}; effective sample for null is GROUPS not edges")

# (3b) compare to frac_ref/co_junc on same independent label + null
for f in ['frac_ref','co_junc']:
    s2=s.dropna(subset=[f]); sc2=s2[f].values
    yv2=s2.compara_pos.astype(int).values; grp2=s2.group.values
    gpos2=s2.groupby('group').compara_pos.max(); gids2=gpos2.index.values
    npg2=int(gpos2.values.sum())
    o=roc_auc_score(yv2,sc2)
    nl=[]
    for _ in range(5000):
        perm=RNG.permutation(len(gids2)); pg=set(gids2[perm[:npg2]])
        yp=np.array([1 if g in pg else 0 for g in grp2])
        if yp.sum()==0 or yp.sum()==len(yp): continue
        nl.append(roc_auc_score(yp,sc2))
    nl=np.array(nl); pp=(np.sum(nl>=o)+1)/(len(nl)+1)
    print(f"  {f:9s} AUC = {o:.3f}  grouped-perm p={pp:.4f}  (n pos-groups={npg2})")

print("\n"+"#"*78)
print("# CHECK 4: the AGREEMENT-with-frac_ref grouped null, on clust_min (claim p=0.0005)")
print("#  but reframed: is this agreement anything beyond reconstructing co-threading?")
print("#"*78)
sub=big[['clust_min','frac_ref','deg_min','deg_max','group']].dropna().copy()
R=sub[['clust_min','frac_ref','deg_min','deg_max']].rank()
Xc=R[['deg_min','deg_max']].values
def resid(yv): return yv-LinearRegression().fit(Xc,yv).predict(Xc)
ra=resid(R.clust_min.values); rb=resid(R.frac_ref.values)
obs=spearmanr(ra,rb)[0]
groups=sub.group.values; uq=np.unique(groups)
null=[]
for _ in range(2000):
    pg=RNG.permutation(uq); gmap={g:p for g,p in zip(uq,pg)}
    order=np.argsort([gmap[g] for g in groups],kind='stable')
    null.append(spearmanr(ra,rb[order])[0])
null=np.array(null); p=(np.sum(np.abs(null)>=abs(obs))+1)/(len(null)+1)
print(f"  partial spearman(clust_min,frac_ref|deg) obs={obs:+.3f} grouped-perm p={p:.4f}")
print(f"    n groups in gsize>=100 set = {len(uq)}  (few big over-merged comps => few units)")
print("  VERDICT NOTE: a significant p here only certifies clust_min RECONSTRUCTS frac_ref,")
print("  which is the EXISTING co-threading lever -> NOT a new independent differentiator.")
