#!/usr/bin/env python
"""A4 -- graph-topology deep dive: is a graph BRIDGE signal a NATURAL differentiator
of real multi-copy-family edges, independent of co-threading?

Label-free first (bimodality), then orient/null with sparse independent labels.
Tests whether edge_betw is just a function of degree / component size (confound)
or adds independent bridge signal that AGREES with the co-threading bridge (frac_ref).
Also evaluates clustering-coefficient / k-core as standalone de-bridge signals.
"""
import numpy as np, pandas as pd
from scipy.stats import spearmanr, mannwhitneyu
from sklearn.mixture import GaussianMixture
from sklearn.metrics import roc_auc_score, silhouette_score
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

RNG = np.random.default_rng(0)
F = '/mnt/c/Users/jfris/Desktop/Rustle/bench/family_def_features.tsv'
df = pd.read_csv(F, sep='\t')
gs = df.groupby('group').size(); df['gsize'] = df['group'].map(gs)

TOPO = ['edge_betw','jaccard_nbr','common_nbr','deg_min','deg_max']
COTHREAD = ['frac_ref','frac_mem','aln_cov_small','contig','co_junc']

def bimodality_coeff(x):
    x = x[np.isfinite(x)]
    if len(x) < 8: return np.nan
    from scipy.stats import skew, kurtosis
    n = len(x); g = skew(x); k = kurtosis(x, fisher=True)
    return (g**2 + 1) / (k + 3*(n-1)**2/((n-2)*(n-3)))

def gmm_bic_gap(x):
    """BIC(1-comp) - BIC(2-comp); >0 => 2 components preferred (more bimodal)."""
    x = x[np.isfinite(x)].reshape(-1,1)
    if len(x) < 20: return np.nan, np.nan
    xs = StandardScaler().fit_transform(x)
    b1 = GaussianMixture(1, random_state=0, n_init=2).fit(xs).bic(xs)
    b2 = GaussianMixture(2, random_state=0, n_init=2).fit(xs).bic(xs)
    g2 = GaussianMixture(2, random_state=0, n_init=2).fit(xs)
    lab = g2.predict(xs)
    sil = silhouette_score(xs, lab) if len(np.unique(lab))==2 and len(xs)<=20000 else np.nan
    return b1-b2, sil

print("="*78)
print("(LABEL-FREE) intrinsic bimodality of each feature on full unlabeled population")
print("="*78)
print(f"{'feature':14s} {'bimod_coeff':>11s} {'BIC1-BIC2':>11s} {'silhouette':>10s}  (BC>0.555 => bimodal-ish)")
for f in TOPO + COTHREAD:
    x = pd.to_numeric(df[f], errors='coerce').values
    bc = bimodality_coeff(x); gap, sil = gmm_bic_gap(x)
    print(f"{f:14s} {bc:11.3f} {gap:11.1f} {sil:10.3f}")

# ---- orient with labels: AUC of each topology feature vs label sets ----
print("\n" + "="*78)
print("(ORIENT) AUC of topology features for separating label classes")
print("="*78)
def auc_safe(y, s):
    s = pd.to_numeric(pd.Series(s), errors='coerce').values
    m = np.isfinite(s)
    if y[m].sum()==0 or y[m].sum()==len(y[m]): return np.nan
    return roc_auc_score(y[m], s[m])

# POSITIVES: compara_pos.  NEGATIVES: antisense / overlap (orthogonal artifacts).
labelsets = {
 'compara_pos_vs_antisense': (df.compara_pos==1, df.antisense==1),
 'compara_pos_vs_overlap':   (df.compara_pos==1, df.overlap==1),
 'compara_pos_vs_all_rest':  (df.compara_pos==1, ~((df.compara_pos==1))),
 'panel_real_vs_panel_counter': (df.panel_real==1, df.panel_counter==1),
}
print(f"{'feature':14s} " + " ".join(f"{k[:22]:>22s}" for k in labelsets))
for f in TOPO + ['frac_ref','co_junc']:
    row=[]
    for k,(pos,neg) in labelsets.items():
        mask = pos|neg
        y = pos[mask].astype(int).values
        a = auc_safe(y, df[f][mask].values)
        row.append(f"{a:22.3f}" if np.isfinite(a) else f"{'na':>22s}")
    print(f"{f:14s} " + " ".join(row))

# ---- confound test: edge_betw vs degree / component size ----
print("\n" + "="*78)
print("(CONFOUND) Is edge_betw just a function of degree / component size?")
print("="*78)
for x in ['deg_min','deg_max','common_nbr','gsize']:
    r,_ = spearmanr(df.edge_betw, df[x]); print(f"  spearman(edge_betw, {x:10s}) = {r:+.3f}")
r,_ = spearmanr(df.edge_betw, 1/df.gsize); print(f"  spearman(edge_betw, 1/gsize)    = {r:+.3f}")

# ---- AGREEMENT test: within large over-merged comps, does high edge_betw coincide
#      with low co-threading (frac_ref)?  If yes => topology independently confirms bridge.
print("\n" + "="*78)
print("(AGREEMENT) Within LARGE comps: edge_betw vs frac_ref (low frac_ref = co-thread bridge)")
print("="*78)
co = df.dropna(subset=['frac_ref']).copy()
for thr in [50, 100, 200]:
    big = co[co.gsize>=thr]
    if len(big)<30: continue
    r,_ = spearmanr(big.edge_betw, big.frac_ref)
    # if edge_betw HIGH should mean bridge => frac_ref LOW => negative correlation expected
    print(f"  gsize>={thr:4d} ({len(big):5d} edges): spearman(edge_betw, frac_ref) = {r:+.3f}"
          + ("  <-- AGREE (neg)" if r<-0.1 else "  <-- no agreement"))

# ---- betweenness BEYOND degree: does edge_betw add signal after partialling out degree? ----
print("\n" + "="*78)
print("(ADDED VALUE) Does edge_betw add bridge signal beyond degree? (partial corr w/ frac_ref)")
print("="*78)
def partial_spearman(a,b,controls,d):
    """spearman(a,b) controlling for `controls` via rank-residuals."""
    from sklearn.linear_model import LinearRegression
    sub = d[[a,b]+controls].dropna()
    R = sub.rank()
    def resid(col):
        X = R[controls].values; y=R[col].values
        return y - LinearRegression().fit(X,y).predict(X)
    ra, rb = resid(a), resid(b)
    return spearmanr(ra,rb)[0]
big = co[co.gsize>=100]
r_raw = spearmanr(big.edge_betw, big.frac_ref)[0]
r_par = partial_spearman('edge_betw','frac_ref',['deg_min','deg_max'], big)
print(f"  gsize>=100: raw  spearman(edge_betw, frac_ref)              = {r_raw:+.3f}")
print(f"  gsize>=100: partial (control deg_min,deg_max)               = {r_par:+.3f}")

# logistic: does edge_betw beat degree at predicting low-frac_ref (proxy bridge) within big comps?
big = big.copy(); big['bridge'] = (big.frac_ref < big.frac_ref.median()).astype(int)
def cv_auc(cols):
    from sklearn.model_selection import GroupKFold
    X = StandardScaler().fit_transform(big[cols].values); y=big.bridge.values; g=big.group.values
    aucs=[]
    gkf=GroupKFold(min(5, big.group.nunique()))
    for tr,te in gkf.split(X,y,g):
        if y[tr].sum() in (0,len(tr)) or y[te].sum() in (0,len(te)): continue
        m=LogisticRegression(max_iter=500).fit(X[tr],y[tr])
        aucs.append(roc_auc_score(y[te], m.predict_proba(X[te])[:,1]))
    return np.mean(aucs) if aucs else np.nan
print(f"  grouped-CV AUC predicting bridge(frac_ref<med) within gsize>=100:")
print(f"    degree only        [deg_min,deg_max]            : {cv_auc(['deg_min','deg_max']):.3f}")
print(f"    edge_betw only                                  : {cv_auc(['edge_betw']):.3f}")
print(f"    degree + edge_betw                              : {cv_auc(['deg_min','deg_max','edge_betw']):.3f}")
print(f"    jaccard_nbr only                                : {cv_auc(['jaccard_nbr']):.3f}")
print(f"    deg+betw+jaccard+common                         : {cv_auc(['deg_min','deg_max','edge_betw','jaccard_nbr','common_nbr']):.3f}")
print("  (NOTE: target frac_ref is itself co-threading => this only measures whether topology"
      " RECONSTRUCTS the co-thread bridge with NO alignment.)")
PY=None
