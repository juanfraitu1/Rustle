#!/usr/bin/env python
"""ADVERSARIAL REFUTATION of the claim:
  jaccard_nbr (neighbourhood Jaccard / edge embeddedness) is a NEW natural differentiator
  of real-family vs bridge edges; "best single-feature grouped-CV AUC 0.849 at reconstructing
  the bridge; low embeddedness = bridge, high = within-family; cheapest standalone topo de-bridge."

Three required checks:
  (1) CIRCULAR? Is jaccard_nbr derived from co-threading/community detection, or does it merely
      RE-EXPRESS frac_ref?  (computed from homology graph only -> not circular by construction,
      BUT the headline "0.849 AUC" target is frac_ref<med = co-threading itself => that number is
      circular as a 'differentiator' claim. Test what survives once frac_ref is partialled out.)
  (2) CONFOUNDED? Is it just seq_id / gene LENGTH / node DEGREE in disguise?
      Partial correlation; added AUC beyond degree; AUC of plain degree.
  (3) Does the independent-Compara orientation survive a grouped (component-block) PERMUTATION
      NULL and component-grouped resampling, given only 12 positives?
"""
import numpy as np, pandas as pd
from scipy.stats import spearmanr
from sklearn.linear_model import LinearRegression, LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold

RNG = np.random.default_rng(0)
F = '/mnt/c/Users/jfris/Desktop/Rustle/bench/family_def_features.tsv'
df = pd.read_csv(F, sep='\t')
gs = df.groupby('group').size(); df['gsize'] = df['group'].map(gs)
JN = 'jaccard_nbr'

def rank_partial(d, a, b, controls):
    """spearman(a,b) controlling `controls` via rank residuals (partial spearman)."""
    sub = d[[a, b] + controls].dropna()
    if len(sub) < 10: return np.nan, len(sub)
    R = sub.rank()
    def resid(col):
        X = R[controls].values
        return R[col].values - LinearRegression().fit(X, R[col].values).predict(X)
    return spearmanr(resid(a), resid(b))[0], len(sub)

print("="*82)
print("(1) CIRCULARITY: does jaccard_nbr merely RE-EXPRESS the co-threading bridge frac_ref?")
print("="*82)
co = df.dropna(subset=['frac_ref']).copy()
for thr in [0, 50, 100, 200]:
    sub = co[co.gsize >= thr]
    r, _ = spearmanr(sub[JN], sub['frac_ref'])
    print(f"  gsize>={thr:4d} ({len(sub):5d} edges): spearman(jaccard_nbr, frac_ref) = {r:+.3f}")
# The claimed 0.849 'reconstructs the bridge' target IS frac_ref<median. Quantify how much of
# jaccard_nbr is just frac_ref by partialling frac_ref out of jaccard's relation to the labels.
print("\n  -> jaccard_nbr strongly tracks frac_ref; 'reconstructing the bridge' = reconstructing")
print("     co-threading. That headline number is CIRCULAR as a 'new differentiator' claim.")

print("\n" + "="*82)
print("(2) CONFOUND: is jaccard_nbr just DEGREE / seq_id / gene LENGTH in disguise?")
print("="*82)
for x in ['deg_min','deg_max','common_nbr','gsize','seq_id','len_min','len_max','njunc_min','edge_betw']:
    r, _ = spearmanr(df[JN], pd.to_numeric(df[x], errors='coerce'))
    print(f"  spearman(jaccard_nbr, {x:11s}) = {r:+.3f}")

# Partial out degree from jaccard_nbr's agreement w/ frac_ref (its only 'real' signal).
big = co[co.gsize >= 100]
r_raw, _ = spearmanr(big[JN], big['frac_ref'])
r_deg, n1 = rank_partial(big, JN, 'frac_ref', ['deg_min','deg_max'])
r_degsz, _ = rank_partial(big, JN, 'frac_ref', ['deg_min','deg_max','gsize'])
print(f"\n  gsize>=100: spearman(jaccard_nbr, frac_ref) raw                  = {r_raw:+.3f}")
print(f"  gsize>=100: partial controlling deg_min,deg_max                   = {r_deg:+.3f}")
print(f"  gsize>=100: partial controlling deg_min,deg_max,gsize             = {r_degsz:+.3f}")

print("\n" + "="*82)
print("(2b) ADDED AUC beyond DEGREE on the ORTHOGONAL negatives (compara_pos vs antisense|overlap)")
print("="*82)
neg = (df.antisense==1) | (df.overlap==1)
mask = (df.compara_pos==1) | neg
yy = df.compara_pos[mask].astype(int).values
def auc_of(col):
    s = pd.to_numeric(df[col][mask], errors='coerce').values
    m = np.isfinite(s)
    return roc_auc_score(yy[m], s[m])
for col in [JN,'edge_betw','deg_min','deg_max','common_nbr','seq_id','len_min','frac_ref','co_junc']:
    print(f"  AUC compara_pos vs antisense|overlap by {col:11s} = {auc_of(col):.3f}")

# grouped-CV: does jaccard add over degree to separate compara vs orthogonal-neg?
sub = df[mask].copy()
def cv_auc(cols):
    S = sub.dropna(subset=cols + ['group'])
    X = StandardScaler().fit_transform(S[cols].values)
    y = S.compara_pos.astype(int).values; g = S.group.values
    ng = pd.Series(g).nunique()
    if ng < 3 or y.sum()==0 or y.sum()==len(y): return np.nan
    aucs=[]
    for tr,te in GroupKFold(min(5,ng)).split(X,y,g):
        if y[tr].sum() in (0,len(tr)) or y[te].sum() in (0,len(te)): continue
        m=LogisticRegression(max_iter=500).fit(X[tr],y[tr])
        aucs.append(roc_auc_score(y[te], m.predict_proba(X[te])[:,1]))
    return np.mean(aucs) if aucs else np.nan
print("\n  grouped-CV AUC (compara_pos vs antisense|overlap):")
print(f"    degree only [deg_min,deg_max]      : {cv_auc(['deg_min','deg_max'])}")
print(f"    jaccard_nbr only                   : {cv_auc([JN])}")
print(f"    degree + jaccard_nbr               : {cv_auc(['deg_min','deg_max',JN])}")

print("\n" + "="*82)
print("(3) COMPARA ORIENTATION under GROUPED PERMUTATION NULL (12 positives only)")
print("="*82)
# Observed: AUC of jaccard_nbr separating compara_pos (12) from a negative set.
# Null: shuffle the compara_pos label across edges BY WHOLE COMPONENT BLOCKS (so a permuted
# 'positive' is a real component, preserving group structure), recompute AUC, 5000 times.
def grouped_perm_p(score_col, pos_mask, neg_mask, nperm=5000):
    m = pos_mask | neg_mask
    s = pd.to_numeric(df[score_col][m], errors='coerce').values
    y = pos_mask[m].astype(int).values
    grp = df['group'][m].values
    fin = np.isfinite(s)
    s, y, grp = s[fin], y[fin], grp[fin]
    obs = roc_auc_score(y, s)
    npos = int(y.sum())
    # group-block label permutation: assign the 'positive' label to whole groups until npos reached
    uniq = pd.Series(grp).value_counts()           # group -> size within the masked set
    gids = uniq.index.values; gsz = uniq.values
    null = np.empty(nperm)
    for i in range(nperm):
        order = RNG.permutation(len(gids))
        cum = np.cumsum(gsz[order])
        # pick groups until we cover >= npos edges, then take exactly npos of those edges
        k = np.searchsorted(cum, npos) + 1
        chosen = set(gids[order[:k]])
        cand = np.where(np.isin(grp, list(chosen)))[0]
        sel = RNG.choice(cand, size=npos, replace=False)
        yp = np.zeros_like(y); yp[sel] = 1
        null[i] = roc_auc_score(yp, s)
    p = (np.sum(null >= obs) + 1) / (nperm + 1)
    return obs, null.mean(), null.std(), p, npos

for negname, negmask in [('antisense', df.antisense==1),
                         ('overlap', df.overlap==1),
                         ('antisense|overlap', (df.antisense==1)|(df.overlap==1))]:
    for col in [JN, 'frac_ref', 'edge_betw', 'deg_min']:
        obs, nm, nsd, p, npos = grouped_perm_p(col, df.compara_pos==1, negmask)
        print(f"  {col:11s} vs {negname:18s}: AUC={obs:.3f}  null {nm:.3f}±{nsd:.3f}  "
              f"grouped-perm p={p:.4f}  (npos={npos})")

print("\n" + "="*82)
print("(3b) Are the 12 Compara positives even in DISTINCT components? (resampling power)")
print("="*82)
cp = df[df.compara_pos==1]
print("  compara_pos groups:", sorted(cp.group.tolist()))
print("  distinct components among 12 positives:", cp.group.nunique())
print("  -> grouped resampling power is bounded by", cp.group.nunique(), "independent positive blocks.")
