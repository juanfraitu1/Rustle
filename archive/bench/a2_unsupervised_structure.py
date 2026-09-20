#!/usr/bin/env python
"""
A2 -- UNSUPERVISED MULTIVARIATE STRUCTURE.

Goal: discover NATURAL differentiators that separate real multi-copy gene-family
edges from over-merge / bridge edges, label-light, using sklearn.

Method:
  1. Standardize all numeric features (median impute), excluding label/id/group cols.
  2. PCA -> variance explained + PC1/PC2 loadings.
  3. 2-cluster GaussianMixture and KMeans on standardized matrix.
  4. WITHOUT having used labels to cluster: are the 12 compara_pos edges enriched
     in one cluster (Fisher)? Are antisense/overlap artifact edges enriched in
     the OTHER? Centroid-difference feature ranking. Is any NEW (topology) feature
     among the top separators, or is the split entirely co-threading?
  5. Re-run EXCLUDING the co-threading columns to see what structure remains.

Plus a LABEL-FREE intrinsic bimodality pass per feature (bimodality coefficient,
2-comp vs 1-comp GMM by BIC, silhouette of the 2-cluster split) to rank single
features as natural differentiators, and orient/null-check with the sparse labels.
"""
import numpy as np, pandas as pd
from scipy import stats
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.mixture import GaussianMixture
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score, roc_auc_score

RNG = 0
np.random.seed(RNG)
PATH = '/mnt/c/Users/jfris/Desktop/Rustle/bench/family_def_features.tsv'
df = pd.read_csv(PATH, sep='\t')

ID_COLS    = ['a', 'b']
LABEL_COLS = ['compara_pos', 'panel_real', 'panel_counter', 'group']
COTHREAD   = ['frac_ref', 'frac_mem', 'aln_cov_small', 'contig', 'co_junc']
TOPOLOGY   = ['edge_betw', 'jaccard_nbr', 'common_nbr', 'deg_min', 'deg_max']

# All numeric feature columns (exclude ids + labels + group)
feat_cols = [c for c in df.columns if c not in ID_COLS + LABEL_COLS]
# antisense/overlap/disjoint/same_chrom are independent-negative *signals*; they are
# legitimate features (orthogonal to co-threading/topology). Keep them in the matrix
# but we will ALSO use them as independent negatives for orientation.
print("=== feature columns used (%d) ===" % len(feat_cols))
print(feat_cols)

X_raw = df[feat_cols].copy()
# median impute
med = X_raw.median(numeric_only=True)
X_imp = X_raw.fillna(med)
scaler = StandardScaler()
Xs = scaler.fit_transform(X_imp.values)
Xs_df = pd.DataFrame(Xs, columns=feat_cols, index=df.index)

# ---- independent label vectors ----
y_compara = df['compara_pos'].values.astype(bool)
y_preal   = df['panel_real'].values.astype(bool)
y_pcount  = df['panel_counter'].values.astype(bool)
y_anti    = (df['antisense'] == 1).values
y_over    = (df['overlap'] == 1).values
# combined independent-negative set (artifact edges): antisense OR overlap
y_neg     = y_anti | y_over
print("\nlabel counts: compara=%d panel_real=%d panel_counter=%d antisense=%d overlap=%d neg(any)=%d"
      % (y_compara.sum(), y_preal.sum(), y_pcount.sum(), y_anti.sum(), y_over.sum(), y_neg.sum()))

# ============================================================
# PART A -- LABEL-FREE INTRINSIC BIMODALITY per feature
# ============================================================
def bimodality_coefficient(x):
    x = x[np.isfinite(x)]
    n = len(x)
    if n < 8 or np.std(x) == 0:
        return np.nan
    g = stats.skew(x)
    k = stats.kurtosis(x, fisher=True)  # excess kurtosis
    # BC = (g^2 + 1) / (k + 3*(n-1)^2/((n-2)(n-3)))
    denom = k + 3.0 * (n - 1) ** 2 / ((n - 2) * (n - 3))
    return (g * g + 1.0) / denom

def gmm_bic_split(x):
    """Return (delta_BIC = BIC1 - BIC2, silhouette of 2-cluster split).
       delta_BIC>0 => 2 comps preferred."""
    x = x[np.isfinite(x)].reshape(-1, 1)
    if len(x) < 50 or np.std(x) == 0:
        return np.nan, np.nan
    xs = StandardScaler().fit_transform(x)
    try:
        g1 = GaussianMixture(1, random_state=RNG).fit(xs)
        g2 = GaussianMixture(2, random_state=RNG, n_init=2).fit(xs)
    except Exception:
        return np.nan, np.nan
    bic1, bic2 = g1.bic(xs), g2.bic(xs)
    lab = g2.predict(xs)
    sil = np.nan
    if len(np.unique(lab)) == 2:
        # subsample for silhouette speed
        idx = np.random.choice(len(xs), min(4000, len(xs)), replace=False)
        try:
            sil = silhouette_score(xs[idx], lab[idx])
        except Exception:
            sil = np.nan
    return bic1 - bic2, sil

print("\n" + "=" * 70)
print("PART A: LABEL-FREE INTRINSIC BIMODALITY (single feature)")
print("=" * 70)
rows = []
for c in feat_cols:
    x = X_imp[c].values.astype(float)
    bc = bimodality_coefficient(x)
    dbic, sil = gmm_bic_split(x)
    rows.append((c, bc, dbic, sil))
bim = pd.DataFrame(rows, columns=['feat', 'bimod_coef', 'dBIC_2v1', 'silhouette'])
# bimod_coef > 0.555 is the classic uniform threshold for bimodality
bim = bim.sort_values('bimod_coef', ascending=False)
print(bim.to_string(index=False, float_format=lambda v: f"{v:8.4f}"))
print("\n(BC>0.555 = bimodal-ish; dBIC>0 => 2-comp GMM beats 1-comp; higher silhouette = cleaner split)")

# orientation + permutation null for the top single features
def orient_and_null(feat, y, nperm=2000, name=""):
    x = X_imp[feat].values.astype(float)
    m = np.isfinite(x)
    x = x[m]; yy = y[m]
    if yy.sum() == 0 or yy.sum() == len(yy):
        return None
    try:
        auc = roc_auc_score(yy, x)
    except Exception:
        return None
    # permutation null: shuffle labels, recompute AUC
    obs_dev = abs(auc - 0.5)
    cnt = 0
    npos = int(yy.sum())
    for _ in range(nperm):
        perm = np.zeros(len(yy), bool)
        perm[np.random.choice(len(yy), npos, replace=False)] = True
        a = roc_auc_score(perm, x)
        if abs(a - 0.5) >= obs_dev:
            cnt += 1
    p = (cnt + 1) / (nperm + 1)
    return auc, p

print("\n--- orientation of top features vs independent labels (AUC, perm-p) ---")
top_feats = bim.head(8)['feat'].tolist()
# always include topology features explicitly
for tf in TOPOLOGY + ['edge_betw']:
    if tf not in top_feats:
        top_feats.append(tf)
top_feats = list(dict.fromkeys(top_feats))
for feat in top_feats:
    line = f"{feat:14s}"
    for lname, yv in [('compara', y_compara), ('neg(anti|over)', y_neg)]:
        r = orient_and_null(feat, yv, nperm=1000, name=lname)
        if r is None:
            line += f"  {lname}: n/a"
        else:
            auc, p = r
            line += f"  {lname}: AUC={auc:.3f} p={p:.3f}"
    print(line)

# ============================================================
# PART B -- PCA on full standardized matrix
# ============================================================
print("\n" + "=" * 70)
print("PART B: PCA (full standardized matrix, all features)")
print("=" * 70)
pca = PCA(random_state=RNG).fit(Xs)
ev = pca.explained_variance_ratio_
print("variance explained (first 8 PCs):", np.round(ev[:8], 4))
print("cumulative:", np.round(np.cumsum(ev[:8]), 4))
load = pd.DataFrame(pca.components_[:3].T, index=feat_cols, columns=['PC1', 'PC2', 'PC3'])
print("\nPC1 loadings (sorted):")
print(load['PC1'].sort_values(key=abs, ascending=False).round(3).to_string())
print("\nPC2 loadings (sorted):")
print(load['PC2'].sort_values(key=abs, ascending=False).round(3).to_string())

# ============================================================
# PART C -- 2-cluster GMM + KMeans on full matrix, label-blind
# ============================================================
def cluster_and_eval(Xmat, cols, tag):
    print("\n" + "=" * 70)
    print(f"PART C [{tag}]: 2-cluster GMM + KMeans (label-blind)")
    print("=" * 70)
    out = {}
    for name, model in [('GMM', GaussianMixture(2, random_state=RNG, n_init=5)),
                        ('KMeans', KMeans(2, random_state=RNG, n_init=10))]:
        lab = model.fit_predict(Xmat)
        idx = np.random.choice(len(Xmat), min(4000, len(Xmat)), replace=False)
        sil = silhouette_score(Xmat[idx], lab[idx])
        sizes = np.bincount(lab)
        # centroid difference => which features separate the clusters
        c0 = Xmat[lab == 0].mean(0)
        c1 = Xmat[lab == 1].mean(0)
        cdiff = pd.Series(c1 - c0, index=cols)
        # orient: which cluster is "real-like"? cluster with MORE compara positives
        comp0 = y_compara[lab == 0].sum(); comp1 = y_compara[lab == 1].sum()
        real_cluster = 1 if comp1 / max(sizes[1], 1) >= comp0 / max(sizes[0], 1) else 0
        # Fisher: compara enrichment in real_cluster
        in_rc = (lab == real_cluster)
        def fisher(yv):
            a = (yv & in_rc).sum(); b = (yv & ~in_rc).sum()
            c = (~yv & in_rc).sum(); d = (~yv & ~in_rc).sum()
            _, p = stats.fisher_exact([[a, b], [c, d]])
            return a, b, c, d, p
        print(f"\n[{name}] sizes={sizes.tolist()} silhouette={sil:.4f}  real_cluster={real_cluster}")
        for lname, yv in [('compara', y_compara), ('panel_real', y_preal),
                          ('panel_counter', y_pcount), ('antisense', y_anti),
                          ('overlap', y_over), ('neg(any)', y_neg)]:
            a, b, c, d, p = fisher(yv)
            frac_rc = a / (a + c) if (a + c) else 0
            frac_oc = b / (b + d) if (b + d) else 0
            print(f"   {lname:14s} in_realcluster={a:4d}/{a+c:5d} ({frac_rc:.4f})  "
                  f"other={b:4d}/{b+d:5d} ({frac_oc:.4f})  fisher_p={p:.2e}")
        # top separating features by |centroid diff| (in standardized units)
        cdiff_abs = cdiff.abs().sort_values(ascending=False)
        print(f"   top separating features (|centroid diff|, std units), sign=>real_cluster side:")
        for f in cdiff_abs.head(12).index:
            sign = cdiff[f] if real_cluster == 1 else -cdiff[f]
            kind = ('COTHREAD' if f in COTHREAD else 'TOPOLOGY' if f in TOPOLOGY
                    else 'other')
            print(f"      {f:14s} |Δ|={cdiff_abs[f]:6.3f}  real-side={'+' if sign>0 else '-'}  [{kind}]")
        out[name] = (lab, real_cluster, sil, cdiff)
    return out

res_full = cluster_and_eval(Xs, feat_cols, "ALL FEATURES")

# ============================================================
# PART D -- EXCLUDE co-threading columns; what structure remains?
# ============================================================
nonct_cols = [c for c in feat_cols if c not in COTHREAD]
Xs_nonct = StandardScaler().fit_transform(X_imp[nonct_cols].values)
print("\n" + "#" * 70)
print("PART D: REMOVE co-threading features. Remaining cols:", nonct_cols)
print("#" * 70)
pca2 = PCA(random_state=RNG).fit(Xs_nonct)
print("variance explained (no-cothread, first 6 PCs):", np.round(pca2.explained_variance_ratio_[:6], 4))
load2 = pd.Series(pca2.components_[0], index=nonct_cols)
print("PC1 loadings (no-cothread):")
print(load2.sort_values(key=abs, ascending=False).round(3).to_string())
res_nonct = cluster_and_eval(Xs_nonct, nonct_cols, "NO CO-THREADING")

# ============================================================
# PART E -- TOPOLOGY-ONLY: is edge_betw a NEW differentiator,
# or confounded with degree? (degree controls)
# ============================================================
print("\n" + "#" * 70)
print("PART E: TOPOLOGY-ONLY structure + degree confound check")
print("#" * 70)
topo_cols = TOPOLOGY
Xs_topo = StandardScaler().fit_transform(X_imp[topo_cols].values)
# correlation matrix among topology features
corr = pd.DataFrame(Xs_topo, columns=topo_cols).corr()
print("topology feature correlation matrix:")
print(corr.round(3).to_string())

# edge_betw AUC vs labels, raw and residualized on degree (partial out deg_min,deg_max)
from sklearn.linear_model import LinearRegression
eb = X_imp['edge_betw'].values.astype(float)
degmat = X_imp[['deg_min', 'deg_max']].values.astype(float)
eb_resid = eb - LinearRegression().fit(degmat, eb).predict(degmat)
print("\nedge_betw vs labels (raw AUC | degree-residualized AUC):")
for lname, yv in [('compara', y_compara), ('neg(anti|over)', y_neg),
                  ('antisense', y_anti), ('overlap', y_over)]:
    if yv.sum() == 0:
        continue
    try:
        a_raw = roc_auc_score(yv, eb)
        a_res = roc_auc_score(yv, eb_resid)
        print(f"   {lname:14s} raw_AUC={a_raw:.3f}  resid_AUC={a_res:.3f}")
    except Exception as e:
        print(f"   {lname}: {e}")

# Does edge_betw separate compara from neg directly? (pos-vs-neg discrimination)
mask = y_compara | y_neg
yy = y_compara[mask]
print("\nedge_betw + topology, compara(+) vs artifact(-) discrimination (n_pos=%d n_neg=%d):"
      % (yy.sum(), (~yy).sum()))
for f in topo_cols:
    xf = X_imp[f].values.astype(float)[mask]
    try:
        a = roc_auc_score(yy, xf)
        print(f"   {f:14s} AUC={a:.3f} (1.0=high->real, 0.0=high->artifact)")
    except Exception as e:
        print(f"   {f}: {e}")

# ============================================================
# PART F -- grouped permutation sanity for the BEST single new feature
# (cluster membership shuffled within nothing; here we null the compara
#  enrichment of the cluster split holding cluster sizes fixed)
# ============================================================
print("\n" + "#" * 70)
print("PART F: grouped permutation null for cluster-compara enrichment (full matrix GMM)")
print("#" * 70)
lab, real_cluster, sil, cdiff = res_full['GMM']
in_rc = (lab == real_cluster)
# observed: fraction of compara in real cluster
obs = y_compara[in_rc].sum()
groups = df['group'].values
# permute compara labels by reassigning to random edges but PRESERVE group structure:
# shuffle group->label assignment. compara positives sit in specific groups; null = move
# the positive mass to random groups of similar membership.
uniq_groups = np.unique(groups)
pos_groups = np.unique(groups[y_compara])
nperm = 5000
cnt = 0
n_posgroups = len(pos_groups)
for _ in range(nperm):
    samp = np.random.choice(uniq_groups, n_posgroups, replace=False)
    permmask = np.isin(groups, samp)
    # within those groups, count how many edges land in real cluster, normalize by total compara count drawn
    # simpler: draw same number of edges as compara from these groups
    pool = np.where(permmask)[0]
    if len(pool) < y_compara.sum():
        draw = pool
    else:
        draw = np.random.choice(pool, y_compara.sum(), replace=False)
    if in_rc[draw].sum() >= obs:
        cnt += 1
p_group = (cnt + 1) / (nperm + 1)
print(f"observed compara-in-real-cluster = {obs}/{y_compara.sum()}; "
      f"grouped permutation p = {p_group:.4f}")

print("\nDONE.")
