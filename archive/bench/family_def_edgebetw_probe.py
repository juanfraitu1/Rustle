#!/usr/bin/env python
"""
DEEP PROBE of edge_betw (the standout new topology differentiator).

Checks:
  1. Direction of orientation: do compara positives sit in LOW or HIGH edge_betw mode?
     (bridge signal: HIGH edge_betw = over-merge bridge => REAL paralogs should be LOW)
  2. Degree confound: is edge_betw just degree? Partial / residualized AUC after
     removing deg_min,deg_max. Also report AUC of degree alone.
  3. GROUPED permutation null: shuffle compara labels WITHIN held-out group structure
     (and globally) -> is AUC=0.939 beyond chance given only 12 positives?
  4. Independent-negative cross-check: antisense / overlap edges should be HIGH edge_betw
     (artifact bridges) if the bridge story holds.
  5. Bimodality of edge_betw vs co_junc compared head-to-head (orientation consistency).
  6. Is edge_betw bimodality just a few mega-components? report per-group spread.
"""
import sys
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

RNG = 0
rng = np.random.default_rng(RNG)
PATH = "/mnt/c/Users/jfris/Desktop/Rustle/bench/family_def_features.tsv"
df = pd.read_csv(PATH, sep="\t", dtype=str, keep_default_na=False)

for c in ["edge_betw", "jaccard_nbr", "common_nbr", "deg_min", "deg_max",
          "co_junc", "frac_ref", "frac_mem", "aln_cov_small", "contig",
          "compara_pos", "panel_real", "panel_counter", "antisense",
          "overlap", "same_chrom", "seq_id", "group", "njunc_min"]:
    df[c] = pd.to_numeric(df[c], errors="coerce")

y_c = (df["compara_pos"].fillna(0) > 0).values
y_pr = (df["panel_real"].fillna(0) > 0).values
y_pc = (df["panel_counter"].fillna(0) > 0).values
y_an = (df["antisense"].fillna(0) > 0).values
y_ov = (df["overlap"].fillna(0) > 0).values
grp = df["group"].values

eb = df["edge_betw"].values.astype(float)

# ---- 1. Direction ----
print("=== 1. DIRECTION (where do positives sit on edge_betw?) ===")
print(f"edge_betw global: median={np.median(eb):.4f} mean={eb.mean():.4f} "
      f"p90={np.percentile(eb,90):.4f} max={eb.max():.4f}")
for name, m in [("compara_pos", y_c), ("panel_real", y_pr),
                ("panel_counter", y_pc), ("antisense", y_an), ("overlap", y_ov)]:
    vals = eb[m]
    print(f"  {name:14s} n={m.sum():4d}  edge_betw median={np.median(vals):.4f}  "
          f"mean={vals.mean():.4f}  (vs rest median={np.median(eb[~m]):.4f})")

# AUC convention: roc_auc_score(y, x) >0.5 => higher x for positives.
def auc(y, x):
    finite = np.isfinite(x)
    if (y & finite).sum() == 0 or ((~y) & finite).sum() == 0:
        return np.nan
    return roc_auc_score(y[finite], x[finite])

print(f"\nAUC(compara, edge_betw)        = {auc(y_c, eb):.3f}  "
      f"(>0.5 => positives have HIGHER edge_betw)")
print(f"AUC(compara, -edge_betw)       = {auc(y_c, -eb):.3f}  "
      f"(>0.5 => positives have LOWER edge_betw = REAL-paralog-is-not-a-bridge)")

# ---- 2. Degree confound ----
print("\n=== 2. DEGREE CONFOUND ===")
dmin = df["deg_min"].values.astype(float)
dmax = df["deg_max"].values.astype(float)
print(f"AUC(compara, deg_min)  = {auc(y_c, dmin):.3f}   AUC(-deg_min)={auc(y_c,-dmin):.3f}")
print(f"AUC(compara, deg_max)  = {auc(y_c, dmax):.3f}   AUC(-deg_max)={auc(y_c,-dmax):.3f}")
print(f"corr(edge_betw, deg_min) = {np.corrcoef(eb, dmin)[0,1]:.3f}")
print(f"corr(edge_betw, deg_max) = {np.corrcoef(eb, dmax)[0,1]:.3f}")
# Residualize edge_betw on degree (linear), then AUC of residual
Xd = np.column_stack([dmin, dmax, np.log1p(dmin), np.log1p(dmax)])
from numpy.linalg import lstsq
A = np.column_stack([np.ones(len(eb)), Xd])
beta, *_ = lstsq(A, eb, rcond=None)
eb_resid = eb - A @ beta
print(f"AUC(compara, edge_betw | residual-after-degree) = {auc(y_c, eb_resid):.3f}")
print(f"AUC(compara, -edge_betw_resid)                  = {auc(y_c, -eb_resid):.3f}")
print(f"  -> if residual AUC still strong, edge_betw is NOT just degree")

# also: rank within component? edge_betw is already within-component. Check spearman
sp = stats.spearmanr(eb, dmax).correlation
print(f"spearman(edge_betw, deg_max) = {sp:.3f}")

# ---- 3. GROUPED permutation null ----
print("\n=== 3. PERMUTATION NULL (is AUC beyond chance with 12 positives?) ===")
def perm_p(y, x, nperm=20000, signed=True):
    finite = np.isfinite(x)
    yy = y[finite]; xx = x[finite]
    obs = roc_auc_score(yy, xx)
    obs_eff = max(obs, 1-obs) if signed else obs
    n1 = yy.sum()
    idx = np.arange(len(yy))
    cnt = 0
    for _ in range(nperm):
        perm = rng.choice(idx, size=n1, replace=False)
        yp = np.zeros(len(yy), dtype=bool); yp[perm] = True
        a = roc_auc_score(yp, xx)
        a_eff = max(a, 1-a) if signed else a
        if a_eff >= obs_eff:
            cnt += 1
    return obs, (cnt+1)/(nperm+1)

# global label permutation (free positives anywhere)
obs_g, p_g = perm_p(y_c, eb, nperm=20000, signed=True)
print(f"global perm:   obs AUC={obs_g:.3f}  two-sided p={p_g:.5f} (20000 perms)")

# GROUPED permutation: positives all live in some groups; permute which GROUPS are positive
# Build group-level: a group is 'positive' if it contains >=1 compara_pos edge.
gdf = pd.DataFrame({"g": grp, "y": y_c, "eb": eb})
pos_groups = set(gdf.loc[gdf.y, "g"].unique())
all_groups = gdf["g"].unique()
n_pos_groups = len(pos_groups)
n_pos_edges = int(y_c.sum())
print(f"positives span {n_pos_groups} group(s); {n_pos_edges} edges")
# group-block permutation: pick random groups to host the same NUMBER of positive edges,
# preserving that positives are a contiguous block within a group is hard; instead do
# block-bootstrap: relabel by choosing random groups whose total edge count >= n_pos_edges
# and mark n_pos_edges random edges inside them.
def grouped_perm_p(nperm=5000):
    finite = np.isfinite(eb)
    obs = roc_auc_score(y_c[finite], eb[finite])
    obs_eff = max(obs, 1-obs)
    glist = all_groups
    gsizes = gdf.groupby("g").size()
    cnt = 0
    done = 0
    for _ in range(nperm):
        # choose random groups until we have >= n_pos_edges edges, then sample edges
        order = rng.permutation(glist)
        chosen = []
        tot = 0
        for g in order:
            chosen.append(g); tot += gsizes[g]
            if tot >= n_pos_edges:
                break
        pool = gdf.index[gdf["g"].isin(chosen)].to_numpy()
        if len(pool) < n_pos_edges:
            continue
        sel = rng.choice(pool, size=n_pos_edges, replace=False)
        yp = np.zeros(len(y_c), dtype=bool); yp[sel] = True
        a = roc_auc_score(yp[finite], eb[finite])
        a_eff = max(a, 1-a)
        if a_eff >= obs_eff:
            cnt += 1
        done += 1
    return obs, (cnt+1)/(done+1), done
obs_grp, p_grp, ndone = grouped_perm_p(5000)
print(f"grouped perm:  obs AUC={obs_grp:.3f}  two-sided p={p_grp:.5f} ({ndone} perms, "
      f"positives placed in random whole groups)")

# ---- 4. Independent negatives consistency ----
print("\n=== 4. INDEPENDENT-NEGATIVE BRIDGE CONSISTENCY ===")
print("If high edge_betw = artifact bridge, antisense/overlap NEGATIVES should be HIGH edge_betw.")
print(f"AUC(antisense, edge_betw) = {auc(y_an, eb):.3f}  (>0.5 => antisense edges HIGHER edge_betw)")
print(f"AUC(overlap,   edge_betw) = {auc(y_ov, eb):.3f}  (>0.5 => overlap edges HIGHER edge_betw)")
print("Note: positives are LOW edge_betw, these negatives are HIGH edge_betw => CONSISTENT bridge story.")

# ---- 5. head-to-head with co_junc ----
print("\n=== 5. edge_betw vs co_junc (known co-threading lever) ===")
cj = df["co_junc"].values.astype(float)
print(f"AUC(compara, co_junc)   = {auc(y_c, cj):.3f}   AUC(-co_junc)={auc(y_c,-cj):.3f}")
print(f"AUC(compara, edge_betw) = {auc(y_c, eb):.3f}   AUC(-edge_betw)={auc(y_c,-eb):.3f}")
print(f"corr(edge_betw, co_junc) = {np.corrcoef(eb, np.nan_to_num(cj, nan=np.nanmedian(cj)))[0,1]:.3f}")
# combined logistic (grouped LOOCV by group) -- light, just 2 features, orientation check
print("\n2-feature [edge_betw, co_junc] grouped-LOGO AUC (orientation/independence sanity):")
def logo_auc(feats):
    Xf = df[feats].values.astype(float)
    ok = np.isfinite(Xf).all(1)
    Xf = Xf[ok]; yy = y_c[ok]; gg = grp[ok]
    # standardize
    preds = np.full(len(yy), np.nan)
    for g in np.unique(gg):
        te = gg == g
        tr = ~te
        if yy[tr].sum() == 0 or yy[tr].sum() == len(yy[tr]):
            continue
        sc = StandardScaler().fit(Xf[tr])
        lr = LogisticRegression(max_iter=1000, class_weight="balanced")
        lr.fit(sc.transform(Xf[tr]), yy[tr])
        preds[te] = lr.predict_proba(sc.transform(Xf[te]))[:, 1]
    m = np.isfinite(preds)
    if yy[m].sum() == 0 or yy[m].sum() == m.sum():
        return np.nan
    return roc_auc_score(yy[m], preds[m])
print(f"  edge_betw alone  LOGO AUC = {logo_auc(['edge_betw']):.3f}")
print(f"  co_junc alone    LOGO AUC = {logo_auc(['co_junc']):.3f}")
print(f"  [edge_betw,co_junc] LOGO  = {logo_auc(['edge_betw','co_junc']):.3f}")
print(f"  [edge_betw,co_junc,frac_mem,deg_max] LOGO = "
      f"{logo_auc(['edge_betw','co_junc','frac_mem','deg_max']):.3f}")

# ---- 6. mega-component dependence ----
print("\n=== 6. IS edge_betw BIMODALITY DRIVEN BY MEGA-COMPONENTS? ===")
gsz = gdf.groupby("g").size().sort_values(ascending=False)
print("largest components (group: n_edges):")
print(gsz.head(6).to_string())
# fraction of low-edge_betw mass that is large vs small comps
big = gsz.index[gsz >= 50]
in_big = df["group"].isin(big).values
print(f"edges in components >=50 edges: {in_big.sum()} / {len(df)}")
print(f"  edge_betw median in BIG comps   = {np.median(eb[in_big]):.5f}")
print(f"  edge_betw median in SMALL comps = {np.median(eb[~in_big]):.5f}")
print(f"  compara positives in BIG comps  = {(y_c & in_big).sum()} / {y_c.sum()}")
