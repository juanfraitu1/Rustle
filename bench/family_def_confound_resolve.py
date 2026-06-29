#!/usr/bin/env python
"""
RESOLVE the edge_betw confound: is its signal REAL bridge detection, or just a
proxy for COMPONENT SIZE (small comp -> high betweenness mechanically)?

Decisive tests:
  A. Component-size stratification of edge_betw and of the positives.
  B. AUC of edge_betw vs compara WITHIN big components only (where the over-merge
     blob lives) -- if it collapses, the signal was component-size, not bridge.
  C. The over-merge target: comp-238-like big same-strand blobs. Do the
     antisense/overlap NEGATIVES that live INSIDE big comps get separated by
     edge_betw? (This is the actual use-case: rank bridges WITHIN a blob.)
  D. Compare to a pure component-size feature: AUC(compara, -comp_size).
  E. njunc_min / len_min sanity: are they ALSO just component-size proxies?
"""
import sys
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

PATH = "/mnt/c/Users/jfris/Desktop/Rustle/bench/family_def_features.tsv"
df = pd.read_csv(PATH, sep="\t", dtype=str, keep_default_na=False)
for c in ["edge_betw", "deg_min", "deg_max", "common_nbr", "jaccard_nbr",
          "co_junc", "frac_mem", "njunc_min", "len_min", "len_max", "seq_id",
          "compara_pos", "panel_real", "panel_counter", "antisense", "overlap",
          "group"]:
    df[c] = pd.to_numeric(df[c], errors="coerce")

y_c = (df["compara_pos"].fillna(0) > 0).values
y_pc = (df["panel_counter"].fillna(0) > 0).values
y_an = (df["antisense"].fillna(0) > 0).values
y_ov = (df["overlap"].fillna(0) > 0).values
grp = df["group"].values
eb = df["edge_betw"].values.astype(float)

# component size in EDGES and component max degree (proxy for n nodes)
csize = df.groupby("group")["group"].transform("size").values.astype(float)
df["csize"] = csize

def auc(y, x, mask=None):
    if mask is None:
        mask = np.ones(len(x), bool)
    f = np.isfinite(x) & mask
    if (y & f).sum() == 0 or ((~y) & f).sum() == 0:
        return np.nan, int((y & f).sum()), int(((~y) & f).sum())
    return roc_auc_score(y[f], x[f]), int((y & f).sum()), int(((~y) & f).sum())

print("=== A. COMPONENT SIZE vs edge_betw (mechanical confound) ===")
print(f"corr(edge_betw, log comp_size)   = {np.corrcoef(eb, np.log(csize))[0,1]:.3f}")
print(f"spearman is captured; small comps -> high betweenness mechanically.")
for lo, hi, name in [(2,3,"size 2-3"),(4,9,"size 4-9"),(10,49,"size 10-49"),
                     (50,1e9,"size>=50")]:
    m = (csize>=lo)&(csize<=hi)
    print(f"  {name:10s} edges={m.sum():5d}  edge_betw median={np.median(eb[m]):.4f}  "
          f"compara={int((y_c&m).sum())}  panelC={int((y_pc&m).sum())}  "
          f"anti={int((y_an&m).sum())}  over={int((y_ov&m).sum())}")

print("\n=== B. edge_betw AUC WITHIN big components (signal survives size control?) ===")
big = csize >= 50
a_big, np1, nn1 = auc(y_c, eb, big)
print(f"compara in big comps: {int((y_c&big).sum())}  -> cannot AUC if 0 positives")
# Since 0 compara in big comps, use antisense/overlap negatives that DO live in big comps
a_an_big, p, n = auc(y_an, eb, big); print(f"AUC(antisense, edge_betw | big comps)  = {a_an_big}  (pos={p} neg={n})")
a_ov_big, p, n = auc(y_ov, eb, big); print(f"AUC(overlap,   edge_betw | big comps)  = {a_ov_big}  (pos={p} neg={n})")
# also small comps
small = csize < 50
a_an_sm, p, n = auc(y_an, eb, small); print(f"AUC(antisense, edge_betw | small comps)= {a_an_sm}  (pos={p} neg={n})")
a_ov_sm, p, n = auc(y_ov, eb, small); print(f"AUC(overlap,   edge_betw | small comps)= {a_ov_sm}  (pos={p} neg={n})")

print("\n=== C. THE ACTUAL USE-CASE: separate artifact bridges WITHIN big blobs ===")
print("Within big comps, do antisense/overlap (artifact) edges have HIGHER edge_betw")
print("than the same-strand co-located edges? (edge_betw>0.5 => >50% of big-comp shortest paths)")
for nm, y in [("antisense", y_an), ("overlap", y_ov)]:
    m = big
    a, p, n = auc(y, eb, m)
    pos_med = np.median(eb[m & y]) if (m&y).sum()>0 else np.nan
    neg_med = np.median(eb[m & ~y & np.isfinite(eb)])
    print(f"  {nm:10s} within big: AUC={a if a is not None else 'NA'}  "
          f"pos_med_eb={pos_med:.4f} vs neg_med_eb={neg_med:.4f}  (pos={p} neg={n})")

print("\n=== D. Pure component-size feature vs compara (is edge_betw just -comp_size?) ===")
a_sz, p, n = auc(y_c, -csize); print(f"AUC(compara, -comp_size)  = {a_sz:.3f}  (smaller comp => positive)")
a_eb, p, n = auc(y_c, eb);     print(f"AUC(compara,  edge_betw)  = {a_eb:.3f}")
# residualize edge_betw on log comp_size, re-AUC
from numpy.linalg import lstsq
A = np.column_stack([np.ones(len(eb)), np.log(csize)])
beta, *_ = lstsq(A, eb, rcond=None)
eb_r = eb - A@beta
a_ebr, p, n = auc(y_c, eb_r); print(f"AUC(compara,  edge_betw | resid log-size) = {a_ebr:.3f}")
a_dmaxr, _, _ = auc(y_c, -df['deg_max'].values.astype(float))
print(f"AUC(compara, -deg_max) = {a_dmaxr:.3f}  (deg_max also ~ comp size/density)")

print("\n=== E. Are njunc_min / len_min ALSO comp-size proxies, or independent? ===")
for f in ["njunc_min", "len_min", "len_max", "seq_id", "co_junc", "frac_mem",
          "common_nbr", "jaccard_nbr"]:
    x = df[f].values.astype(float)
    a, p, n = auc(y_c, x)
    a_neg, _, _ = auc(y_c, -x)
    cs = np.corrcoef(x[np.isfinite(x)], np.log(csize)[np.isfinite(x)])[0,1]
    # residualize on log comp size
    fin = np.isfinite(x)
    Af = np.column_stack([np.ones(fin.sum()), np.log(csize)[fin]])
    b, *_ = lstsq(Af, x[fin], rcond=None)
    xr = np.full(len(x), np.nan); xr[fin] = x[fin] - Af@b
    a_r, _, _ = auc(y_c, xr)
    best = max(a if a==a else 0, a_neg if a_neg==a_neg else 0)
    best_r = max(a_r if a_r==a_r else 0, 1-(a_r if a_r==a_r else 1))
    print(f"  {f:12s} bestAUC={best:.3f}  corr_w_logsize={cs:+.3f}  "
          f"AUC_resid(best)={best_r:.3f}")
