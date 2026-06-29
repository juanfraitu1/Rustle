#!/usr/bin/env python
"""
Drill into edge_betw orientation: it shows AUC 0.939 vs compara AND 0.769 vs neg.
Both can't be 'high=real' AND 'high=artifact'. Resolve the orientation with
DIRECT pos-vs-neg, distribution quantiles, and degree-confound dissection.
"""
import numpy as np, pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score
from sklearn.linear_model import LinearRegression

df = pd.read_csv('/mnt/c/Users/jfris/Desktop/Rustle/bench/family_def_features.tsv', sep='\t')
y_comp = (df.compara_pos == 1).values
y_anti = (df.antisense == 1).values
y_over = (df.overlap == 1).values
y_neg  = y_anti | y_over
y_preal = (df.panel_real == 1).values
y_pcount= (df.panel_counter == 1).values
eb = df.edge_betw.values.astype(float)

def q(x, label):
    print(f"  {label:16s} n={len(x):5d}  min={x.min():.4g} q25={np.percentile(x,25):.4g} "
          f"med={np.median(x):.4g} q75={np.percentile(x,75):.4g} max={x.max():.4g} mean={x.mean():.4g}")

print("=== edge_betw distribution by group ===")
q(eb, "ALL")
q(eb[y_comp], "compara(+)")
q(eb[y_preal], "panel_real(+)")
q(eb[y_pcount], "panel_counter")
q(eb[y_neg], "neg(anti|over)")
q(eb[y_anti], "antisense")
q(eb[y_over], "overlap")
q(eb[~(y_comp|y_neg)], "unlabeled rest")

# The AUC vs compara uses ALL 12200 others as the negative class.
# The AUC vs neg uses ALL 11733 others as negative. These mix differently.
# DIRECT pos-vs-neg:
print("\n=== DIRECT compara(+) vs artifact(-) (only labeled rows) ===")
m = y_comp | y_neg
yy = y_comp[m]
a = roc_auc_score(yy, eb[m])
print(f"edge_betw AUC compara-vs-neg = {a:.3f}  (>0.5 => compara has HIGHER edge_betw than artifacts)")
print(f"  compara median edge_betw={np.median(eb[y_comp]):.4g}  neg median={np.median(eb[y_neg]):.4g}")

# Mann-Whitney for each contrast
def mw(pos, neg, name):
    try:
        u, p = stats.mannwhitneyu(eb[pos], eb[neg], alternative='two-sided')
        print(f"  {name:30s} MW p={p:.2e}  pos_med={np.median(eb[pos]):.4g} neg_med={np.median(eb[neg]):.4g}")
    except Exception as e:
        print(name, e)
print("\n=== Mann-Whitney contrasts ===")
mw(y_comp, ~(y_comp), "compara vs all-others")
mw(y_comp, y_neg, "compara vs neg")
mw(y_neg, ~(y_neg|y_comp), "neg vs unlabeled-rest")
mw(y_comp, ~(y_comp|y_neg), "compara vs unlabeled-rest")

# Now: what does HIGH edge_betw actually capture? Correlate with degree + group size
deg_min = df.deg_min.values; deg_max = df.deg_max.values
gsize = df.group.map(df.group.value_counts()).values
print("\n=== edge_betw vs structural confounds (Spearman) ===")
for nm, v in [('deg_min',deg_min),('deg_max',deg_max),('group_size',gsize),
              ('common_nbr',df.common_nbr.values),('jaccard_nbr',df.jaccard_nbr.values),
              ('same_chrom',df.same_chrom.values)]:
    r,p = stats.spearmanr(eb, v)
    print(f"  edge_betw ~ {nm:12s} rho={r:+.3f} p={p:.1e}")

# Degree-residualized edge_betw: does residual still separate compara vs neg?
X = np.column_stack([deg_min, deg_max, np.log1p(gsize)])
eb_res = eb - LinearRegression().fit(X, eb).predict(X)
print("\n=== degree+group-size residualized edge_betw ===")
m = y_comp | y_neg
a = roc_auc_score(y_comp[m], eb_res[m])
print(f"  resid edge_betw AUC compara-vs-neg = {a:.3f}")
a2 = roc_auc_score(y_comp, eb_res)
print(f"  resid edge_betw AUC compara-vs-all  = {a2:.3f}")
a3 = roc_auc_score(y_neg, eb_res)
print(f"  resid edge_betw AUC neg-vs-all      = {a3:.3f}")

# Key question: do compara positives sit on LOW edge_betw (real, within-clique) or HIGH (bridge)?
# Prior says high edge_betw = bridge = OVER-MERGE. Check the sign explicitly.
print("\n=== SIGN CHECK: is high edge_betw bridge(artifact) or real? ===")
print(f"  frac of compara with edge_betw==0 (no bridge): {(eb[y_comp]==0).mean():.3f}")
print(f"  frac of neg    with edge_betw==0            : {(eb[y_neg]==0).mean():.3f}")
print(f"  frac of all    with edge_betw==0            : {(eb==0).mean():.3f}")
# top-decile edge_betw enrichment
hi = eb >= np.percentile(eb, 90)
print(f"  compara in top-decile edge_betw: {(hi&y_comp).sum()}/{y_comp.sum()}")
print(f"  neg     in top-decile edge_betw: {(hi&y_neg).sum()}/{y_neg.sum()}")
lo = eb <= np.percentile(eb, 10)
print(f"  compara in bottom-decile edge_betw: {(lo&y_comp).sum()}/{y_comp.sum()}")
print(f"  neg     in bottom-decile edge_betw: {(lo&y_neg).sum()}/{y_neg.sum()}")
