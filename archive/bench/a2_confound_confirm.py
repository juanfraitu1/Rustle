#!/usr/bin/env python
"""
Confirm: edge_betw 'signal' is the small-component confound. Compara pairs live
in tiny components where edge_betw is ~1 by construction. Degree-matched test:
within group-size strata, does edge_betw still separate compara from artifacts?
Also: are jaccard_nbr / common_nbr / deg any better, or all the same confound?
And: does ANYTHING non-cothread survive degree-matching?
"""
import numpy as np, pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score
from sklearn.linear_model import LinearRegression

df = pd.read_csv('/mnt/c/Users/jfris/Desktop/Rustle/bench/family_def_features.tsv', sep='\t')
y_comp = (df.compara_pos == 1).values
y_neg  = ((df.antisense==1)|(df.overlap==1)).values
gsize = df.group.map(df.group.value_counts()).values
df['gsize'] = gsize

print("=== component size of labeled edges ===")
print("compara group sizes:", sorted(gsize[y_comp]))
print("neg group sizes (quantiles):", np.percentile(gsize[y_neg],[0,25,50,75,100]).astype(int))
print("all  group sizes (quantiles):", np.percentile(gsize,[0,25,50,75,100]).astype(int))
print(f"compara in 2-node components: {(gsize[y_comp]==1).sum()}/{y_comp.sum()} "
      "(gsize=1 means edge_count=1 i.e. a 2-node component)")
# gsize here = number of EDGES sharing the group id; small => small component

# Degree-matched / size-matched: restrict to edges in SMALL components (where compara lives),
# then ask if edge_betw still discriminates compara vs neg.
TOPO = ['edge_betw','jaccard_nbr','common_nbr','deg_min','deg_max']
COTH = ['frac_ref','frac_mem','aln_cov_small','contig','co_junc']

for cap in [3, 5, 10, 50]:
    sel = (gsize <= cap)
    m = sel & (y_comp | y_neg)
    if (y_comp&sel).sum()<2 or (y_neg&sel).sum()<2:
        print(f"\n[component edge-count <= {cap}] too few labeled"); continue
    yy = y_comp[m]
    print(f"\n[component edge-count <= {cap}]  n_compara={yy.sum()} n_neg={(~yy).sum()}")
    for f in TOPO:
        try:
            a = roc_auc_score(yy, df[f].values[m])
            print(f"   {f:14s} AUC(compara vs neg | size-matched)={a:.3f}")
        except Exception as e:
            print(f"   {f}: {e}")

# Full degree+size residualization for EVERY non-cothread feature: which retain
# any compara-vs-neg signal AFTER removing the size confound?
print("\n=== residualize-on(deg_min,deg_max,log gsize) then compara-vs-neg AUC ===")
Xc = np.column_stack([df.deg_min, df.deg_max, np.log1p(gsize)])
m = y_comp | y_neg
yy = y_comp[m]
noncoth = [c for c in df.columns if c not in
           ['a','b','compara_pos','panel_real','panel_counter','group','gsize']+COTH]
rowz=[]
for f in noncoth:
    x = df[f].values.astype(float)
    x = np.where(np.isfinite(x), x, np.nanmedian(x[np.isfinite(x)]))
    raw = roc_auc_score(yy, x[m])
    res = x - LinearRegression().fit(Xc, x).predict(Xc)
    rA = roc_auc_score(yy, res[m])
    rowz.append((f, raw, rA, abs(raw-0.5), abs(rA-0.5)))
R = pd.DataFrame(rowz, columns=['feat','raw_AUC','resid_AUC','|raw-.5|','|resid-.5|']).sort_values('|resid-.5|',ascending=False)
print(R.to_string(index=False, float_format=lambda v:f"{v:.3f}"))

# And the co-threading baseline for reference (compara-vs-neg, raw):
print("\n=== co-threading baseline (compara vs neg, raw AUC) ===")
for f in COTH:
    x=df[f].values.astype(float); x=np.where(np.isfinite(x),x,np.nanmedian(x[np.isfinite(x)]))
    a=roc_auc_score(yy, x[m]); print(f"   {f:14s} AUC={a:.3f}")
