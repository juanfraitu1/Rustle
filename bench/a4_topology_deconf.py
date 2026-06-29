#!/usr/bin/env python
"""A4 part 2 -- de-confound the edge_betw 'win', and test clustering-coeff / k-core
as STANDALONE no-alignment de-bridge signals. Also a grouped permutation null."""
import numpy as np, pandas as pd, networkx as nx
from scipy.stats import spearmanr, mannwhitneyu
from sklearn.metrics import roc_auc_score

RNG = np.random.default_rng(0)
F='/mnt/c/Users/jfris/Desktop/Rustle/bench/family_def_features.tsv'
df=pd.read_csv(F,sep='\t')
gs=df.groupby('group').size(); df['gsize']=df['group'].map(gs)

print("="*78)
print("(DE-CONFOUND) edge_betw 'compara_pos AUC=0.94' -- is it just small-component bias?")
print("="*78)
# compara_pos lives in tiny comps; so does anything in a tiny comp score high edge_betw.
# Build a size-matched negative: edges in components of the SAME tiny sizes that are NOT compara.
cp = df[df.compara_pos==1]
print("compara_pos comp sizes:", sorted(cp.gsize.tolist()))
print("median edge_betw  compara_pos = %.3f | all = %.3f | gsize<=11 non-compara = %.3f"
      % (cp.edge_betw.median(), df.edge_betw.median(),
         df[(df.gsize<=11)&(df.compara_pos==0)].edge_betw.median()))
small = df[(df.gsize<=11)]
y = small.compara_pos.values
for f in ['edge_betw','jaccard_nbr','deg_min','common_nbr']:
    s=small[f].values;
    auc = roc_auc_score(y,s) if y.sum() else np.nan
    print(f"  size-matched (gsize<=11) AUC compara_pos by {f:11s} = {auc:.3f}")
print("  -> if edge_betw AUC collapses to ~0.5 here, the 0.94 was pure component-size confound.")

# also: does gsize ITSELF separate compara_pos vs all? (the actual signal)
auc_size = roc_auc_score(df.compara_pos.values, -df.gsize.values)
print(f"  AUC compara_pos by (small) -gsize over ALL edges = {auc_size:.3f}  (the real driver)")

print("\n" + "="*78)
print("(STANDALONE no-alignment) clustering-coeff & k-core as de-bridge signals")
print("="*78)
# Build the graph per component from the edge list (a,b). Compute edge-level topo natively
# to get clustering coeff of endpoints and k-core membership -- features minimap2 doesn't need.
G = nx.Graph()
for _,r in df.iterrows():
    G.add_edge(r['a'], r['b'])
print("graph: %d nodes, %d edges, %d comps" % (G.number_of_nodes(), G.number_of_edges(),
        nx.number_connected_components(G)))
clust = nx.clustering(G)                       # node clustering coeff
core  = nx.core_number(G)                       # node k-core
# edge features
def edge_feat(a,b):
    return dict(
        clust_min=min(clust[a],clust[b]), clust_max=max(clust[a],clust[b]),
        core_min=min(core[a],core[b]), core_max=max(core[a],core[b]),
    )
feats = df.apply(lambda r: edge_feat(r['a'],r['b']), axis=1, result_type='expand')
df = pd.concat([df, feats], axis=1)

# Jaccard of endpoints' triangle support is the classic 'within-community' edge embeddedness.
# Low clustering / low embeddedness = BRIDGE. Test vs co-thread bridge (frac_ref) in big comps.
co = df.dropna(subset=['frac_ref'])
big = co[co.gsize>=100]
print("\n  Within gsize>=100, spearman(feature, frac_ref)  [frac_ref low = co-thread bridge]:")
for f in ['clust_min','clust_max','core_min','core_max','edge_betw','jaccard_nbr','common_nbr']:
    r,_=spearmanr(big[f], big.frac_ref)
    tag = " <-- AGREE" if abs(r)>0.2 else ""
    print(f"    {f:11s} = {r:+.3f}{tag}")

# Does ANY topology feature separate the orthogonal NEGATIVES (antisense/overlap) from compara?
print("\n  AUC compara_pos vs (antisense|overlap) by native topo feature:")
neg = (df.antisense==1)|(df.overlap==1)
mask = (df.compara_pos==1)|neg
yy = df.compara_pos[mask].astype(int).values
for f in ['clust_min','clust_max','core_min','core_max','jaccard_nbr','edge_betw']:
    s=df[f][mask].values
    print(f"    {f:11s} = {roc_auc_score(yy,s):.3f}")

print("\n" + "="*78)
print("(GROUPED PERMUTATION NULL) edge_betw vs frac_ref agreement -- is +0.13 partial real?")
print("="*78)
# shuffle frac_ref WITHIN nothing? we test the AGREEMENT signal. Null: permute frac_ref across
# edges holding group structure -> destroy edge-level link, keep marginal. Group-block shuffle.
from sklearn.linear_model import LinearRegression
sub = big[['edge_betw','frac_ref','deg_min','deg_max','group']].dropna().copy()
R=sub[['edge_betw','frac_ref','deg_min','deg_max']].rank()
def resid(y, X):
    return y - LinearRegression().fit(X, y).predict(X)
Xc=R[['deg_min','deg_max']].values
ra=resid(R.edge_betw.values, Xc); rb=resid(R.frac_ref.values, Xc)
obs=spearmanr(ra,rb)[0]
# group-block permutation: permute frac_ref residual by whole group blocks
groups=sub.group.values
null=[]
for _ in range(2000):
    perm_g = RNG.permutation(np.unique(groups))
    gmap={g:pg for g,pg in zip(np.unique(groups), perm_g)}
    order=np.argsort([gmap[g] for g in groups], kind='stable')
    null.append(spearmanr(ra, rb[order])[0])
null=np.array(null)
p=(np.sum(np.abs(null)>=abs(obs))+1)/(len(null)+1)
print(f"  partial spearman(edge_betw, frac_ref | degree) obs = {obs:+.3f}")
print(f"  grouped-block permutation null: mean {null.mean():+.3f} sd {null.std():.3f}  p={p:.4f}")
print("  -> tiny effect; even if p<0.05 the magnitude (~0.13) is far below the degree confound (-0.78).")
