import pandas as pd, numpy as np
from sklearn.metrics import roc_auc_score
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GroupKFold
import warnings; warnings.filterwarnings('ignore')

rng = np.random.default_rng(0)
df = pd.read_csv('bench/family_def_features.tsv', sep='\t')
gs = df.groupby('group').size()
df['comp_edges'] = df['group'].map(gs)
big = df[df['comp_edges']>=50].copy().reset_index(drop=True)

print('='*70)
print('GROUPED-CV: edge_betw vs degree for overlap/antisense (whole groups held out)')
print('='*70)
degcols=['deg_min','deg_max','common_nbr','jaccard_nbr','comp_edges']
groups = big['group'].values
gkf = GroupKFold(n_splits=min(5, big['group'].nunique()))
for target in ['overlap','antisense']:
    t = big[target].values
    # pooled out-of-fold predictions
    def oof(Xcols):
        X = np.nan_to_num(big[Xcols].values.astype(float))
        pred = np.zeros(len(t))
        for tr,te in gkf.split(X,t,groups):
            if t[tr].sum()==0 or t[tr].sum()==len(tr): 
                pred[te]=t[tr].mean(); continue
            lr=LogisticRegression(max_iter=2000).fit(X[tr],t[tr])
            pred[te]=lr.predict_proba(X[te])[:,1]
        return pred
    auc_eb  = roc_auc_score(t, oof(['edge_betw']))
    auc_deg = roc_auc_score(t, oof(degcols))
    auc_both= roc_auc_score(t, oof(degcols+['edge_betw']))
    print(f'{target:10s}: grouped-CV AUC eb={auc_eb:.3f}  deg={auc_deg:.3f}  deg+eb={auc_both:.3f}  delta={auc_both-auc_deg:+.3f}')

print()
print('='*70)
print('PERMUTATION NULL: shuffle target labels grouped by component, recompute AUC(edge_betw)')
print('='*70)
def grouped_perm_auc(score, target_col, n=2000):
    t = big[target_col].values.copy()
    g = big['group'].values
    obs = roc_auc_score(t, score)
    obs = max(obs, 1-obs)  # two-sided magnitude
    # permute labels WITHIN preserving per-group counts: shuffle group->label assignment? 
    # Labels (overlap/antisense) are edge-level. Grouped shuffle = permute the block of labels by component.
    # Build per-component label vectors and shuffle which component each edge's label-block belongs to is ill-defined;
    # instead shuffle labels but keep them clustered: shuffle group ids -> reassign each group's label-set to another group's edges of same size is hard.
    # Standard grouped null: permute labels at GROUP granularity is impossible (labels are per-edge). 
    # Use: shuffle edge labels but restricted so permutation respects component blocks => permute whole rows' labels within random group order:
    null=[]
    uniq=np.unique(g)
    for _ in range(n):
        # assign each component a shuffled copy of another component's label distribution: simplest valid grouped null
        perm = t.copy()
        order = rng.permutation(len(uniq))
        # map: shuffle the label-block of each group to a permuted group of (closest) — approximate by shuffling group labels of whole edges
        # Simpler robust: shuffle the SCORE within, keep labels; but we want label null. 
        # Use block permutation: shuffle group ids then sort
        gmap = dict(zip(uniq, uniq[order]))
        newg = np.array([gmap[x] for x in g])
        idx = np.argsort(newg, kind='stable')
        permt = t[idx]
        a = roc_auc_score(permt, score)
        null.append(max(a,1-a))
    null=np.array(null)
    p=(np.sum(null>=obs)+1)/(n+1)
    return obs, p, null.mean()
for target in ['overlap','antisense']:
    obs,p,nm = grouped_perm_auc(big['edge_betw'].values, target)
    print(f'{target:10s}: obs|AUC|={obs:.3f}  grouped-perm-null-mean={nm:.3f}  p={p:.4f}')
