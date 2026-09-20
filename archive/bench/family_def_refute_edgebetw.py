import pandas as pd, numpy as np
from sklearn.mixture import GaussianMixture
from sklearn.metrics import roc_auc_score, silhouette_score
from sklearn.linear_model import LogisticRegression
from scipy.stats import skew, kurtosis
import warnings; warnings.filterwarnings('ignore')

rng = np.random.default_rng(0)
df = pd.read_csv('bench/family_def_features.tsv', sep='\t')
gs = df.groupby('group').size()
df['comp_edges'] = df['group'].map(gs)
big = df[df['comp_edges']>=50].copy().reset_index(drop=True)

def bc(v):
    v=v[~np.isnan(v)]; n=len(v); s=skew(v); k=kurtosis(v)
    return (s**2+1)/(k + 3*(n-1)**2/((n-2)*(n-3)))

print('='*70)
print('CONFOUND TEST 1: is edge_betw bimodality just a DEGREE effect?')
print('='*70)
# Residualize edge_betw on degree (deg_min,deg_max,common_nbr) and recompute BC/silhouette
X = big[['deg_min','deg_max','common_nbr','comp_edges']].values.astype(float)
X = np.nan_to_num(X)
from sklearn.linear_model import LinearRegression
y = big['edge_betw'].values
res = y - LinearRegression().fit(X, y).predict(X)
print('BC raw edge_betw|big      =', round(bc(y),3))
print('BC edge_betw ~deg-resid   =', round(bc(res),3))
g=GaussianMixture(2,random_state=0).fit(res.reshape(-1,1))
print('silhouette deg-resid 2clu =', round(silhouette_score(res.reshape(-1,1), g.predict(res.reshape(-1,1))),3))

print()
print('='*70)
print('CONFOUND TEST 2: does edge_betw add AUC for overlap/antisense BEYOND degree?')
print('='*70)
for target in ['overlap','antisense']:
    t = big[target].values
    deg = np.nan_to_num(big[['deg_min','deg_max','common_nbr','jaccard_nbr','comp_edges']].values.astype(float))
    # AUC of edge_betw alone
    auc_eb = roc_auc_score(t, big['edge_betw'])
    # AUC of degree-only logistic model (grouped CV would be ideal; in-sample is generous to degree)
    lr_deg = LogisticRegression(max_iter=1000).fit(deg, t)
    auc_deg = roc_auc_score(t, lr_deg.predict_proba(deg)[:,1])
    # AUC of degree + edge_betw
    both = np.column_stack([deg, big['edge_betw'].values])
    lr_both = LogisticRegression(max_iter=1000).fit(both, t)
    auc_both = roc_auc_score(t, lr_both.predict_proba(both)[:,1])
    print(f'{target:10s}: AUC(edge_betw)={auc_eb:.3f}  AUC(degree-only)={auc_deg:.3f}  AUC(degree+eb)={auc_both:.3f}  delta={auc_both-auc_deg:+.3f}')

print()
print('='*70)
print('CONFOUND TEST 3: residualize edge_betw on degree, then AUC for overlap/antisense')
print('='*70)
for target in ['overlap','antisense']:
    t = big[target].values
    print(f'{target:10s}: AUC(edge_betw_raw)={roc_auc_score(t,y):.3f}  AUC(edge_betw|deg-resid)={roc_auc_score(t,res):.3f}')
