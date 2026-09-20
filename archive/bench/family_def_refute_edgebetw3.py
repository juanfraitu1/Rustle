import pandas as pd, numpy as np
from sklearn.metrics import roc_auc_score, silhouette_score
from sklearn.mixture import GaussianMixture
from scipy.stats import skew, kurtosis, spearmanr
import warnings; warnings.filterwarnings('ignore')

df = pd.read_csv('bench/family_def_features.tsv', sep='\t')
gs = df.groupby('group').size()
df['comp_edges'] = df['group'].map(gs)
big = df[df['comp_edges']>=50].copy().reset_index(drop=True)

print('='*70)
print('ORIENTATION CHECK: can the ONLY trustworthy positive orient the big-comp signal?')
print('='*70)
print('compara_pos in big subset:', int((big['compara_pos']==1).sum()), '/ 12 total')
print('panel_real in big        :', int((big['panel_real']==1).sum()), '/ 5 total')
print('=> within-big claim is oriented ONLY by overlap/antisense (no trustworthy real-family positive present)')

print()
print('='*70)
print('IS THE BIMODALITY ITSELF A DEGREE ARTIFACT? (bimodal in deg too?)')
print('='*70)
def bc(v):
    v=v[~np.isnan(v)]; n=len(v); s=skew(v); k=kurtosis(v)
    return (s**2+1)/(k + 3*(n-1)**2/((n-2)*(n-3)))
for c in ['edge_betw','deg_min','deg_max','common_nbr','comp_edges']:
    v=big[c].values.astype(float); v=v[~np.isnan(v)]
    print(f'  BC({c:11s}|big)={bc(v):.3f}')
print('  Spearman edge_betw~deg_min|big =', round(spearmanr(big["edge_betw"],big["deg_min"],nan_policy="omit")[0],3))

print()
print('='*70)
print('WHAT IS IN THE HIGH-BETWEENNESS TAIL? (degree there)')
print('='*70)
thr=np.percentile(big['edge_betw'],99)
tail=big[big['edge_betw']>=thr]
rest=big[big['edge_betw']<thr]
print(f'tail (p99, n={len(tail)}): median deg_min={tail.deg_min.median():.0f} deg_max={tail.deg_max.median():.0f} comp_edges={tail.comp_edges.median():.0f}')
print(f'rest          (n={len(rest)}): median deg_min={rest.deg_min.median():.0f} deg_max={rest.deg_max.median():.0f} comp_edges={rest.comp_edges.median():.0f}')
print(f'tail antisense rate={tail.antisense.mean():.3f}  overlap rate={tail.overlap.mean():.3f}')
print(f'rest antisense rate={rest.antisense.mean():.3f}  overlap rate={rest.overlap.mean():.3f}')

print()
print('='*70)
print('IS edge_betw|big a RE-EXPRESSION of frac_ref (co-threading)? circularity check')
print('='*70)
m = big[['edge_betw','frac_ref','frac_mem','co_junc','aln_cov_small','contig']].dropna()
for c in ['frac_ref','frac_mem','co_junc','aln_cov_small','contig']:
    print(f'  corr edge_betw ~ {c:13s} = {m[["edge_betw",c]].corr().iloc[0,1]:+.3f}')

print()
print('='*70)
print('DOES edge_betw GENERALIZE on the ONE label that matters, full population?')
print('  (compara_pos AUC raw vs degree-residualized, all 12212 edges)')
print('='*70)
from sklearn.linear_model import LinearRegression
full=df.copy()
y=full['edge_betw'].values
X=np.nan_to_num(full[['deg_min','deg_max','common_nbr','comp_edges']].values.astype(float))
resid=y-LinearRegression().fit(X,y).predict(X)
print('  AUC(compara_pos | edge_betw raw)      =', round(roc_auc_score(full['compara_pos'],y),3))
print('  AUC(compara_pos | edge_betw deg-resid)=', round(roc_auc_score(full['compara_pos'],resid),3))
print('  AUC(compara_pos | deg_max alone)      =', round(roc_auc_score(full['compara_pos'],-X[:,1]),3),'(sign-free:',round(max(roc_auc_score(full['compara_pos'],X[:,1]),roc_auc_score(full['compara_pos'],-X[:,1])),3),')')
