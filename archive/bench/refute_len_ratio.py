#!/usr/bin/env python
"""
Adversarial refutation of the claim:
  "len_ratio (cDNA length ratio): weaker secondary -- significant under the full-12
   null (perm_p=0.0065) but not the conservative 5-clean null (perm_p=0.063).
   Use only as a tie-breaker, not standalone."

Checks:
  (1) CIRCULAR? Does len_ratio re-express co-threading (frac_ref) -> partial corr / added value.
  (2) CONFOUNDED? Is it just seq_id / length / degree / cov_min in disguise?
  (3) Does Compara orientation survive a PROPER grouped permutation null and
      component-grouped resampling, with only 12 positives?
"""
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from scipy.stats import spearmanr

rng = np.random.default_rng(2026)
DF = pd.read_csv('/mnt/c/Users/jfris/Desktop/Rustle/bench/family_def_features.tsv', sep='\t')
N = len(DF)

pos_full  = (DF['compara_pos']==1).values
neg_mask  = ((DF['antisense']==1) | (DF['overlap']==1)).values
pos_clean = pos_full & ~neg_mask
neg_clean = neg_mask & ~pos_full          # negatives excluding any contaminated positives

print(f"n={N}  compara_pos={pos_full.sum()}  clean_pos={pos_clean.sum()}  neg={neg_mask.sum()}")
print(f"clean-null: pos={pos_clean.sum()} vs neg={neg_clean.sum()}")
print(f"full-null : pos={pos_full.sum()}  vs neg(excl contam)={neg_clean.sum()}")

def directed_auc(y, x):
    a = roc_auc_score(y, x)
    return max(a, 1-a), a

# ----------------------------------------------------------------------
# PROPER grouped permutation null:
#   We test directed AUC of len_ratio for separating pos from neg in the
#   labelled subset. Null = permute labels at the GROUP level: each connected
#   component is the unit. We randomly relabel which groups carry the positives,
#   keeping the number of positive *groups* fixed, so no within-component leakage
#   and the family structure is the resampling unit.
# ----------------------------------------------------------------------
def grouped_perm_p(pos, neg, feat, n_perm=5000):
    sel = pos | neg
    y = pos[sel].astype(int)
    x = DF.loc[sel, feat].values.astype(float)
    g = DF.loc[sel, 'group'].values
    # impute non-finite with median of the labelled subset (keeps rows; AUC robust)
    if not np.all(np.isfinite(x)):
        med = np.nanmedian(x)
        x = np.where(np.isfinite(x), x, med)
    obs_dir, obs_raw = directed_auc(y, x)
    obs = abs(obs_raw - 0.5)

    # group -> is it a positive group? (a group counts as pos if it contains >=1 pos row)
    groups = np.unique(g)
    grp_pos = np.array([ (y[g==gg].sum() > 0) for gg in groups ])
    n_pos_grp = grp_pos.sum()

    ge = 0
    valid = 0
    for _ in range(n_perm):
        perm_pos_grp = np.zeros(len(groups), dtype=bool)
        perm_pos_grp[rng.choice(len(groups), n_pos_grp, replace=False)] = True
        pos_set = set(groups[perm_pos_grp])
        py = np.array([1 if gg in pos_set else 0 for gg in g])
        # rows in a positive group are labelled positive; only rows that were originally
        # labelled (pos|neg) are in this subset, so this is a clean group relabel.
        if py.sum()==0 or py.sum()==len(py):
            continue
        valid += 1
        a = roc_auc_score(py, x)
        if abs(a-0.5) >= obs - 1e-12:
            ge += 1
    p = (ge+1)/(valid+1)
    return obs_dir, obs_raw, p, valid, n_pos_grp, len(groups)

print("\n=== (3) PROPER grouped permutation null (group-level relabel, 5000 perms) ===")
for feat in ['len_ratio','cov_min','frac_ref','seq_id']:
    for tag, pos, neg in [('clean', pos_clean, neg_clean), ('full12', pos_full, neg_clean)]:
        d, raw, p, v, npg, ng = grouped_perm_p(pos, neg, feat, n_perm=5000)
        print(f"  {feat:10s} [{tag:6s}] AUC_dir={d:.3f} (raw={raw:.3f}) grouped_perm_p={p:.4f}  (pos_groups={npg}/{ng}, valid_perms={v})")

# ----------------------------------------------------------------------
# (1) CIRCULAR / (2) CONFOUNDED: does len_ratio add AUC beyond confounds?
#   Build directed-AUC for len_ratio alone vs partial (residualised) versions,
#   and a logistic added-value test: cov_min alone vs cov_min+len_ratio.
# ----------------------------------------------------------------------
def added_value(pos, neg, base_feats, add_feat, tag):
    sel = pos | neg
    y = pos[sel].astype(int)
    Xb = DF.loc[sel, base_feats].values.astype(float)
    xa = DF.loc[sel, add_feat].values.astype(float).reshape(-1,1)
    # impute any non-finite with column medians
    def imp(M):
        M = M.copy()
        for j in range(M.shape[1]):
            col = M[:,j]; m = np.nanmedian(col); col[~np.isfinite(col)] = m; M[:,j]=col
        return M
    Xb, xa = imp(Xb), imp(xa)
    sc = StandardScaler()
    Xb_s = sc.fit_transform(Xb)
    Xfull_s = sc.fit_transform(np.hstack([Xb, xa]))
    # in-sample AUC (tiny positives -> use as descriptive, paired with perm null elsewhere)
    lr = LogisticRegression(max_iter=2000, C=1.0)
    auc_base = roc_auc_score(y, lr.fit(Xb_s, y).predict_proba(Xb_s)[:,1]) if len(np.unique(y))>1 else np.nan
    auc_full = roc_auc_score(y, lr.fit(Xfull_s, y).predict_proba(Xfull_s)[:,1]) if len(np.unique(y))>1 else np.nan
    # residualise add_feat on base (full population) then directed AUC of the residual
    from numpy.linalg import lstsq
    Xb_all = DF[base_feats].values.astype(float); xa_all = DF[add_feat].values.astype(float)
    Xb_all = imp(Xb_all); xa_all = imp(xa_all.reshape(-1,1)).ravel()
    A = np.hstack([np.ones((N,1)), Xb_all])
    coef,_,_,_ = lstsq(A, xa_all, rcond=None)
    resid = xa_all - A@coef
    d_resid,_ = directed_auc(y, resid[sel])
    d_raw,_   = directed_auc(y, DF.loc[sel,add_feat].values.astype(float))
    print(f"  [{tag}] base={base_feats} -> AUC_base={auc_base:.3f}  +{add_feat} -> AUC_full={auc_full:.3f}  "
          f"(delta={auc_full-auc_base:+.3f}); {add_feat} raw_AUC_dir={d_raw:.3f} residualised_AUC_dir={d_resid:.3f}")

print("\n=== (1)+(2) Added value of len_ratio beyond confounds (in-sample logistic + residualised AUC) ===")
for tag, pos, neg in [('clean', pos_clean, neg_clean), ('full12', pos_full, neg_clean)]:
    added_value(pos, neg, ['cov_min'],          'len_ratio', f'{tag}:cov_min base')
    added_value(pos, neg, ['frac_ref'],         'len_ratio', f'{tag}:frac_ref base')
    added_value(pos, neg, ['seq_id'],           'len_ratio', f'{tag}:seq_id  base')
    added_value(pos, neg, ['cov_min','frac_ref','seq_id'], 'len_ratio', f'{tag}:cov+frac+seq base')
    print()

# ----------------------------------------------------------------------
# Reverse: does cov_min add value beyond len_ratio? (which dominates which)
# ----------------------------------------------------------------------
print("=== Reverse direction: cov_min beyond len_ratio ===")
for tag, pos, neg in [('clean', pos_clean, neg_clean), ('full12', pos_full, neg_clean)]:
    added_value(pos, neg, ['len_ratio'], 'cov_min', f'{tag}:len_ratio base -> add cov_min')

# ----------------------------------------------------------------------
# Stability: leave-one-positive-out on the (already tiny) clean pos set.
# ----------------------------------------------------------------------
print("\n=== Stability: leave-one-clean-positive-out directed AUC for len_ratio ===")
pos_idx = np.where(pos_clean)[0]
aucs=[]
for drop in pos_idx:
    keep_pos = pos_clean.copy(); keep_pos[drop]=False
    sel = keep_pos | neg_clean
    y = keep_pos[sel].astype(int)
    x = DF.loc[sel,'len_ratio'].values.astype(float)
    d,_ = directed_auc(y,x); aucs.append(d)
print(f"  LOO AUC_dir: min={min(aucs):.3f} median={np.median(aucs):.3f} max={max(aucs):.3f}  (n={len(aucs)})")

# Cross-label check: panel_real vs panel_counter
print("\n=== Cross-label: panel_real(+) vs panel_counter(-) directed AUC for len_ratio & cov_min ===")
pr = (DF['panel_real']==1).values; pc=(DF['panel_counter']==1).values
sel = pr|pc; y=pr[sel].astype(int)
for feat in ['len_ratio','cov_min','frac_ref']:
    x=DF.loc[sel,feat].values.astype(float); d,raw=directed_auc(y,x)
    print(f"  {feat:10s} AUC_dir={d:.3f} (raw={raw:.3f}) on panel_real={pr.sum()} vs panel_counter={pc.sum()}")
