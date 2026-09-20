#!/usr/bin/env python
"""
ADVERSARIAL REFUTATION of the claim that cov_min is a NEW, non-co-threading,
permutation-null-surviving differentiator of real-family vs bridge edges.

Three challenges:
  (1) CIRCULAR? Does cov_min just re-express co-threading (frac_ref/aln_cov_small)?
  (2) CONFOUNDED? Is it seq_id / cDNA length (len_ratio) / degree in disguise?
  (3) Does the Compara orientation survive a GROUPED permutation null and
      component-grouped resampling, with only 5 clean (or 12 full) positives?

Default verdict = circular/confounded/does-not-survive unless numbers force otherwise.
"""
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score
from sklearn.linear_model import LogisticRegression
from scipy.stats import spearmanr

rng = np.random.default_rng(7)
DF = pd.read_csv('/mnt/c/Users/jfris/Desktop/Rustle/bench/family_def_features.tsv', sep='\t')

COTHREAD = ['frac_ref','frac_mem','aln_cov_small','contig','co_junc']
SEQLEN   = ['seq_id','cov_max','len_ratio','len_min','len_max']
DEG      = ['deg_min','deg_max','edge_betw']

pos_full  = (DF['compara_pos']==1).values
neg_mask  = ((DF['antisense']==1) | (DF['overlap']==1)).values
pos_clean = pos_full & ~neg_mask
neg_clean = neg_mask & ~pos_full

def auc_dir(y, x):
    m = np.isfinite(x)
    if m.sum() < len(x):
        x = np.where(m, x, np.nanmedian(x[m]))
    a = roc_auc_score(y, x)
    return a, max(a, 1-a)

# ----------------------------------------------------------------------
print("="*92)
print("PART A -- directed AUC of cov_min vs co-threading & seq/len, on labelled subsets")
print("="*92)
for tag, pos, neg in [("CLEAN (5 pos)", pos_clean, neg_clean),
                      ("FULL  (12 pos)", pos_full, neg_clean)]:
    sel = pos | neg
    y = pos[sel].astype(int)
    print(f"\n[{tag}]  n_pos={y.sum()} n_neg={(1-y).sum()}")
    for c in ['cov_min','frac_ref','aln_cov_small','frac_mem','co_junc',
              'seq_id','cov_max','len_ratio']:
        x = DF.loc[sel, c].values.astype(float)
        a, ad = auc_dir(y, x)
        print(f"   {c:14s} AUC_raw={a:6.3f}  AUC_dir={ad:6.3f}")

# ----------------------------------------------------------------------
print("\n"+"="*92)
print("PART B -- CONFOUND: does cov_min add AUC BEYOND co-threading / seq+len?")
print("Residualise cov_min on a baseline block (full-pop OLS), then AUC of the RESIDUAL.")
print("If residual AUC collapses to ~0.5, cov_min is the baseline block in disguise.")
print("="*92)

def residual_auc(target, baseline_cols, pos, neg, tag):
    # fit OLS of target ~ baseline on FULL population (finite rows), get residual,
    # then measure directed AUC of the residual on the labelled subset.
    X = DF[baseline_cols].values.astype(float)
    t = DF[target].values.astype(float)
    m = np.isfinite(t) & np.all(np.isfinite(X), axis=1)
    Xf = np.column_stack([np.ones(m.sum()), X[m]])
    coef, *_ = np.linalg.lstsq(Xf, t[m], rcond=None)
    pred = np.full(len(t), np.nan)
    pred[m] = Xf @ coef
    resid = t - pred
    r2 = 1 - np.nansum((t-pred)**2)/np.nansum((t-np.nanmean(t))**2)
    sel = (pos | neg) & np.isfinite(resid)
    y = pos[sel].astype(int)
    if y.sum()==0 or y.sum()==len(y):
        print(f"   [{tag}] residual after {baseline_cols}: no labelled rows"); return
    a, ad = auc_dir(y, resid[sel])
    # also raw cov_min AUC on same sel for comparison
    araw, adraw = auc_dir(y, DF.loc[sel, target].values.astype(float))
    print(f"   [{tag}] baseline={baseline_cols}")
    print(f"        R2(target~baseline, fullpop)={r2:5.3f}   raw {target} AUC_dir={adraw:.3f}"
          f"   RESIDUAL AUC_dir={ad:.3f}  (drop={adraw-ad:+.3f})")

for tag, pos, neg in [("CLEAN", pos_clean, neg_clean), ("FULL", pos_full, neg_clean)]:
    residual_auc('cov_min', COTHREAD,            pos, neg, tag+" minus co-thread")
    residual_auc('cov_min', ['len_ratio'],       pos, neg, tag+" minus len_ratio")
    residual_auc('cov_min', ['seq_id'],          pos, neg, tag+" minus seq_id")
    residual_auc('cov_min', ['cov_max'],         pos, neg, tag+" minus cov_max")
    residual_auc('cov_min', SEQLEN,              pos, neg, tag+" minus seq+len block")
    residual_auc('cov_min', COTHREAD+SEQLEN+DEG, pos, neg, tag+" minus ALL(ct+seq+deg)")
    print()

# ----------------------------------------------------------------------
print("="*92)
print("PART C -- GROUPED PERMUTATION NULL (honest): shuffle the POSITIVE label over")
print("connected-component groups, holding whole groups together. With 5/12 positives the")
print("question is whether observed |AUC-0.5| beats a null built by relabelling groups.")
print("="*92)

def grouped_perm_null(pos, neg, feat, n_perm=5000, tag=""):
    sel = pos | neg
    y = pos[sel].astype(int)
    g = DF['group'].values[sel]
    x = DF.loc[sel, feat].values.astype(float)
    mfin = np.isfinite(x)
    if mfin.sum() < len(x):
        x = np.where(mfin, x, np.nanmedian(x[mfin]))
    a_obs = roc_auc_score(y, x)
    obs = abs(a_obs - 0.5)
    # group-level: how many groups carry a positive? permute WHICH groups are positive,
    # preserving the number of positive groups and each group's internal label pattern.
    uniq = np.unique(g)
    # which groups contain >=1 positive
    posgrp = np.array([y[g==gg].sum() > 0 for gg in uniq])
    n_posgrp = int(posgrp.sum())
    ge = 0
    valid = 0
    for _ in range(n_perm):
        permpos = np.zeros(len(uniq), bool)
        permpos[rng.choice(len(uniq), n_posgrp, replace=False)] = True
        # assign: a row is "positive" if its group is chosen positive
        # (collapse to group-level label to keep grouped structure honest)
        gmap = dict(zip(uniq, permpos))
        py = np.array([gmap[gg] for gg in g]).astype(int)
        if py.sum()==0 or py.sum()==len(py):
            continue
        valid += 1
        ap = roc_auc_score(py, x)
        if abs(ap-0.5) >= obs - 1e-12:
            ge += 1
    p = (ge+1)/(valid+1)
    print(f"   [{tag}] {feat:14s} AUC_raw={a_obs:6.3f} |AUC-.5|={obs:5.3f}"
          f"  n_posgrp={n_posgrp}/{len(uniq)}  grouped_perm_p={p:6.4f}  (valid={valid})")
    return p

print("\n--- CLEAN positives ---")
for f in ['cov_min','frac_ref','aln_cov_small','co_junc','seq_id','len_ratio','edge_betw']:
    grouped_perm_null(pos_clean, neg_clean, f, n_perm=5000, tag="clean")
print("\n--- FULL positives ---")
for f in ['cov_min','frac_ref','aln_cov_small','co_junc','seq_id','len_ratio','edge_betw']:
    grouped_perm_null(pos_full, neg_clean, f, n_perm=5000, tag="full")

# ----------------------------------------------------------------------
print("\n"+"="*92)
print("PART D -- how many DISTINCT positive GROUPS? (grouped-resampling power ceiling)")
print("="*92)
for tag, pos in [("clean", pos_clean), ("full", pos_full)]:
    g = DF['group'].values[pos]
    print(f"   {tag}: {pos.sum()} positive edges span {len(np.unique(g))} distinct component-groups: {sorted(set(g))}")
ng = DF['group'].values[neg_clean]
print(f"   negatives span {len(np.unique(ng))} distinct groups")

# ----------------------------------------------------------------------
print("\n"+"="*92)
print("PART E -- LEAVE-ONE-GROUP-OUT: does cov_min AUC hold when each positive group is")
print("dropped? If one component drives it, removing it should crater the AUC.")
print("="*92)
def logo(pos, neg, feat, tag):
    sel = pos | neg
    y = pos[sel].astype(int); g = DF['group'].values[sel]
    x = DF.loc[sel, feat].values.astype(float)
    x = np.where(np.isfinite(x), x, np.nanmedian(x[np.isfinite(x)]))
    base = max(roc_auc_score(y,x), 1-roc_auc_score(y,x))
    posg = sorted(set(g[y==1]))
    drops=[]
    for gg in posg:
        keep = g != gg
        if y[keep].sum()==0 or y[keep].sum()==len(y[keep]):
            drops.append((gg, np.nan)); continue
        a = roc_auc_score(y[keep], x[keep])
        drops.append((gg, max(a,1-a)))
    print(f"   [{tag}] {feat} full AUC_dir={base:.3f}; leave-one-pos-group-out:")
    for gg,a in drops:
        print(f"        drop group {gg:>6}: AUC_dir={a:.3f}")
logo(pos_clean, neg_clean, 'cov_min', "clean")
logo(pos_full,  neg_clean, 'cov_min', "full")

print("\nDONE.")
