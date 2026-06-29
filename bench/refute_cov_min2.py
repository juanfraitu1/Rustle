#!/usr/bin/env python
"""
Follow-up: the killer test. cov_min residualised on len_ratio collapses to AUC~0.50.
Now run the GROUPED PERMUTATION NULL on the RESIDUAL (cov_min minus len_ratio, and
minus the seq+len block). If the residual no longer survives the grouped null, then
the 'survives perm null' claim is carried entirely by the cDNA-length confound, not
by reciprocal coverage per se.

Also: len_ratio itself -- is cov_min just length-ratio re-expressed? And does the
GROUPED-RESAMPLING with only 3 positive groups have any power at all (bootstrap CI
on AUC over groups).
"""
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score

rng = np.random.default_rng(11)
DF = pd.read_csv('/mnt/c/Users/jfris/Desktop/Rustle/bench/family_def_features.tsv', sep='\t')

pos_full  = (DF['compara_pos']==1).values
neg_mask  = ((DF['antisense']==1) | (DF['overlap']==1)).values
pos_clean = pos_full & ~neg_mask
neg_clean = neg_mask & ~pos_full

SEQLEN = ['seq_id','cov_max','len_ratio','len_min','len_max']

def make_residual(target, baseline_cols):
    X = DF[baseline_cols].values.astype(float)
    t = DF[target].values.astype(float)
    m = np.isfinite(t) & np.all(np.isfinite(X), axis=1)
    Xf = np.column_stack([np.ones(m.sum()), X[m]])
    coef, *_ = np.linalg.lstsq(Xf, t[m], rcond=None)
    pred = np.full(len(t), np.nan)
    pred[m] = Xf @ coef
    return t - pred

def grouped_perm_null(vals, pos, neg, n_perm=8000, tag=""):
    sel = (pos | neg) & np.isfinite(vals)
    y = pos[sel].astype(int)
    g = DF['group'].values[sel]
    x = vals[sel]
    a_obs = roc_auc_score(y, x)
    obs = abs(a_obs - 0.5)
    uniq = np.unique(g)
    posgrp = np.array([y[g==gg].sum() > 0 for gg in uniq])
    n_posgrp = int(posgrp.sum())
    ge = valid = 0
    for _ in range(n_perm):
        permpos = np.zeros(len(uniq), bool)
        permpos[rng.choice(len(uniq), n_posgrp, replace=False)] = True
        gmap = dict(zip(uniq, permpos))
        py = np.array([gmap[gg] for gg in g]).astype(int)
        if py.sum()==0 or py.sum()==len(py): continue
        valid += 1
        if abs(roc_auc_score(py, x)-0.5) >= obs - 1e-12: ge += 1
    p = (ge+1)/(valid+1)
    print(f"   [{tag}] AUC_raw={a_obs:6.3f} |AUC-.5|={obs:5.3f} n_posgrp={n_posgrp} grouped_perm_p={p:6.4f}")
    return p

print("="*90)
print("KILLER TEST -- grouped permutation null on cov_min RESIDUALS")
print("="*90)
print("\n-- raw cov_min (baseline reference) --")
grouped_perm_null(DF['cov_min'].values.astype(float), pos_clean, neg_clean, tag="clean raw cov_min")
grouped_perm_null(DF['cov_min'].values.astype(float), pos_full,  neg_clean, tag="full  raw cov_min")

for blk_name, blk in [("len_ratio", ['len_ratio']), ("seq+len block", SEQLEN)]:
    res = make_residual('cov_min', blk)
    print(f"\n-- cov_min residualised on {blk_name} --")
    grouped_perm_null(res, pos_clean, neg_clean, tag=f"clean resid|{blk_name}")
    grouped_perm_null(res, pos_full,  neg_clean, tag=f"full  resid|{blk_name}")

print("\n"+"="*90)
print("Is len_ratio ALONE as good a differentiator as cov_min? (if yes, cov_min adds nothing)")
print("="*90)
grouped_perm_null(DF['len_ratio'].values.astype(float), pos_clean, neg_clean, tag="clean len_ratio")
grouped_perm_null(DF['len_ratio'].values.astype(float), pos_full,  neg_clean, tag="full  len_ratio")

print("\n"+"="*90)
print("GROUPED BOOTSTRAP power ceiling: resample WHOLE component-groups with replacement,")
print("recompute cov_min AUC. With only 3 (clean) / 10 (full) positive groups, how wide is the CI?")
print("="*90)
def grouped_bootstrap(pos, neg, feat, n_boot=2000, tag=""):
    sel = pos | neg
    y_all = pos[sel].astype(int)
    g_all = DF['group'].values[sel]
    x_all = DF.loc[sel, feat].values.astype(float)
    x_all = np.where(np.isfinite(x_all), x_all, np.nanmedian(x_all[np.isfinite(x_all)]))
    uniq = np.unique(g_all)
    aucs=[]
    for _ in range(n_boot):
        samp = rng.choice(uniq, len(uniq), replace=True)
        idx = np.concatenate([np.where(g_all==gg)[0] for gg in samp])
        yy, xx = y_all[idx], x_all[idx]
        if yy.sum()==0 or yy.sum()==len(yy): continue
        a = roc_auc_score(yy, xx); aucs.append(max(a,1-a))
    aucs=np.array(aucs)
    print(f"   [{tag}] {feat}: AUC_dir median={np.median(aucs):.3f}  2.5%={np.percentile(aucs,2.5):.3f}  97.5%={np.percentile(aucs,97.5):.3f}  frac<=0.6={np.mean(aucs<=0.6):.3f}")
grouped_bootstrap(pos_clean, neg_clean, 'cov_min', tag="clean")
grouped_bootstrap(pos_full,  neg_clean, 'cov_min', tag="full")
grouped_bootstrap(pos_clean, neg_clean, 'len_ratio', tag="clean")

print("\n"+"="*90)
print("ORIENTATION reality check: the claim says antisense/overlap negs are 'at/near 0'.")
print("But is cov_min~0 SPECIFIC to negatives, or is it just that MOST edges have cov_min~0?")
print("="*90)
x=DF['cov_min'].values.astype(float)
print(f"   full pop:      frac cov_min<0.05 = {np.mean(x<0.05):.3f}  median={np.median(x):.3f}")
print(f"   negatives:     frac cov_min<0.05 = {np.mean(x[neg_clean]<0.05):.3f}  median={np.median(x[neg_clean]):.3f}")
print(f"   clean pos:     cov_min values = {sorted(np.round(x[pos_clean],3))}")
print(f"   full pos:      cov_min values = {sorted(np.round(x[pos_full],3))}")
# base rate: among ALL edges with cov_min>=0.4, what frac are compara-pos? (precision ceiling)
hi = x>=0.4
print(f"\n   edges with cov_min>=0.4: {hi.sum()}  of which compara_pos: {(pos_full&hi).sum()}  -> precision={ (pos_full&hi).sum()/max(hi.sum(),1):.5f}")
print(f"   (i.e. thresholding cov_min flags {hi.sum()} edges to recover {(pos_full&hi).sum()} known paralog pairs)")
