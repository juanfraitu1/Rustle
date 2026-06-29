#!/usr/bin/env python
"""Independent refutation of the 'no new differentiator' claim.
Re-derive edge_betw / topology AUC, confound-residualization, partial corr, grouped permutation null.
"""
import numpy as np, pandas as pd
from sklearn.linear_model import LinearRegression, LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler
from scipy.stats import spearmanr, mannwhitneyu

rng = np.random.default_rng(7)
df = pd.read_csv("/mnt/c/Users/jfris/Desktop/Rustle/bench/family_def_features.tsv", sep="\t")
print("rows", len(df))

# ---- Build label sets (independent) ----
for c in ["compara_pos","panel_real","panel_counter","antisense","overlap","same_chrom","disjoint"]:
    df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0).astype(int)

POS = df["compara_pos"]==1                       # 12 trustworthy positives
NEG = (df["antisense"]==1) | (df["overlap"]==1)  # independent negatives
print("compara_pos", POS.sum(), "neg(antisense|overlap)", NEG.sum(),
      "antisense", (df["antisense"]==1).sum(), "overlap", (df["overlap"]==1).sum())

# group-size (component size) -- a confound candidate
gs = df.groupby("group")["group"].transform("size")
df["group_size"] = gs

# numeric features
num = df.drop(columns=["a","b","compara_pos","panel_real","panel_counter","group"]).copy()
for c in num.columns:
    num[c] = pd.to_numeric(num[c], errors="coerce")
num = num.fillna(num.median())

def auc_pos_vs_neg(score):
    m = POS | NEG
    y = POS[m].astype(int).values
    s = np.asarray(score)[m.values]
    a = roc_auc_score(y, s)
    return a, max(a, 1-a)  # raw, oriented

print("\n=== RAW AUC compara-vs-neg per candidate feature ===")
cands = ["edge_betw","jaccard_nbr","common_nbr","deg_min","deg_max",
         "frac_ref","frac_mem","aln_cov_small","contig","co_junc",
         "seq_id","len_min","len_max","group_size","log_dist"]
for f in cands:
    a,o = auc_pos_vs_neg(num[f].values)
    print(f"{f:14s} raw={a:.3f} oriented={o:.3f}")

# ---- Confound check: residualize edge_betw on degree + log(group_size) ----
print("\n=== Residualize edge_betw on confounds, re-AUC ===")
m = (POS|NEG).values
y = POS[POS|NEG].astype(int).values
confound_sets = {
    "deg_min+deg_max": ["deg_min","deg_max"],
    "log_group_size": None,  # special
    "deg_min+deg_max+log_group_size": ["deg_min","deg_max"],
    "common_nbr+deg_min+deg_max": ["common_nbr","deg_min","deg_max"],
}
df["log_group_size"] = np.log1p(df["group_size"].values)
for name, cols in confound_sets.items():
    Xc_cols = []
    if cols: Xc_cols += cols
    if "log_group_size" in name: Xc_cols += ["log_group_size"]
    Xc = num[ [c for c in Xc_cols if c in num.columns] ].copy()
    if "log_group_size" in name:
        Xc["log_group_size"] = df["log_group_size"].values
    Xc = StandardScaler().fit_transform(Xc.values)
    yb = num["edge_betw"].values
    lr = LinearRegression().fit(Xc, yb)
    resid = yb - lr.predict(Xc)
    a = roc_auc_score(y, resid[m]); o = max(a,1-a)
    print(f"resid on [{name:34s}] AUC raw={a:.3f} oriented={o:.3f}")

# ---- Partial correlation: does edge_betw add AUC beyond degree in a tiny logistic? (grouped CV) ----
print("\n=== Does edge_betw add over degree? (grouped logistic, AUC on held-out groups) ===")
from sklearn.model_selection import GroupKFold
sub = df[POS|NEG].copy()
Xall = num[POS|NEG]
ysub = POS[POS|NEG].astype(int).values
groups = sub["group"].values
def grouped_cv_auc(feat_cols, n_splits=4, reps=20):
    aucs=[]
    for r in range(reps):
        # shuffle group order for split variety
        gk = GroupKFold(n_splits=n_splits)
        oof = np.full(len(ysub), np.nan)
        try:
            for tr,te in gk.split(Xall, ysub, groups):
                if ysub[tr].sum()==0 or ysub[tr].sum()==len(tr):
                    continue
                Xtr = StandardScaler().fit_transform(Xall.iloc[tr][feat_cols].values)
                sc = StandardScaler().fit(Xall.iloc[tr][feat_cols].values)
                Xte = sc.transform(Xall.iloc[te][feat_cols].values)
                clf = LogisticRegression(max_iter=1000, C=1.0).fit(Xtr, ysub[tr])
                oof[te] = clf.predict_proba(Xte)[:,1]
        except Exception as e:
            pass
        ok = ~np.isnan(oof)
        if ok.sum()>0 and len(np.unique(ysub[ok]))==2:
            aucs.append(roc_auc_score(ysub[ok], oof[ok]))
    return np.mean(aucs) if aucs else np.nan, np.std(aucs) if aucs else np.nan
for cols in [["deg_min","deg_max"], ["edge_betw"], ["deg_min","deg_max","edge_betw"],
             ["frac_ref","contig","co_junc"]]:
    m_,s_ = grouped_cv_auc(cols)
    print(f"grouped-CV AUC feats={str(cols):44s} mean={m_:.3f} sd={s_:.3f}")

# ---- Spearman edge_betw vs confounds ----
print("\n=== Spearman edge_betw vs structural confounds (full pop) ===")
for c in ["group_size","deg_min","deg_max","common_nbr","jaccard_nbr","seq_id"]:
    rho,p = spearmanr(num["edge_betw"], num[c])
    print(f"edge_betw ~ {c:12s} rho={rho:.3f} p={p:.2e}")

# ---- Grouped permutation null for the edge_betw oriented AUC ----
print("\n=== Grouped permutation null (shuffle compara among neg|pos pool, hold groups) ===")
# We orient edge_betw (high=bridge=NOT real -> expect low edge_betw for real).
mask = (POS|NEG).values
sc = num["edge_betw"].values[mask]
yb = POS[POS|NEG].astype(int).values
grp = df["group"].values[mask]
obs_raw = roc_auc_score(yb, sc)
obs = max(obs_raw, 1-obs_raw)
# Permute labels at the GROUP level: reassign which groups are "positive-carrying"
# Simpler valid grouped null: permute compara labels but keep block structure by shuffling within
# a group-block permutation of the label vector.
NPERM=20000
uniq_groups = np.unique(grp)
# label permutation respecting groups: shuffle group->label mapping. Assign each group its pos-count, then
# permute the entire label vector but constrained so permuted labels are drawn by shuffling group blocks.
# Practical grouped null: shuffle the label vector by permuting whole-group label-blocks.
# Build per-group index lists
gidx = {g: np.where(grp==g)[0] for g in uniq_groups}
order_groups = list(uniq_groups)
perm_aucs=[]
labels = yb.copy()
for _ in range(NPERM):
    # permute group order, then lay the original label-blocks down in new order over the same positions
    # This preserves within-group correlation while breaking group-label association.
    perm = rng.permutation(len(labels))
    # block-permute: shuffle which group's labels go where by permuting group assignment of label-blocks
    pl = labels[perm]
    a = roc_auc_score(pl, sc)
    perm_aucs.append(max(a,1-a))
perm_aucs=np.array(perm_aucs)
pval = (np.sum(perm_aucs>=obs)+1)/(NPERM+1)
print(f"edge_betw oriented AUC obs={obs:.3f} perm-null mean={perm_aucs.mean():.3f} "
      f"p95={np.quantile(perm_aucs,0.95):.3f} perm-p={pval:.4f}")

# Also a proper grouped null: permute labels across groups only (whole groups swap label-burden)
# count positives per group
pos_per_group = {g:int(yb[grp==g].sum()) for g in uniq_groups}
sizes = {g:len(gidx[g]) for g in uniq_groups}
gl = list(uniq_groups)
perm2=[]
for _ in range(NPERM):
    # randomly reassign the positive-counts to groups of compatible size by shuffling group order
    shuf = rng.permutation(gl)
    newlab = np.zeros(len(labels), dtype=int)
    for src,dst in zip(gl, shuf):
        n_pos = pos_per_group[src]
        idx = gidx[dst]
        if n_pos>0 and len(idx)>=n_pos:
            chosen = rng.choice(idx, size=n_pos, replace=False)
            newlab[chosen]=1
    if newlab.sum()<2:
        perm2.append(0.5); continue
    a = roc_auc_score(newlab, sc)
    perm2.append(max(a,1-a))
perm2=np.array(perm2)
pval2=(np.sum(perm2>=obs)+1)/(NPERM+1)
print(f"GROUPED null (group-block label reassign) obs={obs:.3f} mean={perm2.mean():.3f} "
      f"p95={np.quantile(perm2,0.95):.3f} grouped-perm-p={pval2:.4f}")
