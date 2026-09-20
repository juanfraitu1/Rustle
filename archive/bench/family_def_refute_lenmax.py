#!/usr/bin/env python
"""Adversarial verification of the len_max 'natural differentiator' claim.

Claim under test:
  "len_max: single most intrinsically bimodal feature (BC=0.852, silhouette=0.973),
   low size-correlation (corr +0.099), survives residualization (BC_resid=0.848).
   NEW non-co-threading natural two-mode feature, but NOT a real-vs-artifact
   differentiator (AUC vs compara only 0.630, positives split across both modes) --
   the two modes are short-cDNA vs long-cDNA, not real-vs-overmerge."

We refute / confirm by checking:
  1. CIRCULAR? Is len_max derived from co-threading / community-detection output?
     (test correlation/partial-corr with frac_ref, frac_mem, aln_cov_small, contig, co_junc)
  2. CONFOUNDED? Is it just seq_id, gene length (it IS a length), or degree in disguise?
  3. Does compara orientation survive permutation null (grouped) + grouped resampling,
     given only 12 positives?
  4. Reproduce the reported numbers (BC, silhouette, AUC, residualized).
"""
import numpy as np
import pandas as pd
from sklearn.mixture import GaussianMixture
from sklearn.metrics import roc_auc_score, silhouette_score
from sklearn.linear_model import LinearRegression
from scipy import stats

rng = np.random.default_rng(0)
df = pd.read_csv("/mnt/c/Users/jfris/Desktop/Rustle/bench/family_def_features.tsv",
                 sep="\t", na_values=[""])

print(f"rows={len(df)}  compara_pos={int(df.compara_pos.fillna(0).sum())}")

# ---- helper: sample-corrected bimodality coefficient (SAS) ----
def bimod_coef(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    n = len(x)
    if n < 4:
        return np.nan
    g = stats.skew(x, bias=False)
    k = stats.kurtosis(x, fisher=True, bias=False)  # excess kurtosis
    num = g**2 + 1.0
    den = k + 3.0 * (n - 1)**2 / ((n - 2) * (n - 3))
    return num / den

def gmm2_silhouette(x):
    x = np.asarray(x, float).reshape(-1, 1)
    x = x[np.isfinite(x[:, 0])]
    # standardize for stable GMM
    xs = (x - x.mean()) / (x.std() + 1e-12)
    g1 = GaussianMixture(1, random_state=0).fit(xs)
    g2 = GaussianMixture(2, random_state=0, n_init=3).fit(xs)
    bic_improve = g1.bic(xs) - g2.bic(xs)
    lab = g2.predict(xs)
    sil = np.nan
    if len(np.unique(lab)) == 2:
        # subsample for silhouette speed
        idx = rng.choice(len(xs), size=min(4000, len(xs)), replace=False)
        sil = silhouette_score(xs[idx], lab[idx])
    return bic_improve, sil, lab, g2

feat = "len_max"
x = df[feat].astype(float).values

print("\n=== 1. REPRODUCE reported intrinsic-bimodality numbers ===")
bc = bimod_coef(x)
bic_imp, sil, lab2, gmm = gmm2_silhouette(x)
print(f"  len_max  BC={bc:.3f}  BIC_improve(1->2)={bic_imp:.0f}  silhouette(2clust)={sil:.3f}")

print("\n=== 2. CIRCULARITY: corr of len_max with co-threading features ===")
cothread = ["frac_ref", "frac_mem", "aln_cov_small", "contig", "co_junc"]
for c in cothread:
    m = np.isfinite(x) & np.isfinite(df[c].astype(float).values)
    r, _ = stats.spearmanr(x[m], df[c].astype(float).values[m])
    rp, _ = stats.pearsonr(x[m], df[c].astype(float).values[m])
    print(f"  len_max ~ {c:14s} spearman={r:+.3f}  pearson={rp:+.3f}")

print("\n=== 3. CONFOUND: is len_max just seq_id / a length / degree? ===")
for c in ["seq_id", "len_min", "len_ratio", "deg_min", "deg_max", "common_nbr", "edge_betw"]:
    if c not in df.columns:
        continue
    y = df[c].astype(float).values
    m = np.isfinite(x) & np.isfinite(y)
    r, _ = stats.spearmanr(x[m], y[m])
    print(f"  len_max ~ {c:12s} spearman={r:+.3f}")

# component size confound
csz = df.groupby("group")["group"].transform("size").astype(float).values
m = np.isfinite(x)
r_sz, _ = stats.spearmanr(x[m], csz[m])
print(f"  len_max ~ log(comp_size) spearman={r_sz:+.3f}")

print("\n=== 4. RESIDUALIZE len_max on log(comp_size); recompute BC ===")
logsz = np.log(csz)
mm = np.isfinite(x) & np.isfinite(logsz)
reg = LinearRegression().fit(logsz[mm].reshape(-1, 1), x[mm])
resid = x[mm] - reg.predict(logsz[mm].reshape(-1, 1))
print(f"  BC(resid)={bimod_coef(resid):.3f}   (raw BC={bc:.3f})  pearson(len_max,logsz)={stats.pearsonr(x[mm],logsz[mm])[0]:+.3f}")

print("\n=== 5. ORIENTATION vs Compara (does len_max DISCRIMINATE real?) ===")
pos = df.compara_pos.fillna(0).values.astype(int)
# negatives: antisense or overlap edges (orthogonal artifact signals)
neg = ((df.antisense.fillna(0).values.astype(int) == 1) |
       (df.overlap.fillna(0).values.astype(int) == 1)).astype(int)
print(f"  positives={pos.sum()}  orthogonal-negatives={neg.sum()}")

# AUC vs compara (all non-positive treated as background)
mfin = np.isfinite(x)
auc_bg = roc_auc_score(pos[mfin], x[mfin])
print(f"  AUC(len_max, compara_pos vs background) = {auc_bg:.3f}")
# AUC vs compara among labelled only (pos vs orthogonal neg)
lab_mask = (pos == 1) | (neg == 1)
lm = lab_mask & mfin
if (pos[lm].sum() > 0) and ((neg[lm]).sum() > 0):
    auc_lab = roc_auc_score(pos[lm], x[lm])
    print(f"  AUC(len_max, compara_pos vs orthogonal-neg) = {auc_lab:.3f}  (n={lm.sum()})")

print("\n  Where do the 12 positives sit in the 2-mode GMM split?")
# mode means
order = np.argsort([gmm.means_[0,0], gmm.means_[1,0]])
xfin = x[mfin]
labfin = lab2  # predicted labels on finite rows (same order as gmm2_silhouette used finite x)
# recompute labels on full finite to align with pos
xs_full = ((xfin - xfin.mean())/(xfin.std()+1e-12)).reshape(-1,1)
lab_full = gmm.predict(xs_full)
hi_mode = int(np.argmax([gmm.means_[0,0], gmm.means_[1,0]]))
pos_fin = pos[mfin]
n_pos_hi = int(((lab_full==hi_mode) & (pos_fin==1)).sum())
n_pos_lo = int(((lab_full!=hi_mode) & (pos_fin==1)).sum())
print(f"  positives in HIGH-len_max mode={n_pos_hi}  in LOW mode={n_pos_lo}  (split => not real-vs-artifact)")

print("\n=== 6. PERMUTATION NULL (grouped by component) for AUC orientation ===")
groups = df.group.values
def grouped_perm_auc(score, label, groups, n=2000):
    obs = roc_auc_score(label, score)
    obs = max(obs, 1-obs)  # orientation-free
    # shuffle labels at GROUP level: permute which groups carry positives
    g_unique = np.unique(groups)
    # label is sparse; permute by reassigning each positive's group label to random group rows
    npos = int(label.sum())
    null = np.empty(n)
    gidx = {g: np.where(groups==g)[0] for g in g_unique}
    glist = list(g_unique)
    for i in range(n):
        perm = np.zeros_like(label)
        # pick npos random rows but grouped: choose random groups, mark one row each
        chosen = rng.choice(len(glist), size=npos, replace=True)
        for gi in chosen:
            rows = gidx[glist[gi]]
            perm[rng.choice(rows)] = 1
        a = roc_auc_score(perm, score)
        null[i] = max(a, 1-a)
    p = (np.sum(null >= obs) + 1) / (n + 1)
    return obs, p, null.mean(), null.std()

obs, p, nm, ns = grouped_perm_auc(x[mfin], pos[mfin], groups[mfin], n=2000)
print(f"  observed |AUC|={obs:.3f}  grouped-perm null mean={nm:.3f}+/-{ns:.3f}  p={p:.4f}")

print("\n=== 7. Compare to co_junc (the KNOWN baseline lever) on same labels ===")
cj = df.co_junc.astype(float).values
mcj = np.isfinite(cj)
auc_cj = roc_auc_score(pos[mcj], cj[mcj]); auc_cj = max(auc_cj, 1-auc_cj)
obs_cj, p_cj, _, _ = grouped_perm_auc(cj[mcj], pos[mcj], groups[mcj], n=2000)
print(f"  co_junc  |AUC vs compara|={auc_cj:.3f}  grouped-perm p={p_cj:.4f}")
print(f"  len_max  |AUC vs compara|={obs:.3f}  grouped-perm p={p:.4f}")
