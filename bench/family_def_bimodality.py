#!/usr/bin/env python
"""
A1 -- LABEL-FREE BIMODALITY RANKING of edge features.

For every numeric feature, on the full 12212-edge population:
  (i)   bimodality coefficient BC = (skew^2 + 1) / kurtosis   (BC > 0.555 suggests bimodal)
  (ii)  fit 1- vs 2-component GaussianMixture; report BIC improvement of 2 over 1
  (iii) silhouette of the 2-cluster GMM assignment

Rank by intrinsic bimodality (NO labels). Then ORIENT each top feature with the
12 compara_pos positives (and panel/antisense/overlap negatives). Flag NEW (topology /
len_ratio / njunc) vs co-threading vs strand/genomic.

Missing ('') handled by per-feature drop.
"""
import sys
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.mixture import GaussianMixture
from sklearn.metrics import silhouette_score, roc_auc_score

RNG = 0
np.random.seed(RNG)

PATH = "/mnt/c/Users/jfris/Desktop/Rustle/bench/family_def_features.tsv"
df = pd.read_csv(PATH, sep="\t", dtype=str, keep_default_na=False)
print(f"Loaded {len(df)} edges, {df.shape[1]} columns", file=sys.stderr)

# --- column taxonomy ---
COTHREAD = ["frac_ref", "frac_mem", "aln_cov_small", "contig", "co_junc"]
TOPO     = ["edge_betw", "jaccard_nbr", "common_nbr", "deg_min", "deg_max"]
LEN_NJUNC = ["len_min", "len_max", "len_ratio",
             "njunc_a", "njunc_b", "njunc_min", "njunc_diff", "njunc_ratio"]
SEQ      = ["seq_id", "cov_min", "cov_max"]
STRAND_GEN = ["antisense", "same_chrom", "log_dist", "overlap", "disjoint"]

NUMERIC = SEQ + LEN_NJUNC + COTHREAD + STRAND_GEN + TOPO

def kind_of(f):
    if f in COTHREAD:   return "known-cothread"
    if f in TOPO:       return "new-topology"
    if f in LEN_NJUNC:  return "new-other"
    if f in SEQ:        return "seq"
    if f in STRAND_GEN: return "strand-genomic"
    return "other"

# coerce numeric, '' -> NaN
for c in NUMERIC + ["compara_pos", "panel_real", "panel_counter"]:
    df[c] = pd.to_numeric(df[c], errors="coerce")

# label vectors (for orientation only)
y_compara = (df["compara_pos"].fillna(0) > 0).values
y_panel_r = (df["panel_real"].fillna(0) > 0).values
y_panel_c = (df["panel_counter"].fillna(0) > 0).values
y_anti    = (df["antisense"].fillna(0) > 0).values
y_over    = (df["overlap"].fillna(0) > 0).values
print(f"compara_pos={y_compara.sum()} panel_real={y_panel_r.sum()} "
      f"panel_counter={y_panel_c.sum()} antisense={y_anti.sum()} overlap={y_over.sum()}",
      file=sys.stderr)


def bimodality_coeff(x):
    n = len(x)
    if n < 4:
        return np.nan
    # sample skew/kurtosis (Fisher kurtosis -> add 3 for Pearson)
    g1 = stats.skew(x, bias=False)
    g2 = stats.kurtosis(x, fisher=True, bias=False)  # excess kurtosis
    kurt_pearson = g2 + 3.0
    # SAS bimodality coefficient with small-sample correction
    num = g1**2 + 1.0
    denom = kurt_pearson + (3.0 * (n - 1) ** 2) / ((n - 2) * (n - 3))
    if denom <= 0:
        return np.nan
    return num / denom


def gmm_split(x):
    """Return (bic_improve, silhouette, hard_labels, gmm2) for 2- vs 1-comp."""
    X = x.reshape(-1, 1)
    g1 = GaussianMixture(n_components=1, covariance_type="full",
                         random_state=RNG, n_init=2).fit(X)
    g2 = GaussianMixture(n_components=2, covariance_type="full",
                         random_state=RNG, n_init=5).fit(X)
    bic1, bic2 = g1.bic(X), g2.bic(X)
    bic_improve = bic1 - bic2  # positive => 2-comp better
    lab = g2.predict(X)
    if len(np.unique(lab)) < 2:
        sil = np.nan
    else:
        # subsample for silhouette if huge
        if len(X) > 6000:
            idx = np.random.choice(len(X), 6000, replace=False)
            try:
                sil = silhouette_score(X[idx], lab[idx])
            except Exception:
                sil = np.nan
        else:
            try:
                sil = silhouette_score(X, lab)
            except Exception:
                sil = np.nan
    return bic_improve, sil, lab, g2


def orient(x_full, mask_full, lab_full, gmm2):
    """Which GMM mode do the positives prefer? Return frac of positives in
    the 'higher-mean' mode and in the mode the positives concentrate in."""
    pos = mask_full & ~np.isnan(x_full)
    npos = pos.sum()
    if npos == 0:
        return np.nan, np.nan, 0
    means = gmm2.means_.ravel()
    hi_comp = int(np.argmax(means))
    labs_pos = lab_full[pos]
    frac_hi = np.mean(labs_pos == hi_comp)
    # dominant mode among positives
    vals, cnts = np.unique(labs_pos, return_counts=True)
    dom = vals[np.argmax(cnts)]
    frac_dom = np.max(cnts) / npos
    return frac_hi, frac_dom, npos


rows = []
feat_cache = {}
for f in NUMERIC:
    x = df[f].values.astype(float)
    finite = np.isfinite(x)
    xv = x[finite]
    n = len(xv)
    nuniq = len(np.unique(xv))
    if n < 50 or nuniq < 3:
        rows.append(dict(feature=f, kind=kind_of(f), n=n, nuniq=nuniq,
                         BC=np.nan, bic_improve=np.nan, silhouette=np.nan,
                         note="too-few-uniq" if nuniq < 3 else "too-few-n"))
        continue
    bc = bimodality_coeff(xv)
    bic_improve, sil, lab_sub, gmm2 = gmm_split(xv)
    # map sub labels back to full index space
    lab_full = np.full(len(x), -1)
    lab_full[finite] = lab_sub
    feat_cache[f] = (x, finite, lab_full, gmm2)
    rows.append(dict(feature=f, kind=kind_of(f), n=n, nuniq=nuniq,
                     BC=bc, bic_improve=bic_improve, silhouette=sil, note=""))

res = pd.DataFrame(rows)
# rank: intrinsically bimodal = high BC, high silhouette, positive BIC improve
res["rank_score"] = (
    res["BC"].fillna(0) * 1.0
    + res["silhouette"].fillna(0) * 1.0
    + (res["bic_improve"].fillna(0) > 0).astype(float) * 0.0  # used as gate, not score
)
res = res.sort_values(["BC", "silhouette"], ascending=False).reset_index(drop=True)

pd.set_option("display.width", 200)
pd.set_option("display.max_columns", 30)
print("\n=== A1 LABEL-FREE BIMODALITY RANKING (all numeric features) ===")
print(res[["feature", "kind", "n", "nuniq", "BC", "bic_improve", "silhouette", "note"]]
      .to_string(index=False, float_format=lambda v: f"{v:.4f}"))

# --- ORIENT top bimodal features with sparse labels ---
print("\n=== ORIENTATION of bimodal features (BC>0.45 OR silhouette>0.5) ===")
print("frac_hi = frac of positives in the higher-mean GMM mode")
print("frac_dom = frac of positives in their dominant mode")
print("AUC vs each independent label (orient direction; >0.5 => higher feature = positive)\n")

orient_rows = []
top = res[(res["BC"] > 0.45) | (res["silhouette"] > 0.5)]["feature"].tolist()
# also always include the topology + len/njunc set so we report them explicitly
for f in TOPO + LEN_NJUNC:
    if f not in top and f in feat_cache:
        top.append(f)

for f in top:
    if f not in feat_cache:
        continue
    x, finite, lab_full, gmm2 = feat_cache[f]
    def auc_for(mask):
        m = mask & finite
        if m.sum() == 0:
            return np.nan
        neg = (~mask) & finite
        if neg.sum() == 0:
            return np.nan
        yy = np.concatenate([np.ones(m.sum()), np.zeros(neg.sum())])
        xx = np.concatenate([x[m], x[neg]])
        try:
            return roc_auc_score(yy, xx)
        except Exception:
            return np.nan
    fhi_c, fdom_c, npc = orient(x, y_compara, lab_full, gmm2)
    auc_c  = auc_for(y_compara)
    auc_pr = auc_for(y_panel_r)
    auc_pc = auc_for(y_panel_c)
    auc_an = auc_for(y_anti)
    auc_ov = auc_for(y_over)
    bc = res.loc[res.feature == f, "BC"].values[0]
    sil = res.loc[res.feature == f, "silhouette"].values[0]
    orient_rows.append(dict(
        feature=f, kind=kind_of(f), BC=bc, sil=sil,
        compara_frac_hi=fhi_c, compara_frac_dom=fdom_c, n_compara=npc,
        auc_compara=auc_c, auc_panelR=auc_pr, auc_panelC=auc_pc,
        auc_antisense=auc_an, auc_overlap=auc_ov))

od = pd.DataFrame(orient_rows)
print(od.to_string(index=False, float_format=lambda v: f"{v:.3f}"))

od.to_csv("/mnt/c/Users/jfris/Desktop/Rustle/bench/family_def_bimodality_orient.tsv",
          sep="\t", index=False)
res.to_csv("/mnt/c/Users/jfris/Desktop/Rustle/bench/family_def_bimodality_rank.tsv",
           sep="\t", index=False)
print("\nWrote rank + orient TSVs to bench/", file=sys.stderr)
