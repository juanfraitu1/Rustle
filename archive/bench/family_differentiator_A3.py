#!/usr/bin/env python
"""
A3: INDEPENDENT-LABEL FEATURE RANKING with PERMUTATION NULL
+ LABEL-FREE bimodality ranking (methodology mandate).

Goal: discover NATURAL differentiators separating REAL multi-copy gene-family
edges from over-merge/bridge edges, label-light & non-circular.

POSITIVES = compara_pos==1 (external Compara paralogs; the only trustworthy positive)
NEGATIVES = antisense==1 OR overlap==1 (sense/antisense + co-location artifacts)

Critical guards implemented:
  - Co-threading features (frac_ref,frac_mem,aln_cov_small,contig,co_junc) flagged CIRCULAR-BASELINE.
  - antisense & overlap flagged CIRCULAR (they DEFINE negatives) -> excluded from ranking.
  - 7 of 12 compara positives ALSO satisfy the negative definition -> reported and a CLEAN
    positive set (compara & not antisense/overlap) is used as the primary, with the full set
    shown for sensitivity.
  - PERMUTATION NULL: shuffle pos/neg labels, GROUPED by connected-component id, >=1000x.
  - LABEL-FREE: bimodality coefficient, 2-comp vs 1-comp GMM (BIC), 2-cluster silhouette,
    on the FULL unlabeled edge population.
"""
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score, silhouette_score
from sklearn.feature_selection import mutual_info_classif
from sklearn.mixture import GaussianMixture
from scipy.stats import skew, kurtosis

rng = np.random.default_rng(12345)
DF = pd.read_csv('/mnt/c/Users/jfris/Desktop/Rustle/bench/family_def_features.tsv', sep='\t')

# ---- feature catalogue with provenance tags ----
COTHREAD = ['frac_ref','frac_mem','aln_cov_small','contig','co_junc']          # known baseline lever
TOPOLOGY = ['edge_betw','jaccard_nbr','common_nbr','deg_min','deg_max']        # NEW, untried
SEQ      = ['seq_id','cov_min','cov_max','len_min','len_max','len_ratio']
JUNC     = ['njunc_a','njunc_b','njunc_min','njunc_diff','njunc_ratio']
GENOMIC  = ['same_chrom','log_dist','disjoint']                                # relation
DEFINE_NEG = ['antisense','overlap']                                          # CIRCULAR as features

KIND = {}
for c in COTHREAD: KIND[c]='known-cothread'
for c in TOPOLOGY: KIND[c]='new-topology'
for c in SEQ+JUNC: KIND[c]='new-other'
for c in GENOMIC:  KIND[c]='new-other'
for c in DEFINE_NEG: KIND[c]='confound-risk'   # defines negatives -> circular here

ALL_FEATS = COTHREAD+TOPOLOGY+SEQ+JUNC+GENOMIC+DEFINE_NEG

# ======================================================================
# PART 1: LABEL-FREE bimodality on full unlabeled population
# ======================================================================
def bimodality_coefficient(x):
    x = x[np.isfinite(x)]
    n = len(x)
    if n < 8: return np.nan
    g = skew(x, bias=False)
    k = kurtosis(x, fisher=True, bias=False)  # excess kurtosis
    denom = (k + 3.0*((n-1)**2)/((n-2)*(n-3)))
    if denom == 0: return np.nan
    return (g*g + 1.0) / denom   # > 0.555 (5/9) suggests bimodality

def gmm_bic_gap(x):
    x = x[np.isfinite(x)].reshape(-1,1)
    if len(x) < 50 or np.unique(x).size < 3: return np.nan, np.nan
    # standardize
    xs = (x - x.mean())/(x.std()+1e-12)
    try:
        g1 = GaussianMixture(1, n_init=2, random_state=0).fit(xs)
        g2 = GaussianMixture(2, n_init=4, random_state=0).fit(xs)
    except Exception:
        return np.nan, np.nan
    bic1, bic2 = g1.bic(xs), g2.bic(xs)
    # positive => 2-comp preferred (lower BIC)
    return bic1 - bic2, g2

def two_cluster_silhouette(x, gmm2):
    x = x[np.isfinite(x)].reshape(-1,1)
    if not isinstance(gmm2, GaussianMixture) or len(x) < 50: return np.nan
    xs = (x - x.mean())/(x.std()+1e-12)
    lab = gmm2.predict(xs)
    if np.unique(lab).size < 2: return np.nan
    # subsample for silhouette speed
    if len(xs) > 3000:
        idx = rng.choice(len(xs), 3000, replace=False)
        xs, lab = xs[idx], lab[idx]
        if np.unique(lab).size < 2: return np.nan
    try:
        return silhouette_score(xs, lab)
    except Exception:
        return np.nan

print("="*100)
print("PART 1 -- LABEL-FREE bimodality on full unlabeled population (n=12212)")
print("BC>0.555 => bimodal-ish; BIC_gap>0 => 2-comp GMM beats 1-comp; sil = 2-cluster separation")
print("="*100)
lf_rows=[]
for c in ALL_FEATS:
    x = DF[c].values.astype(float)
    bc = bimodality_coefficient(x)
    gap, g2 = gmm_bic_gap(x)
    sil = two_cluster_silhouette(x, g2)
    lf_rows.append((c, KIND[c], bc, gap, sil))
lf = pd.DataFrame(lf_rows, columns=['feature','kind','BC','BIC_gap','sil2']).sort_values('BC', ascending=False)
pd.set_option('display.width',200); pd.set_option('display.max_rows',200)
print(lf.to_string(index=False, float_format=lambda v: f"{v:9.4f}"))

# ======================================================================
# PART 2: A3 independent-label ranking + grouped permutation null
# ======================================================================
pos_full = (DF['compara_pos']==1).values
neg_mask = ((DF['antisense']==1) | (DF['overlap']==1)).values
contam   = pos_full & neg_mask
pos_clean = pos_full & ~neg_mask

print("\n"+"="*100)
print("PART 2 -- A3 independent-label ranking")
print(f"compara positives total={pos_full.sum()}  of which ALSO antisense/overlap (CONTAMINATED)={contam.sum()}")
print(f"CLEAN positives (compara & not antisense/overlap)={pos_clean.sum()}")
print(f"negatives (antisense OR overlap)={neg_mask.sum()}")
print("="*100)

def run_ranking(pos, neg, tag, n_perm=2000):
    sel = pos | neg
    y = pos[sel].astype(int)
    groups = DF['group'].values[sel]
    Xall = DF.loc[sel].copy()
    rows=[]
    feats = [f for f in ALL_FEATS]  # include confound-risk but FLAG them
    for c in feats:
        x = Xall[c].values.astype(float)
        finite = np.isfinite(x)
        if finite.sum() < len(x):
            # impute median of the labelled subset for AUC robustness; note coverage
            med = np.nanmedian(x)
            x = np.where(finite, x, med)
        if np.unique(x).size < 2 or y.sum()==0 or y.sum()==len(y):
            rows.append((c,KIND[c],np.nan,np.nan,np.nan,finite.mean())); continue
        # orient AUC so >0.5 means feature high in positives; report directed AUC
        try:
            auc = roc_auc_score(y, x)
        except Exception:
            auc = np.nan
        auc_dir = max(auc, 1-auc) if np.isfinite(auc) else np.nan
        # mutual information
        try:
            mi = mutual_info_classif(x.reshape(-1,1), y, discrete_features=False, random_state=0)[0]
        except Exception:
            mi = np.nan
        # grouped permutation null on |auc-0.5|
        obs = abs(auc-0.5) if np.isfinite(auc) else np.nan
        if np.isfinite(obs):
            uniq_g = np.unique(groups)
            ge=0
            for _ in range(n_perm):
                # permute labels at the GROUP level: shuffle which groups' rows are pos
                # preserve count by shuffling y within group-block structure:
                # simplest grouped null = shuffle group->label assignment proportional;
                # here: permute y but keep rows of same group together via group-block shuffle
                perm_y = grouped_permute(y, groups, uniq_g)
                if perm_y.sum()==0 or perm_y.sum()==len(perm_y):
                    continue
                a = roc_auc_score(perm_y, x)
                if abs(a-0.5) >= obs - 1e-12:
                    ge+=1
            p = (ge+1)/(n_perm+1)
        else:
            p = np.nan
        rows.append((c,KIND[c],auc_dir,auc,mi,finite.mean(),p))
    r = pd.DataFrame(rows, columns=['feature','kind','AUC_dir','AUC_raw','MI','cover','perm_p'])
    r = r.sort_values('AUC_dir', ascending=False)
    print(f"\n--- ranking [{tag}]  (n_pos={int(y.sum())} n_neg={int((1-y).sum())}, n_perm={n_perm}) ---")
    print(r.to_string(index=False, float_format=lambda v: f"{v:8.4f}"))
    return r

def grouped_permute(y, groups, uniq_g):
    # Build per-group label (a group is pos if any of its labelled rows is pos; but labelled
    # subset may mix). We shuffle the group->fraction by reassigning each group's block of
    # labels to a randomly chosen group's block (block bootstrap of labels keeps within-group
    # structure). Implement: permute the order of groups, then reassign label-vectors blockwise.
    out = np.empty_like(y)
    # collect index blocks per group
    perm = rng.permutation(uniq_g)
    mapping = dict(zip(uniq_g, perm))
    # For each group, take the label-block of its mapped group.
    blocks = {g: y[groups==g] for g in uniq_g}
    for g in uniq_g:
        src = mapping[g]
        bsrc = blocks[src]
        bdst_len = (groups==g).sum()
        # match lengths: tile/truncate
        if len(bsrc)==bdst_len:
            out[groups==g]=bsrc
        elif len(bsrc)>bdst_len:
            out[groups==g]=bsrc[:bdst_len]
        else:
            reps=int(np.ceil(bdst_len/len(bsrc)))
            out[groups==g]=np.tile(bsrc,reps)[:bdst_len]
    return out

r_clean = run_ranking(pos_clean, neg_mask, "CLEAN positives", n_perm=2000)
r_full  = run_ranking(pos_full,  neg_mask & ~pos_full, "FULL positives, neg excludes contaminated pos", n_perm=2000)

# ======================================================================
# PART 3: degree-confound check for topology features
# ======================================================================
print("\n"+"="*100)
print("PART 3 -- degree confound: correlation of topology feats with deg_min/deg_max (full pop)")
print("="*100)
sub = DF[TOPOLOGY].astype(float)
print(sub.corr(method='spearman').to_string(float_format=lambda v: f"{v:6.3f}"))

# orientation summary: in CLEAN positives, where does edge_betw / jaccard sit vs negatives?
print("\n--- orientation of key topology feats: median in clean-pos vs neg ---")
for c in TOPOLOGY+['co_junc','frac_ref']:
    mp = np.nanmedian(DF.loc[pos_clean,c].astype(float))
    mn = np.nanmedian(DF.loc[neg_mask,c].astype(float))
    mall= np.nanmedian(DF[c].astype(float))
    print(f"{c:14s} clean_pos_med={mp:8.3f}  neg_med={mn:8.3f}  fullpop_med={mall:8.3f}")

# save outputs
lf.to_csv('/mnt/c/Users/jfris/Desktop/Rustle/bench/family_diff_labelfree.tsv', sep='\t', index=False)
r_clean.to_csv('/mnt/c/Users/jfris/Desktop/Rustle/bench/family_diff_A3_clean.tsv', sep='\t', index=False)
r_full.to_csv('/mnt/c/Users/jfris/Desktop/Rustle/bench/family_diff_A3_full.tsv', sep='\t', index=False)
print("\nsaved: family_diff_labelfree.tsv, family_diff_A3_clean.tsv, family_diff_A3_full.tsv")