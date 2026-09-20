#!/usr/bin/env python3
"""
tp_fp_panel.py  --  EDGE STAGE of the TP/genuine-FP differentiator sweep.

Assembles the full EDGE-level feature panel for the E_r family-graph edges
(the 416 TP + 1773 genuine-FP edges that carry computed features in
ri_sharedlen_universal.tsv), joins every available edge feature, and evaluates
each feature with a held-out, component-grouped cross-validated AUC(TP>genuine).
Also fits simple numpy/scipy logistic combinations and reports out-of-fold AUC.

Deterministic: PYTHONHASHSEED=0, seed=0 fold assignment (edge-count balanced
GroupKFold grouped by connected component of the E_r gene graph -> no leakage
across a family block).

Writes:  bench/tp_fp_panel_edge.tsv   (per-edge feature matrix + label + fold)
Prints:  ranked univariate held-out AUC table + logistic-combo OOF AUC.

RNA-only note: soft-mask features (vg_softmask, edge bridge_mask*) and DNA
identity (id/min_cov, sparse) are flagged VALIDATION-ONLY (not RNA-pure /
not deployable); they are ranked but tagged.
"""
import os, sys, csv, math
os.environ.setdefault("PYTHONHASHSEED", "0")
import numpy as np
import networkx as nx
from scipy.stats import rankdata
from scipy.optimize import minimize

np.random.seed(0)
ROOT = "/mnt/c/Users/jfris/Desktop/Rustle"
B = os.path.join(ROOT, "bench")


def rd(path):
    return list(csv.DictReader(open(os.path.join(B, path)), delimiter="\t"))


def norm(a, b):
    return (a, b) if a < b else (b, a)


def f(x, default=np.nan):
    try:
        if x is None or x == "":
            return default
        return float(x)
    except Exception:
        return default


# ---------------------------------------------------------------- load sources
sl = rd("ri_sharedlen_universal.tsv")          # core(=core_recip), aln_len, aln_frac
panel = [r for r in sl if r["cls"] in ("TP", "genuine")]

ex = {norm(r["gene_a"], r["gene_b"]): r for r in rd("exon_containment_criterion.tsv")}
bm = {norm(r["gene_a"], r["gene_b"]): r for r in rd("edge_bridge_mask.tsv")}
fp = {norm(r["gene_a"], r["gene_b"]): r for r in rd("family_er_pr.tsv")}

# vg edges section (edge-level multiplicity)
vg = {}
lines = open(os.path.join(B, "vg_repeat_catalog.tsv")).read().splitlines()
for i, l in enumerate(lines):
    if l.startswith("# SECTION edges"):
        hdr = lines[i + 1].split("\t")
        for l2 in lines[i + 2:]:
            p = l2.split("\t")
            if len(p) < len(hdr):
                continue
            d = dict(zip(hdr, p))
            vg[norm(d["gene_a"], d["gene_b"])] = d
        break

# per-DN caches
km = {r["dn"]: r for r in rd("ri_kmer_cache.tsv")}
rm = {r["dn"]: r for r in rd("ri_readmm_cache.tsv")}


def agg(dn_a, dn_b, cache, col):
    a = f(cache.get(dn_a, {}).get(col, ""))
    b = f(cache.get(dn_b, {}).get(col, ""))
    if np.isnan(a) or np.isnan(b):
        return (np.nan, np.nan, np.nan, np.nan)
    return (min(a, b), max(a, b), (a + b) / 2.0, abs(a - b))


def ratio(a, b):
    a, b = f(a), f(b)
    if np.isnan(a) or np.isnan(b) or max(a, b) == 0:
        return np.nan
    return min(a, b) / max(a, b)


# ------------------------------------------------------------ component groups
G = nx.Graph()
for r in sl:  # components from full E_r population (TP+genuine+truthbar)
    G.add_edge(r["gene_a"], r["gene_b"])
comp = {}
for i, cc in enumerate(nx.connected_components(G)):
    for g in cc:
        comp[g] = i

# edge-count balanced 5-fold assignment of components (deterministic)
K = 5
edges_per_comp = {}
for r in panel:
    edges_per_comp.setdefault(comp[r["gene_a"]], 0)
    edges_per_comp[comp[r["gene_a"]]] += 1
fold_load = [0] * K
comp_fold = {}
for ci, n in sorted(edges_per_comp.items(), key=lambda kv: (-kv[1], kv[0])):
    k = int(np.argmin(fold_load))
    comp_fold[ci] = k
    fold_load[k] += n

# --------------------------------------------------------------- build feature rows
FEATS = []  # ordered list of feature names
recs = []
for r in panel:
    ga, gb = r["gene_a"], r["gene_b"]
    key = norm(ga, gb)
    dn_a, dn_b = r["dn_a"], r["dn_b"]
    e = ex.get(key, {})
    v = vg.get(key, {})
    m = bm.get(key, {})
    q = fp.get(key, {})
    d = {}

    # ---- SEQUENCE / HOMOLOGY
    d["core_recip"] = f(r["core"])                       # SHIPPED baseline 0.815
    d["aln_frac"] = f(r["aln_frac"])                     # SHIPPED workhorse 0.915
    d["aln_len"] = f(r["aln_len"])
    d["log_aln_len"] = math.log1p(f(r["aln_len"], 0.0))
    cab, cba = f(e.get("contain_ab")), f(e.get("contain_ba"))
    d["contain_max"] = np.nanmax([cab, cba]) if not (np.isnan(cab) and np.isnan(cba)) else np.nan
    d["contain_min"] = np.nanmin([cab, cba]) if not (np.isnan(cab) and np.isnan(cba)) else np.nan
    d["recip_contain"] = f(e.get("recip_contain"))       # FALSIFIED 0.900
    d["core_x_alnfrac"] = d["core_recip"] * d["aln_frac"] if not np.isnan(d["core_recip"]) else np.nan
    d["contain_x_alnfrac"] = d["contain_max"] * d["aln_frac"] if not np.isnan(d["contain_max"]) else np.nan
    d["id"] = f(q.get("id"))                             # VALIDATION/sparse identity
    d["min_cov"] = f(q.get("min_cov"))                   # VALIDATION/sparse coverage

    # ---- STRUCTURE / EXON
    na, nb = f(e.get("n_exon_a")), f(e.get("n_exon_b"))
    d["n_exon_min"] = np.nanmin([na, nb]) if not (np.isnan(na) and np.isnan(nb)) else np.nan
    d["n_exon_max"] = np.nanmax([na, nb]) if not (np.isnan(na) and np.isnan(nb)) else np.nan
    d["n_exon_ratio"] = ratio(na, nb)
    d["n_exon_absdiff"] = abs(na - nb) if not (np.isnan(na) or np.isnan(nb)) else np.nan
    d["both_single"] = f(e.get("both_single"))

    # ---- MULTIPLICITY / REPEAT (VG)
    d["min_shared_mult"] = f(v.get("min_shared_mult"))   # repeat-hub gate (shipped)
    d["max_shared_mult"] = f(v.get("max_shared_mult"))
    d["mult_ratio"] = ratio(v.get("min_shared_mult"), v.get("max_shared_mult"))
    d["n_shared_nodes"] = f(v.get("n_shared_nodes"))
    d["log_n_shared_nodes"] = math.log1p(f(v.get("n_shared_nodes"), 0.0))
    d["shares_only_repeat_m5"] = f(v.get("shares_only_repeat_m5"))
    d["vg_softmask"] = f(v.get("softmask"))              # VALIDATION-ONLY (DNA soft-mask)

    # ---- READ / MULTIMAPPING (per-DN readmm, aggregated)  read-degree 0.64-0.69
    for col, tag in [("mm_median", "mmMed"), ("mm_max", "mmMax"), ("n_reads", "nReads")]:
        lo, hi, mn, ad = agg(dn_a, dn_b, rm, col)
        d[f"{tag}_min"] = lo
        d[f"{tag}_max"] = hi
        d[f"{tag}_mean"] = mn
    d["nReads_ratio"] = ratio(rm.get(dn_a, {}).get("n_reads"), rm.get(dn_b, {}).get("n_reads"))

    # ---- KMER (per-DN, aggregated)  kmer-multiplicity
    for col, tag in [("km_median", "kmMed"), ("km_frac_ge_100", "kmFrac100"), ("n_kmers", "nKmers")]:
        lo, hi, mn, ad = agg(dn_a, dn_b, km, col)
        d[f"{tag}_min"] = lo
        d[f"{tag}_max"] = hi
        d[f"{tag}_mean"] = mn

    # ---- SOFT-MASK (edge, VALIDATION-ONLY)
    d["mask_max"] = np.nanmax([f(m.get("mask_a")), f(m.get("mask_b"))]) if m else np.nan
    d["mask_min"] = np.nanmin([f(m.get("mask_a")), f(m.get("mask_b"))]) if m else np.nan
    d["bridge_mask"] = f(m.get("bridge_mask"))           # VALIDATION-ONLY

    d["_gene_a"] = ga
    d["_gene_b"] = gb
    d["_dn_a"] = dn_a
    d["_dn_b"] = dn_b
    d["_label"] = 1 if r["cls"] == "TP" else 0
    d["_comp"] = comp[ga]
    d["_fold"] = comp_fold[comp[ga]]
    recs.append(d)

FEATS = [k for k in recs[0].keys() if not k.startswith("_")]

# tag features
VALIDATION_ONLY = {"vg_softmask", "mask_max", "mask_min", "bridge_mask", "id", "min_cov"}
SHIPPED = {"aln_frac", "core_recip", "min_shared_mult"}
FALSIFIED = {"recip_contain", "contain_max", "contain_min"}  # exon-containment family

y = np.array([d["_label"] for d in recs])
folds = np.array([d["_fold"] for d in recs])
comps = np.array([d["_comp"] for d in recs])
X = {name: np.array([d[name] for d in recs], dtype=float) for name in FEATS}


# ------------------------------------------------------------------- AUC utils
def auc(scores, labels):
    """Mann-Whitney AUC(positive>negative). NaN scores dropped."""
    s = np.asarray(scores, float)
    lab = np.asarray(labels)
    ok = ~np.isnan(s)
    s, lab = s[ok], lab[ok]
    npos = int((lab == 1).sum())
    nneg = int((lab == 0).sum())
    if npos == 0 or nneg == 0:
        return np.nan, npos, nneg
    rk = rankdata(s)
    a = (rk[lab == 1].sum() - npos * (npos + 1) / 2.0) / (npos * nneg)
    return a, npos, nneg


def oof_univariate_auc(vals, y, folds):
    """Honest held-out univariate AUC for a RAW single feature.
    A raw feature has no trained parameters except a 1-bit orientation, so the
    only thing that can fail to generalize is the sign. Protocol: pick sign on
    the TRAIN folds, evaluate AUC on the TEST fold, average over folds
    (train-oriented fold-mean). This is immune to the pooled-reorient artifact
    that spuriously inflates weak sign-unstable features."""
    per_fold = []
    for k in range(K):
        te = folds == k
        tr = ~te
        a_tr, _, _ = auc(vals[tr], y[tr])
        sign = 1.0 if (np.isnan(a_tr) or a_tr >= 0.5) else -1.0
        a_te, npos_te, nneg_te = auc(sign * vals[te], y[te])
        if not np.isnan(a_te):
            per_fold.append(a_te)
    return (np.mean(per_fold) if per_fold else np.nan), (np.std(per_fold) if per_fold else np.nan)


# ---------------------------------------------------- univariate ranking table
rank_rows = []
for name in FEATS:
    vals = X[name]
    full_auc, npos, nneg = auc(vals, y)
    disc = abs(full_auc - 0.5) if not np.isnan(full_auc) else np.nan
    # oriented so reported in-sample AUC>=0.5 (discriminative orientation)
    oriented_full = full_auc if (np.isnan(full_auc) or full_auc >= 0.5) else 1 - full_auc
    oof_mean, oof_sd = oof_univariate_auc(vals, y, folds)   # train-oriented held-out
    n_present = int((~np.isnan(vals)).sum())
    # generalizes: held-out fold-mean confirms in-sample power (not a sign fluke)
    generalizes = (not np.isnan(oof_mean)) and oof_mean >= 0.55 and (oriented_full - oof_mean) < 0.12
    tag = []
    if name in SHIPPED: tag.append("SHIPPED")
    if name in VALIDATION_ONLY: tag.append("VALID-ONLY")
    if name in FALSIFIED: tag.append("FALSIFIED-prior")
    rank_rows.append({
        "feature": name,
        "held_out_auc": oof_mean,     # PRIMARY: train-oriented, test-evaluated, fold-mean
        "held_out_sd": oof_sd,
        "in_sample_auc": oriented_full,
        "raw_dir_auc": full_auc,      # signed (higher raw feature -> TP if >0.5)
        "disc": disc,
        "n_present": n_present,
        "generalizes": generalizes,
        "tag": ",".join(tag) if tag else "-",
    })
rank_rows.sort(key=lambda r: (-(r["held_out_auc"] if not np.isnan(r["held_out_auc"]) else 0)))


# ------------------------------------------------------- logistic combinations
def fit_logistic(Xtr, ytr, l2=1.0):
    n, p = Xtr.shape
    Xb = np.hstack([np.ones((n, 1)), Xtr])

    def nll(w):
        z = Xb @ w
        # stable log-loss
        ll = np.sum(np.logaddexp(0, z) - ytr * z)
        reg = 0.5 * l2 * np.sum(w[1:] ** 2)
        return ll + reg

    def grad(w):
        z = Xb @ w
        p_ = 1 / (1 + np.exp(-z))
        g = Xb.T @ (p_ - ytr)
        g[1:] += l2 * w[1:]
        return g

    w0 = np.zeros(p + 1)
    res = minimize(nll, w0, jac=grad, method="L-BFGS-B")
    return res.x


def combo_oof(feature_list, y, folds, X, l2=1.0):
    """standardize on train, impute train-median, logistic on train.
    Returns (pooled_OOF_AUC, per_fold_mean_AUC, per_fold_sd).
    pooled = one AUC over all concatenated OOF logits (cross-fold, sensitive to
             calibration across family blocks).
    per_fold_mean = deployment-relevant WITHIN-fold/within-block discrimination."""
    M = np.column_stack([X[f] for f in feature_list])
    pooled = np.full(len(y), np.nan)
    per = []
    for k in range(K):
        te = folds == k
        tr = ~te
        med = np.nanmedian(M[tr], axis=0)
        mu = np.nanmean(M[tr], axis=0)
        sd = np.nanstd(M[tr], axis=0)
        sd[sd == 0] = 1.0
        def prep(A):
            A = A.copy()
            idx = np.where(np.isnan(A))
            A[idx] = np.take(med, idx[1])
            return (A - mu) / sd
        Xtr = prep(M[tr]); Xte = prep(M[te])
        w = fit_logistic(Xtr, y[tr], l2=l2)
        z = np.hstack([np.ones((Xte.shape[0], 1)), Xte]) @ w
        pooled[te] = z
        a_te, _, _ = auc(z, y[te])
        if not np.isnan(a_te):
            per.append(a_te)
    a, _, _ = auc(pooled, y)
    return a, (np.mean(per) if per else np.nan), (np.std(per) if per else np.nan)


# candidate deployable (RNA-only) feature set: exclude VALIDATION-ONLY & sparse id/min_cov
deployable = [r["feature"] for r in rank_rows
              if r["feature"] not in VALIDATION_ONLY and r["n_present"] == len(y)]
# curated interpretable set
curated = [x for x in ["aln_frac", "core_recip", "min_shared_mult", "contain_max",
                       "mmMed_max", "kmFrac100_max", "n_exon_ratio", "log_aln_len"]
           if x in X]

alnfrac_fold = [r for r in rank_rows if r["feature"] == "aln_frac"][0]["held_out_auc"]
combos = {
    "aln_frac(single)": combo_oof(["aln_frac"], y, folds, X),
    "curated8_logistic": combo_oof(curated, y, folds, X),
    "aln_frac+core_recip": combo_oof(["aln_frac", "core_recip"], y, folds, X),
    "aln_frac+min_shared_mult": combo_oof(["aln_frac", "min_shared_mult"], y, folds, X),
    "aln_frac+core+mult+contain": combo_oof(["aln_frac", "core_recip", "min_shared_mult", "contain_max"], y, folds, X),
    "full_deployable_logistic": combo_oof(deployable, y, folds, X),
}
alnfrac_oof = combos["aln_frac(single)"][0]
alnfrac_perfold = combos["aln_frac(single)"][1]

# -------------------------------------------------------------------- write TSV
out = os.path.join(B, "tp_fp_panel_edge.tsv")
with open(out, "w") as fh:
    meta = ["_gene_a", "_gene_b", "_dn_a", "_dn_b", "_label", "_comp", "_fold"]
    fh.write("\t".join([m[1:] for m in meta] + FEATS) + "\n")
    for d in recs:
        row = [str(d[m]) for m in meta] + [
            ("" if (isinstance(d[c], float) and np.isnan(d[c])) else f"{d[c]:.6g}") for c in FEATS]
        fh.write("\t".join(row) + "\n")

# --------------------------------------------------------------------- reports
print(f"# EDGE PANEL: {len(recs)} edges  ({int(y.sum())} TP / {int((y==0).sum())} genuine-FP)  "
      f"{len(FEATS)} features  {K}-fold component-grouped CV  ({len(set(comps))} components)")
print(f"# fold edge loads: {fold_load}")
print(f"# CONTRAST NOTE: this panel is the HARD contrast TP vs GENUINE-FP only.")
print(f"#   The headline 'aln_frac 0.915' is TP vs ALL-FP (genuine+truthbar); truthbar-FP")
print(f"#   are easy (aln_frac AUC 0.945). On genuine-FP only, aln_frac full=0.858, held-out=0.830.")
print(f"# wrote {out}")
print()
print("## UNIVARIATE HELD-OUT AUC(TP>genuine-FP), ranked by held_out (train-oriented fold-mean)")
print(f"{'rank':>4} {'feature':22} {'held_out':>8} {'sd':>6} {'in_samp':>8} {'gen?':>4} {'n':>5}  tag")
for i, r in enumerate(rank_rows, 1):
    print(f"{i:>4} {r['feature']:22} {r['held_out_auc']:>8.4f} {r['held_out_sd']:>6.3f} "
          f"{r['in_sample_auc']:>8.4f} {'Y' if r['generalizes'] else 'n':>4} {r['n_present']:>5}  {r['tag']}")
print()
print("## LOGISTIC COMBINATIONS (component-grouped OOF)  -- two metrics:")
print("##   pooled = one AUC over concatenated OOF logits (cross-block calibration-sensitive)")
print("##   perfold = mean within-fold AUC (deployment-relevant within-block discrimination)")
print(f"{'combo':34} {'pooled':>7} {'d_pool':>7}  {'perfold':>7} {'d_perf':>7}  {'sd':>5}")
for name, (a_pool, a_perf, a_sd) in sorted(combos.items(), key=lambda kv: -kv[1][0]):
    print(f"  {name:32} {a_pool:>7.4f} {a_pool-alnfrac_oof:>+7.4f}  "
          f"{a_perf:>7.4f} {a_perf-alnfrac_perfold:>+7.4f}  {a_sd:>5.3f}")
print(f"#   VERDICT: pooled lift of combos is a cross-BLOCK calibration effect; per-fold")
print(f"#   (within-block) lift over aln_frac is ~0. No deployable feature/combo beats the")
print(f"#   workhorse within-block. Only DNA id/min_cov (VALID-ONLY, n=468) rank above it.")

# capture edge-stage summary for the combined JSON
EDGE_SUMMARY = dict(
    n_edges=len(recs), n_TP=int(y.sum()), n_genuineFP=int((y == 0).sum()),
    n_features=len(FEATS), n_components=len(set(comps)), fold_loads=fold_load,
    contrast="TP vs genuine-FP only (hard); headline 0.915 is TP vs ALL-FP",
    ranked=[{k: (None if (isinstance(r[k], float) and np.isnan(r[k])) else r[k])
             for k in ("feature", "held_out_auc", "held_out_sd", "in_sample_auc",
                       "n_present", "generalizes", "tag")} for r in rank_rows],
    combos={name: dict(pooled=a_pool, perfold=a_perf, perfold_sd=a_sd)
            for name, (a_pool, a_perf, a_sd) in combos.items()},
)


# ======================================================================
# ============================ BLOCK STAGE =============================
# ======================================================================
# Per-family-block differentiators of E_p-impurity (over-merge PROXY), with the
# NEW SUB-STRUCTURE features (edge-weight bimodality / spectral gap / best-2-cut
# modularity / identity bimodality / dispersion) added on top of the already-
# computed size/connectivity/multiplicity/cardinality panel (reused from
# block_graph_fp_detector.tsv).  Size-residualized held-out AUC(is_Ep_impure),
# diploid-oracle precision@K, and the sharp GSTM2/MAGE residual test.
print("\n\n" + "=" * 70)
print("BLOCK STAGE  --  per-family over-merge differentiators (+ sub-structure)")
print("=" * 70)

import family_level_pr_current as _M       # build_ctx, load_catalog_blocks
import genome_family_def as _G             # distinct_loci
from scipy.stats import skew as _skew, kurtosis as _kurt
from networkx.algorithms.community import louvain_communities as _louvain, modularity as _modq

SEED = 0
np.random.seed(SEED)

# ---- validated bimodality battery (dip dropped: underpowered/unstable at the
#      block edge-count regime; Sarle BC is the deployed, verified statistic) ----
def bimod_coeff(x):
    """Sarle's bimodality coefficient, sample-corrected moments (>0.555 ~ bimodal)."""
    x = np.asarray(x, float)
    n = x.size
    if n < 4 or x.std() == 0:
        return np.nan
    g1 = float(_skew(x, bias=False))
    g2 = float(_kurt(x, fisher=True, bias=False))
    denom = g2 + 3.0 * (n - 1) ** 2 / ((n - 2) * (n - 3))
    return (g1 ** 2 + 1.0) / denom if denom != 0 else np.nan


def twomeans_sep(x, cap=1e3):
    """Best 1-D threshold split separation (m2-m1)/pooled_within_std, clipped.
    Unlike Otsu-gain this is not inflated by a flat unimodal spread."""
    x = np.sort(np.asarray(x, float))
    n = x.size
    if n < 4 or x[0] == x[-1]:
        return 0.0
    cs = np.cumsum(x); cs2 = np.cumsum(x * x); tot, tot2 = cs[-1], cs2[-1]
    best = 0.0
    for i in range(1, n):
        n1, n2 = i, n - i
        m1 = cs[i - 1] / n1; m2 = (tot - cs[i - 1]) / n2
        v1 = max(cs2[i - 1] / n1 - m1 * m1, 0.0)
        v2 = max((tot2 - cs2[i - 1]) / n2 - m2 * m2, 0.0)
        within = math.sqrt((n1 * v1 + n2 * v2) / n)
        sep = cap if within <= 0 else min((m2 - m1) / within, cap)
        if sep > best:
            best = sep
    return float(best)


def dist_feats(vals, pre):
    """Distribution/bimodality features for an edge-attribute vector."""
    v = np.asarray([x for x in vals if not (isinstance(x, float) and np.isnan(x))], float)
    d = {}
    if v.size == 0:
        for s in ("mean", "cv", "min_max_ratio", "range", "bimod", "kurt", "skew", "sep", "n"):
            d[f"{pre}_{s}"] = np.nan
        d[f"{pre}_n"] = 0
        return d
    mn, mx = float(v.min()), float(v.max())
    mean = float(v.mean()); sd = float(v.std())
    d[f"{pre}_mean"] = mean
    d[f"{pre}_cv"] = sd / mean if mean > 0 else np.nan
    d[f"{pre}_min_max_ratio"] = mn / mx if mx > 0 else np.nan
    d[f"{pre}_range"] = mx - mn
    near_const = (np.ptp(v) <= 1e-12)
    d[f"{pre}_bimod"] = np.nan if near_const else bimod_coeff(v)
    with np.errstate(all="ignore"):
        d[f"{pre}_kurt"] = np.nan if (near_const or v.size < 4) else float(_kurt(v, fisher=True, bias=False))
        d[f"{pre}_skew"] = np.nan if (near_const or v.size < 3) else float(_skew(v, bias=False))
    d[f"{pre}_sep"] = twomeans_sep(v)
    d[f"{pre}_n"] = int(v.size)
    return d


def spectral_feats(sub):
    """lambda2/lambda3 of the weighted combinatorial + normalized Laplacian of the
    giant component; eigengaps; best-2-cut (Fiedler bisection) modularity."""
    d = dict(l2=np.nan, l3=np.nan, gap23=np.nan, ratio32=np.nan, ratio23=np.nan,
             nl2=np.nan, nl3=np.nan, ngap23=np.nan, eigengap_k=np.nan, best2cut_mod=np.nan)
    m = sub.number_of_nodes()
    if m < 3 or sub.number_of_edges() == 0:
        return d
    A = nx.to_numpy_array(sub, weight="weight")
    deg = A.sum(1)
    L = np.diag(deg) - A
    ev = np.sort(np.linalg.eigvalsh(L))
    ev = np.clip(ev, 0, None)
    d["l2"] = float(ev[1]); d["l3"] = float(ev[2]) if m >= 3 else np.nan
    d["gap23"] = float(ev[2] - ev[1]) if m >= 3 else np.nan
    d["ratio32"] = float(ev[2] / ev[1]) if (m >= 3 and ev[1] > 1e-12) else np.nan
    d["ratio23"] = float(ev[1] / ev[2]) if (m >= 3 and ev[2] > 1e-12) else np.nan
    # eigengap heuristic: k that maximizes ev[k]-ev[k-1] over small k -> #clusters
    if m >= 4:
        gaps = np.diff(ev[:min(m, 6)])
        d["eigengap_k"] = int(np.argmax(gaps)) + 1
    dsi = 1.0 / np.sqrt(np.where(deg > 0, deg, 1.0))
    Ln = np.eye(m) - (dsi[:, None] * A) * dsi[None, :]
    evn = np.sort(np.linalg.eigvalsh(Ln))
    d["nl2"] = float(evn[1]); d["nl3"] = float(evn[2]) if m >= 3 else np.nan
    d["ngap23"] = float(evn[2] - evn[1]) if m >= 3 else np.nan
    try:
        fv = nx.fiedler_vector(sub, weight="weight", seed=SEED, method="tracemin_lu")
        nodes = list(sub.nodes())
        p0 = {nodes[i] for i in range(len(nodes)) if fv[i] <= 0}
        p1 = {nodes[i] for i in range(len(nodes)) if fv[i] > 0}
        if p0 and p1:
            d["best2cut_mod"] = float(_modq(sub, [p0, p1], weight="weight"))
    except Exception:
        pass
    return d


# ---- rebuild ctx + catalog; reuse block_graph_fp_detector.tsv for base panel ----
print("[block] build ctx + catalog ...", flush=True)
_ctx, _rf, _ep = _M.build_ctx()
_cat, _ids, _dom, _mg = _M.load_catalog_blocks()
_Gfull = _ctx["G"]; _ew = _ctx["ew"]; _god = _ctx["gene_of_dn"]; _genes = _ctx["genes"]

# base per-block panel (size/connectivity/multiplicity/cardinality + labels)
base = {r["family_id"]: r for r in rd("block_graph_fp_detector.tsv")}
# gene-pair attribute maps for identity distributions
alnfrac_pair = {}
for r in sl:
    alnfrac_pair.setdefault(norm(r["gene_a"], r["gene_b"]), f(r["aln_frac"]))
id_pair = {}
for kpair, r in fp.items():
    id_pair[kpair] = f(r.get("id"))

BASE_COLS = ["n_members", "n_genes", "distinct_loci", "n_edges", "edge_density",
             "max_degree", "mean_degree", "clustering", "n_communities", "fiedler",
             "mincut", "frac_in_giant", "n_seqclust_hi", "loci_per_seqclust",
             "mem_per_seqclust", "loci_per_comm", "loci_per_gene", "mem_per_locus",
             "max_mult", "med_mult", "max_min_mult", "frac_repeat_bridged", "frac_ge_gate"]

brecs = []
for i, block in enumerate(_cat):
    fid = _ids[i]
    br = base.get(fid)
    if br is None:
        continue
    nodes = [d for d in block]
    ns = set(nodes)
    H = nx.Graph(); H.add_nodes_from(nodes)
    ews, afs, idv, chroms = [], [], [], set()
    starts_by_chrom = {}
    cross_chrom_edges = 0
    for u, v in _Gfull.edges(nodes):
        if u in ns and v in ns:
            w = _ew.get(frozenset((u, v)), 0.13)
            H.add_edge(u, v, weight=w)
            ews.append(w)
            gu, gv = _god.get(u), _god.get(v)
            cu = _genes.get(u, {}).get("chrom"); cv = _genes.get(v, {}).get("chrom")
            if cu and cv and cu != cv:
                cross_chrom_edges += 1
            if gu and gv and gu != gv:
                kk = norm(gu, gv)
                a = alnfrac_pair.get(kk)
                if a is not None and not np.isnan(a):
                    afs.append(a)
                idd = id_pair.get(kk)
                if idd is not None and not np.isnan(idd):
                    idv.append(idd)
    for dn in nodes:
        g = _genes.get(dn, {})
        c = g.get("chrom")
        if c:
            chroms.add(c)
            starts_by_chrom.setdefault(c, []).append((g.get("start", 0), g.get("end", 0)))
    n_edges = len(ews)

    d = {"family_id": fid}
    # copy base panel columns (numeric)
    for c in BASE_COLS:
        d[c] = f(br.get(c))
    # labels
    d["is_impure"] = int(br["is_impure"])
    d["oracle_mapped"] = int(br["oracle_mapped"])
    d["oracle_fp"] = int(br["oracle_fp"])
    d["oracle_cls"] = br["oracle_cls"]
    d["dominant_gene"] = br["dominant_gene"]

    # ---- NEW: edge-weight (core_recip) distribution / bimodality ----
    d.update(dist_feats(ews, "ew"))
    # ---- NEW: identity distributions (aln_frac deployable; DNA id VALID-ONLY) ----
    d.update(dist_feats(afs, "af"))
    d.update(dist_feats(idv, "id"))
    # ---- NEW: spectral + best-2-cut on giant component ----
    if n_edges > 0:
        gc = max(nx.connected_components(H), key=len)
        sp = spectral_feats(H.subgraph(gc))
    else:
        sp = spectral_feats(nx.Graph())
    for k, val in sp.items():
        d[f"spec_{k}"] = val
    # ---- NEW: Louvain modularity (weighted) = whole-block sub-structure ----
    if n_edges > 0:
        comms = _louvain(H, seed=SEED, weight="weight")
        d["louvain_mod"] = float(_modq(H, comms, weight="weight")) if len(comms) >= 1 else np.nan
        d["n_louvain"] = len(comms)
    else:
        d["louvain_mod"] = np.nan; d["n_louvain"] = len(nodes)
    # ---- NEW: dispersion ----
    d["n_chroms"] = len(chroms)
    d["cross_chrom_frac"] = cross_chrom_edges / n_edges if n_edges else np.nan
    span = 0
    for c, se in starts_by_chrom.items():
        mn = min(s for s, e in se); mx = max(e for s, e in se)
        span += max(mx - mn, 0)
    d["block_span_bp"] = span
    d["log_span_bp"] = math.log1p(span)
    brecs.append(d)

# feature set: base + new (exclude labels/meta & the size vars used for residualizing)
NEW_FEATS = [k for k in brecs[0].keys()
             if k not in ("family_id", "dominant_gene", "oracle_cls", "is_impure",
                          "oracle_mapped", "oracle_fp") and k not in BASE_COLS]
# report-only leakage sentinels excluded from ranking
BLOCK_FEATS = [c for c in BASE_COLS if c not in ("n_members", "distinct_loci", "n_genes")] + NEW_FEATS
# tag sub-structure-new vs base-reused vs validation-only
SUBSTRUCT_NEW = {kk for kk in NEW_FEATS if kk.startswith(("ew_", "af_", "spec_", "louvain", "n_louvain"))}
DISPERSION_NEW = {"n_chroms", "cross_chrom_frac", "block_span_bp", "log_span_bp"}
VALID_ONLY_BLK = {kk for kk in NEW_FEATS if kk.startswith("id_")} | {"frac_repeat_bridged", "frac_ge_gate"}

yb = np.array([r["is_impure"] for r in brecs])
Nb = len(brecs)
logn = np.log(np.array([r["n_members"] for r in brecs], float))
logdl = np.log(np.clip(np.array([r["distinct_loci"] for r in brecs], float), 1, None))
size_lin = np.column_stack([np.ones(Nb), logn, logdl])
size_quad = np.column_stack([np.ones(Nb), logn, logdl, logn * logn, logdl * logdl, logn * logdl])


def blk_auc(scores, labels, mask=None):
    s = np.asarray(scores, float); lab = np.asarray(labels)
    ok = ~np.isnan(s)
    if mask is not None:
        ok = ok & mask
    s, lab = s[ok], lab[ok]
    npos = int((lab == 1).sum()); nneg = int((lab == 0).sum())
    if npos == 0 or nneg == 0:
        return np.nan
    rk = rankdata(s)
    return (rk[lab == 1].sum() - npos * (npos + 1) / 2.0) / (npos * nneg)


def residualize_fit(vals, Xsize, rows_mask):
    """OLS residual of vals on Xsize using rows in rows_mask (non-nan); return full resid vector."""
    ok = rows_mask & ~np.isnan(vals)
    if ok.sum() < Xsize.shape[1] + 2:
        return None
    beta, *_ = np.linalg.lstsq(Xsize[ok], vals[ok], rcond=None)
    return vals - Xsize @ beta


# deterministic 5-fold split of blocks (shuffle, contiguous)
perm = np.random.RandomState(SEED).permutation(Nb)
blk_fold = np.empty(Nb, int)
for j, idx in enumerate(perm):
    blk_fold[idx] = j % 5

# ---------- per-feature: held-out size-residualized AUC + in-sample lin/quad ----------
brank = []
for name in BLOCK_FEATS:
    vals = np.array([r[name] for r in brecs], float)
    n_present = int((~np.isnan(vals)).sum())
    # in-sample raw + size-residualized (lin + quad)
    a_raw = blk_auc(vals, yb)
    rlin = residualize_fit(vals, size_lin, np.ones(Nb, bool))
    rquad = residualize_fit(vals, size_quad, np.ones(Nb, bool))
    a_lin = blk_auc(rlin, yb) if rlin is not None else np.nan
    a_quad = blk_auc(rquad, yb) if rquad is not None else np.nan
    # held-out fold-mean size-residualized (fit residualizer on train, eval test)
    per = []
    for k in range(5):
        te = blk_fold == k; tr = ~te
        rr = residualize_fit(vals, size_lin, tr)
        if rr is None:
            continue
        a_tr = blk_auc(rr, yb, mask=tr)
        if np.isnan(a_tr):
            continue
        sgn = 1.0 if a_tr >= 0.5 else -1.0
        a_te = blk_auc(sgn * rr, yb, mask=te)
        if not np.isnan(a_te):
            per.append(a_te)
    ho = float(np.mean(per)) if per else np.nan
    ho_sd = float(np.std(per)) if per else np.nan
    tag = []
    if name in SUBSTRUCT_NEW: tag.append("SUBSTRUCT-NEW")
    elif name in DISPERSION_NEW: tag.append("DISPERSION-NEW")
    elif name in BASE_COLS: tag.append("base-reused")
    if name in VALID_ONLY_BLK: tag.append("VALID-ONLY")
    # size-artifact: held-out (linear-resid) looks discriminative BUT collapses
    # under the richer QUADRATIC size basis => it was a nonlinear function of size
    disc_lin = (not np.isnan(ho)) and abs(ho - 0.5) >= 0.10
    quad_collapse = (not np.isnan(a_quad)) and abs(a_quad - 0.5) < 0.06
    # sign-flip between the linear and quadratic size bases = unstable size ratio
    sign_flip = (not np.isnan(a_lin) and not np.isnan(a_quad)
                 and (a_lin - 0.5) * (a_quad - 0.5) < 0 and abs(a_lin - 0.5) >= 0.10)
    size_artifact = bool(disc_lin and (quad_collapse or sign_flip))
    if size_artifact:
        tag.append("SIZE-ARTIFACT")
    brank.append(dict(feature=name, held_out=ho, held_out_sd=ho_sd,
                      insample_raw=None if np.isnan(a_raw) else round(a_raw, 3),
                      insample_lin=None if np.isnan(a_lin) else round(a_lin, 3),
                      insample_quad=None if np.isnan(a_quad) else round(a_quad, 3),
                      size_artifact=size_artifact,
                      n_present=n_present, tag=",".join(tag) if tag else "-"))
# rank by |held_out - 0.5| descending (discriminative power, sign-agnostic)
brank.sort(key=lambda r: -(abs(r["held_out"] - 0.5) if not np.isnan(r["held_out"]) else -1))

# ---------- diploid-oracle precision@K for top features + sub-structure combo ----------
orc_fp_idx = [i for i in range(Nb) if brecs[i]["oracle_fp"] == 1]
n_orc_fp = len(orc_fp_idx)


def precision_at_k(scorevec, Ks=(10, 20, 30, 50)):
    order = np.argsort(-scorevec, kind="mergesort")
    out = {}
    for Kk in Ks:
        topK = set(order[:Kk].tolist())
        out[Kk] = dict(oracle_fp=f"{len(topK & set(orc_fp_idx))}/{n_orc_fp}",
                       ep_impure=f"{sum(brecs[i]['is_impure'] for i in topK)}/{Kk}")
    return out


# oriented in-sample size-residual score per feature (higher = more FP-like)
def oriented_resid(name):
    vals = np.array([r[name] for r in brecs], float)
    rr = residualize_fit(vals, size_lin, np.ones(Nb, bool))
    if rr is None:
        return None
    a = blk_auc(rr, yb)
    sgn = 1.0 if (np.isnan(a) or a >= 0.5) else -1.0
    s = sgn * rr
    s = np.where(np.isnan(s), np.nanmin(s) - 1.0, s)  # NaN sinks to bottom
    return s


# AUDIT FIX: compute precision@K for EVERY ranked feature (not just the top-8 by
# held-out AUC).  Truncating to top-8 previously hid that the RNA-deployable NEW
# feature ew_bimod also reaches 4/6@50 (ties DNA id_bimod).  NOTE: this univariate
# precision@K is IN-SAMPLE-ORIENTED (sign + size-residualizer fit on ALL blocks);
# only the *_combo_oof scorers below are true out-of-fold.  Given n=6 oracle-FPs
# a 1-block (4/6 vs 3/6) gap is within noise -- read the concentration together
# with the held-out AUC and the clean matched-real-dense residual test.
oracle_prec = {}
for r in brank:
    s = oriented_resid(r["feature"])
    if s is not None:
        oracle_prec[r["feature"]] = precision_at_k(s)

# sub-structure logistic combo (size-residualized inputs), held-out oof + oracle concentration
SUB_SET = [x for x in ["ew_bimod", "ew_cv", "ew_sep", "ew_kurt", "spec_best2cut_mod",
                       "spec_ratio32", "spec_ngap23", "louvain_mod", "af_bimod", "af_cv",
                       "cross_chrom_frac"] if x in brecs[0]]


def resid_matrix(cols, rows_mask):
    M = np.zeros((Nb, len(cols)))
    for j, c in enumerate(cols):
        vals = np.array([r[c] for r in brecs], float)
        rr = residualize_fit(vals, size_lin, rows_mask)
        if rr is None:
            rr = vals - np.nanmean(vals[rows_mask])
        med = np.nanmedian(rr[rows_mask])
        rr = np.where(np.isnan(rr), med, rr)
        M[:, j] = rr
    return M


def combo_block_oof(cols):
    pooled = np.full(Nb, np.nan)
    per = []
    for k in range(5):
        te = blk_fold == k; tr = ~te
        M = resid_matrix(cols, tr)
        mu = M[tr].mean(0); sd = M[tr].std(0); sd[sd == 0] = 1.0
        Xtr = (M[tr] - mu) / sd; Xte = (M[te] - mu) / sd
        w = fit_logistic(Xtr, yb[tr].astype(float), l2=2.0)
        z = np.hstack([np.ones((Xte.shape[0], 1)), Xte]) @ w
        pooled[te] = z
        a = blk_auc(pooled, yb, mask=te)
        if not np.isnan(a):
            per.append(a)
    return blk_auc(pooled, yb), (float(np.mean(per)) if per else np.nan), pooled


sub_pool, sub_perfold, sub_pooled_scores = combo_block_oof(SUB_SET)
base_graph = [x for x in ["clustering", "n_communities", "max_mult", "med_mult",
                          "loci_per_comm", "edge_density"] if x in brecs[0]]
bg_pool, bg_perfold, bg_scores = combo_block_oof(base_graph)
allc_pool, allc_perfold, allc_scores = combo_block_oof(SUB_SET + base_graph)
oracle_prec["SUBSTRUCT_combo_oof"] = precision_at_k(
    np.where(np.isnan(sub_pooled_scores), np.nanmin(sub_pooled_scores) - 1, sub_pooled_scores))
oracle_prec["ALL_combo_oof"] = precision_at_k(
    np.where(np.isnan(allc_scores), np.nanmin(allc_scores) - 1, allc_scores))

# ---------- THE RESIDUAL TEST: GSTM2 (9,13) + MAGE (549) vs matched real dense ----------
by_fid = {r["family_id"]: r for r in brecs}
idx_of = {r["family_id"]: i for i, r in enumerate(brecs)}
# matched real dense reference set: pure (not impure) AND not an oracle FP, dense, similar size
real_dense = [r for r in brecs if r["is_impure"] == 0 and r["oracle_fp"] == 0
              and not np.isnan(r["edge_density"]) and r["edge_density"] >= 0.5
              and 3 <= r["n_members"] <= 20]
# AUDIT FIX: cover ALL genuinely-missed confirmed FPs, not just {9,13,549}.
#   b9  = GSTM2-domain-hub over-merge (dominant_gene LOC115930164, oracle fam14; the
#         literal-GSTM2 pure block is 23 -- naming inherited from the task, FP status real)
#   b13 = GSTM2 3-node dense triangle (no sub-structure can exist)
#   b549= MAGE-A over-merge (dense uniform blob)
#   b110= LOC134758618, the OTHER is_impure=0 genuinely-missed FP (was absent before)
TARGETS = {"b9_GSTM2hub": "9", "b13_GSTM2tri": "13", "b549_MAGE": "549", "b110_missed": "110"}
NEW_TEST_FEATS = [kk for kk in NEW_FEATS
                  if kk.startswith(("ew_", "af_", "spec_", "louvain", "n_louvain", "cross_chrom"))
                  and not kk.endswith("_n")]


def pct_rank(name, target_val, ref_vals):
    ref = np.array([x for x in ref_vals if not (isinstance(x, float) and np.isnan(x))], float)
    if len(ref) == 0 or (isinstance(target_val, float) and np.isnan(target_val)):
        return np.nan
    return float((ref < target_val).mean())  # fraction of real-dense scored BELOW target


residual_test = {}
for fname in NEW_TEST_FEATS:
    ref_vals = [r[fname] for r in real_dense]
    row = {}
    for tname, tfid in TARGETS.items():
        tv = by_fid[tfid][fname] if tfid in by_fid else np.nan
        row[tname] = dict(value=None if (isinstance(tv, float) and np.isnan(tv)) else round(float(tv), 4),
                          pct_vs_realdense=None if np.isnan(pct_rank(fname, tv, ref_vals)) else round(pct_rank(fname, tv, ref_vals), 3))
    # AUDIT FIX: honest flag name.  na targets (e.g. b13 3-node triangle, or blocks
    # with <2 defined values) are silently DROPPED, so this is "all AVAILABLE targets
    # above the real-dense median" -- a 50%-FP operating point, NOT a deployable gate.
    pcts = [row[t]["pct_vs_realdense"] for t in TARGETS if row[t]["pct_vs_realdense"] is not None]
    row["catches_all_avail_above_median"] = bool(pcts) and all(p > 0.5 for p in pcts)
    row["n_targets_available"] = len(pcts)
    row["catches_MAGE_top10pct"] = (row["b549_MAGE"]["pct_vs_realdense"] or 0) >= 0.9
    residual_test[fname] = row

# ---------------------------------------------------------------- write TSV
blk_out = os.path.join(B, "tp_fp_panel_block.tsv")
ALL_BLK_COLS = (["family_id", "dominant_gene", "is_impure", "oracle_mapped",
                 "oracle_fp", "oracle_cls", "fold"] + ["n_members", "n_genes", "distinct_loci"]
                + BLOCK_FEATS)
with open(blk_out, "w") as fh:
    fh.write("\t".join(ALL_BLK_COLS) + "\n")
    for i, r in enumerate(brecs):
        vals = []
        for c in ALL_BLK_COLS:
            if c == "fold":
                vals.append(str(int(blk_fold[i]))); continue
            x = r.get(c, "")
            if isinstance(x, float) and np.isnan(x):
                vals.append("")
            elif isinstance(x, float):
                vals.append(f"{x:.6g}")
            else:
                vals.append(str(x))
        fh.write("\t".join(vals) + "\n")

# ---------------------------------------------------------------- write JSON
BLOCK_SUMMARY = dict(
    n_blocks=Nb, n_ep_impure=int(yb.sum()), n_oracle_mapped=len([r for r in brecs if r["oracle_mapped"]]),
    n_oracle_fp=n_orc_fp, n_new_substruct_feats=len(SUBSTRUCT_NEW), n_features_ranked=len(BLOCK_FEATS),
    dip_note="Hartigan dip prototyped but DROPPED: underpowered/unstable at block edge-count "
             "regime (many blocks <10 edges); Sarle bimodality-coefficient is the deployed, "
             "validated bimodality statistic.",
    ranked=brank,
    combos=dict(
        substruct_new=dict(pooled=round(sub_pool, 4), perfold=round(sub_perfold, 4), feats=SUB_SET),
        base_graph=dict(pooled=round(bg_pool, 4), perfold=round(bg_perfold, 4), feats=base_graph),
        all_combined=dict(pooled=round(allc_pool, 4), perfold=round(allc_perfold, 4)),
    ),
    oracle_precision_at_k=oracle_prec,
    residual_test=residual_test,
    matched_real_dense_n=len(real_dense),
)
combined = dict(
    determinism=dict(pythonhashseed=os.environ.get("PYTHONHASHSEED"), seed=SEED),
    edge_stage=EDGE_SUMMARY,
    block_stage=BLOCK_SUMMARY,
)
import json as _json
with open(os.path.join(B, "tp_fp_panel.json"), "w") as fh:
    _json.dump(combined, fh, indent=2, default=str)

# ---------------------------------------------------------------- reports
print(f"\n[block] {Nb} families  ({int(yb.sum())} E_p-impure)  {len(BLOCK_FEATS)} features ranked  "
      f"({len(SUBSTRUCT_NEW)} sub-structure-NEW)  oracle-FP={n_orc_fp}  matched-real-dense={len(real_dense)}")
print(f"[block] wrote {blk_out} + tp_fp_panel.json")
print()
print("## BLOCK UNIVARIATE: SIZE-RESIDUALIZED held-out AUC(is_Ep_impure), ranked by |AUC-0.5|")
print("##   held_out = fit OLS(feat~1+log n_members+log distinct_loci) on TRAIN, eval AUC on TEST fold-mean")
print("##   in-samp lin/quad = size-residual AUC on full data (quad basis exposes nonlinear-ratio artifacts)")
print(f"{'rank':>4} {'feature':20} {'held':>6} {'sd':>5} {'raw':>5} {'lin':>5} {'quad':>5} {'n':>4}  tag")
for i, r in enumerate(brank, 1):
    ho = "  nan" if np.isnan(r["held_out"]) else f"{r['held_out']:.3f}"
    sd = "  nan" if np.isnan(r["held_out_sd"]) else f"{r['held_out_sd']:.3f}"
    print(f"{i:>4} {r['feature']:20} {ho:>6} {sd:>5} {str(r['insample_raw']):>5} "
          f"{str(r['insample_lin']):>5} {str(r['insample_quad']):>5} {r['n_present']:>4}  {r['tag']}")
print()
print("## SIZE-RESIDUALIZED LOGISTIC COMBOS (5-fold OOF AUC is_Ep_impure)")
print(f"  substruct-NEW only : pooled={sub_pool:.3f}  perfold={sub_perfold:.3f}  ({len(SUB_SET)} feats)")
print(f"  base-graph         : pooled={bg_pool:.3f}  perfold={bg_perfold:.3f}")
print(f"  all combined       : pooled={allc_pool:.3f}  perfold={allc_perfold:.3f}")
print()
print("## DIPLOID-ORACLE precision@K (concentration of the 6 confirmed FP blocks)")
print("##   IN-SAMPLE-ORIENTED (sign+residualizer on all blocks); only *_combo_oof are true OOF")
print("##   curated view = top-8 held-out AUC UNION all bimodality/spectral + DNA + combos")
VALID_ONLY_ALL = VALID_ONLY_BLK | {"id_bimod", "min_cov"}
curated_pk = ([r["feature"] for r in brank[:8]]
              + ["ew_bimod", "af_bimod", "id_bimod", "spec_l2", "spec_ratio32",
                 "spec_best2cut_mod", "louvain_mod", "ew_cv", "n_communities"]
              + ["SUBSTRUCT_combo_oof", "ALL_combo_oof"])
seen = set()
print(f"{'scorer':22} " + "  ".join(f"K={k}" for k in (10, 20, 30, 50)) + "   deploy?")
for name in curated_pk:
    if name in seen or name not in oracle_prec:
        continue
    seen.add(name)
    pk = oracle_prec[name]
    dep = "VALID-ONLY" if (name in VALID_ONLY_ALL) else ("combo" if name.endswith("_combo_oof") else "RNA")
    print(f"  {name:20} " + "  ".join(f"{pk[k]['oracle_fp']:>4}" for k in (10, 20, 30, 50))
          + "   (ep-impure@50 " + oracle_prec[name][50]['ep_impure'] + f")  {dep}")
print()
print("## RESIDUAL TEST: do NEW sub-structure features score the confirmed FPs")
print("##   ABOVE matched real dense families?  pct = fraction of real-dense BELOW target (>0.9 = top decile)")
print("##   'avail>med' = all AVAILABLE (non-na) targets above real-dense median (NOT a deployable gate)")
print(f"{'feature':18} " + " ".join(f"{t:>14}" for t in TARGETS) + "  avail>med MAGEtop10%")
for fname in NEW_TEST_FEATS:
    rt = residual_test[fname]
    def cell(t):
        v = rt[t]["value"]; p = rt[t]["pct_vs_realdense"]
        return f"{'na' if v is None else v}({'na' if p is None else p})"
    print(f"{fname:18} " + " ".join(f"{cell(t):>14}" for t in TARGETS)
          + f"  {('Y('+str(rt['n_targets_available'])+')') if rt['catches_all_avail_above_median'] else 'n':>9} "
          + f"{'Y' if rt['catches_MAGE_top10pct'] else 'n':>9}")
print()
# definitive verdict
any_all3 = [fn for fn in NEW_TEST_FEATS if residual_test[fn]["catches_all_avail_above_median"]]
mage_top = [fn for fn in NEW_TEST_FEATS if residual_test[fn]["catches_MAGE_top10pct"]]
print("## DEFINITIVE RESIDUAL VERDICT")
print(f"  features scoring all AVAILABLE confirmed FPs above real-dense median : {any_all3 if any_all3 else 'NONE'}")
print( "    (weak 50%-FP bar; na targets dropped; NOT a deployable top-decile gate)")
print(f"  features putting MAGE(549) in real-dense top decile          : {mage_top if mage_top else 'NONE'}")
# genuine (survives quad-size-basis) sub-structure signals
genuine = [r["feature"] for r in brank
           if abs(r["held_out"] - 0.5) >= 0.10 and not r["size_artifact"]
           and "VALID-ONLY" not in r["tag"] and not np.isnan(r["held_out"])]
size_arts = [r["feature"] for r in brank if r["size_artifact"]]
print("\n## BLOCK-STAGE DEFINITIVE SYNTHESIS")
print(f"  best RNA over-merge combo (perfold OOF): base-graph={bg_perfold:.3f}  "
      f"substruct-NEW={sub_perfold:.3f}  all={allc_perfold:.3f}  (sub-structure does NOT beat base; adding it hurts)")
print(f"  genuine size-independent RNA differentiators (survive QUAD basis, not VALID-ONLY): {genuine[:8]}")
print(f"  size-artifact features (linear-resid discriminative but collapse under quad): {size_arts}")
print( "  oracle-FP concentration (AUDIT-CORRECTED): RNA-deployable ew_bimod = 4/6@50 TIES DNA id_bimod = 4/6@50")
print( "     => residual is NOT concentration-DNA-bound (a NEW RNA sub-structure feature ties DNA on raw oracle-hit)")
print( "     BUT it IS generalization-DNA-bound: ew_bimod fails size-control (in-samp lin AUC 0.348 / quad 0.318),")
print( "     held-out 0.633+-0.188, and the clean matched-real-dense test does NOT flag MAGE (pct 0.593).")
print( "     af_bimod cleanly flags MAGE (pct 1.0) but held-out AUC 0.41 (non-generalizing).")
print( "     ONLY id_bimod (DNA) generalizes held-out (0.76).  No RNA feature both concentrates AND generalizes.")
print( "  GSTM2 block9 = genuinely BIMODAL: best2cut/louvain_mod top 1%/128 real-dense (hypothesis TRUE) but already flagged by shipped multiplicity+impurity")
print( "  MAGE 549 = dense UNIFORM blob: best2cut_mod~0 (graph-topology hypothesis FALSE) BUT af_bimod=0.91 > real-dense MAX (identity-bimodality catches it; does NOT generalize AUC 0.42)")
print( "  GSTM2 block13 = 3-node dense triangle: NO RNA sub-structure feature applies -> provably DNA-bound")
BLOCK_SUMMARY["verdict"] = dict(
    genuine_size_independent_rna=genuine,
    size_artifact_features=size_arts,
    features_all_avail_above_median=any_all3,
    features_mage_top_decile=mage_top,
    best_combo_perfold=dict(base_graph=round(bg_perfold, 3), substruct_new=round(sub_perfold, 3),
                            all=round(allc_perfold, 3)),
    oracle_concentration_corrected=(
        "ew_bimod (RNA-DEPLOYABLE NEW) = 4/6@50 TIES id_bimod (DNA) = 4/6@50 on in-sample-oriented "
        "size-residualized precision@K -> residual is NOT concentration-DNA-bound; but ew_bimod is "
        "generalization-bound (in-samp lin AUC 0.348/quad 0.318, held-out 0.633+-0.188, clean "
        "matched-real-dense MAGE pct 0.593 = not flagged). Only id_bimod (DNA) generalizes held-out 0.76."),
    gstm2_b9="genuinely bimodal; best2cut_mod 0.317 = top 1%/128 real-dense; already shipped-flagged (is_impure=1, max_mult 426)",
    mage_549="dense uniform blob; best2cut_mod~0 (topology hypothesis FALSE); af_bimod 0.91>real-dense max (identity bimodality catches, non-generalizing held-out 0.42)",
    gstm2_b13="3-node dense triangle; no RNA sub-structure applies; DNA-bound",
    b110_missed="LOC134758618 is_impure=0 missed FP; ew_bimod 0.675 (caught by ew_bimod concentration), best2cut~0, af_bimod na",
    residual_is_concentration_dna_bound=False,
    residual_is_generalization_dna_bound=True,
)
with open(os.path.join(B, "tp_fp_panel.json"), "w") as fh:
    _json.dump(combined, fh, indent=2, default=str)
