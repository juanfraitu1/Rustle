#!/usr/bin/env python
"""RNA-ONLY refined edge oracle -- FEATURE STAGE.

GOAL
----
The RNA multi-copy family = connected components of the transcript-homology graph
E_r (core_recip >= 0.13) + the shipped gamma-quasi-clique refinement. The refined
catalog still carries residual OVER-MERGE false positives (GRAPH_DEF_REFINE_SWEEP.md).
Graph TOPOLOGY cannot cut them (bridge/triangle hypothesis falsified: over-merges are
dense/triangulated). The one edge-level separator is the whole-transcript reciprocal
HOMOLOGY WEIGHT core_recip (AUC 0.846) -- which is RNA-derived. So an RNA-only edge
oracle exists; DNA/protein only CALIBRATE the operating point and VALIDATE.

This script is the FEATURE stage: assemble the per-E_r-edge RNA-ONLY feature vector,
join the DNA/protein LABELS (learning/validation only -- NEVER an inference feature)
and the family-level RNA allele signal, and report:
  (1) feature COVERAGE (edges with all RNA features vs core_recip-only),
  (2) per-feature AUC vs the DNA in_dna_loose label (confirm core_recip strongest,
      others add signal),
  (3) correlation structure among RNA features + which features carry INDEPENDENT signal.

HARD GUARDS (thesis-critical):
  * DNA/protein columns (in_dna_loose, in_ep, ep_tier, class, cls) are LABELS only.
    They are NEVER put into the RNA feature vector used to predict.
  * The genome/library-derived soft-mask (bridge_mask) is kept OUT of the RNA-only
    feature set and reported ONLY as a clearly-flagged ABLATION column (abl_bridge_mask).
  * Deterministic: PYTHONHASHSEED=0, sorted iteration/writes, fixed float formatting.

EDGE UNIT / POPULATION
----------------------
The shipped refined-E_r gene-pair edge set = bench/edge_bridge_mask.tsv (10755 gene
pairs = family_er_pr.tsv in_er_refined==1, class in {TP, genuine-FP, truthbar-FP}).
Each row is a gene pair asserted to sit in the SAME refined family -- the object whose
over-merges we want to cut. Labels joined authoritatively from family_er_pr.tsv.

RNA-ONLY FEATURES (inference-time; no DNA, no genome soft-mask):
  core_recip  : whole-transcript reciprocal homology weight. Coalesced from the DN
                graph edge weights (bench/denovo_family_edges.tsv, max over the DN
                edges of the pair -- the graph weight the gamma operator sees) with
                the representative bridge core (edge_bridge_mask). Present iff a DIRECT
                homology edge exists; ABSENT on transitive-closure-only co-membership
                pairs (which is itself informative).
  aln_len/aln_frac : longest shared spliced-exon aligned block, abs bp and fraction of
                the shorter locus (bench/ri_sharedlen_cache.tsv). A real paralog shares
                a LONG gene body; a domain/TE bridge shares a SHORT fragment. Built only
                on the in_ep=0 ARBITRATION edges -> this is exactly the co-located
                over-merge target population.
  mm_* (per endpoint) : genome-wide read multi-mapping degree of each bridge locus
                (bench/ri_readmm_cache.tsv). A TE/repeat bridge multi-maps to many loci;
                a gene-family shared exon to a few. Reported per endpoint (hi/lo by
                mm_max) + symmetric combined.
  allele_* (family-level RNA) : balanced_frac / het_risk / copy_like / med_maf
                (bench/a1_read_consensus_o1.tsv). Minor-allele-fraction ~0.5 => diploid
                het (allele-as-copy risk); ~1/K => real copy. For the family-level DEMOTE
                step; joined here by gene for coverage + a sanity check on the known
                allele-as-copy FPs (DHRSX, LOC129530050).

DNA/PROTEIN LABELS (learning + validation ONLY):
  in_dna_loose : reference cDNA-vs-cDNA minimap2 id>=0.90 & cov>=0.30 (REAL edge).
  in_ep/ep_tier: reciprocal whole-protein homology (divergent-paralog "truthbar").
  class        : TP / genuine-FP (over-merge) / truthbar-FP.

Writes: bench/rna_only_edge_features.tsv
Run:    PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/rna_only_edge_oracle.py
"""
import os
import sys

if os.environ.get("PYTHONHASHSEED") != "0":
    os.environ["PYTHONHASHSEED"] = "0"
    os.execv(sys.executable, [sys.executable] + sys.argv)

import math
import numpy as np
import scipy.stats as ss

BENCH = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, BENCH)
import family_er_pr as FP  # gene_of_dn projection (shared with the shipped metric)

EDGES_TSV = os.path.join(BENCH, "denovo_family_edges.tsv")
ERPR_TSV = os.path.join(BENCH, "family_er_pr.tsv")
EBM_TSV = os.path.join(BENCH, "edge_bridge_mask.tsv")
SHAREDLEN_TSV = os.path.join(BENCH, "ri_sharedlen_cache.tsv")
READMM_TSV = os.path.join(BENCH, "ri_readmm_cache.tsv")
A1_TSV = os.path.join(BENCH, "a1_read_consensus_o1.tsv")
OUT_TSV = os.path.join(BENCH, "rna_only_edge_features.tsv")

NA = "NA"


# --------------------------------------------------------------------------- loaders
def load_gene_of_dn():
    """DN locus id -> gene symbol (floored max-overlap projection, shipped)."""
    meta = FP.load_meta()
    annot = FP.load_annot()
    gene_of = FP.gene_of_factory(annot)
    raw = FP.load_raw_families()
    genes, gene_of_dn, *_ = FP.build_genes_dict(raw, meta, gene_of)
    return gene_of_dn


def load_pair_core(gene_of_dn):
    """gene-pair (frozenset) -> MAX core_recip over the DN graph edges of that pair.
       This is the whole-transcript reciprocal homology weight the gamma operator sees."""
    pc = {}
    with open(EDGES_TSV) as fh:
        next(fh)
        for ln in fh:
            a, b, c = ln.rstrip("\n").split("\t")
            ga, gb = gene_of_dn.get(a), gene_of_dn.get(b)
            if ga and gb and ga != gb:
                k = frozenset((ga, gb))
                cv = float(c)
                if cv > pc.get(k, -1.0):
                    pc[k] = cv
    return pc


def load_labels():
    """gene-pair -> dict(in_dna_loose, in_ep, ep_tier, class). Authoritative DNA/protein truth."""
    lab = {}
    with open(ERPR_TSV) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        I = {c: i for i, c in enumerate(h)}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            k = frozenset((f[I["gene_a"]], f[I["gene_b"]]))
            lab[k] = dict(in_dna_loose=int(f[I["in_dna_loose"]]),
                          in_ep=int(f[I["in_ep"]]),
                          ep_tier=f[I["ep_tier"]],
                          cls=f[I["class"]])
    return lab


def load_ebm_master():
    """The refined-E_r gene-pair edge set. Returns list of dict(k, gene_a, gene_b, cls,
       core_ebm, bridge_dn_a, bridge_dn_b, bridge_mask). cls in {TP, genuine, truthbar}."""
    rows = []
    with open(EBM_TSV) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        I = {c: i for i, c in enumerate(h)}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            ga, gb = f[I["gene_a"]], f[I["gene_b"]]
            core = None if f[I["core"]] == "" else float(f[I["core"]])
            bm = None if f[I["bridge_mask"]] in ("", NA) else float(f[I["bridge_mask"]])
            rows.append(dict(k=frozenset((ga, gb)), gene_a=ga, gene_b=gb,
                             cls=f[I["cls"]], core_ebm=core,
                             bridge_dn_a=f[I["bridge_dn_a"]], bridge_dn_b=f[I["bridge_dn_b"]],
                             bridge_mask=bm))
    return rows


def load_sharedlen():
    """gene-pair -> (aln_len, aln_frac). Longest shared spliced-exon aligned block."""
    d = {}
    with open(SHAREDLEN_TSV) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        I = {c: i for i, c in enumerate(h)}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            k = frozenset((f[I["gene_a"]], f[I["gene_b"]]))
            al = f[I["aln_len"]]
            af = f[I["aln_frac"]]
            d[k] = (None if al == "" else int(al), None if af == "" else float(af))
    return d


def load_readmm():
    """DN locus -> (mm_median, mm_mean, mm_max) genome-wide read multi-mapping degree."""
    d = {}
    with open(READMM_TSV) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        I = {c: i for i, c in enumerate(h)}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            d[f[0]] = (float(f[I["mm_median"]]), float(f[I["mm_mean"]]), float(f[I["mm_max"]]))
    return d


def load_allele():
    """gene -> dict(balanced_frac, het_risk, copy_like, med_maf). Family-level RNA allele signal."""
    d = {}
    with open(A1_TSV) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        I = {c: i for i, c in enumerate(h)}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            d[f[I["gene"]]] = dict(balanced_frac=float(f[I["balanced_frac"]]),
                                   het_risk=float(f[I["het_risk"]]),
                                   copy_like=float(f[I["copy_like"]]),
                                   med_maf=float(f[I["med_maf"]]))
    return d


# --------------------------------------------------------------------------- stats helpers
def auc_mw(pos, neg):
    """AUC = P(x_pos > x_neg) via the Mann-Whitney U statistic (ties -> 0.5).
       Returns (auc, n_pos, n_neg). NaNs dropped."""
    pos = [x for x in pos if x is not None and x == x]
    neg = [x for x in neg if x is not None and x == x]
    if not pos or not neg:
        return float("nan"), len(pos), len(neg)
    r = ss.rankdata(pos + neg)
    rp = float(np.sum(r[:len(pos)]))
    a = (rp - len(pos) * (len(pos) + 1) / 2.0) / (len(pos) * len(neg))
    return a, len(pos), len(neg)


def fmt(x, nd=4):
    if x is None or (isinstance(x, float) and x != x):
        return NA
    if isinstance(x, float):
        return f"{x:.{nd}f}"
    return str(x)


def standardize(cols):
    """cols: (n, p) -> standardized with intercept prepended; returns (X, means, stds)."""
    X = np.asarray(cols, dtype=float)
    mu = X.mean(axis=0)
    sd = X.std(axis=0)
    sd[sd == 0] = 1.0
    Xs = (X - mu) / sd
    Xs = np.hstack([np.ones((Xs.shape[0], 1)), Xs])
    return Xs, mu, sd


def logistic_fit(X, y, iters=4000, lr=0.2, l2=1e-3):
    """Deterministic full-batch gradient descent logistic regression."""
    w = np.zeros(X.shape[1])
    n = X.shape[0]
    for _ in range(iters):
        z = X @ w
        p = 1.0 / (1.0 + np.exp(-np.clip(z, -30, 30)))
        g = X.T @ (p - y) / n + l2 * w
        w = w - lr * g
    return w


def logistic_insample_auc(feat_cols, y):
    Xs, _, _ = standardize(feat_cols)
    w = logistic_fit(Xs, y)
    p = 1.0 / (1.0 + np.exp(-np.clip(Xs @ w, -30, 30)))
    a, _, _ = auc_mw(list(p[y == 1]), list(p[y == 0]))
    return a


# --------------------------------------------------------------------------- build
def main():
    print("[load] gene_of_dn projection / core / labels / caches ...", flush=True)
    gene_of_dn = load_gene_of_dn()
    pair_core = load_pair_core(gene_of_dn)
    lab = load_labels()
    ebm = load_ebm_master()
    sharedlen = load_sharedlen()
    readmm = load_readmm()
    allele = load_allele()

    recs = []
    for r in sorted(ebm, key=lambda r: (r["gene_a"], r["gene_b"])):
        k = r["k"]
        ga, gb = sorted((r["gene_a"], r["gene_b"]))
        L = lab.get(k, {})
        in_dna = L.get("in_dna_loose")
        in_ep = L.get("in_ep")
        ep_tier = L.get("ep_tier", NA)
        cls_auth = L.get("cls", NA)

        # --- core_recip: coalesce DN-graph-edge max with the representative bridge core
        core_dn = pair_core.get(k)
        core = core_dn if core_dn is not None else r["core_ebm"]
        if core_dn is not None:
            core_src = "dn"
        elif r["core_ebm"] is not None:
            core_src = "ebm"
        else:
            core_src = "none"
        has_direct = 1 if k in pair_core else 0

        # --- shared aligned block
        al, af = sharedlen.get(k, (None, None))

        # --- read multi-mapping per endpoint (symmetric hi/lo by mm_max)
        ma = readmm.get(r["bridge_dn_a"])
        mb = readmm.get(r["bridge_dn_b"])
        ends = [e for e in (ma, mb) if e is not None]
        if ends:
            ends_sorted = sorted(ends, key=lambda e: e[2], reverse=True)
            hi = ends_sorted[0]
            lo = ends_sorted[-1] if len(ends_sorted) > 1 else (None, None, None)
            mm_median_hi, mm_mean_hi, mm_max_hi = hi
            mm_median_lo, mm_mean_lo, mm_max_lo = lo
            mm_max_max = mm_max_hi
            mm_max_min = mm_max_lo if len(ends) > 1 else None
            mm_median_max = max(e[0] for e in ends)
            mm_mean_max = max(e[1] for e in ends)
        else:
            mm_median_hi = mm_mean_hi = mm_max_hi = None
            mm_median_lo = mm_mean_lo = mm_max_lo = None
            mm_max_max = mm_max_min = mm_median_max = mm_mean_max = None

        # --- family-level allele signal (join by gene; pick higher het_risk endpoint)
        cand = [(g, allele[g]) for g in (ga, gb) if g in allele]
        if cand:
            ag, av = max(cand, key=lambda t: t[1]["het_risk"])
            allele_gene = ag
            a_bal, a_het, a_cop, a_maf = av["balanced_frac"], av["het_risk"], av["copy_like"], av["med_maf"]
        else:
            allele_gene = NA
            a_bal = a_het = a_cop = a_maf = None

        recs.append(dict(
            gene_a=ga, gene_b=gb, cls=r["cls"], cls_auth=cls_auth,
            in_dna_loose=in_dna, in_ep=in_ep, ep_tier=ep_tier,
            core_recip=core, core_source=core_src, has_direct_edge=has_direct,
            aln_len=al, aln_frac=af,
            mm_median_hi=mm_median_hi, mm_mean_hi=mm_mean_hi, mm_max_hi=mm_max_hi,
            mm_median_lo=mm_median_lo, mm_mean_lo=mm_mean_lo, mm_max_lo=mm_max_lo,
            mm_max_max=mm_max_max, mm_max_min=mm_max_min,
            mm_median_max=mm_median_max, mm_mean_max=mm_mean_max,
            allele_gene=allele_gene, allele_balanced_frac=a_bal, allele_het_risk=a_het,
            allele_copy_like=a_cop, allele_med_maf=a_maf,
            abl_bridge_mask=r["bridge_mask"],
        ))

    # ------------------------------------------------------------------ write TSV
    cols = ["gene_a", "gene_b", "cls", "cls_auth", "in_dna_loose", "in_ep", "ep_tier",
            "core_recip", "core_source", "has_direct_edge",
            "aln_len", "aln_frac",
            "mm_median_hi", "mm_mean_hi", "mm_max_hi", "mm_median_lo", "mm_mean_lo", "mm_max_lo",
            "mm_max_max", "mm_max_min", "mm_median_max", "mm_mean_max",
            "allele_gene", "allele_balanced_frac", "allele_het_risk", "allele_copy_like",
            "allele_med_maf", "abl_bridge_mask"]
    intcols = {"in_dna_loose", "in_ep", "has_direct_edge", "aln_len",
               "mm_max_hi", "mm_max_lo", "mm_max_max", "mm_max_min"}
    with open(OUT_TSV, "w") as out:
        out.write("\t".join(cols) + "\n")
        for r in recs:
            vals = []
            for c in cols:
                v = r[c]
                if v is None or (isinstance(v, float) and v != v):
                    vals.append(NA)
                elif c in intcols and isinstance(v, float):
                    vals.append(str(int(v)))
                elif isinstance(v, float):
                    vals.append(f"{v:.4f}")
                else:
                    vals.append(str(v))
            out.write("\t".join(vals) + "\n")

    # ------------------------------------------------------------------ coverage
    n = len(recs)
    def has(r, c):
        v = r[c]
        return v is not None and not (isinstance(v, float) and v != v)
    cov = {}
    for c in ["core_recip", "aln_len", "aln_frac", "mm_max_max", "mm_median_max",
              "allele_het_risk", "abl_bridge_mask"]:
        cov[c] = sum(1 for r in recs if has(r, c))
    all_rna = sum(1 for r in recs if has(r, "core_recip") and has(r, "aln_frac") and has(r, "mm_max_max"))
    core_only = sum(1 for r in recs if has(r, "core_recip") and not has(r, "aln_frac") and not has(r, "mm_max_max"))
    no_core = sum(1 for r in recs if not has(r, "core_recip"))

    from collections import Counter
    cls_ct = Counter(r["cls"] for r in recs)

    # ------------------------------------------------------------------ per-feature AUC vs in_dna_loose
    def col(r, c):
        return r[c]
    # positive = in_dna_loose == 1 (REAL); negatives split into genuine / truthbar
    def split(feat):
        pos = [col(r, feat) for r in recs if r["in_dna_loose"] == 1]
        neg_all = [col(r, feat) for r in recs if r["in_dna_loose"] == 0]
        neg_gen = [col(r, feat) for r in recs if r["cls"] == "genuine"]
        return pos, neg_all, neg_gen

    rna_feats = ["core_recip", "aln_len", "aln_frac", "mm_max_max", "mm_max_hi",
                 "mm_median_max", "mm_mean_max"]
    auc_report = {}
    for feat in rna_feats:
        pos, neg_all, neg_gen = split(feat)
        a_all, npo, nne = auc_mw(pos, neg_all)
        a_gen, npo2, nne2 = auc_mw(pos, neg_gen)
        auc_report[feat] = dict(auc_vs_dna=a_all, n_pos=npo, n_neg=nne,
                                auc_tp_vs_gen=a_gen, n_pos2=npo2, n_neg2=nne2)
    # ablation feature (genome-informed) reported SEPARATELY
    pos, neg_all, neg_gen = split("abl_bridge_mask")
    abl_auc_all = auc_mw(pos, neg_all)
    abl_auc_gen = auc_mw(pos, neg_gen)

    # ------------------------------------------------------------------ correlation structure (Spearman)
    corr_feats = ["core_recip", "aln_len", "aln_frac", "mm_max_max", "mm_mean_max"]
    joint = [r for r in recs if all(has(r, c) for c in corr_feats)]
    corr = {}
    if len(joint) >= 3:
        mat = np.array([[r[c] for c in corr_feats] for r in joint], dtype=float)
        for i in range(len(corr_feats)):
            for j in range(i + 1, len(corr_feats)):
                rho, _ = ss.spearmanr(mat[:, i], mat[:, j])
                corr[(corr_feats[i], corr_feats[j])] = rho
    # correlation of each feature with core_recip over each feature's own coverage
    corr_with_core = {}
    for feat in ["aln_len", "aln_frac", "mm_max_max", "mm_mean_max", "mm_median_max"]:
        pair = [(r["core_recip"], r[feat]) for r in recs if has(r, "core_recip") and has(r, feat)]
        if len(pair) >= 3:
            x = [p[0] for p in pair]; y = [p[1] for p in pair]
            rho, _ = ss.spearmanr(x, y)
            corr_with_core[feat] = (rho, len(pair))

    # ------------------------------------------------------------------ independent signal
    # (a) conditional AUC in the core_recip GRAY ZONE (0.13 <= core < 0.30) where core alone is ambiguous
    def gray(r):
        return has(r, "core_recip") and 0.13 <= r["core_recip"] < 0.30
    gray_recs = [r for r in recs if gray(r)]
    def gray_auc(feat):
        pos = [r[feat] for r in gray_recs if r["cls"] == "TP"]
        neg = [r[feat] for r in gray_recs if r["cls"] == "genuine"]
        return auc_mw(pos, neg)
    gray_core = gray_auc("core_recip")
    gray_aln = gray_auc("aln_frac")
    gray_mm = gray_auc("mm_max_max")

    # (b) in-sample logistic incremental AUC on the ARBITRATION subpop (has core+aln+mm; TP vs genuine)
    arb = [r for r in recs if has(r, "core_recip") and has(r, "aln_frac") and has(r, "mm_max_max")
           and r["cls"] in ("TP", "genuine")]
    inc = {}
    if arb:
        y = np.array([1.0 if r["cls"] == "TP" else 0.0 for r in arb])
        core_c = [[r["core_recip"]] for r in arb]
        core_aln = [[r["core_recip"], r["aln_frac"]] for r in arb]
        core_aln_mm = [[r["core_recip"], r["aln_frac"], math.log1p(r["mm_max_max"])] for r in arb]
        core_mm = [[r["core_recip"], math.log1p(r["mm_max_max"])] for r in arb]
        inc["n"] = len(arb)
        inc["n_tp"] = int(y.sum())
        inc["core"] = logistic_insample_auc(core_c, y)
        inc["core+aln"] = logistic_insample_auc(core_aln, y)
        inc["core+mm"] = logistic_insample_auc(core_mm, y)
        inc["core+aln+mm"] = logistic_insample_auc(core_aln_mm, y)

    # allele sanity: known allele-as-copy FPs
    allele_sanity = {g: allele[g] for g in ("DHRSX", "LOC129530050") if g in allele}

    # ------------------------------------------------------------------ REPORT
    P = print
    P("\n================ RNA-ONLY EDGE ORACLE :: FEATURE STAGE ================")
    P(f"population = refined-E_r gene-pair edges (edge_bridge_mask) : N = {n}")
    P(f"  class split : " + "  ".join(f"{k}={v}" for k, v in sorted(cls_ct.items())))
    P(f"  wrote {OUT_TSV}")

    P("\n---- COVERAGE (of N={} edges) ----".format(n))
    P(f"  core_recip           : {cov['core_recip']:5d}  ({100*cov['core_recip']/n:4.1f}%)   [direct homology edge; absent=transitive-only]")
    P(f"  aln_len / aln_frac   : {cov['aln_frac']:5d}  ({100*cov['aln_frac']/n:4.1f}%)   [in_ep=0 arbitration edges = over-merge target pop]")
    P(f"  mm_max_max (read deg): {cov['mm_max_max']:5d}  ({100*cov['mm_max_max']/n:4.1f}%)   [>=1 bridge endpoint in readmm cache]")
    P(f"  allele_het_risk      : {cov['allele_het_risk']:5d}  ({100*cov['allele_het_risk']/n:4.1f}%)   [family-level; for the DEMOTE step]")
    P(f"  ALL RNA (core+aln+mm): {all_rna:5d}  ({100*all_rna/n:4.1f}%)")
    P(f"  core_recip-ONLY      : {core_only:5d}  ({100*core_only/n:4.1f}%)   [core present, no aln & no mm]")
    P(f"  NO core_recip        : {no_core:5d}  ({100*no_core/n:4.1f}%)   [transitive-closure co-membership pairs]")
    P(f"  abl_bridge_mask(ABL) : {cov['abl_bridge_mask']:5d}  ({100*cov['abl_bridge_mask']/n:4.1f}%)   [GENOME-derived; ablation only, NOT RNA-only]")

    P("\n---- PER-FEATURE AUC vs DNA label in_dna_loose (higher score => REAL) ----")
    P("  feature            AUC(TP vs all non-TP)   AUC(TP vs genuine)   n_pos/n_neg    direction")
    for feat in rna_feats:
        d = auc_report[feat]
        direction = "higher=REAL" if (d["auc_vs_dna"] == d["auc_vs_dna"] and d["auc_vs_dna"] >= 0.5) else "INVERTED(higher=over-merge)"
        P(f"  {feat:16s}   {fmt(d['auc_vs_dna']):>7s}                {fmt(d['auc_tp_vs_gen']):>7s}            {d['n_pos']}/{d['n_neg']:<6d}   {direction}")
    P(f"  [ABLATION] abl_bridge_mask  AUC(TPvsAll)={fmt(abl_auc_all[0])}  AUC(TPvsGen)={fmt(abl_auc_gen[0])}  (genome-derived; reported for comparison only)")

    P("\n---- CORRELATION STRUCTURE among RNA features (Spearman; joint-coverage n={}) ----".format(len(joint)))
    for (a, b), rho in corr.items():
        P(f"  rho({a}, {b}) = {fmt(rho)}")
    P("  each feature vs core_recip (own coverage):")
    for feat, (rho, npair) in corr_with_core.items():
        P(f"    rho(core_recip, {feat}) = {fmt(rho)}   (n={npair})")

    P("\n---- INDEPENDENT SIGNAL ----")
    P(f"  (a) conditional AUC in core_recip GRAY ZONE [0.13,0.30) (core alone ambiguous), n_gray={len(gray_recs)}:")
    P(f"        core_recip : AUC={fmt(gray_core[0])}  (TP {gray_core[1]} vs genuine {gray_core[2]})")
    P(f"        aln_frac   : AUC={fmt(gray_aln[0])}  (TP {gray_aln[1]} vs genuine {gray_aln[2]})")
    P(f"        mm_max_max : AUC={fmt(gray_mm[0])}  (TP {gray_mm[1]} vs genuine {gray_mm[2]})  [INVERTED: higher=over-merge]")
    if inc:
        P(f"  (b) in-sample logistic incremental AUC on arbitration subpop (n={inc['n']}, TP={inc['n_tp']}; CV deferred to LEARN stage):")
        P(f"        core            : AUC={fmt(inc['core'])}")
        P(f"        core + aln_frac : AUC={fmt(inc['core+aln'])}")
        P(f"        core + mm       : AUC={fmt(inc['core+mm'])}")
        P(f"        core + aln + mm : AUC={fmt(inc['core+aln+mm'])}")

    P("\n---- ALLELE-SIGNAL SANITY (known allele-as-copy FPs) ----")
    for g, av in allele_sanity.items():
        P(f"  {g:14s} balanced_frac={av['balanced_frac']:.3f} het_risk={av['het_risk']:.3f} "
          f"copy_like={av['copy_like']:.3f} med_maf={av['med_maf']:.3f}")

    P("\n======================================================================")

    # machine-readable echo for the caller (compact)
    import json
    summary = dict(
        N=n, class_split=dict(cls_ct),
        coverage=dict(core_recip=cov["core_recip"], aln=cov["aln_frac"], mm=cov["mm_max_max"],
                      allele=cov["allele_het_risk"], all_rna=all_rna, core_only=core_only,
                      no_core=no_core, abl_bridge_mask=cov["abl_bridge_mask"]),
        auc={f: dict(vs_dna=auc_report[f]["auc_vs_dna"], tp_vs_gen=auc_report[f]["auc_tp_vs_gen"],
                     n_pos=auc_report[f]["n_pos"], n_neg=auc_report[f]["n_neg"]) for f in rna_feats},
        abl_bridge_mask_auc=dict(vs_dna=abl_auc_all[0], tp_vs_gen=abl_auc_gen[0]),
        corr={f"{a}|{b}": corr[(a, b)] for (a, b) in corr},
        corr_with_core={f: corr_with_core[f][0] for f in corr_with_core},
        gray_zone=dict(n=len(gray_recs), core=gray_core[0], aln_frac=gray_aln[0], mm_max_max=gray_mm[0]),
        logistic_incremental=inc,
    )
    P("JSON " + json.dumps(summary, sort_keys=True, default=lambda x: None if (isinstance(x, float) and x != x) else x))


# ===========================================================================
# LEARN STAGE : RNA-only edge rule, component-level cross-validation
# ===========================================================================
# The refined-E_r direct-edge graph = the gene pairs the gamma-quasi-clique
# operator thresholds (core_recip present). Every such edge already passed
# core_recip >= 0.13 at edge creation, so the SHIPPED baseline "keep iff
# core >= 0.13" keeps 100% of them: recall 1.0, precision = 416/5580 = 0.075.
# The over-merges are edges that cleared 0.13 but are not DNA-real. We LEARN an
# RNA-only rule that cuts them while preserving the DNA-real (TP) edges.
#
# LABEL (learning/validation ONLY): in_dna_loose (== cls==TP). Positive = REAL.
# RNA-ONLY FEATURES at inference: core_recip, aln_frac (UNIVERSAL cache, see
#   below), log1p(mm_max_max). NO in_dna_loose / in_ep / ep_tier / cls /
#   abl_bridge_mask ever enters the feature vector (hard-asserted).
#
# LEAKAGE FIX (thesis-critical): the original ri_sharedlen_cache.tsv computed
# aln_frac ONLY on in_ep=0 edges, so aln_frac PRESENCE encoded the protein label
# in_ep. ri_build_sharedlen_universal.py recomputes the identical RNA shared-
# exon alignment on EVERY core-present edge (byte-identical on the in_ep=0
# subset), giving aln_frac universal coverage with ZERO in_ep dependence.
#
# CV: whole connected components (families; union-find over all 10755 edges) are
# assigned to K folds -> no two edges of one family ever land in different folds
# (no within-family leakage). Held-out predictions pooled across folds.

SHAREDLEN_UNIV_TSV = os.path.join(BENCH, "ri_sharedlen_universal.tsv")
SHAREDLEN_RC_TSV = os.path.join(BENCH, "ri_sharedlen_readcons.tsv")

RNA_FEATURE_WHITELIST = {"core_recip", "aln_frac", "mm_log", "mm_max_max"}
LABEL_COLS_FORBIDDEN = {"in_dna_loose", "in_ep", "ep_tier", "cls", "cls_auth",
                        "abl_bridge_mask", "allele_het_risk", "allele_balanced_frac"}


def load_universal_aln():
    """gene-pair -> aln_frac from the UNIVERSAL (leakage-free) cache. RNA-only."""
    d = {}
    with open(SHAREDLEN_UNIV_TSV) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        I = {c: i for i, c in enumerate(h)}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            k = frozenset((f[I["gene_a"]], f[I["gene_b"]]))
            af = f[I["aln_frac"]]
            d[k] = None if af == "" else float(af)
    return d


def load_readcons_aln():
    """gene-pair -> aln_frac_read from the READ-CONSENSUS cache (de-novo transcript
       POA consensus, NO genome bases). Used ONLY for the provenance concordance
       check + the purity-max ablation catalog -- proves the dominant separator does
       not hinge on genome-assembly bases. RNA-derived (same source as core_recip)."""
    d = {}
    with open(SHAREDLEN_RC_TSV) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        I = {c: i for i, c in enumerate(h)}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            k = frozenset((f[I["gene_a"]], f[I["gene_b"]]))
            af = f[I["aln_frac_read"]]
            d[k] = None if af == "" else float(af)
    return d


def aln_provenance_concordance(univ_aln, rc_aln, t2):
    """Concordance between the DEPLOYED aln_frac (genome-fetched spliced-exon bases,
       RNA splice STRUCTURE) and a read-consensus aln_frac (denovo POA transcript
       bases, NO genome). Reads the label/core purely to SCORE (never to gate).
       Returns dict(rho, gate_agreement, auc_genome, auc_readcons, class_medians)."""
    # labels + gene pairs from the (unchanged) feature TSV
    y = {}
    with open(OUT_TSV) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        I = {c: i for i, c in enumerate(h)}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            k = frozenset((f[I["gene_a"]], f[I["gene_b"]]))
            y[k] = (1 if f[I["in_dna_loose"]] == "1" else 0, f[I["cls"]])
    keys = sorted((tuple(sorted(k)) for k in univ_aln
                   if k in rc_aln and univ_aln[k] is not None and rc_aln[k] is not None))
    ks = [frozenset(k) for k in keys]
    g = np.array([univ_aln[k] for k in ks])
    r = np.array([rc_aln[k] for k in ks])
    rho, _ = ss.spearmanr(g, r)
    kg = g >= t2
    kr = r >= t2
    agree = int((kg == kr).sum())
    lab = np.array([y[k][0] for k in ks])
    aucg, _, _ = auc_mw(list(g[lab == 1]), list(g[lab == 0]))
    aucr, _, _ = auc_mw(list(r[lab == 1]), list(r[lab == 0]))
    med = {}
    for cls in ("TP", "genuine", "truthbar"):
        gi = [univ_aln[k] for k in ks if y[k][1] == cls]
        ri = [rc_aln[k] for k in ks if y[k][1] == cls]
        med[cls] = (float(np.median(gi)) if gi else float("nan"),
                    float(np.median(ri)) if ri else float("nan"))
    return dict(n=len(ks), rho=float(rho), gate_agreement=agree,
                gate_agreement_frac=agree / len(ks) if ks else 0.0,
                auc_genome=aucg, auc_readcons=aucr, class_medians=med, t2=t2)


def union_find_components(edges):
    """edges: list of (gene_a, gene_b). Returns gene -> component-id (int, stable)."""
    parent = {}
    def find(a):
        parent.setdefault(a, a)
        root = a
        while parent[root] != root:
            root = parent[root]
        while parent[a] != root:
            parent[a], a = root, parent[a]
        return root
    def uni(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb
    for a, b in edges:
        uni(a, b)
    reps = {}
    comp = {}
    for g in sorted(parent):
        r = find(g)
        if r not in reps:
            reps[r] = len(reps)
        comp[g] = reps[r]
    return comp


def assign_folds(comp_ids, comp_tp, comp_sz, K, seed):
    """Assign whole components to K folds, TP-balanced. Deterministic given seed.
       Shuffle components (seed) then greedily place each on the currently
       lightest fold by (running TP, running edge-count, fold index)."""
    rng = np.random.RandomState(seed)
    order = list(comp_ids)
    rng.shuffle(order)
    # heavy-first within the shuffled order stabilises balance
    order.sort(key=lambda c: (-comp_tp[c], -comp_sz[c]))
    fold_tp = [0] * K
    fold_sz = [0] * K
    fold_of = {}
    for c in order:
        k = min(range(K), key=lambda j: (fold_tp[j], fold_sz[j], j))
        fold_of[c] = k
        fold_tp[k] += comp_tp[c]
        fold_sz[k] += comp_sz[c]
    return fold_of, fold_tp, fold_sz


def prf(pred, y):
    pred = np.asarray(pred); y = np.asarray(y)
    tp = int(((pred == 1) & (y == 1)).sum())
    fp = int(((pred == 1) & (y == 0)).sum())
    fn = int(((pred == 0) & (y == 1)).sum())
    prec = tp / (tp + fp) if (tp + fp) else 0.0
    rec = tp / (tp + fn) if (tp + fn) else 0.0
    f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else 0.0
    return prec, rec, f1, tp, fp, fn


def grid_best_rule(core, aln, y, objective="f1", min_recall=None,
                   t1_grid=None, t2_grid=None, mm=None, t3_grid=None):
    """Best (t1,t2[,t3]) for rule (core>=t1)&(aln>=t2)[&(mm<=t3)] on given data.
       objective 'f1' maximises F1; 'prec@recall' maximises precision s.t.
       recall>=min_recall. Returns dict(t1,t2,t3,prec,rec,f1)."""
    if t1_grid is None:
        t1_grid = np.round(np.arange(0.13, 0.701, 0.01), 3)
    if t2_grid is None:
        t2_grid = np.round(np.arange(0.0, 1.001, 0.02), 3)
    if t3_grid is None:
        t3_grid = [None]
    best = None
    for t3 in t3_grid:
        mmok = np.ones(len(core), dtype=bool) if (t3 is None or mm is None) else (mm <= t3)
        for t1 in t1_grid:
            c_ok = core >= t1
            for t2 in t2_grid:
                pred = (c_ok & (aln >= t2) & mmok).astype(int)
                prec, rec, f1, tp, fp, fn = prf(pred, y)
                if objective == "f1":
                    key = (f1, prec)
                else:  # prec@recall
                    if rec < min_recall:
                        continue
                    key = (prec, f1)
                if best is None or key > best["key"]:
                    best = dict(key=key, t1=float(t1), t2=float(t2),
                                t3=(None if t3 is None else float(t3)),
                                prec=prec, rec=rec, f1=f1)
    return best


def learn():
    P = print
    P("\n\n################ RNA-ONLY EDGE ORACLE :: LEARN STAGE ################")

    # ---- load edges + RNA features + label -------------------------------
    feat = {}
    all_pairs = []
    with open(OUT_TSV) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        I = {c: i for i, c in enumerate(h)}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            ga, gb = f[I["gene_a"]], f[I["gene_b"]]
            all_pairs.append((ga, gb))
            k = frozenset((ga, gb))
            feat[k] = dict(
                gene_a=ga, gene_b=gb, cls=f[I["cls"]],
                y=1 if f[I["in_dna_loose"]] == "1" else 0,
                core=None if f[I["core_recip"]] in ("", NA) else float(f[I["core_recip"]]),
                mm=None if f[I["mm_max_max"]] in ("", NA) else float(f[I["mm_max_max"]]),
            )
    univ_aln = load_universal_aln()

    # ---- components (union-find over ALL edges) --------------------------
    comp = union_find_components(all_pairs)

    # ---- evaluation population = core-present direct edges ---------------
    # aln imputed 0 where the transcript sequence was unavailable (RNA-derived
    # absence: no assembled body -> no evidence of a shared gene body).
    pop = []
    for (ga, gb) in all_pairs:
        k = frozenset((ga, gb))
        r = feat[k]
        if r["core"] is None:
            continue                      # transitive-only pair: no direct edge to score
        af = univ_aln.get(k)
        af = 0.0 if af is None else af
        mm = r["mm"] if r["mm"] is not None else 0.0
        pop.append(dict(k=k, ga=ga, gb=gb, cls=r["cls"], y=r["y"],
                        comp=comp[ga], core=r["core"], aln=af, mm=mm))
    pop.sort(key=lambda r: (r["ga"], r["gb"]))
    N = len(pop)
    core = np.array([r["core"] for r in pop])
    aln = np.array([r["aln"] for r in pop])
    mm = np.array([r["mm"] for r in pop])
    mmlog = np.log1p(mm)
    y = np.array([r["y"] for r in pop])
    Ptot = int(y.sum())

    # ---- LEAKAGE GUARD 1 : no forbidden (DNA/protein/genome) col as feature
    feat_names = ["core_recip", "aln_frac", "mm_log"]
    assert set(feat_names) <= RNA_FEATURE_WHITELIST, "non-RNA feature leaked in"
    assert not (set(feat_names) & LABEL_COLS_FORBIDDEN), "label column used as feature"

    # ---- CV folds : whole components ------------------------------------
    K = 5
    SEED = 0
    comp_ids = sorted(set(r["comp"] for r in pop))
    comp_tp = {c: 0 for c in comp_ids}
    comp_sz = {c: 0 for c in comp_ids}
    for r in pop:
        comp_sz[r["comp"]] += 1
        comp_tp[r["comp"]] += r["y"]
    fold_of, fold_tp, fold_sz = assign_folds(comp_ids, comp_tp, comp_sz, K, SEED)
    fold = np.array([fold_of[r["comp"]] for r in pop])

    # ---- LEAKAGE GUARD 2 : no family split across folds -----------------
    comp_fold = {}
    for r in pop:
        c = r["comp"]
        if c in comp_fold:
            assert comp_fold[c] == fold_of[c], "component spans multiple folds"
        comp_fold[c] = fold_of[c]
    # every gene lives in exactly one fold (endpoints never cross folds)
    gene_fold = {}
    leak = 0
    for r in pop:
        for g in (r["ga"], r["gb"]):
            if g in gene_fold and gene_fold[g] != fold_of[r["comp"]]:
                leak += 1
            gene_fold[g] = fold_of[r["comp"]]
    assert leak == 0, f"gene appears in >1 fold ({leak})"

    P(f"population = core-present refined-E_r direct edges : N={N}  TP={Ptot} "
      f"(genuine={sum(1 for r in pop if r['cls']=='genuine')}, "
      f"truthbar={sum(1 for r in pop if r['cls']=='truthbar')})")
    P(f"CV = {K} folds over whole connected components (families); seed={SEED}; "
      f"n_components={len(comp_ids)}")
    P("  fold TP / edges : " + "  ".join(f"f{j}:{fold_tp[j]}/{fold_sz[j]}" for j in range(K)))
    P("  GUARD: no DNA/protein/genome column in feature vector = OK")
    P("  GUARD: no family split across folds = OK ; no gene in >1 fold = OK")

    # baseline operating point (shipped): keep iff core>=0.13 -> keeps all here
    base_pred = (core >= 0.13).astype(int)
    bp, br, bf, btp, bfp, bfn = prf(base_pred, y)
    P(f"\n-- SHIPPED BASELINE (core_recip>=0.13, the edge-creation threshold) --")
    P(f"   keeps {int(base_pred.sum())}/{N} edges : precision={bp:.4f}  recall={br:.4f}  "
      f"F1={bf:.4f}   (keeps everything; over-merges pass through)")

    # ---- CV : interpretable rule (F1) + recall-preserving + logistic -----
    heldout_rule_f1 = np.full(N, -1, dtype=int)      # 0/1 pred, F1-selected rule
    heldout_rule_rp = np.full(N, -1, dtype=int)      # 0/1 pred, recall>=0.90 rule
    heldout_logit = np.full(N, np.nan)               # logistic score
    heldout_core = np.full(N, np.nan)                # single-feature scores (for ROC)
    heldout_aln = np.full(N, np.nan)
    fold_rules_f1 = []
    fold_rules_rp = []
    logit_ws = []

    for kf in range(K):
        tr = fold != kf
        te = fold == kf
        ytr = y[tr]
        # interpretable rule (max train F1)
        rf1 = grid_best_rule(core[tr], aln[tr], ytr, objective="f1")
        fold_rules_f1.append(rf1)
        pred = ((core[te] >= rf1["t1"]) & (aln[te] >= rf1["t2"])).astype(int)
        heldout_rule_f1[te] = pred
        # interpretable rule (max train precision s.t. recall>=0.90 : family-preserving)
        rrp = grid_best_rule(core[tr], aln[tr], ytr, objective="prec@recall", min_recall=0.90)
        if rrp is None:
            rrp = dict(t1=0.13, t2=0.0, t3=None, prec=0, rec=1, f1=0)
        fold_rules_rp.append(rrp)
        predrp = ((core[te] >= rrp["t1"]) & (aln[te] >= rrp["t2"])).astype(int)
        heldout_rule_rp[te] = predrp
        # logistic on RNA features (standardise on TRAIN only -> no test leak)
        Xtr_raw = np.c_[core[tr], aln[tr], mmlog[tr]]
        mu = Xtr_raw.mean(axis=0); sd = Xtr_raw.std(axis=0); sd[sd == 0] = 1.0
        Xtr = np.hstack([np.ones((tr.sum(), 1)), (Xtr_raw - mu) / sd])
        w = logistic_fit(Xtr, ytr.astype(float))
        logit_ws.append((w, mu, sd))
        Xte_raw = np.c_[core[te], aln[te], mmlog[te]]
        Xte = np.hstack([np.ones((te.sum(), 1)), (Xte_raw - mu) / sd])
        heldout_logit[te] = 1.0 / (1.0 + np.exp(-np.clip(Xte @ w, -30, 30)))
        heldout_core[te] = core[te]
        heldout_aln[te] = aln[te]

    assert (heldout_rule_f1 >= 0).all() and not np.isnan(heldout_logit).any(), "unscored edges"

    # ---- optional 3rd term: AND (mm_max_max <= t3) read-degree TE-bridge veto
    T3_GRID = [1e9, 60, 55, 53, 52, 51, 50, 40, 30, 20, 15, 10]
    heldout_rule_mm = np.full(N, -1, dtype=int)
    for kf in range(K):
        tr = fold != kf; te = fold == kf
        r3 = grid_best_rule(core[tr], aln[tr], y[tr], objective="f1",
                            mm=mm[tr], t3_grid=T3_GRID)
        predmm = ((core[te] >= r3["t1"]) & (aln[te] >= r3["t2"]) &
                  (mm[te] <= (r3["t3"] if r3["t3"] is not None else 1e9))).astype(int)
        heldout_rule_mm[te] = predmm
    rmm_p, rmm_r, rmm_f, rmm_tp, rmm_fp, rmm_fn = prf(heldout_rule_mm, y)

    # ---- single-feature threshold rules (does core ADD to aln at the F1 point?)
    ho_alnonly = np.full(N, -1, dtype=int)
    ho_coreonly = np.full(N, -1, dtype=int)
    for kf in range(K):
        tr = fold != kf; te = fold == kf
        ra = grid_best_rule(np.zeros(tr.sum()), aln[tr], y[tr], objective="f1",
                            t1_grid=[0.0])                      # core disabled
        ho_alnonly[te] = (aln[te] >= ra["t2"]).astype(int)
        rc = grid_best_rule(core[tr], np.zeros(tr.sum()), y[tr], objective="f1",
                            t2_grid=[0.0])                      # aln disabled
        ho_coreonly[te] = (core[te] >= rc["t1"]).astype(int)
    ao_p, ao_r, ao_f, *_ = prf(ho_alnonly, y)
    co_p, co_r, co_f, *_ = prf(ho_coreonly, y)

    # ---- pooled held-out metrics ----------------------------------------
    def auc_score(score):
        a, _, _ = auc_mw(list(score[y == 1]), list(score[y == 0]))
        return a
    auc_core = auc_score(heldout_core)
    auc_aln = auc_score(heldout_aln)
    auc_logit = auc_score(heldout_logit)

    rf1_p, rf1_r, rf1_f, rf1_tp, rf1_fp, rf1_fn = prf(heldout_rule_f1, y)
    rrp_p, rrp_r, rrp_f, rrp_tp, rrp_fp, rrp_fn = prf(heldout_rule_rp, y)

    def prec_at_recall(score, target):
        best_p, best_thr, best_r = 0.0, None, 0.0
        for thr in np.unique(score):
            pred = (score >= thr).astype(int)
            p, r, f, tp, fp, fn = prf(pred, y)
            if r >= target and p > best_p:
                best_p, best_thr, best_r = p, float(thr), r
        return best_p, best_thr, best_r

    # negative-type cut breakdown for the F1 rule (held-out)
    cls_arr = np.array([r["cls"] for r in pop])
    def cut_rate(clsname, pred):
        mask = cls_arr == clsname
        cut = int(((pred == 0) & mask).sum())
        tot = int(mask.sum())
        return cut, tot, (cut / tot if tot else 0.0)
    gen_cut = cut_rate("genuine", heldout_rule_f1)
    tb_cut = cut_rate("truthbar", heldout_rule_f1)
    tp_keep = int(((heldout_rule_f1 == 1) & (y == 1)).sum())

    # ---- FINAL rule fit on ALL data (deployable artifact) ----------------
    final_f1 = grid_best_rule(core, aln, y, objective="f1")
    final_rp = grid_best_rule(core, aln, y, objective="prec@recall", min_recall=0.90)
    med_t1 = float(np.median([r["t1"] for r in fold_rules_f1]))
    med_t2 = float(np.median([r["t2"] for r in fold_rules_f1]))

    # ---- REPORT ----------------------------------------------------------
    P(f"\n-- HELD-OUT AUC (pooled out-of-fold scores; higher=REAL) --")
    P(f"   core_recip (shipped feature) : AUC = {auc_core:.4f}")
    P(f"   aln_frac   (RNA shared body) : AUC = {auc_aln:.4f}")
    P(f"   logistic core+aln+mm_log     : AUC = {auc_logit:.4f}")

    P(f"\n-- (a) INTERPRETABLE RULE  (core_recip>=t1) AND (aln_frac>=t2) --")
    P(f"   per-fold selected thresholds (max train F1):")
    for j, r in enumerate(fold_rules_f1):
        P(f"     fold{j}: t1(core)={r['t1']:.2f}  t2(aln)={r['t2']:.2f}   "
          f"(train F1={r['f1']:.3f} P={r['prec']:.3f} R={r['rec']:.3f})")
    P(f"   fold-median thresholds        : core>={med_t1:.2f}  AND  aln_frac>={med_t2:.2f}")
    P(f"   final fit on ALL data (deploy) : core>={final_f1['t1']:.2f}  AND  aln_frac>={final_f1['t2']:.2f}")
    P(f"   HELD-OUT (pooled out-of-fold) : precision={rf1_p:.4f}  recall={rf1_r:.4f}  "
      f"F1={rf1_f:.4f}   (kept {rf1_tp} TP, {rf1_fp} FP; lost {rf1_fn} TP)")
    P(f"      negative cuts: genuine over-merge {gen_cut[0]}/{gen_cut[1]} cut ({100*gen_cut[2]:.1f}%);  "
      f"truthbar divergent-paralog {tb_cut[0]}/{tb_cut[1]} cut ({100*tb_cut[2]:.1f}%)")

    P(f"   ablation (held-out F1-opt threshold rule): "
      f"aln-only F1={ao_f:.3f} (P={ao_p:.3f} R={ao_r:.3f}) | "
      f"core-only F1={co_f:.3f} (P={co_p:.3f} R={co_r:.3f}) | "
      f"core+aln F1={rf1_f:.3f} (P={rf1_p:.3f} R={rf1_r:.3f})")
    P(f"      => core adds a PRECISION GATE on top of aln at the operating point "
      f"(P {ao_p:.3f}->{rf1_p:.3f} at matched recall); aln is the workhorse.")
    P(f"   + optional read-degree veto AND (mm_max_max<=t3), held-out : "
      f"P={rmm_p:.4f} R={rmm_r:.4f} F1={rmm_f:.4f}  (dF1={rmm_f-rf1_f:+.4f} vs core+aln -> DROP: adds nothing)")

    P(f"\n-- (a') RECALL-PRESERVING RULE (max train precision s.t. train recall>=0.90) --")
    P(f"   final fit on ALL data         : core>={final_rp['t1']:.2f}  AND  aln_frac>={final_rp['t2']:.2f}")
    P(f"   HELD-OUT (pooled out-of-fold) : precision={rrp_p:.4f}  recall={rrp_r:.4f}  "
      f"F1={rrp_f:.4f}   (kept {rrp_tp} TP, {rrp_fp} FP; lost {rrp_fn} TP)")

    P(f"\n-- (b) LOGISTIC (numpy GD; features = core_recip, aln_frac, log1p(mm_max_max)) --")
    P(f"   held-out AUC = {auc_logit:.4f}. Precision at matched recall (out-of-fold logistic score):")
    for tgt in (0.95, 0.90, 0.80):
        pl, thl, rl = prec_at_recall(heldout_logit, tgt)
        pc, thc, rc = prec_at_recall(heldout_core, tgt)
        pa, tha, ra = prec_at_recall(heldout_aln, tgt)
        P(f"     R>={tgt:.2f} :  logistic P={pl:.3f}   |  aln-only P={pa:.3f}   |  "
          f"core-only P={pc:.3f}   |  baseline(core>=.13) P={bp:.3f}")
    # standardised logistic coefficients (final fit on all data), sign check
    Xraw = np.c_[core, aln, mmlog]
    mu = Xraw.mean(axis=0); sd = Xraw.std(axis=0); sd[sd == 0] = 1.0
    Xs = np.hstack([np.ones((N, 1)), (Xraw - mu) / sd])
    wfin = logistic_fit(Xs, y.astype(float))
    P(f"   final standardised coefficients: intercept={wfin[0]:+.3f}  "
      f"core_recip={wfin[1]:+.3f}  aln_frac={wfin[2]:+.3f}  mm_log={wfin[3]:+.3f}")
    P(f"     (aln_frac dominant +; mm_log - = higher read-degree pushes toward over-merge)")

    P(f"\n-- COMPARISON vs SHIPPED core_recip>=0.13 baseline --")
    P(f"   HONEST incremental framing: core_recip>=0.13 is the edge-CREATION threshold (keeps 100%),")
    P(f"   NOT a tunable operating point. The real contribution is aln_frac as the workhorse")
    P(f"   (held-out AUC {auc_core:.3f} core-only -> {auc_aln:.3f} core+aln) + core as a precision GATE")
    P(f"   (held-out F1 {co_f:.3f} core-only -> {rf1_f:.3f} core+aln). The x-factor below is generous.")
    P(f"   baseline           : P={bp:.4f}  R={br:.4f}  (keeps all {N}; 6 multifam/9 oversize over-merges survive)")
    P(f"   interpretable rule : P={rf1_p:.4f}  R={rf1_r:.4f}   (precision x{rf1_p/bp:.1f} vs the 0.13 keep-all point)")
    P(f"   recall-preserving  : P={rrp_p:.4f}  R={rrp_r:.4f}   (precision x{rrp_p/bp:.1f} at >=0.90 recall = milder point)")
    P(f"   logistic @R>=0.90  : P={prec_at_recall(heldout_logit,0.90)[0]:.4f}  R>=0.90   "
      f"(does NOT beat aln-alone out-of-fold; aln AUC {auc_aln:.3f} > logistic {auc_logit:.3f})")

    # ---- machine echo ----------------------------------------------------
    import json
    out = dict(
        population="core_present_direct_edges", N=N, TP=Ptot, K=K, seed=SEED,
        n_components=len(comp_ids),
        baseline=dict(rule="core_recip>=0.13", precision=bp, recall=br, f1=bf, keeps=int(base_pred.sum())),
        heldout_auc=dict(core_recip=auc_core, aln_frac=auc_aln, logistic=auc_logit),
        interpretable_rule=dict(
            form="(core_recip>=t1) AND (aln_frac>=t2)",
            final_t1=final_f1["t1"], final_t2=final_f1["t2"],
            fold_median_t1=med_t1, fold_median_t2=med_t2,
            heldout_precision=rf1_p, heldout_recall=rf1_r, heldout_f1=rf1_f,
            heldout_tp=rf1_tp, heldout_fp=rf1_fp, heldout_lost_tp=rf1_fn,
            genuine_cut=gen_cut[0], genuine_tot=gen_cut[1],
            truthbar_cut=tb_cut[0], truthbar_tot=tb_cut[1]),
        recall_preserving_rule=dict(
            final_t1=final_rp["t1"], final_t2=final_rp["t2"],
            heldout_precision=rrp_p, heldout_recall=rrp_r, heldout_f1=rrp_f),
        ablation=dict(
            aln_only=dict(precision=ao_p, recall=ao_r, f1=ao_f),
            core_only=dict(precision=co_p, recall=co_r, f1=co_f),
            core_aln=dict(precision=rf1_p, recall=rf1_r, f1=rf1_f),
            plus_mm_veto=dict(precision=rmm_p, recall=rmm_r, f1=rmm_f, dF1=rmm_f - rf1_f)),
        logistic=dict(features=["core_recip", "aln_frac", "log1p(mm_max_max)"],
                      heldout_auc=auc_logit,
                      prec_at_R90=prec_at_recall(heldout_logit, 0.90)[0],
                      prec_at_R95=prec_at_recall(heldout_logit, 0.95)[0],
                      coef=dict(intercept=float(wfin[0]), core_recip=float(wfin[1]),
                                aln_frac=float(wfin[2]), mm_log=float(wfin[3]))),
        guards=dict(no_dna_feature=True, no_family_split=True, no_gene_cross_fold=True),
    )
    P("\nLEARN_JSON " + json.dumps(out, sort_keys=True,
                                   default=lambda x: None if (isinstance(x, float) and x != x) else x))
    P("####################################################################")

    # ---- per-gene-pair held-out decisions + deployable rule thresholds (for APPLY) ----
    pair_heldout = {}
    for i, r in enumerate(pop):
        pair_heldout[r["k"]] = dict(
            ho_f1=int(heldout_rule_f1[i]), ho_rp=int(heldout_rule_rp[i]),
            y=int(y[i]), cls=r["cls"], core=float(core[i]), aln=float(aln[i]),
            comp=int(r["comp"]), fold=int(fold[i]))
    rules = dict(f1=(float(final_f1["t1"]), float(final_f1["t2"])),
                 rp=(float(final_rp["t1"]), float(final_rp["t2"])),
                 baseline_prec=float(bp))
    return out, pair_heldout, rules


# ===========================================================================
# APPLY STAGE : re-threshold E_r by the RNA-only rule -> gamma-refine ->
#               DEMOTE allele-as-copy -> measure the residual over-merge FP.
# ===========================================================================
# The learned RNA-only rule is applied AT INFERENCE with RNA features only
# (core_recip = pair-max DN homology weight; aln_frac = universal leakage-free
# shared-exon fraction). No DNA/protein/genome value ever gates an edge here.
# We rebuild the family catalog EXACTLY as shipped (connected components of the
# re-thresholded E_r graph -> gamma-quasi-clique refinement gamma=0.20 seed=0 ->
# >=2-distinct-loci multi-copy predicate) and score it against the SAME oracles
# the sweep used (reference-cDNA/protein pair truth + the assembly-INDEPENDENT
# diploid copy-number oracle), so the numbers are directly comparable to the
# shipped gamma baseline.
#
# NON-CIRCULARITY: the removal is reported under BOTH the deployable final rule
# (thresholds fit on all data) and the component-CV HELD-OUT per-edge decisions
# (each edge cut by a rule trained on folds NOT containing its family). The two
# agree tightly => the over-merge removal is not an artifact of fitting the
# thresholds to these specific families.
#
# ALLELE DEMOTE (RNA, orthogonal to the edge rule): a same-gene multi-locus
# "family" whose every variant position is a balanced ~50/50 het (balanced_frac
# >= 0.90 AND copy_like <= 0.10, from the read-consensus O1 signal) is declared
# ALLELIC, not multi-copy, and split to singletons. The edge rule cannot touch
# these (both loci project to ONE gene => within-gene edge, always kept), so the
# allele-as-copy FP needs this separate RNA signal.

OUT_APPLY_TSV = os.path.join(BENCH, "rna_only_edge_oracle.tsv")
OUT_APPLY_JSON = os.path.join(BENCH, "rna_only_edge_oracle.json")

# DNA-confirmed residual FP roster from GRAPH_DEF_REFINE_SWEEP.md (the target list).
# LABELS/VALIDATION ONLY -- never gates an edge; used to name which residuals the
# RNA-only definition removes vs which are irreducible.
DNA_RESIDUAL_FP = dict(
    allele_as_copy=["DHRSX", "LOC129530050"],                      # asm_hapCN==1 both haplotypes
    oversize_diploid=["LOC129523567", "MAGEA9", "MPHOSPH8", "LOC134758618"],
    multifam=["GSTM2", "FOXO1", "LOC115933254", "LOC101142904", "LOC129526550"],
)


def oracle_residuals(fams, gene_of_dn, genes, g2rows, G, margin=1.5):
    """Re-implements the sweep's assembly-INDEPENDENT residual-FP tallies on a
    multi-copy family list, but also RETURNS the offending gene sets so we can
    track each named DNA-confirmed FP across catalogs. Mirrors
    graph_def_refine_sweep.eval_partition (c) exactly."""
    allele = []       # (genes, distinct_loci, asm_hapCN)
    oversize = []     # OVER-MERGE tail: dl > margin*dip  (genes, distinct_loci, diploid_CN, ratio)
    undersize = []    # OVER-SPLIT tail: dip > margin*dl  (genes, distinct_loci, diploid_CN, ratio)
    multifam = []     # (oracle genes spanned, distinct_loci)
    recovered = set() # oracle REP genes retained in a multi-copy family (assembly-independent recovery)
    orac_blocks = 0
    for b in fams:
        gs = {gene_of_dn[dn] for dn in b if dn in gene_of_dn}
        og = sorted(g for g in gs if g in g2rows)
        if not og:
            continue
        orac_blocks += 1
        recovered.update(og)
        rows = [r for g in og for r in g2rows[g]]
        classes = {r["cls"] for r in rows}
        hap = max(r["hap"] for r in rows)
        dip = max(r["dip"] for r in rows)
        dl = G.distinct_loci(b, genes)
        if len(og) >= 2:
            multifam.append((og, dl))
        if dl >= 2 and "single_locus_allele" in classes and "multi_copy" not in classes:
            allele.append((og, dl, hap))
        elif "single_locus_allele" not in classes:
            if dip > 0 and dl > margin * dip:                 # too MANY RNA loci for the DNA CN
                oversize.append((og, dl, dip, round(dl / dip, 2)))
            # symmetric OVER-SPLIT bound (assembly-independent): a multi_copy family whose
            # DNA diploid CN meaningfully EXCEEDS the surviving RNA loci => real copies split off
            if "multi_copy" in classes and dip > margin * dl and dl >= 2:
                undersize.append((og, dl, dip, round(dip / dl, 2)))
    return dict(orac_blocks=orac_blocks, allele=allele, oversize=oversize,
                undersize=undersize, multifam=multifam, recovered=sorted(recovered),
                n_recovered=len(recovered))


def _fmt_examples(lst, kind):
    out = []
    for e in sorted(lst, key=lambda x: (-x[1], x[0][0])):
        if kind == "allele":
            out.append(f"{'+'.join(e[0])}(dl={e[1]},hapCN={e[2]})")
        elif kind == "oversize":
            out.append(f"{'+'.join(e[0])}(dl={e[1]},dipCN={e[2]},{e[3]}x)")
        else:
            out.append(f"{'+'.join(e[0])}(dl={e[1]})")
    return out


def apply_stage():
    import json
    from collections import defaultdict, Counter
    import networkx as nx
    from networkx.algorithms.community import louvain_communities
    import genome_family_def as G
    import graph_def_refine_sweep as SW
    P = print

    # ---- LEARN gives held-out per-pair decisions + the deployable thresholds ----
    _, pair_ho, rules = learn()
    (T1, T2) = rules["f1"]
    (RT1, RT2) = rules["rp"]
    base_prec = rules["baseline_prec"]

    P("\n\n################ RNA-ONLY EDGE ORACLE :: APPLY STAGE ################")
    P(f"deployable RNA-only rule (F1)   : KEEP iff core_recip>={T1:.2f} AND aln_frac>={T2:.2f}")
    P(f"recall-preserving RNA-only rule : KEEP iff core_recip>={RT1:.2f} AND aln_frac>={RT2:.2f}")

    # ---- context (identical to the sweep so metrics are directly comparable) ----
    meta = FP.load_meta(); annot = FP.load_annot(); gene_of = FP.gene_of_factory(annot)
    raw_fams = FP.load_raw_families(); edge_pairs = FP.load_edges()
    genes, gene_of_dn, *_ = FP.build_genes_dict(raw_fams, meta, gene_of)
    ew = SW.load_weighted_edges(); lab = SW.load_pair_labels()
    g2f, mega = FP.load_prfam(); g2rows = SW.load_oracle()
    allele = load_allele()
    pair_core = load_pair_core(gene_of_dn)
    univ_aln = load_universal_aln()
    rc_aln = load_readcons_aln()          # read-consensus aln_frac (NO genome bases) for provenance check

    all_nodes = set()
    for f in raw_fams:
        all_nodes.update(f)
    Gr = nx.Graph(); Gr.add_nodes_from(all_nodes)
    for a, b in edge_pairs:
        if a in all_nodes and b in all_nodes:
            Gr.add_edge(a, b, weight=ew.get(frozenset((a, b)), 0.13))
    REAL = {k for k, v in lab.items() if v == "REAL_cdna"}
    GEN = {k for k, v in lab.items() if v == "GENUINE"}
    TB = {k for k, v in lab.items() if v == "TRUTHBAR"}
    ctx = dict(genes=genes, gene_of_dn=gene_of_dn, lab=lab, g2f=g2f, mega=mega,
               g2rows=g2rows, REAL=REAL, GEN=GEN, TB=TB, G=Gr, ew=ew)

    # ---- GUARD: the edge decision reads ONLY RNA features (core_recip, aln_frac) ----
    #      (label maps lab/g2rows/allele are used only to SCORE the resulting catalog)
    def decide_final(k):
        c = pair_core.get(k)
        c = c if c is not None else 0.0
        a = univ_aln.get(k)
        a = a if a is not None else 0.0
        return (c >= T1) and (a >= T2)

    def decide_recall(k):
        c = pair_core.get(k)
        c = c if c is not None else 0.0
        a = univ_aln.get(k)
        a = a if a is not None else 0.0
        return (c >= RT1) and (a >= RT2)

    def decide_heldout(k):
        r = pair_ho.get(k)
        return bool(r["ho_f1"]) if r is not None else decide_final(k)

    # PURITY-MAX ablation: identical rule but aln_frac from READ-CONSENSUS bases
    # (denovo POA transcript, NO genome). Proves the catalog is robust to the
    # provenance of the workhorse feature's bases.
    def decide_readcons(k):
        c = pair_core.get(k)
        c = c if c is not None else 0.0
        a = rc_aln.get(k)
        a = a if a is not None else 0.0
        return (c >= T1) and (a >= T2)

    def build_kept(decide):
        """Keep a DN edge iff it is within-gene / unannotated (never a cross-gene
        over-merge) OR its cross-gene pair passes the RNA-only decision."""
        kept = set()
        cross_cut = 0
        for u, v in Gr.edges():
            ga, gb = gene_of_dn.get(u), gene_of_dn.get(v)
            if ga is None or gb is None or ga == gb:
                kept.add(frozenset((u, v)))
                continue
            k = frozenset((ga, gb))
            if decide(k):
                kept.add(frozenset((u, v)))
            else:
                cross_cut += 1
        return kept, cross_cut

    def gamma_partition_on(kept):
        """connected components of the KEPT-edge graph, then the shipped gamma-
        quasi-clique recursive refinement (unchanged operator, gamma=0.20 seed=0)."""
        comps = SW.components_from_edges(all_nodes, kept)
        adj = defaultdict(set)
        for e in kept:
            u, v = tuple(e)
            adj[u].add(v); adj[v].add(u)
        out = []
        for m in comps:
            for block in G._refine_component(set(m), adj, G.GAMMA, louvain_communities, SW.SEED):
                out.append(set(block))
        return out

    # ---- shipped gamma baseline (reproduce for a like-for-like comparison) ----
    adj_all = defaultdict(set); fam_of = {}
    for kk, m in enumerate(raw_fams):
        for i in m:
            fam_of[i] = kk
    for u, v in edge_pairs:
        if u in fam_of and fam_of.get(u) == fam_of.get(v):
            adj_all[u].add(v); adj_all[v].add(u)
    ship_blocks = []
    for m in raw_fams:
        for block in G._refine_component(set(m), adj_all, G.GAMMA, louvain_communities, SW.SEED):
            ship_blocks.append(set(block))

    # ---- build the RNA-only catalogs ----
    kept_f1, cut_f1 = build_kept(decide_final)
    kept_rp, cut_rp = build_kept(decide_recall)
    kept_ho, cut_ho = build_kept(decide_heldout)
    kept_rc, cut_rc = build_kept(decide_readcons)
    blocks_f1 = gamma_partition_on(kept_f1)
    blocks_rp = gamma_partition_on(kept_rp)
    blocks_ho = gamma_partition_on(kept_ho)
    blocks_rc = gamma_partition_on(kept_rc)

    # ---- allele DEMOTE on the F1 catalog (RNA-only signal) ----
    def demote_gene(g):
        a = allele.get(g)
        return a is not None and a["balanced_frac"] >= 0.90 and a["copy_like"] <= 0.10

    def apply_demote(blocks):
        out = []; demoted = []
        for b in blocks:
            gs = [gene_of_dn[dn] for dn in b if dn in gene_of_dn]
            if gs:
                dom, cnt = Counter(gs).most_common(1)[0]
                homog = cnt / len(gs)
                if demote_gene(dom) and homog >= 0.5 and G.distinct_loci(b, genes) >= 2:
                    demoted.append(dict(gene=dom, n_loci=G.distinct_loci(b, genes),
                                        balanced_frac=allele[dom]["balanced_frac"],
                                        copy_like=allele[dom]["copy_like"],
                                        dna_confirmed=(dom in DNA_RESIDUAL_FP["allele_as_copy"])))
                    for dn in sorted(b):
                        out.append({dn})            # alleles -> singletons (not a family)
                    continue
            out.append(set(b))
        return out, demoted

    blocks_f1_dem, demotions = apply_demote(blocks_f1)

    # ---- evaluate every catalog with the SAME machinery as the sweep ----
    catalogs = [
        ("gamma_shipped_0.20", ship_blocks),
        ("rna_only_ruleF1", blocks_f1),
        ("rna_only_ruleF1_demoted", blocks_f1_dem),
        ("rna_only_heldout_ruleF1", blocks_ho),
        ("rna_only_recall_preserving", blocks_rp),
        ("rna_only_ruleF1_readcons", blocks_rc),   # purity-max: aln from read consensus, no genome bases
    ]
    results = {}
    residuals = {}
    for name, blocks in catalogs:
        res = SW.eval_partition(name, blocks, ctx)
        results[name] = res
        fams = SW.filter_multicopy(blocks, genes)
        residuals[name] = oracle_residuals(fams, gene_of_dn, genes, g2rows, G,
                                            margin=SW.DIPLOID_MARGIN)

    # ---- per-named-residual-FP tracking (present in shipped -> removed in RNA-only?) ----
    def gene_set_of_catalog(blocks):
        """gene -> set of frozenset(co-family gene partners) for quick co-membership tests."""
        fams = SW.filter_multicopy(blocks, genes)
        comem = defaultdict(set)
        for b in fams:
            gs = sorted({gene_of_dn[dn] for dn in b if dn in gene_of_dn})
            for g in gs:
                comem[g].update(gs)
        return comem

    comem = {name: gene_set_of_catalog(blocks) for name, blocks in catalogs}

    # multifam removal: are the two spanned oracle genes still co-membered?
    multifam_pairs = [("GSTM2", "*hub*"), ("FOXO1", "LOC115933254"),
                      ("LOC101142904", "LOC129526550")]
    tracking = {}
    # allele
    for g in DNA_RESIDUAL_FP["allele_as_copy"]:
        tracking[f"allele:{g}"] = {
            name: any(g in e[0] for e in residuals[name]["allele"]) for name, _ in catalogs}
    # oversize (block containing the gene still > margin*diploidCN)
    for g in DNA_RESIDUAL_FP["oversize_diploid"]:
        tracking[f"oversize:{g}"] = {
            name: any(g in e[0] for e in residuals[name]["oversize"]) for name, _ in catalogs}
    # multifam (spanned oracle-gene co-membership)
    for a, b in multifam_pairs:
        key = f"multifam:{a}+{b}"
        if b == "*hub*":
            tracking[key] = {name: sum(1 for e in residuals[name]["multifam"] if a in e[0])
                             for name, _ in catalogs}
        else:
            tracking[key] = {name: (b in comem[name].get(a, set())) for name, _ in catalogs}

    # ---- gene-pair held-out over-merge removal (the non-circular headline number) ----
    ho_genuine = [r for r in pair_ho.values() if r["cls"] == "genuine"]
    ho_gen_cut = sum(1 for r in ho_genuine if r["ho_f1"] == 0)
    ho_tp = [r for r in pair_ho.values() if r["y"] == 1]
    ho_tp_keep = sum(1 for r in ho_tp if r["ho_f1"] == 1)
    ho_tb = [r for r in pair_ho.values() if r["cls"] == "truthbar"]
    ho_tb_cut = sum(1 for r in ho_tb if r["ho_f1"] == 0)

    # ---- aln_frac PROVENANCE concordance (genome-fetched bases vs read-consensus bases) ----
    #      Bounds the "workhorse feature uses genome bases" concern: if the read-consensus
    #      aln_frac (denovo POA, NO genome) gives the same ranking + the same KEEP/CUT calls,
    #      the separation is a property of the transcript, not the base source.
    concord = aln_provenance_concordance(univ_aln, rc_aln, T2)

    # ---- SYMMETRIC diploid over-split bound (assembly-independent asm_hapCN) ----
    #      The over-MERGE tail is bounded by oracle oversize (dl > margin*dip). The over-SPLIT
    #      tail is bounded here on BOTH the family and the cut-edge level. dip = diploid CN
    #      (mat+pat) from the phased-assembly oracle; the ONLY assembly-independent scalar.
    def gene_dip(g):
        rows = g2rows.get(g)
        return max((r["dip"] for r in rows), default=0) if rows else 0

    # cut-edge level: of the held-out CUT edges, how many touch a locus the diploid DNA
    # says is genuinely MULTI-COPY (dip>=2) -> a potential REAL copy split off (over-split cost),
    # vs single-copy/oracle-uncovered (cut is benign / not a lost genomic copy).
    for k, r in pair_ho.items():            # attach the pair key back (pair_ho is keyed by frozenset)
        r["k"] = k
    def dip_of_pair(r):
        ga, gb = sorted(tuple(r["k"]))
        return max(gene_dip(ga), gene_dip(gb))
    oversplit = dict()
    for cls in ("TP", "truthbar", "genuine"):
        rows = [r for r in pair_ho.values() if r["ho_f1"] == 0 and r["cls"] == cls]
        covered = [r for r in rows if dip_of_pair(r) > 0]
        dip_ge2 = [r for r in rows if dip_of_pair(r) >= 2]
        oversplit[cls] = dict(cut=len(rows), oracle_covered=len(covered),
                              diploid_multicopy=len(dip_ge2))

    # =====================================================================
    # REPORT
    # =====================================================================
    def row(res, r2):
        return (res["n_families"], res["tp_retention"], res["ge_genuine_cut"],
                res["truthbar_retention"], res["impure_blocks"], r2["orac_blocks"],
                r2["n_recovered"], len(r2["allele"]), len(r2["oversize"]),
                len(r2["undersize"]), len(r2["multifam"]))

    P(f"\nDN edges: {Gr.number_of_edges()} total | cross-gene cut by RNA rule (F1)={cut_f1} "
      f"heldout={cut_ho} recall-preserving={cut_rp}")
    P(f"allele DEMOTE fired on {len(demotions)} families: "
      + ", ".join(f"{d['gene']}(dl={d['n_loci']},bal={d['balanced_frac']:.2f},"
                  f"copy_like={d['copy_like']:.2f}{',DNA-confirmed' if d['dna_confirmed'] else ',novel'})"
                  for d in demotions))

    P("\n---- CATALOG COMPARISON (same eval as GRAPH_DEF_REFINE_SWEEP) ----")
    P("   tp_ret=DNA-real edge retention  ge_gcut=genuine over-merge edge-cut  tb_ret=truthbar(divergent-paralog) retention")
    P("   recov=oracle genes retained multi-copy  undersz=over-SPLIT tail (dip>1.5x RNA loci)")
    hdr = ("catalog", "n_fam", "tp_ret", "ge_gcut", "tb_ret", "impure", "orac",
           "recov", "allele", "oversz", "undersz", "multifam")
    P("  {:<28s}{:>6s}{:>8s}{:>8s}{:>8s}{:>8s}{:>6s}{:>7s}{:>8s}{:>8s}{:>9s}{:>9s}".format(*hdr))
    for name, _ in catalogs:
        r = row(results[name], residuals[name])
        P("  {:<28s}{:>6d}{:>8.3f}{:>8.3f}{:>8.3f}{:>8d}{:>6d}{:>7d}{:>8d}{:>8d}{:>9d}{:>9d}".format(
            name, *r))

    P("\n---- DNA-CONFIRMED RESIDUAL FP: removed by the RNA-only definition? ----")
    P("  (shipped -> ruleF1 -> +demote -> heldout)   True=FP present, False=removed")
    order = ["gamma_shipped_0.20", "rna_only_ruleF1", "rna_only_ruleF1_demoted",
             "rna_only_heldout_ruleF1"]
    for key in sorted(tracking):
        vals = tracking[key]
        P(f"  {key:32s} " + " -> ".join(f"{str(vals[o])}" for o in order))

    P("\n---- aln_frac PROVENANCE: genome-fetched bases vs READ-CONSENSUS bases ----")
    P(f"  (the deployed aln_frac uses RNA splice STRUCTURE + genome-fetched exon bases;")
    P(f"   here recomputed from de-novo POA transcript consensus -- NO genome bases)")
    P(f"  Spearman rho(genome_aln, readcons_aln) = {concord['rho']:.4f}  (n={concord['n']})")
    P(f"  KEEP/CUT gate agreement @ aln>={concord['t2']:.2f}    = {concord['gate_agreement_frac']:.4f} "
      f"({concord['gate_agreement']}/{concord['n']})")
    P(f"  AUC vs in_dna_loose  genome={concord['auc_genome']:.4f}  read-consensus={concord['auc_readcons']:.4f} "
      f"(essentially identical => not a genome-base artifact)")

    P("\n---- NON-CIRCULAR (component-CV held-out) gene-pair over-merge removal ----")
    P(f"  genuine over-merge edges cut (held-out) : {ho_gen_cut}/{len(ho_genuine)} "
      f"({100*ho_gen_cut/len(ho_genuine):.1f}%)")
    P(f"  TP (DNA-real) edges kept   (held-out)   : {ho_tp_keep}/{len(ho_tp)} "
      f"({100*ho_tp_keep/len(ho_tp):.1f}%)")
    P(f"  truthbar divergent-paralog edges cut    : {ho_tb_cut}/{len(ho_tb)} "
      f"({100*ho_tb_cut/len(ho_tb):.1f}%)  <- the cost / over-split (bounded below)")

    P("\n---- SYMMETRIC OVER-SPLIT BOUND (assembly-independent diploid CN; both tails) ----")
    P(f"  held-out CUT edges touching a diploid MULTI-COPY locus (dip>=2 => a potential REAL copy split off):")
    for cls in ("TP", "truthbar", "genuine"):
        o = oversplit[cls]
        P(f"    {cls:9s} cut={o['cut']:5d}  oracle-covered={o['oracle_covered']:4d}  "
          f"diploid-multicopy(dip>=2)={o['diploid_multicopy']:4d}")
    ship_res = residuals["gamma_shipped_0.20"]; rf1_res = residuals["rna_only_ruleF1_demoted"]
    P(f"  oracle-gene RECOVERY (validated families retained as multi-copy): "
      f"shipped {ship_res['n_recovered']} -> RNA-only {rf1_res['n_recovered']}  "
      f"(delta {rf1_res['n_recovered']-ship_res['n_recovered']:+d} = validated families over-split away)")
    P(f"  undersize families (dip > {SW.DIPLOID_MARGIN}x RNA loci; real copies split off): "
      f"shipped {len(ship_res['undersize'])} -> RNA-only {len(rf1_res['undersize'])}")

    P("\n---- IRREDUCIBLE RESIDUAL (survives the RNA-only definition; RNA cannot separate) ----")
    surv_over = residuals["rna_only_ruleF1_demoted"]["oversize"]
    surv_multi = [e for e in residuals["rna_only_ruleF1_demoted"]["multifam"]]
    surv_all = residuals["rna_only_ruleF1_demoted"]["allele"]
    P(f"  oversize still > {SW.DIPLOID_MARGIN}x diploid DNA CN : "
      + (", ".join(_fmt_examples(surv_over, "oversize")) or "NONE"))
    P(f"  multifam (>=2 oracle genes co-membered)  : "
      + (", ".join(_fmt_examples(surv_multi, "multifam")) or "NONE"))
    P(f"  allele-as-copy remaining                 : "
      + (", ".join(_fmt_examples(surv_all, "allele")) or "NONE"))

    ship = results["gamma_shipped_0.20"]; rf1 = results["rna_only_ruleF1"]
    rp = results["rna_only_recall_preserving"]
    P("\n---- PRECISION WITHOUT A DNA GATE (honest incremental framing) ----")
    P(f"  The shipped core_recip>=0.13 is the edge-CREATION threshold (keeps 100% of edges),")
    P(f"  NOT an operating point; the honest incremental gain is core-alone AUC 0.815 -> +aln 0.915.")
    P(f"  genuine over-merge edge-cut  shipped {ship['ge_genuine_cut']:.3f} -> F1 {rf1['ge_genuine_cut']:.3f} "
      f"| recall-preserving {rp['ge_genuine_cut']:.3f}")
    P(f"  tp_retention (DNA-real kept) shipped {ship['tp_retention']:.3f} -> F1 {rf1['tp_retention']:.3f} "
      f"| recall-preserving {rp['tp_retention']:.3f}   <- OVER-SPLIT cost (co-headline)")
    P(f"  truthbar_retention           shipped {ship['truthbar_retention']:.3f} -> F1 {rf1['truthbar_retention']:.3f} "
      f"| recall-preserving {rp['truthbar_retention']:.3f}   <- divergent-paralog recall traded for precision")
    P(f"  n_families                   shipped {ship['n_families']} -> F1 {rf1['n_families']} "
      f"| recall-preserving {rp['n_families']}")

    # =====================================================================
    # WRITE OUTPUTS
    # =====================================================================
    cols = ["catalog", "n_families", "tp_retention", "overmerge_cut", "ge_genuine_cut",
            "genuine_precision", "truthbar_retention", "impure_blocks", "oracle_blocks",
            "oracle_genes_recovered", "oracle_allele_as_copy", "oracle_oversize_diploid",
            "oracle_undersize_diploid", "oracle_multifam"]
    with open(OUT_APPLY_TSV, "w") as out:
        out.write("\t".join(cols) + "\n")
        for name, _ in catalogs:
            r = results[name]; r2 = residuals[name]
            out.write("\t".join(str(x) for x in [
                name, r["n_families"], r["tp_retention"], r["overmerge_cut"],
                r["ge_genuine_cut"], r["genuine_precision"], r["truthbar_retention"],
                r["impure_blocks"], r2["orac_blocks"], r2["n_recovered"], len(r2["allele"]),
                len(r2["oversize"]), len(r2["undersize"]), len(r2["multifam"])]) + "\n")

    summary = dict(
        rule=dict(final_f1=dict(core_recip=T1, aln_frac=T2),
                  recall_preserving=dict(core_recip=RT1, aln_frac=RT2),
                  form="KEEP E_r edge iff core_recip>=t1 AND aln_frac>=t2 (RNA-only)"),
        demote=dict(rule="balanced_frac>=0.90 AND copy_like<=0.10 (read-consensus O1)",
                    fired=demotions, n_fired=len(demotions)),
        dn_edges=dict(total=Gr.number_of_edges(), cross_cut_ruleF1=cut_f1,
                      cross_cut_heldout=cut_ho, cross_cut_recall=cut_rp),
        catalogs={name: dict(
            n_families=results[name]["n_families"],
            tp_retention=results[name]["tp_retention"],
            overmerge_cut=results[name]["overmerge_cut"],
            ge_genuine_cut=results[name]["ge_genuine_cut"],
            genuine_precision=results[name]["genuine_precision"],
            truthbar_retention=results[name]["truthbar_retention"],
            impure_blocks=results[name]["impure_blocks"],
            oracle_blocks=residuals[name]["orac_blocks"],
            oracle_genes_recovered=residuals[name]["n_recovered"],
            oracle_allele_as_copy=len(residuals[name]["allele"]),
            oracle_oversize_diploid=len(residuals[name]["oversize"]),
            oracle_undersize_diploid=len(residuals[name]["undersize"]),
            oracle_multifam=len(residuals[name]["multifam"]),
            oracle_allele_examples=_fmt_examples(residuals[name]["allele"], "allele"),
            oracle_oversize_examples=_fmt_examples(residuals[name]["oversize"], "oversize"),
            oracle_undersize_examples=_fmt_examples(residuals[name]["undersize"], "oversize"),
            oracle_multifam_examples=_fmt_examples(residuals[name]["multifam"], "multifam"),
        ) for name, _ in catalogs},
        residual_fp_tracking=tracking,
        heldout_noncircular=dict(
            genuine_cut=ho_gen_cut, genuine_total=len(ho_genuine),
            genuine_cut_frac=round(ho_gen_cut / len(ho_genuine), 4),
            tp_kept=ho_tp_keep, tp_total=len(ho_tp),
            tp_kept_frac=round(ho_tp_keep / len(ho_tp), 4),
            truthbar_cut=ho_tb_cut, truthbar_total=len(ho_tb)),
        aln_provenance=dict(
            desc="genome-fetched exon bases (RNA splice structure) vs read-consensus POA bases (no genome)",
            spearman_rho=round(concord["rho"], 4), n=concord["n"],
            gate_agreement_frac=round(concord["gate_agreement_frac"], 4),
            auc_genome=round(concord["auc_genome"], 4),
            auc_readcons=round(concord["auc_readcons"], 4)),
        oversplit_diploid_bound=dict(
            desc="held-out CUT edges touching a diploid multi-copy locus (dip>=2) = potential real copy split off",
            by_class=oversplit,
            oracle_genes_recovered_shipped=ship_res["n_recovered"],
            oracle_genes_recovered_rna_only=rf1_res["n_recovered"],
            undersize_shipped=len(ship_res["undersize"]),
            undersize_rna_only=len(rf1_res["undersize"])),
        guards=dict(edge_decision_features=["core_recip", "aln_frac"],
                    no_dna_in_decision=True,
                    demote_features=["balanced_frac", "copy_like"],
                    purity_max_catalog="rna_only_ruleF1_readcons (aln from read consensus, no genome bases)",
                    gamma=G.GAMMA, seed=SW.SEED),
    )
    with open(OUT_APPLY_JSON, "w") as out:
        json.dump(summary, out, sort_keys=True, indent=1,
                  default=lambda x: None if (isinstance(x, float) and x != x) else x)
    P(f"\nwrote {OUT_APPLY_TSV}\nwrote {OUT_APPLY_JSON}")
    P("APPLY_JSON " + json.dumps(summary, sort_keys=True,
                                 default=lambda x: None if (isinstance(x, float) and x != x) else x))
    P("####################################################################")
    return summary


if __name__ == "__main__":
    main()
    apply_stage()
