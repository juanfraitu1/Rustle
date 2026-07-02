#!/usr/bin/env python3
"""
block_graph_fp_detector.py  --  PER-FAMILY-BLOCK graph-signature FP triage.

QUESTION (untested angle): can a family's WHOLE-BLOCK graph signature
(multiplicity + connectivity + density + CARDINALITY, controlling for size)
predict that the block is an OVER-MERGE (false positive)?  This is a per-BLOCK
CONFIDENCE / TRIAGE score (flag families for DNA review), NOT an edge cutter.

Prior work that this does NOT repeat:
  (a) per-EDGE node-multiplicity catches the REPEAT-HUB class -> shipped fam17
      repeat gate (family_rna_refine.py, min_shared_mult>=20).
  (b) per-EDGE connectivity/topology is FALSIFIED for the general over-merge
      (over-merges are DENSE not sparse; triangle-support AUC INVERTED 0.28) --
      bench/graph_def_refine_sweep.py.

TRUTH:
  * TRAIN truth  = E_p-impurity: block spans >=2 distinct NON-mega protein
    families (partly reference-derived; large enough to fit -- 80/606 positives).
    Reproduced EXACTLY from family_level_pr_current.impure_block_roster logic.
  * GOLD truth   = diploid_cn_oracle: block is a DNA-confirmed over-merge
    (multifam / oversize / allele), reproduced from rna_only_edge_oracle
    .oracle_residuals logic.  Only ~6 FPs over ~50 oracle-mapped families ->
    TOO SMALL TO TRAIN; used ONLY to VALIDATE the triage flags.

RNA-only features (no DNA in the detector; DNA = validation only).
Deterministic: PYTHONHASHSEED=0, seed=0 folds/Louvain/Fiedler.
"""
import os, sys, json, csv
from collections import defaultdict
from statistics import median

import numpy as np
import networkx as nx
from networkx.algorithms.community import louvain_communities

BENCH = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, BENCH)
import family_level_pr_current as M     # build_ctx, load_catalog_blocks
import genome_family_def as G           # distinct_loci
import graph_def_refine_sweep as SW     # DIPLOID_MARGIN

SEED = 0
np.random.seed(SEED)
HI_CORE = 0.50          # "near-identical" edge threshold for sequence-cluster cardinality
REPEAT_MULT_MIN = 20    # the shipped fam17 repeat-hub gate cut (VG min_shared_mult)
VG_TSV = os.path.join(BENCH, "vg_repeat_catalog.tsv")
OUT_TSV = os.path.join(BENCH, "block_graph_fp_detector.tsv")
OUT_JSON = os.path.join(BENCH, "block_graph_fp_detector.json")


# --------------------------------------------------------------------------- VG edge multiplicity catalog
def load_vg_edges():
    """gene-pair (frozenset) -> dict(min_mult,max_mult,shares_only_repeat_m5) from the
    per-edge section of vg_repeat_catalog.tsv.  Library-free canonical-minimizer multiplicity."""
    out = {}
    with open(VG_TSV) as fh:
        in_edges, ix = False, None
        for ln in fh:
            if ln.startswith("# SECTION edges"):
                in_edges = True; continue
            if not in_edges:
                continue
            if ln.startswith("gene_a\t"):
                hdr = ln.rstrip("\n").split("\t"); ix = {h: i for i, h in enumerate(hdr)}
                continue
            if ix is None:
                continue
            f = ln.rstrip("\n").split("\t")
            try:
                mn = int(f[ix["min_shared_mult"]]); mx = int(f[ix["max_shared_mult"]])
            except (ValueError, KeyError):
                continue
            sor = f[ix["shares_only_repeat_m5"]] if "shares_only_repeat_m5" in ix else "0"
            out[frozenset((f[ix["gene_a"]], f[ix["gene_b"]]))] = dict(
                mn=mn, mx=mx, sor=(sor == "1"))
    return out


# --------------------------------------------------------------------------- per-block features
def block_features(block, ctx, vg, Gfull):
    """All RNA/graph-derived per-block features. No DNA, no protein-family COUNT
    (protein-family count IS the E_p label -> excluded from the classifier; kept
    only as report-only leakage sentinel)."""
    genes = ctx["genes"]; god = ctx["gene_of_dn"]
    nodes = [d for d in block]
    n_mem = len(nodes)
    gene_set = sorted({god[d] for d in nodes if d in god})
    n_genes = len(gene_set)
    dloci = G.distinct_loci(block, genes)

    # induced DN subgraph (unweighted) + core_recip weights
    H = nx.Graph(); H.add_nodes_from(nodes)
    Hhi = nx.Graph(); Hhi.add_nodes_from(nodes)          # near-identical (core>=HI) subgraph
    n_edges = 0
    mults_mx, mults_mn = [], []
    n_repeat_bridged = 0           # shares_only_repeat_m5
    n_ge_gate = 0                  # min_shared_mult >= 20  (the shipped gate; should be ~0 post-gate)
    for u, v in Gfull.edges(nodes):
        if u in H and v in H:
            w = ctx["ew"].get(frozenset((u, v)), 0.13)
            H.add_edge(u, v, weight=w)
            n_edges += 1
            if w >= HI_CORE:
                Hhi.add_edge(u, v)
            gu, gv = god.get(u), god.get(v)
            if gu and gv and gu != gv:
                rec = vg.get(frozenset((gu, gv)))
                if rec is not None:
                    mults_mx.append(rec["mx"]); mults_mn.append(rec["mn"])
                    if rec["sor"]: n_repeat_bridged += 1
                    if rec["mn"] >= REPEAT_MULT_MIN: n_ge_gate += 1

    poss = n_mem * (n_mem - 1) / 2.0
    density = (n_edges / poss) if poss > 0 else 0.0
    degs = [d for _, d in H.degree()]
    max_deg = max(degs) if degs else 0
    mean_deg = float(np.mean(degs)) if degs else 0.0
    clustering = nx.average_clustering(H) if n_mem >= 3 else 0.0

    # Louvain sub-communities (deterministic seed)
    if n_edges > 0:
        comms = louvain_communities(H, seed=SEED, weight="weight")
        n_comm = len(comms)
    else:
        n_comm = n_mem

    # algebraic connectivity (Fiedler) + weighted min-cut, on the giant component
    fiedler = 0.0; mincut = 0.0
    if n_mem >= 2 and n_edges > 0:
        comps = list(nx.connected_components(H))
        gc = max(comps, key=len)
        if len(gc) >= 2:
            sub = H.subgraph(gc)
            try:
                fiedler = float(nx.algebraic_connectivity(sub, weight="weight",
                                                          seed=SEED, method="tracemin_lu"))
            except Exception:
                fiedler = 0.0
            if len(gc) >= 2:
                try:
                    mincut = float(nx.stoer_wagner(sub, weight="weight")[0])
                except Exception:
                    mincut = 0.0
        frac_in_gc = len(gc) / n_mem
    else:
        frac_in_gc = 1.0 if n_mem == 1 else (1.0 if n_edges > 0 else 1.0 / max(n_mem, 1))

    # CARDINALITY: near-identical sequence clusters (core>=HI single-linkage components)
    n_seqclust_hi = nx.number_connected_components(Hhi)
    loci_per_seqclust = dloci / n_seqclust_hi if n_seqclust_hi > 0 else dloci
    mem_per_seqclust = n_mem / n_seqclust_hi if n_seqclust_hi > 0 else n_mem
    loci_per_comm = dloci / n_comm if n_comm > 0 else dloci
    loci_per_gene = dloci / n_genes if n_genes > 0 else 0.0
    mem_per_locus = n_mem / dloci if dloci > 0 else 0.0

    # MULTIPLICITY
    max_mult = max(mults_mx) if mults_mx else 0
    med_mult = median(mults_mx) if mults_mx else 0
    max_min_mult = max(mults_mn) if mults_mn else 0
    cross_edges = len(mults_mx)
    frac_repeat_bridged = n_repeat_bridged / cross_edges if cross_edges else 0.0
    frac_ge_gate = n_ge_gate / cross_edges if cross_edges else 0.0

    return dict(
        n_members=n_mem, n_genes=n_genes, distinct_loci=dloci, n_edges=n_edges,
        # connectivity / density
        edge_density=density, max_degree=max_deg, mean_degree=mean_deg,
        clustering=clustering, n_communities=n_comm, fiedler=fiedler,
        mincut=mincut, frac_in_giant=frac_in_gc,
        # cardinality
        n_seqclust_hi=n_seqclust_hi, loci_per_seqclust=loci_per_seqclust,
        mem_per_seqclust=mem_per_seqclust, loci_per_comm=loci_per_comm,
        loci_per_gene=loci_per_gene, mem_per_locus=mem_per_locus,
        # multiplicity
        max_mult=max_mult, med_mult=med_mult, max_min_mult=max_min_mult,
        frac_repeat_bridged=frac_repeat_bridged, frac_ge_gate=frac_ge_gate,
    )


# --------------------------------------------------------------------------- truth labels
def ep_impure(block, ctx):
    """TRAIN truth: block spans >=2 distinct NON-mega protein families.
    Exact reproduction of family_level_pr_current.impure_block_roster."""
    god = ctx["gene_of_dn"]; g2f = ctx["g2f"]; mega = ctx["mega"]
    gs = {god[d] for d in block if d in god}
    fams = set()
    for g in gs:
        if g in g2f and g2f[g] not in mega:
            fams.add(g2f[g])
    return (len(fams) >= 2, len(fams))


def oracle_fp(block, ctx, margin):
    """GOLD truth: DNA-confirmed over-merge (multifam / oversize / allele).
    Exact per-block reproduction of rna_only_edge_oracle.oracle_residuals."""
    god = ctx["gene_of_dn"]; g2rows = ctx["g2rows"]; genes = ctx["genes"]
    gs = {god[d] for d in block if d in god}
    og = sorted(g for g in gs if g in g2rows)
    if not og:
        return dict(mapped=False, fp=False, cls="")
    rows = [r for g in og for r in g2rows[g]]
    classes = {r["cls"] for r in rows}
    hap = max(r["hap"] for r in rows); dip = max(r["dip"] for r in rows)
    dl = G.distinct_loci(block, genes)
    tags = []
    if len(og) >= 2:
        tags.append("multifam")
    if dl >= 2 and "single_locus_allele" in classes and "multi_copy" not in classes:
        tags.append("allele")
    elif "single_locus_allele" not in classes:
        if dip > 0 and dl > margin * dip:
            tags.append("oversize")
    return dict(mapped=True, fp=bool(tags), cls="+".join(tags), og=og, dip=dip, dl=dl)


# --------------------------------------------------------------------------- metrics
def auc(scores, labels):
    """Mann-Whitney AUC. scores/labels are arrays; ties handled by rank average."""
    s = np.asarray(scores, float); y = np.asarray(labels, int)
    pos = s[y == 1]; neg = s[y == 0]
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    order = np.argsort(s, kind="mergesort")
    ranks = np.empty(len(s), float)
    sr = s[order]
    i = 0
    while i < len(sr):
        j = i
        while j + 1 < len(sr) and sr[j + 1] == sr[i]:
            j += 1
        ranks[order[i:j + 1]] = (i + j) / 2.0 + 1.0
        i = j + 1
    r_pos = ranks[y == 1].sum()
    n1 = len(pos); n0 = len(neg)
    return (r_pos - n1 * (n1 + 1) / 2.0) / (n1 * n0)


def residualize(feat, size_cols):
    """OLS residual of feat on [1, size_cols...]; removes the linear size trend."""
    X = np.column_stack([np.ones(len(feat))] + [np.asarray(c, float) for c in size_cols])
    beta, *_ = np.linalg.lstsq(X, np.asarray(feat, float), rcond=None)
    return np.asarray(feat, float) - X @ beta


class Logit:
    """L2-regularised logistic regression (numpy Newton-IRLS). No sklearn."""
    def __init__(self, lam=1.0, iters=100):
        self.lam = lam; self.iters = iters; self.mu = None; self.sd = None; self.w = None

    def fit(self, X, y):
        X = np.asarray(X, float); y = np.asarray(y, float)
        self.mu = X.mean(0); self.sd = X.std(0); self.sd[self.sd == 0] = 1.0
        Xs = np.column_stack([np.ones(len(X)), (X - self.mu) / self.sd])
        w = np.zeros(Xs.shape[1])
        R = self.lam * np.eye(Xs.shape[1]); R[0, 0] = 0.0
        for _ in range(self.iters):
            p = 1.0 / (1.0 + np.exp(-Xs @ w))
            Wd = p * (1 - p) + 1e-9
            grad = Xs.T @ (p - y) + R @ w
            Hess = Xs.T @ (Xs * Wd[:, None]) + R
            step = np.linalg.solve(Hess, grad)
            w = w - step
            if np.max(np.abs(step)) < 1e-8:
                break
        self.w = w; return self

    def predict(self, X):
        X = np.asarray(X, float)
        Xs = np.column_stack([np.ones(len(X)), (X - self.mu) / self.sd])
        return 1.0 / (1.0 + np.exp(-Xs @ self.w))


# --------------------------------------------------------------------------- main
def main():
    print("[load] ctx + catalog ...", flush=True)
    ctx, raw_fams, edge_pairs = M.build_ctx()
    cat, ids, dom, member_gene = M.load_catalog_blocks()
    Gfull = ctx["G"]
    vg = load_vg_edges()
    margin = SW.DIPLOID_MARGIN
    print(f"       {len(cat)} families; VG edge rows {len(vg)}; margin {margin}", flush=True)

    rows = []
    for i, b in enumerate(cat):
        fid = ids[i]
        feat = block_features(b, ctx, vg, Gfull)
        imp, n_pf = ep_impure(b, ctx)
        orc = oracle_fp(b, ctx, margin)
        gs = sorted({ctx["gene_of_dn"][d] for d in b if d in ctx["gene_of_dn"]})
        rows.append(dict(family_id=fid, dominant_gene=dom.get(fid, ""),
                         genes=";".join(gs), n_protein_families=n_pf,
                         is_impure=int(imp), oracle_mapped=int(orc["mapped"]),
                         oracle_fp=int(orc["fp"]), oracle_cls=orc["cls"], **feat))

    N = len(rows)
    y = np.array([r["is_impure"] for r in rows])
    print(f"       E_p-impure positives: {y.sum()}/{N}", flush=True)

    # -- E_p impurity rate by size (the confound) --
    def rate_by(pred):
        idx = [i for i in range(N) if pred(rows[i])]
        return (len(idx), sum(rows[i]["is_impure"] for i in idx))
    size_bins = {
        "n_loci=2": lambda r: r["distinct_loci"] <= 2,
        "n_loci 3-4": lambda r: 3 <= r["distinct_loci"] <= 4,
        "n_loci 5-9": lambda r: 5 <= r["distinct_loci"] <= 9,
        "n_loci>=10": lambda r: r["distinct_loci"] >= 10,
    }
    conf = {k: rate_by(f) for k, f in size_bins.items()}

    # ---------------- per-feature AUC (raw + size-residualized) ----------------
    FEATS = ["n_members", "n_genes", "distinct_loci", "n_edges",
             "edge_density", "max_degree", "mean_degree", "clustering",
             "n_communities", "fiedler", "mincut", "frac_in_giant",
             "n_seqclust_hi", "loci_per_seqclust", "mem_per_seqclust",
             "loci_per_comm", "loci_per_gene", "mem_per_locus",
             "max_mult", "med_mult", "max_min_mult",
             "frac_repeat_bridged", "frac_ge_gate"]
    logn = np.log(np.array([r["n_members"] for r in rows], float))
    logdl = np.log(np.array([max(r["distinct_loci"], 1) for r in rows], float))
    # A ratio feature (e.g. mem_per_locus = n_members/distinct_loci, or
    # frac_in_giant) is a NONLINEAR function of the size variables, so a LINEAR
    # residualizer on [logn, logdl] cannot remove it -> spurious residual-AUC.
    # Second (quadratic) basis exposes these artifacts: genuine size-independent
    # signal survives the richer basis; residualization artifacts collapse to ~0.5.
    quad = [logn, logdl, logn * logn, logdl * logdl, logn * logdl]
    feat_auc = {}
    for f in FEATS:
        v = np.array([r[f] for r in rows], float)
        a_raw = auc(v, y)
        a_res = auc(residualize(v, [logn, logdl]), y)
        a_res_q = auc(residualize(v, quad), y)
        feat_auc[f] = dict(auc_raw=round(a_raw, 3), auc_size_resid=round(a_res, 3),
                           auc_size_resid_quad=round(a_res_q, 3))
    # size baseline AUC (best single size feature)
    size_auc = max(feat_auc["n_members"]["auc_raw"], feat_auc["distinct_loci"]["auc_raw"],
                   feat_auc["n_genes"]["auc_raw"])

    # ---------------- CV logistic: size-only vs size+graph ----------------
    # graph feature set EXCLUDES protein-family count (=label) and raw size proxies
    GRAPH_FEATS = ["edge_density", "clustering", "n_communities", "fiedler", "mincut",
                   "loci_per_seqclust", "mem_per_seqclust", "loci_per_comm",
                   "loci_per_gene", "mem_per_locus", "max_mult", "med_mult",
                   "max_min_mult", "frac_repeat_bridged", "n_seqclust_hi", "max_degree"]
    SIZE_FEATS = ["n_members", "distinct_loci", "n_genes"]

    def matrix(cols):
        return np.array([[r[c] for c in cols] for r in rows], float)

    rng = np.random.RandomState(SEED)
    fold = rng.randint(0, 5, size=N)
    def cv_scores(cols):
        oof = np.zeros(N)
        for k in range(5):
            tr = fold != k; te = fold == k
            m = Logit(lam=2.0).fit(matrix(cols)[tr], y[tr])
            oof[te] = m.predict(matrix(cols)[te])
        return oof
    oof_size = cv_scores(SIZE_FEATS)
    oof_full = cv_scores(SIZE_FEATS + GRAPH_FEATS)
    oof_graph = cv_scores(GRAPH_FEATS)
    auc_size = auc(oof_size, y)
    auc_full = auc(oof_full, y)
    auc_graph = auc(oof_graph, y)

    # final full-fit FP-score (in-sample, for the shipped tsv) + held-out oof for honesty
    final = Logit(lam=2.0).fit(matrix(SIZE_FEATS + GRAPH_FEATS), y)
    fp_score = final.predict(matrix(SIZE_FEATS + GRAPH_FEATS))
    for i, r in enumerate(rows):
        r["fp_score"] = round(float(fp_score[i]), 5)
        r["fp_score_oof"] = round(float(oof_full[i]), 5)

    # ---------------- diploid-oracle triage validation (precision@K) ----------------
    orc_fp_idx = [i for i in range(N) if rows[i]["oracle_fp"]]
    orc_map_idx = [i for i in range(N) if rows[i]["oracle_mapped"]]
    n_orc_fp = len(orc_fp_idx)
    def precision_at_k(scorevec, Ks):
        order = np.argsort(-scorevec, kind="mergesort")
        res = {}
        for K in Ks:
            topK = set(order[:K].tolist())
            n_fp_in = len(topK & set(orc_fp_idx))
            n_imp_in = sum(rows[i]["is_impure"] for i in topK)
            res[K] = dict(oracle_fp_captured=f"{n_fp_in}/{n_orc_fp}",
                          recall_oracle_fp=round(n_fp_in / n_orc_fp, 3) if n_orc_fp else None,
                          ep_impure_in_topK=f"{n_imp_in}/{K}",
                          precision_ep=round(n_imp_in / K, 3))
        return res
    Ks = [10, 20, 30, 50]
    triage_insample = precision_at_k(fp_score, Ks)
    triage_heldout = precision_at_k(oof_full, Ks)
    # rank of each oracle FP under the held-out score
    order_oof = np.argsort(-oof_full, kind="mergesort")
    rank_of = {int(order_oof[p]): p + 1 for p in range(N)}
    oracle_fp_ranks = [dict(family_id=rows[i]["family_id"],
                            dominant=rows[i]["dominant_gene"],
                            cls=rows[i]["oracle_cls"], distinct_loci=rows[i]["distinct_loci"],
                            fp_score=rows[i]["fp_score"],
                            rank_heldout=rank_of[i], rank_of_N=N)
                       for i in sorted(orc_fp_idx, key=lambda i: rank_of[i])]

    # ---------------- class-by-class flagging ----------------
    def find_block(pred):
        return [rows[i] for i in range(N) if pred(rows[i])]
    def named(family_id):
        for r in rows:
            if r["family_id"] == str(family_id):
                return r
        return None
    # MAGE cardinality over-merge = block 549 (MAGEA1/4/12 + LOC129529978/86)
    mage = named(549)
    # fam17 post-gate remnant = block 17 (20 genes / 10 protein families)
    fam17 = named(17)
    # GSTM2 domain-hub multifam (spans GSTM2 + LOCs)
    gstm2 = [r for r in rows if "GSTM2" in r["genes"].split(";")]
    foxo1 = [r for r in rows if "FOXO1" in r["genes"].split(";")]

    def brief(r):
        if r is None:
            return None
        return dict(family_id=r["family_id"], dominant=r["dominant_gene"],
                    distinct_loci=r["distinct_loci"], n_genes=r["n_genes"],
                    n_protein_families=r["n_protein_families"], is_impure=r["is_impure"],
                    oracle_fp=r["oracle_fp"], oracle_cls=r["oracle_cls"],
                    loci_per_seqclust=round(r["loci_per_seqclust"], 2),
                    n_seqclust_hi=r["n_seqclust_hi"], max_min_mult=r["max_min_mult"],
                    max_mult=r["max_mult"], frac_repeat_bridged=round(r["frac_repeat_bridged"], 2),
                    edge_density=round(r["edge_density"], 3), clustering=round(r["clustering"], 3),
                    fp_score=r["fp_score"], rank_heldout=rank_of[rows.index(r)])

    # unsupervised cardinality flag: does loci_per_seqclust single out MAGE?
    lps = np.array([r["loci_per_seqclust"] for r in rows])
    card_order = np.argsort(-lps, kind="mergesort")
    top_card = [dict(family_id=rows[i]["family_id"], dominant=rows[i]["dominant_gene"],
                     loci_per_seqclust=round(rows[i]["loci_per_seqclust"], 2),
                     distinct_loci=rows[i]["distinct_loci"], oracle_fp=rows[i]["oracle_fp"],
                     is_impure=rows[i]["is_impure"]) for i in card_order[:10]]
    # multiplicity flag: any block with high max_min_mult (repeat-hub survivors)?
    mm_order = np.argsort(-np.array([r["max_min_mult"] for r in rows]), kind="mergesort")
    top_mult = [dict(family_id=rows[i]["family_id"], dominant=rows[i]["dominant_gene"],
                     max_min_mult=rows[i]["max_min_mult"], max_mult=rows[i]["max_mult"],
                     is_impure=rows[i]["is_impure"], oracle_fp=rows[i]["oracle_fp"])
                for i in mm_order[:10]]

    # ---------------- write outputs ----------------
    cols = (["family_id", "dominant_gene", "n_protein_families", "is_impure",
             "oracle_mapped", "oracle_fp", "oracle_cls", "fp_score", "fp_score_oof"]
            + FEATS + ["genes"])
    with open(OUT_TSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)

    out = dict(
        determinism=dict(pythonhashseed=os.environ.get("PYTHONHASHSEED"), seed=SEED,
                         hi_core=HI_CORE, diploid_margin=margin),
        n_families=N, n_ep_impure=int(y.sum()),
        n_oracle_mapped=len(orc_map_idx), n_oracle_fp=n_orc_fp,
        size_confound_ep_rate={k: dict(n=v[0], impure=v[1],
                                       rate=round(v[1] / v[0], 3) if v[0] else None)
                               for k, v in conf.items()},
        per_feature_auc=feat_auc,
        size_baseline_auc=round(size_auc, 3),
        cv_logistic=dict(auc_size_only=round(auc_size, 3),
                         auc_graph_only=round(auc_graph, 3),
                         auc_size_plus_graph=round(auc_full, 3),
                         delta_graph_over_size=round(auc_full - auc_size, 3),
                         note="5-fold random split seed 0, pooled out-of-fold AUC(is_Ep_impure)"),
        triage_precision_at_K=dict(heldout_oof=triage_heldout, in_sample=triage_insample),
        oracle_fp_ranks_heldout=oracle_fp_ranks,
        class_flagging=dict(
            MAGE_block549=brief(mage),
            fam17_block17=brief(fam17),
            GSTM2_blocks=[brief(r) for r in gstm2],
            FOXO1_blocks=[brief(r) for r in foxo1],
            top10_by_cardinality_loci_per_seqclust=top_card,
            top10_by_multiplicity_max_min_mult=top_mult),
    )
    with open(OUT_JSON, "w") as fh:
        json.dump(out, fh, indent=2, default=str)

    # ---------------- console summary ----------------
    P = print
    P("\n==================== BLOCK-GRAPH FP DETECTOR ====================")
    P(f"families={N}  E_p-impure={int(y.sum())}  oracle-mapped={len(orc_map_idx)}  oracle-FP={n_orc_fp}")
    P("\n-- SIZE CONFOUND (E_p-impurity rate by block size) --")
    for k, v in conf.items():
        P(f"  {k:12s}: {v[1]:3d}/{v[0]:3d} = {(v[1]/v[0] if v[0] else 0):.3f}")
    P("\n-- PER-FEATURE AUC(is_Ep_impure)  [raw | lin-size-resid | QUAD-size-resid] --")
    P("   (quad basis removes nonlinear size ratios; artifact features collapse to ~0.5)")
    for f in FEATS:
        fa = feat_auc[f]
        star = " *" if abs(fa["auc_size_resid_quad"] - 0.5) >= 0.08 else "  "
        P(f"  {f:22s} raw={fa['auc_raw']:.3f}  lin={fa['auc_size_resid']:.3f}  "
          f"quad={fa['auc_size_resid_quad']:.3f}{star}")
    P(f"\n  size-only baseline AUC (best size feature)     = {size_auc:.3f}")
    P("\n-- CV LOGISTIC (held-out AUC, 5-fold seed 0) --")
    P(f"  size-only          AUC = {auc_size:.3f}")
    P(f"  graph-only         AUC = {auc_graph:.3f}")
    P(f"  size + graph       AUC = {auc_full:.3f}   (delta over size = {auc_full-auc_size:+.3f})")
    P("\n-- TRIAGE precision@K vs DIPLOID ORACLE (held-out OOF score) --")
    for K in Ks:
        t = triage_heldout[K]
        P(f"  top-{K:2d}: oracle-FP captured {t['oracle_fp_captured']}"
          f"  (recall {t['recall_oracle_fp']}) | E_p-impure {t['ep_impure_in_topK']}")
    P("\n-- ORACLE-FP held-out ranks (of {} families) --".format(N))
    for e in oracle_fp_ranks:
        P(f"  fam {e['family_id']:>4} {e['dominant']:16s} cls={e['cls']:16s} "
          f"dl={e['distinct_loci']:2d}  rank {e['rank_heldout']}/{e['rank_of_N']}  fp_score={e['fp_score']}")
    P("\n-- CLASS FLAGGING --")
    P(f"  MAGE (block 549)  : {brief(mage)}")
    P(f"  fam17 (block 17)  : {brief(fam17)}")
    for r in gstm2:
        P(f"  GSTM2 block {r['family_id']:>4}: is_impure={r['is_impure']} oracle_fp={r['oracle_fp']} "
          f"loci/seqclust={r['loci_per_seqclust']:.2f} fp_score={r['fp_score']} rank={rank_of[rows.index(r)]}")
    P("\n  top-10 by CARDINALITY (loci_per_seqclust):")
    for e in top_card:
        P(f"    fam {e['family_id']:>4} {e['dominant']:16s} lps={e['loci_per_seqclust']:.2f} "
          f"dl={e['distinct_loci']:2d} impure={e['is_impure']} oracle_fp={e['oracle_fp']}")
    P("\n  top-10 by MULTIPLICITY (max_min_mult; the fam17-gate feature):")
    for e in top_mult:
        P(f"    fam {e['family_id']:>4} {e['dominant']:16s} max_min_mult={e['max_min_mult']} "
          f"max_mult={e['max_mult']} impure={e['is_impure']} oracle_fp={e['oracle_fp']}")
    P(f"\n[write] {OUT_TSV}")
    P(f"[write] {OUT_JSON}")
    return out


if __name__ == "__main__":
    main()
