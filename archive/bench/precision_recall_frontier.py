#!/usr/bin/env python
"""PRECISION-RECALL FRONTIER for the RNA multi-copy family DEFINITION.

GOAL
----
Build a HIGH-PRECISION family definition that AVOIDS the residual over-merge false
positives (GSTM2 domain-hub, the fam17 repeat hub, FOXO1, the collapsed-array oversize
blobs) at a controlled RECALL cost, and expose the whole PRECISION-RECALL FRONTIER across
operating points so the operating point can be chosen.  DNA (the phased-diploid CN oracle
bench/diploid_cn_oracle.tsv) is VALIDATION ONLY -- it never gates an RNA edge.

KEY FACTS THIS SCRIPT HONORS (do NOT re-derive -- they are why the levers are what they are)
-------------------------------------------------------------------------------------------
  * Tightening the homology thresholds (core_recip / aln_frac) does NOT remove the residual
    2-cluster over-merges: GSTM2 and the MAGE-class blobs SHARE EXTENSIVE sequence, so they
    pass ANY aln_frac; tightening only kills divergent paralogs => pure RECALL loss.  The
    F1-alnfrac point (core>=0.26 AND aln_frac>=0.72) is included precisely to SHOW this.
  * The precision lever for the 2-cluster over-merges is AGGRESSIVE SPLITTING (stricter
    gamma-quasi-clique cohesion, or a best-2-cut modularity split), NOT a threshold.  The
    GSTM2 domain-hub block is genuinely two weakly-joined clusters (best-2-cut modularity
    ~0.33, top-1% of real-dense families; TP_FP_DIFFERENTIATOR_PANEL.md) so best-2-cut
    splits it.
  * MAGE (MAGEA9 and the dense-uniform collapsed MAGE-family / oversize blobs) is a DENSE
    UNIFORM blob: best-2-cut modularity ~0, edge_density ~1.  NO RNA operator (no threshold,
    no gamma, no best-2-cut) can split it -- only DNA copy-number resolves it.  This is the
    HONEST FLOOR: every RNA operating point below must LEAVE MAGE standing.

OPERATING POINTS (each = a full catalog + its P/R vs the diploid oracle + E_p purity)
-------------------------------------------------------------------------------------
  1. RECALL-PRESERVING (current default) -- the baseline.
     edge KEEP iff core_recip>=0.19 AND aln_frac>=0.24 AND NOT(min_shared_mult>=20)
     -> gamma-quasi-clique(0.20) -> allele-demote.   (== bench/family_rna_refine.py)
  2. F1-alnfrac -- tighter homology (core>=0.26 AND aln_frac>=0.72), else identical.
  3. AGGRESSIVE-GAMMA -- default edge rule, gamma=0.30 and gamma=0.40 (stricter cohesion).
  4. BEST-2-CUT SPLIT -- default catalog (gamma 0.20), then recursively split any family
     whose best-2-cut (spectral Fiedler) modularity Q >= t; sweep t in {0.15..0.35}.
  5. COMBINED high-precision -- default edge rule + aggressive gamma(0.30) + best-2-cut(t*).

Each catalog is scored with the SAME shipped eval machinery the committed catalogs use, so
the numbers line up with bench/family_level_pr_current.json:
    graph_def_refine_sweep.eval_partition   (E_p impurity + retention)
    rna_only_edge_oracle.oracle_residuals   (named diploid-oracle residual-FP roster)
family-level PRECISION  = 1 - distinct_FP_blocks / oracle_mapped_families (gold oracle, dedup)
                        (+ task-formula 1 - flag_instances/mapped, to match the 0.85 headline)
family-level RECALL     = oracle multi_copy genes recovered / 57
E_p BLOCK PURITY        = 1 - impure_blocks / n_families

Deterministic: PYTHONHASHSEED=0 re-exec, gamma seed=0, spectral bisection via numpy.linalg.eigh
(sign-invariant), sorted writes.  Re-runs are byte-identical.
Writes: bench/precision_recall_frontier.tsv + bench/precision_recall_frontier.json
Run:    /home/juanfra/miniforge3/bin/python bench/precision_recall_frontier.py
"""
import os
import sys

if os.environ.get("PYTHONHASHSEED") != "0":
    os.environ["PYTHONHASHSEED"] = "0"
    os.execv(sys.executable, [sys.executable] + sys.argv)

import json
from collections import Counter, defaultdict

import numpy as np
import networkx as nx
from networkx.algorithms.community import modularity as nx_modularity

BENCH = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, BENCH)
import family_er_pr as FP            # noqa: E402
import genome_family_def as G        # noqa: E402
import graph_def_refine_sweep as SW  # noqa: E402
import rna_only_edge_oracle as RN    # noqa: E402
import family_rna_refine as FR       # noqa: E402
import family_level_pr_current as PRC  # noqa: E402

OUT_TSV = os.path.join(BENCH, "precision_recall_frontier.tsv")
OUT_JSON = os.path.join(BENCH, "precision_recall_frontier.json")

SEED = SW.SEED                        # 0
MARGIN = SW.DIPLOID_MARGIN            # 1.5 (oversize_diploid margin above the diploid CN)

# named diploid-DNA-confirmed FP roster (rna_only_edge_oracle.DNA_RESIDUAL_FP), the FP we track
GSTM2 = "GSTM2"
FOXO1 = "FOXO1"
MAGE = "MAGEA9"                       # the MAGE-family oracle representative (dense-uniform floor)


# =========================================================================== builder
def make_kept_edges(ctx, raw_fams, edge_pairs, pair_core, univ_aln, pair_repeat_mult,
                    core_min, aln_min, repeat_gate):
    """Reproduce family_rna_refine.build_catalog's KEEP/CUT edge decision, parametrized by the
    (core_min, aln_min) homology thresholds and the repeat-hub gate.  Returns (all_nodes, kept)."""
    gene_of_dn = ctx["gene_of_dn"]
    all_nodes = set()
    for f in raw_fams:
        all_nodes.update(f)

    def core_aln_keep(k):
        c = pair_core.get(k)
        c = c if c is not None else 0.0
        a = univ_aln.get(k)
        a = a if a is not None else 0.0
        return (c >= core_min) and (a >= aln_min)

    def repeat_hub(k):
        m = pair_repeat_mult.get(k)
        return (m is not None) and (m >= FR.REPEAT_MULT_MIN)

    Gr = nx.Graph()
    Gr.add_nodes_from(all_nodes)
    for a, b in edge_pairs:
        if a in all_nodes and b in all_nodes:
            Gr.add_edge(a, b)
    kept = set()
    for u, v in Gr.edges():
        ga, gb = gene_of_dn.get(u), gene_of_dn.get(v)
        if ga is None or gb is None or ga == gb:
            kept.add(frozenset((u, v)))            # within-gene / unannotated: never an over-merge
            continue
        k = frozenset((ga, gb))
        if core_aln_keep(k) and not (repeat_gate and repeat_hub(k)):
            kept.add(frozenset((u, v)))
    return all_nodes, kept


# --------------------------------------------------------------------------- best-2-cut split
def _induced(block, adj, ew):
    """weighted induced subgraph of `block` from the kept-edge adjacency `adj`; weight = core_recip."""
    sub = nx.Graph()
    sub.add_nodes_from(block)
    bset = set(block)
    for u in block:
        for v in adj.get(u, ()):
            if v in bset and u < v:
                sub.add_edge(u, v, weight=ew.get(frozenset((u, v)), 0.13))
    return sub


def best2cut_modularity(block, adj, ew):
    """Best-2-cut modularity Q of a family block = modularity of the SPECTRAL (Fiedler) bipartition
    of its weighted induced subgraph.  Q ~ 0.3+ => genuinely two weakly-joined clusters (GSTM2-hub);
    Q ~ 0 => dense uniform blob (MAGE).  Deterministic (numpy.linalg.eigh; eigenvector-sign-invariant).
    Returns (Q, sideA, sideB) or (None, None, None) if the block is too small / edgeless to bisect."""
    nodes = sorted(block)
    n = len(nodes)
    if n < 4:
        return None, None, None
    idx = {u: i for i, u in enumerate(nodes)}
    W = np.zeros((n, n))
    for u in nodes:
        for v in adj.get(u, ()):
            if v in idx and idx[u] < idx[v]:
                w = ew.get(frozenset((u, v)), 0.13)
                W[idx[u], idx[v]] = w
                W[idx[v], idx[u]] = w
    if W.sum() == 0:
        return None, None, None
    L = np.diag(W.sum(1)) - W                      # unnormalized Laplacian
    _, vecs = np.linalg.eigh(L)
    fied = vecs[:, 1]                              # Fiedler vector (2nd smallest eigenpair)
    A = [nodes[i] for i in range(n) if fied[i] >= 0]
    B = [nodes[i] for i in range(n) if fied[i] < 0]
    if not A or not B:
        return None, None, None
    sub = _induced(block, adj, ew)
    Q = nx_modularity(sub, [set(A), set(B)], weight="weight")
    return Q, A, B


def best2cut_split(blocks, adj, ew, genes, t):
    """Recursively split any family whose best-2-cut modularity Q >= t; keep only parts with
    >=2 distinct loci (a part isolating a single locus is dropped = the recall cost of over-splitting).
    Returns (new_blocks, n_split_events)."""
    out = []
    n_split = 0

    def rec(block):
        nonlocal n_split
        Q, A, B = best2cut_modularity(block, adj, ew)
        if Q is not None and Q >= t:
            n_split += 1
            rec(A)
            rec(B)
        else:
            if G.distinct_loci(set(block), genes) >= 2:
                out.append(sorted(block))
            # else: dropped (fell below the multi-copy predicate) -- over-split recall cost

    for b in blocks:
        rec(list(b))
    return out, n_split


# --------------------------------------------------------------------------- allele demote (reused logic)
def apply_demote(refined, gene_of_dn, genes, allele):
    """family_rna_refine allele-demote: a same-gene multi-locus family whose dominant gene is a
    balanced diploid het (balanced_frac>=0.90 AND copy_like<=0.10) is ALLELIC -> dropped."""
    def demote_gene(g):
        a = allele.get(g)
        return (a is not None and a["balanced_frac"] >= FR.DEMOTE_BAL_MIN
                and a["copy_like"] <= FR.DEMOTE_COPY_MAX)

    catalog, demoted = [], []
    for b in refined:
        gs = [gene_of_dn[dn] for dn in b if dn in gene_of_dn]
        if gs:
            dom, cnt = Counter(gs).most_common(1)[0]
            homog = cnt / len(gs)
            if demote_gene(dom) and homog >= 0.5 and G.distinct_loci(set(b), genes) >= 2:
                demoted.append(dom)
                continue
        catalog.append(sorted(b))
    return catalog, demoted


def build_point(ctx, raw_fams, edge_pairs, feats, core_min, aln_min, gamma,
                repeat_gate=True, best2cut_t=None):
    """Full parametrized catalog: edge KEEP/CUT (core_min,aln_min,repeat_gate) -> gamma-quasi-clique
    refine(gamma,seed=0) -> optional best-2-cut split(best2cut_t) -> allele-demote.
    Returns (catalog list-of-sorted-DN-lists, bookkeeping dict)."""
    genes = ctx["genes"]; gene_of_dn = ctx["gene_of_dn"]; ew = ctx["ew"]
    all_nodes, kept = make_kept_edges(ctx, raw_fams, edge_pairs, feats["pair_core"],
                                      feats["univ_aln"], feats["repeat_mult"],
                                      core_min, aln_min, repeat_gate)
    comps = SW.components_from_edges(all_nodes, kept)
    refined = G.refine_families(comps, [tuple(e) for e in kept], genes, gamma, SEED)
    n_split = 0
    if best2cut_t is not None:
        adj = defaultdict(set)
        for e in kept:
            u, v = tuple(e)
            adj[u].add(v); adj[v].add(u)
        refined, n_split = best2cut_split(refined, adj, ew, genes, best2cut_t)
    catalog, demoted = apply_demote(refined, gene_of_dn, genes, feats["allele"])
    return catalog, dict(n_kept_edges=len(kept), n_split_events=n_split, n_demoted=len(demoted))


# =========================================================================== scoring
def score_catalog(name, catalog, ctx, mc_oracle, adj_default, ew):
    """Family-level P/R vs the diploid oracle + E_p purity + named-FP roster + recall-cost naming.
    Mirrors bench/family_level_pr_current.py exactly so numbers match the committed catalog."""
    genes = ctx["genes"]; gene_of_dn = ctx["gene_of_dn"]; g2rows = ctx["g2rows"]
    blocks = [set(b) for b in catalog]

    ev = SW.eval_partition(name, blocks, ctx)              # n_families(multicopy), impure_blocks, retention
    fams = SW.filter_multicopy(blocks, genes)
    res = RN.oracle_residuals(fams, gene_of_dn, genes, g2rows, G, margin=MARGIN)

    orac = res["orac_blocks"]
    n_allele, n_over, n_multi = len(res["allele"]), len(res["oversize"]), len(res["multifam"])
    flags = n_allele + n_over + n_multi
    fp_sets = {frozenset(e[0]) for e in res["multifam"] + res["oversize"] + res["allele"]}
    distinct_fp = len(fp_sets)
    distinct_fp_named = sorted("+".join(sorted(s)) for s in fp_sets)
    # retained divergent-paralog (TRUTHBAR) co-membership pairs -- the OVER-SPLIT axis the
    # recovered-somewhere recall metric is BLIND to (breaking a TRUTHBAR pair does NOT drop
    # a recovered gene, it splits two real paralogs into separate families).
    retained_TB = SW.project_pairs(fams, gene_of_dn) & ctx["TB"]
    undersize_counter = Counter("+".join(sorted(e[0])) for e in res["undersize"])
    recovered_mc = set(res["recovered"]) & mc_oracle
    P_task = round(1 - flags / orac, 4) if orac else None
    P_dedup = round(1 - distinct_fp / orac, 4) if orac else None
    R = round(len(recovered_mc) / len(mc_oracle), 4) if mc_oracle else None
    n_fam = ev["n_families"]
    P_Ep = round(1 - ev["impure_blocks"] / n_fam, 4) if n_fam else None

    # ---- named FP roster ----
    multifam_named = ["+".join(e[0]) + f"(dl{e[1]})" for e in
                      sorted(res["multifam"], key=lambda x: (-x[1], x[0][0]))]
    oversize_named = ["+".join(e[0]) + f"(dl{e[1]},dipCN{e[2]},{e[3]}x)" for e in
                      sorted(res["oversize"], key=lambda x: -x[3])]
    allele_named = ["+".join(e[0]) + f"(dl{e[1]})" for e in res["allele"]]
    gstm2_blocks = [e for e in res["multifam"] + res["oversize"] if GSTM2 in e[0]]
    foxo1_blocks = [e for e in res["multifam"] + res["oversize"] if FOXO1 in e[0]]
    mage_over = [e for e in res["oversize"] if MAGE in e[0]]
    mage_recovered = MAGE in res["recovered"]

    # dense-uniform floor: best-2-cut modularity of each FP block (Q<0.15 => graph-invisible floor)
    def fp_block_Q(gene_set):
        target = set(gene_set)
        for b in fams:
            gs = {gene_of_dn.get(dn) for dn in b}
            if target <= gs:
                Q, _, _ = best2cut_modularity(list(b), adj_default, ew)
                return (round(Q, 3) if Q is not None else None), len(b), G.distinct_loci(b, genes)
        return None, None, None

    # ---- E_p-impure repeat-hub (fam17-class) roster ----
    impure = PRC.impure_block_roster(fams, ctx)
    hub_blocks = [b for b in impure if b["n_protein_families"] >= 6]     # repeat-mosaic hubs (fam17-class)

    return dict(
        name=name, n_families=n_fam,
        oracle_mapped=orac, n_multifam=n_multi, n_oversize=n_over, n_allele=n_allele,
        distinct_fp_blocks=distinct_fp, distinct_fp_named=distinct_fp_named,
        _distinct_fp_set=fp_sets, _retained_TB=retained_TB, _undersize_counter=undersize_counter,
        P_oracle_task=P_task, P_oracle_dedup=P_dedup, R_oracle=R,
        recovered_mc=len(recovered_mc), oracle_multicopy_total=len(mc_oracle),
        recovered_genes=sorted(recovered_mc),
        impure_blocks=ev["impure_blocks"], P_Ep=P_Ep,
        truthbar_retention=ev["truthbar_retention"], kept_TRUTHBAR=ev["kept_TRUTHBAR"],
        tp_retention=ev["tp_retention"],
        multifam_named=multifam_named, oversize_named=oversize_named, allele_named=allele_named,
        gstm2_multifam_count=len(gstm2_blocks),
        gstm2_named=["+".join(e[0]) + f"(dl{e[1]})" for e in gstm2_blocks],
        foxo1_count=len(foxo1_blocks),
        mage_oversize_count=len(mage_over), mage_recovered=mage_recovered,
        # over-split (recall cost) roster: DNA CN >> surviving RNA loci
        undersize_named=["+".join(e[0]) + f"(dl{e[1]},dipCN{e[2]},{e[3]}x)" for e in
                         sorted(res["undersize"], key=lambda x: -x[3])],
        n_undersize=len(res["undersize"]),
        n_repeat_hubs=len(hub_blocks),
        _res=res,
    )


# =========================================================================== main
def main():
    print("[load] ctx (graph / projection / truth / diploid oracle) ...", flush=True)
    ctx, raw_fams, edge_pairs = PRC.build_ctx()
    genes = ctx["genes"]; gene_of_dn = ctx["gene_of_dn"]; ew = ctx["ew"]
    mc_oracle = PRC.oracle_multicopy_genes(ctx)

    feats = dict(
        pair_core=RN.load_pair_core(gene_of_dn),
        univ_aln=RN.load_universal_aln(),
        allele=RN.load_allele(),
        repeat_mult=FR.load_repeat_mult(),
    )
    print(f"       oracle multi_copy genes (recall denom) = {len(mc_oracle)}", flush=True)

    # default kept-edge adjacency for best-2-cut modularity of FP blocks (scoring only)
    _, kept_def = make_kept_edges(ctx, raw_fams, edge_pairs, feats["pair_core"], feats["univ_aln"],
                                  feats["repeat_mult"], FR.CORE_MIN, FR.ALN_MIN, True)
    adj_default = defaultdict(set)
    for e in kept_def:
        u, v = tuple(e)
        adj_default[u].add(v); adj_default[v].add(u)

    # ---------- define the operating points ----------
    POINTS = [
        dict(key="1_recall_preserving_DEFAULT", core=FR.CORE_MIN, aln=FR.ALN_MIN, gamma=G.GAMMA,
             repeat_gate=True, best2cut_t=None,
             desc="core>=0.19 AND aln>=0.24, gamma 0.20 (== family_rna_refine DEFAULT)"),
        dict(key="2_f1_alnfrac", core=0.26, aln=0.72, gamma=G.GAMMA, repeat_gate=True, best2cut_t=None,
             desc="TIGHT homology core>=0.26 AND aln_frac>=0.72, gamma 0.20"),
        dict(key="3_aggressive_gamma_0.30", core=FR.CORE_MIN, aln=FR.ALN_MIN, gamma=0.30,
             repeat_gate=True, best2cut_t=None, desc="default edge rule, gamma 0.30"),
        dict(key="3_aggressive_gamma_0.40", core=FR.CORE_MIN, aln=FR.ALN_MIN, gamma=0.40,
             repeat_gate=True, best2cut_t=None, desc="default edge rule, gamma 0.40"),
        dict(key="4_best2cut_t0.15", core=FR.CORE_MIN, aln=FR.ALN_MIN, gamma=G.GAMMA,
             repeat_gate=True, best2cut_t=0.15, desc="default + best-2-cut split t>=0.15"),
        dict(key="4_best2cut_t0.20", core=FR.CORE_MIN, aln=FR.ALN_MIN, gamma=G.GAMMA,
             repeat_gate=True, best2cut_t=0.20, desc="default + best-2-cut split t>=0.20"),
        dict(key="4_best2cut_t0.25", core=FR.CORE_MIN, aln=FR.ALN_MIN, gamma=G.GAMMA,
             repeat_gate=True, best2cut_t=0.25, desc="default + best-2-cut split t>=0.25"),
        dict(key="4_best2cut_t0.30", core=FR.CORE_MIN, aln=FR.ALN_MIN, gamma=G.GAMMA,
             repeat_gate=True, best2cut_t=0.30, desc="default + best-2-cut split t>=0.30"),
        dict(key="4_best2cut_t0.35", core=FR.CORE_MIN, aln=FR.ALN_MIN, gamma=G.GAMMA,
             repeat_gate=True, best2cut_t=0.35, desc="default + best-2-cut split t>=0.35"),
        dict(key="5_combined_g0.30_best2cut_t0.25", core=FR.CORE_MIN, aln=FR.ALN_MIN, gamma=0.30,
             repeat_gate=True, best2cut_t=0.25,
             desc="COMBINED: default edge rule + gamma 0.30 + best-2-cut t>=0.25"),
        dict(key="5_combined_g0.30_best2cut_t0.30", core=FR.CORE_MIN, aln=FR.ALN_MIN, gamma=0.30,
             repeat_gate=True, best2cut_t=0.30,
             desc="COMBINED: default edge rule + gamma 0.30 + best-2-cut t>=0.30"),
        dict(key="5_combined_g0.40_best2cut_t0.15", core=FR.CORE_MIN, aln=FR.ALN_MIN, gamma=0.40,
             repeat_gate=True, best2cut_t=0.15,
             desc="COMBINED: default edge rule + gamma 0.40 + best-2-cut t>=0.15 (max-precision)"),
        dict(key="5_combined_g0.40_best2cut_t0.25", core=FR.CORE_MIN, aln=FR.ALN_MIN, gamma=0.40,
             repeat_gate=True, best2cut_t=0.25,
             desc="COMBINED: default edge rule + gamma 0.40 + best-2-cut t>=0.25"),
    ]

    rows = []
    default_recovered = None
    for pt in POINTS:
        print(f"[build] {pt['key']} :: {pt['desc']}", flush=True)
        cat, book = build_point(ctx, raw_fams, edge_pairs, feats, pt["core"], pt["aln"], pt["gamma"],
                                repeat_gate=pt["repeat_gate"], best2cut_t=pt["best2cut_t"])
        sc = score_catalog(pt["key"], cat, ctx, mc_oracle, adj_default, ew)
        sc.update(desc=pt["desc"], core_min=pt["core"], aln_min=pt["aln"], gamma=pt["gamma"],
                  repeat_gate=pt["repeat_gate"], best2cut_t=pt["best2cut_t"],
                  n_kept_edges=book["n_kept_edges"], n_split_events=book["n_split_events"])
        if pt["key"].startswith("1_"):
            default_recovered = set(sc["recovered_genes"])
        rows.append(sc)

    # ---- recall cost = which oracle multi_copy genes were LOST vs the default ----
    # HONEST CAVEAT: R_oracle is a "recovered-SOMEWHERE" metric (a gene still counts as recovered
    # as long as it survives in ANY multi-copy family). It is BLIND to OVER-SPLITTING: splitting a
    # real family into two multi-copy fragments, or shaving real copies off, does NOT drop R_oracle.
    # We therefore ALSO surface two over-split surrogates that DO move: (i) the undersize tail
    # (DNA diploid-CN >> surviving RNA loci) and (ii) the retained divergent-paralog (TRUTHBAR)
    # co-membership pairs.  Both DEGRADE at every split point and are reported next to R.
    default = rows[0]
    default_orac = default["oracle_mapped"]                # FIXED precision denominator (=48)
    default_fp = default["_distinct_fp_set"]
    default_TB = default["_retained_TB"]
    default_und = default["_undersize_counter"]
    default_ktb = default["kept_TRUTHBAR"]
    for sc in rows:
        lost = sorted(default_recovered - set(sc["recovered_genes"])) if default_recovered else []
        gained = sorted(set(sc["recovered_genes"]) - default_recovered) if default_recovered else []
        sc["lost_oracle_genes_vs_default"] = lost
        sc["gained_oracle_genes_vs_default"] = gained

        # ---- FIXED-DENOMINATOR precision (disclose the MOVING oracle_mapped denominator) ----
        # P_oracle_dedup uses oracle_mapped, which RISES as splitting mints more oracle-mapped
        # families (48 default -> 50 gamma0.40 -> 51 best-2-cut): part of any P gain is pure
        # denominator inflation, not FP removal. Recomputing on the FIXED default denom (48)
        # isolates the REAL FP-removal signal and makes the points directly comparable.
        flags = sc["n_multifam"] + sc["n_oversize"] + sc["n_allele"]
        sc["P_dedup_fixed48"] = round(1 - sc["distinct_fp_blocks"] / default_orac, 4)
        sc["P_task_fixed48"] = round(1 - flags / default_orac, 4)
        sc["oracle_mapped_delta_vs_default"] = sc["oracle_mapped"] - default_orac

        # ---- WHICH distinct FP blocks were removed / added vs default (precision ATTRIBUTION) ----
        sc["distinct_fp_removed_vs_default"] = sorted("+".join(sorted(s))
                                                      for s in (default_fp - sc["_distinct_fp_set"]))
        sc["distinct_fp_added_vs_default"] = sorted("+".join(sorted(s))
                                                    for s in (sc["_distinct_fp_set"] - default_fp))

        # ---- OVER-SPLIT cost (the axis R is blind to) ----
        sc["kept_TRUTHBAR_delta_vs_default"] = sc["kept_TRUTHBAR"] - default_ktb
        sc["truthbar_pairs_cut_vs_default"] = default_ktb - sc["kept_TRUTHBAR"]
        broken = sorted("|".join(sorted(p)) for p in (default_TB - sc["_retained_TB"]))
        sc["truthbar_pairs_broken_named_vs_default"] = broken
        newu = sorted({k for k, c in sc["_undersize_counter"].items() if c > default_und.get(k, 0)})
        sc["undersize_genesets_new_vs_default"] = newu
        sc["undersize_delta_vs_default"] = sc["n_undersize"] - default["n_undersize"]

    # ---- best-2-cut Q DICHOTOMY of the named FP blocks (why best-2-cut splits GSTM2 but not MAGE) ----
    cat_def, _ = build_point(ctx, raw_fams, edge_pairs, feats, FR.CORE_MIN, FR.ALN_MIN, G.GAMMA,
                             repeat_gate=True, best2cut_t=None)
    fams_def = [set(b) for b in cat_def]
    named_fp_targets = [
        ({GSTM2, "LOC115930164", "LOC115930576"}, "GSTM2 domain-hub (multifam)"),
        ({GSTM2, "LOC101129940"}, "GSTM2 triangle (multifam)"),
        ({FOXO1, "LOC115933254"}, "FOXO1 (multifam)"),
        ({"LOC129529978", "LOC129529986"}, "X-array (multifam+oversize)"),
        ({"LOC134758618"}, "LOC134758618 (oversize)"),
        ({MAGE}, "MAGEA9 (MAGE floor)"),
    ]
    best2cut_Q = {}
    for tset, label in named_fp_targets:
        q = None
        for b in fams_def:
            if tset <= {gene_of_dn.get(dn) for dn in b}:
                Qv, _, _ = best2cut_modularity(list(b), adj_default, ew)
                q = dict(n_dn=len(b), distinct_loci=G.distinct_loci(b, genes),
                         best2cut_Q=(round(Qv, 3) if Qv is not None else None),
                         splittable=(Qv is not None and Qv >= 0.30))
                break
        best2cut_Q[label] = q

    # ---- determinism check: rebuild the combined point, compare catalogs ----
    cat_a, _ = build_point(ctx, raw_fams, edge_pairs, feats, FR.CORE_MIN, FR.ALN_MIN, 0.30,
                           repeat_gate=True, best2cut_t=0.30)
    cat_b, _ = build_point(ctx, raw_fams, edge_pairs, feats, FR.CORE_MIN, FR.ALN_MIN, 0.30,
                           repeat_gate=True, best2cut_t=0.30)
    determinism = (sorted(map(tuple, cat_a)) == sorted(map(tuple, cat_b)))
    print(f"[determ] combined-point rebuild byte-identical = {determinism}", flush=True)

    # ---- write TSV ----
    # P_dedup/P_task use the MOVING oracle_mapped denom; P_*_fixed48 use the FIXED default denom (48)
    # and are the comparable numbers.  undersize / truthbar_cut are the OVER-SPLIT surrogates R is blind to.
    cols = ["point", "desc", "core_min", "aln_min", "gamma", "best2cut_t", "n_families",
            "oracle_mapped", "P_oracle_task", "P_oracle_dedup",
            "P_task_fixed48", "P_dedup_fixed48", "R_oracle", "recovered_mc",
            "P_Ep", "impure_blocks", "truthbar_retention",
            "n_multifam", "n_oversize", "distinct_fp_blocks",
            "distinct_fp_removed_vs_default",
            "gstm2_multifam_count", "foxo1_count", "mage_oversize_count", "mage_recovered",
            "n_repeat_hubs", "n_split_events",
            "n_undersize", "undersize_delta_vs_default", "truthbar_pairs_cut_vs_default",
            "n_lost_oracle_genes", "lost_oracle_genes"]
    with open(OUT_TSV, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for sc in rows:
            fh.write("\t".join(str(x) for x in [
                sc["name"], sc["desc"], sc["core_min"], sc["aln_min"], sc["gamma"],
                sc["best2cut_t"], sc["n_families"], sc["oracle_mapped"],
                sc["P_oracle_task"], sc["P_oracle_dedup"],
                sc["P_task_fixed48"], sc["P_dedup_fixed48"], sc["R_oracle"], sc["recovered_mc"],
                sc["P_Ep"], sc["impure_blocks"], sc["truthbar_retention"],
                sc["n_multifam"], sc["n_oversize"], sc["distinct_fp_blocks"],
                ",".join(sc["distinct_fp_removed_vs_default"]) or "-",
                sc["gstm2_multifam_count"], sc["foxo1_count"], sc["mage_oversize_count"],
                sc["mage_recovered"], sc["n_repeat_hubs"], sc["n_split_events"],
                sc["n_undersize"], sc["undersize_delta_vs_default"], sc["truthbar_pairs_cut_vs_default"],
                len(sc["lost_oracle_genes_vs_default"]),
                ",".join(sc["lost_oracle_genes_vs_default"]) or "-",
            ]) + "\n")

    # ---- DOMINANCE: is the COMBINED "max-precision" point actually better than gamma0.40-ALONE? ----
    by_key = {sc["name"]: sc for sc in rows}
    g040 = by_key["3_aggressive_gamma_0.40"]
    comb = by_key["5_combined_g0.40_best2cut_t0.15"]
    b2c15 = by_key["4_best2cut_t0.15"]
    dominance = dict(
        claim=("The COMBINED gamma0.40 + best-2-cut t0.15 point is DOMINATED by gamma0.40-ALONE on "
               "the real metrics: identical distinct_fp, its only P_dedup edge is denominator "
               "inflation, and it OVER-SPLITS strictly more."),
        gamma040_alone=dict(distinct_fp=g040["distinct_fp_blocks"], oracle_mapped=g040["oracle_mapped"],
                            P_dedup_moving=g040["P_oracle_dedup"], P_dedup_fixed48=g040["P_dedup_fixed48"],
                            undersize=g040["n_undersize"], kept_TRUTHBAR=g040["kept_TRUTHBAR"],
                            distinct_fp_removed=g040["distinct_fp_removed_vs_default"]),
        combined_g040_best2cut=dict(distinct_fp=comb["distinct_fp_blocks"], oracle_mapped=comb["oracle_mapped"],
                                    P_dedup_moving=comb["P_oracle_dedup"], P_dedup_fixed48=comb["P_dedup_fixed48"],
                                    undersize=comb["n_undersize"], kept_TRUTHBAR=comb["kept_TRUTHBAR"],
                                    distinct_fp_removed=comb["distinct_fp_removed_vs_default"]),
        best2cut_adds_zero_distinct_fp=(g040["distinct_fp_blocks"] == comb["distinct_fp_blocks"]),
        best2cut_over_split_penalty=dict(
            extra_undersize=comb["n_undersize"] - g040["n_undersize"],
            extra_truthbar_pairs_cut=g040["kept_TRUTHBAR"] - comb["kept_TRUTHBAR"]),
        p_dedup_edge_is_denominator_only=(comb["distinct_fp_blocks"] == g040["distinct_fp_blocks"]
                                          and comb["oracle_mapped"] != g040["oracle_mapped"]),
    )

    # ---- PRECISION ATTRIBUTION: distinct_fp 6->4 is driven by gamma removing OVERSIZE blobs,
    #      best-2-cut adds NO distinct-fp reduction (it only RELABELS the GSTM2 hub) ----
    attribution = dict(
        default_distinct_fp_named=default["distinct_fp_named"],
        gamma040_removed_vs_default=g040["distinct_fp_removed_vs_default"],
        best2cut_only_t015_removed_vs_default=b2c15["distinct_fp_removed_vs_default"],
        combined_removed_vs_default=comb["distinct_fp_removed_vs_default"],
        combined_added_vs_default=comb["distinct_fp_added_vs_default"],
        best2cut_on_top_of_gamma040_gross_removes=sorted(
            set(comb["distinct_fp_removed_vs_default"]) - set(g040["distinct_fp_removed_vs_default"])),
        best2cut_on_top_of_gamma040_NET_distinct_fp_change=(comb["distinct_fp_blocks"]
                                                            - g040["distinct_fp_blocks"]),
        note=("distinct_fp 6->4 comes ENTIRELY from gamma removing the two OVERSIZE collapsed-array "
              "blobs (MPHOSPH8, LOC134758618). best-2-cut then swaps the GSTM2 domain-hub LABEL "
              "(GSTM2+LOC115930164+LOC115930576, gstm2_multifam_count 2->1) for its REMNANT block "
              "(LOC115930164+LOC115930576), which is STILL a multifam FP -- so the NET distinct_fp "
              "change from adding best-2-cut on top of gamma0.40 is 0 (4->4). best-2-cut buys a "
              "cosmetic GSTM2-name split, not a real FP removal."),
    )

    # ---- RESIDUAL FP that SURVIVES at the recommended gamma0.40-ALONE point (corrected floor list) ----
    residual_recommended = dict(
        recommended_point="3_aggressive_gamma_0.40",
        surviving_distinct_fp=g040["distinct_fp_named"],
        interpretation=("At gamma0.40-alone the surviving over-merge FPs are: (1) the GSTM2 domain-hub "
                        "GSTM2+LOC115930164+LOC115930576 (2-cluster, best-2-cut Q=0.335 -- splittable "
                        "only by best-2-cut, but the split leaves the hub REMNANT as an FP so it is not "
                        "a real removal); (2) the GSTM2 triangle GSTM2+LOC101129940 (n=3, Q=None, DNA-only "
                        "floor); (3) the X-array LOC129529978+LOC129529986 (dense-uniform, Q=-0.001, "
                        "edge_density~1 = the actual MAGE-CLASS floor FP present here); (4) FOXO1+LOC115933254 "
                        "(never removed at any point; foxo1_count=1 everywhere). MAGEA9 itself is a 2-node "
                        "block, NOT a counted FP at any operating point (mage_oversize_count=0 throughout)."),
    )

    # ---- write JSON ----
    def clean(sc):
        return {k: v for k, v in sc.items() if not k.startswith("_")}
    summary = dict(
        note=("PRECISION-RECALL FRONTIER for the RNA multi-copy family definition. DNA (diploid CN "
              "oracle) = VALIDATION only. P_oracle_dedup=1-distinct_FP/oracle_mapped (MOVING denom); "
              "P_dedup_fixed48=1-distinct_FP/48 (FIXED denom, comparable across points); "
              "R_oracle=oracle multi_copy genes recovered/57; P_Ep=1-impure/n_families."),
        oracle_multicopy_total=len(mc_oracle),
        determinism_combined_point_byte_identical=determinism,
        disclosure_moving_precision_denominator=(
            "P_oracle_dedup and P_oracle_task both divide by oracle_mapped, which is NOT fixed: it RISES "
            "with splitting (48 default -> 50 gamma0.40 -> 51 best-2-cut/combined) because splitting mints "
            "additional oracle-mapped families. Part of any P gain across points is denominator inflation, "
            "NOT FP removal. P_dedup_fixed48 / P_task_fixed48 recompute on the FIXED default denom (48) and "
            "are the comparable numbers. Worked example: gamma0.40-alone P_dedup 0.920 (distinct_fp 4 / 50) "
            "vs COMBINED 0.9216 (distinct_fp 4 / 51) -- IDENTICAL distinct_fp=4, the 0.0016 edge is purely "
            "50->51 denominator inflation; on the fixed-48 denom BOTH are 0.9167."),
        recall_metric_caveat=(
            "R_oracle is a recovered-SOMEWHERE metric (gene counts as recovered if it survives in ANY "
            "multi-copy family); it is BLIND to OVER-SPLITTING. Every split point holds R_oracle=0.8421 "
            "(48/57) yet DEGRADES the two over-split surrogates: undersize (DNA CN >> surviving RNA loci) "
            "and kept_TRUTHBAR (retained divergent-paralog co-membership pairs). Read n_undersize / "
            "undersize_delta_vs_default and kept_TRUTHBAR / truthbar_pairs_cut_vs_default NEXT TO R."),
        honest_floor=("The DNA-only cardinality floor is the DENSE-UNIFORM blob: the MAGE-A array "
                      "representative MAGEA9 (a 2-node block, best-2-cut Q=None, NOT a counted FP at any "
                      "operating point) and the X-array LOC129529978+LOC129529986 (best-2-cut Q=-0.001, "
                      "edge_density~1) -- the X-array is the dense-uniform FP that ACTUALLY survives as a "
                      "counted over-merge. NO RNA operator (no threshold, no gamma, no best-2-cut) splits "
                      "these; only DNA copy-number resolves the MAGE-class cardinality. mage_recovered "
                      "stays True and mage_oversize_count stays 0 at every point (MAGEA9 is never even "
                      "flagged; it is the recall-side floor, and the X-array is the precision-side floor)."),
        precision_attribution=attribution,
        dominance_analysis=dominance,
        residual_fp_at_recommended=residual_recommended,
        best2cut_Q_of_named_FP_blocks=best2cut_Q,
        gamma_offoracle_caveat=("gamma>=0.27 trips the ZNF716/KRAB-ZNF knife-edge (density 0.261; "
                                "genome_family_def GAMMA note) -> aggressive gamma 0.30/0.40 splits "
                                "divergent KRAB-ZNF paralog families that the diploid oracle (sparse, "
                                "high-CN) cannot see = an OFF-ORACLE recall cost NOT captured by R_oracle. "
                                "gamma 0.20 (default / best-2-cut-only points) preserves them."),
        points=[clean(sc) for sc in rows],
    )
    with open(OUT_JSON, "w") as fh:
        json.dump(summary, fh, indent=1, sort_keys=True,
                  default=lambda x: None if (isinstance(x, float) and x != x) else x)

    # ---- console frontier table ----
    P = print
    P("\n==================== PRECISION-RECALL FRONTIER ====================")
    P("P_ded uses the MOVING oracle_mapped denom; Pfx48 = FIXED-48 denom (comparable);")
    P("R is recovered-SOMEWHERE (blind to over-split) -> read undsz(=undersize) & tbCut next to it.")
    P(f"{'point':<34}{'nFam':>5}{'orac':>5}{'P_ded':>7}{'Pfx48':>7}{'R':>7}{'dFP':>4}"
      f"{'GSTM2':>6}{'FOX':>4}{'MAGE':>5}{'hub':>4}{'undsz':>6}{'tbCut':>6}{'lost':>5}")
    for sc in rows:
        P(f"{sc['name']:<34}{sc['n_families']:>5}{sc['oracle_mapped']:>5}{sc['P_oracle_dedup']:>7}"
          f"{sc['P_dedup_fixed48']:>7}{sc['R_oracle']:>7}{sc['distinct_fp_blocks']:>4}"
          f"{sc['gstm2_multifam_count']:>6}{sc['foxo1_count']:>4}"
          f"{('yes' if sc['mage_recovered'] else 'NO'):>5}{sc['n_repeat_hubs']:>4}"
          f"{sc['n_undersize']:>6}{sc['truthbar_pairs_cut_vs_default']:>6}"
          f"{len(sc['lost_oracle_genes_vs_default']):>5}")
    P("\nlegend: orac=oracle_mapped (MOVING P denom); dFP=distinct FP blocks; GSTM2=multifam blocks carrying GSTM2;")
    P("        FOX=multifam carrying FOXO1 (never 0); MAGE=MAGEA9 still recovered (recall floor -> stays 'yes');")
    P("        hub=E_p-impure repeat-mosaic hubs (fam17-class); undsz=undersize over-split tail;")
    P("        tbCut=divergent-paralog TRUTHBAR pairs cut vs default (OVER-SPLIT cost); lost=oracle genes lost vs default.")
    P("\nper-point residual FP + FULL cost (over-merge FP removed AND over-split cost):")
    for sc in rows:
        P(f"\n  {sc['name']}  [{sc['desc']}]")
        P(f"     P_dedup(moving)={sc['P_oracle_dedup']} on {sc['oracle_mapped']} mapped "
          f"(delta {sc['oracle_mapped_delta_vs_default']:+d})  |  P_dedup_fixed48={sc['P_dedup_fixed48']}  "
          f"R={sc['R_oracle']} ({sc['recovered_mc']}/{sc['oracle_multicopy_total']})  P_Ep={sc['P_Ep']}")
        P(f"     distinct FP blocks={sc['distinct_fp_blocks']}  removed vs default: "
          f"{', '.join(sc['distinct_fp_removed_vs_default']) or 'NONE'}"
          f"{('  added: '+', '.join(sc['distinct_fp_added_vs_default'])) if sc['distinct_fp_added_vs_default'] else ''}")
        P(f"     multifam FP : {sc['multifam_named'] or 'NONE'}")
        P(f"     oversize FP : {sc['oversize_named'] or 'NONE'}")
        P(f"     GSTM2={'name split off hub (1 triangle remains)' if sc['gstm2_multifam_count']==1 else str(sc['gstm2_multifam_count'])+' multifam block(s): '+', '.join(sc['gstm2_named'])}"
          f"   MAGE recovered={sc['mage_recovered']} (counted-FP x{sc['mage_oversize_count']})")
        P(f"     repeat-hubs(fam17-class)={sc['n_repeat_hubs']}  split_events={sc['n_split_events']}")
        P(f"     RECALL COST (on-oracle, recovered-somewhere): lost genes vs default "
          f"({len(sc['lost_oracle_genes_vs_default'])}): {', '.join(sc['lost_oracle_genes_vs_default']) or 'none'}")
        P(f"     OVER-SPLIT COST (R is blind to this): undersize {default['n_undersize']}->{sc['n_undersize']} "
          f"({sc['undersize_delta_vs_default']:+d})  |  TRUTHBAR pairs cut vs default: {sc['truthbar_pairs_cut_vs_default']}")
        if sc["undersize_genesets_new_vs_default"]:
            P(f"        newly over-split gene-sets: {', '.join(sc['undersize_genesets_new_vs_default'][:10])}"
              f"{' ...' if len(sc['undersize_genesets_new_vs_default'])>10 else ''}")
        if sc["truthbar_pairs_broken_named_vs_default"]:
            P(f"        broken divergent-paralog pairs: {', '.join(sc['truthbar_pairs_broken_named_vs_default'][:10])}"
              f"{' ...' if len(sc['truthbar_pairs_broken_named_vs_default'])>10 else ''}")
    P("\n---- DOMINANCE: COMBINED g0.40+best2cut vs gamma0.40-ALONE ----")
    P(f"  gamma0.40-alone      : distinct_fp={g040['distinct_fp_blocks']}  P_dedup(moving)={g040['P_oracle_dedup']} "
      f"P_dedup_fixed48={g040['P_dedup_fixed48']}  undersize={g040['n_undersize']}  kept_TRUTHBAR={g040['kept_TRUTHBAR']}")
    P(f"  combined g0.40+b2c15 : distinct_fp={comb['distinct_fp_blocks']}  P_dedup(moving)={comb['P_oracle_dedup']} "
      f"P_dedup_fixed48={comb['P_dedup_fixed48']}  undersize={comb['n_undersize']}  kept_TRUTHBAR={comb['kept_TRUTHBAR']}")
    P(f"  => best-2-cut adds NET ZERO distinct-fp on top of gamma0.40 (swaps GSTM2-hub LABEL for its "
      f"remnant LOC115930164+LOC115930576, still an FP; distinct_fp 4->4); the P_dedup edge "
      f"{g040['P_oracle_dedup']}->{comb['P_oracle_dedup']} is pure denom inflation "
      f"({g040['oracle_mapped']}->{comb['oracle_mapped']}, fixed48 IDENTICAL {g040['P_dedup_fixed48']}); "
      f"combined OVER-SPLITS more (undersize +{comb['n_undersize']-g040['n_undersize']}, "
      f"TRUTHBAR pairs cut +{(g040['kept_TRUTHBAR']-comb['kept_TRUTHBAR'])}).")
    P("  RECOMMENDED HIGH-PRECISION POINT = gamma0.40-ALONE (dominates combined on the real metrics).")
    P(f"\ndeterminism (combined-point rebuild byte-identical): {determinism}")
    P(f"\nwrote {OUT_TSV}\nwrote {OUT_JSON}")
    return summary


if __name__ == "__main__":
    main()
