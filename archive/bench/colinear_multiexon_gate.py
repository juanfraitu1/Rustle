#!/usr/bin/env python3
"""bench/colinear_multiexon_gate.py -- CHARACTERIZE stage.

Question (user hypothesis)
--------------------------
The shipped E_r edge rule (core_recip>=0.19 AND aln_frac>=0.24) uses aln_frac = the SINGLE
longest shared spliced-exon block / shorter transcript. It does NOT require >=2 shared exons
NOR colinear order, so a substantial shared DOMAIN (one big conserved exon) glues functionally
unrelated genes into a family = the DOMAIN-SHARE FP residual (ABITRAM+CTNNAL1, MOV10+RHOC, ...).

USER FIX: require >=2 shared exons in COLINEAR ORDER. Model of reality = a recent paralog
diverged from a whole shared gene retains >=2 exons in colinear order; a domain-share is ONE
isolated exon in different flanking context. GATE: cut a cross-gene edge iff
  colinear_shared_exons <= K   (sweep K in {1,2}; K=1 == 'require >=2 colinear exons').

THE KILLER TEST (where plain exon-CONTAINMENT failed, cutting real divergent paralogs ~1:1):
do REAL DIVERGENT paralogs retain >=2 COLINEAR shared exons, or erode to a single domain?
If they keep >=2 the gate is recall-safe; if many keep only 1 the gate repeats the containment
failure. Measured, not assumed, on the real TP/truthbar-divergent-paralog set.

colinear_shared_exons
---------------------
For a pair (dn_a, dn_b): reconstruct per-exon spliced sequences (denovo_skeletons introns ->
exon complement; transcript 5'->3' order per strand -- IDENTICAL machinery to
recombination_bridge_detector.family_exons/exon_match_tensor), align every exon of one against
every exon of the other (edlib HW infix, fwd+revcomp, id>=0.70 == ID_THRESH), take each exon's
best target, then LIS = the longest STRICTLY-order-preserving chain of DISTINCT matched target
exons. colinear_shared_exons = max over the two directions (recall-safe). A single isolated
shared domain -> 1; a shared multi-exon gene body -> >=2. (Non-decreasing LIS also reported.)

Substrate = bench/exon_containment_criterion.tsv (5571 MULTI-exon core-present direct E_r edges
with representative loci dn_a/dn_b already selected from ri_sharedlen_universal, the shipped
oracle's LEARN population): 416 TP (DNA-confirmed, near-identical) + 3382 truthbar (E_p
protein-confirmed DIVERGENT real paralogs) + 1773 genuine (over-merge FP == the domain-share
residual). REAL paralogs = TP + truthbar; the DIVERGENT reals = truthbar. Named ep_impure
exemplars + GSTM2/MAGE controls reported separately.

GATE stage (this file, second half)
-----------------------------------
APPLY the gate to the shipped catalog (family_rna_refine.tsv, 573 multi-copy families) and
RE-MEASURE the family-level truths, sweeping K in {0,1,2}. Operationalization (RNA-only, faithful
to the shipped family_level_pr_current machinery):
  * per catalog family: locus subgraph of the shipped E_r edge graph (graph_def_refine_sweep) +
    a SAME-GENE BACKBONE chain (the same-gene guard -> two loci of one gene are NEVER separated);
  * CUT every cross-gene edge (gene_a != gene_b) whose gene-pair colinear_shared_exons <= K;
  * re-derive connected components; DROP components with < 2 distinct loci;
  * re-eval with the SHIPPED eval_partition + oracle_residuals + gene-projection relabel.
Baseline for the gate delta = the SAME operationalization with NO cut (K=None) so the delta is
PURELY the gate (the catalog reference 573/E_p 0.8918/R_oracle 50/57 reproduces to 572/0.8916/50 --
one quasi-clique-only-connected block / 2 loci are lost by re-derivation, NOT by the gate).
Per K we report: (a) DOMAIN-SHARE FP gene-pairs SPLIT (named MOV10+RHOC ...); (b) REAL paralog
pairs wrongly split, DIVERGENT (truthbar / mean-exon-id<0.90) broken out & named; (c) R_oracle
held 50/57?; (d) E_p purity delta; (e) GSTM2/MAGE spared?  A pair is "split" iff co-membered in the
baseline partition and NOT co-membered after the gate (all connecting paths cut) -- the honest
family-level cost, weaker than the raw per-edge cut because same-gene / other-gene paths survive.

Deterministic (PYTHONHASHSEED=0). RNA-only: exon structure + genome-base exon homology; DNA
(in_dna_loose) and E_p classes are VALIDATION labels only, never a gate.
Out: bench/colinear_multiexon_gate.tsv (per pair) + .json (characterize + gate_stage) + stdout.
"""
import os
import sys

if os.environ.get("PYTHONHASHSEED") != "0":
    os.environ["PYTHONHASHSEED"] = "0"
    os.execv(sys.executable, [sys.executable] + sys.argv)

import bisect
import csv
import json
import collections

import pysam

BENCH = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, BENCH)
import recombination_bridge_detector as R   # reuse family_exons / exon_match_tensor / ID_THRESH

CONT_TSV = os.path.join(BENCH, "exon_containment_criterion.tsv")
PR_JSON = os.path.join(BENCH, "family_level_pr_current.json")
OUT_TSV = os.path.join(BENCH, "colinear_multiexon_gate.tsv")
OUT_JSON = os.path.join(BENCH, "colinear_multiexon_gate.json")
CACHE_TSV = os.path.join(BENCH, "colinear_multiexon_cache.tsv")

ID_THRESH = R.ID_THRESH   # 0.70 -- exon-match identity floor for a colinear link

# MAGE / GSTM2 control gene sets (NAMED_FP roster in exon_containment_criterion.py)
GSTM2_GENES = {"GSTM2", "LOC115930164", "LOC115930576", "LOC101129940"}
MAGE_GENES = {"MAGEA9", "LOC129529978", "LOC129529986"}


# --------------------------------------------------------------------------- LIS colinear count
def _lis(seq, strict):
    """length of longest (strictly / non-strictly) increasing subsequence of seq."""
    tails = []
    for x in seq:
        p = bisect.bisect_left(tails, x) if strict else bisect.bisect_right(tails, x)
        if p == len(tails):
            tails.append(x)
        else:
            tails[p] = x
    return len(tails)


def colinear_count(best, b, n, nexb, strict=True, thresh=ID_THRESH):
    """longest order-preserving chain of DISTINCT matched exons of b covered by n (LIS on target j).
       strict=True  -> each target exon used once (distinct shared exons, the honest count).
       strict=False -> non-decreasing (== recombination_bridge_detector.colinear_cov numerator).
       thresh       -> exon-match identity floor (default ID_THRESH=0.70; swept for sensitivity)."""
    pairs = []
    for i in range(nexb):
        idv, j = best.get((b, i, n), (0.0, -1))
        if idv >= thresh and j >= 0:
            pairs.append((i, j))
    if not pairs:
        return 0
    js = [j for _, j in sorted(pairs)]
    return _lis(js, strict)


# identity floors for the recall-vs-threshold sensitivity (0.70 == the specified gate)
SENS_THRESH = [0.50, 0.60, 0.65, 0.70, 0.75, 0.80]


def matched_ids(best, b, n, nexb):
    """list of best identities of b's exons that match some exon of n at id>=ID_THRESH."""
    out = []
    for i in range(nexb):
        idv, j = best.get((b, i, n), (0.0, -1))
        if idv >= ID_THRESH and j >= 0:
            out.append(idv)
    return out


# --------------------------------------------------------------------------- per-pair analysis
def analyze_pair(recs):
    """recs = 2-member exon record list -> colinear metrics."""
    best = R.exon_match_tensor(recs)
    n0 = len(recs[0]["exons"])
    n1 = len(recs[1]["exons"])
    c01 = colinear_count(best, 0, 1, n0, strict=True)
    c10 = colinear_count(best, 1, 0, n1, strict=True)
    d01 = colinear_count(best, 0, 1, n0, strict=False)
    d10 = colinear_count(best, 1, 0, n1, strict=False)
    m01 = matched_ids(best, 0, 1, n0)
    m10 = matched_ids(best, 1, 0, n1)
    n_matched = max(len(m01), len(m10))
    # mean identity over matched exons (divergence proxy for real paralogs); best single exon id
    prim = m01 if len(m01) >= len(m10) else m10
    mean_matched_id = round(sum(prim) / len(prim), 4) if prim else 0.0
    best_exon_id = round(max(m01 + m10), 4) if (m01 or m10) else 0.0
    colinear = max(c01, c10)
    colinear_nd = max(d01, d10)
    n_exon_min = min(n0, n1)
    exon_match_frac = round(n_matched / n_exon_min, 4) if n_exon_min else 0.0
    # colinear count at each identity floor (sensitivity: how much erosion is a 0.70-floor artifact)
    sens = {}
    for t in SENS_THRESH:
        sens[f"colinear_{int(t*100)}"] = max(colinear_count(best, 0, 1, n0, thresh=t),
                                             colinear_count(best, 1, 0, n1, thresh=t))
    return dict(n_exon_a=n0, n_exon_b=n1, n_exon_min=n_exon_min,
                colinear=colinear, colinear_nd=colinear_nd, n_matched=n_matched,
                exon_match_frac=exon_match_frac, mean_matched_id=mean_matched_id,
                best_exon_id=best_exon_id,
                isolated_domain=int(colinear == 1), **sens)


# --------------------------------------------------------------------------- cache
SENS_COLS = [f"colinear_{int(t*100)}" for t in SENS_THRESH]
_CACHE_COLS = (["dn_a", "dn_b", "n_exon_a", "n_exon_b", "n_exon_min", "colinear", "colinear_nd",
                "n_matched", "exon_match_frac", "mean_matched_id", "best_exon_id",
                "isolated_domain"] + SENS_COLS)


def _empty_metrics():
    m = dict(n_exon_a=0, n_exon_b=0, n_exon_min=0, colinear=0, colinear_nd=0, n_matched=0,
             exon_match_frac=0.0, mean_matched_id=0.0, best_exon_id=0.0, isolated_domain=0)
    for c in SENS_COLS:
        m[c] = 0
    return m


def load_cache():
    d = {}
    if not os.path.exists(CACHE_TSV):
        return d
    with open(CACHE_TSV) as f:
        h = f.readline().rstrip("\n").split("\t")
        if h != _CACHE_COLS:          # schema drift -> recompute from scratch
            return {}
        I = {c: i for i, c in enumerate(h)}
        for ln in f:
            t = ln.rstrip("\n").split("\t")
            rec = dict(
                n_exon_a=int(t[I["n_exon_a"]]), n_exon_b=int(t[I["n_exon_b"]]),
                n_exon_min=int(t[I["n_exon_min"]]), colinear=int(t[I["colinear"]]),
                colinear_nd=int(t[I["colinear_nd"]]), n_matched=int(t[I["n_matched"]]),
                exon_match_frac=float(t[I["exon_match_frac"]]),
                mean_matched_id=float(t[I["mean_matched_id"]]),
                best_exon_id=float(t[I["best_exon_id"]]),
                isolated_domain=int(t[I["isolated_domain"]]))
            for c in SENS_COLS:
                rec[c] = int(t[I[c]])
            d[(t[I["dn_a"]], t[I["dn_b"]])] = rec
    return d


def write_cache(cache):
    with open(CACHE_TSV, "w") as out:
        out.write("\t".join(_CACHE_COLS) + "\n")
        for k in sorted(cache):
            v = cache[k]
            row = [k[0], k[1], v["n_exon_a"], v["n_exon_b"], v["n_exon_min"], v["colinear"],
                   v["colinear_nd"], v["n_matched"], v["exon_match_frac"], v["mean_matched_id"],
                   v["best_exon_id"], v["isolated_domain"]] + [v[c] for c in SENS_COLS]
            out.write("\t".join(str(x) for x in row) + "\n")


# --------------------------------------------------------------------------- ep_impure named pairs
def ep_impure_pairs():
    d = json.load(open(PR_JSON))
    roster = d["residual_fp_roster"]["ep_impure_not_dna_named"]
    pairs = set()
    for b in roster:
        pf = b.get("protein_families", {})
        g2pf = {}
        for p, gs in pf.items():
            for g in gs:
                g2pf[g] = p
        named = [g for g in b["genes"] if g in g2pf]
        for i in range(len(named)):
            for j in range(i + 1, len(named)):
                if g2pf[named[i]] != g2pf[named[j]]:
                    pairs.add(frozenset((named[i], named[j])))
    return pairs


# --------------------------------------------------------------------------- distributions
def dist(vals):
    c = collections.Counter(vals)
    return {str(k): c[k] for k in sorted(c)}


def frac_ge2(vals):
    if not vals:
        return None
    return round(sum(1 for v in vals if v >= 2) / len(vals), 4)


def div_bin(mid):
    if mid >= 0.95:
        return "id>=0.95"
    if mid >= 0.90:
        return "0.90-0.95"
    if mid >= 0.80:
        return "0.80-0.90"
    if mid >= 0.70:
        return "0.70-0.80"
    return "<0.70"


# =========================================================================== GATE STAGE
# named exemplar domain-share FPs + control gene sets (task roster)
DOMAIN_SHARE_EXEMPLARS = [("MOV10", "RHOC"), ("RBBP4", "SYNC"), ("RHD", "SDHD"),
                          ("DEDD", "NIT1"), ("KMT2C", "TMEM128"),
                          ("LOC129532935", "SORT1"), ("ABITRAM", "CTNNAL1")]


def _build_block_graph(block, gene_of_dn, Gr, nx):
    """Locus graph for one catalog family: same-gene backbone chain (the guard) + the
       shipped E_r cross-gene edges (tagged with their gene pair)."""
    H = nx.Graph()
    H.add_nodes_from(block)
    byg = collections.defaultdict(list)
    for dn in sorted(block):
        byg[gene_of_dn.get(dn)].append(dn)
    for g, ls in byg.items():
        for i in range(len(ls) - 1):
            H.add_edge(ls[i], ls[i + 1], kind="same")          # same-gene guard
    for u, v in Gr.subgraph(block).edges():
        gu, gv = gene_of_dn.get(u), gene_of_dn.get(v)
        if gu != gv:
            H.add_edge(u, v, kind="cross", gp=frozenset((gu, gv)))
    return H


def _partition(cat_blocks, gene_of_dn, Gr, genes, col, K, nx, G):
    """Apply the gate at cut-threshold K (None == ungated baseline); return multi-copy blocks."""
    out = []
    for b in cat_blocks:
        H = _build_block_graph(b, gene_of_dn, Gr, nx)
        if K is not None:
            rm = [(u, v) for u, v, d in H.edges(data=True)
                  if d["kind"] == "cross" and col.get(d["gp"]) is not None and col[d["gp"]] <= K]
            H.remove_edges_from(rm)
        for comp in nx.connected_components(H):
            if G.distinct_loci(comp, genes) >= 2:
                out.append(set(comp))
    return out


def _gene_families(blocks, gene_of_dn):
    gf = collections.defaultdict(set)
    for i, b in enumerate(blocks):
        for dn in b:
            g = gene_of_dn.get(dn)
            if g:
                gf[g].add(i)
    return gf


def _eval_family(blocks, ctx, C, SW, RN, G, mc_oracle):
    """E_p purity + diploid-oracle recall (with the shipped gene-projection relabel credit)."""
    genes = ctx["genes"]; gene_of_dn = ctx["gene_of_dn"]; g2rows = ctx["g2rows"]
    ev = SW.eval_partition("gate", blocks, ctx)
    ep = round(1 - ev["impure_blocks"] / ev["n_families"], 4) if ev["n_families"] else None
    fams = SW.filter_multicopy(blocks, genes)
    res = RN.oracle_residuals(fams, gene_of_dn, genes, g2rows, G, margin=SW.DIPLOID_MARGIN)
    recov = set(res["recovered"])
    credit, _, _ = C.relabel_recovered_genes(blocks, genes)
    recov_full = (recov | (credit & mc_oracle)) & mc_oracle
    return dict(n_families=ev["n_families"], impure_blocks=ev["impure_blocks"], E_p=ep,
                R_oracle=len(recov_full), recovered_mc=recov_full)


def gate_stage(out_rows):
    """Apply the colinear gate to the shipped catalog, sweep K in {0,1,2}, re-measure the
       family-level truths (R_oracle, E_p) + name the domain-share FP removed and real
       (esp. divergent) paralogs wrongly cut."""
    import networkx as nx
    import family_level_pr_current as C
    import graph_def_refine_sweep as SW
    import rna_only_edge_oracle as RN
    import genome_family_def as G

    # per-gene-pair metrics from the characterize pass (already in memory)
    sub = {}
    for r in out_rows:
        sub[frozenset((r["gene_a"], r["gene_b"]))] = r
    col = {gp: r["colinear"] for gp, r in sub.items()}

    ctx, raw_fams, edge_pairs = C.build_ctx()
    genes = ctx["genes"]; gene_of_dn = ctx["gene_of_dn"]
    Gr = ctx["G"]
    cat_blocks, fam_ids, fam_dom, member_gene = C.load_catalog_blocks()
    mc_oracle = C.oracle_multicopy_genes(ctx)
    n_mc = len(mc_oracle)

    # colinear coverage of the intra-family cross-gene edges the gate can act on
    blkof = {}
    for bi, b in enumerate(cat_blocks):
        for dn in b:
            blkof[dn] = bi
    n_cross = 0; gpset = set(); covered_gp = set()
    for u, v in edge_pairs:
        if u in blkof and v in blkof and blkof[u] == blkof[v]:
            gu, gv = gene_of_dn.get(u), gene_of_dn.get(v)
            if gu and gv and gu != gv:
                n_cross += 1
                gp = frozenset((gu, gv)); gpset.add(gp)
                if gp in col:
                    covered_gp.add(gp)

    # ungated baseline (same operationalization -> isolate the gate delta)
    base_blocks = _partition(cat_blocks, gene_of_dn, Gr, genes, col, None, nx, G)
    base_gf = _gene_families(base_blocks, gene_of_dn)
    base_ev = _eval_family(base_blocks, ctx, C, SW, RN, G, mc_oracle)
    base_R = base_ev["R_oracle"]; base_recov = base_ev["recovered_mc"]

    def comem(gf, ga, gb):
        return bool(gf.get(ga, set()) & gf.get(gb, set()))

    base_co = {gp: comem(base_gf, r["gene_a"], r["gene_b"]) for gp, r in sub.items()}
    n_base_co = sum(base_co.values())

    exemplar_set = {frozenset(p) for p in DOMAIN_SHARE_EXEMPLARS}

    per_K = {}
    for K in (0, 1, 2):
        blocks = _partition(cat_blocks, gene_of_dn, Gr, genes, col, K, nx, G)
        gf = _gene_families(blocks, gene_of_dn)
        ev = _eval_family(blocks, ctx, C, SW, RN, G, mc_oracle)

        # pair-level SPLITS: co-membered in baseline, NOT co-membered after gate
        split = collections.defaultdict(list)
        for gp, r in sub.items():
            if base_co[gp] and not comem(gf, r["gene_a"], r["gene_b"]):
                split[r["cls"]].append(r)
        fp = split["genuine"]; tp = split["TP"]; tb = split["truthbar"]
        ep_fp = [r for r in fp if r["ep_impure"]]
        div_tb = [r for r in tb if r["mean_matched_id"] < 0.90]
        real_split = len(tp) + len(tb)
        total_split = len(fp) + real_split

        # named domain-share FP removed (exemplars that were co-membered and are now split)
        named_fp_removed = []
        for p in DOMAIN_SHARE_EXEMPLARS:
            gp = frozenset(p)
            if gp in sub and base_co[gp] and not comem(gf, *p):
                named_fp_removed.append(f"{p[0]}+{p[1]}(col={sub[gp]['colinear']})")
        exemplars_status = {}
        for p in DOMAIN_SHARE_EXEMPLARS:
            gp = frozenset(p)
            if gp not in sub:
                exemplars_status[f"{p[0]}+{p[1]}"] = "absent"
            elif not base_co[gp]:
                exemplars_status[f"{p[0]}+{p[1]}"] = f"not-comembered(col={sub[gp]['colinear']})"
            else:
                st = "SPLIT" if not comem(gf, *p) else "kept"
                exemplars_status[f"{p[0]}+{p[1]}"] = f"{st}(col={sub[gp]['colinear']},{sub[gp]['cls']})"
        # all ep_impure FP gene-pairs removed (named)
        ep_fp_named = sorted(f"{r['gene_a']}+{r['gene_b']}(col={r['colinear']})" for r in ep_fp)
        # divergent real paralogs wrongly cut (named, sorted by identity)
        div_named = [f"{r['gene_a']}+{r['gene_b']}(id{r['mean_matched_id']:.2f},col{r['colinear']})"
                     for r in sorted(tb, key=lambda r: (r["mean_matched_id"], r["gene_a"], r["gene_b"]))]
        tp_named = sorted(f"{r['gene_a']}+{r['gene_b']}(id{r['mean_matched_id']:.2f},col{r['colinear']})"
                          for r in tp)

        # oracle-gene recall casualties (named)
        lost = sorted(base_recov - ev["recovered_mc"])
        gained = sorted(ev["recovered_mc"] - base_recov)

        # GSTM2 / MAGE spared?  (same # families containing their genes as baseline)
        def ctrl_fams(gf_, gset):
            fams = set()
            for g in gset:
                fams |= gf_.get(g, set())
            present = sorted(g for g in gset if g in gf_)
            return dict(n_families=len(fams), genes_present=present)
        gstm2 = ctrl_fams(gf, GSTM2_GENES)
        mage = ctrl_fams(gf, MAGE_GENES)
        gstm2_base = ctrl_fams(base_gf, GSTM2_GENES)
        mage_base = ctrl_fams(base_gf, MAGE_GENES)
        gstm2_spared = (gstm2["genes_present"] == gstm2_base["genes_present"]
                        and gstm2["n_families"] == gstm2_base["n_families"])
        mage_spared = (mage["genes_present"] == mage_base["genes_present"]
                       and mage["n_families"] == mage_base["n_families"])

        per_K[str(K)] = dict(
            K=K,
            gate=f"cut cross-gene edge iff colinear_shared_exons <= {K} "
                 f"(== require >= {K+1} colinear exons to KEEP)",
            n_families=ev["n_families"], impure_blocks=ev["impure_blocks"],
            E_p=ev["E_p"], E_p_delta_vs_rederive_baseline=round(ev["E_p"] - base_ev["E_p"], 4),
            E_p_delta_vs_catalog_0_8918=round(ev["E_p"] - 0.8918, 4),
            R_oracle=f"{ev['R_oracle']}/{n_mc}",
            R_oracle_held_50_57=(ev["R_oracle"] >= base_R),
            oracle_genes_lost=lost, oracle_genes_gained=gained,
            # (a) domain-share FP removed
            FP_genuine_split=len(fp), FP_ep_impure_split=len(ep_fp),
            named_domain_share_FP_removed=named_fp_removed,
            exemplar_status=exemplars_status,
            ep_impure_FP_removed_named=ep_fp_named,
            # (b) real paralog cost, divergent broken out
            REAL_split_total=real_split, REAL_TP_split=len(tp),
            REAL_truthbar_divergent_split=len(tb),
            REAL_divergent_meanid_lt_0_90_split=len(div_tb),
            split_FP_fraction_pct=(round(100 * len(fp) / total_split, 1) if total_split else None),
            divergent_real_cut_named=div_named,
            TP_real_cut_named=tp_named,
            # (e) controls
            GSTM2_spared=gstm2_spared, MAGE_spared=mage_spared,
            GSTM2=gstm2, MAGE=mage,
        )

    # pick best K: highest E_p among the K that HOLD R_oracle >= baseline; else honest 'none'
    holding = [K for K in (0, 1, 2) if per_K[str(K)]["R_oracle_held_50_57"]]
    if holding:
        best_K = max(holding, key=lambda K: per_K[str(K)]["E_p"])
        best_note = (f"K={best_K}: highest E_p among K that HOLD R_oracle 50/57. "
                     "Note the recall-safe K only cut ZERO-shared-colinear-exon edges; it does NOT "
                     "implement the user's 'require >=2 colinear' proposal (that is K=1, which fails "
                     "the R_oracle bar).")
    else:
        best_K = None
        best_note = "no K holds R_oracle 50/57"

    return dict(
        operationalization=("locus subgraph per catalog family + same-gene backbone (guard); "
                            "cut cross-gene edge iff gene-pair colinear<=K; re-derive components; "
                            "drop <2 distinct loci; re-eval via shipped eval_partition + "
                            "oracle_residuals + gene-projection relabel"),
        catalog_reference=dict(n_families=573, E_p=0.8918, R_oracle="50/57"),
        rederive_baseline_K_none=dict(
            n_families=base_ev["n_families"], E_p=base_ev["E_p"], R_oracle=f"{base_R}/{n_mc}",
            impure_blocks=base_ev["impure_blocks"],
            note="1 quasi-clique-only-connected block / 2 loci lost by re-derivation, NOT by the "
                 "gate; ALL gate deltas are measured against THIS baseline"),
        colinear_coverage=dict(intra_family_crossgene_edges=n_cross,
                               distinct_gene_pairs=len(gpset),
                               covered_by_colinear_substrate=len(covered_gp)),
        substrate_pairs_comembered_in_baseline=n_base_co,
        per_K=per_K,
        best_K=best_K, best_K_note=best_note,
    )


def main():
    skel = R.load_skeletons()
    strand = R.load_strand()
    fa = pysam.FastaFile(R.GENOME)
    meta = {}
    with open(R.META_TSV) as f:
        next(f)
        for line in f:
            t = line.rstrip("\n").split("\t")
            meta[t[0]] = (t[1], int(t[2]), int(t[3]))

    # substrate
    rows = []
    with open(CONT_TSV) as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            rows.append(row)

    ep_pairs = ep_impure_pairs()
    cache = load_cache()
    n_new = 0

    out_rows = []
    for row in rows:
        da, db = row["dn_a"], row["dn_b"]
        key = (da, db)
        if key not in cache:
            if da not in meta or db not in meta:
                cache[key] = _empty_metrics()
            else:
                cha, sa, ea = meta[da]
                chb, sb, eb = meta[db]
                mA = dict(dn=da, gene=row["gene_a"], chrom=cha, start=sa, end=ea)
                mB = dict(dn=db, gene=row["gene_b"], chrom=chb, start=sb, end=eb)
                recs = R.family_exons([mA, mB], skel, strand, fa)
                cache[key] = analyze_pair(recs) if len(recs) >= 2 else _empty_metrics()
            n_new += 1
            if n_new % 500 == 0:
                print(f"[colinear] computed {n_new} new pairs", flush=True)
        m = cache[key]
        gp = frozenset((row["gene_a"], row["gene_b"]))
        rec = dict(gene_a=row["gene_a"], gene_b=row["gene_b"], dn_a=da, dn_b=db,
                   cls=row["cls"], in_dna_loose=int(row["in_dna_loose"]),
                   aln_frac=float(row["aln_frac"]), core_recip=float(row["core_recip"]),
                   ep_impure=int(gp in ep_pairs),
                   gstm2=int(bool(gp & GSTM2_GENES)), mage=int(bool(gp & MAGE_GENES)),
                   **m)
        rec["cut_K1"] = int(rec["colinear"] <= 1)
        rec["cut_K2"] = int(rec["colinear"] <= 2)
        rec["div_bin"] = div_bin(rec["mean_matched_id"])
        out_rows.append(rec)
    write_cache(cache)

    # per-pair TSV
    cols = (["gene_a", "gene_b", "dn_a", "dn_b", "cls", "in_dna_loose", "aln_frac", "core_recip",
             "n_exon_a", "n_exon_b", "n_exon_min", "n_matched", "exon_match_frac",
             "mean_matched_id", "best_exon_id", "colinear", "colinear_nd", "isolated_domain",
             "cut_K1", "cut_K2", "div_bin", "ep_impure", "gstm2", "mage"] + SENS_COLS)
    out_rows.sort(key=lambda r: (r["cls"], r["gene_a"], r["gene_b"], r["dn_a"], r["dn_b"]))
    with open(OUT_TSV, "w") as out:
        out.write("\t".join(cols) + "\n")
        for r in out_rows:
            out.write("\t".join(str(r[c]) for c in cols) + "\n")

    # ---------------- distributions per class -------------------------------------------------
    by_cls = collections.defaultdict(list)
    for r in out_rows:
        by_cls[r["cls"]].append(r)

    def clip(v):
        # bucket colinear counts 0,1,2,3,4,>=5 for a readable histogram
        return str(v) if v <= 5 else "6+"

    class_summary = {}
    for cls, rs in sorted(by_cls.items()):
        cols_v = [r["colinear"] for r in rs]
        class_summary[cls] = dict(
            n=len(rs),
            colinear_hist=dist([clip(v) for v in cols_v]),
            colinear_eq1_pct=round(100 * sum(1 for v in cols_v if v == 1) / len(rs), 2),
            colinear_ge2_pct=round(100 * sum(1 for v in cols_v if v >= 2) / len(rs), 2),
            colinear_ge3_pct=round(100 * sum(1 for v in cols_v if v >= 3) / len(rs), 2),
            median_colinear=sorted(cols_v)[len(cols_v) // 2],
            median_n_exon_min=sorted([r["n_exon_min"] for r in rs])[len(rs) // 2],
        )

    # REAL paralogs = TP + truthbar ; DIVERGENT reals = truthbar
    real = by_cls["TP"] + by_cls["truthbar"]
    fp = by_cls["genuine"]

    # gate outcomes at K in {1,2}
    def gate_block(rs):
        cut1 = [r for r in rs if r["cut_K1"]]
        cut2 = [r for r in rs if r["cut_K2"]]
        return dict(n=len(rs),
                    cut_K1=len(cut1), cut_K1_pct=round(100 * len(cut1) / len(rs), 2) if rs else None,
                    cut_K2=len(cut2), cut_K2_pct=round(100 * len(cut2) / len(rs), 2) if rs else None)
    gate = dict(
        FP_genuine=gate_block(fp),                       # want HIGH cut (precision gain)
        REAL_TP_plus_truthbar=gate_block(real),          # want LOW cut (recall cost)
        REAL_TP=gate_block(by_cls["TP"]),
        REAL_truthbar_DIVERGENT=gate_block(by_cls["truthbar"]),
    )

    # DECISIVE: divergent-real-paralog colinear distribution (truthbar), stratified by identity
    tb = by_cls["truthbar"]
    tb_by_id = collections.defaultdict(list)
    for r in tb:
        tb_by_id[r["div_bin"]].append(r)
    divergent_stratified = {}
    order = ["id>=0.95", "0.90-0.95", "0.80-0.90", "0.70-0.80", "<0.70"]
    for b in order:
        rs = tb_by_id.get(b, [])
        if not rs:
            continue
        cv = [r["colinear"] for r in rs]
        divergent_stratified[b] = dict(
            n=len(rs),
            colinear_ge2_pct=round(100 * sum(1 for v in cv if v >= 2) / len(rs), 2),
            colinear_eq1_pct=round(100 * sum(1 for v in cv if v == 1) / len(rs), 2),
            median_colinear=sorted(cv)[len(cv) // 2],
            median_n_exon_min=sorted([r["n_exon_min"] for r in rs])[len(rs) // 2],
        )

    # the single decisive number: divergent real paralogs (mean_matched_id<0.90) retaining >=2
    div_real = [r for r in real if r["mean_matched_id"] < 0.90]
    div_real_ge2 = round(100 * sum(1 for r in div_real if r["colinear"] >= 2) / len(div_real), 2) if div_real else None
    div_real_lt80 = [r for r in real if r["mean_matched_id"] < 0.80]
    div_real_lt80_ge2 = round(100 * sum(1 for r in div_real_lt80 if r["colinear"] >= 2) / len(div_real_lt80), 2) if div_real_lt80 else None

    # TP paralogs cut at K1 (recall casualties) -- inspect: are they DNA-loose-impure segment shares?
    tp_cut_K1 = sorted([dict(gene_a=r["gene_a"], gene_b=r["gene_b"], colinear=r["colinear"],
                             n_exon_min=r["n_exon_min"], aln_frac=r["aln_frac"],
                             best_exon_id=r["best_exon_id"], mean_matched_id=r["mean_matched_id"])
                        for r in by_cls["TP"] if r["cut_K1"]],
                       key=lambda x: (x["gene_a"], x["gene_b"]))
    truthbar_cut_K1 = [r for r in by_cls["truthbar"] if r["cut_K1"]]

    # named ep_impure exemplars (from task): report cls + colinear
    exemplars = [("MOV10", "RHOC"), ("RBBP4", "SYNC"), ("RHD", "SDHD"), ("DEDD", "NIT1"),
                 ("KMT2C", "TMEM128"), ("LOC129532935", "SORT1"), ("ABITRAM", "CTNNAL1")]
    by_gp = {}
    for r in out_rows:
        by_gp.setdefault(frozenset((r["gene_a"], r["gene_b"])), []).append(r)
    exemplar_out = []
    for ga, gb in exemplars:
        rs = by_gp.get(frozenset((ga, gb)), [])
        if rs:
            r = min(rs, key=lambda x: x["colinear"])  # strongest-cut representative
            exemplar_out.append(dict(pair=f"{ga}+{gb}", cls=r["cls"], colinear=r["colinear"],
                                     aln_frac=r["aln_frac"], core_recip=r["core_recip"],
                                     n_exon_min=r["n_exon_min"], cut_K1=r["cut_K1"],
                                     in_dna_loose=r["in_dna_loose"]))
        else:
            exemplar_out.append(dict(pair=f"{ga}+{gb}", note="no direct multi-exon E_r edge in substrate"))

    # ep_impure edges present in substrate: colinear + cut rates (roster is IMPURE -> report cls mix)
    ep_rows = [r for r in out_rows if r["ep_impure"]]
    ep_summary = dict(
        n=len(ep_rows),
        cls_mix=dist([r["cls"] for r in ep_rows]),
        colinear_eq1_pct=round(100 * sum(1 for r in ep_rows if r["colinear"] == 1) / len(ep_rows), 2) if ep_rows else None,
        cut_K1_pct=round(100 * sum(r["cut_K1"] for r in ep_rows) / len(ep_rows), 2) if ep_rows else None,
        note="ep_impure roster is IMPURE: includes DNA-confirmed TP pairs (e.g. ABCD1/LOC...); "
             "cls_mix shows the genuine-FP vs TP/truthbar split",
    )

    # controls: GSTM2 / MAGE (share many exons -> should NOT be cut)
    def ctrl(flag):
        rs = [r for r in out_rows if r[flag]]
        cv = [r["colinear"] for r in rs]
        nem = [r["n_exon_min"] for r in rs]
        return dict(n=len(rs), cls_mix=dist([r["cls"] for r in rs]),
                    median_colinear=sorted(cv)[len(cv) // 2] if cv else None,
                    colinear_ge2_pct=round(100 * sum(1 for v in cv if v >= 2) / len(rs), 2) if rs else None,
                    cut_K1_pct=round(100 * sum(r["cut_K1"] for r in rs) / len(rs), 2) if rs else None,
                    median_n_exon_min=sorted(nem)[len(nem) // 2] if nem else None)
    controls = dict(GSTM2=ctrl("gstm2"), MAGE=ctrl("mage"))

    # R_oracle: the shipped rule keeps aln_frac>=0.24 & core>=0.19; how many FP does colinear cut
    # AMONG the ones the shipped rule KEEPS (the residual)? and does it hold real recall?
    shipped_keep = lambda r: (r["aln_frac"] >= 0.24 and r["core_recip"] >= 0.19)
    fp_kept = [r for r in fp if shipped_keep(r)]
    real_kept = [r for r in real if shipped_keep(r)]
    residual = dict(
        genuine_FP_shipped_keeps=len(fp_kept),
        of_those_colinear_cut_K1=sum(r["cut_K1"] for r in fp_kept),
        of_those_colinear_cut_K1_pct=round(100 * sum(r["cut_K1"] for r in fp_kept) / len(fp_kept), 2) if fp_kept else None,
        real_shipped_keeps=len(real_kept),
        real_of_those_colinear_cut_K1=sum(r["cut_K1"] for r in real_kept),
        real_of_those_colinear_cut_K1_pct=round(100 * sum(r["cut_K1"] for r in real_kept) / len(real_kept), 2) if real_kept else None,
    )

    # identity-threshold sensitivity: is the id>=0.70 erosion a NUCLEOTIDE-FLOOR artifact?
    # (divergent paralog colinear exons are real but sit at 0.55-0.70 nt id -> hidden at 0.70)
    sens_summary = {}
    for cls, rs in sorted(by_cls.items()):
        row = {}
        for t, col in zip(SENS_THRESH, SENS_COLS):
            cv = [r[col] for r in rs]
            row[f"id>={t:.2f}"] = round(100 * sum(1 for v in cv if v >= 2) / len(rs), 2) if rs else None
        sens_summary[cls] = dict(n=len(rs), colinear_ge2_pct_at_identity_floor=row)

    summary = dict(
        substrate=dict(n_pairs=len(out_rows), class_counts=dist([r["cls"] for r in out_rows]),
                       note="MULTI-exon core-present direct E_r edges; TP=DNA-confirmed near-identical, "
                            "truthbar=E_p-protein-confirmed DIVERGENT real paralogs, genuine=over-merge FP"),
        class_colinear_distribution=class_summary,
        gate_outcomes_by_K=gate,
        DECISIVE_divergent_real_paralogs=dict(
            definition="truthbar = E_p-confirmed real paralogs DNA did NOT confirm (the divergent reals)",
            truthbar_colinear_ge2_pct=class_summary["truthbar"]["colinear_ge2_pct"],
            truthbar_colinear_eq1_pct=class_summary["truthbar"]["colinear_eq1_pct"],
            truthbar_stratified_by_exon_identity=divergent_stratified,
            real_paralogs_meanid_lt_0_90_n=len(div_real),
            real_paralogs_meanid_lt_0_90_colinear_ge2_pct=div_real_ge2,
            real_paralogs_meanid_lt_0_80_n=len(div_real_lt80),
            real_paralogs_meanid_lt_0_80_colinear_ge2_pct=div_real_lt80_ge2,
        ),
        recall_casualties=dict(
            TP_cut_at_K1_n=len(tp_cut_K1),
            TP_cut_at_K1=tp_cut_K1,
            truthbar_cut_at_K1_n=len(truthbar_cut_K1),
        ),
        identity_floor_sensitivity=dict(
            note="colinear>=2 %% per class as the exon-match identity floor varies; 0.70 == the "
                 "specified gate. If divergent-real recovery jumps as the floor drops, the 0.70 "
                 "'erosion to a single domain' is a NUCLEOTIDE-FLOOR artifact, not lost structure.",
            per_class=sens_summary),
        named_exemplars=exemplar_out,
        ep_impure_edges=ep_summary,
        controls_GSTM2_MAGE=controls,
        residual_vs_shipped_rule=residual,
        params=dict(ID_THRESH=ID_THRESH, MIN_EXON=R.MIN_EXON,
                    PYTHONHASHSEED=os.environ.get("PYTHONHASHSEED")),
        out_tsv=OUT_TSV,
    )

    # ---------------- GATE STAGE: apply gate to the catalog, re-measure R_oracle / E_p ----------
    print("[gate] applying colinear gate to family_rna_refine.tsv (K in 0,1,2) ...", flush=True)
    summary["gate_stage"] = gate_stage(out_rows)

    with open(OUT_JSON, "w") as jf:
        json.dump(summary, jf, indent=1, sort_keys=True)
        jf.write("\n")
    print(json.dumps(summary, indent=1, sort_keys=True))


if __name__ == "__main__":
    main()
