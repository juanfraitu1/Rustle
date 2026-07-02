#!/usr/bin/env python
"""Graph-theory refinement sweep for the RNA family E_r graph.

QUESTION
--------
(1) Does the E_r family GRAPH definition (raw transcript-homology components +
    gamma-quasi-clique refinement R) still carry FALSE POSITIVES (over-merges)?
(2) Can MORE GRAPH THEORY -- k-truss / triangle support, edge embeddedness
    (Jaccard neighbourhood overlap), community detection (Louvain / greedy
    modularity), k-core, graph bridges, combined -- cut the over-merge BRIDGE
    edges better than the shipped single gamma threshold, at fixed real-edge
    retention?

HYPOTHESIS
----------
An over-merge edge is a BRIDGE: it joins two dense family blobs but sits in FEW
triangles (low edge embeddedness / low triangle support); a real family is a
DENSE cluster whose edges are highly embedded. So triangle-support (k-truss) and
embeddedness should CUT bridges while KEEPING dense families -- more principled
than one tuned gamma.

GRAPH
-----
Nodes = de-novo loci (DN_<contig>_<start>_<nreads>), edges = transcript-homology
(core_recip >= 0.13), weight = core_recip. bench/denovo_family_edges.tsv.
1130 raw connected components == the raw families; largest = the 728-locus blob
that holds 98.9% of the over-merge (GENUINE) gene-pair mass (see below).

TWO GROUND-TRUTH MODALITIES (honest labelling -- read carefully)
----------------------------------------------------------------
(A) EDGE truth = REFERENCE-cDNA / PROTEIN homology, NOT the DNA assembly. Every
    raw E_r gene pair is labelled from bench/family_er_pr.tsv:
      REAL_cdna := in_dna_loose==1                      (reference cDNA-vs-cDNA
                   minimap2 identity>=0.90 & max-cov>=0.30; see family_def_dna_pr.py)
      TRUTHBAR  := !dna & in_ep & id>=.50 & mincov>=.60 (reciprocal whole-protein
                   homology divergent paralog -> KEEP-ish)
      GENUINE   := else                                 (over-merge -> CUT)  <- the FP
    This is *sequence-homology* truth in the SAME MODALITY as E_r (cDNA/protein,
    reference-derived). It is NON-CIRCULAR w.r.t. the graph OPERATORS (they never
    see cDNA-id, only core_recip), so it is a fair RANKING oracle for the
    operators -- but it is NOT assembly-independent and must NOT be called
    "independent DNA".

(B) The ONLY assembly-independent quantity is asm_hapCN in the phased-DIPLOID DNA
    oracle (bench/diploid_cn_oracle.tsv): the number of distinct HAPLOTYPE loci in
    the DNA assembly (a DIFFERENT gorilla) that carry a full-length copy of the
    family's transcript, per maternal/paternal haplotype. Its family PARTITION,
    however, comes from validated_families.tsv (a backbone-reinforced E_r RNA
    variant) -- so the assembly supplies only the SCALAR copy number, not the
    family boundaries. Coverage is small: only ~50/858 refined families map to an
    oracle gene (~6%). The genuinely assembly-independent over-merge signals are:
       - allele_as_copy : block has >=2 loci but DNA class==single_locus_allele
                          (asm_hapCN==1 on BOTH haplotypes) => the RNA "copies"
                          are ALLELES.  (threshold-free: DNA CN literally = 1)
       - oversize_diploid: block distinct-loci > 1.5 x DIPLOID DNA CN (mat+pat),
                          a margin above the pure-family baseline offset (the pure
                          GSTM2 family already sits at ~2.5x HAPLOID CN, so the
                          naive >1.5x-HAPLOID flag over-counts and is reported only
                          as a caveated upper bound).
    'multifam' (block spans >=2 distinct oracle REPRESENTATIVE genes) is
    RNA-PARTITION-anchored (validated_families), i.e. an RNA-vs-RNA comparison, so
    it is reported as corroboration, NOT as assembly-family evidence.

EVALUATION per refinement
-------------------------
 (a) EDGE, two honest levels:
       * PAIR-CLOSURE (Rand-like): over-merge CUT rate / REAL RETENTION over the
         TRANSITIVE CLOSURE of gene pairs. GENUINE denom is 129,059 pairs, 98.9%
         inside the ONE 728-node megablob -> this is a quadratic pair-counting
         (Rand) metric dominated by how finely the megablob is subdivided; a high
         overmerge_cut is NOT "96.8% of bridge edges removed".
       * GRAPH-EDGE: over-merge CUT / REAL KEPT over the 11,400 REAL DN edges
         (mapped to labels). This is the honest per-edge number (gamma cuts only
         ~32% of GENUINE graph edges while keeping 99.6% of REAL graph edges).
     Plus per-edge SCORE AUC (triangle support / embeddedness / core_recip weight /
     betweenness): does the score separate GENUINE from REAL edges -- the direct
     bridge-hypothesis test.
 (b) BLOCK (PRIMARY honest number): E_p-impure block count / rate (a block with
     >=2 distinct non-mega protein families = an over-merged family).
 (c) FAMILY-vs-DNA (assembly-independent): allele_as_copy, oversize_diploid,
     ratio distribution vs asm_hapCN; multifam reported as RNA-vs-RNA corroboration.

BEATS analysis uses TWO criteria: (i) the TASK criterion (max PAIR overmerge_cut at
tp_retention >= baseline) and (ii) a multi-criterion PARETO test that ADDS
truthbar_retention and impure_blocks (gamma destroys ~85% of TRUTHBAR; some
community operators keep more).

Deterministic: PYTHONHASHSEED=0, seed=0 community detection, sorted iteration/writes.
Outputs: bench/graph_def_refine_sweep.tsv + .json.
Run: /home/juanfra/miniforge3/bin/python bench/graph_def_refine_sweep.py
"""
import json
import math
import os
import statistics
import sys
from collections import defaultdict
from itertools import combinations

if os.environ.get("PYTHONHASHSEED") != "0":
    os.environ["PYTHONHASHSEED"] = "0"
    os.execv(sys.executable, [sys.executable] + sys.argv)

import networkx as nx
from networkx.algorithms.community import louvain_communities, greedy_modularity_communities

BENCH = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, BENCH)
import family_er_pr as FP           # loaders + gene projection (shared with the shipped metric)
import genome_family_def as G       # refine_families / _refine_component / distinct_loci / GAMMA

GAMMA = G.GAMMA                     # 0.20 (shipped)
SEED = 0
DIPLOID_MARGIN = 1.5               # oversize_diploid: distinct_loci > 1.5 * (mat+pat), a margin above the
                                   # pure-family reference/assembly cross-individual offset (GSTM2 ~2.5x HAPLOID)
EDGES_TSV = os.path.join(BENCH, "denovo_family_edges.tsv")
ERPR_TSV = os.path.join(BENCH, "family_er_pr.tsv")
ORACLE_TSV = os.path.join(BENCH, "diploid_cn_oracle.tsv")


# ------------------------------------------------------------------ load graph + projection + truth
def load_weighted_edges():
    ew = {}
    with open(EDGES_TSV) as fh:
        next(fh)
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            ew[frozenset((f[0], f[1]))] = float(f[2])
    return ew


def okf(v, thr):
    return True if v in ("", "NA", "None") else (float(v) >= thr)


def load_pair_labels():
    """gene-pair -> label in {REAL_cdna, TRUTHBAR, GENUINE} for every raw-E_r pair (in_er_raw==1).
       Exact reconstruction of the shipped `class` (verified 0 mismatch).
       NB: REAL_cdna is REFERENCE cDNA-homology (in_dna_loose), NOT the DNA assembly."""
    lab = {}
    with open(ERPR_TSV) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        I = {c: i for i, c in enumerate(h)}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if f[I["in_er_raw"]] != "1":
                continue
            key = frozenset((f[I["gene_a"]], f[I["gene_b"]]))
            if f[I["in_dna_loose"]] == "1":
                lab[key] = "REAL_cdna"
            elif f[I["in_ep"]] == "1" and okf(f[I["id"]], 0.50) and okf(f[I["min_cov"]], 0.60):
                lab[key] = "TRUTHBAR"
            else:
                lab[key] = "GENUINE"
    return lab


def load_oracle():
    """gene -> LIST of oracle rows (a gene may represent >=2 validated families).
       asm_hapCN = HAPLOID assembly CN; diploid = hap_CN_mat + hap_CN_pat (assembly-independent scalars)."""
    g2rows = defaultdict(list)
    with open(ORACLE_TSV) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        I = {c: i for i, c in enumerate(h)}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            g = f[I["gene"]]
            if g in ("NA", "", "None"):
                continue
            try:
                hap = int(f[I["asm_hapCN"]])
                mat = int(f[I["hap_CN_mat"]]); pat = int(f[I["hap_CN_pat"]])
            except ValueError:
                continue
            g2rows[g].append(dict(fam=f[I["family"]], hap=hap, dip=mat + pat, cls=f[I["class"]]))
    return dict(g2rows)


# ------------------------------------------------------------------ helpers on the DN graph
def components_from_edges(all_nodes, kept_edges):
    """partition all_nodes into connected components using ONLY kept_edges (a set of frozenset pairs)."""
    H = nx.Graph()
    H.add_nodes_from(all_nodes)
    for e in kept_edges:
        u, v = tuple(e)
        H.add_edge(u, v)
    return [set(c) for c in nx.connected_components(H)]


def filter_multicopy(blocks, genes):
    """apply the shipped >=2-distinct-loci multi-copy predicate (same object across all operators)."""
    return [b for b in blocks if G.distinct_loci(b, genes) >= 2]


def project_pairs(blocks, gene_of_dn):
    pairs, _ = FP.blocks_to_gene_pairs(blocks, gene_of_dn)
    return pairs


def block_of(blocks):
    bof = {}
    for i, b in enumerate(blocks):
        for n in b:
            bof[n] = i
    return bof


# ------------------------------------------------------------------ evaluation of ONE partition
def eval_partition(name, blocks_full, ctx):
    genes = ctx["genes"]; gene_of_dn = ctx["gene_of_dn"]; lab = ctx["lab"]
    g2f = ctx["g2f"]; mega = ctx["mega"]; g2rows = ctx["g2rows"]; Gr = ctx["G"]
    REAL = ctx["REAL"]; GEN = ctx["GEN"]; TB = ctx["TB"]

    fams = filter_multicopy(blocks_full, genes)          # multi-copy families (the shipped object)
    kept = project_pairs(fams, gene_of_dn)

    # ---- (a1) PAIR-CLOSURE (Rand-like, megablob-dominated) ----
    kR = len(kept & REAL); kG = len(kept & GEN); kT = len(kept & TB)
    tp_ret = kR / len(REAL) if REAL else 0.0
    om_cut = 1 - kG / len(GEN) if GEN else 0.0
    tb_ret = kT / len(TB) if TB else 0.0
    labeled_kept = kR + kG + kT
    genuine_prec = 1 - kG / labeled_kept if labeled_kept else 1.0

    # ---- (a2) GRAPH-EDGE level (the honest per-edge number) ----
    bof = block_of(blocks_full)
    ge_cutG = ge_keptG = ge_cutR = ge_keptR = 0
    for u, v in Gr.edges():
        ga, gb = gene_of_dn.get(u), gene_of_dn.get(v)
        if not ga or not gb or ga == gb:
            continue
        L = lab.get(frozenset((ga, gb)))
        kept_edge = bof.get(u) is not None and bof.get(u) == bof.get(v)
        if L == "GENUINE":
            ge_keptG += kept_edge; ge_cutG += (not kept_edge)
        elif L == "REAL_cdna":
            ge_keptR += kept_edge; ge_cutR += (not kept_edge)
    ge_genuine_cut = ge_cutG / (ge_cutG + ge_keptG) if (ge_cutG + ge_keptG) else 0.0
    ge_real_kept = ge_keptR / (ge_cutR + ge_keptR) if (ge_cutR + ge_keptR) else 0.0

    # ---- (b) BLOCK Ep-impurity (>=2 distinct NON-mega protein families) = over-merged block ----
    impure = 0
    for b in fams:
        nonmega = set()
        for dn in b:
            g = gene_of_dn.get(dn)
            if g and g in g2f and g2f[g] not in mega:
                nonmega.add(g2f[g])
        if len(nonmega) >= 2:
            impure += 1
    n_fam = len(fams)
    impure_rate = impure / n_fam if n_fam else 0.0

    # ---- (c) FAMILY-vs-DIPLOID-DNA-ORACLE (assembly-independent scalars) ----
    orac_blocks = 0
    allele_as_copy = 0; oversize_diploid = 0; oversize_haploid_naive = 0; multifam = 0
    logr_hap = []; logr_dip = []
    recovered = set()
    oversize_examples = []; allele_examples = []
    for b in fams:
        gs = {gene_of_dn[dn] for dn in b if dn in gene_of_dn}
        og = sorted(g for g in gs if g in g2rows)            # oracle REPRESENTATIVE genes present
        if not og:
            continue
        orac_blocks += 1
        recovered.update(og)
        rows = [r for g in og for r in g2rows[g]]
        classes = {r["cls"] for r in rows}
        hap = max(r["hap"] for r in rows)                    # DNA-asserted haploid CN
        dip = max(r["dip"] for r in rows)                    # DNA-asserted diploid CN
        dl = G.distinct_loci(b, genes)
        if len(og) >= 2:                                     # spans >=2 oracle REP genes (RNA-partition-anchored)
            multifam += 1
        if hap > 0:
            logr_hap.append(math.log2(dl / hap))
        if dip > 0:
            logr_dip.append(math.log2(dl / dip))
        # allele_as_copy: RNA multi-loci but DNA says single locus on both haplotypes
        if dl >= 2 and "single_locus_allele" in classes and "multi_copy" not in classes:
            allele_as_copy += 1
            allele_examples.append(dict(genes=og, distinct_loci=dl, asm_hapCN=hap))
        elif "single_locus_allele" not in classes:          # multi_copy-only blocks: size-vs-CN
            if dip > 0 and dl > DIPLOID_MARGIN * dip:
                oversize_diploid += 1
                oversize_examples.append(dict(genes=og, distinct_loci=dl, diploid_CN=dip, asm_hapCN=hap,
                                              ratio_diploid=round(dl / dip, 2)))
            if hap > 0 and dl > 1.5 * hap:                   # naive haploid flag = OVER-COUNTS (caveated)
                oversize_haploid_naive += 1
    oversize_examples.sort(key=lambda e: -e["ratio_diploid"])
    med_hap = round(statistics.median(logr_hap), 3) if logr_hap else 0.0
    med_dip = round(statistics.median(logr_dip), 3) if logr_dip else 0.0

    return dict(
        name=name, n_families=n_fam, kept_pairs=len(kept),
        kept_REAL=kR, kept_GENUINE=kG, kept_TRUTHBAR=kT,
        tp_retention=round(tp_ret, 4), overmerge_cut=round(om_cut, 4),
        truthbar_retention=round(tb_ret, 4), genuine_precision=round(genuine_prec, 4),
        ge_genuine_cut=round(ge_genuine_cut, 4), ge_real_kept=round(ge_real_kept, 4),
        ge_cut_GENUINE=ge_cutG, ge_kept_GENUINE=ge_keptG, ge_cut_REAL=ge_cutR, ge_kept_REAL=ge_keptR,
        impure_blocks=impure, impure_rate=round(impure_rate, 4),
        oracle_blocks=orac_blocks, oracle_recovered_genes=len(recovered),
        oracle_allele_as_copy=allele_as_copy, oracle_oversize_diploid=oversize_diploid,
        oracle_oversize_haploid_naive=oversize_haploid_naive, oracle_multifam=multifam,
        oracle_median_log2_ratio_haploid=med_hap, oracle_median_log2_ratio_diploid=med_dip,
        _oversize_examples=oversize_examples, _allele_examples=allele_examples,
    )


# ------------------------------------------------------------------ per-edge SCORE discriminability (hypothesis)
def edge_score_auc(ctx):
    """For every DN edge whose two loci map to a LABELLED gene pair (REAL vs GENUINE), score the edge by
       triangle support / Jaccard embeddedness / core_recip weight / edge-betweenness (control) and ask:
       does a HIGH score => REAL family edge, LOW => over-merge bridge?  AUC = P(score_REAL > score_GENUINE).
       This is the DIRECT bridge-hypothesis test at the GRAPH-EDGE level (474 REAL, 3136 GENUINE edges)."""
    from scipy.stats import mannwhitneyu
    Gr = ctx["G"]; gene_of_dn = ctx["gene_of_dn"]; lab = ctx["lab"]; ew = ctx["ew"]
    adj = {u: set(Gr[u]) for u in Gr.nodes()}
    print("[edge-score] edge_betweenness (per-component) ...", flush=True)
    # size-confounded control -- compute per component so big-blob edges are not globally inflated
    ebc = {}
    for comp in nx.connected_components(Gr):
        if len(comp) < 3:
            for u, v in Gr.subgraph(comp).edges():
                ebc[frozenset((u, v))] = 0.0
            continue
        sub = Gr.subgraph(comp)
        for (u, v), val in nx.edge_betweenness_centrality(sub, normalized=True).items():
            ebc[frozenset((u, v))] = val

    scores = {"triangle_support": {}, "jaccard_embeddedness": {}, "core_recip_weight": {}, "edge_betweenness_ctrl": {}}
    labels = {}
    tri_by_class = {"REAL": [], "GENUINE": []}
    for u, v in Gr.edges():
        ga, gb = gene_of_dn.get(u), gene_of_dn.get(v)
        if not ga or not gb or ga == gb:
            continue
        key = frozenset((ga, gb))
        L = lab.get(key)
        if L not in ("REAL_cdna", "GENUINE"):
            continue
        Nu = adj[u] - {v}; Nv = adj[v] - {u}
        inter = len(Nu & Nv); union = len(Nu | Nv)
        ekey = frozenset((u, v))
        scores["triangle_support"][ekey] = inter
        scores["jaccard_embeddedness"][ekey] = (inter / union) if union else 0.0
        scores["core_recip_weight"][ekey] = ew.get(ekey, 0.0)
        scores["edge_betweenness_ctrl"][ekey] = ebc.get(ekey, 0.0)
        labels[ekey] = 1 if L == "REAL_cdna" else 0
        tri_by_class["REAL" if L == "REAL_cdna" else "GENUINE"].append(inter)

    out = {}
    pos_keys = [k for k, y in labels.items() if y == 1]
    neg_keys = [k for k, y in labels.items() if y == 0]
    out["n_real_edges"] = len(pos_keys)
    out["n_genuine_edges"] = len(neg_keys)
    for sc, d in scores.items():
        pos = [d[k] for k in pos_keys]
        neg = [d[k] for k in neg_keys]
        if not pos or not neg:
            out[sc] = None
            continue
        U = mannwhitneyu(pos, neg, alternative="two-sided").statistic
        out[sc] = round(U / (len(pos) * len(neg)), 4)
    # concrete triangle-support summary (bridge hypothesis says REAL few, GENUINE many)
    for cl in ("REAL", "GENUINE"):
        t = tri_by_class[cl]
        out[f"triangles_{cl}_median"] = statistics.median(t) if t else None
        out[f"triangles_{cl}_mean"] = round(statistics.fmean(t), 2) if t else None
        out[f"triangles_{cl}_frac_zero"] = round(sum(x == 0 for x in t) / len(t), 3) if t else None
    return out


# ------------------------------------------------------------------ ROC over a threshold sweep
def roc_auc(points):
    """points = list of (tp_loss, overmerge_cut). Add (0,0) keep-all and (1,1) cut-all; trapezoid AUC.
       NB: uses the PAIR-CLOSURE overmerge_cut (Rand-like) -- a ranking device, magnitude megablob-inflated."""
    pts = sorted(set(points) | {(0.0, 0.0), (1.0, 1.0)})
    auc = 0.0
    for (x0, y0), (x1, y1) in zip(pts, pts[1:]):
        auc += (x1 - x0) * (y0 + y1) / 2.0
    return round(auc, 4)


def best_at_retention(rows_with_thr, target_ret):
    """among sweep operating points with tp_retention >= target_ret, the max overmerge_cut (and its thr)."""
    ok = [r for r in rows_with_thr if r["tp_retention"] >= target_ret]
    if not ok:
        return None
    b = max(ok, key=lambda r: r["overmerge_cut"])
    return dict(thr=b["thr"], tp_retention=b["tp_retention"], overmerge_cut=b["overmerge_cut"],
                impure_rate=b["impure_rate"], n_families=b["n_families"])


# ------------------------------------------------------------------ main
def main():
    print("[load] graph / projection / truth / oracle ...", flush=True)
    meta = FP.load_meta(); annot = FP.load_annot(); gene_of = FP.gene_of_factory(annot)
    raw_fams = FP.load_raw_families(); edge_pairs = FP.load_edges()
    genes, gene_of_dn, gene_of_dn_raw, miss, ovstats = FP.build_genes_dict(raw_fams, meta, gene_of)
    ew = load_weighted_edges()
    lab = load_pair_labels()
    g2f, mega = FP.load_prfam()
    g2rows = load_oracle()

    all_nodes = set()
    for f in raw_fams:
        all_nodes.update(f)
    Gr = nx.Graph()
    Gr.add_nodes_from(all_nodes)
    for a, b in edge_pairs:
        if a in all_nodes and b in all_nodes:
            Gr.add_edge(a, b, weight=ew.get(frozenset((a, b)), 0.13))

    REAL = {k for k, v in lab.items() if v == "REAL_cdna"}
    GEN = {k for k, v in lab.items() if v == "GENUINE"}
    TB = {k for k, v in lab.items() if v == "TRUTHBAR"}

    # megablob share of the GENUINE pair mass (justifies "pair-closure = Rand-like, megablob-dominated")
    comps = sorted((set(c) for c in nx.connected_components(Gr)), key=len, reverse=True)
    megablob = comps[0]
    mega_genes = sorted({gene_of_dn[dn] for dn in megablob if dn in gene_of_dn and gene_of_dn[dn]})
    mega_pairs = {frozenset(p) for p in combinations(mega_genes, 2)}
    mega_genuine_frac = round(len(mega_pairs & GEN) / len(GEN), 4) if GEN else 0.0

    print(f"       nodes={Gr.number_of_nodes()} edges={Gr.number_of_edges()} raw_comps={len(raw_fams)} "
          f"labelled_pairs={len(lab)} (REAL_cdna={len(REAL)} GENUINE={len(GEN)} TRUTHBAR={len(TB)}) "
          f"oracle_genes={len(g2rows)} megablob_loci={len(megablob)} "
          f"megablob_GENUINE_pair_frac={mega_genuine_frac}", flush=True)

    ctx = dict(genes=genes, gene_of_dn=gene_of_dn, lab=lab, g2f=g2f, mega=mega, g2rows=g2rows,
               REAL=REAL, GEN=GEN, TB=TB, G=Gr, ew=ew)

    # ---- gamma partition (full, incl singletons) for baseline + combined ops ----
    adj = defaultdict(set)
    fam_of = {}
    for k, m in enumerate(raw_fams):
        for i in m:
            fam_of[i] = k
    for u, v in edge_pairs:
        if u in fam_of and fam_of.get(u) == fam_of.get(v):
            adj[u].add(v); adj[v].add(u)

    def gamma_partition(gamma):
        out = []
        for m in raw_fams:
            for block in G._refine_component(set(m), adj, gamma, louvain_communities, SEED):
                out.append(set(block))
        return out

    def gamma_kept_edges(gamma):
        blocks = gamma_partition(gamma)
        bof = block_of(blocks)
        return {frozenset((u, v)) for u, v in Gr.edges() if bof.get(u) == bof.get(v)}, blocks

    all_edge_set = {frozenset((u, v)) for u, v in Gr.edges()}

    results = []
    sweeps = {}

    # ---- 0. RAW (no refinement) ----
    raw_blocks = [set(f) for f in raw_fams]
    results.append(eval_partition("raw_components", raw_blocks, ctx))

    # ---- 1. BASELINE shipped gamma-quasi-clique(0.20) ----
    gk20, gblocks20 = gamma_kept_edges(GAMMA)
    baseline = eval_partition("gamma_quasiclique_0.20_SHIPPED", gblocks20, ctx)
    results.append(baseline)

    # ---- 1b. PRINCIPLED VARIANT: gamma with a WEIGHT-AWARE splitter ----
    # The shipped _split_once builds an UNWEIGHTED subgraph and calls Louvain without weights, so it
    # discards core_recip -- the ONLY signal that separates over-merge from real edges (edge-AUC 0.846).
    # This variant keeps the EXACT unweighted gamma-quasi-clique STOP certificate (rho_in >= gamma) but
    # makes the SPLIT DIRECTION weighted by core_recip, so it preferentially severs low-homology (bridge)
    # edges. It folds the one informative signal into the operator WITHOUT adding a new tuned threshold.
    def _refine_weighted(nodes, gamma):
        sub = nx.Graph(); sub.add_nodes_from(nodes)
        for u in nodes:
            for v in adj.get(u, ()):
                if v in nodes and u < v:
                    sub.add_edge(u, v, weight=ew.get(frozenset((u, v)), 0.13))
        n = sub.number_of_nodes()
        if n <= 2 or G._induced_density(sub) >= gamma:        # SAME unweighted stop certificate as shipped
            return [set(nodes)]
        parts = None
        for res in (1.0, 2.0, 4.0, 8.0):
            p = louvain_communities(sub, weight="weight", resolution=res, seed=SEED)  # WEIGHTED split
            if len(p) >= 2:
                parts = p; break
        if parts is None:
            if not nx.is_connected(sub):
                parts = [set(c) for c in nx.connected_components(sub)]
            else:
                nl = sorted(sub.nodes()); parts = [set(nl[:n // 2]), set(nl[n // 2:])]
        out = []
        for p in parts:
            p = set(p)
            if len(p) == n:
                return [set(nodes)]
            out.extend(_refine_weighted(p, gamma))
        return out

    def gamma_weighted_partition(gamma):
        out = []
        for m in raw_fams:
            out.extend(_refine_weighted(set(m), gamma))
        return out

    gw_blocks = gamma_weighted_partition(GAMMA)
    results.append(eval_partition("gamma_weightedsplit_0.20_PRINCIPLED", gw_blocks, ctx))

    # gamma sweep (is 0.20 the sweet spot?)
    gamma_rows = []
    for gm in (0.05, 0.10, 0.15, 0.20, 0.26, 0.30, 0.40):
        blk = gamma_partition(gm)
        r = eval_partition(f"gamma_{gm:.2f}", blk, ctx)
        r["thr"] = gm
        gamma_rows.append(r)
    sweeps["gamma"] = gamma_rows

    # ---- 2. k-TRUSS (triangle support >= k-2) ----
    truss_rows = []
    for k in (3, 4, 5, 6, 7):
        H = nx.k_truss(Gr, k)
        kept_edges = {frozenset((u, v)) for u, v in H.edges()}
        blk = components_from_edges(all_nodes, kept_edges)
        r = eval_partition(f"ktruss_k{k}", blk, ctx)
        r["thr"] = k
        truss_rows.append(r)
        results.append(r)
    sweeps["ktruss"] = truss_rows

    # ---- 3. EDGE EMBEDDEDNESS (Jaccard neighbourhood overlap >= t) ----
    adjsets = {u: set(Gr[u]) for u in Gr.nodes()}
    jac = {}
    for u, v in Gr.edges():
        Nu = adjsets[u] - {v}; Nv = adjsets[v] - {u}
        union = len(Nu | Nv)
        jac[frozenset((u, v))] = (len(Nu & Nv) / union) if union else 0.0
    embed_rows = []
    for t in (0.001, 0.01, 0.02, 0.05, 0.08, 0.10, 0.15, 0.20, 0.30, 0.40):
        kept_edges = {e for e, j in jac.items() if j >= t}
        blk = components_from_edges(all_nodes, kept_edges)
        r = eval_partition(f"embed_t{t:.3f}", blk, ctx)
        r["thr"] = t
        embed_rows.append(r)
        results.append(r)
    sweeps["embeddedness"] = embed_rows

    # ---- 4. COMMUNITY DETECTION (Louvain / greedy modularity, weighted) -- EXTENDED resolution ----
    louvain_rows = []
    for res in (0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 16.0, 24.0, 32.0, 64.0, 128.0, 256.0, 512.0):
        parts = louvain_communities(Gr, weight="weight", resolution=res, seed=SEED)
        r = eval_partition(f"louvain_res{res:g}", [set(p) for p in parts], ctx)
        r["thr"] = res
        results.append(r)
        louvain_rows.append(r)
    sweeps["louvain"] = louvain_rows
    greedy_rows = []
    for res in (0.5, 1.0, 2.0, 4.0, 8.0):
        try:
            parts = greedy_modularity_communities(Gr, weight="weight", resolution=res)
            r = eval_partition(f"greedy_mod_res{res:g}", [set(p) for p in parts], ctx)
            r["thr"] = res
            results.append(r)
            greedy_rows.append(r)
        except Exception as e:
            print(f"[greedy] res={res} failed: {e}", flush=True)
    sweeps["greedy"] = greedy_rows

    # ---- 5. k-CORE control (peel low-degree nodes) ----
    kcore_rows = []
    for k in (2, 3, 4, 5, 6):
        core = nx.k_core(Gr, k)
        kept_edges = {frozenset((u, v)) for u, v in core.edges()}
        blk = components_from_edges(all_nodes, kept_edges)
        r = eval_partition(f"kcore_k{k}", blk, ctx)
        r["thr"] = k
        kcore_rows.append(r)
        results.append(r)
    sweeps["kcore"] = kcore_rows

    # ---- 6. GRAPH BRIDGES (2-edge-connected: remove all cut edges) ----
    bridge_edges = set()
    for comp in nx.connected_components(Gr):
        sub = Gr.subgraph(comp)
        if sub.number_of_edges() == 0:
            continue
        for u, v in nx.bridges(sub):
            bridge_edges.add(frozenset((u, v)))
    kept_edges = all_edge_set - bridge_edges
    blk = components_from_edges(all_nodes, kept_edges)
    results.append(eval_partition("bridge_removal_2ec", blk, ctx))

    # ---- 7. COMBINED: gamma AND truss(3); gamma AND embed(t) ----
    truss3_kept = {frozenset((u, v)) for u, v in nx.k_truss(Gr, 3).edges()}
    kept_edges = gk20 & truss3_kept
    blk = components_from_edges(all_nodes, kept_edges)
    results.append(eval_partition("combined_gamma_AND_truss3", blk, ctx))
    for t in (0.02, 0.05, 0.10):
        kept_edges = gk20 & {e for e, j in jac.items() if j >= t}
        blk = components_from_edges(all_nodes, kept_edges)
        results.append(eval_partition(f"combined_gamma_AND_embed{t:.2f}", blk, ctx))

    # ---- per-edge score discriminability (the direct hypothesis test) ----
    escore = edge_score_auc(ctx)
    print(f"[edge-score] AUC(high-score=>REAL): {escore}", flush=True)

    # ---- ROC AUC (PAIR-CLOSURE, ranking device) + best-at-fixed-retention ----
    def sweep_points(rows):
        return [(1 - r["tp_retention"], r["overmerge_cut"]) for r in rows]

    roc = {"_note": "PAIR-CLOSURE (Rand-like) overmerge_cut; magnitude megablob-inflated; ranking device only",
           "ktruss": roc_auc(sweep_points(truss_rows)),
           "embeddedness": roc_auc(sweep_points(embed_rows)),
           "kcore": roc_auc(sweep_points(kcore_rows)),
           "gamma": roc_auc(sweep_points(gamma_rows)),
           "louvain": roc_auc(sweep_points(louvain_rows))}

    base_ret = baseline["tp_retention"]
    base_cut = baseline["overmerge_cut"]
    oppoints = {}
    for target in (base_ret, 0.95, 0.90):
        oppoints[f"ret>={target:.3f}"] = {
            "ktruss": best_at_retention(truss_rows, target),
            "embeddedness": best_at_retention(embed_rows, target),
            "kcore": best_at_retention(kcore_rows, target),
            "gamma": best_at_retention(gamma_rows, target),
            "louvain": best_at_retention(louvain_rows, target),
            "baseline_gamma0.20": dict(thr=0.20, tp_retention=base_ret, overmerge_cut=base_cut,
                                       impure_rate=baseline["impure_rate"], n_families=baseline["n_families"]),
        }

    # ---- BEATS analysis, TWO criteria ----
    # (i) TASK criterion: higher PAIR overmerge_cut at >= baseline retention
    beats_task = []
    for r in results:
        if "tp_retention" not in r:
            continue
        if r["tp_retention"] >= base_ret - 1e-9 and r["overmerge_cut"] > base_cut + 1e-9:
            beats_task.append(dict(name=r["name"], tp_retention=r["tp_retention"], overmerge_cut=r["overmerge_cut"],
                                   impure_rate=r["impure_rate"]))
    beats_task.sort(key=lambda x: -x["overmerge_cut"])

    # (ii) multi-criterion PARETO: dominate gamma on ALL of {tp_ret>=, overmerge_cut>=, truthbar_ret>=, impure<=}
    def dominates(r, b):
        ge = (r["tp_retention"] >= b["tp_retention"] - 1e-9 and
              r["overmerge_cut"] >= b["overmerge_cut"] - 1e-9 and
              r["truthbar_retention"] >= b["truthbar_retention"] - 1e-9 and
              r["impure_blocks"] <= b["impure_blocks"])
        strict = (r["tp_retention"] > b["tp_retention"] + 1e-9 or
                  r["overmerge_cut"] > b["overmerge_cut"] + 1e-9 or
                  r["truthbar_retention"] > b["truthbar_retention"] + 1e-9 or
                  r["impure_blocks"] < b["impure_blocks"])
        return ge and strict
    beats_pareto = [dict(name=r["name"], tp_retention=r["tp_retention"], overmerge_cut=r["overmerge_cut"],
                         truthbar_retention=r["truthbar_retention"], impure_blocks=r["impure_blocks"])
                    for r in results if "tp_retention" in r and dominates(r, baseline)]

    # (iii) GRAPH-EDGE criterion (the task-literal "over-merge edges CUT vs TP edges KEPT", NOT pair-closure):
    #       operators that WEAKLY DOMINATE gamma on {ge_real_kept >=, ge_genuine_cut >, impure_blocks <=}.
    #       This is the honest per-edge test the megablob-inflated pair metric masks.
    beats_graphedge = []
    for r in results:
        if "ge_real_kept" not in r or r["name"] == baseline["name"]:
            continue
        if (r["ge_real_kept"] >= baseline["ge_real_kept"] - 1e-9 and
                r["ge_genuine_cut"] > baseline["ge_genuine_cut"] + 1e-9 and
                r["impure_blocks"] <= baseline["impure_blocks"]):
            beats_graphedge.append(dict(
                name=r["name"], thr=r.get("thr"), tp_retention=r["tp_retention"],
                ge_real_kept=r["ge_real_kept"], ge_genuine_cut=r["ge_genuine_cut"],
                ge_cut_GENUINE=r["ge_cut_GENUINE"], impure_blocks=r["impure_blocks"],
                pair_overmerge_cut=r["overmerge_cut"], truthbar_retention=r["truthbar_retention"]))
    beats_graphedge.sort(key=lambda x: -x["ge_genuine_cut"])

    # 2-of-3 Pareto advantage of community detection over gamma (real-edge, impurity) at fixed retention
    comm_advantage = []
    for r in louvain_rows + greedy_rows:
        adv = []
        if r["kept_REAL"] > baseline["kept_REAL"]:
            adv.append("more_REAL")
        if r["impure_blocks"] < baseline["impure_blocks"]:
            adv.append("fewer_impure")
        if r["truthbar_retention"] > baseline["truthbar_retention"]:
            adv.append("more_TRUTHBAR")
        if len(adv) >= 2 and r["tp_retention"] >= base_ret - 1e-9:
            comm_advantage.append(dict(name=r["name"], thr=r.get("thr"), advantages=adv,
                                       kept_REAL=r["kept_REAL"], impure_blocks=r["impure_blocks"],
                                       truthbar_retention=r["truthbar_retention"],
                                       overmerge_cut=r["overmerge_cut"], tp_retention=r["tp_retention"]))

    # frontier crossing: Louvain op at matched retention with LOWER over-merge-cut (proves gamma on frontier)
    frontier = [dict(name=r["name"], thr=r.get("thr"), tp_retention=r["tp_retention"],
                     overmerge_cut=r["overmerge_cut"])
                for r in louvain_rows if abs(r["tp_retention"] - base_ret) < 1e-9]

    # ---- determinism check ----
    gk20b, gblocks20b = gamma_kept_edges(GAMMA)
    det_gamma = (sorted(map(sorted, gblocks20)) == sorted(map(sorted, gblocks20b)))
    p1 = louvain_communities(Gr, weight="weight", resolution=16.0, seed=SEED)
    p2 = louvain_communities(Gr, weight="weight", resolution=16.0, seed=SEED)
    det_louv = (sorted(map(sorted, (set(x) for x in p1))) == sorted(map(sorted, (set(x) for x in p2))))
    det = det_gamma and det_louv
    print(f"[determ] gamma reproducible={det_gamma}  louvain reproducible={det_louv}", flush=True)

    # ---- write TSV ----
    cols = ["name", "n_families", "kept_pairs", "kept_REAL", "kept_GENUINE", "kept_TRUTHBAR",
            "tp_retention", "overmerge_cut", "truthbar_retention", "genuine_precision",
            "ge_genuine_cut", "ge_real_kept", "impure_blocks", "impure_rate",
            "oracle_blocks", "oracle_recovered_genes", "oracle_allele_as_copy",
            "oracle_oversize_diploid", "oracle_oversize_haploid_naive", "oracle_multifam",
            "oracle_median_log2_ratio_haploid", "oracle_median_log2_ratio_diploid"]
    tsv = os.path.join(BENCH, "graph_def_refine_sweep.tsv")
    with open(tsv, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in sorted(results, key=lambda x: x["name"]):
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")

    def clean(r):
        return {k: v for k, v in r.items() if not k.startswith("_")}

    summary = dict(
        graph=dict(nodes=Gr.number_of_nodes(), edges=Gr.number_of_edges(), raw_components=len(raw_fams),
                   labelled_pairs=len(lab), REAL_cdna=len(REAL), GENUINE=len(GEN), TRUTHBAR=len(TB),
                   oracle_genes=len(g2rows), megablob_loci=len(megablob),
                   megablob_GENUINE_pair_frac=mega_genuine_frac,
                   truth_note="EDGE truth (REAL_cdna/TRUTHBAR/GENUINE) = REFERENCE cDNA/protein homology, "
                              "NOT the DNA assembly; non-circular w.r.t. operators but NOT assembly-independent. "
                              "ONLY asm_hapCN (diploid oracle) is assembly-independent.",
                   oracle_coverage_note=f"only ~{baseline['oracle_blocks']}/{baseline['n_families']} refined "
                                        f"families map to an oracle gene (~"
                                        f"{round(100*baseline['oracle_blocks']/baseline['n_families'],1)}%).",
                   pair_metric_note="overmerge_cut/tp_retention/ROC are PAIR-CLOSURE (Rand-like, quadratic) metrics; "
                                    f"{round(100*mega_genuine_frac,1)}% of GENUINE pair mass is in the ONE "
                                    f"{len(megablob)}-node megablob, so overmerge_cut measures megablob subdivision, "
                                    "NOT fraction of bridge edges removed (see ge_genuine_cut for that)."),
        baseline_gamma_0_20=clean(baseline),
        raw_components=clean(results[0]),
        edge_score_auc=escore,
        roc_auc_pairclosure=roc,
        operating_points_fixed_retention=oppoints,
        beats_shipped_gamma_task_criterion=beats_task,
        beats_shipped_gamma_pareto_4crit=beats_pareto,
        beats_shipped_gamma_graphedge=beats_graphedge,
        community_2of3_pareto_advantage=comm_advantage,
        louvain_frontier_at_matched_retention=frontier,
        residual_fp_examples=dict(
            oversize_diploid=baseline["_oversize_examples"],
            allele_as_copy=baseline["_allele_examples"],
        ),
        sweeps={k: [dict(name=r["name"], thr=r.get("thr"), tp_retention=r["tp_retention"],
                         overmerge_cut=r["overmerge_cut"], ge_genuine_cut=r["ge_genuine_cut"],
                         ge_real_kept=r["ge_real_kept"], truthbar_retention=r["truthbar_retention"],
                         impure_blocks=r["impure_blocks"], impure_rate=r["impure_rate"],
                         n_families=r["n_families"], kept_REAL=r["kept_REAL"],
                         oracle_oversize_diploid=r["oracle_oversize_diploid"],
                         oracle_allele_as_copy=r["oracle_allele_as_copy"],
                         oracle_multifam=r["oracle_multifam"])
                    for r in v] for k, v in sweeps.items()},
        all_results=[clean(r) for r in sorted(results, key=lambda x: x["name"])],
        determinism=dict(gamma=det_gamma, louvain=det_louv, all=det),
    )
    js = os.path.join(BENCH, "graph_def_refine_sweep.json")
    with open(js, "w") as fh:
        json.dump(summary, fh, indent=2, sort_keys=True)

    # ---- console summary ----
    print("\n==== SUMMARY ====")
    print(f"baseline gamma0.20: families={baseline['n_families']} tp_ret={base_ret} "
          f"overmerge_cut(PAIR)={base_cut} ge_genuine_cut(EDGE)={baseline['ge_genuine_cut']} "
          f"ge_real_kept={baseline['ge_real_kept']} impure={baseline['impure_blocks']}/{baseline['n_families']} "
          f"({baseline['impure_rate']}) truthbar_ret={baseline['truthbar_retention']}")
    print(f"megablob: {len(megablob)} loci hold {mega_genuine_frac*100:.1f}% of GENUINE pair mass "
          f"=> PAIR overmerge_cut is Rand-like/megablob-dominated")
    print(f"oracle coverage: {baseline['oracle_blocks']}/{baseline['n_families']} families "
          f"(~{100*baseline['oracle_blocks']/baseline['n_families']:.1f}%) assembly-independent")
    print(f"residual (assembly-independent): allele_as_copy={baseline['oracle_allele_as_copy']} "
          f"oversize_diploid={baseline['oracle_oversize_diploid']} "
          f"(naive-haploid over-count={baseline['oracle_oversize_haploid_naive']}) "
          f"multifam(RNA-anchored)={baseline['oracle_multifam']}")
    print(f"edge-score AUC (high=>REAL): triangle_support={escore['triangle_support']} "
          f"jaccard={escore['jaccard_embeddedness']} core_recip={escore['core_recip_weight']} "
          f"betweenness={escore['edge_betweenness_ctrl']}")
    print(f"  triangles REAL median={escore['triangles_REAL_median']} GENUINE median={escore['triangles_GENUINE_median']}")
    print(f"ROC-AUC (pair-closure) ops: {dict((k,v) for k,v in roc.items() if k!='_note')}")
    print(f"BEATS task-criterion (>=ret & >PAIR-cut): {[b['name'] for b in beats_task] or 'NONE'}")
    print(f"BEATS 4-crit Pareto (dominate gamma incl TRUTHBAR): {[b['name'] for b in beats_pareto] or 'NONE'}")
    print(f"BEATS GRAPH-EDGE (ge_real_kept>= & ge_genuine_cut> & impure<=): "
          f"{[(b['name'],b['ge_genuine_cut'],b['impure_blocks']) for b in beats_graphedge] or 'NONE'}")
    print(f"community 2-of-3 advantage over gamma: {[c['name'] for c in comm_advantage] or 'NONE'}")
    print(f"louvain frontier at matched retention {base_ret}: {frontier}")
    print(f"\nwrote {tsv}\nwrote {js}")


if __name__ == "__main__":
    main()
