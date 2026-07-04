#!/usr/bin/env python
"""BYTE-PARITY golden-fixture generator for the Rust port of the O1 recombinant-SPLIT
gate `bench/recombinant_split.py` -> `vg_family/recombinant_split.rs`
(the O1 over-merge-gate migration, RETIREMENT_AND_MIGRATION.md Tier-2 #6, full-parity path).

It imports the SHIPPED `recombinant_split` (RS) + the shared exon-graph library
`recombination_bridge_detector` (RBD) and reconstructs the EXACT (refined, genes,
gene_of_dn2) that `family_rna_refine.build_catalog` feeds to `RS.split_families`
(the same reconstruction prefix as `gen_refine_families_fixture.py`, proven to
reproduce the shipped catalog).  For each selected family block it runs the REAL
`RS.split_block` and captures a SELF-CONTAINED record (block + the genes/gene_of_dn/
strand/skel subset + the RAW case-carrying genome slices via a mock fetcher, so no
3.4G genome is needed at test time) -> the Python split_block PARTITION.

Two parts:
  (a) "articulation_graphs": >= 20 hand-built graphs (chains, cycles, trees, star,
      shared-vertex triangles, barbell, K4, multi-component, single node, bridge
      edges) -> the networkx articulation-point SET.  PROVES the Rust Tarjan
      cut-vertex SET == nx.articulation_points EXACTLY.
  (b) "blocks": the recall-safe SPLIT case (fid-210 GALNT17|LOC101126070) + fid-187
      + representative NON-split blocks + <3 / >MAX_LOCI skips + colocated-chimera
      skips + NA-dominant skips + a multi-articulation-point family -> the REAL
      RS.split_block partition for each.

Run: PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/gen_recombinant_split_fixture.py
Deterministic (sorted selection, PYTHONHASHSEED=0).
"""
import os
import sys

# determinism: pin the hash seed BEFORE the heavy imports (re-exec preserves argv)
if os.environ.get("PYTHONHASHSEED") != "0":
    os.environ["PYTHONHASHSEED"] = "0"
    os.execv(sys.executable, [sys.executable] + sys.argv)

os.chdir(os.path.dirname(os.path.abspath(__file__)) + "/..")   # repo root (RBD uses relative paths)

import collections
import json

sys.path.insert(0, "bench")
import family_er_pr as FP              # shipped graph context / genes dict
import genome_family_def as G          # refinement operator R
import graph_def_refine_sweep as SW    # components_from_edges + SEED
import rna_only_edge_oracle as RO      # RNA feature loaders
import family_rna_refine as FR         # DEFAULT RNA operating point (gates/constants)
import recombination_bridge_detector as RBD   # shared exon-graph library (ground truth)
import recombinant_split as RS         # the port target (ground truth)
import pysam
import networkx as nx

OUT = os.path.join("src", "rustle", "vg_family", "testdata", "recombinant_split_fixture.json")


# ------------------------------------------------------------------ reconstruction (== build_catalog prefix)
def reconstruct():
    """EXACT copy of the family_rna_refine.build_catalog prefix that produces the
    (refined, genes, gene_of_dn2) fed to RS.split_families -- default gates."""
    FR.rna_only_guard()
    gene_of_dn = RO.load_gene_of_dn()
    pair_core = RO.load_pair_core(gene_of_dn)
    univ_aln = RO.load_universal_aln()
    pair_repeat_mult = FR.load_repeat_mult()

    meta = FP.load_meta()
    annot = FP.load_annot()
    gene_of = FP.gene_of_factory(annot)
    raw_fams = FP.load_raw_families()
    edge_pairs = FP.load_edges()
    genes, gene_of_dn2, *_ = FP.build_genes_dict(raw_fams, meta, gene_of)
    all_nodes = set()
    for f in raw_fams:
        all_nodes.update(f)

    def core_aln_keep(k):
        c = pair_core.get(k)
        c = c if c is not None else 0.0
        a = univ_aln.get(k)
        a = a if a is not None else 0.0
        return (c >= FR.CORE_MIN) and (a >= FR.ALN_MIN)

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
        ga, gb = gene_of_dn2.get(u), gene_of_dn2.get(v)
        if ga is None or gb is None or ga == gb:
            kept.add(frozenset((u, v)))
            continue
        k = frozenset((ga, gb))
        if core_aln_keep(k) and not repeat_hub(k):
            kept.add(frozenset((u, v)))

    comps = SW.components_from_edges(all_nodes, kept)
    refined = G.refine_families(comps, [tuple(e) for e in kept], genes, FR.GAMMA, FR.SEED)
    return refined, genes, gene_of_dn2


# ------------------------------------------------------------------ trace (classify a block's split path)
def trace(block, genes, gene_of_dn, skel, strand, fa):
    """Mirror split_block + _recombinant_bridges guards to classify the block under
    split_block's OWN coord derivation (via genes).  Returns (tag, n_art, reasons)."""
    blk = list(block)
    if len(blk) < 3:
        return ("skip_block_lt3", 0, [])
    members = RS._members(blk, genes, gene_of_dn)
    if len(members) < 3:
        return ("skip_members_lt3", 0, [])
    if len(members) > RS.MAX_LOCI:
        return ("skip_gt_maxloci", 0, [])
    recs = RBD.family_exons(members, skel, strand, fa)
    if len(recs) < 3:
        return ("skip_recs_lt3", 0, [])
    best = RBD.exon_match_tensor(recs)
    Gr, _ = RBD.build_graph(recs, best)
    if Gr.number_of_edges() == 0:
        return ("no_edges", 0, [])
    arts = sorted(set(nx.articulation_points(Gr)))
    n_art = len(arts)
    reasons = []
    n_pass = 0
    for b in arts:
        Gb = Gr.copy(); Gb.remove_node(b)
        nbrs = set(Gr.neighbors(b))
        bridged = sorted((sorted(c) for c in nx.connected_components(Gb) if set(c) & nbrs),
                         key=lambda c: c[0])
        if len(bridged) < 2:
            reasons.append("lt2_bridged"); continue
        cid_of = {m: ci for ci, comp in enumerate(bridged) for m in comp}
        single = collections.Counter()
        for i in range(len(recs[b]["exons"])):
            hit = set()
            for m in nbrs:
                idv, _ = best.get((b, i, m), (0.0, -1))
                if idv >= RBD.ID_THRESH and m in cid_of:
                    hit.add(cid_of[m])
            if len(hit) == 1:
                single[next(iter(hit))] += 1
        if len(single) < 2:
            reasons.append("lt2_single"); continue
        best_single = max((RBD.colinear_cov(recs, best, b, m) for m in nbrs), default=0.0)
        if best_single >= RBD.COLINEAR_COV:
            reasons.append("colinear_ge"); continue
        if RBD.subfams_colocated(recs, bridged):
            reasons.append("colocated"); continue
        involved = sorted(single)
        minor = min(single[c] for c in involved)
        genes_inv = [RBD.subfam_dom_gene(recs, bridged[c]) for c in involved]
        both = all(g not in ("NA", "") for g in genes_inv)
        if not (both and minor >= 2):
            reasons.append("low_conf"); continue
        reasons.append("PASS"); n_pass += 1
    tag = "split" if n_pass else "nonsplit_graph"
    return (tag, n_art, reasons)


# ------------------------------------------------------------------ HW-DP cost proxy (keep the Rust test snappy)
def dp_cost(block, genes, gene_of_dn, skel, strand, fa):
    """Estimate the exon_match_tensor HW-DP cell-ops for a block (== the Rust test cost)."""
    members = RS._members(list(block), genes, gene_of_dn)
    if not (len(members) >= 3 and len(members) <= RS.MAX_LOCI):
        return 0
    recs = RBD.family_exons(members, skel, strand, fa)
    if len(recs) < 3:
        return 0
    n = len(recs)
    cost = 0
    for a in range(n):
        for b in range(n):
            if a == b:
                continue
            for qa in recs[a]["exons"]:
                if len(qa) < RBD.MIN_EXON:
                    continue
                for tb in recs[b]["exons"]:
                    if len(tb) < RBD.MIN_EXON:
                        continue
                    lo, hi = (len(qa), len(tb)) if len(qa) <= len(tb) else (len(tb), len(qa))
                    if lo / hi < 0.5:
                        continue
                    cost += 2 * len(qa) * len(tb)     # x2 for fwd+revcomp
    return cost


# ------------------------------------------------------------------ capture one block (self-contained)
def capture_block(label, cls, block, n_art, genes, gene_of_dn, skel, strand, fa):
    """Run the REAL RS.split_block and ship its exact inputs (genes/gene_of_dn/strand/
    skel subset + raw genome slices) + partition.  Blocks that split_block EARLY-RETURNS
    on (before family_exons: <3 loci / <3 members / >MAX_LOCI) ship NO skel/slices -- the
    Rust test reproduces the identical early return with an empty mock fetcher."""
    blk = list(block)
    result = RS.split_block(blk, genes, gene_of_dn, skel, strand, fa)

    # genes subset (chrom/start/end) for every DN in the block that IS in genes
    genes_out = {}
    for dn in blk:
        g = genes.get(dn)
        if g is not None:
            genes_out[dn] = [g["chrom"], g["start"], g["end"]]
    # gene_of_dn subset (floored) for DNs present
    gof_out = {dn: gene_of_dn[dn] for dn in blk if dn in gene_of_dn}
    # members that split_block will actually build (via _members) -> their skel keys + strand + slices
    members = RS._members(blk, genes, gene_of_dn)
    reaches_family_exons = (len(blk) >= 3 and 3 <= len(members) <= RS.MAX_LOCI)
    strand_out = {}
    skel_out = []
    slices = {}
    seen_keys = set()
    for m in members:
        strand_out[m["dn"]] = strand.get(m["dn"], "+")
        if not reaches_family_exons:
            continue                                   # early return -> no exon graph -> no slices
        key = (m["chrom"], m["start"], m["end"])
        if key in seen_keys:
            continue
        seen_keys.add(key)
        exons = skel.get(key)
        if exons is None:
            continue
        skel_out.append(dict(chrom=m["chrom"], start=m["start"], end=m["end"],
                             exons=[list(e) for e in exons]))
        for (s, e) in exons:
            slk = (m["chrom"], s - 1, e)
            if slk in slices:
                continue
            slices[slk] = fa.fetch(m["chrom"], s - 1, e)   # case-carrying, BEFORE .upper()
    return dict(
        label=label, cls=cls, n_art=n_art,
        block=list(blk),
        genes=genes_out,
        gene_of_dn=gof_out,
        strand=strand_out,
        skel=skel_out,
        slices=[dict(chrom=k[0], start0=k[1], end0=k[2], bytes=v) for k, v in slices.items()],
        result=[list(p) for p in result],
    )


# ------------------------------------------------------------------ articulation-point graph fixtures
def build_articulation_graphs():
    cases = []

    def add(label, n, edges):
        Gr = nx.Graph()
        Gr.add_nodes_from(range(n))
        Gr.add_edges_from(edges)
        art = sorted(set(nx.articulation_points(Gr)))
        cases.append(dict(label=label, n=n, edges=[list(e) for e in edges], art=art))

    add("single_node", 1, [])
    add("two_isolated", 2, [])
    add("dyad_edge", 2, [(0, 1)])                                   # a dyad is biconnected -> no AP
    add("path_p3", 3, [(0, 1), (1, 2)])                            # 1
    add("path_p5", 5, [(0, 1), (1, 2), (2, 3), (3, 4)])            # 1,2,3
    add("cycle_c3", 3, [(0, 1), (1, 2), (2, 0)])                   # none
    add("cycle_c4", 4, [(0, 1), (1, 2), (2, 3), (3, 0)])          # none
    add("cycle_c5", 5, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)])  # none
    add("star_s5", 6, [(0, 1), (0, 2), (0, 3), (0, 4), (0, 5)])   # 0 (center)
    add("k4", 4, [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)])  # none
    add("two_triangles_shared_vertex", 5,
        [(0, 1), (1, 2), (2, 0), (2, 3), (3, 4), (4, 2)])         # 2 (multiply-yielded, deduped)
    add("triangle_with_tail", 4, [(0, 1), (1, 2), (2, 0), (2, 3)])  # 2
    add("barbell", 6, [(0, 1), (0, 2), (1, 2), (2, 3), (3, 4), (3, 5), (4, 5)])  # 2,3
    add("bridge_two_cliques", 8,
        [(0, 1), (0, 2), (1, 2), (3, 4), (3, 5), (4, 5), (2, 6), (6, 3)])  # 2,3,6
    add("tree_binary", 7,
        [(0, 1), (0, 2), (1, 3), (1, 4), (2, 5), (2, 6)])         # 0,1,2
    add("multi_component_two_paths", 6,
        [(0, 1), (1, 2), (3, 4), (4, 5)])                        # 1,4
    add("multi_component_cycle_and_path", 7,
        [(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 6)])        # 4,5
    add("wheel_w5", 6,
        [(0, 1), (0, 2), (0, 3), (0, 4), (0, 5),
         (1, 2), (2, 3), (3, 4), (4, 5), (5, 1)])                 # none
    add("two_cycles_shared_vertex", 7,
        [(0, 1), (1, 2), (2, 0), (2, 3), (3, 4), (4, 5), (5, 2)])  # 2
    add("path_with_leaf", 5, [(0, 1), (1, 2), (2, 3), (1, 4)])    # 1,2
    add("complete_bipartite_k23", 5,
        [(0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4)])        # none
    add("long_chain", 8, [(i, i + 1) for i in range(7)])          # 1..6
    add("cycle_with_two_tails", 6,
        [(0, 1), (1, 2), (2, 0), (0, 3), (3, 4), (2, 5)])        # 0,2,3
    add("isolated_plus_edge", 4, [(1, 2)])                        # none (0,3 isolated; 1-2 dyad)
    add("double_bridge_chain", 9,
        [(0, 1), (1, 2), (2, 0), (2, 3), (3, 4), (4, 5), (5, 3), (5, 6), (6, 7), (7, 8), (8, 6)])
    return cases


# ------------------------------------------------------------------ main
def main():
    print("[load] reconstruct refined / genes / gene_of_dn ...", flush=True)
    refined, genes, gene_of_dn = reconstruct()
    print(f"       refined blocks: {len(refined)}")
    skel = RBD.load_skeletons()
    strand = RBD.load_strand()
    fa = pysam.FastaFile(RBD.GENOME)

    # ---- (b) select representative blocks ---------------------------------------------------------
    # multi-gene >=3-member blocks are the ONLY blocks that can split (a single-gene family is
    # rejected by the same-gene guard); they also host the colocated / NA-dominant / multi-art guards.
    def floored_genes(block):
        return {gene_of_dn.get(dn) for dn in block if gene_of_dn.get(dn) is not None}

    multigene = [b for b in refined if len(b) >= 3 and len(floored_genes(b)) >= 2]
    print(f"       multi-gene >=3-member candidates: {len(multigene)}")

    emitted = []
    seen_blocks = set()

    def emit(label, cls, block, n_art):
        key = tuple(sorted(block))
        if key in seen_blocks:
            return False
        seen_blocks.add(key)
        emitted.append(capture_block(label, cls, block, n_art, genes, gene_of_dn, skel, strand, fa))
        return True

    # trace every candidate; bucket (with a HW-DP cost proxy so we pick the CHEAPEST representative
    # per required guard-class -> the Rust byte-parity test stays snappy).
    buckets = collections.defaultdict(list)   # tag -> [(cost, block, n_art, reasons)]
    for b in multigene:
        tag, n_art, reasons = trace(b, genes, gene_of_dn, skel, strand, fa)
        cost = dp_cost(b, genes, gene_of_dn, skel, strand, fa)
        buckets[tag].append((cost, b, n_art, reasons))
    print(f"       traced: split={len(buckets['split'])} nonsplit_graph={len(buckets['nonsplit_graph'])} "
          f"no_edges={len(buckets['no_edges'])} recs_lt3={len(buckets['skip_recs_lt3'])}")

    # (1) EVERY split block (the 2 recall-safe HIGH-confidence mosaics: fid 210 + fid 187)
    for i, (cost, b, n_art, reasons) in enumerate(sorted(buckets["split"], key=lambda x: sorted(x[1]))):
        parts = RS.split_block(list(b), genes, gene_of_dn, skel, strand, fa)
        doms = []
        for p in parts:
            gs = [gene_of_dn.get(dn) for dn in p if gene_of_dn.get(dn) not in (None, "NA")]
            doms.append(collections.Counter(gs).most_common(1)[0][0] if gs else "NA")
        lab = "flagship_210" if set(g for g in doms) >= {"GALNT17", "LOC101126070"} else f"split_{i}"
        emit(lab, "split", b, n_art)
        print(f"       SPLIT {lab}: {len(parts)} parts, subfam_genes={doms} (dp_cost={cost:,})")

    # (2) colocated-chimera skip (READTHROUGH) -- CHEAPEST block that reaches the subfams_colocated guard
    for cost, b, n_art, reasons in sorted(buckets["nonsplit_graph"], key=lambda x: x[0]):
        if "colocated" in reasons and emit("colocated_skip", "colocated_skip", b, n_art):
            print(f"       colocated_skip dp_cost={cost:,}"); break
    # (3) NA-dominant / low-confidence skip -- CHEAPEST
    for cost, b, n_art, reasons in sorted(buckets["nonsplit_graph"], key=lambda x: x[0]):
        if "low_conf" in reasons and emit("na_dominant_skip", "na_dominant_skip", b, n_art):
            print(f"       na_dominant_skip dp_cost={cost:,}"); break
    # (4) multi-articulation-point family (>=2 cut vertices) -- CHEAPEST
    for cost, b, n_art, reasons in sorted(buckets["nonsplit_graph"] + buckets["split"], key=lambda x: x[0]):
        if n_art >= 2 and emit("multi_articulation", "multi_articulation", b, n_art):
            print(f"       multi_articulation n_art={n_art} dp_cost={cost:,}"); break
    # (5) the CHEAPEST representative nonsplit-with-graph blocks (keeps the Rust HW-DP test snappy)
    for cost, b, n_art, reasons in sorted(buckets["nonsplit_graph"], key=lambda x: x[0]):
        if sum(1 for e in emitted if e["cls"] == "nonsplit_graph") >= 14:
            break
        emit(f"nonsplit_graph_{len(emitted)}", "nonsplit_graph", b, n_art)

    # (6) <3-member skip blocks (real 2-member refined families; cost-free -- no exon graph)
    small = [b for b in refined if len(b) == 2]
    for b in small[:2]:
        emit(f"skip_lt3_{len(emitted)}", "skip_block_lt3", b, 0)

    # (8) SYNTHETIC >MAX_LOCI skip (guaranteed early return; needs NO slices -> tiny).  Build a block
    #     of MAX_LOCI+5 real DN ids that ARE in genes so _members keeps them (>60 -> early return).
    all_dn_in_genes = [dn for b in refined for dn in b if dn in genes]
    big = sorted(set(all_dn_in_genes))[: RS.MAX_LOCI + 5]
    if len(big) > RS.MAX_LOCI:
        emit("synthetic_gt_maxloci", "skip_gt_maxloci", big, 0)

    fa.close()

    # ---- (a) articulation graphs -----------------------------------------------------------------
    art_graphs = build_articulation_graphs()
    print(f"       articulation graphs: {len(art_graphs)}")

    # coverage report
    classes = collections.Counter(e["cls"] for e in emitted)
    print(f"[emit] {len(emitted)} blocks; classes={dict(classes)}")

    fixture = dict(
        MAX_LOCI=RS.MAX_LOCI, COLINEAR_COV=RS.COLINEAR_COV, ID_THRESH=RS.ID_THRESH,
        articulation_graphs=art_graphs,
        blocks=emitted,
    )
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as fh:
        json.dump(fixture, fh, indent=1, sort_keys=True)
    print(f"[fixture] wrote {os.path.normpath(OUT)} ({os.path.getsize(OUT):,} bytes)")


if __name__ == "__main__":
    main()
