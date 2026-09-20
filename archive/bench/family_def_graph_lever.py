#!/usr/bin/env python3
"""family_def_graph_lever.py — unbiased lever search in GRAPH / higher-dimensional space.

Pairwise sequence/coverage features are all soft. New idea: the TOPOLOGY of the cross-mapping
graph. A REAL family is a CLIQUE (all copies mutually cross-map -> edges sit in triangles, share
neighbors); a BRIDGE is a STAR/hub (a shared domain links genes that don't connect to each other
-> the bridge edge has no common neighbors, no triangles). Build the ~R gene-gene graph, compute
graph features per edge (common neighbors, neighborhood Jaccard, triangles, degrees, Adamic-Adar,
k-truss support), and run the unbiased discovery (label-free bimodality + EXTERNAL DNA AUC +
permutation null) -- does graph structure separate bridge from paralog where pairwise could not?
Run: python bench/family_def_graph_lever.py NC_073240.2 NC_073244.2 NC_073248.2 [more]
"""
import collections
import json
import math
import os
import sys

import numpy as np
import networkx as nx

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from family_def_genomewide import pair_evidence, components, GENES_BED, best_gene, DELTA, DE_MAX, MIN_READS
from family_def_newbam_validate import scan_region, build_model, recip_cov, GENOME, NEW, COV_MIN
from family_def_feature_discovery import auc, bimodality
from family_def_read_filters import dna_homology
import pysam

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    chroms = set(sys.argv[1:] or ["NC_073240.2", "NC_073244.2", "NC_073248.2"])
    by = collections.defaultdict(list)
    with open(GENES_BED) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 4 and p[0] in chroms:
                by[p[0]].append((int(p[1]), int(p[2]), p[3]))
    for c in by:
        by[c].sort()
    # ~R scan: region-restricted over the new BAM per chrom, merged (captures cross-chrom edges too)
    mm = collections.defaultdict(dict)
    for c in chroms:
        mc, _ = scan_region(by, NEW, c)
        for q, gd in mc.items():
            for g, de in gd.items():
                if g not in mm[q] or de < mm[q][g]:
                    mm[q][g] = de
    ev = pair_evidence(mm)
    edges, _ = components(ev, DELTA, DE_MAX, MIN_READS)   # (ga, gb, n)
    G = nx.Graph()
    for ga, gb, n in edges:
        G.add_edge(ga, gb, w=n)
    deg = dict(G.degree())
    Hd, _ = dna_homology()

    def dna(ga, gb):
        r = Hd.get((ga, gb) if ga < gb else (gb, ga))
        if r is None or r.get("id", 0) == 0:
            return 1
        return 0 if (r["id"] >= 0.90 and max(r["cov_a"], r["cov_b"]) >= 0.30) else None

    # gene coords + copy models for ~B (so we can test graph-lever vs / with ~B)
    coord = {}
    by2 = collections.defaultdict(list)
    with open(GENES_BED) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 4 and p[0] in chroms:
                coord[p[3]] = (p[0], int(p[1]), int(p[2])); by2[p[0]].append((int(p[1]), int(p[2]), p[3]))
    bam = pysam.AlignmentFile(NEW, "rb"); genome = pysam.FastaFile(GENOME)
    mcache = {}
    def rcov(ga, gb):
        for g in (ga, gb):
            if g not in mcache and g in coord:
                mcache[g] = build_model(bam, genome, *coord[g])
        return recip_cov(mcache.get(ga, ""), mcache.get(gb, "")) if ga in coord and gb in coord else 0.0

    rows, labels, jacc_recip = [], [], []
    for ga, gb, n in edges:
        lab = dna(ga, gb)
        if lab is None:
            continue
        jacc_recip.append((set(G[ga]) - {gb}, set(G[gb]) - {ga}, ga, gb, lab))
        Na, Nb = set(G[ga]) - {gb}, set(G[gb]) - {ga}
        common = Na & Nb
        uni = Na | Nb
        # k-truss support = # common neighbors (triangles through the edge)
        aa = sum(1.0 / math.log(deg[w]) for w in common if deg[w] > 1)         # Adamic-Adar
        ra = sum(1.0 / deg[w] for w in common if deg[w] > 0)                   # resource allocation
        feat = dict(
            common_neighbors=len(common),
            jaccard_nbr=len(common) / max(len(uni), 1),
            triangles=len(common),                       # = support
            deg_min=min(deg[ga], deg[gb]),
            deg_max=max(deg[ga], deg[gb]),
            deg_prod=deg[ga] * deg[gb],
            adamic_adar=aa,
            resource_alloc=ra,
            in_triangle=1.0 if common else 0.0,
        )
        rows.append(feat); labels.append(lab)

    labels = np.array(labels)
    feats = list(rows[0].keys())
    n_par, n_brg = int((labels == 0).sum()), int((labels == 1).sum())
    print(f"=== GRAPH-topology lever search: {len(rows)} ~R edges "
          f"({n_par} DNA-paralog / {n_brg} bridge), graph |V|={G.number_of_nodes()} |E|={G.number_of_edges()} ===")
    print(f"  {'graph feature':18} {'bimod':>6} {'AUC_vs_DNA':>11}  paralog-med / bridge-med")
    res = []
    for fk in feats:
        vals = [r[fk] for r in rows]
        a = auc(vals, labels); sep, _ = bimodality(vals)
        pm = np.median([r[fk] for r, l in zip(rows, labels) if l == 0])
        bm = np.median([r[fk] for r, l in zip(rows, labels) if l == 1])
        res.append((fk, sep, a))
        print(f"  {fk:18} {sep:>6.2f} {a:>11.3f}  {pm:.2f} / {bm:.2f}"
              f"{'  <== separates' if a >= 0.70 else ''}")
    res.sort(key=lambda t: -t[2])

    rng = np.random.default_rng(1729)
    real_max = max(a for _, _, a in res)
    null = [max(auc([r[fk] for r in rows], rng.permutation(labels)) for fk in feats) for _ in range(2000)]
    p = float((np.array(null) >= real_max).mean())
    print(f"\n  best graph lever: {res[0][0]} AUC={res[0][2]:.3f}; permutation p={p:.4f} -> "
          f"{'real lever' if p < 0.05 else 'not past chance'}")
    # operating point of the best lever (does it cleanly cut bridges?)
    bfk = res[0][0]
    pv = sorted(r[bfk] for r, l in zip(rows, labels) if l == 0)
    bv = sorted(r[bfk] for r, l in zip(rows, labels) if l == 1)
    print(f"  TRUE paralog {bfk}: p10={pv[len(pv)//10]:.2f} median={pv[len(pv)//2]:.2f}")
    lo = np.mean(bv) < np.mean(pv)
    for q in [0.25, 0.5, 0.75]:
        thr = np.quantile(bv, q) if lo else np.quantile(bv, 1 - q)
        bc = sum(1 for x in bv if (x <= thr if lo else x >= thr))
        pc = sum(1 for x in pv if (x <= thr if lo else x >= thr))
        print(f"    thr={thr:.2f}: cuts {bc}/{n_brg} bridges, costs {pc}/{n_par} paralogs")
    # === does the graph lever (jaccard_nbr) close the ~B residual? combine with ~B ===
    print("\n  === graph lever vs / with ~B (does topology close what sequence-homology can't?) ===")
    comb = []
    for Na, Nb, ga, gb, lab in jacc_recip:
        j = len(Na & Nb) / max(len(Na | Nb), 1)
        comb.append((j, rcov(ga, gb), lab))
    bam.close(); genome.close()
    Bsurv = [(j, l) for j, rc, l in comb if rc >= COV_MIN]      # ~B survivors
    bt, bf = sum(1 for _, l in Bsurv if l == 0), sum(1 for _, l in Bsurv if l == 1)
    print(f"  ~B survivors: {len(Bsurv)} ({bt} paralog / {bf} bridge = the residual ~B can't remove)")
    if bf:
        a_res = auc([j for j, _ in Bsurv], np.array([l for _, l in Bsurv]))
        print(f"  jaccard_nbr AUC on the ~B-RESIDUAL (paralog vs surviving bridge): {a_res:.3f}")
        for thr in [0.20, 0.33, 0.50]:
            bc = sum(1 for j, l in Bsurv if l == 1 and j < thr)
            pc = sum(1 for j, l in Bsurv if l == 0 and j < thr)
            print(f"    ~B ∩ jaccard>={thr}: removes {bc}/{bf} residual bridges, costs {pc}/{bt} paralogs")
    json.dump(dict(n_paralog=n_par, n_bridge=n_brg, nodes=G.number_of_nodes(), edges=G.number_of_edges(),
                   ranking=[dict(feature=f, bimod=s, dna_auc=a) for f, s, a in res], perm_p=p,
                   B_survivors=len(Bsurv), B_residual_bridges=bf),
              open(os.path.join(HERE, "family_def_graph_lever.json"), "w"), indent=2)
    print("\n[+] wrote family_def_graph_lever.json")


if __name__ == "__main__":
    main()
