#!/usr/bin/env python3
"""family_def_community_debridge.py — de-bridge over-merged components by co-threading-WEIGHTED
community detection (the lever a global frac_ref threshold could not provide).

A global threshold failed because real partial-paralogs (GGT1~GGTLC2, frac_ref 0.20) overlap the
bridge-edge range, and because the dense comp-238 blob co-threads a shared repeat so NO single
threshold splits it. Community detection sidesteps both: within each large component it finds
sub-groups that co-thread MORE junctions internally than across the bridge, weighting each edge by
its junction-position co-threading. Small components (real families) are left untouched.

Pipeline:
  1. post-retrocopy graph (optionally also drop antisense '-' edges: --strand).
  2. edge weight = frac_ref (junction-position co-threading) from family_def_cothread_genomewide.tsv;
     edges with no measured weight get a small epsilon (weak link).
  3. for each component with >= MIN_COMP members: louvain_communities(weight, seed=0, resolution=R).
     final families = small components (unchanged) + the communities of each large component.
  4. validate: (a) named real families stay inside ONE community (no real-family fragmentation);
     (b) the mutually-UNRELATED named genes that anchor each over-merge (e.g. comp-238:
     FOXD4L1/NDUFA8/RHBDF1/SLC35F5) land in DIFFERENT communities (the bridge is cut);
     (c) family-count / max-component change.

Run: /home/juanfra/miniforge3/bin/python bench/family_def_community_debridge.py [--strand] [--res R] [--weight frac_ref|co]
Out: bench/family_def_community_debridge.json
"""
import collections
import json
import os
import re
import sys

import networkx as nx
from networkx.algorithms import community as nxcom

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from family_def_read_filters import dna_homology
from family_def_retrocopy_filter import build_family_graph, load_intron_index
from family_def_strand_probe import edge_strand

COTHREAD_TSV = os.path.join(HERE, "family_def_cothread_genomewide.tsv")
MIN_COMP = 20
EPS = 0.02                       # weak-link weight for edges with no co-threading measurement

NAMED_REAL = [("RABL2A", "RABL2B"), ("APOBEC3D", "APOBEC3F"),
              ("RFPL1", "RFPL2"), ("RFPL2", "RFPL3"), ("GGT1", "GGTLC2")]
# mutually-unrelated named anchors of the big over-merges (should end up SEPARATED)
UNRELATED_ANCHORS = {
    "comp-238": ["FOXD4L1", "NDUFA8", "RHBDF1", "SLC35F5"],
    "comp-201": ["ACTR3B", "AGAP5", "ARPC2", "BCO2", "CASR", "CCDC62"],
    "comp-140": ["ATG10", "BCR", "CPO", "EPPK1", "FUT3", "IFNE"],
}


def _cothread_weights(key):
    w = {}
    if not os.path.exists(COTHREAD_TSV):
        return w
    with open(COTHREAD_TSV) as f:
        h = f.readline().rstrip("\n").split("\t"); ix = {c: i for i, c in enumerate(h)}
        for line in f:
            t = line.rstrip("\n").split("\t")
            e = tuple(sorted((t[ix["ref"]], t[ix["mem"]])))
            try:
                w[e] = float(t[ix[key]])
            except (ValueError, KeyError):
                pass
    if key == "co":   # normalise raw counts to [0,1]-ish for modularity stability
        mx = max(w.values()) if w else 1
        w = {e: v / mx for e, v in w.items()}
    return w


def _covmin_weights(Hd):
    """cov_min = min(cov_a, cov_b) straight from the all-vs-all -- FREE, no alignment.
    Best de-bridge weight (max comp 238->44, reals safe): real paralogs cover each other
    reciprocally; domain bridges cover asymmetrically -> low cov_min."""
    return {e: min(r["cov_a"], r["cov_b"]) for e, r in Hd.items()}


def load_weights(key, Hd):
    """key in {cov_min (default), frac_ref/frac_mem/aln_cov_small/contig/co, product}.
    'product' = cov_min * frac_ref (mechanism-grounded + alignment refinement)."""
    if key == "cov_min":
        return _covmin_weights(Hd)
    if key == "product":
        cm = _covmin_weights(Hd); fr = _cothread_weights("frac_ref")
        return {e: cm[e] * fr.get(e, 1.0) for e in cm}
    return _cothread_weights(key)


def root(g):
    return re.sub(r"[0-9]+[A-Z]*$", "", g)


def main():
    use_strand = "--strand" in sys.argv
    res = float(sys.argv[sys.argv.index("--res") + 1]) if "--res" in sys.argv else 1.0
    # frac_ref (co-threading) is the de-bridge weight VALIDATED against protein orthogroups:
    # it recovers the protein-family partition at F1 0.70-0.72 vs cov_min 0.62-0.64 (cov_min's
    # smaller max-comp was OVER-FRAGMENTATION -- it splits real protein families; see
    # family_def_protein_validate.py). For coding genes prefer protein orthogroups directly;
    # this de-bridge is the cDNA-level proxy where protein is unavailable.
    wkey = sys.argv[sys.argv.index("--weight") + 1] if "--weight" in sys.argv else "frac_ref"

    Hd, _ = dna_homology()
    idx = load_intron_index()
    G, prov, _ = build_family_graph(Hd, idx, filter_retrocopy=True)
    if use_strand:
        strand = edge_strand()
        drop = [e for e in G.edges() if strand.get(tuple(sorted(e))) == "-"]
        G.remove_edges_from(drop)
        G.remove_nodes_from([n for n in list(G.nodes()) if G.degree(n) == 0])
        print(f"[strand] dropped {len(drop)} antisense edges")

    W = load_weights(wkey, Hd)
    for a, b in G.edges():
        G[a][b]["w"] = W.get(tuple(sorted((a, b))), EPS)

    comps = sorted([c for c in nx.connected_components(G) if len(c) >= 2], key=len, reverse=True)
    base_n = len(comps)
    base_max = len(comps[0])

    # de-bridge each large component via weighted louvain; keep small ones intact
    final = []
    big_report = []
    comp_anchor_idx = {}     # map our hardcoded anchors to the actual largest comps in order
    for i, c in enumerate(comps):
        if len(c) < MIN_COMP:
            final.append(set(c)); continue
        sub = G.subgraph(c)
        coms = nxcom.louvain_communities(sub, weight="w", resolution=res, seed=0)
        coms = [set(x) for x in coms]
        final.extend([x for x in coms if len(x) >= 2])
        # name composition of the resulting communities
        named = [g for g in c if not g.startswith("LOC")]
        comp_anchor_idx[f"comp-{len(c)}"] = (c, coms)
        big_report.append(dict(orig_size=len(c), n_named=len(named),
                               communities=sorted((len(x) for x in coms), reverse=True)[:8],
                               n_comm=len(coms)))

    final_sizes = sorted((len(c) for c in final), reverse=True)
    # node -> final community id
    node2fam = {}
    for fi, c in enumerate(final):
        for g in c:
            node2fam[g] = fi

    print(f"=== co-threading community de-bridge (weight={wkey}, res={res}, strand={use_strand}) ===")
    print(f"  base : {base_n} families, max comp {base_max}, top {[len(c) for c in comps[:5]]}")
    print(f"  final: {len(final)} families, max comp {final_sizes[0]}, top {final_sizes[:6]}")

    print("\n=== big-component breakup ===")
    for r in big_report:
        print(f"  orig n={r['orig_size']:3d} (named {r['n_named']:3d}) -> {r['n_comm']} communities, "
              f"sizes {r['communities']}")

    # (a) named real families intact?
    print("\n=== (a) real families NOT fragmented? ===")
    real_ok = 0
    for a, b in NAMED_REAL:
        fa, fb = node2fam.get(a), node2fam.get(b)
        ok = fa is not None and fa == fb
        real_ok += ok
        print(f"  {a}~{b}: {'TOGETHER' if ok else 'SPLIT'}  (fam {fa} / {fb})")

    # (b) unrelated anchors separated?
    print("\n=== (b) unrelated anchors SEPARATED (bridge cut)? ===")
    for comp_label, anchors in UNRELATED_ANCHORS.items():
        present = [(g, node2fam.get(g)) for g in anchors if g in node2fam]
        fams = {f for _, f in present if f is not None}
        print(f"  {comp_label}: anchors {[g for g,_ in present]} -> "
              f"{len(fams)} distinct families "
              f"({'SEPARATED' if len(fams) == len([g for g,_ in present]) else 'some still merged'})")

    # (c) coherence: among resulting families with >=2 named genes, root-name homogeneity
    homo = []
    for c in final:
        named = [g for g in c if not g.startswith("LOC")]
        if len(named) >= 2:
            rc = collections.Counter(root(g) for g in named)
            homo.append(rc.most_common(1)[0][1] / len(named))
    print(f"\n=== (c) coherence: mean root-name homogeneity over {len(homo)} multi-named families: "
          f"{sum(homo)/max(len(homo),1):.2f} (1.0 = each family one root) ===")

    json.dump(dict(weight=wkey, resolution=res, strand=use_strand,
                   base_families=base_n, base_max=base_max,
                   final_families=len(final), final_max=final_sizes[0],
                   final_top=final_sizes[:10], big_report=big_report,
                   real_intact=real_ok, real_total=len(NAMED_REAL),
                   mean_homogeneity=round(sum(homo)/max(len(homo),1), 3)),
              open(os.path.join(HERE, "family_def_community_debridge.json"), "w"), indent=2)
    print(f"\n[+] wrote family_def_community_debridge.json")


if __name__ == "__main__":
    main()
