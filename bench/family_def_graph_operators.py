#!/usr/bin/env python3
"""
family_def_graph_operators.py — is there a CLEANER family object than the connected
component of the read-conflict graph?

The one impurity in "family = connected component" is that read-indistinguishability is
NOT transitive: one spurious cross-locus edge (a retrocopy/repeat bridge) fuses two
real clusters. Connected components force transitivity; a real multi-copy array is a
DENSE block, a bridge is a sparse cut-edge. We compare standard, (near) parameter-free
graph objects on the genome-wide gene-vertex read-conflict graph and score each against
the independent DNA paralogy truth:

  - CC                : connected components (the current definition)
  - biconnected       : biconnected components (+ isolated edges as size-2 families)
  - cut-edge-pruned   : drop bridge edges whose BOTH endpoints have degree>=2 (so size-2
                        isolated families survive), then connected components
  - modularity        : greedy-modularity communities (networkx)

Metrics (higher purity + fewer cross-chrom bridges + arrays retained = cleaner):
  - families, size-2 count, cross-chrom-bridge families (>=3 chromosomes = artifact)
  - DNA-purity = within-family gene-pairs that are DNA paralog edges / all within-family pairs
  - retention of the panel families (RABL2A~RABL2B; MAGEA arrays via LOC co-location)
  - genes covered (a recall proxy)

Reads bench/family_def_dna_pr_edges.tsv (in_rna / in_dna_loose / chrom / n_tied).
Run: /home/juanfra/miniforge3/bin/python bench/family_def_graph_operators.py
"""
import collections
import itertools
import json
import os
import networkx as nx

HERE = os.path.dirname(os.path.abspath(__file__))
DUMP = os.path.join(HERE, "family_def_dna_pr_edges.tsv")


def load():
    rna = nx.Graph()
    dna_loose = set()
    g2chrom = {}
    with open(DUMP) as f:
        h = f.readline().rstrip("\n").split("\t"); idx = {k: i for i, k in enumerate(h)}
        for line in f:
            c = line.rstrip("\n").split("\t")
            a, b = c[idx["gene_a"]], c[idx["gene_b"]]
            g2chrom[a] = c[idx["chrom_a"]]; g2chrom[b] = c[idx["chrom_b"]]
            if c[idx["in_dna_loose"]] == "1":
                dna_loose.add(frozenset((a, b)))
            if c[idx["in_rna"]] == "1":
                rna.add_edge(a, b, w=int(c[idx["n_tied"]]))
    return rna, dna_loose, g2chrom


def families_CC(G):
    return [c for c in nx.connected_components(G) if len(c) >= 2]


def families_biconnected(G):
    fams = []
    for block in nx.biconnected_components(G):   # sets of >=2 nodes; isolated edges included
        if len(block) >= 2:
            fams.append(set(block))
    return fams


def families_cutedge_pruned(G):
    H = G.copy()
    # graph-theoretic bridges (edges in no cycle) whose endpoints are BOTH non-trivial
    bridges = [(u, v) for u, v in nx.bridges(G) if G.degree[u] >= 2 and G.degree[v] >= 2]
    H.remove_edges_from(bridges)
    return [c for c in nx.connected_components(H) if len(c) >= 2], len(bridges)


def families_modularity(G):
    fams = []
    for comp in nx.connected_components(G):
        sub = G.subgraph(comp)
        if len(comp) <= 3:
            fams.append(set(comp)); continue
        comms = nx.community.greedy_modularity_communities(sub, weight="w")
        for cm in comms:
            if len(cm) >= 2:
                fams.append(set(cm))
    return fams


def score(name, fams, dna_loose, g2chrom, extra=""):
    n = len(fams)
    size2 = sum(1 for c in fams if len(c) == 2)
    def nchrom(c): return len({g2chrom.get(g) for g in c})
    bridges = sum(1 for c in fams if nchrom(c) >= 3)
    genes = sum(len(c) for c in fams)
    # DNA purity: within-family pairs that are DNA paralog edges / all within-family pairs
    wp = dp = 0
    for c in fams:
        for u, v in itertools.combinations(c, 2):
            wp += 1
            if frozenset((u, v)) in dna_loose:
                dp += 1
    purity = dp / wp if wp else 0.0
    rabl2 = any({"RABL2A", "RABL2B"} <= c for c in fams)
    print(f"  {name:18} fams={n:4d} size2={size2:4d} xchrom-bridge={bridges:3d} "
          f"genes={genes:5d} DNA-purity={purity:.3f} RABL2={'Y' if rabl2 else 'n'} {extra}")
    return dict(name=name, families=n, size2=size2, xchrom_bridges=bridges,
                genes=genes, dna_purity=purity, rabl2=rabl2)


def main():
    G, dna_loose, g2chrom = load()
    print(f"genome-wide read-conflict graph: {G.number_of_edges()} edges, "
          f"{G.number_of_nodes()} vertices\n")
    print("operator            families   size2  xchrom-bridges  genes  DNA-purity")
    res = []
    res.append(score("CC (current)", families_CC(G), dna_loose, g2chrom))
    res.append(score("biconnected", families_biconnected(G), dna_loose, g2chrom))
    cut_fams, n_bridges = families_cutedge_pruned(G)
    res.append(score("cut-edge-pruned", cut_fams, dna_loose, g2chrom,
                     extra=f"(removed {n_bridges} inter-cluster bridges)"))
    res.append(score("modularity", families_modularity(G), dna_loose, g2chrom))

    # ---- confound-free test: are the CUT edges disproportionately NON-paralog? ----
    # (granularity raises pairwise purity mechanically; this measures whether the
    #  operator cuts the RIGHT edges, independent of family sizes.)
    def edge_cut_analysis(name, fams):
        node2fam = {}
        for i, c in enumerate(fams):
            for g in c:
                node2fam.setdefault(g, i)  # first family wins (biconnected shares cut-vertices)
        kept_par = kept_tot = cut_par = cut_tot = 0
        for u, v, d in G.edges(data=True):
            par = frozenset((u, v)) in dna_loose
            if node2fam.get(u) is not None and node2fam.get(u) == node2fam.get(v):
                kept_tot += 1; kept_par += par
            else:
                cut_tot += 1; cut_par += par
        kr = kept_par / kept_tot if kept_tot else 0
        cr = cut_par / cut_tot if cut_tot else 0
        print(f"  {name:18} KEPT edges DNA-paralog={kr:.3f} ({kept_par}/{kept_tot})   "
              f"CUT edges DNA-paralog={cr:.3f} ({cut_par}/{cut_tot})   ratio={kr/cr if cr else float('inf'):.2f}")
    print("\n--- confound-free: DNA-paralog rate of KEPT (intra-family) vs CUT edges ---")
    print("  (a correct cleanup KEEPS paralog edges and CUTS non-paralog ones -> ratio >> 1)")
    edge_cut_analysis("biconnected", families_biconnected(G))
    edge_cut_analysis("cut-edge-pruned", cut_fams)
    edge_cut_analysis("modularity", families_modularity(G))

    # what cut-edge pruning specifically removes (sanity on the bridges)
    print("\n--- inter-cluster bridge edges removed by cut-edge pruning (sample) ---")
    bridges = [(u, v) for u, v in nx.bridges(G) if G.degree[u] >= 2 and G.degree[v] >= 2]
    xchrom = [(u, v) for u, v in bridges if g2chrom.get(u) != g2chrom.get(v)]
    print(f"  total inter-cluster bridges removed: {len(bridges)}  "
          f"({len(xchrom)} cross-chromosome)")
    for u, v in bridges[:12]:
        dna = "DNA-paralog" if frozenset((u, v)) in dna_loose else "no-DNA"
        sc = "same-chrom" if g2chrom.get(u) == g2chrom.get(v) else "CROSS-chrom"
        print(f"     {u} ~ {v}   {sc}  {dna}")

    json.dump(res, open(os.path.join(HERE, "family_def_graph_operators.json"), "w"), indent=2)
    print(f"\n[+] wrote {os.path.join(HERE,'family_def_graph_operators.json')}")


if __name__ == "__main__":
    main()
