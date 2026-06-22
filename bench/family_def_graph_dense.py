#!/usr/bin/env python3
"""
family_def_graph_dense.py — do CLIQUE / DENSE-CLUSTER family objects beat connected
components? (k-truss, k-core, maximal cliques, vs CC.)

Two things to settle:
  (1) size-2 problem: a single paralog pair has no triangle/density, so any
      density/clique/truss object discards it. How many families is that?
  (2) on the LARGE components (where the over-merge bridges live), does a dense object
      cut the right edges? Confound-free metric = DNA-paralog rate of CUT vs KEPT edges
      (granularity inflates within-family purity, so we score the edges that are cut).

Reads bench/family_def_dna_pr_edges.tsv (in_rna / in_dna_loose / chrom / n_tied).
Run: /home/juanfra/miniforge3/bin/python bench/family_def_graph_dense.py
"""
import collections
import os
import networkx as nx

HERE = os.path.dirname(os.path.abspath(__file__))
DUMP = os.path.join(HERE, "family_def_dna_pr_edges.tsv")


def load():
    G = nx.Graph(); dna = set(); g2chrom = {}
    with open(DUMP) as f:
        h = f.readline().rstrip("\n").split("\t"); idx = {k: i for i, k in enumerate(h)}
        for line in f:
            c = line.rstrip("\n").split("\t")
            a, b = c[idx["gene_a"]], c[idx["gene_b"]]
            g2chrom[a] = c[idx["chrom_a"]]; g2chrom[b] = c[idx["chrom_b"]]
            if c[idx["in_dna_loose"]] == "1":
                dna.add(frozenset((a, b)))
            if c[idx["in_rna"]] == "1":
                G.add_edge(a, b, w=int(c[idx["n_tied"]]))
    return G, dna, g2chrom


def edge_cut_ratio(G, node2fam, dna):
    kp = kt = cp = ct = 0
    for u, v in G.edges():
        par = frozenset((u, v)) in dna
        if node2fam.get(u) is not None and node2fam.get(u) == node2fam.get(v):
            kt += 1; kp += par
        else:
            ct += 1; cp += par
    kr = kp / kt if kt else 0; cr = cp / ct if ct else 0
    return kr, cr, (kr / cr if cr else float("inf")), kt, ct


def fams_from_subgraph(H):
    return [set(c) for c in nx.connected_components(H) if len(c) >= 2]


def report(name, G, fams, dna, g2chrom):
    node2fam = {}
    for i, c in enumerate(fams):
        for g in c:
            node2fam.setdefault(g, i)
    size2 = sum(1 for c in fams if len(c) == 2)
    def nch(c): return len({g2chrom.get(g) for g in c})
    bridges = sum(1 for c in fams if nch(c) >= 3)
    genes = sum(len(c) for c in fams)
    rabl2 = any({"RABL2A", "RABL2B"} <= c for c in fams)
    kr, cr, ratio, kt, ct = edge_cut_ratio(G, node2fam, dna)
    print(f"  {name:20} fams={len(fams):4d} size2={size2:4d} xbridge={bridges:3d} genes={genes:5d} "
          f"RABL2={'Y' if rabl2 else 'n'}  KEPTpar={kr:.2f} CUTpar={cr:.2f} ratio={ratio:.2f} cut={ct}")


def main():
    G, dna, g2chrom = load()
    comps = list(nx.connected_components(G))
    sizes = collections.Counter(len(c) for c in comps if len(c) >= 2)
    n2 = sizes[2]; ntot = sum(sizes.values())
    print(f"graph: {G.number_of_edges()} edges, {G.number_of_nodes()} nodes, {ntot} families")
    print(f"** size-2 families = {n2}/{ntot} = {100*n2/ntot:.0f}% (a dense/clique/truss object discards ALL of these) **\n")

    print(f"  {'method':20} {'families':>8}  ... KEPT/CUT DNA-paralog rate (ratio>>1 = cuts the right edges)")
    report("CC (current)", G, fams_from_subgraph(G), dna, g2chrom)
    report("k-core>=2", G, fams_from_subgraph(nx.k_core(G, 2)), dna, g2chrom)
    report("k-truss>=3 (triangle)", G, fams_from_subgraph(nx.k_truss(G, 3)), dna, g2chrom)
    report("k-truss>=4", G, fams_from_subgraph(nx.k_truss(G, 4)), dna, g2chrom)

    # HYBRID: keep size-2 / size-3 as-is, dense-decompose only the LARGE components
    print("\n--- HYBRID: keep small families (<=4), apply k-truss>=3 only to large components ---")
    hybrid = []
    cut_big = 0
    for c in comps:
        if len(c) < 2:
            continue
        if len(c) <= 4:
            hybrid.append(set(c)); continue
        sub = G.subgraph(c)
        H = nx.k_truss(sub, 3)
        before = sub.number_of_edges(); after = H.number_of_edges()
        cut_big += before - after
        cores = fams_from_subgraph(H)
        if cores:
            hybrid.extend(cores)
        # nodes dropped from the truss (no triangle) -> their surviving edges form smaller fams
        leftover = set(c) - set().union(*cores) if cores else set(c)
        if leftover:
            lg = G.subgraph(leftover)
            hybrid.extend(fams_from_subgraph(lg))
    report(f"hybrid (cut {cut_big} big-comp edges)", G, hybrid, dna, g2chrom)

    # maximal cliques: how fragmented, how pure
    cliques = [c for c in nx.find_cliques(G) if len(c) >= 2]
    csz = collections.Counter(len(c) for c in cliques)
    pure = 0; wp = 0
    import itertools
    for c in cliques:
        for u, v in itertools.combinations(c, 2):
            wp += 1; pure += frozenset((u, v)) in dna
    print(f"\n--- maximal cliques: {len(cliques)} cliques (vs {ntot} CC families); "
          f"size dist {dict(sorted(csz.items())[:6])}; within-clique DNA-paralog purity={pure/wp if wp else 0:.2f}")
    print("  (a single edge is its own 2-clique, so cliques don't discard size-2 — but they over-count "
          "dense components into many overlapping cliques)")


if __name__ == "__main__":
    main()
