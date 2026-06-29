#!/usr/bin/env python3
"""family_def_strand_probe.py — does cDNA-alignment STRAND fix the over-merge?

rep_cdna.fa is in transcribed (mRNA-sense) orientation, so a homology row with strand '-' means
gene A's sense cDNA matches gene B's ANTISENSE -- a repeat/transposon inserted in opposite
orientation, or a sense/antisense relationship, NOT a co-linear paralog duplication. dna_homology()
currently ignores strand. Test: are the over-merge components (esp. the dense comp-238 LOC blob)
bridged by opposite-strand edges, and does dropping them break the over-merge while keeping real
paralogs connected?

Per edge we record, among rows that meet the homology bar (id>=0.90 & cov>=0.30), whether the
qualifying alignment is same-strand ('+') or opposite ('-') (taking the strand of the max-id
qualifying row). An edge is 'antisense' if its qualifying homology is '-'.
Run: /home/juanfra/miniforge3/bin/python bench/family_def_strand_probe.py
"""
import collections
import os
import sys

import networkx as nx

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from family_def_read_filters import dna_homology
from family_def_retrocopy_filter import build_family_graph, load_intron_index

AVA = "/home/juanfra/winloci_scratch/rep_ava.tsv"
ID_MIN, COV_MIN = 0.90, 0.30


def edge_strand():
    """edge (a<b) -> '+' or '-' : strand of the max-id qualifying alignment."""
    best = {}     # edge -> (id, strand)
    with open(AVA) as f:
        for line in f:
            q, t, ql, tl, m, bl, qs, qe, strand, mapq, de = line.rstrip("\n").split("\t")
            if q == t:
                continue
            ql = int(ql); m = int(m); bl = int(bl); qs = int(qs); qe = int(qe)
            if bl == 0 or ql == 0:
                continue
            ident = m / bl; cov = (qe - qs) / ql
            if ident < ID_MIN or cov < COV_MIN:
                continue
            e = (q, t) if q < t else (t, q)
            if e not in best or ident > best[e][0]:
                best[e] = (ident, strand)
    return {e: s for e, (i, s) in best.items()}


def main():
    Hd, _ = dna_homology()
    idx = load_intron_index()
    strand = edge_strand()
    G, prov, part = build_family_graph(Hd, idx, filter_retrocopy=True)

    comps = sorted([c for c in nx.connected_components(G) if len(c) >= 2], key=len, reverse=True)
    n_minus = sum(1 for e in G.edges() if strand.get(tuple(sorted(e))) == "-")
    print(f"=== strand of post-retrocopy graph edges ===")
    print(f"  total edges {G.number_of_edges()} ; antisense('-') edges: {n_minus} "
          f"({100*n_minus/G.number_of_edges():.1f}%)")

    print("\n=== antisense-edge enrichment per big component ===")
    print(f"  {'comp':10} {'n':>4} {'edges':>6} {'minusE':>7} {'%minus':>7}")
    for i, c in enumerate(comps[:8]):
        sub = G.subgraph(c)
        me = sum(1 for e in sub.edges() if strand.get(tuple(sorted(e))) == "-")
        print(f"  comp-{i:<5} {len(c):>4} {sub.number_of_edges():>6} {me:>7} "
              f"{100*me/max(sub.number_of_edges(),1):>6.1f}%")

    # effect of dropping antisense edges
    named = [("RABL2A", "RABL2B"), ("APOBEC3D", "APOBEC3F"), ("RFPL1", "RFPL2"),
             ("RFPL2", "RFPL3"), ("GGT1", "GGTLC2")]
    def stats(drop_minus):
        H = nx.Graph()
        H.add_nodes_from(G.nodes())
        for a, b in G.edges():
            if drop_minus and strand.get(tuple(sorted((a, b)))) == "-":
                continue
            H.add_edge(a, b)
        comps2 = [c for c in nx.connected_components(H) if len(c) >= 2]
        sizes = sorted((len(c) for c in comps2), reverse=True)
        okn = sum(1 for a, b in named if a in H and b in H and nx.has_path(H, a, b))
        return len(comps2), sizes[:5], okn

    b = stats(False); d = stats(True)
    print("\n=== drop antisense('-') edges ===")
    print(f"  {'config':18} {'families':>9} {'maxcomp':>8} {'reals':>6}  top_sizes")
    print(f"  {'base':18} {b[0]:>9} {b[1][0]:>8} {b[2]:>4}/5  {b[1]}")
    print(f"  {'drop antisense':18} {d[0]:>9} {d[1][0]:>8} {d[2]:>4}/5  {d[1]}")
    print(f"\n  named real edges strand: " +
          ", ".join(f"{a}~{b}:{strand.get(tuple(sorted((a,b))),'NA')}" for a, b in named))


if __name__ == "__main__":
    main()
