#!/usr/bin/env python3
"""family_def_immune_and_protein.py — two refinements to the protein-confirmed family definition.

TASK 1 (segment-aware recovery): immune-receptor segments (IG/TCR V/D/J/C) have NO representative
protein (not translated until VDJ recombination), so the protein bar drops them. They ARE cDNA-
homologous, but the shared V-DOMAIN over-merges segments ACROSS loci (one cDNA component spanned 6
chromosomes). Real immune families are TANDEM ARRAYS at ONE locus -> recover them from segment-biotype
cDNA homology with a SAME-LOCUS constraint (same chrom, within LOCUS_WIN), the inverse of the
disjoint-loci rule used for dispersed paralogs.

TASK 2 (protein-level community detection): check whether any large PROTEIN orthogroup is itself
domain-over-merged. Weight protein edges by reciprocal protein coverage and run weighted Louvain on
orthogroups >= MIN_OG. High-density orthogroups (cliques) should stay whole; only genuinely
domain-bridged ones split.

Run: /home/juanfra/miniforge3/bin/python bench/family_def_immune_and_protein.py
"""
import collections
import os
import re
import sys

import networkx as nx
from networkx.algorithms import community as nxcom

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from family_def_read_filters import dna_homology
from family_def_retrocopy_filter import build_family_graph, load_intron_index
from family_def_protein_validate import gene_biotype

PROT_AVA = "/home/juanfra/winloci_scratch/prot_ava.m8"
GENES_BED = "/home/juanfra/winloci_scratch/unmapped_poc/genes.bed"
LOCUS_WIN = 3_000_000           # tandem immune locus span
SEG_BT = ("V_segment", "C_region", "D_segment", "J_segment")
MIN_OG = 20
PID, PCOV, EVAL = 0.30, 0.50, 1e-5


def load_coord():
    coord = {}
    for line in open(GENES_BED):
        p = line.rstrip("\n").split("\t")
        if len(p) >= 4:
            coord[p[3]] = (p[0], int(p[1]), int(p[2]))
    return coord


def protein_edges_weighted():
    """gene-pair -> reciprocal-coverage weight (min(qcov,tcov)), passing the family bar."""
    coding = set()
    pe = {}
    with open(PROT_AVA) as f:
        for line in f:
            q, t, fid, qc, tc, ev, bits, al = line.rstrip("\n").split("\t")
            coding.add(q)
            if q == t:
                continue
            fid, qc, tc, ev = float(fid), float(qc), float(tc), float(ev)
            if fid >= PID and min(qc, tc) >= PCOV and ev <= EVAL:
                e = (q, t) if q < t else (t, q)
                w = min(qc, tc)
                if e not in pe or w > pe[e]:
                    pe[e] = w
    return coding, pe


def main():
    Hd, _ = dna_homology()
    idx = load_intron_index()
    G, _, _ = build_family_graph(Hd, idx, filter_retrocopy=True)
    bt = gene_biotype()
    coord = load_coord()
    coding, pe = protein_edges_weighted()

    # ---------------- TASK 1: locus-constrained immune-segment families ----------------
    print("=== TASK 1: immune-receptor families (segment cDNA homology + SAME-LOCUS constraint) ===")
    segset = {g for g, b in bt.items() if b in SEG_BT and g in G}
    # naive: cDNA components among segments (over-merge across loci)
    naive = [c for c in nx.connected_components(G.subgraph(segset)) if len(c) >= 2]
    # locus-constrained: keep an edge only if same chrom AND within LOCUS_WIN
    Gi = nx.Graph(); Gi.add_nodes_from(segset)
    for a, b in G.subgraph(segset).edges():
        ca, cb = coord.get(a), coord.get(b)
        if ca and cb and ca[0] == cb[0] and abs(ca[1] - cb[1]) <= LOCUS_WIN:
            Gi.add_edge(a, b)
    fams = [c for c in nx.connected_components(Gi) if len(c) >= 2]
    fams.sort(key=len, reverse=True)
    print(f"  segment genes in graph: {len(segset)}; naive cDNA families: {len(naive)} "
          f"(sizes {sorted((len(c) for c in naive), reverse=True)[:6]}, cross-locus over-merged)")
    print(f"  LOCUS-CONSTRAINED immune families: {len(fams)} "
          f"(sizes {[len(c) for c in fams][:8]})")
    # validate: each family single-chrom & compact + which locus
    single = sum(1 for c in fams if len({coord[g][0] for g in c if g in coord}) == 1)
    print(f"  families that are single-chromosome (compact tandem arrays): {single}/{len(fams)}")
    for c in fams[:6]:
        ch = collections.Counter(coord[g][0] for g in c if g in coord).most_common(1)[0][0]
        pos = [coord[g][1] for g in c if coord.get(g, (None,))[0] == ch]
        print(f"    n={len(c):2d}  {ch}:{min(pos)//1000}k-{max(pos)//1000}k  span={(max(pos)-min(pos))/1e6:.2f}Mb")

    # ---------------- TASK 2: protein-level community detection ----------------
    print("\n=== TASK 2: protein-level community detection on large protein orthogroups ===")
    Gp = nx.Graph()
    for (a, b), w in pe.items():
        Gp.add_edge(a, b, w=w)
    ogs = sorted([c for c in nx.connected_components(Gp) if len(c) >= 2], key=len, reverse=True)
    big = [c for c in ogs if len(c) >= MIN_OG]
    print(f"  protein orthogroups: {len(ogs)} (>= {MIN_OG} members: {len(big)})")
    split_count = 0
    for c in big[:10]:
        sub = Gp.subgraph(c)
        dens = nx.density(sub)
        coms = nxcom.louvain_communities(sub, weight="w", resolution=1.0, seed=0)
        coms = sorted((len(x) for x in coms), reverse=True)
        is_split = len([s for s in coms if s >= 2]) > 1
        split_count += is_split
        named = [g for g in c if not g.startswith("LOC")]
        root = collections.Counter(re.sub(r"[0-9]+[A-Z]*$", "", g) for g in named).most_common(1)
        tag = root[0][0] if root else "(all-LOC)"
        print(f"    OG n={len(c):2d} density={dens:.2f} [{tag}] -> communities {coms}  "
              f"{'SPLIT (domain-bridge?)' if is_split else 'intact (real family)'}")
    print(f"  => {len(big)-split_count}/{len(big)} large orthogroups are single real families; "
          f"{split_count} show internal sub-structure (candidate domain over-merge)")


if __name__ == "__main__":
    main()
