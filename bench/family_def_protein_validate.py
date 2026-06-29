#!/usr/bin/env python3
"""family_def_protein_validate.py — protein-level GROUND TRUTH for the cDNA family de-bridge.

The cDNA-homology over-merge could never be adjudicated from inside (0 trustworthy positives in the
big blobs). Protein homology is the orthogonal validator the project flagged as needed:
  - a real protein-coding paralog FAMILY requires members that are PROTEIN-CODING and PROTEIN-
    HOMOLOGOUS. Non-coding (lncRNA) cDNA-homology bridges have NO protein -> excluded by construction.
  - protein RECIPROCAL coverage (whole-protein homology) separates real paralogs from shared-domain
    bridges more robustly than cDNA (protein is under selection).

Inputs: proteins.fa (1 representative protein per gene, gene-name keyed) + mmseqs all-vs-all
(prot_ava.m8: query target fident qcov tcov evalue bits alnlen). A protein edge (gene pair) passes
if fident>=PID and min(qcov,tcov)>=PCOV and evalue<=EVAL.

Validations:
  (a) named real families are protein-homologous (sanity);
  (b) coding composition + #protein-orthogroups per big cDNA component (over-merge confirmation);
  (c) within-community vs cross-(cut) cDNA edges: protein-confirmed rate (does cov_min de-bridge
      track protein truth?);
  (d) protein-confirmed family graph: keep only cDNA edges that are coding+protein-homologous ->
      family count / max comp / real-family preservation (the "real claim").

Run: /home/juanfra/miniforge3/bin/python bench/family_def_protein_validate.py
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
from family_def_community_debridge import load_weights

PROT_FA = "/home/juanfra/winloci_scratch/proteins.fa"
PROT_AVA = "/home/juanfra/winloci_scratch/prot_ava.m8"
GFF = "/mnt/c/Users/jfris/Desktop/GGO_genomic.gff"
PID, PCOV, EVAL = 0.30, 0.50, 1e-5

NAMED_REAL = [("RABL2A", "RABL2B"), ("APOBEC3D", "APOBEC3F"),
              ("RFPL1", "RFPL2"), ("RFPL2", "RFPL3"), ("GGT1", "GGTLC2")]


def load_protein_edges():
    coding = set()
    for line in open(PROT_FA):
        if line.startswith(">"):
            coding.add(line[1:].strip().split()[0])
    pe = {}
    with open(PROT_AVA) as f:
        for line in f:
            q, t, fident, qcov, tcov, ev, bits, aln = line.rstrip("\n").split("\t")
            if q == t:
                continue
            fident, qcov, tcov, ev = float(fident), float(qcov), float(tcov), float(ev)
            if fident >= PID and min(qcov, tcov) >= PCOV and ev <= EVAL:
                e = (q, t) if q < t else (t, q)
                key = (fident, min(qcov, tcov))
                if e not in pe or key > pe[e]:
                    pe[e] = key
    return coding, pe


def gene_biotype():
    bt = {}
    for line in open(GFF):
        if "\tgene\t" in line or "\tpseudogene\t" in line:
            m = re.search(r"Name=([^;]+)", line)
            if m:
                b = re.search(r"gene_biotype=([^;]+)", line)
                bt[m.group(1)] = b.group(1).strip() if b else "?"
    return bt


def main():
    if not os.path.exists(PROT_AVA):
        print(f"[!] {PROT_AVA} not ready yet"); sys.exit(1)
    Hd, _ = dna_homology()
    idx = load_intron_index()
    G, _, _ = build_family_graph(Hd, idx, filter_retrocopy=True)
    coding, pe = load_protein_edges()
    bt = gene_biotype()
    print(f"=== protein truth: {len(coding)} coding genes, {len(pe)} protein-homology gene-pairs "
          f"(fident>={PID}, min-cov>={PCOV}, e<={EVAL}) ===")

    def pedge(a, b):
        return tuple(sorted((a, b))) in pe

    # (a) named real families protein-homologous?
    print("\n(a) named real families protein-homologous?")
    for a, b in NAMED_REAL:
        print(f"    {a}~{b}: {'PROTEIN-HOMOLOGOUS' if pedge(a, b) else ('one/both non-coding' if a not in coding or b not in coding else 'NOT (sub-bar)')}")

    # (b) coding composition + #protein-orthogroups per big cDNA component
    comps = sorted([c for c in nx.connected_components(G) if len(c) >= 2], key=len, reverse=True)
    print("\n(b) big cDNA components: coding% and protein-orthogroup count (over-merge confirmation)")
    for i, c in enumerate([c for c in comps if len(c) >= 20][:8]):
        pc = sum(1 for g in c if bt.get(g) == "protein_coding")
        # protein orthogroups WITHIN this component
        sub = nx.Graph()
        sub.add_nodes_from(g for g in c if g in coding)
        for a in c:
            for b in c:
                if a < b and pedge(a, b):
                    sub.add_edge(a, b)
        ogs = [x for x in nx.connected_components(sub) if len(x) >= 2]
        print(f"    comp-{i} n={len(c):3d}: coding={pc:3d}({100*pc/len(c):3.0f}%)  "
              f"protein-orthogroups(>=2)={len(ogs):3d}  largest={max((len(x) for x in ogs), default=0)}")

    # (c) does cov_min de-bridge track protein truth? within-community vs cut edges
    W = load_weights("cov_min", Hd)
    final = []
    for c in comps:
        if len(c) < 20:
            final.append(set(c)); continue
        H = G.subgraph(c).copy()
        for a, b in H.edges():
            H[a][b]["w"] = W.get(tuple(sorted((a, b))), 0.02)
        final.extend([set(x) for x in nxcom.louvain_communities(H, weight="w", resolution=1.0, seed=0) if len(x) >= 2])
    node2fam = {g: i for i, c in enumerate(final) for g in c}
    within = cut = within_pc = cut_pc = 0
    for a, b in G.edges():
        if a not in coding or b not in coding:
            continue                                   # protein truth only defined for coding pairs
        same = node2fam.get(a) is not None and node2fam.get(a) == node2fam.get(b)
        ok = pedge(a, b)
        if same:
            within += 1; within_pc += ok
        else:
            cut += 1; cut_pc += ok
    print("\n(c) cov_min de-bridge vs protein truth (coding-coding cDNA edges only):")
    print(f"    WITHIN-community edges: {within_pc}/{within} protein-confirmed ({100*within_pc/max(within,1):.0f}%)")
    print(f"    CUT (cross-community)  : {cut_pc}/{cut} protein-confirmed ({100*cut_pc/max(cut,1):.0f}%)")
    print(f"    => de-bridge {'AGREES with' if within_pc/max(within,1) > 2*cut_pc/max(cut,1) else 'is weakly aligned to'} protein truth")

    # (d) protein-confirmed family graph (the 'real claim'): keep cDNA edges that are coding+protein-homologous
    Gp = nx.Graph()
    Gp.add_nodes_from(G.nodes())
    kept = 0
    for a, b in G.edges():
        if pedge(a, b):
            Gp.add_edge(a, b); kept += 1
    Gp.remove_nodes_from([n for n in list(Gp.nodes()) if Gp.degree(n) == 0])
    pc_comps = [c for c in nx.connected_components(Gp) if len(c) >= 2]
    sizes = sorted((len(c) for c in pc_comps), reverse=True)
    okn = sum(1 for a, b in NAMED_REAL if a in Gp and b in Gp and nx.has_path(Gp, a, b))
    print("\n(d) PROTEIN-CONFIRMED family graph (coding + protein-homologous cDNA edges only):")
    print(f"    cDNA edges {G.number_of_edges()} -> protein-confirmed {kept}; "
          f"families {len(pc_comps)}, max comp {sizes[0] if sizes else 0}, top {sizes[:5]}")
    print(f"    named real families connected: {okn}/{len(NAMED_REAL)}")
    print(f"    (vs cov_min cDNA de-bridge: 1283 families, max comp 43)")


if __name__ == "__main__":
    main()
