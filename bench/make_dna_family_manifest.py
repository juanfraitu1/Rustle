#!/usr/bin/env python3
"""make_dna_family_manifest.py — DNA-family manifest for rustle --family-manifest.
Builds the family-manifest TSV (family_id, gene_id, chrom, start, end, strand) from the
genome cDNA all-vs-all paralog families. NOTE: the cDNA components are OVER-MERGED (max 389
members) -- de-bridge before production use (see family_def_scaffold_wiring.md).
Run: python bench/make_dna_family_manifest.py > /home/juanfra/winloci_scratch/dna_family_manifest.tsv
"""
import sys
import networkx as nx
sys.path.insert(0, "bench")
from family_def_read_filters import dna_homology

coord = {}
for line in open("/home/juanfra/winloci_scratch/unmapped_poc/genes.bed"):
    p = line.rstrip("\n").split("\t")
    if len(p) >= 4:
        coord[p[3]] = (p[0], int(p[1]), int(p[2]))
Hd, _ = dna_homology()
DG = nx.Graph()
for (ga, gb), r in Hd.items():
    if r.get("id", 0) >= 0.90 and max(r["cov_a"], r["cov_b"]) >= 0.30:
        DG.add_edge(ga, gb)
print("family_id\tgene_id\tchrom\tstart\tend\tstrand")
for i, fam in enumerate(sorted((c for c in nx.connected_components(DG) if len(c) >= 2), key=lambda c: -len(c))):
    members = [g for g in sorted(fam) if g in coord]
    if len(members) >= 2:
        for g in members:
            c, s, e = coord[g]
            print(f"DNAFAM_{i}\t{g}\t{c}\t{s}\t{e}\t.")
