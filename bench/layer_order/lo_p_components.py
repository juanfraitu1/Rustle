#!/usr/bin/env python3
"""P interpretation check: is each member P cluster an exact connected component of the §6ko graph over all saved blastp
searches (light/work/P/blastp.tsv), or an MCL cut inside a larger component? Also: the shipped §6ko dev table family that
holds the TBC1D3 genes (o1_falsemerge/lit/pfam_dev/r2_refseq.families.tsv)."""
import collections, csv, sys
sys.path.insert(0, "/mnt/c/Users/jfris/Desktop/Rustle/bench")
from protein_families import edges_from, pair_hsps
L = "/mnt/linuxdisk/home/juanfraitu/layer_order/npip_tbc1d3/light"
idx = {r["pid"]: r for r in csv.DictReader(open(f"{L}/work/P/proteins.index.tsv"), delimiter="\t")}
plen = {k: int(v["length_aa"]) for k, v in idx.items()}
E = edges_from(pair_hsps(f"{L}/work/P/blastp.tsv", plen), plen, 0.0)
searched = {x.strip() for x in open(f"{L}/work/P/searched.txt") if x.strip()}
adj = collections.defaultdict(set)
for u, v in E:
    adj[u].add(v); adj[v].add(u)
g2p = {r["gene_id"]: p for p, r in idx.items()}
groups = collections.defaultdict(set)
for r in csv.DictReader(open(f"{L}/P.groups.corrected.tsv"), delimiter="\t"):
    if r["p_status"] == "in_P":  # member clusters only (P|other rows are coding U genes outside them)
        groups[r["group_id"]].add(g2p[r["gene_id"]])
print(f"§6ko edges over the saved searches (bench/protein_families.pair_hsps + edges_from, the shipped greedy-by-bitscore "
      f"HSP order): {len(E)}; searched proteins {len(searched)}")
for gid, G in sorted(groups.items()):
    seen, stack = set(G), list(G)
    while stack:
        x = stack.pop()
        for y in adj[x]:
            if y not in seen:
                seen.add(y); stack.append(y)
    out = sorted({idx[y]["name"] for x in G for y in adj[x]} - {idx[x]["name"] for x in G})
    print(f"{gid}: {len(G)} genes; connected component containing it: {len(seen)} genes ({sum(1 for x in seen if x not in searched)} "
          f"never searched); direct §6ko neighbours outside the group: {len(out)} {out[:15]}")
fam = collections.defaultdict(list)
for r in csv.DictReader(open("/mnt/linuxdisk/home/juanfraitu/o1_falsemerge/lit/pfam_dev/r2_refseq.families.tsv"), delimiter="\t"):
    fam[r["family_id"]].append(r["name"])
for f, ns in fam.items():
    if any(n.startswith("TBC1D3") for n in ns) or any(n.startswith("NPIP") for n in ns):
        print(f"pfam_dev r2_refseq {f}: {sorted(ns)}")
