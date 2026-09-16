#!/usr/bin/env python3
"""PC2 component and outside neighbours under two equal-bitscore HSP orders: shipped (protein_families.edges_from:
sorted(rows, reverse=True), ties broken by the remaining tuple fields) vs stable file order (verify_slim_protein/v_core.py)."""
import collections, csv, sys
sys.path.insert(0, "/mnt/c/Users/jfris/Desktop/Rustle/bench")
from protein_families import edges_from, pair_hsps
L = "/mnt/linuxdisk/home/juanfraitu/layer_order/npip_tbc1d3/light"
idx = {r["pid"]: r for r in csv.DictReader(open(f"{L}/work/P/proteins.index.tsv"), delimiter="\t")}
plen = {k: int(v["length_aa"]) for k, v in idx.items()}
searched = {x.strip() for x in open(f"{L}/work/P/searched.txt") if x.strip()}
hs = pair_hsps(f"{L}/work/P/blastp.tsv", plen)
E_ship = edges_from(hs, plen, 0.0)
E_stab = {}
for (q, s), rows in hs.items():
    use_q = plen[q] >= plen[s]
    acc, nid, aln = [], 0, 0
    for bits, n, l, q0, q1, s0, s1 in sorted(rows, key=lambda r: -r[0]):
        a, b = (q0, q1) if use_q else (s0, s1)
        if any(a < y and x < b for x, y in acc):
            continue
        acc.append((a, b)); nid += n; aln += l
    cov = sum(y - x for x, y in acc) / max(plen[q], plen[s])
    if cov >= 0.30:
        k = (min(q, s, key=int), max(q, s, key=int)); w = (nid / aln) * min(cov, 1.0)
        if w > E_stab.get(k, (0,))[0]:
            E_stab[k] = (w,)
pc2 = {p for p, v in idx.items() if v["name"] in {"TBC1D26", "USP6NL"} or (v["name"].startswith("TBC1D3") and v["name"][6:7] in ("", "B", "D", "E", "F", "G", "H", "I", "K") and len(v["name"]) <= 7)}
print("PC2 genes", sorted(idx[p]["name"] for p in pc2))
for tag, E in (("shipped edges_from", E_ship), ("stable file order (verifier)", E_stab)):
    adj = collections.defaultdict(set)
    for u, v in E:
        adj[u].add(v); adj[v].add(u)
    seen, st = set(pc2), list(pc2)
    while st:
        x = st.pop()
        for y in adj[x]:
            if y not in seen:
                seen.add(y); st.append(y)
    out = sorted({idx[y]["name"] for x in pc2 for y in adj[x]} - {idx[x]["name"] for x in pc2})
    print(f"{tag}: edges {len(E)}; PC2 component {len(seen)} genes ({len(seen - searched)} never searched); direct outside neighbours {len(out)} {out}")
