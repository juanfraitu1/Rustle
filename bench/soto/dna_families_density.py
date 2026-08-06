#!/usr/bin/env python3
"""DNA-only gene families: self-alignment edges + DENSITY-based clustering (no annotation, no reads).

Usage: dna_families_density.py <self_alignment.paf> <region_offset> [min_bp]

Three attempts have now been made to define families from genome self-alignment, and the first two failed
in the AGGREGATION rather than in the evidence:

  1. merge overlapping intervals into blocks   -> 86 mega-blocks at 6.3x gene size, best NPIP component
                                                  12% pure. Any-overlap merging chains transitively.
  2. reciprocal-overlap nodes + average linkage -> nodes are gene-scale (median 7.5 kb) but the cluster
                                                  holding NPIP still carries 138 other genes at every k.

The reason the second failed is measurable and is the point of this third attempt: 74% of alignments
starting inside an NPIP gene land on another NPIP gene, and 26% land elsewhere. Linkage merges on the
EXISTENCE of a connection, so a quarter of edges leaking outward collapses the family. The signal lives in
the PROPORTION of an node's connections, which only a weighted, density-aware method can use -- which is
also why the RNA pipeline already partitions E_r by gamma-quasi-clique rather than by components.

Edge weight is the number of supporting alignment records, so a node pair joined by many independent
duplications outweighs one joined by a single bridge. Greedy modularity is used as the density method
(the same objective family as the quasi-clique, and available without extra dependencies).

Genes are read only at the end, to score what the blind procedure produced.
"""
import gzip
import re
import sys
from collections import defaultdict

import networkx as nx

PAF = sys.argv[1]
OFF = int(sys.argv[2])
MIN_BP = int(sys.argv[3]) if len(sys.argv) > 3 else 2000
RECIP = 0.50

recs = []
for ln in open(PAF):
    f = ln.split("\t")
    qs, qe, ts, te = int(f[2]), int(f[3]), int(f[7]), int(f[8])
    if qe - qs < MIN_BP or te - ts < MIN_BP:
        continue
    de = next((float(x[5:]) for x in f[12:] if x.startswith("de:f:")), None)
    if de is None:
        continue
    recs.append((qs, qe, ts, te, 1.0 - de))

iv = sorted({(a, b) for (qs, qe, ts, te, _) in recs for (a, b) in ((qs, qe), (ts, te))})
parent = list(range(len(iv)))


def find(x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x


for i in range(len(iv)):
    a, b = iv[i]
    for j in range(i + 1, len(iv)):
        c, d = iv[j]
        if c >= b:
            break
        o = min(b, d) - max(a, c)
        if o >= RECIP * (b - a) and o >= RECIP * (d - c):
            ri, rj = find(i), find(j)
            if ri != rj:
                parent[ri] = rj
grp = defaultdict(list)
for i in range(len(iv)):
    grp[find(i)].append(i)
nodes = sorted((min(iv[i][0] for i in m), max(iv[i][1] for i in m)) for m in grp.values())
print(f"{len(recs)} records >= {MIN_BP} bp -> {len(nodes)} reciprocal-overlap nodes")


def node_of(a, b):
    best, bo = None, 0
    for k, (x, y) in enumerate(nodes):
        o = min(b, y) - max(a, x)
        if o > bo:
            best, bo = k, o
    return best


wt = defaultdict(float)
cnt = defaultdict(int)
for (qs, qe, ts, te, ident) in recs:
    u, v = node_of(qs, qe), node_of(ts, te)
    if u is None or v is None or u == v:
        continue
    k = (min(u, v), max(u, v))
    cnt[k] += 1
    wt[k] += ident * min(qe - qs, te - ts) / 1000.0   # weight by support AND aligned length

G = nx.Graph()
G.add_nodes_from(range(len(nodes)))
for (u, v), w in wt.items():
    G.add_edge(u, v, weight=w, n=cnt[(u, v)])
print(f"graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} weighted edges\n")

genes = []
with gzip.open("/mnt/c/Users/jfris/Desktop/Reference/chm13v2.0_RefSeq_full.gff.gz", "rt") as fh:
    for ln in fh:
        if ln[0] == "#":
            continue
        f = ln.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "gene" or f[0] != "chr16":
            continue
        m = re.search(r"Name=([^;]+)", f[8])
        s, e = int(f[3]) - 1 - OFF, int(f[4]) - OFF
        if m and 0 <= s <= 17_000_000:
            genes.append((s, e, m.group(1)))


def gat(a, b):
    return {g for (x, y, g) in genes if a < y and b > x}


npip = {g for _, _, g in genes if g.upper().startswith("NPIP")}


def score(comms, label):
    best = None
    for c in comms:
        gs = set()
        for i in c:
            gs |= gat(*nodes[i])
        hit = len(gs & npip)
        if best is None or hit > best[0]:
            best = (hit, gs, c)
    hit, gs, c = best
    print(f"{label:<38}{len(comms):>7}{f'{hit}/{len(npip)}':>10}"
          f"{len(gs)-hit:>13}{hit/max(len(gs),1):>9.2f}")
    return best


print(f"{'method':<38}{'#comm':>7}{'NPIP in':>10}{'other genes':>13}{'purity':>9}")
score(list(nx.connected_components(G)), "connected components (the old way)")
for res in (0.5, 1.0, 2.0, 4.0, 8.0):
    comms = nx.community.greedy_modularity_communities(G, weight="weight", resolution=res)
    best = score(comms, f"greedy modularity, resolution {res}")
    if res == 4.0:
        hit, gs, c = best
        named = sorted(g for g in gs if not g.startswith("LOC") and g not in npip)
        print(f"      NPIP found: {', '.join(sorted(gs & npip))}")
        print(f"      named non-NPIP in the same community: {', '.join(named[:12]) if named else '(none)'}")
