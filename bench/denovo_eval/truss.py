#!/usr/bin/env python3
"""k-truss families from a sketchfam edge list (2026-10-29, user: "is there a more lenient graph theory concept" between cliques and
components): keep an edge only while it lies in >= k-2 triangles (iterated to a fixed point), families = connected components of the
surviving edges (>= 2 loci). A lone edge (adjacent gene, retrocopy in a host, a single linker exon) has no triangle and is dropped;
nested families keep their edges because every member has common neighbours.
  truss.py EDGES_TSV LOCI_GFF3 OUT [--k 3] [--mode truss|bridge]   (EDGES_TSV = sketchfam OUT.edges.tsv) -> OUT.fam.copies.tsv
--mode bridge: instead of the truss, drop every BRIDGE (edge whose removal disconnects its component) whose two endpoints both have
other neighbours, i.e. a lone edge joining two groups; a standalone pair (most Simiiformes truth pairs are two-gene families, which a
triangle rule destroys) is kept because its endpoints have no other neighbour. Iterated to a fixed point."""
import sys, collections
edges_f, loci_f, out = sys.argv[1:4]; K = int(sys.argv[sys.argv.index('--k') + 1]) if '--k' in sys.argv else 3
nbr = collections.defaultdict(set)
for l in open(edges_f):
    x, y = l.split('\t')[:2]; nbr[x].add(y); nbr[y].add(x)
MODE = sys.argv[sys.argv.index('--mode') + 1] if '--mode' in sys.argv else 'truss'
sys.setrecursionlimit(1000000)
def bridges(adj):
    disc = {}; low = {}; out = []; t = [0]
    def dfs(u, p):
        disc[u] = low[u] = t[0]; t[0] += 1
        for v in adj[u]:
            if v == p: continue
            if v not in disc:
                dfs(v, u); low[u] = min(low[u], low[v])
                if low[v] > disc[u]: out.append((u, v))
            else: low[u] = min(low[u], disc[v])
    for u in list(adj):
        if u not in disc and adj[u]: dfs(u, None)
    return out
changed = True; dropped = 0
while changed and MODE == 'bridge':
    changed = False
    for x, y in bridges(nbr):
        if len(nbr[x]) > 1 and len(nbr[y]) > 1:
            nbr[x].discard(y); nbr[y].discard(x); changed = True; dropped += 1
while changed:
    changed = False
    for x in list(nbr):
        for y in list(nbr[x]):
            if x < y and len(nbr[x] & nbr[y]) < K - 2:
                nbr[x].discard(y); nbr[y].discard(x); changed = True; dropped += 1
parent = {}
def find(a):
    parent.setdefault(a, a)
    while parent[a] != a: parent[a] = parent[parent[a]]; a = parent[a]
    return a
for x in nbr:
    for y in nbr[x]: parent[find(x)] = find(y)
comp = collections.defaultdict(list)
for x in nbr:
    if nbr[x]: comp[find(x)].append(x)
fams = sorted([v for v in comp.values() if len(v) >= 2], key=lambda v: -len(v))
loci = {}; exons = collections.defaultdict(list)
for l in open(loci_f):
    if l[0] == '#': continue
    f = l.rstrip('\n').split('\t')
    if len(f) < 9 or f[2] != 'exon': continue
    a = dict(x.split('=', 1) for x in f[8].split(';') if '=' in x); t = a.get('gene') or a.get('Parent', '').replace('gene-', '')
    loci[t] = (f[0], f[6]); exons[t].append((int(f[3]), int(f[4])))
with open(out + '.fam.copies.tsv', 'w') as w:
    w.write("family_id\tcopy_idx\ttid\tchrom\tstart\tend\tn_exon\tstrand\tn_reads\texons\n")
    for i, v in enumerate(fams):
        for j, t in enumerate(sorted(v)):
            c, s = loci[t]; ex = sorted(exons[t])
            w.write(f"TR{i}\t{j}\t{t}\t{c}\t{ex[0][0]}\t{ex[-1][1]}\t{len(ex)}\t{s}\t0\t{','.join(f'{a}-{b}' for a, b in ex)}\n")
print(f"k {K}: edges dropped {dropped}, families {len(fams)}, members {sum(len(v) for v in fams)}, largest {[len(v) for v in fams[:5]]}", file=sys.stderr)
