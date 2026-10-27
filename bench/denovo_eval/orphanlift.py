"""Orphan lift edges for PREREG DENOVOFAM-L5 (docs/LEDGER_O1.md).

  orphanlift.py GRAPH_TSV CLUSTERS_TSV LIFT_EDGES OUT_E6 OUT_E7

GRAPH_TSV = baseline `mcl_families --dump-graph`; CLUSTERS_TSV = baseline `.clusters.tsv`; LIFT_EDGES = liftchain.py output.
E6 = lifts with both loci orphans (no baseline edge). E7 = E6 + orphan-to-family lifts with w >= the family's weakest internal edge.
"""
import collections
import sys

graph, clusters, lift, out6, out7 = sys.argv[1:6]
cl = {}
for line in open(clusters):
    f = line.rstrip('\n').split('\t')
    if f[0] != 'cluster_id':
        cl[f'{f[5]}:{f[6]}-{f[7]}'] = f[0]
node = set()
weakest = {}
for line in open(graph):
    a, b, w = line.split('\t')
    node.update((a, b))
    if a in cl and cl.get(a) == cl.get(b):
        weakest[cl[a]] = min(weakest.get(cl[a], 1.0), float(w))
n = collections.Counter()
with open(out6, 'w') as o6, open(out7, 'w') as o7:
    for line in open(lift):
        a, b, w = line.rstrip('\n').split('\t')
        oa, ob = a not in node, b not in node
        if oa and ob:
            o6.write(line); o7.write(line); n['orphan-orphan'] += 1
        elif oa or ob:
            fam = cl.get(b if oa else a)
            if fam is not None and fam in weakest and float(w) >= weakest[fam]:
                o7.write(line); n['orphan-family kept'] += 1
            else:
                n['orphan-family dropped'] += 1
print(dict(n), file=sys.stderr)
