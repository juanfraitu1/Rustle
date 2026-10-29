#!/usr/bin/env python3
"""EDGE-HUMAN-2 readout: verdict of every U1 lift edge (keyed to locus coordinates) and what the family floor keeps / drops.

  edge_verdicts.py LOCI_GFF3 FLOORS_JSON EDGES_U1_TSV EDGES_FLOOR_TSV [--species gorilla]
"""
import collections, json, sys
from exon_rule_judge import Judge

gff, fl, u1, flo = sys.argv[1:5]
J = Judge(sys.argv[6] if len(sys.argv) > 6 else 'gorilla', gff)
F = json.load(open(fl))
kept = {tuple(l.split('\t')[:2]) for l in open(flo)}
order = ['true-recent', 'true-simiiformes', 'older', 'hgnc', 'same-symbol', 'WRONG', 'same-gene', 'unjudged']
c = collections.defaultdict(collections.Counter)
dropped_true, true_added = [], []
for l in open(u1):
    a, b, w, ident = l.rstrip('\n').split('\t')[:4]
    v = J.verdict(a, b)
    k = 'kept' if (a, b) in kept else 'dropped'
    c[k][v] += 1
    fa, fb = F['member'].get(a), F['member'].get(b)
    fs = [F['family'][f] for f in (fa, fb) if f]
    floor = min(fs) if fs else F['global_floor']
    if v.startswith('true'):
        true_added.append((J.name(a), J.name(b), v, float(ident), floor, k))
        if k == 'dropped':
            dropped_true.append((J.name(a), J.name(b), v, float(ident), floor, fa, fb))
print('set\tn\t' + '\t'.join(order))
for k in ('kept', 'dropped'):
    print(f'{k}\t{sum(c[k].values())}\t' + '\t'.join(str(c[k][v]) for v in order))
print('\n# true edges (either level) with identity, floor applied, kept/dropped:')
for r in sorted(true_added, key=lambda r: r[3]):
    print('\t'.join(str(x) if not isinstance(x, float) else f'{x:.3f}' for x in r))
print('\n# families whose floor dropped a true edge:')
for r in dropped_true:
    print(f'{r[0]}/{r[1]}\t{r[2]}\tid {r[3]:.3f} < floor {r[4]:.3f}\tfamilies {r[5]} {r[6]}')
