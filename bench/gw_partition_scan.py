#!/usr/bin/env python3
"""Partition permutation test over every family in a gw pairs dump (from gw_subfamily_scan.py).
usage: gw_partition_scan.py <pairs.tsv> <out.tsv> [--min-pairs 10] [--B 200]"""
import sys, os, csv, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import modality as M
pairs_p, out_p = sys.argv[1:3]
MINP = int(sys.argv[sys.argv.index('--min-pairs') + 1]) if '--min-pairs' in sys.argv else 10
B = int(sys.argv[sys.argv.index('--B') + 1]) if '--B' in sys.argv else 200
fam = collections.defaultdict(list)
for r in csv.DictReader(open(pairs_p), delimiter='\t'):
    fam[r['family']].append((r['a'], r['b'], float(r['identity'])))
with open(out_p, 'w') as fh:
    fh.write('family\tn_pairs\tcontrast\tpartition_p\n')
    done = 0
    for f, pr in sorted(fam.items()):
        if len(pr) < MINP:
            continue
        p, obs = M.partition_perm_p(pr, B=B)
        fh.write(f'{f}\t{len(pr)}\t{obs:.4f}\t{p:.4f}\n'); fh.flush()
        done += 1
        if done % 25 == 0:
            print(f'  {done} families', file=sys.stderr)
print(f'wrote {out_p}: {done} families', file=sys.stderr)
