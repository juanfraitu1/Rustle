#!/usr/bin/env python3
"""r1524 bar 1: of the assemblable misses (assemblable.tsv rows with ours = 0 and maxP_chain >= 2), how many gain a model on the
gene's strand overlapping its exons (any exon count) in each GTF.  usage: score34.py ASSEMBLABLE.tsv NAME=GTF ..."""
import sys, os, re, csv, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'denovo_eval'))
import famtruth as F
G = F.load_cat()
one = '--one-chain' in sys.argv  # r1526: the 97 members with reads but no two identical chains
rows = [r for r in csv.DictReader(open(sys.argv[1]), delimiter='\t') if r['ours'] == '0' and
        ((int(r['maxP_chain']) < 2 and int(r['maxT_chain']) < 2 and int(r['n_local_spliced_P']) >= 1) if one else int(r['maxP_chain']) >= 2)]
bar = 25 if one else 20
for arg in sys.argv[2:]:
    if arg == '--one-chain': continue
    name, path = arg.split('=', 1)
    ex = collections.defaultdict(list)
    for l in open(path):
        f = l.split('\t')
        if len(f) > 8 and f[2] == 'exon': ex[(f[0], f[6])].append((int(f[3]) - 1, int(f[4])))
    byc = {k: sorted(v) for k, v in ex.items()}
    hit = []
    for r in rows:
        v = G[r['gene']]; xs = byc.get((v['chrom'], v['strand']), [])
        if any(a < y and x < b for a, b in xs for x, y in v['exons']): hit.append(r['name'])
    print(f"[bar1] {name}: {len(hit)} of {len(rows)} assemblable misses gain a same-strand model (bar >= {bar}): {' '.join(hit)}")
