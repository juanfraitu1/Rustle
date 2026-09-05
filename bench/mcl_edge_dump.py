#!/usr/bin/env python3
"""Write `<arm>/dump/e.edges.tsv` (chrom_i,start_i,end_i,chrom_j,start_j,end_j,identity) for an adapted MCL
arm so `soto_adjudicate.py`'s band table stratifies its asserted pairs by the best direct PAF identity between
the two copies' GENES.  usage: mcl_edge_dump.py <prefix> <paf> <arm_dir>"""
import sys, os, csv, collections
prefix, paf, arm = sys.argv[1:4]
members = [r for r in csv.DictReader(open(prefix + '.clusters.tsv'), delimiter='\t')]
gkey = lambda r: f"{r['chrom']}:{r['start']}-{r['end']}"
fam = {gkey(r): r['cluster_id'] for r in members}
ident = collections.defaultdict(float)
for l in open(paf):
    f = l.split('\t')
    if f[0] == f[5] or f[0] not in fam or fam.get(f[5]) != fam[f[0]]: continue
    p = tuple(sorted((f[0], f[5]))); ident[p] = max(ident[p], int(f[9]) / max(1, int(f[10])))
copies = list(csv.DictReader(open(f'{arm}/cat.copies.tsv'), delimiter='\t'))
# copy row -> gene key: a member of the same family whose span overlaps the copy
by_fam = collections.defaultdict(list)
for r in members: by_fam[r['cluster_id']].append(r)
gene_of = {}
for c in copies:
    s, e = int(c['start']), int(c['end'])
    for m in by_fam[c['family_id']]:
        if m['chrom'] == c['chrom'] and int(m['start']) - 1 < e and s < int(m['end']):
            gene_of[(c['chrom'], s, e)] = gkey(m); break
os.makedirs(f'{arm}/dump', exist_ok=True); n = 0
with open(f'{arm}/dump/e.edges.tsv', 'w') as o:
    o.write('chrom_i\tstart_i\tend_i\tchrom_j\tstart_j\tend_j\tidentity\n')
    keys = list(gene_of.items())
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            (ci, gi), (cj, gj) = keys[i], keys[j]
            p = tuple(sorted((gi, gj)))
            if p in ident:
                o.write(f'{ci[0]}\t{ci[1]}\t{ci[2]}\t{cj[0]}\t{cj[1]}\t{cj[2]}\t{ident[p]:.6f}\n'); n += 1
print(f'{arm}: {len(gene_of)}/{len(copies)} copies mapped to a gene, {n} direct edges')
