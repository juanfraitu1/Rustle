#!/usr/bin/env python3
"""Build junction probes for the cross-species replication test (PREREG md5 b0794d03).
A probe = last 150 bp of the upstream exon + first 150 bp of the downstream exon: the sequence that exists
only if the junction is spliced. Three classes: test (read-through junctions), positive control (ordinary
introns of the same units), negative control (scrambled halves).
usage: xspecies_junction_probes.py <readthrough.tsv> <units.tsv> <genome.fa> <out.fa>
"""
import sys, csv, subprocess, random
rt_p, units_p, fa, out_p = sys.argv[1:5]
W = 150
random.seed(1337)
def seq(c, s, e):
    if s < 1: return ''
    r = subprocess.run(f"samtools faidx {fa} {c}:{s}-{e}", shell=True, capture_output=True, text=True).stdout
    return ''.join(r.split('\n')[1:]).upper()
def probe(c, d, a):
    """d = last base of the upstream exon (1-based), a = first base of the downstream exon"""
    up, dn = seq(c, d - W + 1, d), seq(c, a, a + W - 1)
    return up, dn
test = []
for r in csv.DictReader(open(rt_p), delimiter='\t'):
    c, js, je = r['chrom'], int(r['intron_start']), int(r['intron_end'])
    up, dn = probe(c, js, je + 1)          # intron is 0-based half-open: last exon base js, next exon base je+1
    if len(up) == W and len(dn) == W: test.append((f"TEST_{r['from_family']}_{r['to_family']}_{js}", up, dn))
pos = []
for u in csv.DictReader(open(units_p), delimiter='\t'):
    if u['member_status'] not in ('kept_full', 'kept_trimmed'): continue
    ex = [tuple(map(int, x.split('-'))) for x in u['exons'].split(',')]
    for i in range(len(ex) - 1):
        d, a = ex[i][1], ex[i + 1][0]
        if a - d < 200: continue
        up, dn = probe(u['chrom'], d, a + 1)
        if len(up) == W and len(dn) == W: pos.append((f"POS_{u['family_id']}_{u['copy_idx']}_{d}", up, dn))
        if len(pos) >= 120: break
    if len(pos) >= 120: break
neg = []
for i, (n, up, _) in enumerate(test):
    j = (i + 7) % len(test)
    neg.append((f"NEG_{i}_{j}", up, test[j][2]))
with open(out_p, 'w') as o:
    for name, up, dn in test + pos + neg:
        o.write(f">{name}\n{up}{dn}\n")
print(f"probes: {len(test)} test, {len(pos)} positive control, {len(neg)} negative control -> {out_p}")
