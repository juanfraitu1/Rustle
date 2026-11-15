"""Exons of every catalog copy (bench/sdcn_units.py, exon part): the copy's self-hit in ref.paf split at introns, >= 100 bp,
autosomal, deduplicated by reference coordinates. usage: ex_units.py COPIES.tsv REF.paf SEX.tsv GENOME.fa OUTDIR"""
import re, sys, pysam
from common import sexset
copies, paf, sexf, fa, out = sys.argv[1:6]
SEX = {t for g, t in sexset(sexf) if g == 'ref'}
span = {}
for l in open(copies):
    if l.startswith('family_id'): continue
    f = l.split('\t'); span[f"{f[0]}|{f[1]}"] = (f[3], int(f[4]), int(f[5]))
ex = {}; seen = set()
for l in open(paf):
    f = l.rstrip('\n').split('\t'); q = f[0].split('|'); key = f"{q[0]}|{q[1]}"
    c, s, e = span[key]
    if f[5] != c or int(f[7]) >= e or int(f[8]) <= s: continue
    seen.add(key)
    cg = [x for x in f[12:] if x.startswith('cg:Z:')][0][5:]
    p = int(f[7]); a = p
    for n, o in re.findall(r'(\d+)([MIDN=X])', cg):
        n = int(n)
        if o == 'N':
            if p > a: ex.setdefault((c, a, p), q[0])
            p += n; a = p
        elif o in 'MD=X': p += n
    if p > a: ex.setdefault((c, a, p), q[0])
ex = {k: v for k, v in ex.items() if k[2] - k[1] >= 100 and k[0] not in SEX}
g = pysam.FastaFile(fa)
with open(f'{out}/exons.bed', 'w') as o, open(f'{out}/exons.fa', 'w') as fo:
    for (c, a, b), fam in sorted(ex.items()):
        o.write(f"{c}\t{a}\t{b}\tEX|{fam}|{c}:{a}-{b}\n"); fo.write(f">EX|{fam}|{c}:{a}-{b}\n{g.fetch(c, a, b).upper()}\n")
print('copies', len(span), 'with a self-hit', len(seen), 'exons', len(ex), 'families', len(set(ex.values())))
