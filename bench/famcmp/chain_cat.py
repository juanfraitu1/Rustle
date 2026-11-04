#!/usr/bin/env python3
"""r1524/r1526 chain bars: base spliced chains kept in the arm, and the arm's NEW spliced chains' exact match to a CAT Liftoff
transcript chain.  usage: chain_cat.py BASE.gtf ARM.gtf CAT.gff3[.gz]   (chain = (chrom, strand, introns))"""
import sys, re, gzip, collections
def chains(path, feat='exon', tid=r'transcript_id "([^"]+)"', parent=None):
    ex = collections.defaultdict(list); strand = {}
    op = gzip.open if path.endswith('.gz') else open
    for l in op(path, 'rt'):
        if l.startswith('#'): continue
        f = l.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != feat: continue
        m = re.search(parent or tid, f[8])
        if not m: continue
        ex[m.group(1)].append((int(f[3]) - 1, int(f[4]))); strand[m.group(1)] = (f[0], f[6])
    out = {}
    for t, e in ex.items():
        e.sort(); ch = tuple((a[1], b[0]) for a, b in zip(e, e[1:]))
        if ch: out[t] = strand[t] + (ch,)
    return out
base, arm, cat = sys.argv[1:4]
B = chains(base); A = chains(arm)
C = {v for v in chains(cat, parent=r'Parent=([^;]+)').values()}
bset = set(B.values()); aset = set(A.values())
kept = len(bset & aset); new = aset - bset
m = sum(1 for c in new if c in C)
print(f'transcripts base {len(B)} arm {len(A)}; base spliced chains {len(bset)} kept {kept} = {kept/max(1,len(bset)):.4f} (bar >= 0.99)')
print(f'new chains {len(new)}, exact CAT match {m} = {m/max(1,len(new)):.3f} (bar >= 0.05)')
