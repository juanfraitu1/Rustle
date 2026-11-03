#!/usr/bin/env python3
"""Cut audit (r1516/r1517 bar 3): of the cuts in a units list, how many have BOTH flanking exon runs overlapping exons of the SAME
RefSeq coding gene and of no other (a cut inside one gene). usage: cut_audit.py UNITS.tsv ASSEMBLY.gtf REFSEQ.gtf.gz"""
import sys, re, gzip, collections, bisect
units, gtf, rs = sys.argv[1:4]
A = lambda k, s: (re.search(k + r' "([^"]+)"', s) or [None, ''])[1]
ex = collections.defaultdict(list)
for l in open(gtf):
    f = l.rstrip('\n').split('\t')
    if len(f) >= 9 and f[2] == 'exon': ex[A('transcript_id', f[8])].append((int(f[3]), int(f[4])))
genes = collections.defaultdict(list); coding = set()
for l in gzip.open(rs, 'rt'):
    f = l.rstrip('\n').split('\t')
    if len(f) < 9 or f[2] != 'exon': continue
    g = A('gene_name', f[8]) or A('gene_id', f[8])
    if A('transcript_id', f[8]).startswith(('rna-NM_', 'rna-XM_')): coding.add(g)
    genes[(f[0], f[6])].append((int(f[3]), int(f[4]), g))
for k in genes: genes[k].sort()
starts = {k: [v[0] for v in lst] for k, lst in genes.items()}
def hits(chrom, strand, a, b):
    lst = genes.get((chrom, strand), []); i = bisect.bisect_right(starts.get((chrom, strand), []), b)
    return {g for s, e, g in lst[max(0, i - 5000):i] if e >= a and g in coding}
n = inside = 0; by = collections.Counter()
for l in open(units):
    f = l.rstrip('\n').split('\t')
    if f[0] == 'tid' or len(f) < 4: continue
    tid, chrom, strand = f[0], f[1], f[2]; exs = sorted(ex[tid])
    for j in f[3].split(';'):
        s, e = map(int, j.split('-')); n += 1
        left = set().union(*[hits(chrom, strand, a, b) for a, b in exs if b < s] or [set()])
        right = set().union(*[hits(chrom, strand, a, b) for a, b in exs if a > e] or [set()])
        if left and left == right and len(left) == 1: inside += 1; by['inside_one_gene'] += 1
        elif left and right and left.isdisjoint(right): by['between_two_genes'] += 1
        elif left and right: by['overlapping_gene_sets'] += 1
        else: by['a_side_without_coding_gene'] += 1
print(f"cuts {n}; inside one RefSeq coding gene {inside} ({100 * inside / max(1, n):.1f}%); {dict(by)}")
