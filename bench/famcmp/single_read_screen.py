#!/usr/bin/env python3
"""Screen (exploratory): which extra requirement makes single-read chains worth accepting?

  single_read_screen.py REF_GTF BAM FASTA DEFAULT_GTF FLOOR1_GTF OUT_TSV
One row per chain the floor-1 run adds over the default run: RefSeq match and candidate features.
"""
import sys, re, gzip, collections, pysam
def models(path):
    ex = collections.defaultdict(list)
    for l in (gzip.open(path, 'rt') if path.endswith('.gz') else open(path)):
        f = l.split('\t')
        if len(f) > 8 and f[2] == 'exon': ex[(f[0], f[6], re.search(r'transcript_id "([^"]+)"', f[8]).group(1))].append((int(f[3]), int(f[4])))
    out = {}
    for (c, s, _), e in ex.items():
        e.sort()
        if len(e) > 1: out.setdefault((c, s, tuple((e[i][1], e[i + 1][0] - 1) for i in range(len(e) - 1))), []).append((e[0][0], e[-1][1]))
    return out
ref_p, bam_p, fa_p, base_p, f1_p, out_p = sys.argv[1:7]
ref = models(ref_p); base = models(base_p); f1 = models(f1_p); add = [k for k in f1 if k not in base]
two = set()  # default chains carried by exactly 2 reads: the comparison class
for l in open(base_p):
    f = l.split('\t')
    if len(f) > 8 and f[2] == 'transcript' and 'reads "2"' in f[8]: two.add((f[0], f[6], re.search(r'transcript_id "([^"]+)"', f[8]).group(1)))
ex2 = collections.defaultdict(list)
for l in open(base_p):
    f = l.split('\t')
    if len(f) > 8 and f[2] == 'exon':
        k = (f[0], f[6], re.search(r'transcript_id "([^"]+)"', f[8]).group(1))
        if k in two: ex2[k].append((int(f[3]), int(f[4])))
two = set()
for (c, s, _), e in ex2.items():
    e.sort()
    if len(e) > 1: two.add((c, s, tuple((e[i][1], e[i + 1][0] - 1) for i in range(len(e) - 1))))
rows = [('added', k) for k in add] + [('default2', k) for k in two]
known = {(c, i) for c, s, ch in base for i in ch}
ends = collections.defaultdict(lambda: (set(), set()))  # (chrom, strand) -> left ends, right ends of default models, 50 bp bins
for (c, s, ch), se in base.items():
    for a, b in se: ends[(c, s)][0].add(a // 50); ends[(c, s)][1].add(b // 50)
near = lambda S, x: any(x // 50 + d in S for d in (-1, 0, 1))
wi = {(c, i) for _, (c, s, ch) in rows for i in ch}; wp = {(c, ch[k], ch[k + 1]) for _, (c, s, ch) in rows for k in range(len(ch) - 1)}
isup = collections.Counter(); psup = collections.Counter()
for a in pysam.AlignmentFile(bam_p):
    if a.is_unmapped or a.is_supplementary or a.is_secondary or 'N' not in a.cigarstring: continue
    c = a.reference_name; ch = []; p = a.reference_start
    for op, n in a.cigartuples:
        if op == 3: ch.append((p, p + n))
        if op in (0, 2, 3, 7, 8): p += n
    for k, i in enumerate(ch):
        if (c, i) in wi: isup[(c, i)] += 1
        if k and (c, ch[k - 1], i) in wp: psup[(c, ch[k - 1], i)] += 1
fa = pysam.FastaFile(fa_p); canon = {'GTAG', 'CTAC', 'GCAG', 'CTGC', 'ATAC', 'GTAT'}
with open(out_p, 'w') as w:
    w.write('set\tchrom\tstrand\tn_introns\tin_ref\tnew_introns\tmin_intron_reads\tmin_pair_reads\tcanonical\tstart_known\tend_known\n')
    for st, (c, s, ch) in rows:
        mot = all((fa.fetch(c, a, a + 2) + fa.fetch(c, b - 2, b)).upper() in canon for a, b in ch)
        a, b = (f1 if st == 'added' else base)[(c, s, ch)][0]; L, R = ends[(c, s)]
        w.write(f"{st}\t{c}\t{s}\t{len(ch)}\t{int((c, s, ch) in ref)}\t{sum((c, i) not in known for i in ch)}\t{min(isup[(c, i)] for i in ch)}\t"
                f"{min((psup[(c, ch[k], ch[k + 1])] for k in range(len(ch) - 1)), default=-1)}\t{int(mot)}\t{int(near(L, a))}\t{int(near(R, b))}\n")
