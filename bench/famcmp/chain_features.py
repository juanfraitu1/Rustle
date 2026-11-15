#!/usr/bin/env python3
"""PREREG r1619: per-chain features of our own models, with the RefSeq exact-chain label.

  chain_features.py REF_GTF OURS_GTF FASTA OUT_TSV
"""
import sys, re, gzip, bisect, collections, pysam
ref_p, ours_p, fa_p, out = sys.argv[1:5]
def load(path, attrs=False):
    ex = collections.defaultdict(list); at = {}
    for l in (gzip.open(path, 'rt') if path.endswith('.gz') else open(path)):
        f = l.split('\t')
        if len(f) < 9: continue
        k = (f[0], f[6], re.search(r'transcript_id "([^"]+)"', f[8]).group(1)) if f[2] in ('exon', 'transcript') else None
        if f[2] == 'exon': ex[k].append((int(f[3]) - 1, int(f[4])))
        elif f[2] == 'transcript' and attrs:
            at[k] = (float(re.search(r'reads "([^"]+)"', f[8]).group(1)), re.search(r'gene_id "([^"]+)"', f[8]).group(1))
    for e in ex.values(): e.sort()
    return ex, at
chain = lambda e: tuple((e[i][1], e[i + 1][0]) for i in range(len(e) - 1))  # 0-based half-open introns
ref = {(k[0], k[1], chain(e)) for k, e in load(ref_p)[0].items() if len(e) > 1}
ex, at = load(ours_p, True)
T = {}  # (chrom, strand, chain) -> best-supported model
for k, e in ex.items():
    if len(e) < 2: continue
    c = (k[0], k[1], chain(e)); r, g = at[k]
    if c not in T or r > T[c][0]: T[c] = (r, g, e)
isup = collections.Counter(); gtop = collections.Counter(); by_intron = collections.defaultdict(list)
for (c, s, ch), (r, g, e) in T.items():
    gtop[g] = max(gtop[g], r)
    for i in ch: isup[(c, s, i)] += r; by_intron[(c, s, i)].append((r, ch))
gitop = collections.Counter()
for (c, s, ch), (r, g, e) in T.items():
    for i in ch: gitop[g] = max(gitop[g], isup[(c, s, i)])
starts = collections.defaultdict(list)
for (c, s, i) in isup: starts[(c, s)].append(i)
for v in starts.values(): v.sort()
fa = pysam.FastaFile(fa_p); comp = str.maketrans('ACGT', 'TGCA')
CANON = {('GT', 'AG'), ('GC', 'AG'), ('AT', 'AC')}
def motif(c, s, i):
    d, a = fa.fetch(c, i[0], i[0] + 2).upper(), fa.fetch(c, i[1] - 2, i[1]).upper()
    return (d, a) if s == '+' else (a.translate(comp)[::-1], d.translate(comp)[::-1])
mcache = {}
sub = lambda a, b: any(b[k:k + len(a)] == a for k in range(len(b) - len(a) + 1))
with open(out, 'w') as w:
    w.write('chrom\tstrand\tin_ref\treads\tn_introns\tlocus_frac\tnoncanon\twobble\tsubchain\tretention\tmin_isup\tmin_isup_rel\ta_rich\tlength\n')
    for (c, s, ch), (r, g, e) in T.items():
        non = 0; wob = 0
        st = starts[(c, s)]
        for i in ch:
            if (c, s, i) not in mcache: mcache[(c, s, i)] = motif(c, s, i) in CANON
            non |= not mcache[(c, s, i)]
            for j in st[bisect.bisect_left(st, (i[0] - 10, 0)):bisect.bisect_right(st, (i[0] + 10, 1 << 60))]:
                if j != i and abs(j[1] - i[1]) <= 10 and isup[(c, s, j)] >= 5 * isup[(c, s, i)]: wob = 1
        sc = int(any(r2 > r and len(c2) > len(ch) and sub(ch, c2) for r2, c2 in by_intron[(c, s, ch[0])]))
        ret = 0
        for a, b in e:
            for j in st[bisect.bisect_left(st, (a, 0)):bisect.bisect_left(st, (b, 0))]:
                if j[1] <= b and isup[(c, s, j)] >= 5 * r: ret = 1
        mi = min(isup[(c, s, i)] for i in ch)
        L = fa.get_reference_length(c)
        tail = fa.fetch(c, e[-1][1], min(L, e[-1][1] + 20)).upper().count('A') if s == '+' else fa.fetch(c, max(0, e[0][0] - 20), e[0][0]).upper().count('T')
        w.write(f"{c}\t{s}\t{int((c, s, ch) in ref)}\t{r:g}\t{len(ch)}\t{r / gtop[g]:.3f}\t{int(non)}\t{wob}\t{sc}\t{ret}\t{mi:g}\t{mi / gitop[g]:.3f}\t{tail}\t{sum(b - a for a, b in e)}\n")
