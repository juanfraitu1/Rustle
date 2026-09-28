#!/usr/bin/env python3
"""Readthrough junctions: what minimap2 -G 50k would remove vs cost, and an annotation-free 'polyA site inside
the intron' rule. Report: docs/READTHROUGH_G50K_AND_LAST_EXON_2026-09-25.md
usage: readthrough_rules.py TAG BAM GFF CHROM OUTDIR [FASTA]"""
import sys, os, re, bisect, collections, pickle
import pysam
tag, bam, gff, chrom, out = sys.argv[1:6]
fasta = sys.argv[6] if len(sys.argv) > 6 else None
P = lambda s: os.path.join(out, f'{tag}.{s}')

# ---- annotation: genes (minus readthrough-described records), exon unions, annotated introns ----
genes, par, tx_ex, rt_names = {}, {}, collections.defaultdict(list), []
for ln in open(gff):
    if ln[0] == '#':
        continue
    r = ln.rstrip('\n').split('\t')
    if len(r) < 9 or r[0] != chrom:
        continue
    at = dict(kv.split('=', 1) for kv in r[8].split(';') if '=' in kv)
    if r[2] in ('gene', 'pseudogene'):
        if 'readthrough' in at.get('description', '').lower():
            rt_names.append(at.get('Name')); continue
        genes[at['ID']] = (int(r[3]), int(r[4]), r[6], at.get('Name', at['ID']))
    elif r[2] == 'exon' and 'Parent' in at:
        tx_ex[at['Parent']].append((int(r[3]), int(r[4])))
    elif 'Parent' in at and 'ID' in at:
        par[at['ID']] = at['Parent']
exon_iv, ann_intron = [], {}
for t, xs in tx_ex.items():
    g = t
    while g not in genes and g in par:
        g = par[g]
    if g not in genes:
        continue
    xs.sort()
    for a, b in xs:
        exon_iv.append((a, b, g))
    for (a1, b1), (a2, b2) in zip(xs, xs[1:]):
        ann_intron[(b1 + 1, a2 - 1)] = g
BW = 10000
EXB = collections.defaultdict(list)
for iv in exon_iv:
    for k in range(iv[0] // BW, iv[1] // BW + 1):
        EXB[k].append(iv)


def genes_at(p):
    return {g for s, e, g in EXB.get(p // BW, ()) if s <= p <= e}


# ---- reads ----
junc = collections.Counter()           # (s, e, strand) intron 1-based inclusive
first_donor = {}                       # (strand, 5' end, 3' end, read) -> donor base of the read's FIRST intron
ends = collections.defaultdict(list)   # strand -> 3' ends of spliced reads
jreads = collections.defaultdict(list)
for rd in pysam.AlignmentFile(bam).fetch(chrom):
    if rd.flag & 2308:
        continue
    ts = rd.get_tag('ts') if rd.has_tag('ts') else '+'
    st = '+' if (ts == '+') != rd.is_reverse else '-'
    pos, introns = rd.reference_start, []
    for op, L in rd.cigartuples:
        if op == 3:
            introns.append((pos + 1, pos + L)); pos += L
        elif op in (0, 2, 7, 8):
            pos += L
    if not introns:
        continue
    e3, e5 = (rd.reference_end, rd.reference_start + 1) if st == '+' else (rd.reference_start + 1, rd.reference_end)
    ends[st].append((e3, e5))
    first_donor[(st, e5, e3, rd.query_name)] = introns[0][0] - 1 if st == '+' else introns[-1][1] + 1
    for s, e in introns:
        junc[(s, e, st)] += 1
        jreads[(s, e, st)].append(rd.query_name)
for st in ends:
    ends[st].sort()

# ---- 5' START density: clusters of spliced-read 5' ends (gap > 100 bp opens a new one; >= 3 reads = real start)
starts = {}
for st, E in ends.items():
    S5 = sorted((e5, e3) for e3, e5 in E)
    FD = collections.defaultdict(list)
    for (st2, e5, e3, _), d in first_donor.items():
        if st2 == st:
            FD[(e5, e3)].append(d)
    lab, grp = [], []
    def close(g):
        for i in g:
            lab.append(len(g) >= 3)
    for i, (e5, _) in enumerate(S5):
        if grp and e5 - S5[grp[-1]][0] > 100:
            close(grp); grp = []
        grp.append(i)
    close(grp)
    starts[st] = (S5, lab, FD)


def motif_ok(s, e, st):
    if FA is None:
        return 'NA'
    d, a = FA.fetch(chrom, s - 1, s + 1).upper(), FA.fetch(chrom, e - 2, e).upper()
    pairs = {('GT', 'AG'), ('GC', 'AG'), ('AT', 'AC')} if st == '+' else {('CT', 'AC'), ('CT', 'GC'), ('GT', 'AT')}
    return (d, a) in pairs


# ---- 3' end DENSITY: clusters of spliced-read 3' ends (gap > 25 bp opens a new one) ----
FA = pysam.FastaFile(fasta) if fasta else None
COMP = str.maketrans('ACGTN', 'TGCAN')
real_end = {}      # (strand, end position) -> (is_real_end, cluster mode, has_PAS)
clusters = {}
for st, E in ends.items():
    grp = []
    for i, (e3, _) in enumerate(E):
        if grp and e3 - E[grp[-1]][0] > 25:
            clusters.setdefault(st, []).append(grp); grp = []
        grp.append(i)
    if grp:
        clusters.setdefault(st, []).append(grp)
    for grp in clusters.get(st, []):
        mode = collections.Counter(E[i][0] for i in grp).most_common(1)[0][0]
        ok, pas = len(grp) >= 3, False
        if FA is not None:
            if st == '+':
                down = FA.fetch(chrom, mode, mode + 20).upper(); up = FA.fetch(chrom, max(0, mode - 50), mode).upper()
            else:
                down = FA.fetch(chrom, max(0, mode - 21), mode - 1).upper().translate(COMP)[::-1]
                up = FA.fetch(chrom, mode - 1, mode + 49).upper().translate(COMP)[::-1]
            ok = ok and down.count('A') < 12
            pas = 'AATAAA' in up or 'ATTAAA' in up
        for i in grp:
            real_end[(st, E[i][0])] = (ok, mode, pas, len(grp))

rows = []
for (s, e, st), n in junc.items():
    if n < 2:
        continue
    L = e - s + 1
    if (s, e) in ann_intron and genes[ann_intron[(s, e)]][2] == st:
        cls = 'ANN_ge50k' if L >= 50000 else 'ANN_lt50k'; ga = gb = genes[ann_intron[(s, e)]][3]
    else:
        don, acc = (s - 1, e + 1) if st == '+' else (e + 1, s - 1)
        A = {g for g in genes_at(don) if genes[g][2] == st}
        B = {g for g in genes_at(acc) if genes[g][2] == st}
        cls = None
        for a in A:
            for b in B:
                if a != b and (genes[a][1] < genes[b][0] or genes[b][1] < genes[a][0]) and not (A & B):
                    cls = 'RT'; ga, gb = genes[a][3], genes[b][3]
        if cls is None:
            continue
    E = ends[st]
    lo, hi = bisect.bisect_left(E, (s, 0)), bisect.bisect_left(E, (e + 1, 0))
    T = hi - lo                                                   # R_any: any spliced read's 3' end inside the intron
    # R_up (refined on dev, frozen before held-out): only reads that START upstream of the donor, i.e. transcripts
    # of the donor's gene that terminate inside the intron (they polyadenylate where J keeps going)
    up = [e3 for e3, e5 in E[lo:hi] if (e5 < s if st == '+' else e5 > e)]
    U = len(up)
    Up = [x for x in up if real_end[(st, x)][0]]
    pk = collections.Counter(real_end[(st, x)][1] for x in Up)
    top = pk.most_common(1)[0][0] if pk else None
    top_pas = next((real_end[(st, x)][2] for x in Up if real_end[(st, x)][1] == top), 'NA') if top else 'NA'
    S5, lab, FD = starts[st]
    if st == '+':      # 5' inside the intron, 3' beyond the acceptor (e)
        a, b = bisect.bisect_left(S5, (s, 0)), bisect.bisect_left(S5, (e + 1, 0))
        idx = [i for i in range(a, b) if S5[i][1] > e]
    else:              # minus: 5' = high coordinate inside the intron, 3' below the acceptor (s)
        a, b = bisect.bisect_left(S5, (s, 0)), bisect.bisect_left(S5, (e + 1, 0))
        idx = [i for i in range(a, b) if S5[i][1] < s]
    V, Vp = len(idx), sum(lab[i] for i in idx)
    # strict: the read's own first exon lies inside J's intron (its first donor is inside the intron)
    Vs = sum(1 for i in idx if lab[i] and any(s <= d <= e for d in FD[S5[i]]))
    rows.append((cls, s, e, st, L, n, T, U / (U + n), ga, gb, len(Up) / (len(Up) + n), top_pas,
                 V, Vp / (Vp + n), motif_ok(s, e, st), Vs / (Vs + n)))
with open(P('junctions.tsv'), 'w') as f:
    f.write('class\tstart\tend\tstrand\tintron_len\treads\tends_inside\tR\tgene_a\tgene_b\tR_peak\ttop_peak_PAS\tV\tQ\tcanonical\tQ1\n')
    for r in rows:
        f.write('\t'.join(map(str, r[:7])) + f'\t{r[7]:.4f}\t{r[8]}\t{r[9]}\t{r[10]:.4f}\t{r[11]}\t{r[12]}\t{r[13]:.4f}\t{r[14]}\t{r[15]:.4f}\n')
pickle.dump({(r[1], r[2], r[3]): jreads[(r[1], r[2], r[3])] for r in rows if r[0] in ('RT', 'ANN_ge50k')},
            open(P('jreads.pkl'), 'wb'))
c = collections.Counter(r[0] for r in rows)
print(f'[{tag}] junctions scored: {dict(c)}; readthrough-described records removed: {len(rt_names)}', file=sys.stderr)
