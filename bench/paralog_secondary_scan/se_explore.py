#!/usr/bin/env python3
"""Quick exploration: are SINGLE-EXON multi-copy families worth including?
The main scan skipped single-exon genes (no splice chain). Here we run an
overlap-based (no chain) variant of the dNM copy-anchor on single-exon genes
that sit in multi-mapping hotspots, to see how many are expressed and resolvable.

Usage: se_explore.py GENE_INTRONS.tsv GENE_META.tsv SECONDARY_POS.tsv MM_PLACEMENTS.tsv
"""
import sys, re, bisect, collections, statistics as st
MIN_SEC = 10; T = 2; KEXPR = 3; EXTENT = 0.7
def parse(c): return re.findall(r'(\d+)([MIDNSH=X])', c)
def alen(c): return sum(int(l) for l, o in parse(c) if o in 'M=XD') or 1

gi, meta, secf, mm = sys.argv[1:5]
genes = {}; se = []
for line in open(gi):
    if line.startswith('gene_id'): continue
    g, chrom, s, e, strand, ntx, chain = line.rstrip('\n').split('\t')
    if chain.strip(): continue            # single-exon only
    genes[g] = dict(chrom=chrom, start=int(s), end=int(e), strand=strand)
    se.append(g)
print(f"single-exon genes: {len(se)}", file=sys.stderr)
desc = {}
for line in open(meta):
    f = line.rstrip('\n').split('\t')
    if len(f) >= 4: desc[f[0]] = (f[2], f[3])
# secondary density candidate gen
sec = collections.defaultdict(list)
for line in open(secf):
    c, p = line.rstrip('\n').split('\t'); sec[c].append(int(p))
for c in sec: sec[c].sort()
cands = []
for g in se:
    G = genes[g]; arr = sec.get(G['chrom'])
    if not arr: continue
    n = bisect.bisect_right(arr, G['end']) - bisect.bisect_left(arr, G['start'])
    if n >= MIN_SEC: cands.append((n, g))
print(f"single-exon multi-mapping candidates (>={MIN_SEC} secondaries): {len(cands)}", file=sys.stderr)
# placements
reads = collections.defaultdict(list); idx = collections.defaultdict(list)
for line in open(mm):
    f = line.rstrip('\n').split('\t')
    if len(f) < 6 or f[5] == '' or f[4] == '*': continue
    qn, flag, chrom, pos, cig, nm = f
    reads[qn].append((chrom, int(pos), cig, int(nm), int(flag))); idx[chrom].append((int(pos), qn))
for c in idx: idx[c].sort()
def rin(chrom, a, b):
    arr = idx.get(chrom, []); lo = bisect.bisect_left(arr, (a - 50, '')); hi = bisect.bisect_right(arr, (b + 50, chr(0x10ffff)))
    return set(q for _, q in arr[lo:hi])

tax = collections.Counter(); examples = collections.defaultdict(list)
for n, g in cands:
    G = genes[g]; glen = G['end'] - G['start']; ga = sa = amb = 0
    for qn in rin(G['chrom'], G['start'], G['end']):
        gnm = None; galen = None
        for (chrom, pos, cig, nm, flag) in reads[qn]:
            if chrom != G['chrom']: continue
            st_ = '-' if flag & 16 else '+'
            if st_ != G['strand']: continue
            al = alen(cig); ov = min(pos + al, G['end']) - max(pos, G['start'])
            if ov < 0.5 * min(al, glen): continue          # must cover a real chunk
            if gnm is None or nm < gnm: gnm = nm; galen = al
        if gnm is None: continue
        best = None
        for (chrom, pos, cig, nm, flag) in reads[qn]:
            if chrom == G['chrom'] and G['start'] - 50 <= pos <= G['end']: continue
            if alen(cig) < EXTENT * galen: continue
            if best is None or nm < best: best = nm
        if best is None: ga += 1; continue
        d = best - gnm
        if d >= T: ga += 1
        elif d <= -T: sa += 1
        else: amb += 1
    tot = ga + sa + amb
    if tot < KEXPR: t = 'not_expressed'
    elif ga >= KEXPR: t = 'expressed_real_copy'
    elif ga == 0 and sa > 0: t = 'pure_spillover'
    elif amb >= max(sa, ga) and amb >= 3: t = 'non_identifiable'
    else: t = 'weak'
    tax[t] += 1
    if len(examples[t]) < 6:
        nm = desc.get(g, ('', ''))
        examples[t].append(f"{nm[0] or g}: {(nm[1] or '')[:38]} (G={ga},sib={sa},tie={amb})")
print("\n=== single-exon taxonomy ===")
for k in ['expressed_real_copy', 'non_identifiable', 'pure_spillover', 'weak', 'not_expressed']:
    print(f"  {k:20s} {tax.get(k,0)}")
print("\n=== examples ===")
for k in ['expressed_real_copy', 'non_identifiable', 'pure_spillover']:
    print(f"-- {k}:")
    for ex in examples[k]: print(f"     {ex}")
