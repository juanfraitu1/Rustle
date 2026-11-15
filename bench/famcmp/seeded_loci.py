#!/usr/bin/env python3
"""Loci the default assembly (tied secondary alignments seeded) has and the primary-only assembly, StringTie and FLAIR do not.

  seeded_loci.py REF_GTF DEFAULT_GTF PRIMARY_GTF STRINGTIE_GTF FLAIR_GTF OUT_TSV
A locus = merged same-strand overlapping transcripts of DEFAULT_GTF. "Not there" = no exon of the other GTF overlaps an exon of
the locus on the same strand ('.' strands match either). Descriptive.
"""
import sys, re, gzip, collections, bisect
def tx(path):
    ex = collections.defaultdict(list); nm = {}
    for l in (gzip.open(path, 'rt') if path.endswith('.gz') else open(path)):
        f = l.rstrip('\n').split('\t')
        if len(f) > 8 and f[2] == 'exon':
            k = (f[0], f[6], re.search(r'transcript_id "([^"]+)"', f[8]).group(1)); ex[k].append((int(f[3]) - 1, int(f[4])))
            g = re.search(r'gene "([^"]+)"', f[8]) or re.search(r'gene_id "([^"]+)"', f[8]); nm[k] = g.group(1) if g else ''
    for e in ex.values(): e.sort()
    return ex, nm
def index(ex):  # (chrom, strand) -> sorted exon list + starts
    d = collections.defaultdict(list)
    for (c, s, t), e in ex.items():
        for a, b in e: d[(c, s)].append((a, b, t))
    for v in d.values(): v.sort()
    return {k: (v, [x[0] for x in v], max(b - a for a, b, _ in v)) for k, v in d.items()}
def hits(ix, c, s, exons):
    out = set()
    for st in ((s,) if s in '+-' else ('+', '-')) + (('.',) if s in '+-' else ('.',)):
        if (c, st) not in ix: continue
        v, starts, mx = ix[(c, st)]
        for a, b in exons:
            i = bisect.bisect_left(starts, a - mx)
            while i < len(v) and v[i][0] < b:
                if v[i][1] > a: out.add(v[i][2])
                i += 1
    return out
ref, rn = tx(sys.argv[1]); dflt, _ = tx(sys.argv[2]); ix = [index(tx(p)[0]) for p in sys.argv[3:6]]; rix = index(ref)
chain = lambda e: tuple((e[i][1], e[i + 1][0]) for i in range(len(e) - 1))
refc = {(c, chain(e)) for (c, s, t), e in ref.items() if len(e) > 1}
by = collections.defaultdict(list)
for (c, s, t), e in dflt.items(): by[(c, s)].append(e)
loci = []
for (c, s), ts in by.items():
    ts.sort(key=lambda e: e[0][0]); cur = [ts[0]]; end = ts[0][-1][1]
    for e in ts[1:]:
        if e[0][0] < end: cur.append(e); end = max(end, e[-1][1])
        else: loci.append((c, s, cur)); cur = [e]; end = e[-1][1]
    loci.append((c, s, cur))
n = collections.Counter(); rows = []
for c, s, ts in loci:
    exons = sorted({x for e in ts for x in e}); p, st, fl = (bool(hits(i, c, s, exons)) for i in ix)
    n['loci'] += 1; n['not in StringTie and not in FLAIR'] += not st and not fl
    if p or st or fl: n['in primary-only run'] += p; continue
    rh = hits(rix, c, s, exons); genes = sorted({rn[(c, x, t)] for t in rh for x in '+-.' if (c, x, t) in rn})
    exact = any(len(e) > 1 and (c, chain(e)) in refc for e in ts)
    n['new with seeding, absent from StringTie and FLAIR'] += 1; n['  of those on an annotated gene'] += bool(rh); n['  with an exact RefSeq intron chain'] += exact; n['  spliced'] += any(len(e) > 1 for e in ts)
    rows.append((c, exons[0][0], exons[-1][1], s, len(ts), max(len(e) for e in ts), ','.join(genes)[:80], int(exact)))
with open(sys.argv[6], 'w') as w:
    w.write('chrom\tstart\tend\tstrand\ttranscripts\tmax_exons\trefseq_genes\texact_chain\n')
    for r in sorted(rows): w.write('\t'.join(map(str, r)) + '\n')
for k, v in n.items(): print(f'{k}\t{v}')
