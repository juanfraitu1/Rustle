#!/usr/bin/env python3
"""r1518 bars 4-5 (docs/PREREG_seed_unspliced_primary_2026-11-03.md): A.gtf B.gtf CAT.gff3.gz [contig,...]
bar 4: CAT single-exon genes (every transcript one exon) overlapped by >= 1 same-strand single-exon model, per arm.
bar 5: >= 2-exon transcripts identical in (chrom, start, end, introns) between the arms, as a share of A's."""
import gzip, re, sys
from collections import defaultdict

def models(path, contigs):
    ex = defaultdict(list); strand = {}
    for l in open(path):
        f = l.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'exon' or (contigs and f[0] not in contigs): continue
        t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
        ex[(f[0], t)].append((int(f[3]), int(f[4]))); strand[(f[0], t)] = f[6]
    mono, spliced = [], set()
    for (c, t), e in ex.items():
        e.sort()
        if len(e) == 1: mono.append((c, strand[(c, t)], e[0][0], e[0][1]))
        else: spliced.add((c, strand[(c, t)], e[0][0], e[-1][1], tuple((a[1], b[0]) for a, b in zip(e, e[1:]))))
    return mono, spliced

def cat_mono(path, contigs):
    tex = defaultdict(int); tg = {}; gpos = {}
    for l in gzip.open(path, 'rt'):
        if l.startswith('#'): continue
        f = l.rstrip('\n').split('\t')
        if len(f) < 9 or (contigs and f[0] not in contigs): continue
        a = dict(kv.split('=', 1) for kv in f[8].split(';') if '=' in kv)
        if f[2] == 'exon': tex[a['Parent']] += 1
        elif f[2] == 'transcript': tg[a['ID']] = a.get('Parent') or a.get('gene')
        elif f[2] == 'gene': gpos[a['ID']] = (f[0], f[6], int(f[3]), int(f[4]))
    bad = {g for t, g in tg.items() if tex.get(t, 0) != 1}
    return {g: p for g, p in gpos.items() if g in set(tg.values()) and g not in bad}

def hit(genes, mono, tight=False):
    # tight (Amendment 1): the model's span must be <= 2x the gene's, so one array-wide blob does not "hit" every gene it covers
    by = defaultdict(list)
    for c, s, a, b in mono: by[(c, s)].append((a, b))
    n = 0
    for g, (c, s, a, b) in genes.items():
        if any(x <= b and a <= y and (not tight or y - x <= 2 * (b - a)) for x, y in by[(c, s)]): n += 1
    return n

A, B, CAT = sys.argv[1:4]; contigs = set(sys.argv[4].split(',')) if len(sys.argv) > 4 else None
ma, sa = models(A, contigs); mb, sb = models(B, contigs); genes = cat_mono(CAT, contigs)
ha, hb = hit(genes, ma), hit(genes, mb)
print(f"bar4 as registered: CAT single-exon genes {len(genes)}: hit A {ha} B {hb} ratio {hb/max(ha,1):.4f} (bar >= 0.98); single-exon models A {len(ma)} B {len(mb)}")
ta, tb = hit(genes, ma, True), hit(genes, mb, True)
print(f"bar4 Amendment 1 (model span <= 2x gene): hit A {ta} B {tb} ratio {tb/max(ta,1):.4f} (bar >= 0.98)")
print(f"bar5 spliced identical {len(sa & sb)} of A {len(sa)} = {len(sa & sb)/max(len(sa),1):.4f} (bar >= 0.99); B spliced {len(sb)}")
