#!/usr/bin/env python3
"""PREREG r1602 (docs/PREREG_hard_members_signal_2026-11-10.md): decisive / tied / weaker molecules per hard member.

  hard_signal.py HARD_TSV(r1601) CAT_TSV BAM AS_TABLE OUT_TSV
"""
import sys, collections, pysam

hard_tsv, cat_tsv, bam_path, as_table, out = sys.argv[1:6]
with open(hard_tsv) as fh:
    hdr = fh.readline().rstrip('\n').split('\t'); H = [dict(zip(hdr, l.rstrip('\n').split('\t'))) for l in fh]
G = {}
for l in open(cat_tsv):
    f = l.rstrip('\n').split('\t')
    G[f[0]] = (f[4], f[7], [tuple(map(int, b.split('-'))) for b in f[9].split(',')]) if f[0] != 'gene_id' else None
bam = pysam.AlignmentFile(bam_path)
recs = {}  # gene -> {molecule: best AS of its records here}
for r in H:
    chrom, strand, ex = G[r['gene']]; m = {}
    for a in bam.fetch(chrom, ex[0][0], ex[-1][1]):
        if a.is_supplementary or a.is_unmapped or a.is_reverse != (strand == '-'): continue
        if not any(x < e and s < y for x, y in a.get_blocks() for s, e in ex): continue
        s = a.get_tag('AS') if a.has_tag('AS') else 0
        m[a.query_name] = max(m.get(a.query_name, s), s)
    recs[r['gene']] = m
need = set(n for m in recs.values() for n in m); best = {}
with open(as_table) as fh:
    for l in fh:
        if l[0] == '#': continue
        f = l.split('\t', 4)
        if f[0] in need: best[f[0]] = (int(f[1]), int(f[2]))
cnt = {}
for r in H:
    c = collections.Counter()
    for n, s in recs[r['gene']].items():
        b, sec = best.get(n, (s, -1))
        c['decisive' if s == b and sec < b else 'tied' if s == b else 'weaker'] += 1
    cnt[r['gene']] = c
grp = lambda r: 'L3+' if int(r.get('a1', 0)) >= 3 else 'L1-2' if int(r.get('a1', 0)) >= 1 else 'missed'  # no levels column: all 'missed'
cls = lambda c: 'distinguishable' if c['decisive'] >= 2 else 'weak' if c['decisive'] == 1 else 'indistinguishable'
with open(out, 'w') as w:
    w.write('family\tgene\tname\tclass\tr1601_group\tdecisive\ttied\tweaker\tsignal\n')
    for r in H:
        c = cnt[r['gene']]
        w.write(f"{r['family']}\t{r['gene']}\t{r['name']}\t{r['class']}\t{grp(r)}\t{c['decisive']}\t{c['tied']}\t{c['weaker']}\t{cls(c)}\n")
tab = collections.Counter((grp(r), cls(cnt[r['gene']])) for r in H)
for g in ('L3+', 'L1-2', 'missed'):
    print(g, {k: tab[(g, k)] for k in ('distinguishable', 'weak', 'indistinguishable')})
for g in ('L1-2', 'missed'):
    for k in ('distinguishable', 'weak', 'indistinguishable'):
        print(g, k, collections.Counter((r['name'].rstrip('0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ') or r['name'], r['class']) for r in H if grp(r) == g and cls(cnt[r['gene']]) == k).most_common(12))
