#!/usr/bin/env python3
"""PREREG r1631 (docs/PREREG_tie_groups_transfer_2026-11-11.md): do tie-group models add RefSeq intron chains the default lacks?

  tie_chains.py REFSEQ_GTF[.gz] TIE_GROUPS_GTF DEFAULT_GTF [DEFAULT_GTF ...]     (selftest: tie_chains.py selftest)

New chain = a multi-exon tie-group model's (contig, strand, introns) in no DEFAULT_GTF. Correct = a RefSeq chain."""
import sys, re, gzip, collections

def chains(lines):
    ex = collections.defaultdict(list); loci = {}
    for l in lines:
        f = l.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'exon': continue
        k = (f[0], f[6], re.search(r'transcript_id "([^"]*)"', f[8]).group(1))
        ex[k].append((int(f[3]) - 1, int(f[4])))
        m = re.search(r'group_loci "(\d+)"', f[8])
        if m: loci[k] = int(m.group(1))
    out = {}
    for k, e in ex.items():
        e.sort(); c = tuple((a[1], b[0]) for a, b in zip(e, e[1:]))
        if c: out[(k[0], k[1], c)] = max(out.get((k[0], k[1], c), 0), loci.get(k, 0))
    return out

def score(ref, tg, default):
    rows = []
    for name, sel in (('all', tg), ('groups of >= 2 loci', {c: n for c, n in tg.items() if n >= 2})):
        new = [c for c in sel if c not in default]; ok = sum(c in ref for c in new)
        old = [c for c in sel if c in default]
        rows.append((name, len(sel), len(new), ok, len(old), sum(c in ref for c in old)))
    return rows

def read(p): return chains((gzip.open if p.endswith('.gz') else open)(p, 'rt'))

if sys.argv[1] == 'selftest':
    g = lambda c, s, t, exons, n=None: [f'{c}\tx\texon\t{a + 1}\t{b}\t.\t{s}\t.\ttranscript_id "{t}";' + (f' group_loci "{n}";' if n else '') for a, b in exons]
    ref = chains(g('c', '+', 'r1', [(0, 10), (20, 30)]) + g('c', '+', 'r2', [(0, 10), (50, 60)]))
    tg = chains(g('c', '+', 'a', [(5, 10), (20, 40)], 3) + g('c', '+', 'b', [(0, 10), (50, 60)], 1) + g('c', '-', 'd', [(0, 10), (20, 30)], 2) + g('c', '+', 'm', [(0, 99)], 2))
    de = chains(g('c', '+', 'x', [(0, 10), (50, 60)]))
    assert score(ref, tg, de) == [('all', 3, 2, 1, 1, 1), ('groups of >= 2 loci', 2, 2, 1, 0, 0)], score(ref, tg, de)
    print('ok'); sys.exit()
ref = read(sys.argv[1]); tg = read(sys.argv[2]); default = {}
for p in sys.argv[3:]: default.update(read(p))
print('set\tchains\tnew\tnew_correct\tshare\tin_default\tin_default_correct')
rows = score(ref, tg, default)
for n, a, b, c, d, e in rows: print(f'{n}\t{a}\t{b}\t{c}\t{c / max(b, 1):.3f}\t{d}\t{e}')
_, _, b, c, _, _ = rows[0]
print(f'A1 new correct {c} (bar >= 15): {"PASS" if c >= 15 else "FAIL"}; A2 share {c / max(b, 1):.3f} (bar >= 0.070): {"PASS" if b and c / b >= 0.070 else "FAIL"}')
