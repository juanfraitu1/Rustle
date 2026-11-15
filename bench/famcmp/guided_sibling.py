#!/usr/bin/env python3
"""PREREG r1632 (docs/PREREG_sibling_guide_2026-11-11.md): guided_sibling.py REFSEQ DEFAULT_GTF GUIDED_GTF PROJECTED_GTF | selftest"""
import sys, re, gzip, collections

def models(lines):
    ex = collections.defaultdict(list)
    for l in lines:
        f = l.rstrip('\n').split('\t')
        if len(f) > 8 and f[2] == 'exon': ex[(f[0], f[6], re.search(r'transcript_id "([^"]*)"', f[8]).group(1))].append((int(f[3]) - 1, int(f[4])))
    return [(c, s, sorted(e)) for (c, s, _), e in ex.items()]
def chains(ms): return {(c, s, tuple((a[1], b[0]) for a, b in zip(e, e[1:]))) for c, s, e in ms if len(e) > 1}
def inside(ch, spans): return any(a < ch[2][-1][1] and ch[2][0][0] < b for a, b in spans.get(ch[:2], ()))
def score(ref, default, guided, proj):
    spans = collections.defaultdict(list)
    for c, s, e in proj: spans[(c, s)].append((e[0][0], e[-1][1]))
    d, g = chains(default), chains(guided); out = {}
    for name, x in (('new', g - d), ('lost', d - g)):
        i = [c for c in x if inside(c, spans)]
        out[name] = (len(i), sum(c in ref for c in i), len(x) - len(i))
    return out
def read(p): return models((gzip.open if p.endswith('.gz') else open)(p, 'rt'))

if sys.argv[1] == 'selftest':
    m = lambda c, t, exons: [f'{c}\tx\texon\t{a + 1}\t{b}\t.\t+\t.\ttranscript_id "{t}";' for a, b in exons]
    A, B, C = [(0, 10), (20, 30)], [(0, 10), (40, 50)], [(900, 910), (920, 930)]
    r = score(chains(models(m('c', 'r', A))), models(m('c', 'a', B) + m('c', 'z', C)), models(m('c', 'a', A)), models(m('c', 'p', [(0, 60)])))
    assert r == {'new': (1, 1, 0), 'lost': (1, 0, 1)}, r
    print('ok'); sys.exit()
r = score(chains(read(sys.argv[1])), read(sys.argv[2]), read(sys.argv[3]), read(sys.argv[4]))
(n, nc, no), (l, lc, lo) = r['new'], r['lost']
print(f'new chains at family loci {n}, correct {nc} = {nc / max(n, 1):.3f}; elsewhere {no}')
print(f'lost chains at family loci {l}, correct {lc} = {lc / max(l, 1):.3f}; elsewhere (partial-guide artefact) {lo}')
print(f'N1 {nc} >= 15: {"PASS" if nc >= 15 else "FAIL"}; N2 {nc / max(n, 1):.3f} >= 0.070: {"PASS" if n and nc / n >= 0.07 else "FAIL"}; L1 net {nc - lc} > 0 and lost share {lc / max(l, 1):.3f} <= 0.070: {"PASS" if nc > lc and lc <= 0.07 * l else "FAIL"}')
