#!/usr/bin/env python3
"""Pool `diff_model.py` outputs over a species' chromosomes and check claims D1-D6 of
docs/archive/2026-10/merged/PREREG_tool_difference_model_2026-10-11.md.
usage: score_claims.py SPECIES DIR [DIR ...]   (each DIR holds diff_StringTie.txt and diff_FLAIR.txt)"""
import sys, re, collections

def parse(path):
    """-> {section: {'n': total, 'ref': ref-exact count, cls: (n, ref-exact count)}}"""
    out = {}; sec = None
    for ln in open(path):
        m = re.match(r'(\S+)-only: ([\d,]+) chains, reference-exact ([\d.]+)', ln)
        if m:
            sec = m.group(1); n = int(m.group(2).replace(',', ''))
            out[sec] = {'n': n, 'ref': round(n * float(m.group(3)))}
            continue
        m = re.match(r'\s+(\S+)\s+([\d,]+)\s+share [\d.]+\s+reference-exact ([\d.]+)', ln)
        if m and sec:
            n = int(m.group(2).replace(',', '')); out[sec][m.group(1)] = (n, round(n * float(m.group(3))))
    return out

species, dirs = sys.argv[1], sys.argv[2:]
pool = {}
for t in ('StringTie', 'FLAIR'):
    acc = {}
    for d in dirs:
        for sec, v in parse(f"{d}/diff_{t}.txt").items():
            a = acc.setdefault(sec, collections.defaultdict(lambda: [0, 0]))
            a['_all'][0] += v['n']; a['_all'][1] += v['ref']
            for k, x in v.items():
                if isinstance(x, tuple):
                    a[k][0] += x[0]; a[k][1] += x[1]
    pool[t] = acc
rate = lambda p: p[1] / max(1, p[0])
st_only, st_ours = pool['StringTie']['StringTie'], pool['StringTie']['ours']
fl_only, fl_ours = pool['FLAIR']['FLAIR'], pool['FLAIR']['ours']
for t, (oth, ours) in (('StringTie', (st_only, st_ours)), ('FLAIR', (fl_only, fl_ours))):
    print(f"[{species}] {t}-only {oth['_all'][0]:,} (ref-exact {rate(oth['_all']):.3f}); ours-only {ours['_all'][0]:,} ({rate(ours['_all']):.3f})")
    for side, a in ((f'{t}-only', oth), ('ours-only', ours)):
        tot = a['_all'][0]
        print('   ' + side + ': ' + ', '.join(f"{k} {v[0]:,} ({v[0] / max(1, tot):.3f}; ref {rate(v):.3f})" for k, v in a.items() if k != '_all'))
st_cls = {k: v for k, v in st_only.items() if k != '_all'}
claims = {
    'D1 single_shared is the largest StringTie-only class': max(st_cls, key=lambda k: st_cls[k][0]) == 'single_shared',
    'D2 StringTie-only single_shared ref-exact > stitched': rate(st_only['single_shared']) > rate(st_only['stitched']),
    'D3 ours-only near_stronger ref-exact > StringTie-only overall': rate(st_ours['near_stronger']) > rate(st_only['_all']),
    'D4 FLAIR-only all read-carried; single_novel ref-exact < single_shared':
        fl_only['corrected'][0] + fl_only['stitched'][0] == 0 and rate(fl_only['single_novel']) < rate(fl_only['single_shared']),
    'D5 ours-only two_reads ref-exact > FLAIR-only overall': rate(fl_ours['two_reads']) > rate(fl_only['_all']),
    'D6 ours-only ref-exact > other-only, both tools':
        rate(st_ours['_all']) > rate(st_only['_all']) and rate(fl_ours['_all']) > rate(fl_only['_all']),
}
for k, v in claims.items():
    print(f"  {'PASS' if v else 'FAIL'}  {k}")
