#!/usr/bin/env python3
"""Row 1682: SEDEF identity of the holding SD pair, inside the SD-held low-SHARE cell (docs/PREREG_sd_identity_cell_2026-11-13.md).

  edgesd.py REL.tsv SEDEF.bed [REL.tsv SEDEF.bed ...]      edgesd.py --selftest
"""
import sys, collections
BIN, THR = 100_000, (0.90, 0.95, 0.98, 0.99)


def sd_index(lines):
    ix = collections.defaultdict(list); col = None
    for l in lines:
        f = l.split('\t')
        if len(f) < 22 or not f[1].isdigit(): continue
        if col is None: col = 20 if '.' in f[20] else 21   # fracMatch: col 21, col 22 in the HG002 build; fixed per file from its first row
        fm = float(f[col]); assert 0.4 <= fm <= 1, l[:200]
        for c, s, e, c2, s2, e2 in ((f[0], int(f[1]), int(f[2]), f[3], int(f[4]), int(f[5])), (f[3], int(f[4]), int(f[5]), f[0], int(f[1]), int(f[2]))):
            for b in range(s // BIN, e // BIN + 1): ix[(c, b)].append((s, e, c2, s2, e2, fm))
    return ix


def sd_id(ix, a, b, minov=300):
    """highest fracMatch among SD pairs with a mate over each locus; None when none holds them"""
    v = [fm for k in range(a[1] // BIN, a[2] // BIN + 1) for s, e, c2, s2, e2, fm in ix.get((a[0], k), ())
         if min(e, a[2]) - max(s, a[1]) >= minov and c2 == b[0] and min(e2, b[2]) - max(s2, b[1]) >= minov]
    return max(v) if v else None


def loc(n):
    c, r = n.rsplit(':', 1); s, e = r.split('-'); return c, int(s), int(e)


def selftest():
    row = lambda fm: '\t'.join(['c1', '1000', '5000', 'c2', '100', '4100'] + ['x'] * 14 + [fm, '0.5'])
    ix = sd_index([row('0.93'), row('0.99'), 'c1\tx\n'])
    assert sd_id(ix, ('c1', 2000, 3000), ('c2', 500, 900)) == 0.99 and sd_id(ix, ('c2', 500, 900), ('c1', 2000, 3000)) == 0.99
    assert sd_id(ix, ('c1', 2000, 3000), ('c3', 500, 900)) is None and sd_id(ix, ('c1', 4800, 9000), ('c2', 500, 900)) is None
    print('selftest ok')


def main(a):
    for rel, bed in zip(a[::2], a[1::2]):
        ix = sd_index(open(bed)); rows = [l.rstrip('\n').split('\t') for l in open(rel)]; h = {x: i for i, x in enumerate(rows[0])}
        cell = [(sd_id(ix, loc(r[h['locus_a']]), loc(r[h['locus_b']])), float(r[h['identity']]), r[h['truth']] == '1', r[h['rel_kind']])
                for r in rows[1:] if r[h['linked']] == '0' and r[h['stage']] == 'low_share' and r[h['sd_held']] == '1']
        lost = sum(c[0] is None for c in cell); cell = [c for c in cell if c[0] is not None]
        print(f'{rel}: cell {len(cell)} pairs, truth {sum(c[2] for c in cell)} (no holding pair found again: {lost})')
        byid = sorted(cell, key=lambda c: -c[1])
        for t in THR:
            p = [c for c in cell if c[0] >= t]; tr = sum(c[2] for c in p); ct = sum(c[2] for c in byid[:len(p)])
            kinds = collections.Counter(c[3] for c in p if not c[2])
            print(f'  sd_id >= {t}: picked {len(p)} true {tr} precision {tr / max(1, len(p)):.3f} | control (PAF identity, same picks) true {ct} precision {ct / max(1, len(p)):.3f} | false by kind {dict(kinds)}')


if __name__ == '__main__':
    selftest() if sys.argv[1:] == ['--selftest'] else main(sys.argv[1:])
