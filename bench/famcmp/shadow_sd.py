#!/usr/bin/env python3
"""Row 1683: SD mate with primary support as the mark of a shadow locus (docs/PREREG_shadow_sd_mate_2026-11-13.md).

  shadow_sd.py LOCI.tsv ASSEMBLY.gtf SEDEF.bed      shadow_sd.py --selftest
"""
import sys, os, collections, bisect
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'denovo_eval'))
import edgesd as E


def marked(ix, loc, starts, ivs):
    """loc (chrom, start, end); starts/ivs: chrom -> sorted transcript starts / (start, end) with primary support"""
    c, s, e = loc; need = min(300, e - s)
    for k in range(s // E.BIN, e // E.BIN + 1):
        for a, b, c2, a2, b2, _ in ix.get((c, k), ()):
            if min(b, e) - max(a, s) < need: continue
            i = bisect.bisect_right(starts.get(c2, []), b2)
            if any(y > a2 for _, y in ivs.get(c2, [])[max(0, i - 2000):i]): return True   # shortcut: looks back 2000 transcripts, enough unless one locus nests more
    return False


def selftest():
    row = '\t'.join(['c1', '1000', '5000', 'c2', '100', '4100'] + ['x'] * 14 + ['0.95', '0.5'])
    ix = E.sd_index([row]); ivs = {'c2': [(50, 90), (3000, 3500)]}; st = {'c2': [50, 3000]}
    assert marked(ix, ('c1', 2000, 2600), st, ivs) and not marked(ix, ('c1', 6000, 7000), st, ivs)
    assert not marked(ix, ('c1', 2000, 2600), {'c2': [50]}, {'c2': [(50, 90)]}) and marked(ix, ('c1', 4900, 5050), st, ivs) is False
    print('selftest ok')


def main(a):
    rows = [l.rstrip('\n').split('\t') for l in open(a[0])]; h = {x: i for i, x in enumerate(rows[0])}; rows = rows[1:]
    fl = collections.defaultdict(list)
    for r in rows: fl[r[0]].append((int(r[1]), int(r[2])))
    ivs = collections.defaultdict(list)
    for l in open(a[1]):
        f = l.split('\t')
        if len(f) > 4 and f[2] == 'transcript':
            s, e = int(f[3]), int(f[4])
            if not any(x <= s and e <= y for x, y in fl.get(f[0], ())): ivs[f[0]].append((s, e))
    for c in ivs: ivs[c].sort()
    starts = {c: [s for s, _ in v] for c, v in ivs.items()}; ix = E.sd_index(open(a[2]))
    n = collections.Counter((marked(ix, (r[0], int(r[1]), int(r[2])), starts, ivs), r[h['truth']]) for r in rows)
    for m in (True, False):
        d = {k: v for (mm, k), v in n.items() if mm == m}; j = sum(v for k, v in d.items() if k != 'thin')
        print(f'{"marked" if m else "unmarked"}: {sum(d.values())} loci, judged {j}: ' + ' '.join(f'{k} {v} ({v / max(1, j):.3f})' for k, v in sorted(d.items()) if k != 'thin') + f' thin {d.get("thin", 0)}')
    p = n[(True, 'P')] + n[(False, 'P')]; print(f'P loci kept by the mark: {n[(True, "P")]} of {p}')


if __name__ == '__main__':
    selftest() if sys.argv[1:] == ['--selftest'] else main(sys.argv[1:])
