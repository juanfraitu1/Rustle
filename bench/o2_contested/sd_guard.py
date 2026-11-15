#!/usr/bin/env python3
"""Row 1685 (measurement only): SD identity between sibling copies as a guard on O2 assignments
(docs/PREREG_o2_sd_guard_2026-11-13.md).   sd_guard.py CATALOG.fam.copies.tsv SEDEF.bed DEV_DUMP.tsv TEST_DUMP.tsv"""
import sys, os, csv, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'denovo_eval'))
import edgesd as E
ARMS = (0.98, 0.99, 0.995)


def main(cat, bed, dev, test):
    fam = collections.defaultdict(dict)
    for r in csv.DictReader(open(cat), delimiter='\t'): fam[r['family_id']][r['copy_idx']] = (r['chrom'], int(r['start']), int(r['end']))
    ix = E.sd_index(open(bed)); memo = {}
    def attr(f, i):
        if (f, i) not in memo:
            memo[(f, i)] = max([E.sd_id(ix, fam[f][i], loc) or 0 for j, loc in fam[f].items() if j != i] or [0])
        return memo[(f, i)]
    chosen = None
    for name, path in (('development', dev), ('test', test)):
        rows = [l.rstrip('\n').split('\t') for l in open(path)]; rows = [(r[1], r[5], r[3]) for r in rows if r[2] == 'assigned' and r[3] in ('right', 'wrong')]
        R, Wr = sum(v == 'right' for *_, v in rows), sum(v == 'wrong' for *_, v in rows)
        print(f'{name}: right {R} wrong {Wr}; no SD pair with a sibling: right {sum(v == "right" and attr(f, i) == 0 for f, i, v in rows)} wrong {sum(v == "wrong" and attr(f, i) == 0 for f, i, v in rows)}')
        for t in ARMS:
            r_, w_ = sum(v == 'right' and attr(f, i) >= t for f, i, v in rows), sum(v == 'wrong' and attr(f, i) >= t for f, i, v in rows)
            ok = w_ >= 0.5 * Wr and r_ <= 0.2 * R
            fams = collections.Counter(f for f, i, v in rows if v == 'wrong' and attr(f, i) >= t)
            print(f'  T {t}: removes wrong {w_} ({w_ / max(1, Wr):.2f}) right {r_} ({r_ / max(1, R):.2f}) -> {"meets" if ok else "no"}; wrong rows by family {dict(fams)}')
            if name == 'development' and ok and chosen is None: chosen = t
        if name == 'development':
            print('chosen arm:', chosen)
            if chosen is None: print('STOP'); return


if __name__ == '__main__':
    main(*sys.argv[1:5])
