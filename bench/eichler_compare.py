#!/usr/bin/env python3
"""Eichler-style AS-margin assignment, computed alongside ours and compared.

The method the advisor cites: a multi-mapping read is assigned to its best alignment iff no other
alignment scores within T units of it (T = 10 by convention); otherwise the read is discarded as
ambiguous. It is a MARGIN rule on the aligner's own score.

No pipeline change is needed to compute it: `copy_assign` already emits `as_best`, `as_second` and
`as_margin` per read alongside our `status`, so both calls come from one file.

    EICHLER(T):  as_margin >= T          -> assign to the best-AS copy
                 as_margin <  T          -> discard (ambiguous)

    OURS:        status in {assigned, tied, ambiguous}, assign-or-abstain, never 1/k

⚠ The two rules do not have the same SUBJECT, and the comparison is meaningless unless that is stated:
our AS-tied gate deliberately selects the reads where the aligner is indifferent (margin ~ 0), which is
exactly the population Eichler's rule discards by construction. So "agreement" is not the interesting
number; the interesting number is what each rule decides on the population the other keeps.

Usage: eichler_compare.py --assignments A.tsv [--threshold 10] [--out report.tsv]
"""
import argparse
import collections
import csv


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--assignments', required=True)
    ap.add_argument('--threshold', type=float, default=10.0)
    ap.add_argument('--out')
    a = ap.parse_args()

    rows = list(csv.DictReader(open(a.assignments), delimiter='\t'))
    if not rows:
        raise SystemExit('no rows')

    def num(r, k):
        v = (r.get(k) or '').strip()
        try:
            return float(v)
        except ValueError:
            return None

    joint = collections.Counter()
    margins = collections.Counter()
    n = 0
    for r in rows:
        m = num(r, 'as_margin')
        if m is None:
            continue
        n += 1
        eich = 'assign' if m >= a.threshold else 'discard'
        ours = (r.get('status') or '').strip()
        joint[(ours, eich)] += 1
        margins['>=T' if m >= a.threshold else ('0' if m == 0 else '0<m<T')] += 1

    print(f"reads with an AS margin: {n}   (Eichler threshold T = {a.threshold:g})\n")
    print("AS-margin distribution")
    for k in ('0', '0<m<T', '>=T'):
        print(f"  margin {k:6s} {margins[k]:>7}  {100*margins[k]/n if n else 0:>5.1f}%")

    ours_vals = sorted({k[0] for k in joint})
    print(f"\njoint decision table  (rows = OURS, cols = EICHLER T={a.threshold:g})")
    print(f"  {'':12s} {'assign':>9} {'discard':>9} {'total':>8}")
    for o in ours_vals:
        aa, dd = joint[(o, 'assign')], joint[(o, 'discard')]
        print(f"  {o:12s} {aa:>9} {dd:>9} {aa+dd:>8}")
    ta = sum(joint[(o, 'assign')] for o in ours_vals)
    td = sum(joint[(o, 'discard')] for o in ours_vals)
    print(f"  {'TOTAL':12s} {ta:>9} {td:>9} {ta+td:>8}")

    ours_assign = sum(v for k, v in joint.items() if k[0] == 'assigned')
    print(f"\n  Eichler assigns  {ta:>7} / {n} = {100*ta/n if n else 0:.1f}%")
    print(f"  we assign        {ours_assign:>7} / {n} = {100*ours_assign/n if n else 0:.1f}%")
    both = joint[('assigned', 'assign')]
    print(f"  both assign      {both:>7}")
    print(f"  we assign where Eichler discards: {joint[('assigned','discard')]}")
    print(f"  Eichler assigns where we abstain: "
          f"{sum(joint[(o,'assign')] for o in ours_vals if o != 'assigned')}")

    if a.out:
        with open(a.out, 'w') as fh:
            fh.write('ours\teichler\tn\n')
            for (o, e), v in sorted(joint.items()):
                fh.write(f'{o}\t{e}\t{v}\n')
        print(f"\n  wrote {a.out}")


if __name__ == '__main__':
    main()
