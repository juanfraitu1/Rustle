#!/usr/bin/env python3
"""Does a variation-graph signal separate real containments from repeat/domain ones?
Per `docs/PREREG_vg_containment_2026-09-20.md`.

§6t7 left the containment problem with no pairwise signal: among pairs the shipped rule rejects at
containment >= 0.90, there are 19 TRUE and 57 FALSE, alignment identity is 0.987 vs 0.989, and the
containment precision curve is flat at ~0.25 from 0.30 to 0.99.

A variation graph's distinctive information is MULTI-WAY: how many other sequences traverse a segment.
That is node multiplicity, and the existing all-vs-all PAF already determines it — so this tests the
VG hypothesis WITHOUT building a VG. If multiplicity does not separate the classes, a VG cannot help.

    mult(pair) = number of DISTINCT other genes whose PAF alignment overlaps the long gene's aligned
                 interval by >= 50% of that interval (the two genes of the pair excluded)

Usage: vg_multiplicity_containment.py --gff G --soto S1C --pafs chr2=P1,chr8=P2,... --graphs DIR
"""
import argparse
import collections
import csv
import re


def gene_names(gff, chrom):
    n = {}
    for line in open(gff):
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[0] != chrom or f[2] not in ('gene', 'pseudogene'):
            continue
        m = re.search(r'Name=([^;]+)', f[8])
        if m:
            n[f'{f[0]}:{f[3]}-{f[4]}'] = m.group(1)
    return n


def soto(s1c):
    ids = collections.defaultdict(set)
    for r in csv.DictReader(open(s1c), delimiter='\t'):
        fid = (r.get('Family ID') or '').strip(); nm = (r.get('Gene Name') or '').strip()
        if fid and fid != 'N/A' and nm:
            ids[nm].add(fid)
    return {nm: next(iter(v)) for nm, v in ids.items() if len(v) == 1}


def load(paf):
    """best record per pair, plus every alignment interval per gene (for multiplicity)."""
    best = {}
    ivs = collections.defaultdict(list)          # gene -> [(start, end, partner)]
    for line in open(paf):
        f = line.rstrip('\n').split('\t')
        if len(f) < 12:
            continue
        a, b = f[0], f[5]
        if a == b:
            continue
        la, lb, m = int(f[1]), int(f[6]), int(f[9])
        qs, qe, ts, te = int(f[2]), int(f[3]), int(f[7]), int(f[8])
        ivs[a].append((qs, qe, b))
        ivs[b].append((ts, te, a))
        k = tuple(sorted((a, b)))
        rec = (m, la, lb, qs, qe, ts, te) if k[0] == a else (m, lb, la, ts, te, qs, qe)
        if k not in best or m > best[k][0]:
            best[k] = rec
    return best, ivs


def main():
    ap = argparse.ArgumentParser()
    for x in ('--gff', '--soto', '--pafs', '--graphs'):
        ap.add_argument(x, required=True)
    a = ap.parse_args()
    fam = soto(a.soto)

    rows = []
    for spec in a.pafs.split(','):
        chrom, paf = spec.split('=')
        names = gene_names(a.gff, chrom)
        shipped = set()
        for line in open(f'{a.graphs}/{chrom}.graph.tsv'):
            f = line.rstrip('\n').split('\t')
            if len(f) == 3 and f[0] != f[1]:
                shipped.add(tuple(sorted((f[0], f[1]))))
        best, ivs = load(paf)
        for k, (m, la, lb, qs, qe, ts, te) in best.items():
            if k in shipped:
                continue
            ga, gb = names.get(k[0]), names.get(k[1])
            if not ga or not gb:
                continue
            fa, fb = fam.get(ga), fam.get(gb)
            if not fa or not fb:
                continue
            lo = min(la, lb)
            if not lo or m / lo < 0.90:
                continue
            # the LONG gene and its aligned interval
            if la >= lb:
                longg, s, e = k[0], qs, qe
            else:
                longg, s, e = k[1], ts, te
            span = max(1, e - s)
            others = set()
            for (xs, xe, partner) in ivs.get(longg, []):
                if partner in (k[0], k[1]):
                    continue
                ov = min(e, xe) - max(s, xs)
                if ov >= 0.5 * span:
                    others.add(partner)
            rows.append((fa == fb, len(others), ga, gb))

    T = [r for r in rows if r[0]]
    F = [r for r in rows if not r[0]]
    print(f"population: {len(rows)} rejected containment>=0.90 pairs — {len(T)} TRUE, {len(F)} FALSE\n")
    import statistics
    print(f"  multiplicity  TRUE median {statistics.median([r[1] for r in T]) if T else 0:.1f}"
          f"   FALSE median {statistics.median([r[1] for r in F]) if F else 0:.1f}\n")
    print(f"  {'mult <= k':>10} {'TRUE kept':>10} {'FALSE kept':>11} {'precision':>10} {'recall':>8}")
    best_p = 0
    for k in (0, 1, 2, 3, 5, 10, 20, 50, 10**9):
        t = sum(1 for r in T if r[1] <= k); f_ = sum(1 for r in F if r[1] <= k)
        p = t / (t + f_) if t + f_ else 0
        if t >= 10:
            best_p = max(best_p, p)
        lbl = 'any' if k == 10**9 else str(k)
        print(f"  {lbl:>10} {t:>10} {f_:>11} {p:>10.3f} {t/len(T) if T else 0:>8.3f}")
    # AUC (higher multiplicity should mean FALSE, so score = -mult)
    import itertools
    pairs = list(itertools.product([r[1] for r in T], [r[1] for r in F]))
    wins = sum(1 for x, y in pairs if x < y) + 0.5 * sum(1 for x, y in pairs if x == y)
    print(f"\n  AUC (mult separating TRUE from FALSE) = {wins/len(pairs) if pairs else 0:.3f}"
          f"   [r384's general-separator AUC was 0.686]")
    print(f"  baseline precision (containment alone) = 0.250")
    print(f"  best precision at >=10 TRUE retained   = {best_p:.3f}")


if __name__ == '__main__':
    main()
