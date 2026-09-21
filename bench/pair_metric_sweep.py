#!/usr/bin/env python3
"""Pair-metric bakeoff for `docs/PREREG_asymmetry_metric_2026-09-20.md`.

r359 and r913 are the two horns of the same problem:

    jaccard      m / (la + lb - m)     PENALISES length asymmetry (EEF1A1 retrocopies J ~ 0.07 at
                                       core identity 0.92)
    containment  m / min(la, lb)       SATURATES on it -- a short gene inside a long one scores ~1.0,
                                       so long genes become hubs (864 chr8 pairs >= 1.0)

    ochiai       m / sqrt(la * lb)     the GEOMETRIC MEAN of the two containments: equal to containment
                                       when la == lb, but it cannot reach 1.0 unless the alignment
                                       covers BOTH genes. Parameter-free.
    dice         2m / (la + lb)        symmetric but milder than Jaccard
    guarded      m / min(la, lb), but only if min/max >= --guard

Families = connected components of pairs scoring >= t. Identical machinery to §6t3 arm B, so the metric
is the only thing that varies. Writes `.clusters.tsv` in `mcl_families`' column order.

Usage: pair_metric_sweep.py --paf P --chrom chrN --out PREFIX --metric ochiai --t 0.30
"""
import argparse
import collections
import math


def best_records(paf):
    """(a, b) -> (best m, la, lb) over the pair's PAF records."""
    best = {}
    for line in open(paf):
        f = line.rstrip('\n').split('\t')
        if len(f) < 11:
            continue
        a, b = f[0], f[5]
        if a == b:
            continue
        la, lb, m = int(f[1]), int(f[6]), int(f[9])
        k = tuple(sorted((a, b)))
        if k[0] != a:
            la, lb = lb, la
        if k not in best or m > best[k][0]:
            best[k] = (m, la, lb)
    return best


def score(metric, m, la, lb, guard):
    if metric == 'jaccard':
        d = la + lb - m
        return m / d if d > 0 else 0.0
    if metric == 'dice':
        d = la + lb
        return 2 * m / d if d > 0 else 0.0
    if metric == 'ochiai':
        d = math.sqrt(float(la) * float(lb))
        return m / d if d > 0 else 0.0
    if metric == 'containment':
        d = min(la, lb)
        return m / d if d > 0 else 0.0
    if metric == 'union':
        # POST-HOC (§6t5). The measurement says one threshold cannot serve two populations: symmetric
        # metrics get ~96% of similar-length true pairs and 0% of asymmetric ones, while containment
        # gets both but builds hubs. So accept a pair if EITHER criterion fires, with containment
        # guarded by a length-ratio floor so it cannot hub.
        d = la + lb - m
        jac = m / d if d > 0 else 0.0
        lo, hi = min(la, lb), max(la, lb)
        con = (m / lo) if (lo > 0 and hi > 0 and lo / hi >= guard) else 0.0
        return max(jac, con)
    if metric == 'guarded':
        lo, hi = min(la, lb), max(la, lb)
        if hi == 0 or lo / hi < guard:
            return 0.0
        return m / lo if lo > 0 else 0.0
    raise SystemExit(f'unknown metric {metric}')


def main():
    ap = argparse.ArgumentParser()
    for x in ('--paf', '--chrom', '--out', '--metric'):
        ap.add_argument(x, required=True)
    ap.add_argument('--t', type=float, required=True)
    ap.add_argument('--guard', type=float, default=0.25)
    a = ap.parse_args()

    best = best_records(a.paf)
    parent = {}

    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x

    kept = 0
    for (x, y), (m, la, lb) in best.items():
        if score(a.metric, m, la, lb, a.guard) >= a.t:
            rx, ry = find(x), find(y)
            if rx != ry:
                parent[rx] = ry
            kept += 1
    comp = collections.defaultdict(list)
    for n in list(parent):
        comp[find(n)].append(n)
    clusters = {k: v for k, v in comp.items() if len(v) >= 2}

    out = a.out + '.clusters.tsv'
    with open(out, 'w') as fh:
        fh.write('cluster_id\tsize\tdensity\tfrac_in\tcorroborated\tchrom\tstart\tend\n')
        for i, (_, mem) in enumerate(sorted(clusters.items(), key=lambda kv: -len(kv[1]))):
            for name in sorted(mem):
                try:
                    ch, rng = name.rsplit(':', 1); s, e = rng.split('-')
                except ValueError:
                    continue
                fh.write(f'M{i}\t{len(mem)}\tNA\tNA\tNA\t{ch}\t{int(s)-1}\t{int(e)}\n')
    sizes = sorted((len(v) for v in clusters.values()), reverse=True)
    print(f'{a.chrom} {a.metric} t={a.t}: {kept} pairs | {len(clusters)} components | '
          f'{sum(sizes)} members | largest {sizes[0] if sizes else 0}')


if __name__ == '__main__':
    main()
