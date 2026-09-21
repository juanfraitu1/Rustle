#!/usr/bin/env python3
"""Arm B of `docs/PREREG_jaccard_sweep_2026-09-20.md`: the simplest possible family rule.

Connected components of the pair graph, keeping every pair whose alignment Jaccard

    J = nmatch / (len_a + len_b - nmatch)

is >= t. No MCL, no shared-exon conjunct, no core step, no read corroboration — the whole point is that
this is the rule an advisor would write on a whiteboard. Emits a `.clusters.tsv` in `mcl_families`'
own column order so the SAME scorer (`bench/heldout_family_score.py`) can score both arms.

Usage: jaccard_sweep_families.py --paf P --chrom chrN --out PREFIX --t 0.05
"""
import argparse
import collections


def parse(paf, chrom):
    """(chrom,start,end) keyed spans -> best Jaccard per unordered pair."""
    best = {}
    span = {}
    for line in open(paf):
        f = line.rstrip('\n').split('\t')
        if len(f) < 11:
            continue
        a, b = f[0], f[5]
        if a == b:
            continue
        la, lb, nm = int(f[1]), int(f[6]), int(f[9])
        den = la + lb - nm
        if den <= 0:
            continue
        j = nm / den
        k = tuple(sorted((a, b)))
        if j > best.get(k, 0.0):
            best[k] = j
        span[a] = la; span[b] = lb
    return best, span


def main():
    ap = argparse.ArgumentParser()
    for x in ('--paf', '--chrom', '--out'):
        ap.add_argument(x, required=True)
    ap.add_argument('--t', type=float, required=True)
    a = ap.parse_args()

    best, span = parse(a.paf, a.chrom)
    parent = {}

    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x

    def union(x, y):
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[rx] = ry

    kept = 0
    for (x, y), j in best.items():
        if j >= a.t:
            union(x, y); kept += 1

    comp = collections.defaultdict(list)
    for n in list(parent):
        comp[find(n)].append(n)
    clusters = {k: v for k, v in comp.items() if len(v) >= 2}

    out = a.out + '.clusters.tsv'
    with open(out, 'w') as fh:
        fh.write('cluster_id\tsize\tdensity\tfrac_in\tcorroborated\tchrom\tstart\tend\n')
        for i, (_, mem) in enumerate(sorted(clusters.items(), key=lambda kv: -len(kv[1]))):
            for name in sorted(mem):
                # faidx names are chrom:start-end, 1-based inclusive; mcl_families writes start0
                try:
                    ch, rng = name.rsplit(':', 1); s, e = rng.split('-')
                    s0 = int(s) - 1; e0 = int(e)
                except ValueError:
                    continue
                fh.write(f'JAC{i}\t{len(mem)}\tNA\tNA\tNA\t{ch}\t{s0}\t{e0}\n')
    sizes = sorted((len(v) for v in clusters.values()), reverse=True)
    print(f'{a.chrom} t={a.t}: kept {kept} pairs | {len(clusters)} components >=2 | '
          f'{sum(sizes)} members | largest {sizes[0] if sizes else 0} -> {out}')


if __name__ == '__main__':
    main()
