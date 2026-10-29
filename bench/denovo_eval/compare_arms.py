#!/usr/bin/env python3
"""Compare two mcl_families runs (baseline vs --shared-exon-columns) on the same loci: partition differences, pair
agreement, families that change, and the graph counters from the logs."""
import collections, itertools, re, sys

def read_clusters(path):
    fam = {}
    with open(path) as fh:
        header = fh.readline().rstrip('\n').split('\t')
        c = {h: i for i, h in enumerate(header)}
        for line in fh:
            r = line.rstrip('\n').split('\t')
            if len(r) < len(header):
                continue
            fam[(r[c['chrom']], int(r[c['start']]), int(r[c['end']]))] = r[c['cluster_id']]
    return fam

def members(fam):
    m = collections.defaultdict(set)
    for k, f in fam.items():
        m[f].add(k)
    return m

def pairs(m):
    out = set()
    for f, ks in m.items():
        if len(ks) >= 2:
            for a, b in itertools.combinations(sorted(ks), 2):
                out.add((a, b))
    return out

def log_counts(path):
    keys = ('rejected_low_shared_exon', 'rejected_no_exonic', 'pair(s) dropped for sharing too small', 'pair(s) dropped for resting on no exonic',
            'edges', 'nodes', 'gate', 'all-vs-all', 'records')
    out = []
    try:
        with open(path) as fh:
            for line in fh:
                if any(k in line for k in keys) and len(line) < 400:
                    out.append(line.rstrip())
    except OSError:
        pass
    return out

def main(base, cols, name_a='baseline', name_b='columns'):
    fa, fb = read_clusters(base), read_clusters(cols)
    ma, mb = members(fa), members(fb)
    print(f'{name_a}: {len(fa)} loci in {len(ma)} clusters ({sum(1 for v in ma.values() if len(v) >= 2)} with >= 2 members, '
          f'{sum(len(v) for v in ma.values() if len(v) >= 2)} members)')
    print(f'{name_b}: {len(fb)} loci in {len(mb)} clusters ({sum(1 for v in mb.values() if len(v) >= 2)} with >= 2 members, '
          f'{sum(len(v) for v in mb.values() if len(v) >= 2)} members)')
    pa, pb = pairs(ma), pairs(mb)
    print(f'co-clustered pairs: {name_a} {len(pa)}, {name_b} {len(pb)}, shared {len(pa & pb)}, '
          f'only {name_a} {len(pa - pb)}, only {name_b} {len(pb - pa)}, Jaccard {len(pa & pb) / max(1, len(pa | pb)):.4f}')
    # families of A that are not reproduced as a whole in B
    inv_b = {k: f for f, ks in mb.items() for k in ks}
    changed = 0
    rows = []
    for f, ks in sorted(ma.items(), key=lambda x: -len(x[1])):
        if len(ks) < 2:
            continue
        tgt = collections.Counter(inv_b.get(k, 'NONE') for k in ks)
        whole = len(tgt) == 1 and len(mb.get(next(iter(tgt)), ())) == len(ks)
        if not whole:
            changed += 1
            if len(rows) < 25:
                rows.append((f, len(ks), dict(tgt), {t: len(mb.get(t, ())) for t in tgt}))
    print(f'{name_a} families (>= 2 members) not reproduced as a whole in {name_b}: {changed}')
    for f, n, tgt, sizes in rows:
        print(f'  {f} ({n}) -> ' + ', '.join(f'{t}:{c}/{sizes[t]}' for t, c in tgt.items()))

if __name__ == '__main__':
    main(*sys.argv[1:])
