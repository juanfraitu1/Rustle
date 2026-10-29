#!/usr/bin/env python3
"""EDGE-HUMAN-2 floors (docs/LEDGER_O1.md): per shipped family, the LOWEST best-record identity among its member pairs in the
all-vs-all PAF; global floor = the lowest family floor. Then filter liftmode edges.tsv (4th column = placement identity):
an edge A-B is kept iff identity >= min(floor(family(A)), floor(family(B))); a singleton has no floor; two singletons -> global.

  liftfloor.py floors CLUSTERS_TSV PAF > floors.json
  liftfloor.py filter floors.json EDGES_TSV KEYMAP_TSV > kept_edges.tsv    (KEYMAP: locus_id<TAB>chrom:start-end; stderr: counts)
"""
import collections, itertools, json, sys


def floors(clusters, paf):
    fam = {}
    for line in open(clusters):
        f = line.rstrip('\n').split('\t')
        if f[0] == 'cluster_id':
            continue
        fam[f'{f[5]}:{f[6]}-{f[7]}'] = f[0]
    best = {}
    for line in open(paf):
        f = line.split('\t', 12)
        q, t = f[0], f[5]
        if q == t or fam.get(q) is None or fam.get(q) != fam.get(t):
            continue
        k = (min(q, t), max(q, t))
        best[k] = max(best.get(k, 0.0), int(f[9]) / max(1, int(f[10])))
    members = collections.defaultdict(list)
    for g, c in fam.items():
        members[c].append(g)
    out, nopair = {}, 0
    for c, gs in members.items():
        ids = [best[(a, b)] for a, b in itertools.combinations(sorted(gs), 2) if (a, b) in best]
        if ids:
            out[c] = min(ids)
        else:
            nopair += 1
    g = min(out.values()) if out else 0.0
    vals = sorted(out.values())
    print(f'families {len(members)}, with a PAF pair {len(out)} (none: {nopair}); floor median {vals[len(vals)//2]:.3f} global {g:.3f}', file=sys.stderr)
    return dict(family=out, member=fam, global_floor=g)


def filt(fl, edges, keymap):
    F = json.load(open(fl))
    key = dict(l.rstrip('\n').split('\t')[:2] for l in open(keymap))
    kept = dropped = 0
    for line in open(edges):
        a, b, w, ident = line.rstrip('\n').split('\t')[:4]
        ka, kb = key.get(a, a), key.get(b, b)
        fa, fb = F['member'].get(ka), F['member'].get(kb)
        fs = [F['family'][f] for f in (fa, fb) if f is not None and f in F['family']]
        floor = min(fs) if fs else F['global_floor']
        if float(ident) >= floor:
            kept += 1
            print(f'{ka}\t{kb}\t{w}\t{ident}')
        else:
            dropped += 1
    print(f'edges kept {kept}, dropped {dropped}', file=sys.stderr)


if __name__ == '__main__':
    if sys.argv[1] == 'floors':
        json.dump(floors(sys.argv[2], sys.argv[3]), sys.stdout)
    else:
        filt(*sys.argv[2:5])
