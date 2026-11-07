#!/usr/bin/env python3
"""Row 1570: build the four family definitions on universe U2 (docs/PREREG_structure_vs_sequence_2026-11-07.md).

usage: structdefs.py RUN_PREFIX OUT_DIR   -> OUT_DIR/{SEQ,STRUCT_EXACT,STRUCT_TOL,AND}.copies.tsv + U2.gtf
       structdefs.py --selftest
"""
import collections, os, re, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import structfam as S

MIN_INT, TOL = 2, 3


def components(items, edges):
    """Union-find over item indices -> list of index lists (size >= 2 only)."""
    parent = list(range(len(items)))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for a, b in edges:
        parent[find(a)] = find(b)
    groups = collections.defaultdict(list)
    for i in range(len(items)):
        groups[find(i)].append(i)
    return [g for g in groups.values() if len(g) >= 2]


def tol_edges(loci):
    """Pairs with the same n_int and every internal exon within TOL bp (sorted window on the first exon)."""
    by_n = collections.defaultdict(list)
    for i, x in enumerate(loci):
        by_n[len(x['ints'])].append(i)
    for idx in by_n.values():
        idx.sort(key=lambda i: loci[i]['ints'])
        for p, i in enumerate(idx):
            a = loci[i]['ints']
            for j in idx[p + 1:]:
                b = loci[j]['ints']
                if b[0] > a[0] + TOL:
                    break
                if all(abs(x - y) <= TOL for x, y in zip(a, b)):
                    yield i, j


def write(path, name, groups, loci):
    with open(path, 'w') as fh:
        fh.write('family_id\ttid\tchrom\tstrand\texons\n')
        for k, g in enumerate(groups):
            for i in g:
                x = loci[i]
                fh.write(f"{name}{k}\t{x['tid']}\t{x['chrom']}\t{x['strand']}\t{x['exons']}\n")
    return sum(len(g) for g in groups), len(groups)


def main(prefix, out):
    os.makedirs(out, exist_ok=True)
    members, null, _ = S.load(prefix)
    loci = [x for x in members + null if x['ints'] and len(x['ints']) >= MIN_INT]
    log = [f"U2 loci={len(loci)} (members {sum(x['fam'] is not None for x in loci)})"]

    seq = collections.defaultdict(list)
    for i, x in enumerate(loci):
        if x['fam'] is not None:
            seq[x['fam']].append(i)
    arms = {'SEQ': [g for g in seq.values() if len(g) >= 2]}
    exact = collections.defaultdict(list)
    for i, x in enumerate(loci):
        exact[x['ints']].append(i)
    arms['STRUCT_EXACT'] = [g for g in exact.values() if len(g) >= 2]
    arms['STRUCT_TOL'] = components(loci, tol_edges(loci))
    and_groups = []
    for g in arms['SEQ']:
        sub = [loci[i] for i in g]
        e = [(a, b) for a in range(len(sub)) for b in range(a + 1, len(sub)) if S.alike(sub[a]['ints'], sub[b]['ints'])]
        and_groups += [[g[i] for i in c] for c in components(sub, e)]
    arms['AND'] = and_groups
    for name, groups in arms.items():
        n, k = write(f'{out}/{name}.copies.tsv', name, groups, loci)
        sizes = sorted((len(g) for g in groups), reverse=True)
        log.append(f"{name}: families={k} loci={n} largest={sizes[:3]}")

    keep = {x['gene_id'] for x in loci}
    with open(prefix + '.families.gtf') as fh, open(f'{out}/U2.gtf', 'w') as oh:
        for line in fh:
            m = re.search(r'gene_id "([^"]+)"', line)
            if m and m.group(1) in keep:
                oh.write(line)
    print('\n'.join(log))
    open(f'{out}/arms.log', 'w').write('\n'.join(log) + '\n')


def selftest():
    loci = [dict(ints=(100, 200)), dict(ints=(102, 199)), dict(ints=(104, 201)), dict(ints=(100, 200, 50))]
    assert sorted(tol_edges(loci)) == [(0, 1), (1, 2)]               # 0-2 differ by 4 at exon 1, linked through 1
    assert components(loci, tol_edges(loci)) == [[0, 1, 2]]
    print('selftest ok')


if __name__ == '__main__':
    selftest() if sys.argv[1:] == ['--selftest'] else main(*sys.argv[1:3])
