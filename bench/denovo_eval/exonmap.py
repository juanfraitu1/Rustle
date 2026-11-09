#!/usr/bin/env python3
"""Row 1595 (docs/PREREG_family_exon_map_2026-11-09.md): family exon map -- exon columns shared across a family's members.

  exonmap.py [--threads N] FAM_PREFIX OUT_PREFIX     FAM_PREFIX = PREFIX.fam (reads .clusters.tsv, .loci.gff3, .loci.gated.paf)
  exonmap.py split S1|S50 FAM_PREFIX OUT_PREFIX.columns.tsv OUT.copies.tsv   row 1597: split families into components of
                                                  members sharing >= 1 column (S1) or >= half the smaller member's columns (S50)
  exonmap.py --selftest
Member = a clusters.tsv node key CONTIG:START-END (= the PAF name; its loci.gff3 exons = copies.tsv exons).
Exon match = exonshare.py's 50% reciprocal projection through the pair's max-SHARE record (ties: earlier in the PAF).
Writes OUT.columns.tsv (one row per column), OUT.map.tsv (one row per family), OUT.summary.tsv; verdict on stdout.
"""
import sys, os, csv, collections, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from exonshare import blocks_from_cigar, seg_map, project, inter


def matches(ea, eb, cg, qs, qe, ts, strand):
    """(share, [(i, j)], any_a, any_b): exon i of a matches exon j of b; share as record_scores' SHARE;
    any_a / any_b = exon indices with any aligned base onto any exon base of the other locus (container-style core)."""
    bl = blocks_from_cigar(cg, qs, qe, ts, strand)
    sq, st = seg_map(bl, strand, True), seg_map(bl, strand, False)
    la, lb = sum(b - a for a, b in ea), sum(b - a for a, b in eb)
    a_src = (len(ea), la) <= (len(eb), lb)
    src, dst, segs = (ea, eb, sq) if a_src else (eb, ea, st)
    pairs = []
    for i, (a, b) in enumerate(src):
        pr = project(segs, a, b)
        for j, (c, d) in enumerate(dst):
            ov = inter(pr, c, d)
            if ov >= 0.5 * (b - a) and ov >= 0.5 * (d - c):
                pairs.append((i, j) if a_src else (j, i))
    nsrc = len(src)
    share = len({p[0] if a_src else p[1] for p in pairs}) / nsrc
    any_a = {i for i, (a, b) in enumerate(ea) if any(inter(project(sq, a, b), c, d) > 0 for c, d in eb)}
    any_b = {j for j, (c, d) in enumerate(eb) if any(inter(project(st, c, d), a, b) > 0 for a, b in ea)}
    return share, pairs, any_a, any_b


_REL = {}


def _chunk(job):
    out = []
    for idx, line in job:
        f = line.rstrip('\n').split('\t')
        a, b = f[0], f[5]
        cg = next((x[5:] for x in f[12:] if x.startswith('cg:Z:')), None)
        if not cg:
            continue
        sh, pairs, aa, ab = matches(_REL[a], _REL[b], cg, int(f[2]), int(f[3]), int(f[7]), f[4])
        out.append((idx, a, b, sh, pairs, aa, ab))
    return out


class UF:
    def __init__(self):
        self.p = {}

    def find(self, x):
        self.p.setdefault(x, x)
        while self.p[x] != x:
            self.p[x] = self.p[self.p[x]]; x = self.p[x]
        return x

    def union(self, x, y):
        self.p[self.find(x)] = self.find(y)


def build(fam_members, rel, best):
    """fam_members: fam -> [keys]; rel: key -> exons; best: (a, b) sorted -> (pairs as ((a,i),(b,j)))
    returns columns per family [(fam, [(key, i)])], triangles (n, closed)."""
    cols, tri_n, tri_c = [], 0, 0
    for fam, mem in fam_members.items():
        uf = UF()
        for k in mem:
            for i in range(len(rel[k])):
                uf.find((k, i))
        adj = collections.defaultdict(set)        # exon node -> matched exon nodes
        M = {}
        for x, y in itertools.combinations(sorted(mem), 2):
            if (x, y) in best:
                M[(x, y)] = set(best[(x, y)])
                for u, v in best[(x, y)]:
                    uf.union(u, v); adj[u].add(v); adj[v].add(u)
        g = collections.defaultdict(list)
        for k in mem:
            for i in range(len(rel[k])):
                g[uf.find((k, i))].append((k, i))
        cols += [(fam, sorted(v)) for v in g.values()]
        for bnode, nb in adj.items():           # open paths a - b - c, each (a, b, c) once (a < c), A != C
            for u, w in itertools.combinations(sorted(nb), 2):
                if u[0] == w[0] or u[0] == bnode[0] or w[0] == bnode[0]:
                    continue
                key = (min(u[0], w[0]), max(u[0], w[0]))
                if key not in M:
                    continue
                tri_n += 1
                pr = (u, w) if u[0] < w[0] else (w, u)
                tri_c += pr in M[key]
    return cols, tri_n, tri_c


def classify(col, nmem):
    loci = {k for k, _ in col}
    return 'private' if len(loci) == 1 else 'core' if len(loci) == nmem else 'accessory'


def main(prefix, out, threads):
    import protedge as P, multiprocessing as mp
    fam_of, fam_members = {}, collections.defaultdict(list)
    for r in csv.DictReader(open(prefix + '.clusters.tsv'), delimiter='\t'):
        k = f"{r['chrom']}:{r['start']}-{r['end']}"
        fam_of[k] = r['cluster_id']; fam_members[r['cluster_id']].append(k)
    loci = P.read_loci(prefix + '.loci.gff3')
    off = lambda n: int(n.rsplit(':', 1)[1].split('-')[0]) - 1
    rel = {k: sorted((a - off(k), b - off(k)) for a, b in loci[k][2]) for k in fam_of}
    _REL.update(rel)
    jobs, cur = [], []
    for idx, line in enumerate(open(prefix + '.loci.gated.paf')):
        a, _, rest = line.partition('\t')
        b = rest.split('\t', 5)[4]
        if a == b or a not in fam_of or fam_of.get(b) != fam_of[a]:
            continue
        cur.append((idx, line))
        if len(cur) == 5000:
            jobs.append(cur); cur = []
    jobs.append(cur)
    nrec = sum(len(j) for j in jobs)
    best, core_any = {}, set()
    with mp.Pool(threads) as pool:
        for part in pool.imap_unordered(_chunk, jobs):
            for idx, a, b, sh, pairs, aa, ab in part:
                core_any |= {(a, i) for i in aa} | {(b, j) for j in ab}
                key = (min(a, b), max(a, b))
                pp = [((a, i), (b, j)) if a < b else ((b, j), (a, i)) for i, j in pairs]
                if key not in best or (sh, -idx) > (best[key][0], -best[key][1]):
                    best[key] = (sh, idx, pp)
    best = {k: v[2] for k, v in best.items()}
    cols, tri_n, tri_c = build(fam_members, rel, best)
    nmem = {f: len(m) for f, m in fam_members.items()}
    cls_n, ex_n, coll, nonpriv = collections.Counter(), collections.Counter(), 0, 0
    agree = collections.Counter()
    per_fam = collections.defaultdict(lambda: collections.Counter())
    with open(out + '.columns.tsv', 'w') as w:
        w.write('family\tmembers\tcolumn\tclass\tn_loci\tn_exons\tcollision\tlen_ratio\texons\n')
        for ci, (fam, col) in enumerate(cols):
            c = classify(col, nmem[fam]); cls_n[c] += 1; ex_n[c] += len(col); per_fam[fam][c] += 1
            loci_c = collections.Counter(k for k, _ in col)
            col_coll = c != 'private' and max(loci_c.values()) > 1
            nonpriv += c != 'private'; coll += col_coll
            for node in col:
                agree[(node in core_any, c != 'private')] += 1
            lens = [rel[k][i][1] - rel[k][i][0] for k, i in col]
            ex = ','.join(f'{k}#{i}:{off(k) + rel[k][i][0]}-{off(k) + rel[k][i][1]}' for k, i in col)
            w.write(f'{fam}\t{nmem[fam]}\tC{ci}\t{c}\t{len(loci_c)}\t{len(col)}\t{int(col_coll)}\t{max(lens) / max(1, min(lens)):.2f}\t{ex}\n')
    with open(out + '.map.tsv', 'w') as w:
        w.write('family\tmembers\tcore\taccessory\tprivate\n')
        for fam in sorted(fam_members, key=lambda f: (-nmem[f], f)):
            p = per_fam[fam]; w.write(f'{fam}\t{nmem[fam]}\t{p["core"]}\t{p["accessory"]}\t{p["private"]}\n')
    coll_rate = coll / nonpriv if nonpriv else float('nan'); clos = tri_c / tri_n if tri_n else float('nan')
    lines = [f'families\t{len(fam_members)}', f'members\t{len(fam_of)}', f'same_family_records\t{nrec}',
             f'pairs_with_record\t{len(best)}']
    lines += [f'columns_{c}\t{cls_n[c]}\nexons_{c}\t{ex_n[c]}' for c in ('core', 'accessory', 'private')]
    lines += [f'families_without_core\t{sum(1 for f in fam_members if per_fam[f]["core"] == 0)}',
              f'collision\t{coll}/{nonpriv}\t{coll_rate:.4f}', f'triangles\t{tri_c}/{tri_n}\t{clos:.4f}',
              f'containercore_in_nonprivate\t{agree[(True, True)]}/{agree[(True, True)] + agree[(True, False)]}',
              f'containeraccessory_in_private\t{agree[(False, False)]}/{agree[(False, False)] + agree[(False, True)]}']
    for lab, lo, hi in (('2', 2, 2), ('3-5', 3, 5), ('>=6', 6, 10 ** 9)):
        fs = [f for f in fam_members if lo <= nmem[f] <= hi]
        lines.append(f'size_{lab}\tfamilies {len(fs)}\t' + '\t'.join(f'{c} {sum(per_fam[f][c] for f in fs)}'
                                                                     for c in ('core', 'accessory', 'private')))
    if tri_n < 1000:
        verdict = f'STOP (triangles {tri_n} < 1000)'
    else:
        verdict = ('PASS' if coll_rate <= 0.10 and clos >= 0.80 else 'FAIL') + \
            f' (collision {coll_rate:.4f} <= 0.10, closure {clos:.4f} >= 0.80)'
    lines.append(f'verdict\t{verdict}')
    open(out + '.summary.tsv', 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))


def components(mem, cols_of, mode):
    """mem: member keys; cols_of: key -> set of column ids -> list of member groups."""
    uf = UF()
    for k in mem:
        uf.find(k)
    for x, y in itertools.combinations(mem, 2):
        sh = len(cols_of[x] & cols_of[y])
        need = 1 if mode == 'S1' else max(1, 0.5 * min(len(cols_of[x]), len(cols_of[y])))
        if sh >= need:
            uf.union(x, y)
    g = collections.defaultdict(list)
    for k in mem:
        g[uf.find(k)].append(k)
    return sorted(sorted(v) for v in g.values())


def split(mode, prefix, columns, out):
    fam_members = collections.defaultdict(list); span = {}
    for r in csv.DictReader(open(prefix + '.clusters.tsv'), delimiter='\t'):
        k = f"{r['chrom']}:{r['start']}-{r['end']}"
        fam_members[r['cluster_id']].append(k); span[k] = (r['chrom'], int(r['start']) - 1, int(r['end']))
    cols_of = collections.defaultdict(set)
    for r in csv.DictReader(open(columns), delimiter='\t'):
        for e in r['exons'].split(','):
            cols_of[e.split('#')[0]].add(r['column'])
    comp_of, nsplit = {}, 0
    for f, mem in fam_members.items():
        comps = components(mem, cols_of, mode); nsplit += len(comps) > 1
        for i, c in enumerate(comps):
            for k in c:
                comp_of[k] = f'{f}.s{i}'
    rows = list(csv.DictReader(open(prefix + '.copies.tsv'), delimiter='\t'))
    unmapped = 0
    with open(out, 'w') as w:
        w.write('\t'.join(rows[0].keys()) + '\n')
        for r in rows:
            s, e = int(r['start']), int(r['end'])
            ks = [k for k in fam_members[r['family_id']] if span[k][0] == r['chrom'] and span[k][1] <= s and e <= span[k][2]]
            if ks:
                r['family_id'] = comp_of[min(ks, key=lambda k: span[k][2] - span[k][1])]
            else:
                unmapped += 1; r['family_id'] = f"{r['family_id']}.u{r['copy_idx']}"
            w.write('\t'.join(r.values()) + '\n')
    print(f'{mode}: {nsplit} of {len(fam_members)} families split, {len(set(comp_of.values()))} components, '
          f'copies unmapped to a node {unmapped} -> {out}')


def selftest():
    # identical 3-exon loci, '+': every exon matches its counterpart
    ea = [(0, 100), (200, 300), (400, 500)]
    sh, pairs, aa, ab = matches(ea, ea, '500M', 0, 500, 0, '+')
    assert sh == 1.0 and pairs == [(0, 0), (1, 1), (2, 2)] and aa == ab == {0, 1, 2}
    # b lacks the middle exon (intron there): 2-exon b is the source; a's middle exon unmatched
    eb = [(0, 100), (400, 500)]
    sh, pairs, aa, ab = matches(ea, eb, '500M', 0, 500, 0, '+')
    assert sh == 1.0 and pairs == [(0, 0), (2, 1)] and aa == {0, 2}
    # minus strand: q runs backwards; a's exon (0,100) lands on b's (400,500)
    sh, pairs, _, _ = matches(ea, ea, '500M', 0, 500, 0, '-')
    assert pairs == [(0, 2), (1, 1), (2, 0)], pairs
    # columns, collision, triangles
    rel = {'A': [(0, 1)], 'B': [(0, 1)], 'C': [(0, 1), (2, 3)]}
    best = {('A', 'B'): [(('A', 0), ('B', 0))], ('B', 'C'): [(('B', 0), ('C', 0))], ('A', 'C'): [(('A', 0), ('C', 0))]}
    cols, n, c = build({'F': ['A', 'B', 'C']}, rel, best)
    assert sorted(classify(x, 3) for _, x in cols) == ['core', 'private'] and n == 3 and c == 3, (cols, n, c)
    best[('A', 'C')] = [(('A', 0), ('C', 1))]               # A-C record matches a different exon: open, and a collision
    cols, n, c = build({'F': ['A', 'B', 'C']}, rel, best)
    assert c < n and any(len({k for k, _ in x}) < len(x) for _, x in cols)
    cols = {'A': {1, 2, 3, 4}, 'B': {1, 2, 5, 6}, 'C': {7}, 'D': {1, 8, 9, 10}}
    assert components(['A', 'B', 'C', 'D'], cols, 'S1') == [['A', 'B', 'D'], ['C']]
    assert components(['A', 'B', 'C', 'D'], cols, 'S50') == [['A', 'B'], ['C'], ['D']]
    print('selftest ok')


if __name__ == '__main__':
    a = sys.argv[1:]
    if a[0] == '--selftest':
        selftest()
    elif a[0] == 'split':
        split(*a[1:5])
    else:
        th = int(a[a.index('--threads') + 1]) if '--threads' in a else 8
        a = [x for i, x in enumerate(a) if x != '--threads' and (i == 0 or a[i - 1] != '--threads')]
        main(a[0], a[1], th)
