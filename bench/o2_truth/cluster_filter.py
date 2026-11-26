#!/usr/bin/env python3
"""The 10-point lead between clusters of copies on one reference (register 1974).

    cluster_filter.py GENES.tsv GENES_SELF.paf FAM.copies.tsv CHM13.bam TRUTH_READS.tsv

GENES.tsv from eichler_haplotypes.py genes; GENES_SELF.paf = the genes' sequences against themselves (minimap2 -cx asm20 -X);
FAM.copies.tsv = the pipeline's families; TRUTH_READS.tsv = EICHLER_READS_OUT of eichler_haplotypes.py score (read,
verdict, paralog). Clusterings: copy (one gene), 99%, 95% (single linkage on identity over >= half of the shorter gene),
pipeline (families stage), family (annotation). Per read and score (ms, AS): best score per cluster; the top cluster is
retained when it leads every other cluster and every alignment on no family gene by >= 10."""
import sys, collections

def uf():
    p = {}
    def find(x):
        while p.setdefault(x, x) != x: p[x] = p[p[x]]; x = p[x]
        return x
    return p, find

def identity_clusters(genes, paf, thr):
    p, find = uf(); L = {}
    for l in open(paf):
        f = l.split('\t'); a, b = f[0], f[5]; L[a] = int(f[1]); L[b] = int(f[6])
        if a != b and int(f[10]) >= 0.5 * min(int(f[1]), int(f[6])) and int(f[9]) / max(1, int(f[10])) >= thr: p[find(a)] = find(b)
    return {g: find(g) for g in genes}

def decide(scores, lead=10):
    """scores {cluster or None: best score} -> retained cluster or None (ambiguous)"""
    best = max((v, k is not None, k) for k, v in scores.items() if k is not None) if any(k is not None for k in scores) else None
    if best is None: return None
    top, _, c = best
    return c if all(top - v >= lead for k, v in scores.items() if k != c) else None

def main(gtsv, paf, famp, bam, truthp):
    import pysam
    G = {}; byc = collections.defaultdict(list)
    for l in list(open(gtsv))[1:]:
        g, name, fam, c, s, e, st, ex = l.rstrip('\n').split('\t'); exs = [tuple(map(int, x.split('-'))) for x in ex.split(',')]; G[g] = (name, fam); byc[c].append((int(s), int(e), exs, g))
    F = collections.defaultdict(list)
    for l in list(open(famp))[1:]:
        f = l.rstrip('\n').split('\t'); F[f[3]].append((int(f[4]), int(f[5]), [tuple(map(int, x.split('-'))) for x in f[9].split(',')], f[0]))
    pipe = {}
    for c, v in byc.items():
        for s, e, exs, g in v:
            b = (49, None)
            for fs, fe, fex, fid in F.get(c, ()):
                if fs < e and fe > s:
                    o = sum(max(0, min(y, q) - max(x, p)) for x, y in exs for p, q in fex)
                    if o > b[0]: b = (o, fid)
            pipe[g] = b[1] or ('single:' + g)
    CL = {'copy': {g: g for g in G}, '99%': identity_clusters(G, paf, 0.99), '95%': identity_clusters(G, paf, 0.95), 'pipeline': pipe, 'family': {g: G[g][1] for g in G}}
    size = {k: collections.Counter(v.values()) for k, v in CL.items()}
    print('clusterings: ' + '; '.join(f'{k}: {len(size[k])} clusters over {len(G)} genes, largest {max(size[k].values())}' for k in CL))
    truth = {}
    for l in open(truthp):
        q, v, par = l.rstrip('\n').split('\t')
        if v == 'retained' and par in G: truth[q] = par
    R = collections.defaultdict(list)
    for r in pysam.AlignmentFile(bam, check_sq=False):
        if r.is_unmapped or r.is_supplementary: continue
        bl = r.get_blocks(); lab = (49, None)
        for s, e, exs, g in byc.get(r.reference_name, ()):
            if s < bl[-1][1] and e > bl[0][0]:
                o = sum(max(0, min(y, q) - max(x, p)) for x, y in bl for p, q in exs)
                if o > lab[0]: lab = (o, g)
        R[r.query_name].append((lab[1], r.get_tag('ms'), r.get_tag('AS'), not r.is_secondary))
    def table(title, keep):
        print(f'## {title}')
        print('clusters\tscore\treads\tretained\tshare\tretained with a truth\ttruth in the retained cluster\tshare\tcopies in the retained cluster, median / mean\tambiguous')
        for k, cl in CL.items():
            for si, sn in ((1, 'ms'), (2, 'AS')):
                n = ret = wt = ok = 0; sz = []
                for q, al in R.items():
                    if not keep(q, al): continue
                    n += 1; sc = {}
                    for a in al:
                        c = cl[a[0]] if a[0] else None; sc[c] = max(sc.get(c, -10 ** 9), a[si])
                    d = decide(sc)
                    if d is None: continue
                    ret += 1; sz.append(size[k][d])
                    if q in truth: wt += 1; ok += cl[truth[q]] == d
                sz.sort(); print(f"{k}\t{sn}\t{n}\t{ret}\t{ret / max(1, n):.3f}\t{wt}\t{ok}\t{ok / max(1, wt):.3f}\t{sz[len(sz) // 2] if sz else 0} / {sum(sz) / max(1, len(sz)):.1f}\t{n - ret}")
    def near(al):
        top = max(a[1] for a in al); return sum(1 for a in al if a[1] >= top - 10) > 1
    fam_read = lambda al: any(a[0] for a in al if a[3])
    table('reads that are near ties per copy on ms (primary on a family gene)', lambda q, al: fam_read(al) and near(al))
    table('all reads with the primary on a family gene', lambda q, al: fam_read(al))
    for fam in ('TBC1D3', 'NPIP'): table(f'{fam}: reads with the primary on a member', lambda q, al, f=fam: any(a[0] and a[3] and G[a[0]][0].startswith(f) for a in al))

if __name__ == '__main__':
    if sys.argv[1:] == ['--selftest']:
        assert decide({'A': 100, 'B': 90}) == 'A' and decide({'A': 100, 'B': 91}) is None and decide({'A': 100}) == 'A' and decide({'A': 100, None: 95}) is None and decide({None: 5}) is None
        import tempfile, os
        d = tempfile.mkdtemp(); p = os.path.join(d, 'x.paf'); open(p, 'w').write('a\t1000\t0\t1000\t+\tb\t1000\t0\t1000\t995\t1000\t60\na\t1000\t0\t300\t+\tc\t1000\t0\t300\t300\t300\t60\n')
        c = identity_clusters(['a', 'b', 'c'], p, 0.99); assert c['a'] == c['b'] != c['c']; print('cluster_filter selftest ok')
    else: main(*sys.argv[1:])
