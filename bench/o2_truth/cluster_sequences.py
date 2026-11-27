#!/usr/bin/env python3
"""Do recurring differences among tied reads mark where the cluster verdict is wrong? (register 1978)

    cluster_sequences.py SEQUENCES.tsv TIES.tsv.gz CLUSTER_LOCI.tsv ASSEMBLY.gtf TRUTH.tsv GENES.tsv

SEQUENCES.tsv = bench/o2_families/recurrence_test.py's output (alleles = `reference` or the differing sites). Decided
reads with a truth (row 1972; any assembled locus on the truth gene) by whether their sequence is non-reference."""
import sys, os, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'tools'))
import cluster_filter_run as R, cluster_verdict as V

def table(rows):
    """rows [(kind, right)] -> {kind: (n, wrong)}"""
    c = collections.defaultdict(lambda: [0, 0])
    for k, ok in rows: c[k][0] += 1; c[k][1] += not ok
    return {k: tuple(v) for k, v in c.items()}

def main(seqp, ties, cluster_loci, gtf, truthp, genesp):
    by, st, mx = R.load_loci(gtf); cl = {}
    for l in list(open(cluster_loci))[1:]: f = l.split('\t'); cl[f[0]] = f[1]
    gl = {}
    for l in list(open(genesp))[1:]:
        g, name, fam, c, s, e, stn, ex = l.rstrip('\n').split('\t'); gl[g] = (R.loci_under(by, st, mx, c, [tuple(map(int, x.split('-'))) for x in ex.split(',')]), name)
    truth = {}
    for l in open(truthp):
        q, v, par = l.rstrip('\n').split('\t')
        if v == 'retained' and par in gl and gl[par][0]: truth[q] = par
    kind = {}; nsites = {}
    for l in list(open(seqp))[1:]:
        f = l.rstrip('\n').split('\t'); k = 'reference sequence' if f[5] == 'reference' else 'non-reference sequence'
        for q in f[8].split(','):
            if q: kind[q] = k; nsites[q] = 0 if f[5] == 'reference' else len(f[5].split(','))
    rows = []; fam = collections.defaultdict(lambda: [0, 0, 0, 0])
    for q, r in V.read_ties(ties).items():
        v = V.verdict(r, cl, 10)
        if v is None or q not in truth: continue
        ok = any(cl.get(x) == v for x in gl[truth[q]][0]); k = kind.get(q, 'not tested'); rows.append((k, ok))
        if k != 'not tested': x = fam[gl[truth[q]][1].rstrip('0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ') or gl[truth[q]][1]]; x[0 if k.startswith('non') else 2] += 1; x[1 if k.startswith('non') else 3] += not ok
    T = table(rows); n = len(rows); w = sum(1 for _, ok in rows if not ok)
    print(f'decided reads with a truth: {n}, wrong {w} ({w / max(1, n):.3f})')
    print('kind\treads\twrong\tshare wrong')
    for k in ('non-reference sequence', 'reference sequence', 'not tested'):
        a, b = T.get(k, (0, 0)); print(f'{k}\t{a}\t{b}\t{b / max(1, a):.3f}')
    nr, nw = T.get('non-reference sequence', (0, 0)); rr, rw = T.get('reference sequence', (0, 0)); cov = nr + rr
    print(f'covered by the test: {cov} of {n} ({cov / max(1, n):.3f}); wrong among the covered {nw + rw}')
    print(f'flagging non-reference sequences: catches {nw} of {nw + rw} wrong ({nw / max(1, nw + rw):.3f}); flags {nr - nw} of {cov - nw - rw} right ({(nr - nw) / max(1, cov - nw - rw):.3f})')
    print('gene-name stems with the most flagged reads: ' + '; '.join(f'{k}: flagged {v[0]} (wrong {v[1]}), unflagged {v[2]} (wrong {v[3]})' for k, v in sorted(fam.items(), key=lambda z: -z[1][0])[:8]))

if __name__ == '__main__':
    if sys.argv[1:] == ['--selftest']: assert table([('a', True), ('a', False), ('b', True)]) == {'a': (2, 1), 'b': (1, 0)}; print('cluster_sequences selftest ok')
    else: main(*sys.argv[1:])
