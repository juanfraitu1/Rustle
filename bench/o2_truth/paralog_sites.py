#!/usr/bin/env python3
"""Are the recurring differences of a tied sequence the bases of another copy? (register 1979)

    paralog_sites.py SEQUENCES.tsv BAM TIES.tsv.gz CLUSTER_LOCI.tsv ASSEMBLY.gtf TRUTH.tsv GENES.tsv

Per non-reference sequence of recurrence_test.py: each substitution site (reference position at the locus the reads are
read at) is carried to the read's own coordinate, and looked up in the read's alignments to other loci (--eqx: match or
mismatch). paralog-explained = one other locus matches at every site (majority of the reads that reach it); private =
no locus matches at any; mixed otherwise."""
import sys, os, re, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'tools'))
import cluster_filter_run as R, cluster_verdict as V

def ref_to_read(cig, start, rev, qlen):
    """{reference pos: read pos in the read's own orientation} for aligned bases"""
    out = {}; rp = start; qp = 0
    for op, n in cig:
        if op in (0, 7, 8):
            for i in range(n): out[rp + i] = (qlen - 1 - (qp + i)) if rev else qp + i
            rp += n; qp += n
        elif op in (2, 3): rp += n
        elif op in (1, 4, 5): qp += n
    return out

def state_at(cig, rev, qlen, pos):
    """the alignment's column at the read's own position: '=' match, 'X' mismatch, None not aligned there"""
    want = (qlen - 1 - pos) if rev else pos; qp = 0
    for op, n in cig:
        if op in (0, 7, 8, 1, 4, 5):
            if qp <= want < qp + n: return '=' if op == 7 else 'X' if op == 8 else None
            qp += n
    return None

def classify(site_by_locus, nsites):
    """site_by_locus {locus: {site: [matches, reads]}} -> ('paralog-explained', Y) | ('private', None) | ('mixed', best Y)"""
    best = (0, None)
    for y, d in site_by_locus.items():
        k = sum(1 for s, (m, n) in d.items() if n and m * 2 > n)
        if k > best[0]: best = (k, y)
    return ('paralog-explained', best[1]) if best[0] == nsites else ('private', None) if best[0] == 0 else ('mixed', best[1])

def main(seqp, bam, ties, cluster_loci, gtf, truthp, genesp):
    import pysam
    by, st, mx = R.load_loci(gtf); cl = {}
    for l in list(open(cluster_loci))[1:]: f = l.split('\t'); cl[f[0]] = f[1]
    S = []; want = set()
    for l in list(open(seqp))[1:]:
        f = l.rstrip('\n').split('\t')
        if f[5] == 'reference': continue
        sites = [int(m.group(1)) for m in (re.match(r'(\d+)S', x) for x in f[5].split(',')) if m]
        if not sites: continue
        reads = [q for q in f[8].split(',') if q]; S.append((f[2], sites, reads)); want.update(reads)
    A = collections.defaultdict(list)
    for r in pysam.AlignmentFile(bam):
        if r.query_name in want and not (r.is_unmapped or r.is_supplementary):
            A[r.query_name].append((R.locus_of(by, st, mx, r.reference_name, r.get_blocks()), r.cigartuples, r.reference_start, r.is_reverse, r.infer_read_length()))
    kind = {}; expl = {}; nseq = collections.Counter()
    for read_at, sites, reads in S:
        d = collections.defaultdict(lambda: collections.defaultdict(lambda: [0, 0]))
        for q in reads:
            home = next((a for a in A[q] if a[0] == read_at), None)
            if home is None: continue
            m = ref_to_read(home[1], home[2], home[3], home[4]); c0 = cl.get(read_at)
            for y, cig, start, rev, qlen in A[q]:
                if y is None or y == read_at or (c0 is not None and cl.get(y) == c0): continue
                for s in sites:
                    if s in m:
                        x = state_at(cig, rev, qlen, m[s])
                        if x: d[y][s][1] += 1; d[y][s][0] += x == '='
        k, y = classify(d, len(sites)); nseq[k] += 1
        for q in reads: kind[q] = k; expl[q] = y
    print('non-reference sequences with substitution sites: ' + ', '.join(f'{k} {v}' for k, v in nseq.most_common()))
    gl = {}
    for l in list(open(genesp))[1:]:
        g, name, fam, c, s, e, stn, ex = l.rstrip('\n').split('\t'); gl[g] = R.loci_under(by, st, mx, c, [tuple(map(int, x.split('-'))) for x in ex.split(',')])
    truth = {}
    for l in open(truthp):
        q, v, par = l.rstrip('\n').split('\t')
        if v == 'retained' and gl.get(par): truth[q] = gl[par]
    T = collections.defaultdict(lambda: [0, 0]); y_ok = y_n = 0
    for q, r in V.read_ties(ties).items():
        v = V.verdict(r, cl, 10)
        if v is None or q not in truth or q not in kind: continue
        ok = any(cl.get(x) == v for x in truth[q]); T[kind[q]][0] += 1; T[kind[q]][1] += not ok
        if kind[q] == 'paralog-explained' and not ok:
            y_n += 1; tc = {cl.get(x) for x in truth[q]} - {None}; y_ok += (cl.get(expl[q]) in tc) or expl[q] in truth[q]
    print('sequence\tdecided reads with a truth\twrong\tshare wrong')
    for k in ('paralog-explained', 'mixed', 'private'): print(f'{k}\t{T[k][0]}\t{T[k][1]}\t{T[k][1] / max(1, T[k][0]):.3f}')
    n = sum(v[0] for v in T.values()); w = sum(v[1] for v in T.values()); a, b = T['paralog-explained']
    print(f'flagging paralog-explained: catches {b} of {w} wrong ({b / max(1, w):.3f}); flags {a - b} of {n - w} right ({(a - b) / max(1, n - w):.3f})')
    print(f'flagged wrong reads whose explaining locus is at the truth: {y_ok} of {y_n}')

if __name__ == '__main__':
    if sys.argv[1:] == ['--selftest']:
        cig = [(4, 2), (7, 3), (8, 1), (7, 2)]; m = ref_to_read(cig, 100, False, 8); assert m[100] == 2 and m[103] == 5 and state_at(cig, False, 8, 5) == 'X' and state_at(cig, False, 8, 2) == '=' and state_at(cig, False, 8, 0) is None
        mr = ref_to_read(cig, 100, True, 8); assert mr[100] == 5 and mr[103] == 2 and state_at(cig, True, 8, 2) == 'X'
        assert classify({'Y': {1: [3, 3], 2: [2, 3]}}, 2) == ('paralog-explained', 'Y') and classify({'Y': {1: [3, 3], 2: [1, 3]}}, 2) == ('mixed', 'Y') and classify({'Y': {1: [0, 3]}}, 1) == ('private', None) and classify({}, 2) == ('private', None)
        print('paralog_sites selftest ok')
    else: main(*sys.argv[1:])
