#!/usr/bin/env python3
"""PREREG row 1715: read chains that are a terminal sub-run of an emitted chain; read-end features, annotated or not.

  nested_sub.py REF_GTF DEFAULT_GTF FLOOR1_GTF BAM FASTA OUT_TSV [THREADS] [CONTIG,...] [--models GTF (row 1716)]   |   nested_sub.py --selftest"""
import sys, collections, multiprocessing
from alt_events import load, count
from alt_ends import flags
MIN = 5

def features(M, K, k, suf, cnt, lo, hi):
    """M = K[k:] (suf) or K[:k]. -> (inside share, peak share, commonest end, reads continuing past M on its variable side)"""
    ctr = (lo if suf else hi)[M]; n = cnt[M]; b = K[k - 1][1] if suf else K[k][0]
    inside = sum(v for p, v in ctr.items() if (p >= b if suf else p <= b)) / n
    mode = min(ctr, key=lambda p: (-ctr[p], p if suf else -p)); peak = sum(v for p, v in ctr.items() if abs(p - mode) <= 10) / n
    return inside, peak, mode

def cont_counts(cnt, cands):
    c = collections.Counter()
    for R, n in cnt.items():
        for k in range(1, len(R)):
            if (R[k:], True) in cands: c[(R[k:], True)] += n
            if (R[:k], False) in cands: c[(R[:k], False)] += n
    return c

def work(job):
    bam, fa_p, contig, D, ref, fl = job; import pysam
    _, cnt, lo, hi = count((bam, contig)); fa = pysam.FastaFile(fa_p); best = {}
    for K, (strand, attr, ks, ke) in D.items():
        for k in range(1, len(K)):
            for M, suf in ((K[k:], True), (K[:k], False)):
                if cnt.get(M, 0) >= MIN and cnt.get(K, 0) >= best.get((M, suf), (-1,))[0]: best[(M, suf)] = (cnt.get(K, 0), K, k, strand, attr, ks, ke)
    cont = cont_counts(cnt, best); rows = []; models = []
    for (M, suf), (rk, K, k, strand, attr, ks, ke) in best.items():
        ins, peak, mode = features(M, K, k, suf, cnt, lo, hi); plus = strand == '+'; cls = 'NEST5' if suf == plus else 'NEST3'
        sig, arich = flags(fa.fetch, contig, plus, mode) if cls == 'NEST3' else (False, False)
        third = (cont[(M, suf)] < cnt[M]) if cls == 'NEST5' else (sig and not arich)
        rows.append((contig, strand, cls, int(M in D), int(M in ref.get(strand, ())), int(M in fl), cnt[M], rk, attr if attr is not None else -1,
                     cont[(M, suf)], round(ins, 3), round(peak, 3), int(third), len(M), len(K)))
        grp = ('E' if M in D else 'N') + ('A' if M in ref.get(strand, ()) else 'U')
        models.append((f'{cls}_{grp}_{int(peak >= 0.5)}_{cnt[M]}_{len(models)}', strand, M, mode if suf else ks, ke if suf else mode))
    for K in {v[1] for v in best.values()}:   # FULL: annotated containers, 5' end = commonest start of their own reads
        strand, _, ks, ke = D[K]
        if K in ref.get(strand, ()) and cnt.get(K, 0) >= MIN:
            c5 = lo[K] if strand == '+' else hi[K]; m = min(c5, key=lambda p: (-c5[p], p)); pk = sum(v for p, v in c5.items() if abs(p - m) <= 10) / cnt[K]
            models.append((f'FULL_EA_{int(pk >= 0.5)}_{cnt[K]}_{len(models)}', strand, K, m if strand == '+' else ks, ke if strand == '+' else m))
    return rows, [gtf_lines(contig, *m) for m in models]

def gtf_lines(contig, tid, strand, chain, start, end):
    """0-based half-open start/end + introns -> GTF transcript and exon lines; None when an end falls inside the chain."""
    if not (start < chain[0][0] and end > chain[-1][1]): return None
    b = [start] + [i[1] for i in chain]; e = [i[0] for i in chain] + [end]; a = f'gene_id "{contig}_{tid}"; transcript_id "{contig}_{tid}";'
    return ''.join(f'{contig}\tnested\t{f}\t{x + 1}\t{y}\t.\t{strand}\t.\t{a}\n' for f, x, y in [('transcript', start, end)] + [('exon', x, y) for x, y in zip(b, e)])

def selftest():
    K = ((100, 200), (300, 400), (500, 600)); M = K[1:]
    cnt = {K: 10, M: 6}; lo = {M: {250: 5, 150: 1}}; hi = {}
    ins, peak, mode = features(M, K, 1, True, cnt, lo, hi)
    assert (round(ins, 2), round(peak, 2), mode) == (0.83, 0.83, 250)          # one read starts before K's exon 2 (200..300)
    P = K[:2]; ins, peak, mode = features(P, K, 2, False, {P: 4}, {}, {P: {450: 2, 470: 1, 520: 1}})
    assert (ins, peak, mode) == (0.75, 0.5, 450)                                # one read ends past the exon's far boundary (500)
    assert cont_counts({K: 10, M: 6}, {(M, True): 1})[(M, True)] == 10
    g = gtf_lines('c', 't', '-', M, 250, 700).split('\n')
    assert g[0].split('\t')[3:5] == ['251', '700'] and g[1].split('\t')[3:5] == ['251', '300'] and g[2].split('\t')[3:5] == ['401', '500'] and gtf_lines('c', 't', '-', M, 350, 700) is None
    print('selftest ok')

if __name__ == '__main__':
    if sys.argv[1] == '--selftest': selftest(); sys.exit()
    mp = sys.argv.pop(sys.argv.index('--models') + 1) if '--models' in sys.argv else None
    if mp: sys.argv.remove('--models')
    ref_p, def_p, fl_p, bam, fa_p, out = sys.argv[1:7]; thr = int(sys.argv[7]) if len(sys.argv) > 7 else 16
    only = set(sys.argv[8].split(',')) if len(sys.argv) > 8 else None
    import pysam
    D = collections.defaultdict(dict)
    for (c, s), d in load(def_p, 'reads').items():
        for ch, v in d.items(): D[c][ch] = (s, v[0][2], min(t[0] for t in v), max(t[1] for t in v))
    ref = collections.defaultdict(dict); fl = collections.defaultdict(set)
    for (c, s), d in load(ref_p).items(): ref[c][s] = set(d)
    for (c, s), d in load(fl_p).items(): fl[c] |= set(d)
    jobs = [(bam, fa_p, c, D[c], ref[c], fl[c]) for c in sorted(D, key=lambda c: -len(D[c])) if c in set(pysam.AlignmentFile(bam).references) and (not only or c in only)]
    with multiprocessing.Pool(thr) as pool: res = list(pool.imap_unordered(work, jobs))
    rows = [r for rs, _ in res for r in rs]
    if mp:
        with open(mp, 'w') as o: o.write(''.join(g for _, gs in res for g in gs if g))
    H = 'contig strand class emitted annotated floor1 reads k_reads k_attr cont inside peak third n_introns k_introns'.split()
    with open(out, 'w') as o:
        o.write('\t'.join(H) + '\n'); [o.write('\t'.join(map(str, r)) + '\n') for r in sorted(rows)]
    R = [dict(zip(H, r)) for r in rows]; sh = lambda xs, f: f'{sum(map(f, xs)) / max(1, len(xs)):.3f}'
    for cls, third in (('NEST5', 'own'), ('NEST3', 'polyA')):
        X = [r for r in R if r['class'] == cls]; g = lambda e, a: [r for r in X if r['emitted'] == e and r['annotated'] == a]
        print(f"## {cls}: candidates {len(X)}; emitted annotated {len(g(1, 1))}, emitted unannotated {len(g(1, 0))}, not emitted annotated {len(g(0, 1))}, not emitted unannotated {len(g(0, 0))}")
        print(f'group\tn\tinside\tpeaked\t{third}\tfloor1\tpooled>=reads\tmedian reads')
        for name, xs in (('annotated kept', g(1, 1)), ('annotated lost', g(0, 1)), ('unannotated kept', g(1, 0)), ('unannotated not emitted', g(0, 0))):
            print(f"{name}\t{len(xs)}\t{sh(xs, lambda r: r['inside'] >= 0.9)}\t{sh(xs, lambda r: r['peak'] >= 0.5)}\t{sh(xs, lambda r: r['third'])}\t{sh(xs, lambda r: r['floor1'])}\t"
                  f"{sh(xs, lambda r: r['k_attr'] - r['k_reads'] >= r['reads'])}\t{sorted(r['reads'] for r in xs)[len(xs) // 2] if xs else 0}")
        print(f'signature (inside, peaked, {third})\tnot emitted\tannotated\tshare annotated\tshare of annotated lost')
        ne = [r for r in X if not r['emitted']]; lost = max(1, len(g(0, 1)))
        for i in (1, 0):
            for p in (1, 0):
                for t in (1, 0):
                    xs = [r for r in ne if (r['inside'] >= 0.9) == i and (r['peak'] >= 0.5) == p and r['third'] == t]; a = sum(r['annotated'] for r in xs)
                    print(f'{i}{p}{t}\t{len(xs)}\t{a}\t{a / max(1, len(xs)):.3f}\t{a / lost:.3f}')
