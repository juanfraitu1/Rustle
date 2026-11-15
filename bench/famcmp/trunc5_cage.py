#!/usr/bin/env python3
"""PREREG row 1742: CAGE support of the 5' ends of chains the `trunc5` stage drops vs keeps.

  trunc5_cage.py PREFIX.gtf PREFIX.trunc5.gtf CAGE_BED6 REF_GTF   |   trunc5_cage.py --selftest
Multi-exon transcripts only. Background = the position 100 bp downstream in the transcript."""
import sys, re, collections
from sq_nested import tss_and_shift, Peaks

def load(path):
    """-> {transcript_id: (contig, strand, sorted exons 0-based half-open)}"""
    ex = collections.defaultdict(list); meta = {}
    for l in open(path):
        f = l.split('\t')
        if len(f) > 8 and f[2] == 'exon':
            t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1); ex[t].append((int(f[3]) - 1, int(f[4]))); meta[t] = (f[0], f[6])
    return {t: (*meta[t], sorted(e)) for t, e in ex.items()}

def chain(c, s, e): return (c, s, tuple((e[i][1], e[i + 1][0]) for i in range(len(e) - 1)))

def score(models, kept_ids, peaks, ref_chains):
    T = collections.Counter()
    for t, (c, s, e) in models.items():
        if len(e) < 2: continue
        tss, sh = tss_and_shift(s, e); d = peaks.dist(c, s, tss)
        for g in (('kept' if t in kept_ids else 'dropped', 'all'), ('kept' if t in kept_ids else 'dropped', 'refseq' if chain(c, s, e) in ref_chains else 'not_refseq')):
            T[g, 'n'] += 1; T[g, 'in'] += d == 0; T[g, 'w50'] += d <= 50
            if sh is not None: T[g, 'bn'] += 1; T[g, 'bin'] += peaks.dist(c, s, sh) == 0
    return T

def selftest():
    import tempfile, os
    d = tempfile.mkdtemp(); bed = os.path.join(d, 'p.bed'); open(bed, 'w').write('c\t100\t110\tp\t1\t+\n')
    m = {'a': ('c', '+', [(105, 300), (400, 500)]), 'b': ('c', '+', [(5, 300), (400, 500)]), 'm': ('c', '+', [(105, 300)])}
    T = score(m, {'a'}, Peaks(bed), {('c', '+', ((300, 400),))})
    assert T[('kept', 'all'), 'n'] == 1 and T[('kept', 'all'), 'in'] == 1 and T[('dropped', 'all'), 'in'] == 0
    assert T[('dropped', 'all'), 'bin'] == 1 and T[('kept', 'refseq'), 'n'] == 1   # b's shifted position 105 is in the peak
    print('selftest ok')

if __name__ == '__main__':
    if sys.argv[1] == '--selftest': selftest(); sys.exit()
    gtf, tr, bed, ref = sys.argv[1:5]
    models = load(gtf); kept = set(load(tr)); refc = {chain(*v) for v in load(ref).values() if len(v[2]) > 1}
    T = score(models, kept, Peaks(bed), refc)
    print('group\tsplit\tchains\tin peak\twithin 50 bp\tbackground in peak\tratio')
    for g in ('dropped', 'kept'):
        for sp in ('all', 'refseq', 'not_refseq'):
            k = (g, sp); n = max(1, T[k, 'n']); a = T[k, 'in'] / n; b = T[k, 'bin'] / max(1, T[k, 'bn'])
            print(f"{g}\t{sp}\t{T[k, 'n']}\t{a:.4f}\t{T[k, 'w50'] / n:.4f}\t{b:.4f}\t{a / b if b else float('inf'):.2f}")
