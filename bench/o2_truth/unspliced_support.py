#!/usr/bin/env python3
"""Unspliced reads inside the exons of spliced transcripts (register 1989).

    unspliced_support.py GTF BAM OUT.tsv [THREADS]

One line per multi-exon transcript of GTF: unspliced primary reads lying wholly inside one of its end exons / internal exons
(SLACK bp allowed). Last lines (#): where the unspliced primary reads lie against the transcripts WITHOUT a one_read_locus tag."""
import sys, os, bisect, collections, multiprocessing
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../tools'))
import oneread
SLACK = 4
_G = {}

def inside(starts, exons, a, b):
    """indices of the sorted exons that hold [a, b) with SLACK; exons are short, so a bounded walk back is enough"""
    i = bisect.bisect_right(starts, a + SLACK) - 1; out = []
    while i >= 0 and starts[i] >= a - _G['maxlen']:
        if exons[i][1] >= b - SLACK: out.append(i)
        i -= 1
    return out

def _scan(contig):
    import pysam
    ex = _G['exons'].get(contig, []); st = [e[0] for e in ex]; n = collections.Counter(); base = _G['base'].get(contig, ([], [], [])); where = collections.Counter()
    for r in pysam.AlignmentFile(_G['bam']).fetch(contig):
        if r.is_secondary or r.is_supplementary or r.is_unmapped or 'N' in r.cigarstring: continue
        a, b = r.reference_start, r.reference_end
        for i in inside(st, ex, a, b): n[ex[i]] += 1
        bs, be, bx = base
        if inside(bs, be, a, b) if bs else False: where['inside an exon of a spliced transcript'] += 1
        else:
            j = bisect.bisect_right(bx[0], a) - 1 if bx else -1
            where['inside the span of a spliced transcript, not inside an exon' if j >= 0 and bx[1][j] >= b else 'elsewhere'] += 1
    return contig, n, where

def main(gtf, bam, out, thr=32):
    import pysam
    M = {t: v for t, v in oneread.models(gtf).items() if len(v[2]) > 1}
    exons = collections.defaultdict(set); bex = collections.defaultdict(set); span = collections.defaultdict(list)
    for c, s, ex, lines in M.values():
        exons[c].update(ex)
        if 'one_read_locus' not in lines[0]: bex[c].update(ex); span[c].append((ex[0][0], ex[-1][1]))
    base = {}
    for c in bex:
        e = sorted(bex[c]); sp = sorted(span[c]); mx = []; m = 0      # spans: running max of ends, so "inside some span" is one lookup
        for a, b in sp: m = max(m, b); mx.append(m)
        base[c] = ([x[0] for x in e], e, ([x[0] for x in sp], mx))
    _G.update(bam=bam, exons={c: sorted(v) for c, v in exons.items()}, base=base, maxlen=max(b - a for v in exons.values() for a, b in v))
    with multiprocessing.get_context('fork').Pool(int(thr)) as pool: res = pool.map(_scan, [c for c in pysam.AlignmentFile(bam).references], chunksize=1)
    n = {c: k for c, k, _ in res}; where = sum((w for _, _, w in res), collections.Counter())
    with open(out, 'w') as o:
        o.write('tid\tunspliced_end_exons\tunspliced_internal_exons\n')
        for t, (c, s, ex, _) in M.items(): k = n.get(c, {}); o.write(f"{t}\t{k.get(ex[0], 0) + k.get(ex[-1], 0)}\t{sum(k.get(e, 0) for e in ex[1:-1])}\n")
        for w, v in where.most_common(): o.write(f'# unspliced primary reads {w}\t{v}\n')

def selftest():
    _G['maxlen'] = 100; ex = [(0, 50), (10, 100), (200, 300)]; st = [e[0] for e in ex]
    assert sorted(inside(st, ex, 12, 48)) == [0, 1] and inside(st, ex, 60, 90) == [1] and inside(st, ex, 90, 210) == [] and inside(st, ex, 198, 303) == [2]
    print('unspliced_support selftest ok')

if __name__ == '__main__': selftest() if sys.argv[1:] == ['--selftest'] else main(*sys.argv[1:])
