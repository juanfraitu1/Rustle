#!/usr/bin/env python3
"""3' ends: the assembled end against the unspliced reads in the last exon (register 1989, addendum B).

    unspliced_ends.py GTF BAM TRUTH.gtf [THREADS]

Units: multi-exon transcripts of GTF (no one_read_locus tag) whose intron chain is annotated, with >= MIN unspliced primary
reads lying in the last exon (not crossing its splice site by more than SLACK; they may run past the assembled end)."""
import sys, os, bisect, collections, multiprocessing, statistics
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../tools'))
import oneread
SLACK, MIN = 4, 3
_G = {}

def _scan(contig):
    import pysam
    U = _G['units'].get(contig, []); out = []
    if not U: return out
    B = pysam.AlignmentFile(_G['bam'])
    for s, (a, b), asm, ann in U:
        ends = []
        for r in B.fetch(contig, a, b):
            if r.is_secondary or r.is_supplementary or 'N' in r.cigarstring: continue
            if s == '+' and r.reference_start >= a - SLACK: ends.append(r.reference_end)
            if s == '-' and r.reference_end <= b + SLACK: ends.append(r.reference_start)
        if len(ends) >= MIN:
            m = statistics.median(ends); near = lambda x: min(abs(x - e) for e in ann)
            out.append((near(asm), near(m), len(ends)))
    return out

def main(gtf, bam, truth, thr=32):
    import pysam
    ann = collections.defaultdict(set)
    for c, s, ex, _ in oneread.models(truth).values():
        if len(ex) > 1: ann[(c, s, oneread.chain(ex))].add(ex[-1][1] if s == '+' else ex[0][0])
    units = collections.defaultdict(list); n = 0
    for c, s, ex, lines in oneread.models(gtf).values():
        if len(ex) < 2 or 'one_read_locus' in lines[0]: continue
        k = (c, s, oneread.chain(ex))
        if k in ann: n += 1; units[c].append((s, ex[-1] if s == '+' else ex[0], ex[-1][1] if s == '+' else ex[0][0], sorted(ann[k])))
    _G.update(bam=bam, units=units)
    with multiprocessing.get_context('fork').Pool(int(thr)) as pool: R = [x for v in pool.map(_scan, list(units), chunksize=1) for x in v]
    print(f'transcripts with an annotated chain {n}; with >= {MIN} unspliced reads in the last exon {len(R)} ({len(R) / max(1, n):.3f})')
    print("within\tassembled 3' end\tmedian end of the unspliced reads\tunspliced better / same / worse (by > 10 bp)")
    for w in (10, 20, 50, 100):
        print(w, f'{sum(a <= w for a, u, _ in R) / max(1, len(R)):.3f}', f'{sum(u <= w for a, u, _ in R) / max(1, len(R)):.3f}', sep='\t', end='\t' if w == 10 else '\n')
        if w == 10: print(sum(u < a - 10 for a, u, _ in R), sum(abs(u - a) <= 10 for a, u, _ in R), sum(u > a + 10 for a, u, _ in R), sep=' / ')

if __name__ == '__main__': main(*sys.argv[1:])
