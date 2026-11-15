#!/usr/bin/env python3
"""Driver stage `shortexons` (register 1753 / 1754 / 1757): rewrite reads that skip an exon of <= 15 bp through that exon.

  short_exons.py BAM FASTA OUT_BAM OUT_PREFIX [THREADS] [MIN_READS=2]   |   short_exons.py --selftest
Two annotation-free rules give the units (donor, exon, acceptor): `junction` = another primary read of the locus spells the
exon (skip_fix.py), `search` = the read's extra bases at the junction match one splice-flanked site of the intron, on at
least MIN_READS reads (skip_search.py). OUT_BAM = BAM with those reads' CIGARs rewritten (same order, primary alignments
only, unmapped reads left out); OUT_PREFIX.units.tsv and OUT_PREFIX.reads.tsv say what was changed."""
import sys, os, collections, multiprocessing, functools
from skip_fix import fix_all, primary_spliced
import skip_search

def pick(rows, floor):
    """scan rows (contig, d, x, y, a, reads with the exon, reads with the skip, reads the search placed) -> units with their rule"""
    return sorted((*r[:5], 'search' if r[7] >= floor else 'junction', *r[5:]) for r in rows if r[7] >= floor or (r[5] and r[6]))

def rewrite(job):
    bam, fa_p, contig, cands, tmp = job; import pysam
    fa = pysam.FastaFile(fa_p); fetch = functools.lru_cache(200000)(lambda s, e: fa.fetch(contig, s, e).upper())
    src = pysam.AlignmentFile(bam); log = []
    with pysam.AlignmentFile(tmp, 'wb', template=src) as dst:
        for r in src.fetch(contig):
            if cands and primary_spliced(r) and r.query_sequence:
                cig, done = fix_all(r.cigartuples, r.reference_start, r.query_sequence, cands, fetch)
                if done:
                    old = r.cigarstring; r.cigartuples = cig
                    for t in ('NM', 'MD', 'cs'):                # stale after the rewrite (the assembler reads AS only)
                        if r.has_tag(t): r.set_tag(t, None)
                    log.append((r.query_name, contig, r.reference_start, old, r.cigarstring, ';'.join('%d-%d' % (x, y) for _, x, y, _ in done)))
            dst.write(r)
    return log

if __name__ == '__main__':
    if sys.argv[1] == '--selftest':
        rows = [('c', 1, 2, 3, 4, 0, 9, 2), ('c', 1, 5, 6, 9, 3, 7, 0), ('c', 1, 7, 8, 9, 0, 7, 1), ('c', 2, 3, 4, 5, 4, 6, 5)]
        assert [u[5] for u in pick(rows, 2)] == ['search', 'junction', 'search'] and len(pick(rows, 1)) == 4   # one read: under the floor
        print('selftest ok'); sys.exit()
    import pysam
    bam, fa, out, pre = sys.argv[1:5]; thr = int(sys.argv[5]) if len(sys.argv) > 5 else 16; floor = int(sys.argv[6]) if len(sys.argv) > 6 else 2
    refs = pysam.AlignmentFile(bam).references; rows = []
    with multiprocessing.Pool(thr) as pool:
        for r, _ in pool.imap_unordered(skip_search.scan, [(bam, fa, c) for c in refs]): rows += r
        units = pick(rows, floor); cands = collections.defaultdict(lambda: collections.defaultdict(list))
        for c, d, x, y, a, *_ in units: cands[c][d, a].append((x, y))
        tmps = [f'{out}.{i}.tmp.bam' for i in range(len(refs))]
        logs = pool.map(rewrite, [(bam, fa, c, dict(cands[c]) if c in cands else None, t) for c, t in zip(refs, tmps)])
    pysam.cat('-o', out, *tmps); pysam.index(out)
    for t in tmps: os.remove(t)
    with open(pre + '.units.tsv', 'w') as o:
        o.write('contig\tdonor\texon_start\texon_end\tacceptor\trule\treads_with_exon\treads_with_skip\treads_placed_by_search\n')
        for u in units: o.write('\t'.join(map(str, u)) + '\n')
    n = 0
    with open(pre + '.reads.tsv', 'w') as o:
        o.write('read\tcontig\tstart\told_cigar\tnew_cigar\texons\n')
        for l in logs:
            for r in l: o.write('\t'.join(map(str, r)) + '\n'); n += 1
    k = collections.Counter(u[5] for u in units)
    print(f"{len(units)} short-exon units ({k['search']} by intron search, {k['junction']} by another read's junctions), {n} reads rewritten")
