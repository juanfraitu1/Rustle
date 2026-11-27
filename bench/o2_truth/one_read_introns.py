#!/usr/bin/env python3
"""One-read loci judged by their introns, and by the sample's other reads (register 1987, second addendum).

    one_read_introns.py ONELOCI.gtf BAM TRUTH.gtf REP1.bam[,REP2.bam] OUT.tsv [THREADS]

One line per multi-exon transcript of ONELOCI.gtf: class (one_read_locus or the transcript's reads), introns, how many are
in a replicate's primary reads, in the annotation, and in >= 2 primary reads of BAM itself (the lone read + another)."""
import sys, os, re, collections, multiprocessing
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../tools'))
import oneread, one_read_loci as O
_G = {}

def _introns(a):
    import pysam
    bam, contig = a; n = collections.Counter(); want = _G['want'].get(contig, ())
    if not want: return contig, n
    for r in pysam.AlignmentFile(bam).fetch(contig):
        if r.is_secondary or r.is_supplementary or 'N' not in r.cigarstring: continue
        for i in O.read_chain(r):
            if i in want: n[i] += 1
    return contig, n

def counts(bam, pool):
    import pysam
    return dict(pool.imap_unordered(_introns, [(bam, c) for c in pysam.AlignmentFile(bam).references if c in _G['want']]))

def main(gtf, bam, truth, reps, out, thr=32):
    M = {t: v for t, v in oneread.models(gtf).items() if len(v[2]) > 1}
    want = collections.defaultdict(set)
    for c, s, ex, _ in M.values(): want[c].update(oneread.chain(ex))
    _G['want'] = want
    tr = collections.defaultdict(set)
    for c, s, ex, _ in oneread.models(truth).values(): tr[c].update(oneread.chain(ex))
    with multiprocessing.get_context('fork').Pool(int(thr)) as pool:
        own = counts(bam, pool); rep = [counts(b, pool) for b in reps.split(',')]
    with open(out, 'w') as o:
        o.write('tid\tclass\tintrons\tin_replicate\tannotated\tin_two_own_reads\tmin_own_reads\n')
        for t, (c, s, ex, lines) in M.items():
            ch = oneread.chain(ex); m = re.search(r' reads "([0-9.]+)"', lines[0])
            k = 'one' if 'one_read_locus' in lines[0] else m.group(1) if m else 'NA'; ow = [own.get(c, {}).get(i, 0) for i in ch]
            o.write(f"{t}\t{k}\t{len(ch)}\t{sum(any(r.get(c, {}).get(i, 0) for r in rep) for i in ch)}\t{sum(i in tr[c] for i in ch)}\t{sum(x >= 2 for x in ow)}\t{min(ow)}\n")

if __name__ == '__main__': main(*sys.argv[1:])
