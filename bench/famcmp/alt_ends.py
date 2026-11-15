#!/usr/bin/env python3
'''Alternative 3' ends as attributes (opt-in driver stage `altends`; docs/STEP0_tool_match_2026-11-11.md §11).

  alt_ends.py GTF BAM FASTA OUT_GTF [THREADS]        (selftest: alt_ends.py selftest)
OUT_GTF = GTF with, on a multi-exon transcript line whose commonest full-chain primary read 3' end lies more than MIN_DIST bp
from the emitted 3' end and has >= MIN_READS reads: alt_tes (1-based), alt_tes_reads "n/total", alt_tes_polya_signal
(AATAAA/ATTAAA 5-40 bp upstream), alt_tes_a_rich (genome after it: >= 12 A in 20 bp or 6 in a row = internal-priming
signature). No coordinate changes; every other line is copied unchanged.'''
import sys, re, collections
MIN_DIST, MIN_READS = 10, 2
TR = str.maketrans('ACGTN', 'TGCAN')

def alt_end(ends, plus, emitted):
    """ends: 3' end of each full-chain read (0-based: exclusive end on +, start on -) -> (end, reads there, reads) or None."""
    if not ends: return None
    c = collections.Counter(ends); m = min(c, key=lambda x: (-c[x], -x if plus else x))  # ties to the more extended
    return (m, c[m], len(ends)) if abs(m - emitted) > MIN_DIST and c[m] >= MIN_READS else None

def flags(fetch, contig, plus, pos):
    """(polyA signal upstream, A-rich downstream) in transcript sense; fetch(contig, start, end) -> sequence."""
    if plus: up, down = fetch(contig, max(pos - 40, 0), max(pos - 5, 0)).upper(), fetch(contig, pos, pos + 20).upper()
    else: up, down = fetch(contig, pos + 5, pos + 40).upper().translate(TR)[::-1], fetch(contig, max(pos - 20, 0), pos).upper().translate(TR)[::-1]
    return 'AATAAA' in up or 'ATTAAA' in up, down.count('A') >= 12 or 'AAAAAA' in down

def attrs(plus, alt, sig, arich):
    pos, n, tot = alt; b = lambda x: 'true' if x else 'false'
    return f' alt_tes "{pos if plus else pos + 1}"; alt_tes_reads "{n}/{tot}"; alt_tes_polya_signal "{b(sig)}"; alt_tes_a_rich "{b(arich)}";'

def work(job):
    import pysam
    bam, contig, want = job; g = collections.defaultdict(list)
    for r in pysam.AlignmentFile(bam).fetch(contig):
        if r.is_secondary or r.is_supplementary or r.is_unmapped or 'N' not in r.cigarstring: continue
        pos = r.reference_start; ch = []
        for op, ln in r.cigartuples:
            if op == 3: ch.append((pos, pos + ln))
            if op in (0, 2, 3, 7, 8): pos += ln
        ch = tuple(ch)
        if ch in want: g[ch].append((r.reference_start, pos))
    return contig, dict(g)

if __name__ == '__main__':
    if sys.argv[1] == 'selftest':
        assert alt_end([900] * 3 + [1000], True, 1000) == (900, 3, 4)        # commonest end upstream of the emitted one
        assert alt_end([995] * 3 + [1000], True, 1000) is None               # within MIN_DIST
        assert alt_end([900, 950, 1000], True, 1000) is None                 # no end with 2 reads away from the emitted one
        assert alt_end([200] * 2 + [100] * 2, False, 100) is None and alt_end([200] * 3 + [100] * 2, False, 100) == (200, 3, 5)  # - strand: tie to the more extended (smaller)
        seq = 'C' * 60 + 'AATAAA' + 'C' * 14 + 'A' * 20 + 'C' * 100; f = lambda c, a, b: seq[a:b]
        assert flags(f, 'c', True, 80) == (True, True) and flags(f, 'c', True, 150) == (False, False)
        rc = seq.translate(TR)[::-1]; f2 = lambda c, a, b: rc[a:b]; assert flags(f2, 'c', False, len(seq) - 80) == (True, True)
        assert attrs(False, (200, 3, 5), True, False) == ' alt_tes "201"; alt_tes_reads "3/5"; alt_tes_polya_signal "true"; alt_tes_a_rich "false";'
        print('ok'); sys.exit()
    import multiprocessing, pysam
    gtf, bam, fa, out = sys.argv[1:5]; thr = int(sys.argv[5]) if len(sys.argv) > 5 else 8
    ex = collections.defaultdict(list); tid = lambda a: re.search(r'transcript_id "([^"]*)"', a).group(1)
    for l in open(gtf):
        f = l.split('\t')
        if len(f) > 8 and f[2] == 'exon' and f[6] in '+-': ex[(f[0], f[6], tid(f[8]))].append((int(f[3]) - 1, int(f[4])))
    model = {}; want = collections.defaultdict(set)   # transcript -> (contig, plus, chain, emitted 3' end)
    for (c, s, t), e in ex.items():
        e.sort()
        if len(e) > 1: ch = tuple((a[1], b[0]) for a, b in zip(e, e[1:])); model[t] = (c, s == '+', ch, e[-1][1] if s == '+' else e[0][0]); want[c].add(ch)
    have = set(pysam.AlignmentFile(bam).references)
    with multiprocessing.Pool(thr) as p: reads = dict(p.imap_unordered(work, [(bam, c, w) for c, w in want.items() if c in have]))
    g = pysam.FastaFile(fa); n = sig = ar = 0
    with open(out, 'w') as w:
        for l in open(gtf):
            f = l.rstrip('\n').split('\t')
            if len(f) > 8 and f[2] == 'transcript' and tid(f[8]) in model:
                c, plus, ch, em = model[tid(f[8])]; a = alt_end([(b if plus else a_) for a_, b in reads.get(c, {}).get(ch, ())], plus, em)
                if a:
                    s_, r_ = flags(g.fetch, c, plus, a[0]); n += 1; sig += s_; ar += r_; l = l.rstrip('\n') + attrs(plus, a, s_, r_) + '\n'
            w.write(l)
    print(f'{len(model)} multi-exon models; {n} with an alternative 3\' end ({sig} with a polyA signal, {ar} A-rich downstream)')
