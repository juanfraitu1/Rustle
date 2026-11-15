#!/usr/bin/env python3
"""Step 0b (tool-matching goal, 2026-11-11): full-chain primary-record support of every chain ours, StringTie or FLAIR output.

  tool_match_reads.py REFSEQ OURS_GTF STRINGTIE_GTF FLAIR_GTF BAM OUT_TSV [THREADS]

OUT_TSV: contig, strand, n_introns, reads (primary records whose intron chain equals the chain), in_ours, in_stringtie, in_flair,
in_refseq. Stdout: per tool, the tool-only / shared / ours-only chains by read bin, with RefSeq-correct counts."""
import sys, collections, multiprocessing, pysam
sys.path.insert(0, __import__('os').path.dirname(__file__))
from tool_match import read

def count(job):
    bam_path, contig, want = job; n = collections.Counter()
    for r in pysam.AlignmentFile(bam_path).fetch(contig):
        if r.is_secondary or r.is_supplementary or 'N' not in r.cigarstring: continue
        pos = r.reference_start; ch = []
        for op, ln in r.cigartuples:
            if op == 3: ch.append((pos, pos + ln))
            if op in (0, 2, 3, 7, 8): pos += ln
        ch = tuple(ch)
        if ch in want: n[ch] += 1
    return contig, n

if __name__ == '__main__':
    ref, ours, st, fl = (read(p) for p in sys.argv[1:5]); bam, out = sys.argv[5:7]; thr = int(sys.argv[7]) if len(sys.argv) > 7 else 16
    allc = ours | st | fl; by = collections.defaultdict(set)
    for c, s, ch in allc: by[c].add(ch)
    with multiprocessing.Pool(thr) as p: got = dict(p.imap_unordered(count, [(bam, c, w) for c, w in by.items()]))
    bins = lambda n: '0' if n == 0 else '1' if n == 1 else '2' if n == 2 else '3-4' if n < 5 else '5+'
    tab = collections.defaultdict(lambda: [0, 0])
    with open(out, 'w') as w:
        for c, s, ch in sorted(allc):
            n = got[c][ch]; k = (c, s, ch); w.write(f'{c}\t{s}\t{len(ch)}\t{n}\t{int(k in ours)}\t{int(k in st)}\t{int(k in fl)}\t{int(k in ref)}\n')
            for name, t in (('stringtie', st), ('flair', fl)):
                cls = 'shared' if k in ours and k in t else 'tool_only' if k in t else 'ours_only' if k in ours else None
                if cls: v = tab[(name, cls, bins(n))]; v[0] += 1; v[1] += k in ref
    print('tool\tclass\treads\tchains\tcorrect\tshare')
    for (name, cls, b), (n, ok) in sorted(tab.items()): print(f'{name}\t{cls}\t{b}\t{n}\t{ok}\t{ok / n:.3f}')
