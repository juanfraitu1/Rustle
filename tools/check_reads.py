#!/usr/bin/env python3
"""Pre-flight check of the reads (driver: start of `assemble`, warnings only).

  check_reads.py BAM FASTA [N=100000]

Reads the first N primary alignments and warns when the reads look untrimmed or unoriented. The pipeline takes the
transcript strand as the alignment strand flipped by minimap2's `ts:A:-`, which is right only for reads already in
transcript orientation with the polyA tail removed (FLNC). Exit status is always 0."""
import sys, collections
RC = str.maketrans('ACGT', 'TGCA')
SENSE, ANTI = {'GTAG', 'GCAG', 'ATAC'}, {'CTAC', 'CTGC', 'GTAT'}

def read_flags(reverse, seq, motifs, ts):
    """seq: SEQ field; motifs: intron first 2 + last 2 genome bases, + strand. Flags of one primary read."""
    own = seq[::-1].translate(RC) if reverse else seq  # the read as sequenced
    if reverse: motifs = [m[::-1].translate(RC) for m in motifs]
    s = sum(m in SENSE for m in motifs); a = sum(m in ANTI for m in motifs)
    return dict(tail=own[-20:].count('A') >= 16, polyt='T' * 12 in own[:30], spliced=bool(motifs), anti=a > s,
                ts_minus=ts == '-', unspliced_no_ts=not motifs and ts is None)

def warnings(T, n):
    """T: flag counts over n reads."""
    w = []; sp = max(1, T['spliced']); f = lambda k, d=n: f'{T[k]} of {d} ({T[k] / d:.1%})'
    if T['tail'] >= 0.01 * n: w.append(f"{f('tail')} reads end in a polyA run: tails look untrimmed (run isoseq refine / a polyA trimmer before aligning)")
    if T['polyt'] >= 0.01 * n: w.append(f"{f('polyt')} reads start with a polyT run: reads look unoriented (reversed molecules keep their tail at the 5' end)")
    if T['anti'] >= 0.02 * sp: w.append(f"{f('anti', sp)} spliced reads have mostly antisense introns: reads look unoriented but were aligned as oriented (minimap2 -uf); those reads are on the wrong strand")
    if T['ts_minus'] >= 0.02 * sp: w.append(f"{f('ts_minus', sp)} spliced reads carry ts:A:-: reads are not in transcript orientation; the {T['unspliced_no_ts']} unspliced reads have no strand evidence and are taken at read orientation")
    return w

def selftest():
    assert read_flags(False, 'C' * 40 + 'A' * 20, ['GTAG'], '+')['tail'] and not read_flags(False, 'C' * 60, [], None)['tail']
    r = read_flags(True, 'C' * 40 + 'A' * 20, ['CTAC', 'CTAC'], '+')  # reverse alignment: own start is polyT, introns sense
    assert r['polyt'] and not r['anti'] and not r['tail']
    assert read_flags(False, 'T' * 20 + 'C' * 40, ['CTAC'], '+')['anti'] and read_flags(False, 'C' * 60, [], None)['unspliced_no_ts']
    T = collections.Counter(tail=5, spliced=900, anti=30); assert len(warnings(T, 1000)) == 1 and 'antisense' in warnings(T, 1000)[0]
    assert not warnings(collections.Counter(tail=1, polyt=1, spliced=900, anti=1), 1000)
    print('selftest ok')

if __name__ == '__main__':
    if sys.argv[1:] == ['--selftest']: selftest(); sys.exit(0)
    import pysam
    bam, fa_p = sys.argv[1:3]; N = int(sys.argv[3]) if len(sys.argv) > 3 else 100000
    fa = pysam.FastaFile(fa_p); T = collections.Counter(); n = 0
    for r in pysam.AlignmentFile(bam):
        if r.is_secondary or r.is_supplementary or r.is_unmapped or not r.query_sequence: continue
        pos = r.reference_start; motifs = []
        for op, ln in r.cigartuples:
            if op == 3: motifs.append((fa.fetch(r.reference_name, pos, pos + 2) + fa.fetch(r.reference_name, pos + ln - 2, pos + ln)).upper())
            if op in (0, 2, 3, 7, 8): pos += ln
        ts = r.get_tag('ts') if r.has_tag('ts') else None
        T.update(k for k, v in read_flags(r.is_reverse, r.query_sequence, motifs, ts).items() if v); n += 1
        if n >= N: break
    w = warnings(T, n) if n else []
    for x in w: print(f'[check_reads] WARNING: {x}', file=sys.stderr)
    print(f"[check_reads] {n} primary reads: polyA tail {T['tail']}, leading polyT {T['polyt']}, spliced {T['spliced']}, antisense-majority {T['anti']}, ts:- {T['ts_minus']}; {len(w)} warning(s)", file=sys.stderr)
