#!/usr/bin/env python3
"""One-read chains with junction support and a clean reading frame (register 1766, rule S3b; driver stage `oneread`).

    oneread.py DEFAULT.gtf FLOOR1.gtf BAM GENOME.fa OUT.gtf   |   oneread.py --selftest

FLOOR1.gtf = the same assembly at a 1-read floor. OUT.gtf holds the multi-exon models of FLOOR1 whose intron chain DEFAULT lacks and
that pass every test (primary alignments only, intron chains compared without strand):
  one read has exactly the chain, and no read has a longer chain containing it as a contiguous part
  every junction is also in another read
  on every junction, more than 0.05 of the reads are colinear with the chain (their chain is a contiguous part of it)
  the longest ATG..stop ORF is >= 900 nt and does not stop > 50 nt before the last junction (gtf_ptc.py)
Admitted transcript lines get oneread "1". Prints one summary line (the driver logs it)."""
import collections, os, re, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gtf_ptc

MIN_ORF, MIN_SHARE = 900, 0.05

def run(a, b):
    """a is a contiguous part of b"""
    if len(a) > len(b) or a[0] not in b: return False
    k = b.index(a[0]); return b[k:k + len(a)] == a

def junction_ok(C, others):
    """others = {read chain: reads} for the reads sharing an intron with C. -> True when C passes the read tests."""
    if others.get(C, 0) != 1: return False
    cls = {r: 'full' if r == C else 'agree' if run(r, C) else 'longer' if run(C, r) else 'conflict' for r in others}
    if 'longer' in cls.values(): return False
    for i in C:
        on = [r for r in others if i in r]; tot = sum(others[r] for r in on)
        if tot < 2 or sum(others[r] for r in on if cls[r] != 'conflict') / tot <= MIN_SHARE: return False
    return True

def models(path):
    """-> ordered {tid: (chrom, strand, sorted 0-based half-open exons, lines)}"""
    T = collections.OrderedDict()
    for l in open(path):
        f = l.split('\t')
        if len(f) < 9 or f[2] not in ('transcript', 'exon'): continue
        t = T.setdefault(re.search(r'transcript_id "([^"]+)"', f[8]).group(1), [f[0], f[6], [], []])
        t[3].append(l)
        if f[2] == 'exon': t[2].append((int(f[3]) - 1, int(f[4])))
    for t in T.values(): t[2].sort()
    return T

chain = lambda ex: tuple((ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1))

def main(default, floor1, bam, fasta, out):
    import pysam
    have = {(c, chain(ex)) for c, s, ex, _ in models(default).values() if len(ex) > 1}
    cand = collections.OrderedDict()
    for tid, (c, s, ex, lines) in models(floor1).items():
        k = (c, chain(ex))
        if len(ex) > 1 and s in '+-' and k not in have and k not in cand: cand[k] = (s, ex, lines)
    want = {(c, i) for c, ch in cand for i in ch}; prim = collections.Counter()
    for a in pysam.AlignmentFile(bam, threads=4):
        if a.is_unmapped or a.is_secondary or a.is_supplementary: continue
        p = a.reference_start; ch = []
        for op, n in a.cigartuples:
            if op == 3: ch.append((p, p + n))
            if op in (0, 2, 3, 7, 8): p += n
        if ch and any((a.reference_name, i) in want for i in ch): prim[(a.reference_name, tuple(ch))] += 1
    has = collections.defaultdict(set)
    for c, ch in prim:
        for i in ch: has[(c, i)].add(ch)
    fa = pysam.FastaFile(fasta); n = 0
    with open(out, 'w') as w:
        for (c, C), (s, ex, lines) in cand.items():
            if prim.get((c, C), 0) != 1: continue
            if not junction_ok(C, {r: prim[(c, r)] for i in C for r in has[(c, i)]}): continue
            o, p = gtf_ptc.call(ex, s, lambda a, b: fa.fetch(c, a, b))
            if p or o < MIN_ORF: continue
            n += 1
            for l in lines:
                # the 1-read run reuses DN_<contig>_<start>_<n> ids of other default transcripts: gtf_ptc joined their exons
                l = re.sub(r'(transcript_id "[^"]+)', r'\1.oneread', l, 1)
                if l.split('\t')[2] == 'transcript': l = l.rstrip('\n').rstrip() + ('' if l.rstrip().endswith(';') else ';') + ' oneread "1";\n'
                w.write(l)
    print('%d one-read chains admitted of %d chains the default lacks' % (n, len(cand)))

def selftest():
    C = ((10, 20), (30, 40), (50, 60)); main_iso = ((10, 20), (30, 60))
    ok = {C: 1, ((30, 40), (50, 60)): 1, ((10, 20), (30, 40)): 1, main_iso: 5}
    assert run(((30, 40),), C) and not run(main_iso, C)
    assert junction_ok(C, ok)
    assert not junction_ok(C, {C: 2})                                   # two full reads: not this rule's case
    assert not junction_ok(C, {C: 1, ((10, 20),): 3, ((30, 40),): 1})    # last junction only in the full read
    assert not junction_ok(C, {**ok, ((5, 8),) + C: 1})            # a longer read contains it
    assert not junction_ok(C, {C: 1, ((10, 20), (30, 40)): 1, ((30, 40), (50, 60)): 1, main_iso: 60})  # junction 1: 2 of 62 colinear
    print('oneread selftest ok')

if __name__ == '__main__':
    if sys.argv[1:] == ['--selftest']: selftest()
    elif len(sys.argv) == 6: main(*sys.argv[1:])
    else: sys.exit(__doc__)
