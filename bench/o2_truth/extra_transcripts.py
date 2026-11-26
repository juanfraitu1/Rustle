#!/usr/bin/env python3
"""What speaks for the unannotated transcripts an arm adds (register 1965).

    extra_transcripts.py ANNOTATION.gff GENOME.fa BAM TODAY.gtf ARM.gtf OTHER_BAM OTHER_TODAY.gtf [THREADS]

EXTRA = an intron chain of ARM.gtf that TODAY.gtf does not have and the annotation does not have. For comparison, the same
table for the unannotated chains of TODAY.gtf. Per chain: every intron annotated (a new combination of known introns) or
not; every intron GT-AG / GC-AG / AT-AC on its strand; reads of the sample's own with exactly the chain (primary, mapping
quality >= 1); the same in OTHER_BAM (a second individual); the chain in OTHER_TODAY.gtf (assembled there without the arm)."""
import sys, os, collections, multiprocessing
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import thin_isoforms as T
COMP = str.maketrans('ACGTacgt', 'TGCAtgca')

def canonical(fa, c, s, intron):
    a, b = intron; d = fa.fetch(c, a, a + 2).upper(); e = fa.fetch(c, b - 2, b).upper()
    if s == '-': d, e = e[::-1].translate(COMP), d[::-1].translate(COMP)
    return (d, e) in (('GT', 'AG'), ('GC', 'AG'), ('AT', 'AC'))

def own_reads(bam, chains, thr):
    import pysam
    by = collections.defaultdict(set)
    for c, s, ch in chains: by[c].add(ch)
    T.G.update(bam=bam, by=by)
    with multiprocessing.get_context('fork').Pool(thr) as pool: res = pool.map(T.scan, list(pysam.AlignmentFile(bam).references), chunksize=1)
    n = collections.Counter()
    for d in res: n.update(d)
    return n

def main(gff, fasta, bam, today, arm, obam, otoday, thr=32):
    import pysam
    thr = int(thr); ann = T.chains_gff(gff); known = {(c, s, i) for c, s, ch in ann for i in ch}; t0 = T.chains_gtf(today); a = T.chains_gtf(arm); o = T.chains_gtf(otoday)
    extra = a - t0 - ann; base = t0 - ann; both = extra | base
    n1 = own_reads(bam, both, thr); n2 = own_reads(obam, both, thr); fa = pysam.FastaFile(fasta)
    def table(lab, S):
        S = list(S); n = len(S); f = lambda g: sum(1 for x in S if g(x))
        ka = lambda x: all((x[0], x[1], i) in known for i in x[2]); ca = lambda x: all(canonical(fa, x[0], x[1], i) for i in x[2]); r1 = lambda x: n1[(x[0], x[2])]; r2 = lambda x: n2[(x[0], x[2])]
        ev = lambda x: r1(x) >= 1 or r2(x) >= 1 or x in o
        print(f'{lab}\t{n}\t{f(ka)} ({f(ka) / max(1, n):.3f})\t{f(ca)} ({f(ca) / max(1, n):.3f})\t{f(lambda x: r1(x) == 0)} / {f(lambda x: 1 <= r1(x) <= 2)} / {f(lambda x: r1(x) >= 3)}\t{f(lambda x: r2(x) >= 1)} ({f(lambda x: r2(x) >= 1) / max(1, n):.3f})\t{f(lambda x: x in o)} ({f(lambda x: x in o) / max(1, n):.3f})\t{f(ev)} ({f(ev) / max(1, n):.3f})\t{f(lambda x: not ev(x) and not ca(x))}')
    print('set\tchains\tevery intron annotated\tevery intron canonical\town full-length reads: 0 / 1-2 / 3 or more\tfull-length reads in the other individual\tassembled in the other individual today\tany of: own read, other individual\tno read anywhere and a non-canonical intron')
    table('extra transcripts of the arm', extra); table("today's unannotated transcripts", base)

if __name__ == '__main__':
    if sys.argv[1:] == ['--selftest']:
        import tempfile, pysam
        d = tempfile.mkdtemp(); p = os.path.join(d, 'g.fa'); open(p, 'w').write('>c\nAAAAGTCCCCAGAAAACTGGGGACAAAA\n'); pysam.faidx(p); fa = pysam.FastaFile(p)
        assert canonical(fa, 'c', '+', (4, 12)) and not canonical(fa, 'c', '-', (4, 12)) and canonical(fa, 'c', '-', (16, 24)); print('extra_transcripts selftest ok')
    else: main(*sys.argv[1:])
