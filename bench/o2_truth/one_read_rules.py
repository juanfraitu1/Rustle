#!/usr/bin/env python3
"""One-read chains by the reads they EXPLAIN (register 1990).

    one_read_rules.py BASE.gtf FLOOR1.gtf BAM GENOME.fa TRUTH.gtf OUT.tsv [THREADS]

Units: multi-exon chains of FLOOR1 that BASE lacks and exactly one primary read has. Columns: see the header line."""
import sys, os, re, collections, multiprocessing
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../tools'))
import oneread, gtf_ptc, multicopy_rescue as MR

def relations(C, others, base):
    """others = {read chain: reads} sharing an intron with C; base = assembled chains sharing an intron with C.
    -> (sub, sub_only, longer, conflict, min reads on an intron of C incl. its own read, sub reads holding C's last / first intron)"""
    sub = only = longer = conflict = last = first = 0
    for r, n in others.items():
        if r == C: continue
        if oneread.run(r, C):
            sub += n; only += n * (not any(oneread.run(r, b) for b in base)); last += n * (C[-1] in r); first += n * (C[0] in r)
        elif oneread.run(C, r): longer += n
        else: conflict += n
    return sub, only, longer, conflict, min(sum(n for r, n in others.items() if i in r) for i in C), last, first

def main(base, f1, bam, fasta, truth, out, thr=32):
    import pysam
    B = oneread.models(base); have = collections.defaultdict(set); bhas = collections.defaultdict(set)
    for c, s, ex, _ in B.values():
        if len(ex) > 1:
            ch = oneread.chain(ex); have[c].add(ch)
            for i in ch: bhas[(c, i)].add(ch)
    T = oneread.models(truth); tr = set(); tby = collections.defaultdict(list); tin = set()
    for c, s, ex, _ in T.values():
        if len(ex) > 1:
            ch = oneread.chain(ex); tr.add((c, ch))
            for i in ch: tby[(c, i)].append(ch); tin.add((c, i))
    cand = collections.OrderedDict()
    for t, (c, s, ex, lines) in oneread.models(f1).items():
        k = (c, oneread.chain(ex))
        if len(ex) > 1 and s in '+-' and k[1] not in have[c] and k not in cand: cand[k] = (t, s, ex)
    want = {(c, i) for c, ch in cand for i in ch}; prim = collections.Counter()
    contigs = sorted({c for c, _ in want} & set(pysam.AlignmentFile(bam).references))
    oneread._SCAN.update(bam=bam, want=want)
    with multiprocessing.get_context('fork').Pool(int(thr)) as pool:
        for part in pool.imap_unordered(oneread._scan, contigs, chunksize=1): prim.update(part)
    has = collections.defaultdict(set)
    for c, ch in prim:
        for i in ch: has[(c, i)].add(ch)
    sp = MR.spans(base); fa = pysam.FastaFile(fasta)
    with open(out, 'w') as o:
        o.write('tid\texact\tpart\tintrons_annotated\tintrons\tsub\tsub_only\tlonger\tconflict\tmin_intron_reads\tsub_last\tsub_first\tbase\torf\tptc\n')
        for (c, C), (t, s, ex) in cand.items():
            if prim.get((c, C), 0) != 1: continue
            oth = {r: prim[(c, r)] for i in C for r in has[(c, i)]}; bs = {b for i in C for b in bhas.get((c, i), ())}
            rel = relations(C, oth, bs)
            if s == '-': rel = rel[:5] + (rel[6], rel[5])        # last = the 3' most intron of the transcript
            b = 'part' if any(oneread.run(C, x) for x in bs) else 'extends' if any(oneread.run(x, C) for x in bs) else 'conflict' if bs else \
                'near' if MR.occupied(sp, c, s, ex[0][0], ex[-1][1]) else 'empty'
            orf, ptc = gtf_ptc.call(ex, s, lambda a, e: fa.fetch(c, a, e)); tc = tby.get((c, C[0]), ())
            o.write('\t'.join(map(str, (t, int((c, C) in tr), int((c, C) not in tr and any(oneread.run(C, x) for x in tc)), int(all((c, i) in tin for i in C)), len(C)) + rel + (b, orf, int(bool(ptc))))) + '\n')

def selftest():
    C = ((10, 20), (30, 40), (50, 60))
    oth = {C: 1, ((30, 40), (50, 60)): 3, ((10, 20), (30, 40)): 1, ((5, 8),) + C: 2, ((10, 20), (30, 45)): 4}
    assert relations(C, oth, {((30, 40), (50, 60), (70, 80))}) == (4, 1, 2, 4, 6, 3, 1), relations(C, oth, {((30, 40), (50, 60), (70, 80))})
    print('one_read_rules selftest ok')

if __name__ == '__main__': selftest() if sys.argv[1:] == ['--selftest'] else main(*sys.argv[1:])
