#!/usr/bin/env python3
"""Why annotated isoforms with reads of their own are not in the assembly (register 1966).

    thin_reasons.py ANNOTATION.gff BAM ASSEMBLY.gtf POLISH_REASONS.tsv [THREADS]

Own reads of an isoform = primary alignments, mapping quality >= 1, whose intron chain is exactly the annotated chain.
For each annotated isoform with own reads that ASSEMBLY.gtf lacks: the rule of POLISH_REASONS.tsv that removed a candidate
with exactly that chain, else `no candidate`. Also the DISTINCT (start, end) records among the own reads (the read floor
counts distinct coordinate records)."""
import sys, os, collections, multiprocessing
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import thin_isoforms as T
CL = ('1', '2', '3', '4', '5-9', '10 or more')
def kl(n): return str(n) if n < 5 else CL[4] if n < 10 else CL[5]

def scan(contig):
    import pysam
    want = T.G['by'].get(contig)
    if not want: return {}
    out = collections.defaultdict(set); n = collections.Counter()
    for r in pysam.AlignmentFile(T.G['bam']).fetch(contig):
        if r.is_secondary or r.is_supplementary or r.mapping_quality < 1 or 'N' not in r.cigarstring: continue
        c = []; p = r.reference_start
        for o, k in r.cigartuples:
            if o == 3: c.append((p, p + k))
            if o in (0, 2, 3, 7, 8): p += k
        c = tuple(c)
        if c in want: n[c] += 1; out[c].add((r.reference_start, p))
    return {(contig, c): (v, len(out[c])) for c, v in n.items()}

def main(gff, bam, gtf, reasons, thr=32):
    import pysam
    ann = T.chains_gff(gff); by = collections.defaultdict(set)
    for c, s, ch in ann: by[c].add(ch)
    T.G.update(bam=bam, by=by)
    with multiprocessing.get_context('fork').Pool(int(thr)) as pool: res = pool.map(scan, list(pysam.AlignmentFile(bam).references), chunksize=1)
    own = {}
    for d in res: own.update(d)
    asm = T.chains_gtf(gtf); why = {}
    for l in list(open(reasons))[1:]:
        t, c, s, ex, n, rule = l.rstrip('\n').split('\t')[:6]; e = [tuple(map(int, x.split('-'))) for x in ex.split(',')]
        if len(e) > 1: why.setdefault((c, s, tuple((e[i][1], e[i + 1][0] - 1) for i in range(len(e) - 1))), rule)
    C = collections.Counter(); D = collections.Counter(); tot = collections.Counter(); got = collections.Counter()
    for a in ann:
        k = (a[0], a[2])
        if k not in own: continue
        n, dist = own[k]; c = kl(n); tot[c] += 1
        if a in asm: got[c] += 1; continue
        r = why.get(a, 'no candidate'); C[(c, r)] += 1
        if r == 'no candidate': D[(c, 'one distinct record' if dist == 1 else 'two or more distinct records')] += 1
    rules = sorted({r for _, r in C}, key=lambda r: -sum(v for (c, x), v in C.items() if x == r))
    print('own reads\tannotated isoforms\tassembled\tnot assembled\t' + '\t'.join(rules) + '\tno candidate: one distinct record\tno candidate: two or more')
    for c in CL: print(c, tot[c], got[c], tot[c] - got[c], *[C[(c, r)] for r in rules], D[(c, 'one distinct record')], D[(c, 'two or more distinct records')], sep='\t')

if __name__ == '__main__':
    if sys.argv[1:] == ['--selftest']: assert kl(1) == '1' and kl(4) == '4' and kl(7) == '5-9' and kl(10) == '10 or more'; print('thin_reasons selftest ok')
    else: main(*sys.argv[1:])
