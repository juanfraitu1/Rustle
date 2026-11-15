#!/usr/bin/env python3
"""PREREG row 1746: do reads that skip a short annotated exon carry the exon's bases at the junction?

  skip_realign.py REF_GTF BAM FASTA OUT_TSV [THREADS]   |   skip_realign.py --selftest
Unit = annotated chain x internal exon that no annotated chain skips (row 1744's skip table). For up to CAP full-chain
reads and CAP skip reads per unit, the read bases between W bp before the donor and W bp after the acceptor are compared
by edit distance with the spliced sequence with the exon (A) and without it (B)."""
import sys, collections, multiprocessing, random
from alt_events import load, skip_bin
W, CAP, MARGIN, MAXUNITS = 25, 20, 2, 1500

def edit(a, b):
    prev = list(range(len(b) + 1))
    for i, x in enumerate(a, 1):
        cur = [i]
        for j, y in enumerate(b, 1): cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (x != y)))
        prev = cur
    return prev[-1]

def qpos(start, cigar, x):
    """query index of the first read base at reference position >= x (bases inserted before it come after the index)."""
    r, q = start, 0
    for op, n in cigar:
        if op in (0, 7, 8):
            if r + n > x: return q + max(0, x - r)
            r += n; q += n
        elif op in (2, 3):
            if r + n > x: return q
            r += n
        elif op in (1, 4): q += n
    return q

def call(seg, a, b):
    da, db = edit(seg, a), edit(seg, b)
    return 'has_exon' if da + MARGIN <= db else 'skips' if db + MARGIN <= da else 'unclear'

def work(job):
    bam, fa_p, contig, units = job; import pysam   # units: [(chain, i, bin)]
    fa = pysam.FastaFile(fa_p); want = collections.defaultdict(list)
    for u, (ch, i, _) in enumerate(units):
        want[ch].append((u, 'full')); want[ch[:i] + ((ch[i][0], ch[i + 1][1]),) + ch[i + 2:]].append((u, 'skip'))
    seen = collections.Counter(); out = collections.Counter(); per = collections.defaultdict(collections.Counter)
    for r in pysam.AlignmentFile(bam).fetch(contig):
        if r.is_secondary or r.is_supplementary or 'N' not in r.cigarstring: continue
        pos = r.reference_start; ch = []
        for o, n in r.cigartuples:
            if o == 3: ch.append((pos, pos + n))
            if o in (0, 2, 3, 7, 8): pos += n
        for u, kind in want.get(tuple(ch), ()):
            if seen[u, kind] >= CAP: continue
            seen[u, kind] += 1; c, i, b = units[u]; don, acc = c[i][0], c[i + 1][1]
            lo, hi = max(r.reference_start, don - W), min(pos, acc + W)
            seg = r.query_sequence[qpos(r.reference_start, r.cigartuples, lo):qpos(r.reference_start, r.cigartuples, hi)].upper()
            up, dn = fa.fetch(contig, lo, don).upper(), fa.fetch(contig, acc, hi).upper(); ex = fa.fetch(contig, c[i][1], c[i + 1][0]).upper()
            v = call(seg, up + ex + dn, up + dn); out[b, kind, v] += 1; per[u][kind, v] += 1
    rows = [(contig, units[u][0][units[u][1]][1], units[u][0][units[u][1] + 1][0], units[u][2], *(p[k, v] for k in ('full', 'skip') for v in ('has_exon', 'skips', 'unclear'))) for u, p in per.items()]
    return out, rows

def selftest():
    assert edit('ACGT', 'ACGT') == 0 and edit('ACGT', 'AGT') == 1 and edit('', 'AC') == 2
    cig = [(0, 10), (3, 100), (0, 10)]                      # 10 M, 100 N, 10 M from reference 1000
    assert qpos(1000, cig, 1005) == 5 and qpos(1000, cig, 1010) == 10 and qpos(1000, cig, 1110) == 10 and qpos(1000, cig, 1115) == 15
    assert qpos(1000, [(0, 10), (1, 6), (3, 100), (0, 10)], 1110) == 16   # bases inserted at the junction lie inside the segment
    up, ex, dn = 'ACGTTGCAAC', 'GGATCC', 'TTGACCATGA'
    assert call(up + ex + dn, up + ex + dn, up + dn) == 'has_exon' and call(up + dn, up + ex + dn, up + dn) == 'skips'
    print('selftest ok')

if __name__ == '__main__':
    if sys.argv[1] == '--selftest': selftest(); sys.exit()
    ref_p, bam, fa, out_p = sys.argv[1:5]; thr = int(sys.argv[5]) if len(sys.argv) > 5 else 16
    import pysam
    ref = load(ref_p); have = set(pysam.AlignmentFile(bam).references); by = collections.defaultdict(list); nbin = collections.Counter()
    random.seed(1746); keys = sorted(k for k in ref if k[0] in have); random.shuffle(keys)   # shuffled so the unit cap is not one chromosome
    for c, s in keys:
        d = ref[c, s]
        for ch in sorted(d):
            for i in range(len(ch) - 1):
                L = ch[i + 1][0] - ch[i][1]; b = skip_bin(L)
                if b == '>60' or ch[:i] + ((ch[i][0], ch[i + 1][1]),) + ch[i + 2:] in d or nbin[b] >= MAXUNITS: continue
                nbin[b] += 1; by[c].append((ch, i, b))
    T = collections.Counter(); rows = []
    with multiprocessing.Pool(thr) as pool:
        for o, r in pool.imap_unordered(work, [(bam, fa, c, u) for c, u in by.items()]): T.update(o); rows += r
    with open(out_p, 'w') as o:
        o.write('contig\texon_start\texon_end\tbin\tfull_has\tfull_skips\tfull_unclear\tskip_has\tskip_skips\tskip_unclear\n')
        for r in sorted(rows): o.write('\t'.join(map(str, r)) + '\n')
    print('bin\tannotated units taken\tunits with skip reads\tfull-chain reads: has exon / skips / unclear (has share)\tskip reads: has exon / skips / unclear (has share)\tunits where >= half the skip reads have the exon')
    for b in ('<=15', '16-30', '31-60'):
        f = [T[b, 'full', v] for v in ('has_exon', 'skips', 'unclear')]; k = [T[b, 'skip', v] for v in ('has_exon', 'skips', 'unclear')]
        us = [r for r in rows if r[3] == b and sum(r[7:10])]
        print(f"{b}\t{nbin[b]}\t{len(us)}\t{'/'.join(map(str, f))} ({f[0] / max(1, sum(f)):.3f})\t{'/'.join(map(str, k))} ({k[0] / max(1, sum(k)):.3f})\t{sum(2 * r[7] >= sum(r[7:10]) for r in us)}")
