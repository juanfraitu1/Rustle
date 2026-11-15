#!/usr/bin/env python3
"""Row 1687, Rule S: cut transcripts fused across a tandem segmental duplication (docs/PREREG_sd_cut_rule_2026-11-13.md).

  sd_cut.py ASSEMBLY.gtf SEDEF.bed OUT.units.tsv      sd_cut.py --selftest
Output rows are a `--bridge-units-list` (tid chrom strand junctions label gene); no annotation is read.
"""
import sys, re, collections
BIN, MINBP = 100_000, 100


def tandem_pairs(lines):
    """direct, disjoint, same-contig SD pairs -> {(chrom, bin of the gap midpoint): [(l0, l1, r0, r1)]}"""
    ix = collections.defaultdict(list)
    for l in lines:
        f = l.split('\t')
        if len(f) < 10 or not f[1].isdigit() or f[0] != f[3] or f[8] != f[9]: continue
        a, b = sorted(((int(f[1]), int(f[2])), (int(f[4]), int(f[5]))))
        if a[1] <= b[0]: ix[(f[0], (a[1] + b[0]) // 2 // BIN)].append((a[0], a[1], b[0], b[1]))
    return ix


def inside(exons, s, e):
    return sum(max(0, min(b, e) - max(a, s)) for a, b in exons)


def cuts(exons, pairs):
    """exons: sorted 0-based half-open; pairs [(l0, l1, r0, r1)] -> sorted cut introns (exon end, next exon start)"""
    out = set()
    for l0, l1, r0, r1 in pairs:
        p = (l1 + r0) // 2
        for (a, b), (c, d) in zip(exons, exons[1:]):
            if b <= p < c:
                L, R = [x for x in exons if x[1] <= p], [x for x in exons if x[0] > p]
                il, ir = inside(L, l0, l1), inside(R, r0, r1)
                if il >= MINBP and ir >= MINBP and 2 * il >= sum(y - x for x, y in L) and 2 * ir >= sum(y - x for x, y in R): out.add((b, c))
    return sorted(out)


def selftest():
    ex = [(100, 300), (400, 600), (1100, 1300), (1400, 1600)]
    assert cuts(ex, [(50, 700, 1000, 1700)]) == [(600, 1100)]
    assert cuts(ex, [(50, 1150, 1250, 1700)]) == []                    # p inside an exon: no intron contains it
    assert cuts(ex, [(350, 700, 1000, 1700)]) == [(600, 1100)] and cuts(ex, [(550, 700, 1000, 1700)]) == []   # 50 bp in the left mate
    assert cuts([(0, 3000), (5100, 5300), (5400, 5600), (6100, 6300), (6400, 6600)], [(5050, 5700, 6000, 6700)]) == []   # most of the left side lies outside the mate
    row = lambda c2, s1, s2: '\t'.join(['c', '50', '700', c2, '1000', '1700', 'S', '1', s1, s2])
    ix = tandem_pairs([row('c', '+', '+'), row('c', '+', '-'), row('d', '+', '+')]); assert sum(map(len, ix.values())) == 1 and ix[('c', 0)] == [(50, 700, 1000, 1700)]
    print('selftest ok')


def main(gtf, bed, out):
    ix = tandem_pairs(open(bed)); A = lambda k, s: (re.search(k + r' "([^"]+)"', s) or [None, None])[1]; tx = {}
    for l in open(gtf):
        f = l.rstrip('\n').split('\t')
        if len(f) >= 9 and f[2] == 'exon':
            tx.setdefault(A('transcript_id', f[8]), [f[0], f[6], A('gene_id', f[8]), []])[3].append((int(f[3]) - 1, int(f[4])))
    n = 0; genes = set()
    with open(out, 'w') as w:
        w.write('tid\tchrom\tstrand\tjunctions\tlabel\tgene\n')
        for t, (c, st, g, ex) in tx.items():
            ex.sort()
            if len(ex) < 2: continue
            ps = [p for b in range(ex[0][0] // BIN, ex[-1][1] // BIN + 1) for p in ix.get((c, b), ())]
            cs = cuts(ex, ps) if ps else []
            if cs:
                n += 1; genes.add(g); w.write(f"{t}\t{c}\t{st}\t{';'.join(f'{b + 1}-{d}' for b, d in cs)}\tsdcut\t{g}\n")   # exon end + 1 .. next exon start - 1, 1-based
    print(f'transcripts cut {n}; assembly loci with a cut {len(genes)}')


if __name__ == '__main__':
    selftest() if sys.argv[1:] == ['--selftest'] else main(*sys.argv[1:4])
