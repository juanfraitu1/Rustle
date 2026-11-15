#!/usr/bin/env python3
"""PREREG row 1754: find the exon of <= 15 bp a read skipped by searching the intron for the extra bases the read carries.

  skip_search.py search REF_GTF BAM FASTA OUT_PREFIX THREADS R1746_TSV R1751_UNITS   |   skip_search.py --selftest
A read's intron (d, a) with L = 3..15 extra read bases in the junction window (row 1751's window) is searched for exons of
L bp with the intron's own splice signals around them and >= MININ bp of intron on each side. The search is used only when
(number of such sites) x 4^-L <= EMAX, and a site is taken only when it is the single one whose sequence is in the read and
leaves <= 1 edit and none of the extra bases unexplained. Decoy = the same search on the reversed intron.
OUT.units.tsv (these units + row 1751's) is the input of `skip_fix.py fix` and `skip_fix.py score`."""
import sys, collections, multiprocessing, functools
from skip_realign import edit, W
from skip_fix import introns, expand, cut, primary_spliced, read_units, annot, write_units, MAXEX
EMAX, MININ = 0.01, 50
SIG = {('GT', 'AG'): ('AG', 'GT'), ('CT', 'AC'): ('AC', 'CT')}   # intron ends -> bases before / after an exon inside it

def sites(intr, L, pre, post):
    """offsets of the L bp windows of the intron with the splice signals around them, >= MININ bp from each end"""
    out, p = [], intr.find(pre, MININ - 2)
    while p != -1 and p + 2 + L <= len(intr) - MININ:
        if intr.startswith(post, p + 2 + L): out.append(p + 2)
        p = intr.find(pre, p + 1)
    return out

def search(seg, up, dn, intr, L, pos):
    """offset of the one site that explains the read's extra bases; 'weak' when the intron has too many sites for L; else None"""
    if len(pos) * 4.0 ** -L > EMAX: return 'weak'
    lim = min(1, edit(seg, up + dn) - L)
    hits = [q for q in pos if intr[q:q + L] in seg and edit(seg, up + intr[q:q + L] + dn) <= lim]
    return hits[0] if len(hits) == 1 else None

def scan(job):
    bam, fa_p, contig = job; import pysam
    fa = pysam.FastaFile(fa_p); fetch = functools.lru_cache(20000)(lambda s, e: fa.fetch(contig, s, e).upper())
    site = functools.lru_cache(50000)(lambda d, a, L, rev: sites(fetch(d, a)[::-1] if rev else fetch(d, a), L, *SIG[fetch(d, d + 2), fetch(a - 2, a)]))
    n, pair, found, st = collections.Counter(), collections.Counter(), collections.Counter(), collections.Counter()
    for r in pysam.AlignmentFile(bam).fetch(contig):
        if not primary_spliced(r): continue
        cig, start = r.cigartuples, r.reference_start; ivs, end = introns(start, cig); n.update(ivs)
        for (d, x), (y, a) in zip(ivs, ivs[1:]):
            if 0 < y - x <= MAXEX: pair[d, x, y, a] += 1
        ks = [i for i, (o, _) in enumerate(cig) if o == 3]; seq = r.query_sequence; ops = None
        for k, ci in enumerate(ks):
            if not seq or not any(o == 1 and m >= 2 for o, m in cig[max(0, ci - 3):ci + 4]): continue   # no inserted bases near
            d, a = ivs[k]; ops = ops or expand(cig)
            lo = max(ivs[k - 1][1] if k else start, d - W); hi = min(ivs[k + 1][0] if k + 1 < len(ivs) else end, a + W)
            c = cut(ops, start, lo, hi)
            if not c: continue
            seg = seq[c[1]:c[3]].upper(); up, dn = fetch(lo, d), fetch(a, hi); L = len(seg) - len(up) - len(dn)
            if not 3 <= L <= MAXEX or a - d < 2 * MININ + L or (fetch(d, d + 2), fetch(a - 2, a)) not in SIG: continue
            intr = fetch(d, a); st[L, 'tested'] += 1
            q = search(seg, up, dn, intr, L, site(d, a, L, False))
            if q == 'weak': st[L, 'weak'] += 1; continue
            if q is not None: st[L, 'found'] += 1; found[d, d + q, d + q + L, a] += 1
            st[L, 'decoy'] += search(seg, up, dn, intr[::-1], L, site(d, a, L, True)) not in (None, 'weak')
    keys = set(found) | {k for k in pair if n[k[0], k[3]]}
    return [(contig, *k, pair[k], n[k[0], k[3]], found[k]) for k in keys], st

def selftest():
    import random; rng = random.Random(7); body = ''.join(rng.choice('ACGT') for _ in range(400)).replace('AG', 'AC')
    ex = 'TTGACCATGCA'; intr = 'GT' + body[:150] + 'AG' + ex + 'GT' + body[150:] + 'AG'; q = 2 + 150 + 2
    assert q in sites(intr, 11, 'AG', 'GT') and sites('GT' + 'C' * 30 + 'AG' + ex + 'GT' + 'C' * 30 + 'AG', 11, 'AG', 'GT') == []   # too near the ends
    up, dn = 'ACGTTGCAACGGATCCAAGT', 'CCATGATTGACAGGATACCA'; pos = sites(intr, 11, 'AG', 'GT')
    assert search(up + ex + dn, up, dn, intr, 11, pos) == q
    assert search(up + ex[::-1] + dn, up, dn, intr, 11, pos) is None            # other bases at the junction
    assert search(up + ex + dn, up, dn, intr, 11, pos + [q]) is None            # two sites fit: no call
    assert search(up + 'TTG' + dn, up, dn, intr, 3, list(range(700))) == 'weak'  # 700 sites x 4^-3: chance
    print('selftest ok')

if __name__ == '__main__':
    if sys.argv[1] == '--selftest': selftest(); sys.exit()
    import pysam
    ref_p, bam, fa, out, thr, prior, junc = sys.argv[2:9]; klass = annot(ref_p); rows, st = {}, collections.Counter()
    with multiprocessing.Pool(int(thr)) as pool:
        for r, s in pool.imap_unordered(scan, [(bam, fa, c) for c in pysam.AlignmentFile(bam).references]):
            st.update(s); rows.update({x[:5]: x[5:] for x in r})
    jfired = {u[:5] for u in read_units(junc) if u[8]}; sfound = {k for k, v in rows.items() if v[2]}
    # one table for the rewrite: reads_fixed = reads the search placed, or row 1751's count for a unit only it fired on
    old = {u[:5]: u[8] for u in read_units(junc)}
    write_units(sorted((*k, klass(*k), v[0], v[1], v[2] or old.get(k, 0)) for k, v in rows.items()), out)
    print('exon length\treads with that many extra bases at a junction\tintron too large for the length\texon found\tfound in the reversed intron (decoy)')
    for L in range(3, MAXEX + 1): print(f'{L}\t' + '\t'.join(str(st[L, k]) for k in ('tested', 'weak', 'found', 'decoy')))
    f, d = (sum(v for (L, k), v in st.items() if k == x) for x in ('found', 'decoy'))
    print(f'gate: decoy / found reads\t{d}\t{f}\t{d / max(1, f):.4f}')
    print('annotation of the unit\tunits found by the search\tof them not fired by row 1751\treads placed')
    for k in ('incl', 'both', 'skip', 'neither', 'all'):
        u = [x for x in sfound if k in (klass(*x), 'all')]; print(f'{k}\t{len(u)}\t{len(set(u) - jfired)}\t{sum(rows[x][2] for x in u)}')
    short = collections.defaultdict(set)                    # annotated intron pair around an exon of <= 15 bp, by its outer ends
    from alt_events import load
    for (c, s), dd in load(ref_p).items():
        for ch in dd:
            for (d0, x), (y, a) in zip(ch, ch[1:]):
                if 0 < y - x <= MAXEX: short[c, d0, a].add((x, y))
    judged = [x for x in sfound if (x[0], x[1], x[4]) in short]; ok = sum((x[2], x[3]) in short[x[0], x[1], x[4]] for x in judged)
    print(f'precision: found units whose outer junctions flank an annotated exon of <= 15 bp / of them the annotated exon\t{len(judged)}\t{ok}\t{ok / max(1, len(judged)):.3f}')
    want = {(r[0], int(r[1]), int(r[2])) for r in (l.split('\t') for l in list(open(prior))[1:]) if r[3] == '<=15' and sum(map(int, r[7:10])) and 2 * int(r[7]) >= sum(map(int, r[7:10]))}
    s_, j_ = {(x[0], x[2], x[3]) for x in sfound}, {(x[0], x[2], x[3]) for x in jfired}
    print(f'row 1746 exons\t{len(want)}\tsearch\t{len(want & s_)}\trow 1751 rule\t{len(want & j_)}\teither\t{len(want & (s_ | j_))}\treach (either)\t{len(want & (s_ | j_)) / max(1, len(want)):.3f}')
    print('row 1746 exons not reached, by length\t' + ' '.join(f'{L}={n}' for L, n in sorted(collections.Counter(y - x for c, x, y in want - s_ - j_).items())))
