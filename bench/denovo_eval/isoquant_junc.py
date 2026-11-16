#!/usr/bin/env python3
"""PREREG r1765 (docs/PREREG_isoquant_junctions_2026-11-16.md): junction precision of ours vs IsoQuant genome-wide, per RefSeq gene
(gffcompare .tmap), at a matched junction count per gene, then slip classification of the wrong junctions.
  isoquant_junc.py REFSEQ.gtf BAM FASTA OURS.gtf OURS.tmap IQ.gtf IQ.tmap THREADS [WRONG.tsv]"""
import sys, re, collections, multiprocessing, pysam
ref, bam, fa, ours, otmap, iq, itmap, thr = sys.argv[1:9]
def introns(path):  # transcript_id -> (chrom, strand, [(donor0, acceptor0)])
    ex = collections.defaultdict(list); meta = {}
    for line in open(path):
        f = line.split('\t')
        if len(f) < 9 or f[2] != 'exon': continue
        t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1); ex[t].append((int(f[3]) - 1, int(f[4]))); meta[t] = (f[0], f[6])
    return {t: (*meta[t], [(a[1], b[0]) for a, b in zip(sorted(e), sorted(e)[1:])]) for t, e in ex.items() if len(e) > 1}
def by_gene(gtf, tmap):  # (chrom, gene) -> {junction: strand}
    tx = introns(gtf); out = collections.defaultdict(dict)
    for l in list(open(tmap))[1:]:
        f = l.split('\t')
        if f[0] != '-' and f[2] not in 'upxs' and f[4] in tx:
            c, s, js = tx[f[4]]
            for j in js: out[(c, f[0])][j] = s
    return out
def count(c):  # every intron of every non-supplementary record on one contig -> [primary, secondary]
    S = collections.defaultdict(lambda: [0, 0])
    for r in pysam.AlignmentFile(bam).fetch(c):
        if r.is_supplementary or r.is_unmapped: continue
        p = r.reference_start
        for op, n in r.cigartuples:
            if op == 3: S[(p, p + n)][r.is_secondary] += 1
            if op in (0, 2, 3, 7, 8): p += n
    return c, {j: v for j, v in S.items() if j in need[c]}
RC = str.maketrans('ACGT', 'TGCA')
def slip(j, J, S, tr):  # as row 1764
    best = None
    for k in J:
        if k != j and k in tr and abs(k[0] - j[0]) <= 10 and abs(k[1] - j[1]) <= 10 and (best is None or sum(S[k]) > sum(S[best])): best = k
    if best is None: return 'novel'
    return 'slip' if sum(S[j]) <= 0.2 * sum(S[best]) else 'alt_site'
if __name__ == '__main__':
    truth = collections.defaultdict(set)
    for c, s, js in introns(ref).values(): truth[c] |= set(js)
    O, I = by_gene(ours, otmap), by_gene(iq, itmap); both = sorted(g for g in O if g in I)
    need = collections.defaultdict(set)
    for g in both: need[g[0]] |= set(O[g]) | set(I[g])
    with multiprocessing.Pool(int(thr)) as P: sup = dict(P.imap_unordered(count, sorted(need)))  # fork: workers see `need`
    FA = pysam.FastaFile(fa); T = collections.defaultdict(collections.Counter); wrong = []
    for g in both:
        c = g[0]; S = collections.defaultdict(lambda: [0, 0], sup[c]); tr = truth[c]
        cut = set(sorted(O[g], key=lambda j: (-sum(S[j]), j))[:len(I[g])])  # shortcut: ties at the cut broken by coordinate
        for lab, J, st in (('ours', set(O[g]), O[g]), ('ours_cut', cut, O[g]), ('isoquant', set(I[g]), I[g])):
            for j in J:
                ok = j in tr; T[lab]['n'] += 1; T[lab]['ok'] += ok
                if not ok and lab != 'ours':
                    kind = slip(j, J, S, tr); m = (FA.fetch(c, j[0], j[0] + 2) + FA.fetch(c, j[1] - 2, j[1])).upper()
                    if st[j] == '-': m = m.translate(RC)[::-1]
                    canon = 'GTAG' if m == 'GTAG' else 'GCAG/ATAC' if m in ('GCAG', 'ATAC') else 'other'
                    n = sum(S[j]); T[lab][kind] += 1; T[lab]['m_' + canon] += 1; T[lab]['r_' + ('0' if n == 0 else '1-2' if n < 3 else '3+')] += 1
                    wrong.append((lab, g[1], c, j[0], j[1], S[j][0], S[j][1], m, kind))
    print(f"genes with a spliced model in both tools: {len(both)} (ours {len(O)}, isoquant {len(I)}); truth = RefSeq introns")
    print("set\tjunctions\tcorrect\tprecision\twrong\tslip\talt_site\tnovel\tGTAG\tGCAG/ATAC\tother\treads 0\t1-2\t3+")
    for lab in ('ours', 'ours_cut', 'isoquant'):
        c = T[lab]; print(lab, c['n'], c['ok'], f"{c['ok'] / max(1, c['n']):.3f}", c['n'] - c['ok'], c['slip'], c['alt_site'], c['novel'], c['m_GTAG'], c['m_GCAG/ATAC'], c['m_other'], c['r_0'], c['r_1-2'], c['r_3+'], sep='\t')
    d = T['ours_cut']['ok'] / max(1, T['ours_cut']['n']) - T['isoquant']['ok'] / max(1, T['isoquant']['n'])
    print(f"step B: ours at the cut minus isoquant = {d:+.3f} -> {'OPERATING POINT (bar -0.02 met)' if d >= -0.02 else 'bar missed, step C applies'}")
    if len(sys.argv) > 9: open(sys.argv[9], 'w').write(''.join('\t'.join(map(str, w)) + '\n' for w in sorted(wrong)))
