#!/usr/bin/env python3
"""PREREG r1677: read-derived features of every candidate multi-exon chain (the 1-read-floor superset), RefSeq label.

  chain_rules_feat.py REF_GTF DEFAULT_GTF FLOOR1_GTF BAM FASTA OUT_TSV [THREADS] [CONTIG,...]

Reads = primary spliced records. Chains are keyed (contig, intron chain); strand comes from the GTF model.
Features r1619 lacked: full / partial / container read counts, intron support over all reads, read-end
concentration, polyA signal, tied-read share, read divergence, chains per locus."""
import sys, re, collections, multiprocessing, pysam

def load(path):
    """(contig, chain) -> (reads, strand, gene, exons) of the best-supported model; chain = 0-based half-open introns."""
    ex = collections.defaultdict(list); at = {}
    for l in open(path):
        f = l.split('\t')
        if len(f) < 9 or f[2] not in ('exon', 'transcript'): continue
        k = (f[0], f[6], re.search(r'transcript_id "([^"]+)"', f[8]).group(1))
        if f[2] == 'exon': ex[k].append((int(f[3]) - 1, int(f[4])))
        else: at[k] = (float(re.search(r'reads "([^"]+)"', f[8]).group(1)), re.search(r'gene_id "([^"]+)"', f[8]).group(1))
    T = {}
    for k, e in ex.items():
        if len(e) < 2: continue
        e.sort(); c = (k[0], tuple((e[i][1], e[i + 1][0]) for i in range(len(e) - 1))); r, g = at[k]
        if c not in T or r > T[c][0]: T[c] = (r, k[1], g, e)
    return T

def ref_chains(path):
    import gzip
    ex = collections.defaultdict(list)
    for l in (gzip.open(path, 'rt') if path.endswith('.gz') else open(path)):
        f = l.split('\t')
        if len(f) > 8 and f[2] == 'exon': ex[(f[0], re.search(r'transcript_id "([^"]+)"', f[8]).group(1))].append((int(f[3]) - 1, int(f[4])))
    out = set()
    for k, e in ex.items():
        e.sort()
        if len(e) > 1: out.add((k[0], tuple((e[i][1], e[i + 1][0]) for i in range(len(e) - 1))))
    return out

def is_run(a, b, pos):  # a is the contiguous run of b starting at pos
    return b[pos:pos + len(a)] == a

def work(job):
    bam, fa_p, contig, cand = job  # cand: chain -> (reads, strand, gene, exons, in_ref, default_reads)
    cnt = collections.Counter(); isup = collections.Counter()
    lo = collections.defaultdict(collections.Counter); hi = collections.defaultdict(collections.Counter)
    mq0 = collections.Counter(); de = collections.Counter()
    for r in pysam.AlignmentFile(bam).fetch(contig):
        if r.is_secondary or r.is_supplementary or 'N' not in r.cigarstring: continue
        pos = r.reference_start; ch = []
        for op, ln in r.cigartuples:
            if op == 3: ch.append((pos, pos + ln))
            if op in (0, 2, 3, 7, 8): pos += ln
        ch = tuple(ch); cnt[ch] += 1
        for i in ch: isup[i] += 1
        if ch in cand:
            lo[ch][r.reference_start] += 1; hi[ch][pos] += 1; mq0[ch] += r.mapping_quality < 5
            de[ch] += r.get_tag('de') if r.has_tag('de') else 0.0
    by_first = collections.defaultdict(list)  # intron -> candidates containing it, with its position
    for c in cand:
        for p, i in enumerate(c): by_first[i].append((c, p))
    part = collections.Counter(); cl = collections.Counter(); cr = collections.Counter()
    for R, n in cnt.items():
        for c, p in by_first.get(R[0], ()):  # R a strict contiguous run of candidate c -> partial read of c
            if len(R) < len(c) and is_run(R, c, p): part[c] += n
        for p, i in enumerate(R):  # candidate c a strict contiguous run of R -> R is a container of c
            for c, q in by_first.get(i, ()):
                if q == 0 and len(c) < len(R) and is_run(c, R, p):
                    if p > 0: cl[c] += n
                    if p + len(c) < len(R): cr[c] += n
    gtop = collections.Counter(); gn = collections.Counter(); gir = collections.Counter()
    for c, v in cand.items():
        gn[v[2]] += 1; gtop[v[2]] = max(gtop[v[2]], v[0]); gir[v[2]] = max(gir[v[2]], max(isup[i] for i in c))
    starts = sorted(isup)
    import bisect
    fa = pysam.FastaFile(fa_p); L = fa.get_reference_length(contig); comp = str.maketrans('ACGT', 'TGCA')
    rows = []
    for c, (reads, strand, gene, e, in_ref, dreads) in cand.items():
        full = cnt[c]; mi = min(isup[i] for i in c); wob = 0; ret = 0
        for i in c:
            for j in starts[bisect.bisect_left(starts, (i[0] - 10, 0)):bisect.bisect_right(starts, (i[0] + 10, 1 << 60))]:
                if j != i and abs(j[1] - i[1]) <= 10 and isup[j] >= 5 * max(1, isup[i]): wob = 1
        for a, b in e:
            for j in starts[bisect.bisect_left(starts, (a, 0)):bisect.bisect_left(starts, (b, 0))]:
                if j[1] <= b and isup[j] >= 5 * max(1.0, reads): ret = 1
        ends3 = (hi if strand == '+' else lo)[c]; ends5 = (lo if strand == '+' else hi)[c]
        def conc(E):  # share of full-chain reads within 10 bp of the commonest end, and that end
            if not E: return 0.0, -1
            m = max(E, key=E.get); return sum(n for p, n in E.items() if abs(p - m) <= 10) / sum(E.values()), m
        c3, m3 = conc(ends3); c5, m5 = conc(ends5)
        pas = arich = 0
        if m3 >= 0:
            if strand == '+': up = fa.fetch(contig, max(0, m3 - 40), m3).upper(); dn = fa.fetch(contig, m3, min(L, m3 + 20)).upper()
            else: up = fa.fetch(contig, m3, min(L, m3 + 40)).upper().translate(comp)[::-1]; dn = fa.fetch(contig, max(0, m3 - 20), m3).upper().translate(comp)[::-1]
            pas = int('AATAAA' in up[:-5] or 'ATTAAA' in up[:-5]); arich = dn.count('A')
        c5n, c3n = (cl[c], cr[c]) if strand == '+' else (cr[c], cl[c])
        rows.append((contig, strand, in_ref, int(dreads > 0), reads, dreads, len(c), sum(b - a for a, b in e), full, part[c], c5n, c3n,
                     mi, round(mi / max(1, gir[gene]), 4), round(reads / gtop[gene], 4), gn[gene], wob, ret,
                     round(c3, 3), round(c5, 3), pas, arich, round(mq0[c] / full, 3) if full else 0, round(de[c] / full, 5) if full else 0,
                     ','.join(f'{x}-{y}' for x, y in c)))
    return rows

COLS = 'chrom strand in_ref in_default reads default_reads n_introns length full partial cont5 cont3 min_isup min_isup_rel locus_frac locus_chains wobble retention end3_conc end5_conc pas a_rich mq0_share de key'.split()  # key = the intron chain, 0-based half-open (joins other tables)
if __name__ == '__main__':
    ref_p, def_p, fl_p, bam, fa_p, out = sys.argv[1:7]
    thr = int(sys.argv[7]) if len(sys.argv) > 7 else 16; only = set(sys.argv[8].split(',')) if len(sys.argv) > 8 else None
    ref = ref_chains(ref_p); D = load(def_p); F = load(fl_p)
    for k, v in D.items(): F.setdefault(k, v)  # superset = floor-1 chains plus any default chain the floor-1 run dropped
    per = collections.defaultdict(dict)
    for (ctg, c), (r, s, g, e) in F.items():
        if only and ctg not in only: continue
        per[ctg][c] = (r, s, g, e, int((ctg, c) in ref), D[(ctg, c)][0] if (ctg, c) in D else 0)
    jobs = [(bam, fa_p, ctg, cand) for ctg, cand in sorted(per.items(), key=lambda x: -len(x[1]))]
    with multiprocessing.Pool(thr) as pool, open(out, 'w') as w:
        w.write('\t'.join(COLS) + '\n')
        for rows in pool.imap_unordered(work, jobs):
            for r in rows: w.write('\t'.join(f'{x:g}' if isinstance(x, float) else str(x) for x in r) + '\n')
