#!/usr/bin/env python3
"""Member discovery: the places every expressed locus's reads also align to (register 1910; O2, the `available` layer).

    member_discovery.py ASSEMBLY.gtf BAM OUT_PREFIX [--min-reads 5] [--min-bases 100] [--min-ratio 0.3] [--merge 1000]
                        [--all-hits OUT.hits.tsv] [--threads N]

Loci = gene_id of the GTF (rule of eq_classes.py). A read's SOURCE = the locus of its best-AS record. Every other
non-supplementary record of the read that is not on the source locus is a hit; hits on one strand are merged into
places (overlapping or within --merge). A place is a member candidate of source S when >= --min-reads reads of S hit it
with >= --min-bases aligned bases and AS >= --min-ratio x the read's best AS.
Writes OUT.members.tsv (place, chrom, start, end, strand, source, reads, best_ratio, median_ratio, max_bases, locus = the
assembled locus the place coincides with, else `.`) and OUT.members.gtf (one model per place that is no assembled locus:
exons of its best record; o2_status "member_candidate"). It says where relatives are, nothing about their expression.
Defaults (5 reads, ratio 0.3) are row 1910's: the loosest setting at which <= 10% of control genes without relatives got a
place (draw A), confirmed on draw B: 6% of controls, 90% of places on listed relatives or annotated pseudogenes, 72% of
recent paralogues and 48% of named processed pseudogenes reached.
--all-hits also writes every hit (read, source, chrom, start, end, strand, bases, ratio, locus) for threshold studies.
Memory: one entry per record while scanning, then one per read (source and best AS) and one per hit. --selftest runs the pure parts."""
import sys, os, re, argparse, collections, statistics, multiprocessing
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
OUTSIDE = 'outside'

def places(hits, merge=1000):
    """hits [(chrom, strand, start, end, ...)] -> [(chrom, strand, start, end, [hit indices])], merged per chrom + strand"""
    out = []; by = collections.defaultdict(list)
    for i, h in enumerate(hits): by[(h[0], h[1])].append((h[2], h[3], i))
    for (c, st), v in sorted(by.items()):
        v.sort(); cs, ce, idx = v[0][0], v[0][1], [v[0][2]]
        for s, e, i in v[1:]:
            if s <= ce + merge: ce = max(ce, e); idx.append(i)
            else: out.append((c, st, cs, ce, idx)); cs, ce, idx = s, e, [i]
        out.append((c, st, cs, ce, idx))
    return out

def candidates(hits, min_reads=2, min_bases=100, min_ratio=0.0, merge=1000):
    """hits [(chrom, strand, start, end, read, source, bases, ratio, locus)] -> rows (chrom, start, end, strand, source, reads, best, median, max_bases, locus)"""
    hs = [h for h in hits if h[6] >= min_bases and h[7] >= min_ratio]; rows = []
    for c, st, s, e, idx in places(hs, merge):
        by = collections.defaultdict(dict); loc = collections.Counter(); mb = collections.Counter()
        for i in idx:
            h = hs[i]; by[h[5]][h[4]] = max(by[h[5]].get(h[4], 0), h[7]); mb[h[5]] = max(mb[h[5]], h[6])
            if h[8] != OUTSIDE: loc[h[8]] += 1
        l = loc.most_common(1)[0][0] if loc else '.'
        for src, rd in sorted(by.items()):
            if len(rd) >= min_reads: rows.append((c, s, e, st, src, len(rd), max(rd.values()), statistics.median(rd.values()), mb[src], l))
    return rows

G = {}
def scan(c):
    import pysam
    from gtf_names import overlap_genes
    import eq_classes as E
    gex = G['gex']; out = []
    idx = {s: ([a for a, _, _ in gex[(c, s)]], max(b - a for a, b, _ in gex[(c, s)])) for s in '+-' if (c, s) in gex}
    for r in pysam.AlignmentFile(G['bam']).fetch(c):
        if r.is_unmapped or r.is_supplementary or not r.has_tag('AS'): continue
        s = E.tx_strand(r); bl = r.get_blocks()
        ov = overlap_genes([(a + 1, b) for a, b in bl], gex[(c, s)], *idx[s]) if s in idx else None
        g = max(ov, key=lambda k: (ov[k], k)) if ov else OUTSIDE
        out.append((r.query_name, r.get_tag('AS'), g, s, r.reference_start, r.reference_end, sum(b - a for a, b in bl)))
    return out

def selftest():
    h = [('c', '+', 100, 400, 'r1', 'A', 300, 0.9, OUTSIDE), ('c', '+', 150, 500, 'r2', 'A', 350, 0.6, OUTSIDE), ('c', '+', 5000, 5300, 'r3', 'A', 300, 0.9, OUTSIDE),
         ('c', '-', 120, 380, 'r4', 'A', 260, 0.9, OUTSIDE), ('c', '+', 1400, 1700, 'r5', 'B', 300, 0.4, 'L9'), ('c', '+', 1450, 1800, 'r6', 'B', 80, 0.9, 'L9'), ('c', '+', 1500, 1900, 'r7', 'B', 400, 0.5, 'L9')]
    p = places(h); assert [(x[1], x[2], x[3], len(x[4])) for x in p] == [('+', 100, 1900, 5), ('+', 5000, 5300, 1), ('-', 120, 380, 1)], p
    r = candidates(h); assert r == [('c', 100, 1900, '+', 'A', 2, 0.9, 0.75, 350, 'L9'), ('c', 100, 1900, '+', 'B', 2, 0.5, 0.45, 400, 'L9')], r
    assert candidates(h, min_ratio=0.5) == [('c', 100, 1900, '+', 'A', 2, 0.9, 0.75, 350, 'L9')] and candidates(h, min_reads=3) == []
    assert candidates(h, min_ratio=0.5, merge=0) == [('c', 100, 500, '+', 'A', 2, 0.9, 0.75, 350, '.')]
    print('member_discovery selftest ok')

def main():
    if '--selftest' in sys.argv: return selftest()
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('gtf'); ap.add_argument('bam'); ap.add_argument('out'); ap.add_argument('--min-reads', type=int, default=5); ap.add_argument('--min-bases', type=int, default=100)
    ap.add_argument('--min-ratio', type=float, default=0.3); ap.add_argument('--merge', type=int, default=1000); ap.add_argument('--all-hits'); ap.add_argument('--threads', type=int, default=8)
    a = ap.parse_args()
    import pysam, eq_classes as E
    meta, gex = E.load_loci(a.gtf); G.update(gex=gex, bam=a.bam)
    with multiprocessing.get_context('fork').Pool(a.threads) as pool:
        res = pool.map(scan, list(pysam.AlignmentFile(a.bam).references), chunksize=1)
    best = {}                                     # read -> (best AS, its locus)
    for part in res:
        for i, (q, s, g, *_) in enumerate(part):
            b = best.get(q)
            if b is None or s > b[0] or (s == b[0] and g != OUTSIDE and b[1] == OUTSIDE): best[q] = (s, g)
    hits = []
    for ci, part in enumerate(res):
        c = None
        for q, s, g, st, a0, a1, nb in part:
            b, src = best[q]
            if src == OUTSIDE or g == src or b <= 0: continue
            hits.append((G_REF[ci], st, a0, a1, q, src, nb, s / b, g))
    del res
    if a.all_hits:
        with open(a.all_hits, 'w') as w:
            w.write('read\tsource\tchrom\tstart\tend\tstrand\tbases\tratio\tlocus\n')
            for h in hits: w.write(f'{h[4]}\t{h[5]}\t{h[0]}\t{h[2]}\t{h[3]}\t{h[1]}\t{h[6]}\t{h[7]:.3f}\t{h[8]}\n')
    rows = candidates(hits, a.min_reads, a.min_bases, a.min_ratio, a.merge)
    with open(a.out + '.members.tsv', 'w') as w:
        w.write('place\tchrom\tstart\tend\tstrand\tsource\treads\tbest_ratio\tmedian_ratio\tmax_bases\tlocus\n')
        for r in rows: w.write(f'M_{r[0]}_{r[1]}_{r[3]}\t{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{r[4]}\t{r[5]}\t{r[6]:.3f}\t{r[7]:.3f}\t{r[8]}\t{r[9]}\n')
    new = collections.OrderedDict()
    for r in rows:
        if r[9] == '.': new.setdefault((r[0], r[3], r[1], r[2]), []).append(r)
    B = pysam.AlignmentFile(a.bam)   # the blocks of each model's record are read back here: the scan keeps no blocks (195 GB on a gorilla sample when it did)
    with open(a.out + '.members.gtf', 'w') as w:
        byc = collections.defaultdict(list)
        for i, h in enumerate(hits):
            if h[6] >= a.min_bases and h[7] >= a.min_ratio: byc[(h[0], h[1])].append(i)
        for (c, st, s, e), rs in new.items():
            cand = [i for i in byc[(c, st)] if hits[i][2] >= s and hits[i][3] <= e]
            if not cand: continue
            i = max(cand, key=lambda i: (hits[i][7] * hits[i][6], -hits[i][2])); ex = []; h = hits[i]
            bl = next((r.get_blocks() for r in B.fetch(c, h[2], h[2] + 1) if r.query_name == h[4] and r.reference_start == h[2] and not r.is_supplementary), [(h[2], h[3])])
            for x, y in bl:
                if ex and x - ex[-1][1] < 20: ex[-1] = (ex[-1][0], y)
                else: ex.append((x, y))
            at = f'gene_id "M_{c}_{s}_{st}"; transcript_id "M_{c}_{s}_{st}.1"; o2_status "member_candidate"; sources "{",".join(sorted({r[4] for r in rs}))}"; reads "{sum(r[5] for r in rs)}"; best_ratio "{max(r[6] for r in rs):.3f}";'
            w.write(f'{c}\tmember_discovery\ttranscript\t{ex[0][0] + 1}\t{ex[-1][1]}\t.\t{st}\t.\t{at}\n')
            for x, y in ex: w.write(f'{c}\tmember_discovery\texon\t{x + 1}\t{y}\t.\t{st}\t.\t{at}\n')
    print(f"member_discovery: {len(best)} reads, {len(hits)} hits off their source locus, {len(rows)} place x source rows ({len(new)} places on no assembled locus, {len({(r[0], r[1], r[3]) for r in rows if r[9] != '.'})} on one) -> {a.out}.members.tsv, .members.gtf")

G_REF = []
if __name__ == '__main__':
    if '--selftest' not in sys.argv and len(sys.argv) > 2:
        import pysam
        G_REF = list(pysam.AlignmentFile(sys.argv[2]).references)
    main()
