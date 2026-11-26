#!/usr/bin/env python3
"""Member discovery: the places every expressed locus's reads also align to (register 1910; O2, the `available` layer).

    member_discovery.py ASSEMBLY.gtf BAM OUT_PREFIX [--min-reads 5] [--min-bases 100] [--min-ratio 0.3] [--as-table TSV]
                        [--exclude BED] [--merge 0] [--repeats BED [--max-repeat 0.5] [--neighbour-bp 1000000]] [--threads N]

Loci = gene_id of the GTF (rule of eq_classes.py). A read's SOURCE = the locus of its primary record. A hit = a
secondary record of the read that is not on the source locus, with >= --min-bases aligned bases and AS >= --min-ratio x
the read's best AS (from --as-table, the driver's score table; without it a first pass finds it). Hits form a PLACE
when their aligned BLOCKS overlap on one strand: spans are not used, so an intron joins nothing (rows 1913 / 1914: places
merged by span chained neighbouring genes, median 21.5 kb on a human sample). A place is a member candidate of source S
when >= --min-reads reads of S hit it.
Writes OUT.members.tsv (place, chrom, start, end, strand, source, reads, best_ratio, median_ratio, max_bases, locus = the
assembled locus most of its hits lie on, else `.`, footprint_bp, footprint = its merged blocks) and OUT.members.gtf (one
model per place that is no assembled locus: the blocks of its best hit; o2_status "member_candidate").
It says where relatives are, nothing about their expression.
--repeats BED (opt-in; register 1956 / 1957): the table gains a last column `repeat_share` = the share of the place's footprint
inside the BED's intervals, and rows with a share >= --max-repeat (0.5) are left out of both files (2 = keep all, only the
column). On a human sample a third of the places at >= 0.9 of the best score are copies of a transposable element; the
filter also removes about a third of the real relatives (reads that end in a repeat-rich untranslated region).
A place on its source locus's contig within --neighbour-bp (1,000,000; 0 = no exemption) of it is kept whatever its share:
a tandem copy sits beside its source (register 1968). Give a BED of INTERSPERSED repeats only (SINE, LINE, LTR, DNA,
Retroposon): satellite and simple-repeat intervals cover much of a tandem array (rows 1967 / 1968: with both, the share of
places at a relative rises from 0.50 to 0.72 on an unused sample with 99-100% of the paralogues and pseudogenes kept).
--exclude BED: records that start inside an interval are not read (satellite arrays: one contig of A119b ran for hours).
Two passes over the BAM (primaries, then secondaries); the places of a contig are built in its worker.
--selftest runs the pure parts."""
import sys, os, re, argparse, bisect, collections, statistics, multiprocessing
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
OUTSIDE = 'outside'

def places(hits, merge=0):
    """hits [(chrom, strand, start, end, read, source, bases, ratio, locus, blocks)] -> [(chrom, strand, [hit indices])].
    Two hits are in one place when an aligned BLOCK of one overlaps (or lies within `merge` of) a block of the other:
    spans are not used, a spliced hit's introns join nothing (rows 1913 / 1914: span merging chained neighbouring genes)."""
    by = collections.defaultdict(list)
    for i, h in enumerate(hits):
        for a, b in h[9]: by[(h[0], h[1])].append((a, b, i))
    parent = list(range(len(hits)))
    def find(x):
        while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for k, v in by.items():
        v.sort(); ce, ci = v[0][1], v[0][2]
        for a, b, i in v[1:]:
            if a < ce + merge + (1 if merge else 0): parent[find(i)] = find(ci); ce = max(ce, b)
            else: ce, ci = b, i
    comp = collections.defaultdict(list)
    for i in range(len(hits)): comp[find(i)].append(i)
    return [(hits[v[0]][0], hits[v[0]][1], v) for _, v in sorted(comp.items(), key=lambda kv: (hits[kv[1][0]][0], hits[kv[1][0]][1], min(hits[i][2] for i in kv[1])))]

def footprint(hits, idx):
    o = []
    for a, b in sorted(x for i in idx for x in hits[i][9]):
        if o and a <= o[-1][1]: o[-1][1] = max(o[-1][1], b)
        else: o.append([a, b])
    return o

def candidates(hits, min_reads=2, min_bases=100, min_ratio=0.0, merge=0):
    """-> rows (chrom, start, end, strand, source, reads, best, median, max_bases, locus, footprint, best hit index in the filtered list), and that list"""
    hs = [h for h in hits if h[6] >= min_bases and h[7] >= min_ratio]; rows = []
    for c, st, idx in places(hs, merge):
        by = collections.defaultdict(dict); loc = collections.Counter(); mb = collections.Counter()
        for i in idx:
            h = hs[i]; by[h[5]][h[4]] = max(by[h[5]].get(h[4], 0), h[7]); mb[h[5]] = max(mb[h[5]], h[6])
            if h[8] != OUTSIDE: loc[h[8]] += 1
        l = loc.most_common(1)[0][0] if loc else '.'; fp = footprint(hs, idx)
        bi = max(idx, key=lambda i: (hs[i][7] * hs[i][6], -hs[i][2]))
        for src, rd in sorted(by.items()):
            if len(rd) >= min_reads: rows.append((c, fp[0][0], fp[-1][1], st, src, len(rd), max(rd.values()), statistics.median(rd.values()), mb[src], l, fp, bi))
    return rows, hs

G = {}
def _locus(r, c, gex, idx, E, overlap_genes):
    s = E.tx_strand(r); bl = r.get_blocks()
    ov = overlap_genes([(a + 1, b) for a, b in bl], gex[(c, s)], *idx[s]) if s in idx else None
    return (max(ov, key=lambda k: (ov[k], k)) if ov else OUTSIDE), s, bl

def _open(c):
    import pysam, bisect
    from gtf_names import overlap_genes
    import eq_classes as E
    gex = G['gex']
    idx = {s: ([a for a, _, _ in gex[(c, s)]], max(b - a for a, b, _ in gex[(c, s)])) for s in '+-' if (c, s) in gex}
    xs = G.get('excl', {}).get(c, ([], []))
    def skip(r):
        if r.is_unmapped or r.is_supplementary or not r.has_tag('AS'): return True
        if xs[0]:
            i = bisect.bisect_right(xs[0], r.reference_start) - 1
            if i >= 0 and r.reference_start < xs[1][i]: return True
        return False
    return pysam.AlignmentFile(G['bam']).fetch(c), skip, gex, idx, E, overlap_genes

def scan_primaries(c):
    """pass 1: the locus of every primary record (the read's source) and, without a score table, each read's best AS here"""
    it, skip, gex, idx, E, og = _open(c); src = {}; best = {}
    need = G['best'] is None
    for r in it:
        if skip(r): continue
        if need:
            a = r.get_tag('AS')
            if a > best.get(r.query_name, 0): best[r.query_name] = a
        if not r.is_secondary:
            g = _locus(r, c, gex, idx, E, og)[0]
            if g != OUTSIDE: src[r.query_name] = g
    return src, best

def scan_hits(c):
    """pass 2: the secondary records of reads with a source, off that source, that pass the thresholds; places of this contig"""
    it, skip, gex, idx, E, og = _open(c); hits = []
    SRC, BEST, mb, mr = G['src'], G['best'], G['min_bases'], G['min_ratio']
    for r in it:
        if not r.is_secondary or skip(r): continue
        src = SRC.get(r.query_name)
        if src is None: continue
        b = BEST.get(r.query_name, 0)
        if b <= 0: continue
        ratio = r.get_tag('AS') / b
        if ratio < mr: continue
        g, s, bl = _locus(r, c, gex, idx, E, og)
        if g == src: continue
        nb = sum(y - x for x, y in bl)
        if nb < mb: continue
        m = []
        for x, y in bl:
            if m and x - m[-1][1] < 20: m[-1] = (m[-1][0], y)
            else: m.append((x, y))
        hits.append((c, s, r.reference_start, r.reference_end, r.query_name, src, nb, ratio, g, tuple(m)))
    rows, hs = candidates(hits, G['min_reads'], mb, mr, G['merge'])
    return len(hits), [(r[:11], hs[r[11]][9]) for r in rows]

def selftest():
    import tempfile
    d = tempfile.mkdtemp(); open(d + '/r.bed', 'w').write('c\t100\t200\nc\t150\t300\nc\t900\t950\n'); bed = read_bed(d + '/r.bed')
    assert neighbour(('c', 5000, 6000), ['c', 100, 4000, '+', 'x'], 1000) and not neighbour(('c', 5001, 6000), ['c', 100, 4000, '+', 'x'], 1000) and not neighbour(('d', 5000, 6000), ['c', 100, 4000, '+', 'x'], 1000)
    assert neighbour(('c', 10, 90), ['c', 100, 4000, '+', 'x'], 1000) and not neighbour(('c', 5000, 6000), None, 1000) and not neighbour(('c', 5000, 6000), ['c', 100, 4000, '+', 'x'], 0)
    assert bed == {'c': ([100, 900], [300, 950])} and inside([(0, 200), (250, 350)], bed['c']) == 0.5 and inside([(0, 10)], None) == 0.0
    B = lambda *x: tuple(x)
    h = [('c', '+', 100, 400, 'r1', 'A', 300, 0.9, OUTSIDE, B((100, 400))), ('c', '+', 150, 9000, 'r2', 'A', 350, 0.6, OUTSIDE, B((150, 300), (8800, 9000))),
         ('c', '+', 5000, 5300, 'r3', 'A', 300, 0.9, OUTSIDE, B((5000, 5300))),            # inside r2's intron: NOT joined to it
         ('c', '-', 120, 380, 'r4', 'A', 260, 0.9, OUTSIDE, B((120, 380))), ('c', '+', 8900, 9100, 'r5', 'B', 200, 0.4, 'L9', B((8900, 9100))),
         ('c', '+', 8950, 9200, 'r6', 'B', 80, 0.9, 'L9', B((8950, 9200))), ('c', '+', 9050, 9300, 'r7', 'B', 250, 0.5, 'L9', B((9050, 9300)))]
    p = places(h); assert sorted(sorted(x[2]) for x in p) == [[0, 1, 4, 5, 6], [2], [3]], p
    r, hs = candidates(h); assert [x[:10] for x in r] == [('c', 100, 9300, '+', 'A', 2, 0.9, 0.75, 350, 'L9'), ('c', 100, 9300, '+', 'B', 2, 0.5, 0.45, 250, 'L9')], r
    assert r[0][10] == [[100, 400], [8800, 9300]] and footprint(h, [2]) == [[5000, 5300]]
    r, hs = candidates(h, min_ratio=0.5); assert [x[:10] for x in r] == [('c', 100, 9000, '+', 'A', 2, 0.9, 0.75, 350, '.')] and candidates(h, min_reads=3)[0] == []
    r, hs = candidates([x for x in h if x[4] != 'r2'], min_reads=1); assert sorted((x[1], x[2], x[4]) for x in r) == [(100, 400, 'A'), (120, 380, 'A'), (5000, 5300, 'A'), (8900, 9300, 'B')]
    print('member_discovery selftest ok')

def read_bed(path):
    """-> {chrom: (starts, ends)} of the BED's intervals, merged and sorted"""
    iv = collections.defaultdict(list)
    for l in open(path):
        f = l.split('\t')
        if len(f) >= 3 and f[1].isdigit(): iv[f[0]].append((int(f[1]), int(f[2])))
    ex = {}
    for c, v in iv.items():
        v.sort(); m = []
        for x, y in v:
            if m and x <= m[-1][1]: m[-1][1] = max(m[-1][1], y)
            else: m.append([x, y])
        ex[c] = ([x for x, _ in m], [y for _, y in m])
    return ex

def inside(blocks, bed):
    """share of the bases of `blocks` inside (starts, ends) of read_bed"""
    if not bed: return 0.0
    st, en = bed; n = 0
    for a, b in blocks:
        k = max(0, bisect.bisect_right(st, a) - 1)
        while k < len(st) and st[k] < b: n += max(0, min(b, en[k]) - max(a, st[k])); k += 1
    return n / max(1, sum(b - a for a, b in blocks))

def neighbour(row, src, bp):
    """row (chrom, start, end, ...) lies on the contig of its source locus [chrom, start, end, ...] within bp of it"""
    return bool(bp) and src is not None and src[0] == row[0] and max(row[1] - src[2], src[1] - row[2]) <= bp

def main():
    if '--selftest' in sys.argv: return selftest()
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('gtf'); ap.add_argument('bam'); ap.add_argument('out'); ap.add_argument('--min-reads', type=int, default=5); ap.add_argument('--min-bases', type=int, default=100)
    ap.add_argument('--min-ratio', type=float, default=0.3); ap.add_argument('--merge', type=int, default=0); ap.add_argument('--as-table'); ap.add_argument('--exclude'); ap.add_argument('--repeats'); ap.add_argument('--max-repeat', type=float, default=0.5); ap.add_argument('--neighbour-bp', type=int, default=1_000_000); ap.add_argument('--threads', type=int, default=8)
    a = ap.parse_args()
    import pysam, eq_classes as E
    meta, gex = E.load_loci(a.gtf); G.update(gex=gex, bam=a.bam, best=None, min_reads=a.min_reads, min_bases=a.min_bases, min_ratio=a.min_ratio, merge=a.merge)
    if a.exclude: G['excl'] = read_bed(a.exclude)
    if a.as_table:
        best = {}
        for l in open(a.as_table):
            if l[0] == '#':       # the header records the BAM's size: the scores of another BAM give wrong ratios without any error
                m = re.search(r'\tbam_bytes=(\d+)', l)
                if m and int(m.group(1)) != os.path.getsize(a.bam): sys.exit(f'member_discovery: {a.as_table} was made from a BAM of {m.group(1)} bytes; {a.bam} has {os.path.getsize(a.bam)}')
            else:
                f = l.split('\t', 2); best[f[0]] = int(f[1]) if f[1].lstrip('-').isdigit() else 0
        G['best'] = best
    refs = list(pysam.AlignmentFile(a.bam).references)
    with multiprocessing.get_context('fork').Pool(a.threads) as pool: res = pool.map(scan_primaries, refs, chunksize=1)
    src = {}; best = G['best'] if G['best'] is not None else {}
    for s_, b_ in res:
        src.update(s_)
        if G['best'] is None:
            for q, v in b_.items():
                if v > best.get(q, 0): best[q] = v
    del res; G['src'] = src; G['best'] = best
    with multiprocessing.get_context('fork').Pool(a.threads) as pool: res = pool.map(scan_hits, refs, chunksize=1)
    nh = sum(r[0] for r in res); rows = [x for r in res for x in r[1]]; share = {}; nrep = 0
    if a.repeats:
        rp = read_bed(a.repeats); share = {id(x[0]): inside(x[0][10], rp.get(x[0][0])) for x in rows}; nrep = len(rows)
        rows = [x for x in rows if share[id(x[0])] < a.max_repeat or neighbour(x[0], meta.get(x[0][4]), a.neighbour_bp)]; nrep -= len(rows)
    with open(a.out + '.members.tsv', 'w') as w:
        w.write('place\tchrom\tstart\tend\tstrand\tsource\treads\tbest_ratio\tmedian_ratio\tmax_bases\tlocus\tfootprint_bp\tfootprint' + ('\trepeat_share' if a.repeats else '') + '\n')
        for r, _ in rows:
            fp = r[10]; tail = f'\t{share[id(r)]:.3f}' if a.repeats else ''; w.write(f"M_{r[0]}_{r[1]}_{r[3]}\t{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{r[4]}\t{r[5]}\t{r[6]:.3f}\t{r[7]:.3f}\t{r[8]}\t{r[9]}\t{sum(b - a for a, b in fp)}\t{','.join(f'{x}-{y}' for x, y in fp[:400])}{tail}\n")
    done = set(); nnew = 0; sources = collections.defaultdict(set)      # every source of a place, in one pass (was a scan of all rows per place, cut to one source past 200,000 rows)
    for x, _ in rows: sources[(x[0], x[1], x[3])].add(x[4])
    with open(a.out + '.members.gtf', 'w') as w:
        for r, bl in rows:
            k = (r[0], r[1], r[3])
            if r[9] != '.' or k in done: continue
            done.add(k); nnew += 1; srcs = sorted(sources[k])
            at = f'gene_id "M_{r[0]}_{r[1]}_{r[3]}"; transcript_id "M_{r[0]}_{r[1]}_{r[3]}.1"; o2_status "member_candidate"; sources "{",".join(srcs)}"; reads "{r[5]}"; best_ratio "{r[6]:.3f}";'
            w.write(f'{r[0]}\tmember_discovery\ttranscript\t{bl[0][0] + 1}\t{bl[-1][1]}\t.\t{r[3]}\t.\t{at}\n')
            for x, y in bl: w.write(f'{r[0]}\tmember_discovery\texon\t{x + 1}\t{y}\t.\t{r[3]}\t.\t{at}\n')
    if a.repeats: print(f'member_discovery: {nrep} rows left out, repeat share >= {a.max_repeat} and not within {a.neighbour_bp} bp of their source ({a.repeats})')
    print(f"member_discovery: {len(src)} reads with a source locus, {nh} secondary hits off it past the thresholds, {len(rows)} place x source rows ({nnew} places on no assembled locus, {len({(r[0], r[1], r[3]) for r, _ in rows if r[9] != '.'})} on one) -> {a.out}.members.tsv, .members.gtf")

if __name__ == '__main__':
    main()
