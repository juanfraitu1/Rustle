#!/usr/bin/env python3
"""Read-derived gene nodes for the de novo mode (development rule of ledger §6kh): same-strand primary MAPQ >= 1 reads;
strong transcript end sites = clusters (within W bp) of >= 3 read 5' starts or 3' ends where those reads are >= FRAC of
the reads covering the site; a read that crosses a strong site is clipped to its 5' side of the first crossed site;
clipped reads grouped by exon overlap; exons = bases at read depth >= 2; groups split into sub-loci linked by >= 2 reads
(`denovo_shared_def.split_linked`); sub-loci with >= 3 supporting reads and >= 100 exonic bp are nodes.

usage: read_gene_nodes.py BAM CONTIGS(comma) OUT.nodes.tsv [W=50] [FRAC=0.5] [nolink]
"""
import csv, collections, sys, pysam, statistics
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
import denovo_shared_def as d
bam_path, contigs, out_path = sys.argv[1], sys.argv[2].split(","), sys.argv[3]
W = int(sys.argv[4]) if len(sys.argv) > 4 else 50
FRAC = float(sys.argv[5]) if len(sys.argv) > 5 else 0.5
NOLINK = len(sys.argv) > 6 and sys.argv[6] == 'nolink'
bam = pysam.AlignmentFile(bam_path)
reads = []
for c in contigs:
    for r in bam.fetch(c):
        if r.is_unmapped or r.is_secondary or r.is_supplementary or r.mapping_quality < 1: continue
        st = "-" if r.is_reverse else "+"
        if r.has_tag("ts") and r.get_tag("ts") == "-": st = "+" if st == "-" else "-"
        blocks, pos, cur = [], r.reference_start, r.reference_start
        for op, n in r.cigartuples:
            if op in (0, 2, 7, 8): pos += n
            elif op == 3:
                if pos > cur: blocks.append((cur, pos))
                pos += n; cur = pos
        if pos > cur: blocks.append((cur, pos))
        reads.append([c, st, blocks])
# strong end sites per (chrom, strand): cluster 5' starts and 3' ends within W; strong if >= 3 reads and >= FRAC of reads covering the site end/start there
def clusters(pos):
    pos.sort(); out = []
    for p in pos:
        if out and p - out[-1][-1] <= W: out[-1].append(p)
        else: out.append([p])
    return [(int(statistics.median(g)), len(g)) for g in out if len(g) >= 3]
by = collections.defaultdict(list)
for i, (c, st, bl) in enumerate(reads): by[(c, st)].append(i)
cut_sites = collections.defaultdict(list)  # (c, st) -> positions where genes are separated
for (c, st), ids in by.items():
    ends3 = [(reads[i][2][-1][1] if st == "+" else reads[i][2][0][0]) for i in ids]
    ends5 = [(reads[i][2][0][0] if st == "+" else reads[i][2][-1][1]) for i in ids]
    spans = sorted((reads[i][2][0][0], reads[i][2][-1][1]) for i in ids)
    starts = [s for s, e in spans]
    import bisect
    for site, n in clusters(ends3) + clusters(ends5):
        lo = bisect.bisect_left(starts, site - 10**7)
        cover = sum(1 for s, e in spans[:bisect.bisect_right(starts, site)] if e > site)
        if cover and n >= FRAC * cover:
            cut_sites[(c, st)].append(site)
n_cut = 0
for (c, st), ids in by.items():
    sites = sorted(set(cut_sites[(c, st)]))
    if not sites: continue
    for i in ids:
        bl = reads[i][2]
        s5 = bl[0][0] if st == "+" else bl[-1][1]
        s3 = bl[-1][1] if st == "+" else bl[0][0]
        lo, hi = min(s5, s3), max(s5, s3)
        inner = [x for x in sites if lo + W < x < hi - W]
        if not inner: continue
        # keep only the part of the read on its 5' side of the first crossed site (transcript direction)
        if st == "+":
            cutp = min(inner); newb = [(a, min(b, cutp)) for a, b in bl if a < cutp]
        else:
            cutp = max(inner); newb = [(max(a, cutp), b) for a, b in bl if b > cutp]
        newb = [(a, b) for a, b in newb if b > a]
        if newb and newb != bl: reads[i][2] = newb; n_cut += 1
print(f"reads {len(reads)}; strong end sites {sum(len(set(v)) for v in cut_sites.values())}; reads clipped at a crossed site {n_cut}")
# read components (same rule as AF-3 split, alpha 0) on the clipped reads
items = [((c, st), bl) for c, st, bl in reads]
parent = list(range(len(reads)))
def find(x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]; x = parent[x]
    return x
byk = collections.defaultdict(list)
for i, (k, bl) in enumerate(items):
    for s, e in bl: byk[k].append((s, e, i))
for iv in byk.values():
    iv.sort(); end, owner = -1, None
    for s, e, i in iv:
        if owner is not None and s < end: parent[find(i)] = find(owner)
        if e > end: end, owner = e, i
groups = collections.defaultdict(list)
for i in range(len(reads)): groups[find(i)].append(i)
nodes = []
for mem in groups.values():
    if len(mem) < 3: continue
    c, st = reads[mem[0]][0], reads[mem[0]][1]
    ex = d.depth2_exons([reads[i][2] for i in mem])
    if NOLINK:
        if sum(e - s for s, e in ex) >= 100:
            nodes.append((c, ex[0][0], ex[-1][1], st, len(mem), ",".join(f"{a}-{b}" for a, b in ex)))
        continue
    for ks, sup in d.split_linked(ex, [reads[i][2] for i in mem]):
        sub = [ex[k] for k in ks]
        if sup >= 3 and sum(e - s for s, e in sub) >= 100:
            nodes.append((c, sub[0][0], sub[-1][1], st, sup, ",".join(f"{a}-{b}" for a, b in sub)))
nodes.sort()
with open(out_path, "w") as fh:
    fh.write("idx\tchrom\tstart\tend\tstrand\tn_exon\tn_reads\texons\n")
    for i, (c, s, e, st, n, ex) in enumerate(nodes):
        fh.write(f"{i}\t{c}\t{s}\t{e}\t{st}\t{ex.count(',') + 1}\t{n}\t{ex}\n")
print("nodes", len(nodes))
