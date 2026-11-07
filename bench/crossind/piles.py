#!/usr/bin/env python3
"""Cross-individual missing-copy test (2026-11-07). Does a non-Jim gorilla's Iso-Seq on Jim's T2T show more
divergence piles (reads of a copy the reference lacks, forced onto a relative) than Jim's own reads?
  piles.py BAM OUT [KEEP_FRACS=1.0]   (comma list: deterministic read-name subsamples, one pass)
Primary reads, clustered by overlapping landing. Read is multi-copy when a secondary has AS >= 0.9 x primary AS.
Cluster multi-copy when >= 30% of its reads are. Read de thresholds T = 0.005 and 0.01 (paralog-level divergence); the
self-aligned sample (Jim) is the SNP + error baseline; T95 is printed only as a description. Pile = >= 10 reads with de > T and >= 25% of the cluster. Only clusters with >= 10 reads count."""
import sys, zlib, collections
import pysam

bam, out = sys.argv[1:3]
keeps = [float(x) for x in sys.argv[3].split(',')] if len(sys.argv) > 3 else [1.0]
prim, sec = [], collections.defaultdict(int)
for a in pysam.AlignmentFile(bam, threads=8):
    if a.is_unmapped or a.is_supplementary:
        continue
    if a.is_secondary:
        sec[a.query_name] = max(sec[a.query_name], a.get_tag('AS'))
    else:
        prim.append((a.reference_name, a.reference_start, a.reference_end, a.query_name, a.get_tag('de'), a.get_tag('AS')))
prim.sort()
prim_all = prim
for keep in keeps:
  prim = [r for r in prim_all if keep >= 1 or zlib.crc32(r[3].encode()) / 2**32 < keep]
  if True:
    clusters, cur = [], None
    for c, s, e, q, de, AS in prim:
        mc = sec.get(q, -1) >= 0.9 * AS
        if cur and cur[0] == c and s < cur[2]:
            cur[2] = max(cur[2], e); cur[3].append((de, mc))
        else:
            cur = [c, s, e, [(de, mc)]]; clusters.append(cur)
    clusters = [x for x in clusters if len(x[3]) >= 10]
    for x in clusters:
        x.append(sum(m for _, m in x[3]) >= 0.3 * len(x[3]))
    single = sorted(de for x in clusters if not x[4] for de, _ in x[3])
    T95 = single[int(0.95 * (len(single) - 1))]
    def is_pile(rs, T):
        hi = sum(de > T for de, _ in rs)
        return hi, hi >= 10 and hi >= 0.25 * len(rs)
    res = {}
    with open(f'{out}.k{keep}.clusters.tsv', 'w') as w:
        w.write('contig\tstart\tend\treads\tmulticopy\tmedian_de\thigh_005\tpile_005\thigh_01\tpile_01\n')
        for c, s, e, rs, mc in clusters:
            row = [c, s, e, len(rs), int(mc), f'{sorted(d for d, _ in rs)[len(rs)//2]:.5f}']
            for T in (0.005, 0.01):
                hi, pile = is_pile(rs, T)
                res[T, mc, pile] = res.get((T, mc, pile), 0) + 1
                row += [hi, int(pile)]
            w.write('\t'.join(map(str, row)) + '\n')
    print(f'{bam}\tkeep={keep}\tprimary_reads={len(prim)}\tclusters={len(clusters)}\tT95={T95:.5f}')
    for T in (0.005, 0.01):
        for mc in (False, True):
            k, m = res.get((T, mc, True), 0), res.get((T, mc, False), 0)
            print(f'  T={T:.4f}\t{"multi " if mc else "single"}\tpiles {k}/{k+m}\t{k/max(k+m,1):.4f}')
