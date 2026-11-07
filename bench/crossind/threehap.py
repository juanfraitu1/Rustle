#!/usr/bin/env python3
"""Three-haplotype test (2026-11-07): a diploid shows <= 2 haplotypes at any locus; >= 3 read-supported haplotypes at a
pair of variable sites means a copy the reference lacks is collapsed there.
  threehap.py BAM CLUSTERS.tsv OUT [KEEP=1.0] [N_SINGLE=1500]
Reads with de > 0.02 dropped. Loci = piles.py clusters (all multi-copy + a fixed hash sample of single-copy as control). Primary reads only, mismatches
from --eqx X ops (indels ignored: HiFi homopolymer errors). Site = alt base on >= 3 reads and >= 15% of covering reads.
Site pair (within 2 kb) has >= 3 haplotypes when each has >= 3 reads and >= 10% of the reads spanning both sites.
Locus flagged when >= 2 site pairs that share no site have >= 3 haplotypes (one bad column cannot flag a locus)."""
import sys, zlib, collections, itertools
import pysam

bam, cl, out = sys.argv[1:4]
keep = float(sys.argv[4]) if len(sys.argv) > 4 else 1.0
nsingle = int(sys.argv[5]) if len(sys.argv) > 5 else 1500
loci = [l.rstrip('\n').split('\t') for l in open(cl) if not l.startswith('contig')]
multi = [x for x in loci if x[4] == '1']
single = sorted((x for x in loci if x[4] == '0'), key=lambda x: zlib.crc32(f'{x[0]}:{x[1]}'.encode()))[:nsingle]
B = pysam.AlignmentFile(bam)

def three_hap(c, s, e):
    reads = []
    for a in B.fetch(c, s, e):
        if a.is_secondary or a.is_supplementary or a.is_unmapped:
            continue
        if a.get_tag('de') > 0.02:  # junk tail (OR6737 de p95 0.0285) makes random third haplotypes
            continue
        if keep < 1 and zlib.crc32(a.query_name.encode()) / 2**32 >= keep:
            continue
        mm, rp, qp, q = {}, a.reference_start, 0, a.query_sequence
        for op, n in a.cigartuples:
            if op == 8:
                for i in range(n): mm[rp + i] = q[qp + i]
            if op in (0, 7, 8): rp += n; qp += n
            elif op in (1, 4): qp += n
            elif op in (2, 3): rp += n
        reads.append((a.get_blocks(), mm))
    if len(reads) < 10:
        return None
    alt = collections.Counter((p, b) for _, mm in reads for p, b in mm.items())
    cand = {p: b for (p, b), k in alt.items() if k >= 3}
    def covers(blocks, p): return any(x <= p < y for x, y in blocks)
    cov = {p: [i for i, (bl, _) in enumerate(reads) if covers(bl, p)] for p in cand}
    sites = sorted(p for p in cand if alt[p, cand[p]] >= 0.15 * len(cov[p]))
    hits = []
    for p1, p2 in itertools.combinations(sites, 2):
        if p2 - p1 > 2000: continue
        span = set(cov[p1]) & set(cov[p2])
        if len(span) < 9: continue
        hap = collections.Counter((reads[i][1].get(p1) == cand[p1], reads[i][1].get(p2) == cand[p2]) for i in span)
        if sum(k >= 3 and k >= 0.1 * len(span) for k in hap.values()) >= 3:
            hits.append((p1, p2))
    used, indep = set(), 0  # greedy count of site-disjoint flagged pairs
    for p1, p2 in hits:
        if p1 not in used and p2 not in used:
            used |= {p1, p2}; indep += 1
    return len(reads), len(sites), indep

tot = collections.Counter()
with open(out + '.threehap.tsv', 'w') as w:
    w.write('contig\tstart\tend\tmulticopy\treads\tsites\tindep_3hap_pairs\tflag\n')
    for x in multi + single:
        r = three_hap(x[0], int(x[1]), int(x[2]))
        if r is None: continue
        f = r[2] >= 2
        tot[x[4], f] += 1
        w.write('\t'.join(x[:3] + [x[4]] + list(map(str, r)) + [str(int(f))]) + '\n')
for m, name in (('1', 'multi '), ('0', 'single')):
    k, n = tot[m, True], tot[m, True] + tot[m, False]
    print(f'{bam}\tkeep={keep}\t{name}\tflagged {k}/{n}\t{k/max(n,1):.4f}')
