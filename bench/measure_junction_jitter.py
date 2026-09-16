#!/usr/bin/env python3
"""PREREG (docs/PREREG_junction_fuzz_2026-09-15.md): measure real per-junction alignment jitter on chr20,
BEFORE any fuzzy-merge code is written or evaluated against a metric. For every real chr20 read's spliced
junction that sits within CAPTURE_BP of a real annotated RefSeq intron boundary (matched by donor-site
distance), record the signed offset at both the donor and acceptor site. The 90th percentile of the
pooled |offset| becomes RUSTLE_JUNCTION_FUZZ_BP's pre-registered value.

usage: python3 bench/measure_junction_jitter.py <chr20.bam> <chr20_ref.gtf>
"""
import sys, re, collections, pysam

bam_p, gtf_p = sys.argv[1], sys.argv[2]
CAPTURE_BP = 500

# annotated introns: for each ref transcript, gaps between consecutive exons
ex = collections.defaultdict(list)
for line in open(gtf_p):
    if line.startswith('#'):
        continue
    f = line.rstrip('\n').split('\t')
    if len(f) < 9 or f[2] != 'exon':
        continue
    m = re.search(r'transcript_id "([^"]+)"', f[8])
    if not m:
        continue
    ex[m.group(1)].append((int(f[3]) - 1, int(f[4]), f[0]))

ann_introns = collections.defaultdict(set)  # chrom -> set of (start, end)
for t, v in ex.items():
    v.sort()
    for a, b in zip(v, v[1:]):
        ann_introns[a[2]].add((a[1], b[0]))
ann_sorted = {c: sorted(v) for c, v in ann_introns.items()}


def introns_of(pos, cig):
    o, p = [], pos
    for n, op in re.findall(r'(\d+)([MIDNSHP=X])', cig):
        n = int(n)
        if op in 'M=XD':
            p += n
        elif op == 'N':
            o.append((p, p + n))
            p += n
    return o


def nearest(chrom, don):
    """closest annotated intron by donor-site distance, within CAPTURE_BP, else None."""
    best, bd = None, CAPTURE_BP + 1
    for a, b in ann_sorted.get(chrom, ()):
        d = abs(a - don)
        if d < bd:
            best, bd = (a, b), d
    return best if bd <= CAPTURE_BP else None


bam = pysam.AlignmentFile(bam_p, 'rb')
offsets = []
n_reads = n_junctions = n_matched = 0
for rec in bam:
    if rec.is_unmapped or rec.is_secondary or rec.is_supplementary or not rec.cigarstring:
        continue
    n_reads += 1
    chrom = rec.reference_name
    for don, acc in introns_of(rec.reference_start, rec.cigarstring):
        n_junctions += 1
        ref = nearest(chrom, don)
        if ref is None:
            continue
        n_matched += 1
        rd, ra = ref
        offsets.append(abs(don - rd))
        offsets.append(abs(acc - ra))

offsets.sort()
n = len(offsets)
p50 = offsets[int(n * 0.50)]
p90 = offsets[int(n * 0.90)]
p95 = offsets[int(n * 0.95)]
print(f'reads scanned: {n_reads}')
print(f'junctions seen: {n_junctions}')
print(f'junctions matched to an annotated intron within {CAPTURE_BP}bp: {n_matched}')
print(f'pooled donor+acceptor |offset| samples: {n}')
print(f'median |offset|: {p50}')
print(f'90th percentile |offset|: {p90}   <-- this is RUSTLE_JUNCTION_FUZZ_BP\'s pre-registered value')
print(f'95th percentile |offset|: {p95}')
