#!/usr/bin/env python3
"""PREREG (docs/PREREG_junction_fuzz_2026-09-15.md): measure real per-junction alignment jitter on chr20,
BEFORE any fuzzy-merge code is written or evaluated against a metric. For every real chr20 read's spliced
junction that sits within CAPTURE_BP of a real annotated RefSeq intron boundary (matched by donor-site
distance), record the signed offset at both the donor and acceptor site. Since merge_fuzzy_skeletons
requires BOTH donor AND acceptor to be within tolerance (AND condition), the relevant statistic is the
90th percentile of per-junction max(donor_offset, acceptor_offset), not the pooled distribution.

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
donor_offsets = []
acceptor_offsets = []
max_offsets = []  # per-junction max(donor, acceptor)
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
        don_offset = abs(don - rd)
        acc_offset = abs(acc - ra)
        donor_offsets.append(don_offset)
        acceptor_offsets.append(acc_offset)
        max_offsets.append(max(don_offset, acc_offset))

donor_offsets.sort()
acceptor_offsets.sort()
max_offsets.sort()
pooled = donor_offsets + acceptor_offsets
pooled.sort()

# Statistics
n_pooled = len(pooled)
n_max = len(max_offsets)
p50_pooled = pooled[int(n_pooled * 0.50)]
p90_pooled = pooled[int(n_pooled * 0.90)]
p95_pooled = pooled[int(n_pooled * 0.95)]

p50_donor = donor_offsets[int(n_matched * 0.50)]
p90_donor = donor_offsets[int(n_matched * 0.90)]
p95_donor = donor_offsets[int(n_matched * 0.95)]

p50_acceptor = acceptor_offsets[int(n_matched * 0.50)]
p90_acceptor = acceptor_offsets[int(n_matched * 0.90)]
p95_acceptor = acceptor_offsets[int(n_matched * 0.95)]

p50_max = max_offsets[int(n_max * 0.50)]
p90_max = max_offsets[int(n_max * 0.90)]
p95_max = max_offsets[int(n_max * 0.95)]

print(f'reads scanned: {n_reads}')
print(f'junctions seen: {n_junctions}')
print(f'junctions matched to an annotated intron within {CAPTURE_BP}bp: {n_matched}')
print()
print('=== PER-AXIS STATISTICS (for context) ===')
print(f'donor-only |offset| samples: {len(donor_offsets)}')
print(f'  median: {p50_donor}, p90: {p90_donor}, p95: {p95_donor}')
print(f'acceptor-only |offset| samples: {len(acceptor_offsets)}')
print(f'  median: {p50_acceptor}, p90: {p90_acceptor}, p95: {p95_acceptor}')
print()
print('=== POOLED BOTH AXES (methodologically flawed for AND-condition merge test) ===')
print(f'pooled donor+acceptor |offset| samples: {n_pooled}')
print(f'  median: {p50_pooled}, p90: {p90_pooled}, p95: {p95_pooled}')
print()
print('=== PER-JUNCTION MAX(donor, acceptor) (CORRECT for AND-condition merge test) ===')
print(f'per-junction max(|donor_offset|, |acceptor_offset|) samples: {n_max}')
print(f'  median: {p50_max}')
print(f'  90th percentile: {p90_max}   <-- THIS IS RUSTLE_JUNCTION_FUZZ_BP\'s pre-registered value')
print(f'  95th percentile: {p95_max}')
