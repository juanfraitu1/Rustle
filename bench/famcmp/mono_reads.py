#!/usr/bin/env python3
"""Per reference single-exon transcript: unspliced primary alignments that overlap it by >= 50% of both lengths, on its strand
and on any strand, plus spliced primaries covering >= 50% of it. One pass over the BAM. usage: mono_reads.py REF_GTF BAM OUT_TSV"""
import sys, re, gzip, collections, pysam
T = collections.defaultdict(list)
for l in (gzip.open(sys.argv[1], 'rt') if sys.argv[1].endswith('.gz') else open(sys.argv[1])):
    f = l.split('\t')
    if len(f) > 8 and f[2] == 'exon': T[(f[0], f[6], re.search(r'transcript_id "([^"]*)"', f[8]).group(1))].append((int(f[3]) - 1, int(f[4])))
bins = collections.defaultdict(list); cnt = {}
for k, e in T.items():
    if len(e) != 1: continue
    a, b = e[0]; cnt[k] = [a, b, 0, 0, 0]
    for i in range(a // 10000, b // 10000 + 1): bins[(k[0], i)].append(k)
for r in pysam.AlignmentFile(sys.argv[2]):
    if r.is_unmapped or r.is_secondary or r.is_supplementary: continue
    c, s, e = r.reference_name, r.reference_start, r.reference_end; spl = 'N' in r.cigarstring
    for k in {k for i in range(s // 10000, e // 10000 + 1) for k in bins.get((c, i), ())}:
        v = cnt[k]; ov = min(v[1], e) - max(v[0], s)
        if spl: v[4] += ov >= 0.5 * (v[1] - v[0])
        elif ov >= 0.5 * max(v[1] - v[0], e - s): v[3] += 1; v[2] += r.is_reverse == (k[1] == '-')
with open(sys.argv[3], 'w') as w:
    for (c, s, t), v in sorted(cnt.items()): w.write(f'{c}\t{v[0]}\t{v[1]}\t{s}\t{t}\t{v[2]}\t{v[3]}\t{v[4]}\n')
