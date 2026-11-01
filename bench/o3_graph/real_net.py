#!/usr/bin/env python3
"""copy_phasing inputs from a real BAM and a run's copies table (2026-11-01, r1477 port measurement): the families with >= 2 copies,
and a net = every read with a non-supplementary record overlapping a copy of the family (all used), the P16 / P18 convention.
  real_net.py BAM FAM.copies.tsv OUT_PREFIX   -> OUT_PREFIX.copies.tsv (the >= 2-copy families' rows, as given) and OUT_PREFIX.net.tsv
                                                (family, read, used=1, via=B)"""
import sys, csv, collections
import pysam
bam, cp, out = sys.argv[1:4]
rows = list(csv.DictReader(open(cp), delimiter='\t'))
by = collections.defaultdict(list)
for r in rows: by[r['family_id']].append(r)
multi = {f: v for f, v in by.items() if len(v) >= 2}
with open(f"{out}.copies.tsv", 'w') as w:
    w.write('\t'.join(rows[0].keys()) + '\n')
    for f in multi:
        for r in multi[f]: w.write('\t'.join(r.values()) + '\n')
b = pysam.AlignmentFile(bam)
n = 0
with open(f"{out}.net.tsv", 'w') as w:
    w.write("family\tread\tused\tvia\n")
    for f, v in multi.items():
        reads = set()
        for r in v:
            ex = sorted(tuple(map(int, x.split('-'))) for x in r['exons'].split(',') if x)
            for rec in b.fetch(r['chrom'], ex[0][0], ex[-1][1]):
                if rec.is_unmapped or rec.is_supplementary: continue
                reads.add(rec.query_name)
        for q in sorted(reads): w.write(f"{f}\t{q}\t1\tB\n")
        n += len(reads)
print(f"{len(multi)} families with >= 2 copies, {n} net rows", file=sys.stderr)
