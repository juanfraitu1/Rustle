#!/usr/bin/env python3
"""Score a --from-genome copies.tsv against the Soto 83-family benchmark, using the SAME member-recovery
rule as the RNA scorer: a Soto member is RECOVERED iff a >=2-copy DNA family has a copy locus overlapping it.
Usage: score_from_genome.py <dna_mode.copies.tsv> [bed=bench/soto/80_fams.chr.bed]"""
import sys
from collections import defaultdict

copies_tsv = sys.argv[1]
BED = sys.argv[2] if len(sys.argv) > 2 else "bench/soto/80_fams.chr.bed"

# family_id -> list of (chrom, start, end)
fam = defaultdict(list)
for i, ln in enumerate(open(copies_tsv)):
    if i == 0:
        continue
    f = ln.rstrip("\n").split("\t")  # family_id copy_idx tid chrom start end n_exon strand n_reads
    if len(f) >= 6:
        fam[f[0]].append((f[3], int(f[4]), int(f[5])))

# only >=2-copy families count as detections
copyloci = [(c, s, e) for _fid, locs in fam.items() if len(locs) >= 2 for (c, s, e) in locs]

members = []
for ln in open(BED):
    c, s, e, name, *_ = ln.rstrip("\n").split("\t")
    members.append((name.split("|")[1], name.split("|")[0], c, int(s), int(e)))

def hit(c, s, e):
    return any(cc == c and s < ee and e > ss for (cc, ss, ee) in copyloci)

rec = sum(1 for (_f, _g, c, s, e) in members if hit(c, s, e))
n = len(members)
print(f"DNA member recovery: {rec}/{n} = {100*rec/n:.1f}%")
print(f"  ({len(fam)} DNA families, {sum(1 for l in fam.values() if len(l)>=2)} with >=2 copies; {len(copyloci)} scored copy loci)")
