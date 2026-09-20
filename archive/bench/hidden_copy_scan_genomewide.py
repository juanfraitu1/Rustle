#!/usr/bin/env python3
"""Genome-wide hidden-copy scan: the detect_hidden_copy signal at EVERY expressed de-novo locus
(not just the 70 multi-copy families) -- so a single-copy reference gene expressed as >=2 copies
(a CNV the reference collapsed) is now catchable. Reuses bench/hidden_copy_scan.py's detector."""
import sys, json, time
sys.path.insert(0, "bench")
import copy_assign as CA
import pysam
from hidden_copy_scan import alts_in_span, detect_hidden, MIN_DEPTH

READ_CAP = 800  # cap reads/locus (signal is statistical; bounds the few huge collapsed loci)
loci = json.load(open("/home/juanfra/winloci_scratch/refabsent/genomewide_loci.json"))
bam = pysam.AlignmentFile(CA.GGO_BAM, "rb")
flags = []
t0 = time.time()
for i, (chrom, s, e) in enumerate(loci):
    if i % 1000 == 0:
        print(f"  {i}/{len(loci)} loci ({time.time()-t0:.0f}s), {len(flags)} flagged", flush=True)
    reads = []
    for aln in bam.fetch(chrom, s, e):
        if aln.is_secondary or aln.is_supplementary or aln.is_unmapped:
            continue
        reads.append((aln.reference_start, aln.reference_end, alts_in_span(aln, s, e)))
        if len(reads) >= READ_CAP:
            break
    if len(reads) < MIN_DEPTH:
        continue
    ev = detect_hidden(reads)
    if ev["flagged"]:
        flags.append(dict(chrom=chrom, start=s, end=e, **ev))
flags.sort(key=lambda r: -r["n_alt_positions"])
json.dump(flags, open("/home/juanfra/winloci_scratch/refabsent/genomewide_flags.json", "w"), indent=1)
print(f"\nDONE {len(loci)} loci in {time.time()-t0:.0f}s -> {len(flags)} flagged loci (hidden 2nd haplotype)")
# how many are at NEW (non-family) loci vs the 70 families already known?
meta = CA.load_meta(); skel = CA.load_skel(); spliced = CA.load_spliced()
resc = CA.load_rescued(meta, skel, spliced)
fams = CA.colocated_families(meta, rescued_fam=resc)
famspans = []
for fid, dom, copies in fams:
    for (tid, c, s1, en, st) in copies: famspans.append((c, s1 - 1, en))
def in_family(chrom, s, e):
    return any(c == chrom and not (e < fs or s > fe) for c, fs, fe in famspans)
new = [f for f in flags if not in_family(f["chrom"], f["start"], f["end"])]
print(f"  of which NEW (outside the 70 multi-copy families): {len(new)}")
json.dump(new, open("/home/juanfra/winloci_scratch/refabsent/genomewide_flags_new.json", "w"), indent=1)
