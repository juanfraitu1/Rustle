#!/usr/bin/env python3
"""o4_fp_bound.py — a specificity / false-positive BOUND for the reference-absent hidden-copy detector
(closes L15 of LOOSE_ENDS_AUDIT.md). No BAM re-scan: intersects the ALREADY-COMPUTED genome-wide flags
(genomewide_flags.json, from hidden_copy_scan_genomewide.py) with a SINGLE-COPY control gene set.

Control = de-novo expressed loci overlapping a protein_coding gene whose NAME is UNIQUE in the RefSeq
annotation (single-copy by construction) and NOT LOC* (uncharacterized/paralog-prone) and NOT overlapping
any MULTI-copy (name appears >1) gene. Such a locus cannot host a real reference-absent paralog COPY, so
any "hidden second haplotype" flag there is a FALSE POSITIVE for the O4 claim (it is a het / RNA-edit /
alignment artifact). The flag rate on the control = an UPPER BOUND on the detector's FP rate.

Contrast: the flag rate on multi-copy-family loci (where real collapsed copies live) should be HIGHER.
"""
import json
import bisect
from collections import defaultdict

SCRATCH = "/home/juanfra/winloci_scratch"
GENES = "/tmp/ggo_genes.gff"

# --- parse genes: name -> list of (chrom,start,end,biotype) ---
by_name = defaultdict(list)
for line in open(GENES):
    f = line.rstrip("\n").split("\t")
    if len(f) < 9:
        continue
    chrom, s, e, attr = f[0], int(f[3]), int(f[4]), f[8]
    name = bt = None
    for kv in attr.split(";"):
        if kv.startswith("Name="):
            name = kv[5:]
        elif kv.startswith("gene_biotype="):
            bt = kv[13:]
    if name:
        by_name[name].append((chrom, s, e, bt))

# single-copy = unique-name, protein_coding, non-LOC ; multi-copy = name appears >1
single_genes = defaultdict(list)  # chrom -> [(start,end)]
multi_iv = defaultdict(list)
n_single = n_multi = 0
for name, locs in by_name.items():
    if len(locs) > 1:
        for c, s, e, bt in locs:
            multi_iv[c].append((s, e))
        n_multi += 1
    else:
        c, s, e, bt = locs[0]
        if bt == "protein_coding" and not name.startswith("LOC"):
            single_genes[c].append((s, e))
            n_single += 1
for c in single_genes:
    single_genes[c].sort()
for c in multi_iv:
    multi_iv[c].sort()

def overlaps_any(iv_list, s, e):
    if not iv_list:
        return False
    starts = [x[0] for x in iv_list]
    i = bisect.bisect_right(starts, e)
    for j in range(max(0, i - 200), i):       # scan back a bounded window of candidates
        a, b = iv_list[j]
        if a <= e and b >= s:
            return True
    return False

# --- flags + all scanned loci ---
loci = json.load(open(f"{SCRATCH}/refabsent/genomewide_loci.json"))
flags = json.load(open(f"{SCRATCH}/refabsent/genomewide_flags.json"))
flagged = {(f["chrom"], f["start"], f["end"]) for f in flags}

# --- classify every scanned locus ---
ctrl_total = ctrl_flagged = 0
multi_total = multi_flagged = 0
for entry in loci:
    if isinstance(entry, (list, tuple)):
        c, s, e = entry[0], int(entry[1]), int(entry[2])
    else:
        c, s, e = entry["chrom"], int(entry["start"]), int(entry["end"])
    is_flag = (c, s, e) in flagged
    in_multi = overlaps_any(multi_iv.get(c, []), s, e)
    in_single = overlaps_any(single_genes.get(c, []), s, e)
    if in_single and not in_multi:
        ctrl_total += 1
        ctrl_flagged += is_flag
    if in_multi:
        multi_total += 1
        multi_flagged += is_flag

def rate(a, b):
    return f"{a}/{b} = {100*a/b:.2f}%" if b else "n/a"

print("=== O4 hidden-copy detector — false-positive BOUND (no BAM re-scan) ===\n")
print(f"genes: {n_single:,} single-copy (unique-name protein_coding, non-LOC) | {n_multi:,} multi-copy (name appears >1)")
print(f"scanned de-novo loci: {len(loci):,} | total flagged genome-wide: {len(flags):,} "
      f"({100*len(flags)/len(loci):.2f}%)\n")
print(f"FP BOUND  (flag rate on SINGLE-COPY control loci): {rate(ctrl_flagged, ctrl_total)}")
print(f"  -> on genes that CANNOT host a reference-absent copy, this fraction still trips the detector")
print(f"     (het / RNA-edit / alignment artifact) = an upper bound on the O4 false-positive rate.\n")
print(f"CONTRAST (flag rate on MULTI-COPY-family loci): {rate(multi_flagged, multi_total)}")
print(f"  -> where real collapsed copies live; should be HIGHER than the control (enrichment = signal).")
enr = (multi_flagged / multi_total) / (ctrl_flagged / ctrl_total) if ctrl_flagged and multi_total else float('nan')
print(f"\nenrichment (multi / control flag rate): {enr:.1f}x")
