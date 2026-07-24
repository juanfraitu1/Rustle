#!/usr/bin/env python3
"""Score a rebuilt RNA-split catalog against the Soto members (member-detection recall/precision), keeping
the other 3 detection legs. Prints recall vs the committed baseline + which members changed. Used by
recompute_soto.sh for the fast cached-BAM validation loop.

Usage: soto_cache_score.py <catalog.copies.tsv>
"""
import csv, sys, os
from collections import defaultdict
CAT = sys.argv[1] if len(sys.argv) > 1 else "/home/juanfra/winloci_scratch/soto_cache/catalog.copies.tsv"
HERE = os.path.dirname(os.path.abspath(__file__))
D = "/mnt/linuxdisk/home/juanfraitu/winloci_data"
BED = f"{HERE}/80_fams.chr.bed"
MDET = f"{HERE}/soto_member_detection.tsv"   # committed baseline (76.2%)
# 4 legs: RNA-split = the freshly rebuilt catalog; the other 3 unchanged (genome-wide). The 6th tuple
# element is the projection-loci column (semicolon `chrom:start-end@id` list of extra K0-collapsed loci);
# expr-collapse (col 7) MUST be expanded exactly like soto_detection_eval.load_loci, else its members are
# under-counted (they appear only via the projected loci, not the primary span).
LEGS = [("RNA-split", CAT, 3, 4, 5, None),
        ("protein-tail", f"{D}/soto_gw_prot.copies.tsv", 3, 4, 5, None),
        ("projection", f"{D}/soto_pall.allproj.tsv", 1, 2, 3, None),
        ("expr-collapse", f"{D}/soto_ce.expressed_collapsed.tsv", 1, 2, 3, 7)]

def load(path, cc, sc, ec, proj_col):
    idx = defaultdict(list)
    try:
        r = csv.reader(open(path), delimiter="\t"); next(r, None)
        for row in r:
            if len(row) <= max(cc, sc, ec): continue
            try: idx[row[cc]].append((int(row[sc]), int(row[ec])))
            except ValueError: continue
            if proj_col is not None and len(row) > proj_col:
                for seg in row[proj_col].split(";"):
                    if not seg: continue
                    try:
                        ch, rng = seg.split("@")[0].split(":"); s, e = rng.split("-")
                        idx[ch].append((int(s), int(e)))
                    except ValueError: continue
    except FileNotFoundError:
        return None
    return idx

legs = {lab: load(p, cc, sc, ec, pj) for lab, p, cc, sc, ec, pj in LEGS}
for lab, p, *_ in LEGS:
    n = "MISSING" if legs[lab] is None else f"{sum(len(v) for v in legs[lab].values())} copies"
    print(f"    leg {lab:<13}: {n}")

def hit(mc, ms, me):
    for lab, idx in legs.items():
        if idx and any(not (a > me or b < ms) for a, b in idx.get(mc, ())):
            return lab
    return None

members = []
for l in open(BED):
    c = l.rstrip("\n").split("\t")
    members.append((c[3].split("|")[0], c[3].split("|")[-1], c[0], int(c[1]), int(c[2])))
oldY = set()
for l in open(MDET):
    if l.startswith("family_id"): continue
    c = l.rstrip("\n").split("\t")
    if c[5] == "Y": oldY.add((c[0], c[1]))
newY = set(); newly = []; lost = []
for g, fid, mc, ms, me in members:
    lab = hit(mc, ms, me)
    if lab: newY.add((fid, g))
    if lab and (fid, g) not in oldY: newly.append((fid, g, lab))
    if (not lab) and (fid, g) in oldY: lost.append((fid, g))
N = len(members)
print(f"\n  members: {N}   OLD (committed): {len(oldY)}/{N} = {100*len(oldY)/N:.1f}%   "
      f"NEW: {len(newY)}/{N} = {100*len(newY)/N:.1f}%  (delta {len(newY)-len(oldY):+d})")
if newly: print(f"  newly detected ({len(newly)}): {[f'{f}/{g}({l})' for f,g,l in newly]}")
if lost:  print(f"  LOST vs baseline ({len(lost)}): {[f'{f}/{g}' for f,g in lost]}")
print("  VALIDATION: recall near 76.2% => cache+recipe reproduce the genome-wide result; far off => not comparable.")
