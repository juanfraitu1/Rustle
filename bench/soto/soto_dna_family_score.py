#!/usr/bin/env python3
"""Score member detection with the dna_family_fallback re-admissions as a 5th leg, and audit its precision.

The dna_family.tsv from `gw_family_catalog --homology-primary` (RUSTLE_DNA_FAMILY_FALLBACK=1) has columns:
  family_id[0]  chrom[1]  start[2]  end[3]  famCN[4]  min_locus_reads[5]  status[6]  projection_loci[7]
where projection_loci is a `chrom:start-end@id;...` list of the genomic paralogs the orphan locus projects to.
A Soto member is recovered by this leg if EITHER the re-admitted locus (chrom,start,end) OR one of its
projection loci overlaps the member (the projections recover the paralog members — e.g. AC134878.2 projects
to TEKT4P2's locus).

Measures: recall WITHOUT vs WITH the dna_family leg (what the fallback adds), target recovery, and precision
(how many dna_family loci/projections land on NO Soto member = off-benchmark, the over-detection risk).

Usage: soto_dna_family_score.py <homprim.copies.tsv> <homprim.dna_family.tsv>
"""
import csv, sys, os
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
D = "/mnt/linuxdisk/home/juanfraitu/winloci_data"
if not os.path.exists(f"{D}/soto_gw_prot.copies.tsv"):
    D = "/mnt/wsl/PHYSICALDRIVE0p2/home/juanfraitu/winloci_data"  # bind-source fallback when /mnt/linuxdisk is unmounted
BED = f"{HERE}/80_fams.chr.bed"
CAT = sys.argv[1]
DNAFAM = sys.argv[2]

LEGS = [("RNA-split", CAT, 3, 4, 5, None),
        ("protein-tail", f"{D}/soto_gw_prot.copies.tsv", 3, 4, 5, None),
        ("projection", f"{D}/soto_pall.allproj.tsv", 1, 2, 3, None),
        ("expr-collapse", f"{D}/soto_ce.expressed_collapsed.tsv", 1, 2, 3, 7)]
DNA_LEG = ("dna-family", DNAFAM, 1, 2, 3, 7)

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

base_legs = {lab: load(p, cc, sc, ec, pj) for lab, p, cc, sc, ec, pj in LEGS}
dna_idx = load(*DNA_LEG[1:])
for lab in base_legs:
    n = "MISSING" if base_legs[lab] is None else sum(len(v) for v in base_legs[lab].values())
    print(f"  leg {lab:<13}: {n}")
print(f"  leg {'dna-family':<13}: {0 if dna_idx is None else sum(len(v) for v in dna_idx.values())}")

def any_hit(idxs, mc, ms, me):
    for idx in idxs:
        if idx and any(not (a > me or b < ms) for a, b in idx.get(mc, ())):
            return True
    return False

members = []
for l in open(BED):
    c = l.rstrip("\n").split("\t")
    members.append((c[3].split("|")[-1], c[3].split("|")[0], c[0], int(c[1]), int(c[2])))
N = len(members)
base = [v for v in base_legs.values()]
without = sum(1 for f, g, mc, ms, me in members if any_hit(base, mc, ms, me))
withdna = sum(1 for f, g, mc, ms, me in members if any_hit(base + [dna_idx], mc, ms, me))
print(f"\nrecall WITHOUT dna-family: {without}/{N} = {100*without/N:.1f}%")
print(f"recall WITH    dna-family: {withdna}/{N} = {100*withdna/N:.1f}%   (dna-family adds +{withdna-without})")

# which members does dna-family newly recover?
newly = [(f, g, mc) for f, g, mc, ms, me in members
         if any_hit([dna_idx], mc, ms, me) and not any_hit(base, mc, ms, me)]
print(f"\nmembers newly recovered by dna-family ({len(newly)}):")
for f, g, mc in newly: print(f"  {f:<8} {g:<14} {mc}")

# target recovery
tgt = set()
for l in open(f"{HERE}/member_attribution.final.tsv"):
    x = l.split("\t")
    if len(x) > 12 and x[12] in ("MISS:mischain", "MISS:genuine-miss", "MISS:collapse-K0"):
        tgt.add((x[0], x[1]))
tg = sum(1 for f, g, mc in newly if (f, g) in tgt)
print(f"  of which {tg} are mis-chain/genuine-miss/collapse-K0 TARGETS")

# PRECISION: dna_family loci (primary + projections) that overlap NO Soto member = off-benchmark
def overlaps_any_member(ch, s, e):
    for f, g, mc, ms, me in members:
        if mc == ch and not (ms > e or me < s): return True
    return False
onm = offm = 0
r = csv.reader(open(DNAFAM), delimiter="\t"); next(r, None)
for row in r:
    if len(row) < 8: continue
    locs = [(row[1], int(row[2]), int(row[3]))]
    for seg in row[7].split(";"):
        if not seg: continue
        try:
            ch, rng = seg.split("@")[0].split(":"); s, e = rng.split("-"); locs.append((ch, int(s), int(e)))
        except ValueError: continue
    if any(overlaps_any_member(c, s, e) for c, s, e in locs): onm += 1
    else: offm += 1
print(f"\nPRECISION of dna_family re-admissions: {onm} on-member, {offm} off-benchmark (over-detection risk)")
print(f"  on-member rate: {100*onm/max(onm+offm,1):.0f}%")
