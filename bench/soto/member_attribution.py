#!/usr/bin/env python3
"""Attribute EVERY Soto member: HOW it was found (which detection leg) or, if not, WHY not.

Motivation: Soto members are DNA segmental duplications overlapping >=1 exon, so the catalog INCLUDES
pseudogenes and gene fragments that need not be transcribed. An RNA method cannot detect an unexpressed
DNA feature, so those misses are outside the detectable universe -- not method failures. This script
attributes all 362 members so recall can be reported on the RNA-detectable denominator.

Signals joined per member (all by bed row index / coordinate, no name collisions):
  * found_leg  -- which of the 4 scorer legs overlaps it (RNA-split=definitive catalog, projection,
                  expr-collapse, protein-tail), or NONE.  = HOW found.
  * cls        -- annotation class from artifacts/member_annotation.json (coding / transcribed_pseudogene
                  / pseudogene / piece-of-other / lncRNA / fragment / unannotated).
  * mean_total -- primary (-F 2308) depth-sum / bp  = expression incl. multimapped spillover.
  * mean_uniq  -- MAPQ>=1 primary depth-sum / bp    = UNIQUE support (reads that prefer THIS copy).
  * fam_found  -- does any other member of the same family get found (co-location / paralog-resolved).

Buckets for a NOT-FOUND member (why not):
  not-expressed        mean_total < MT_EXPR                      -> no transcript here; DNA-only Soto leftover
  K0-no-unique         expressed, mean_uniq ~ 0, fam has a found -> paralog tie, no unique read: identifiability
  low-unique-orphan    mean_uniq ~ 0, fam has NO found          -> only multimapped/spillover, likely silent
  coverage-limited     some unique support but below floor      -> more depth would recover
  genuine-investigate  substantial unique support yet missed    -> real miss / mis-chain: needs adjudication

Thresholds are printed with the signal distributions so they can be calibrated to the data, not guessed.
Usage: member_attribution.py [catalog.copies.tsv]   (default = cache definitive.copies.tsv)
"""
import csv, sys, os, json
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = "/home/juanfra/winloci_scratch/soto_cache"
D = "/mnt/linuxdisk/home/juanfraitu/winloci_data"
BED = f"{HERE}/80_fams.chr.bed"
ANN = f"{HERE}/artifacts/member_annotation.json"
CAT = sys.argv[1] if len(sys.argv) > 1 else f"{CACHE}/definitive.copies.tsv"
COV_TOTAL = "/tmp/cov_total.tsv"   # bedcov -Q0 -G2308 : name<TAB>depthsum, bed order
COV_UNIQ = "/tmp/cov_uniq.tsv"     # bedcov -Q1 -G2308 : name<TAB>depthsum, bed order

# thresholds (mean per-base primary depth). Calibrated below from the printed distributions.
MT_EXPR = float(os.environ.get("MT_EXPR", 3.0))    # below this mean_total -> not expressed
MU_UNIQ = float(os.environ.get("MU_UNIQ", 1.0))    # below this mean_uniq  -> no real unique support
COV_FLOOR = float(os.environ.get("COV_FLOOR", 8.0))  # unique support present but below this -> coverage-limited

LEGS = [("rna-split", CAT, 3, 4, 5, None),
        ("protein-tail", f"{D}/soto_gw_prot.copies.tsv", 3, 4, 5, None),
        ("projection", f"{D}/soto_pall.allproj.tsv", 1, 2, 3, None),
        ("expr-collapse", f"{D}/soto_ce.expressed_collapsed.tsv", 1, 2, 3, 7)]

def load_leg(path, cc, sc, ec, proj_col):
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

legs = {lab: load_leg(p, cc, sc, ec, pj) for lab, p, cc, sc, ec, pj in LEGS}
for lab, p, *_ in LEGS:
    n = "MISSING" if legs[lab] is None else f"{sum(len(v) for v in legs[lab].values())}"
    print(f"  leg {lab:<13}: {n}", file=sys.stderr)

def found_leg(mc, ms, me):
    # priority order: RNA de-novo first, then copy-number projections, then protein
    for lab in ("rna-split", "projection", "expr-collapse", "protein-tail"):
        idx = legs.get(lab)
        if idx and any(not (a > me or b < ms) for a, b in idx.get(mc, ())):
            return lab
    return None

# annotation class by coordinate
ann = {}
for x in json.load(open(ANN)):
    ann[(x["chrom"], x["start"], x["end"])] = x

# coverage by bed row index
def col2(path):
    return [float(l.split("\t")[1]) for l in open(path)]
tot = col2(COV_TOTAL); uniq = col2(COV_UNIQ)

# members from bed (authoritative), in order
members = []
for l in open(BED):
    c = l.rstrip("\n").split("\t")
    gene, fam = c[3].split("|")[0], c[3].split("|")[-1]
    members.append({"gene": gene, "fam": fam, "chrom": c[0], "start": int(c[1]), "end": int(c[2])})

# pass 1: found status + which members' families have any found member
fam_has_found = defaultdict(bool)
for m in members:
    m["leg"] = found_leg(m["chrom"], m["start"], m["end"])
    if m["leg"]:
        fam_has_found[m["fam"]] = True

# pass 2: attach coverage, class, bucket
rows = []
for i, m in enumerate(members):
    bp = m["end"] - m["start"]
    mt = tot[i] / bp if bp else 0.0
    mu = uniq[i] / bp if bp else 0.0
    a = ann.get((m["chrom"], m["start"], m["end"]), {})
    cls = a.get("cls", "?")
    leg = m["leg"]
    if leg:
        bucket = f"FOUND:{leg}"
    else:
        famf = fam_has_found[m["fam"]]
        if mt < MT_EXPR:
            bucket = "MISS:not-expressed"
        elif mu < MU_UNIQ and famf:
            bucket = "MISS:K0-no-unique"
        elif mu < MU_UNIQ and not famf:
            bucket = "MISS:low-unique-orphan"
        elif mu < COV_FLOOR:
            bucket = "MISS:coverage-limited"
        else:
            bucket = "MISS:genuine-investigate"
    rows.append({**m, "bp": bp, "cls": cls, "mean_total": round(mt, 1),
                 "mean_uniq": round(mu, 1), "fam_found": fam_has_found[m["fam"]], "bucket": bucket})

# write per-member TSV
out = f"{HERE}/member_attribution.tsv"
cols = ["fam", "gene", "chrom", "start", "end", "bp", "cls", "mean_total", "mean_uniq", "fam_found", "bucket"]
with open(out, "w") as f:
    f.write("\t".join(cols) + "\n")
    for r in rows:
        f.write("\t".join(str(r[c]) for c in cols) + "\n")

# ---- summary ----
N = len(rows)
found = [r for r in rows if r["bucket"].startswith("FOUND")]
miss = [r for r in rows if r["bucket"].startswith("MISS")]
print(f"\n=== ATTRIBUTION SUMMARY ({N} members, catalog={os.path.basename(CAT)}) ===")
print(f"thresholds: MT_EXPR={MT_EXPR} MU_UNIQ={MU_UNIQ} COV_FLOOR={COV_FLOOR} (mean per-base primary depth)")

bc = defaultdict(int)
for r in rows: bc[r["bucket"]] += 1
print("\nbucket breakdown:")
for b, n in sorted(bc.items(), key=lambda kv: (not kv[0].startswith("FOUND"), -kv[1])):
    print(f"  {b:<26}{n:>4}")

# reframed recall
not_expr = sum(1 for r in miss if r["bucket"] == "MISS:not-expressed")
orphan = sum(1 for r in miss if r["bucket"] == "MISS:low-unique-orphan")
detectable = N - not_expr - orphan  # RNA-detectable universe = drop the silent DNA-only members
print(f"\nrecall (raw):                 {len(found)}/{N} = {100*len(found)/N:.1f}%")
print(f"not-expressed (DNA-only):     {not_expr}   low-unique-orphan (likely silent): {orphan}")
print(f"recall (RNA-detectable univ): {len(found)}/{detectable} = {100*len(found)/detectable:.1f}%")

# distributions to calibrate thresholds
print("\nmean_total distribution among MISSES (to calibrate MT_EXPR):")
mts = sorted(r["mean_total"] for r in miss)
if mts:
    q = lambda p: mts[min(len(mts)-1, int(p*len(mts)))]
    print(f"  n={len(mts)} min={mts[0]} p25={q(.25)} med={q(.5)} p75={q(.75)} max={mts[-1]}")
print("mean_uniq distribution among MISSES (to calibrate MU_UNIQ/COV_FLOOR):")
mus = sorted(r["mean_uniq"] for r in miss)
if mus:
    q = lambda p: mus[min(len(mus)-1, int(p*len(mus)))]
    print(f"  n={len(mus)} min={mus[0]} p25={q(.25)} med={q(.5)} p75={q(.75)} max={mus[-1]}")

# class x found (current catalog)
print("\nannotation class x detection (current catalog):")
byc = defaultdict(lambda: [0, 0])
for r in rows:
    byc[r["cls"]][0] += 1
    if r["bucket"].startswith("FOUND"): byc[r["cls"]][1] += 1
for c, (n, d) in sorted(byc.items(), key=lambda kv: -kv[1][0]):
    print(f"  {c:<24}{d:>3}/{n:<3} = {100*d/n:>3.0f}%")

print(f"\nwrote {out}")
