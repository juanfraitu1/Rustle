#!/usr/bin/env python3
"""Finalize the per-member attribution by joining three signals onto the coarse buckets:
  1. found-leg               (member_attribution.tsv)  -> HOW found (rna-split/projection/expr-collapse/protein)
  2. top-up recovery         (topup_redetect/topup_catalog.copies.tsv) -> does ideal coverage seed a locus?
  3. mis-chain fraction      (residual_mischain.tsv)   -> do the real reads mis-align across giant introns?

Final WHY-NOT buckets for the residual (members with reads but not found), disambiguated:
  MISCHAIN            mischain_frac >= MC  -> reads mis-align across giant introns (NCF1-style); gate drops them.
                      Real-data limitation, fixable by salvaging the local portion of the read.
  K0-FLOOR           low mis-chain AND top-up did NOT seed a locus -> identical to a resolved sibling; even
                      ideal coverage of the member's own transcript maps MAPQ-0 and never seeds. Fundamental.
  GENUINE-BUG        low mis-chain, top-up seeds, AND ample real unique coverage yet still missed -> real miss.
  COVERAGE-LIMITED   low mis-chain, top-up seeds, but thin real coverage -> more depth recovers it.
Plus the already-clear buckets: FOUND:<leg> and MISS:not-expressed (no transcript / 0-template).

Usage: member_attribution_final.py
"""
import os, csv
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
ATTR = f"{HERE}/member_attribution.tsv"
MISCHAIN = f"{HERE}/residual_mischain.tsv"
TOPUP = "/home/juanfra/winloci_scratch/soto_cache/topup_redetect/topup_catalog.copies.tsv"
OUT = f"{HERE}/member_attribution_final.tsv"

MC = 0.30        # mis-chain fraction threshold
TOPUP = "/home/juanfra/winloci_scratch/soto_cache/topup_flank/cache_topup_catalog.copies.tsv"  # corrected (flank)

# topup catalog: chrom -> list of (start,end); member is topup-recovered if a copy overlaps it
topup = defaultdict(list)
if os.path.exists(TOPUP):
    r = csv.reader(open(TOPUP), delimiter="\t"); next(r, None)
    for row in r:
        if len(row) < 6: continue
        try: topup[row[3]].append((int(row[4]), int(row[5])))
        except ValueError: continue
def topup_seeds(chrom, s, e):
    return any(not (a > e or b < s) for a, b in topup.get(chrom, ()))

# mis-chain fraction per member
mcfrac = {}
if os.path.exists(MISCHAIN):
    for l in open(MISCHAIN):
        c = l.rstrip("\n").split("\t")
        if c[0] == "fam" or len(c) < 9: continue
        mcfrac[(c[0], c[1])] = (float(c[8]), int(c[5]))  # (frac, n_prim)

rows = list(csv.DictReader(open(ATTR), delimiter="\t"))
topup_present = bool(topup)

# POWER CHECK: the topup verdict is only trustworthy if the detection recovers the real-found members on
# these chroms. If control recovery is low, the substrate is under-powered and topup_seeds is unreliable.
RES_CHROMS = set("chr1 chr2 chr3 chr5 chr7 chr8 chr9 chr11 chr14 chr15 chr16 chr17 chr21 chr22".split())
ctrl = [r for r in rows if r["bucket"] == "FOUND:rna-split" and r["chrom"] in RES_CHROMS]
ctrl_rec = sum(1 for r in ctrl if topup_seeds(r["chrom"], int(r["start"]), int(r["end"])))
power = ctrl_rec / max(len(ctrl), 1)
topup_trustworthy = topup_present and power >= 0.75
print(f"TOPUP POWER CHECK: control real-found recovery = {ctrl_rec}/{len(ctrl)} = {100*power:.0f}%  "
      f"-> topup verdicts {'TRUSTED' if topup_trustworthy else 'NOT TRUSTED (under-powered)'}\n")

final = []
for r in rows:
    b = r["bucket"]; key = (r["fam"], r["gene"])
    if b.startswith("FOUND") or b == "MISS:not-expressed":
        r["final_bucket"] = b  # already clear
        final.append(r); continue
    # residual: all have unique support (uniq_frac>0.5), so NOT no-unique K=0. Causes: mis-chain / collapse / bug.
    frac, nprim = mcfrac.get(key, (0.0, 0))
    seeds = topup_seeds(r["chrom"], int(r["start"]), int(r["end"])) if topup_trustworthy else None
    famf = r["fam_found"] == "True"
    if frac >= MC:
        fb = "MISS:mischain"                 # reads mis-align across giant introns; gate drops them (fixable)
    elif seeds is True:
        fb = "MISS:coverage-recoverable"     # lacked reads; ideal depth+flank now seeds a locus
    elif seeds is False and famf:
        fb = "MISS:collapse-K0"              # family IS found (sibling); this copy collapses into it (O2 limit)
    elif seeds is False:
        fb = "MISS:genuine-bug"              # clean reads, no sibling, ideal coverage still nothing -> real miss
    else:  # topup untrusted -> fall back to co-location only
        fb = "MISS:collapse-K0" if famf else "MISS:genuine-bug"
    r["mischain_frac"] = frac
    r["topup_seeds"] = seeds
    r["final_bucket"] = fb
    final.append(r)

# write
cols = ["fam", "gene", "chrom", "start", "end", "bp", "cls", "mean_total", "mean_uniq",
        "fam_found", "mischain_frac", "topup_seeds", "bucket", "final_bucket"]
with open(OUT, "w") as f:
    f.write("\t".join(cols) + "\n")
    for r in final:
        f.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")

# summary
N = len(final)
fb = defaultdict(int)
for r in final: fb[r["final_bucket"]] += 1
print(f"=== FINAL ATTRIBUTION ({N} members) ===")
order = ["FOUND:rna-split", "FOUND:projection", "FOUND:expr-collapse", "FOUND:protein-tail",
         "MISS:mischain", "MISS:coverage-recoverable", "MISS:collapse-K0", "MISS:genuine-bug",
         "MISS:not-expressed"]
for k in order:
    if fb.get(k): print(f"  {k:<28}{fb[k]:>4}")
for k in sorted(fb):
    if k not in order: print(f"  {k:<28}{fb[k]:>4}")
found = sum(v for k, v in fb.items() if k.startswith("FOUND"))
not_expr = fb.get("MISS:not-expressed", 0)
k0 = fb.get("MISS:collapse-K0", 0)
detectable = N - not_expr
print(f"\nraw recall:                 {found}/{N} = {100*found/N:.1f}%")
print(f"RNA-detectable universe:    {found}/{detectable} = {100*found/detectable:.1f}%  (drop {not_expr} not-expressed)")
print(f"copy-resolution ceiling:    {found}/{detectable-k0} = {100*found/(detectable-k0):.1f}%  "
      f"(also drop {k0} collapse-K0: family found, copy not resolved)")
for lab, tag in [("genuine-bug (clean reads, ideal coverage still nothing -> REAL MISS; VERIFY)", "MISS:genuine-bug"),
                 ("mis-chain (reads mis-align across giant introns; fixable by read salvage)", "MISS:mischain"),
                 ("coverage-recoverable (lacked reads; ideal depth seeds it)", "MISS:coverage-recoverable"),
                 ("collapse-K0 (family found via sibling; this copy not resolved)", "MISS:collapse-K0")]:
    ms = [r for r in final if r["final_bucket"] == tag]
    if ms:
        print(f"\n{tag} ({len(ms)}) -- {lab}:")
        for r in ms:
            print(f"  {r['fam']:<8} {r['gene']:<14} {r['chrom']:<6} uniq={r['mean_uniq']:>7} mischain={r.get('mischain_frac','')}")
print(f"\nwrote {OUT}")
