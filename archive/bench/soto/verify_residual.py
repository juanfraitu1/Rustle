#!/usr/bin/env python3
"""Deterministic verification of the residual buckets, WITHOUT the (compromised) top-up test.
Reliable signals only:
  * found-leg          -> FOUND (306) : 4-leg overlap
  * mean coverage      -> not-expressed (22): mean primary depth < 3
  * mis-chain fraction -> mischain (12): real reads carry giant introns
For the remaining residual (unique reads, low mis-chain) split by CATALOG geometry, not top-up:
  * nearest resolved copy (any leg) on the same chrom, and whether the member's FAMILY has any found copy.
    - a found copy overlaps/abuts the member (<=2kb gap)   -> RESOLVED-ELSEWHERE (coordinate/scoring near-miss)
    - family HAS a found member, member itself unresolved  -> COLLAPSE-K0 (family found; copy not separated)
    - family has NO found member, clean unique reads       -> GENUINE-MISS (verify: real bug)
Emits the final per-member TSV and the verified genuine-miss shortlist.
"""
import csv, os
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = "/home/juanfra/winloci_scratch/soto_cache"
D = "/mnt/linuxdisk/home/juanfraitu/winloci_data"
FINAL = f"{HERE}/member_attribution_final.tsv"
OUT = f"{HERE}/member_attribution_verified.tsv"

# all detection legs -> resolved copies per chrom (used to find the nearest resolved copy to a member)
LEGS = [(f"{CACHE}/definitive.copies.tsv", 3, 4, 5),
        (f"{D}/soto_gw_prot.copies.tsv", 3, 4, 5),
        (f"{D}/soto_pall.allproj.tsv", 1, 2, 3),
        (f"{D}/soto_ce.expressed_collapsed.tsv", 1, 2, 3)]
copies = defaultdict(list)
for path, cc, sc, ec in LEGS:
    try:
        r = csv.reader(open(path), delimiter="\t"); next(r, None)
        for row in r:
            if len(row) <= max(cc, sc, ec): continue
            try: copies[row[cc]].append((int(row[sc]), int(row[ec])))
            except ValueError: continue
    except FileNotFoundError: pass

def nearest_copy_bp(chrom, s, e):
    best = None
    for a, b in copies.get(chrom, ()):
        gap = 0 if not (a > e or b < s) else min(abs(a - e), abs(b - s))
        if best is None or gap < best: best = gap
    return best  # None if no copy on chrom; 0 if overlaps

rows = list(csv.DictReader(open(FINAL), delimiter="\t"))
# family -> list of members with found flag
fam_members = defaultdict(list)
for r in rows:
    fam_members[r["fam"]].append(r)

final = []
for r in rows:
    b = r["final_bucket"]
    if b in ("MISS:mischain",) or b.startswith("FOUND") or b == "MISS:not-expressed":
        r["verified_bucket"] = b
        r["nearest_copy_bp"] = nearest_copy_bp(r["chrom"], int(r["start"]), int(r["end"]))
        final.append(r); continue
    # residual collapse-K0 / genuine-bug / coverage: re-derive from catalog geometry
    nc = nearest_copy_bp(r["chrom"], int(r["start"]), int(r["end"]))
    fam_found = any(m["final_bucket"].startswith("FOUND") for m in fam_members[r["fam"]])
    if nc is not None and nc <= 2000:
        vb = "MISS:resolved-elsewhere"   # a resolved copy overlaps/abuts -> coordinate near-miss, not a true miss
    elif fam_found:
        vb = "MISS:collapse-K0"          # family detected via a sibling; this copy not individually resolved
    else:
        vb = "MISS:genuine-miss"         # family entirely missed, clean unique reads -> real bug to verify
    r["verified_bucket"] = vb
    r["nearest_copy_bp"] = nc
    final.append(r)

cols = ["fam", "gene", "chrom", "start", "end", "bp", "cls", "mean_total", "mean_uniq",
        "fam_found", "mischain_frac", "nearest_copy_bp", "final_bucket", "verified_bucket"]
with open(OUT, "w") as f:
    f.write("\t".join(cols) + "\n")
    for r in final:
        f.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")

N = len(final); fb = defaultdict(int)
for r in final: fb[r["verified_bucket"]] += 1
print("=== VERIFIED ATTRIBUTION (top-up discarded; catalog-geometry + reliable signals) ===")
order = ["FOUND:rna-split", "FOUND:projection", "FOUND:expr-collapse", "FOUND:protein-tail",
         "MISS:mischain", "MISS:resolved-elsewhere", "MISS:collapse-K0", "MISS:genuine-miss", "MISS:not-expressed"]
for k in order:
    if fb.get(k): print(f"  {k:<26}{fb[k]:>4}")
found = sum(v for k, v in fb.items() if k.startswith("FOUND"))
ne = fb.get("MISS:not-expressed", 0); k0 = fb.get("MISS:collapse-K0", 0); re_ = fb.get("MISS:resolved-elsewhere", 0)
print(f"\nraw recall:               {found}/{N} = {100*found/N:.1f}%")
print(f"+ resolved-elsewhere:     {found+re_}/{N} = {100*(found+re_)/N:.1f}%  (near-miss coords count)")
print(f"RNA-detectable universe:  {found+re_}/{N-ne} = {100*(found+re_)/(N-ne):.1f}%  (drop {ne} not-expressed)")
print(f"copy-resolution ceiling:  {found+re_}/{N-ne-k0} = {100*(found+re_)/(N-ne-k0):.1f}%  (also drop {k0} collapse-K0)")
print("\ngenuine-miss shortlist (family entirely missed, clean unique reads -> VERIFY as real bugs):")
for r in final:
    if r["verified_bucket"] == "MISS:genuine-miss":
        sibs = [m['gene'] for m in fam_members[r['fam']]]
        print(f"  {r['fam']:<8} {r['gene']:<14} {r['chrom']:<6} uniq={r['mean_uniq']:>7} nearest_copy={r['nearest_copy_bp']}bp  fam={sibs}")
print(f"\nwrote {OUT}")
