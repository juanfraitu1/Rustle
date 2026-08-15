#!/usr/bin/env python3
"""Can the MULTIMAP FOOTPRINT define locus extent better than the annotation does?

THE PROPOSAL (user's, and it is NOT the conflict graph): use multimapping as a PER-BASE property of a
single locus -- "is this base covered by reads that also map elsewhere" -- rather than as a pairwise
relation between loci. A per-base mappability track never names the other locus, so no E_c information
enters O1. Node extent, not edges.

THE FALSIFIABLE PREDICTION: NPIP's duplicated unit is a SIZE-INVARIANT ~16 kb cassette while annotated
spans vary 4.65x. If the multimap footprint tracks the real duplicated unit, its length must be
CONSISTENT ACROSS COPIES -- i.e. its coefficient of variation must beat annotation's. If the footprint's
CV is no better, the proposal is dead and we say so.

⚠⚠ THE SIGNAL IS SECONDARY-ALIGNMENT DEPTH, NOT MAPQ=0. First attempt used MAPQ=0 on primaries and got
  a zero-length footprint at 13/19 loci. Diagnosis: at NPIPB11, 234 of 271 PRIMARIES carry MAPQ 60 --
  minimap2 breaks the tie confidently -- while the same window holds 25,287 SECONDARY alignments. The
  standing "-F 2308 before any per-read statistic" rule guards per-read CIGAR/clip stats; here the
  secondary alignments ARE the object of measurement, so they must be kept.
⚠ Still per-locus: a secondary alignment says "a read placed elsewhere also fits HERE". It never names
  where, so no pairwise E_c information enters. Node extent, not edges.
⚠ Reference-consuming span includes N (introns) -- a locus spans its introns. Read-through chimeras are
  guarded by MAX_SPAN and reported separately, never silently dropped.

usage: multimap_extent.py BAM GENE_BED [FLANK] [out.tsv]
"""
import subprocess
import sys
from collections import defaultdict

bam, gene_bed = sys.argv[1], sys.argv[2]
FLANK = int(sys.argv[3]) if len(sys.argv) > 3 else 25000
out_tsv = sys.argv[4] if len(sys.argv) > 4 else None

# ⚠⚠ AN ABSOLUTE DEPTH CUTOFF DOES NOT EXIST. Measured over a 625 kb window at NPIPB11, secondary depth
#    spans 8 to 21,670 -- background is 8-60, the duplicated block is 1,000-21,670. A cutoff of 5 sat
#    BELOW background and saturated the whole window, producing footprints of exactly (annotated + 2*FLANK)
#    and a CV of 0.134 that measured the window size, not the genome. Same shape as the absolute-block-floor
#    result: the hidden oracle is the assumption that one number exists.
# ⟹ threshold is RELATIVE to each locus's own profile, and saturation is flagged, never hidden.
MIN_DEPTH = 5        # absolute floor, only to exclude the true background
FRAC_OF_MAX = 0.05   # a base counts if it also reaches 5% of THIS locus's peak secondary depth
MAX_SPAN = 200000    # alignments longer than this are read-through, counted separately
REF_CONSUMING = {"M", "=", "X", "D", "N"}


def ref_span(cigar):
    n, tot = "", 0
    for ch in cigar:
        if ch.isdigit():
            n += ch
        else:
            if ch in REF_CONSUMING:
                tot += int(n or 0)
            n = ""
    return tot


genes = []
for line in open(gene_bed):
    f = line.rstrip("\n").split("\t")
    if len(f) >= 4:
        genes.append((f[0], int(f[1]), int(f[2]), f[3]))

rows = []
skipped_total = 0
for c, gs, ge, name in genes:
    lo, hi = max(0, gs - FLANK), ge + FLANK
    # keep secondary (256); drop only unmapped (4) and supplementary (2048)
    p = subprocess.run(["samtools", "view", "-F", "2052", bam, f"{c}:{lo+1}-{hi}"],
                       capture_output=True, text=True)
    L = hi - lo
    d_sec = [0] * (L + 1)   # difference arrays -- O(alignments + L), not O(alignments * span)
    skipped = 0
    for line in p.stdout.splitlines():
        f = line.split("\t", 9)
        if len(f) < 6:
            continue
        flag, pos, cig = int(f[1]), int(f[3]) - 1, f[5]
        if cig == "*" or not (flag & 256):
            continue
        sp = ref_span(cig)
        if sp > MAX_SPAN:
            skipped += 1
            continue
        a, b = max(pos, lo) - lo, min(pos + sp, hi) - lo
        if b <= a:
            continue
        d_sec[a] += 1
        d_sec[b] -= 1
    skipped_total += skipped

    sec, run = [0] * L, 0
    for i in range(L):
        run += d_sec[i]
        sec[i] = run

    # footprint = longest run above a threshold RELATIVE to this locus's own peak
    thr = max(MIN_DEPTH, FRAC_OF_MAX * max(sec)) if sec else MIN_DEPTH
    best_len = best_a = best_b = 0
    cur = None
    for i in range(L):
        ok = sec[i] >= thr
        if ok and cur is None:
            cur = i
        elif not ok and cur is not None:
            if i - cur > best_len:
                best_len, best_a, best_b = i - cur, cur, i
            cur = None
    if cur is not None and L - cur > best_len:
        best_len, best_a, best_b = L - cur, cur, L

    ann = ge - gs
    sat = (best_a == 0) or (best_b == L)   # touching a window edge = the window, not the genome
    rows.append((name, c, ann, best_len, lo + best_a, lo + best_b, skipped, sat))

print(f"multimap footprint vs annotated span   (flank {FLANK//1000} kb, "
      f"min depth {MIN_DEPTH}, secondary-alignment depth)")
print(f"\n{'gene':<10}{'annotated':>11}{'footprint':>11}{'ratio':>8}   footprint interval")
for name, c, ann, fp, a, b, sk, sat in rows:
    r = fp / ann if ann else 0
    print(f"{name:<10}{ann:>11,}{fp:>11,}{r:>8.2f}   {c}:{a+1:,}-{b:,}"
          + ("   ⚠SATURATED" if sat else "")
          + (f"   [{sk} rt]" if sk else ""))


def cv(v):
    v = [x for x in v if x > 0]
    if len(v) < 2:
        return float("nan")
    m = sum(v) / len(v)
    return (sum((x - m) ** 2 for x in v) / (len(v) - 1)) ** 0.5 / m


ann_v = [r[2] for r in rows]
fp_v = [r[3] for r in rows]
print(f"\n{'':<14}{'n':>4}{'median':>10}{'min':>10}{'max':>10}{'spread':>9}{'CV':>8}")
for lbl, v in (("annotated span", ann_v), ("multimap fp", fp_v)):
    w = sorted(x for x in v if x > 0)
    if not w:
        print(f"{lbl:<14}{0:>4}   (none)")
        continue
    print(f"{lbl:<14}{len(w):>4}{w[len(w)//2]:>10,}{w[0]:>10,}{w[-1]:>10,}"
          f"{w[-1]/w[0]:>8.2f}x{cv(v):>8.3f}")

print(f"\n  PREDICTION: footprint CV must BEAT annotation CV, else the proposal is dead.")
c_a, c_f = cv(ann_v), cv(fp_v)
print(f"  annotation CV {c_a:.3f}   footprint CV {c_f:.3f}   -> "
      + ("FOOTPRINT WINS" if c_f < c_a else "FOOTPRINT DOES NOT WIN"))
nsat = sum(1 for r in rows if r[7])
print(f"  saturated (footprint touches a window edge -> measures the WINDOW, not the genome): {nsat}/{len(rows)}")
if nsat:
    print("  ⚠ CV above is NOT trustworthy while any locus is saturated -- widen FLANK.")
if skipped_total:
    print(f"  ⚠ {skipped_total} alignments longer than {MAX_SPAN:,} bp excluded as read-through.")

if out_tsv:
    with open(out_tsv, "w") as fh:
        fh.write("gene\tchrom\tannotated_bp\tfootprint_bp\tratio\tfp_start\tfp_end\treadthrough\tsaturated\n")
        for name, c, ann, fp, a, b, sk, sat in rows:
            fh.write(f"{name}\t{c}\t{ann}\t{fp}\t{fp/ann if ann else 0:.3f}\t{a}\t{b}\t{sk}\t{int(sat)}\n")
    print(f"\nwrote {out_tsv}")
