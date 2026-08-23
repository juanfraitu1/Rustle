#!/usr/bin/env python3
"""Supersedes o1_classify_absent.py: classify SD-predicted catalog-absent intervals on rep EXONS.

WHY THIS REPLACES THE SPAN VERSION. The earlier classification asked "does a rep's genomic SPAN cover
this interval". Rep spans are ~90.8% intron by bp, so "the interval sits inside a gene's intron"
(a CORRECT rejection -- there is no transcript there to miss) was indistinguishable from "real
transcript sequence that no node covers" (a genuine miss). That ambiguity invalidated the 2026-08-21
audit's expression angle. `nodes.tsv` now carries per-rep exon blocks, so the two separate.

TWO FIXES TO THE OLD SCRIPT, both from the audit:
  * EXACT overlap (max-end prefix sweep) instead of a fixed 800-row backscan window. The window never
    actually fired (max required backscan was 4), but it could have silently inflated "never a node".
  * The stub stratum is reported SEPARATELY and marked INDETERMINATE rather than pooled. A single-exon
    rep's exon array is one block == its whole span, so `exon_bp == span` and the array carries ZERO
    information. 33.07% of reps are stubs; pooling them would claim precision that does not exist on a
    third of the data.

⚠ WHAT THIS DOES NOT ANSWER. The rate of "never a node" is NOT a gap: measured against matched
comparators it IS the base rate for catalog-absent SD sequence (0.9086 over 73,324 non-copy SD sides).
This script asks only the conditional question -- GIVEN a rep is there, is the interval exonic or
intronic -- which is not affected by that base rate because it conditions on rep presence.
⚠ n = 3,928 is not an independent sample size (78.77% self-overlap); the merged-block count is also
reported and every rate should be read at both units.
"""
import bisect, collections, csv, sys

R = "/mnt/linuxdisk/home/juanfraitu/o1_reps2"
G = "/mnt/linuxdisk/home/juanfraitu/o1_gw"
SD = "/tmp/claude-1000/-mnt-c-Users-jfris-Desktop/931c208e-8acb-4dd2-aacb-cf92d5ad051f/scratchpad/sd_absent_intervals.bed"


def load_reps():
    by = collections.defaultdict(list)
    for r in csv.DictReader(open(f"{R}/dump/ggo.nodes.tsv"), delimiter="\t"):
        ex = [tuple(map(int, b.split("-"))) for b in r["exons"].split(",") if b]
        by[r["chrom"]].append((int(r["start"]), int(r["end"]), int(r["n_exon"]), int(r["degree"]), ex))
    idx = {}
    for c, v in by.items():
        v.sort()
        starts = [x[0] for x in v]
        mx, run = [], 0                      # max end over v[0..i] -> exact backward-scan terminator
        for x in v:
            run = max(run, x[1]); mx.append(run)
        idx[c] = (v, starts, mx)
    return idx


def overlappers(idx, c, s, e):
    """EXACT: every rep whose span intersects [s,e). No window heuristic."""
    if c not in idx:
        return []
    v, starts, mx = idx[c]
    hi = bisect.bisect_left(starts, e)       # reps starting at/after e cannot overlap
    out = []
    for j in range(hi - 1, -1, -1):
        if mx[j] <= s:                       # no rep in v[0..j] reaches s -> provably done
            break
        if v[j][1] > s:
            out.append(v[j])
    return out


def main():
    idx = load_reps()
    iv = []
    for line in open(SD):
        f = line.split()
        if len(f) >= 3:
            iv.append((f[0], int(f[1]), int(f[2])))

    cls = collections.Counter()
    exonic_bp = []
    for c, s, e in iv:
        ov = overlappers(idx, c, s, e)
        if not ov:
            cls["(a) NEVER A NODE"] += 1
            continue
        spliced = [r for r in ov if r[2] >= 2]
        if not spliced:
            cls["(s) STUB REP ONLY — INDETERMINATE (exon array == span)"] += 1
            continue
        bp = 0
        for (_, _, _, _, ex) in spliced:
            for a, b in ex:
                bp += max(0, min(e, b) - max(s, a))
        if bp > 0:
            cls["(b-EXONIC) overlaps a spliced rep's EXON — real transcript is there"] += 1
            exonic_bp.append(bp)
        else:
            cls["(b-INTRONIC) inside a spliced rep's INTRON only — CORRECT rejection"] += 1

    n = sum(cls.values())
    print(f"SD-predicted catalog-absent intervals: {n}")
    print("\n=== classification on EXONS (supersedes the span version) ===")
    for k, v in cls.most_common():
        print(f"  {k:<62} {v:>5}/{n} = {v/n:.4f}")

    det = n - cls["(s) STUB REP ONLY — INDETERMINATE (exon array == span)"]
    ex_ = cls["(b-EXONIC) overlaps a spliced rep's EXON — real transcript is there"]
    inn = cls["(b-INTRONIC) inside a spliced rep's INTRON only — CORRECT rejection"]
    print(f"\n  DETERMINATE subset (stubs excluded): {det}/{n} = {det/n:.4f}")
    if ex_ + inn:
        print(f"  Of intervals with a SPLICED rep ({ex_+inn}): "
              f"EXONIC {ex_/(ex_+inn):.4f}   INTRONIC {inn/(ex_+inn):.4f}")
        print(f"  ⟹ the span-based classification counted all {ex_+inn} as 'a node is here';"
              f" only {ex_} have transcript sequence in the interval.")
    if exonic_bp:
        import statistics as st
        print(f"  median exonic bp in-interval where exonic: {st.median(exonic_bp):.0f}")

    merged, cur = 0, None
    for c, s, e in sorted(iv):
        if cur and cur[0] == c and s < cur[2]:
            cur = (c, cur[1], max(cur[2], e))
        else:
            merged += 1; cur = (c, s, e)
    print(f"\n  ⚠ effective n after merging overlaps: {merged} blocks (not {n} intervals)")


if __name__ == "__main__":
    main()
