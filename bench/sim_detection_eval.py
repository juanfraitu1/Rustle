#!/usr/bin/env python3
"""Family-DETECTION sensitivity/precision on the fully-simulated ground-truth genome (O1).

NON-CIRCULAR: `sim_genome.py` PLANTS the families (positions, copy number, divergence) and labels the truth
in `simgw_truth.tsv`; we run `gw_family_catalog` (O1) and score the detected `simdet.copies.tsv` against the
plant. Isolates the ALGORITHM's ceiling from real-data confounds (coverage/expression) — every planted copy
is well-covered and expressed, so a miss would be an algorithm gap, not a data gap.

Detection is SPATIAL (a family = >=2 distinct loci), so even the K=0 identical-copy family is DETECTED at
the member level; K=0 is a per-READ ASSIGNMENT limit (O2), not a detection limit.

Outputs: bench/sim_detection.tsv (per-family) + a summary block to stdout.
Run: /home/juanfra/miniforge3/bin/python bench/sim_detection_eval.py
"""
import csv
from collections import defaultdict

OUT = "/home/juanfra/winloci_scratch"
TRUTH = f"{OUT}/simgw_truth.tsv"
DETECTED = f"{OUT}/simdet.copies.tsv"
DEST = "bench/sim_detection.tsv"

# families sim_genome plants as SINGLE-copy / over-merge CONTROLS (must NOT form a >=2-copy family).
CONTROL_PREFIXES = ("single", "domshare")


def load_intervals(path, fam_col, chrom_col, start_col, end_col, skip_header=True):
    fams = defaultdict(list)  # family -> [(chrom, start, end)]
    with open(path) as f:
        r = csv.reader(f, delimiter="\t")
        if skip_header:
            next(r, None)
        for row in r:
            if not row or len(row) <= max(fam_col, chrom_col, start_col, end_col):
                continue
            fams[row[fam_col]].append((row[chrom_col], int(row[start_col]), int(row[end_col])))
    return fams


def overlaps(a, dets):
    ch, s, e = a
    return any(dc == ch and not (ds > e or de < s) for (dc, ds, de) in dets)


def main():
    planted = load_intervals(TRUTH, 0, 2, 3, 4)          # family, chrom, start, end
    detected = load_intervals(DETECTED, 0, 3, 4, 5)      # family_id, chrom, start, end
    all_det = [iv for ivs in detected.values() for iv in ivs]

    planted_mc = {f: c for f, c in planted.items() if len(c) >= 2}     # planted MULTI-copy (the recall target)
    planted_ctrl = {f: c for f, c in planted.items() if not (len(c) >= 2)}  # single-copy controls

    rows = []
    fam_detected = 0
    total_members = 0
    total_members_found = 0
    for fam, copies in sorted(planted_mc.items()):
        found = sum(1 for c in copies if overlaps(c, all_det))
        total_members += len(copies)
        total_members_found += found
        # which detected family(ies) do this family's copies land in?
        hit_dets = sorted({df for df, divs in detected.items() for c in copies if overlaps(c, divs)})
        is_detected = found >= 2  # a family needs >=2 detected member loci
        if is_detected:
            fam_detected += 1
        rows.append([fam, len(copies), found, f"{found}/{len(copies)}",
                     "DETECTED" if is_detected else "MISSED", ";".join(hit_dets)])

    # PRECISION: detected families whose copies do NOT overlap any planted multi-copy family = false;
    # controls that formed a >=2-copy detected family = over-merge FP.
    planted_mc_ivs = [iv for f in planted_mc for iv in planted_mc[f]]
    false_families = []
    for df, divs in detected.items():
        if not any(overlaps(iv, planted_mc_ivs) for iv in divs):
            false_families.append(df)
    ctrl_merged = [f for f in planted_ctrl if any(overlaps(c, all_det) and
                    any(len(detected[df]) >= 2 and overlaps(c, detected[df]) for df in detected)
                    for c in planted_ctrl[f])]

    with open(DEST, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["planted_family", "planted_copies", "members_detected", "sensitivity", "status", "detected_as"])
        w.writerows(rows)

    n_mc = len(planted_mc)
    print(f"\n=== FAMILY-DETECTION on planted ground truth (O1) ===")
    print(f"  multi-copy families:  {fam_detected}/{n_mc} DETECTED  ({100*fam_detected/n_mc:.0f}% family sensitivity)")
    print(f"  family MEMBERS (copies): {total_members_found}/{total_members} detected  ({100*total_members_found/total_members:.0f}% member sensitivity)")
    print(f"  false families (precision): {len(false_families)} spurious  -> precision {100*(len(detected)-len(false_families))/max(1,len(detected)):.0f}%")
    print(f"  single-copy / domshare controls that WRONGLY merged: {len(ctrl_merged)} (expect 0)")
    print(f"  wrote per-family table -> {DEST}")
    for r in rows:
        print("   ", "\t".join(str(x) for x in r))


if __name__ == "__main__":
    main()
