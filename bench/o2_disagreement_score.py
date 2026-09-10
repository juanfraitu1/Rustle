#!/usr/bin/env python3
"""Score PREREG_aligner_disagreement (P1-P4): admitted-by-disagreement molecules
vs their PRIMARY copy, and status identity of the pre-existing contested set.

  python3 bench/o2_disagreement_score.py --bam hsa16.bam --copies copies16.tsv \
      --new ours_dis.assignments.tsv --base ours_odi.assignments.tsv

Primary copy per molecule = the -F 2308 record's unit-span overlap in --copies
(same overlap rule the binary uses: read [start,end) vs unit [start,end)).
"""
import argparse, csv, statistics, sys
from collections import Counter

import pysam


def load_units(path):
    units = []
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            units.append((r["chrom"], int(r["start"]), int(r["end"]), r["copy_idx"]))
    return units


def primary_copy_map(bam, units, region=None):
    out = {}
    with pysam.AlignmentFile(bam) as fh:
        it = fh.fetch(region=region) if region else fh
        for a in it:
            if a.is_secondary or a.is_supplementary or a.is_unmapped:
                continue
            s0, e0 = a.reference_start, a.reference_end
            hit = None
            for c, a0, b0, idx in units:
                if c == a.reference_name and s0 < b0 and e0 > a0:
                    hit = idx
                    break
            out[a.query_name] = hit
    return out


def read_rows(path):
    # `assigned_copy` is the binary's internal unit id; `catalog_copy_idx` is the assigned
    # copy's index in --copies (identical in human copies16, NOT in the gorilla catalogs).
    with open(path) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        rows = {r["read_name"]: r for r in rd}
    for r in rows.values():
        r["assigned_copy"] = r.get("catalog_copy_idx") or r["assigned_copy"]
    return rows


def median(xs):
    return statistics.median(xs) if xs else float("nan")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--copies", required=True)
    ap.add_argument("--new", required=True, help="run WITH --admit-aligner-disagreement")
    ap.add_argument("--base", required=True, help="run WITHOUT (same other flags)")
    ap.add_argument("--region", help="restrict the BAM scan (needs .bai), e.g. NC_073242.2:1-2")
    a = ap.parse_args()

    units = load_units(a.copies)
    prim = primary_copy_map(a.bam, units, a.region)
    new, base = read_rows(a.new), read_rows(a.base)

    dis = [r for r in new.values() if r.get("aligner_disagreement") == "1"]
    n_dis = len(dis)
    st = Counter(r["status"] for r in dis)
    assigned = [r for r in dis if r["status"] == "assigned"]
    to_prim = [r for r in assigned if prim.get(r["read_name"]) is not None
               and r["assigned_copy"] == prim[r["read_name"]]]
    n_dec = [int(r["n_decisive"]) for r in assigned]
    marg = [float(r["margin"]) for r in assigned]

    print(f"P1 admitted_by_disagreement = {n_dis}  (pass iff 2000..2600)")
    print(f"   status: " + ", ".join(f"{k}={v}" for k, v in st.most_common()))
    pa = 100.0 * len(assigned) / n_dis if n_dis else float("nan")
    pp = 100.0 * len(to_prim) / len(assigned) if assigned else float("nan")
    print(f"P2 assigned = {len(assigned)}/{n_dis} = {pa:.1f}%  (pass >=60, refuted <40)")
    print(f"   assigned->PRIMARY copy = {len(to_prim)}/{len(assigned)} = {pp:.1f}%  (pass >=80, refuted <60)")
    print(f"   assigned->non-primary breakdown: "
          + ", ".join(f"{k}={v}" for k, v in Counter(
              (prim.get(r['read_name']), r['assigned_copy']) for r in assigned
              if r['assigned_copy'] != prim.get(r['read_name'])).most_common(8)))
    print(f"P3 n_decisive median = {median(n_dec):.1f} (pass >=8, refuted <4); "
          f"margin median = {median(marg):.1f} (pass >=40, refuted <20)")

    # P4: pre-existing contested set in base = contested==1 rows; status must be identical in new
    base_c = {n for n, r in base.items() if r.get("contested") == "1"}
    changed = [(n, base[n]["status"], new.get(n, {}).get("status", "<missing>"))
               for n in base_c if new.get(n, {}).get("status") != base[n]["status"]]
    print(f"P4 base contested = {len(base_c)}; status changed in new = {len(changed)}  (pass iff 0)")
    for n, s0, s1 in changed[:10]:
        print(f"   {n}: {s0} -> {s1}")
    # sanity: the disagreement rows must not be in the base at all (they were gate-skipped)
    overlap = sum(1 for r in dis if r["read_name"] in base)
    print(f"   sanity: disagreement rows already present in base = {overlap} (expect 0)")
    print(f"   sanity: disagreement rows with no primary unit = "
          f"{sum(1 for r in dis if prim.get(r['read_name']) is None)} (expect 0)")


if __name__ == "__main__":
    main()
