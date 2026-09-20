#!/usr/bin/env python3
"""make_psv_filtered_catalog.py — split family_rna_refine.tsv by PSV evidence.

Creates:
  - bench/family_rna_refine_psv_filtered.tsv  (families with PSV evidence)
  - bench/family_rna_refine_psv_rejected.tsv  (families without PSV evidence)
  - bench/family_psv_evidence.tsv             (already produced by filter_families_by_psv.py)

A family is kept if it has >= min_psv PSV columns from the mafft alignment of
member transcripts. Rejected families are likely isoforms, domain sharers, or
spurious minimizer hits.
"""
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REFINE_TSV = os.path.join(HERE, "family_rna_refine.tsv")
EVIDENCE_TSV = os.path.join(HERE, "family_psv_evidence.tsv")
FILTERED_TSV = os.path.join(HERE, "family_rna_refine_psv_filtered.tsv")
REJECTED_TSV = os.path.join(HERE, "family_rna_refine_psv_rejected.tsv")

MIN_PSV = 1


def main():
    evidence = {}
    if os.path.exists(EVIDENCE_TSV):
        for r in csv.DictReader(open(EVIDENCE_TSV), delimiter="\t"):
            evidence[int(r["family_id"])] = r
    else:
        print(f"[!] {EVIDENCE_TSV} not found; run filter_families_by_psv.py first", file=sys.stderr)
        sys.exit(1)

    kept_rows = []
    rejected_rows = []
    with open(REFINE_TSV) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            fid = int(r["family_id"])
            ev = evidence.get(fid, {})
            n_psv = int(ev.get("n_psv", 0))
            r["psv_status"] = ev.get("status", "missing")
            r["n_psv"] = str(n_psv)
            r["n_aligned_cols"] = ev.get("n_aligned_cols", "0")
            if n_psv >= MIN_PSV:
                kept_rows.append(r)
            else:
                rejected_rows.append(r)

    fieldnames = list(csv.DictReader(open(REFINE_TSV), delimiter="\t").fieldnames) + \
                 ["psv_status", "n_psv", "n_aligned_cols"]

    for path, rows in [(FILTERED_TSV, kept_rows), (REJECTED_TSV, rejected_rows)]:
        with open(path, "w", newline="") as o:
            w = csv.DictWriter(o, fieldnames=fieldnames, delimiter="\t")
            w.writeheader()
            for r in rows:
                w.writerow(r)
        print(f"[+] wrote {path}: {len(rows)} rows", file=sys.stderr)

    print(f"\nSummary:")
    print(f"  kept (>= {MIN_PSV} PSV): {len(kept_rows)}")
    print(f"  rejected (0 PSV):        {len(rejected_rows)}")


if __name__ == "__main__":
    main()
