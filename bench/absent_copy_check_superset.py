#!/usr/bin/env python3
"""Superset checker for --absent-copies two-stage freeze guarantee.

Given an OFF assignments.tsv (no --absent-copies) and an ON assignments.tsv
(--absent-copies enabled), verifies that every read ASSIGNED in OFF is still
ASSIGNED to the SAME copy in ON (the freeze guarantee: Stage-1 results must
not be overwritten for reads that were already resolved).

CONDITION ON RUSTLE_XFAM_RECONCILE (added 2026-08-31). The freeze guarantee this script enforces holds
only when BOTH arms were run with the SAME value of RUSTLE_XFAM_RECONCILE. Under
RUSTLE_XFAM_RECONCILE=abstain a molecule assigned in two families at once whose two claims rest on
DIFFERENT alignment records naming DISJOINT copy intervals is demoted to `ambiguous` in EVERY family it
is assigned in — deliberately, because a molecule has one origin. Those are exactly the flips this script
reports as violations: 221 of them on the 12-family mec batch (221/53,715 assigned rows, 110 molecules).
They are NOT --absent-copies overwriting a Stage-1 result, and comparing an OFF arm run without the flag
against an ON arm run with it will report them as such. Set the env var identically in both arms, or
subtract the rows listed in <out>.xfam_conflicts.tsv with demoted=true before reading this script's exit
code.

CRITICAL GOTCHA: read_name alone is NOT a unique key — the same read appears in
multiple rows when it multimaps to multiple families.  We key on (read_name,
family_id) and compare complete sorted rows set-wise; a naive join on read_name
would cartesian-explode and produce false mismatches.

Usage:
  python3 bench/absent_copy_check_superset.py OFF.assignments.tsv ON.assignments.tsv

Exit 0 if the superset guarantee holds (0 flips), nonzero on any flip.
Also reports the net assigned delta (ON − OFF) so callers can see if the ON
path admitted extra reads.

Testing with identical files produces 0 flips:
  python3 bench/absent_copy_check_superset.py x.tsv x.tsv
"""

import sys
import csv
from collections import defaultdict


def load(path: str) -> dict[tuple[str, str], list[dict]]:
    """Return {(read_name, family_id): [row, ...]} for rows where status==assigned.
    Multiple rows per key are possible (e.g. multimapper in multiple families).
    """
    by_key: dict[tuple[str, str], list[dict]] = defaultdict(list)
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row.get("status") == "assigned":
                key = (row["read_name"], row["family_id"])
                by_key[key].append(row)
    return dict(by_key)


def main() -> int:
    if len(sys.argv) != 3:
        print("Usage: absent_copy_check_superset.py OFF.tsv ON.tsv", file=sys.stderr)
        return 2

    off_path, on_path = sys.argv[1], sys.argv[2]
    off = load(off_path)
    on  = load(on_path)

    flips: list[tuple[str, str, str, str]] = []  # (read_name, family_id, off_copy, on_copy)
    missing: list[tuple[str, str]] = []           # (read_name, family_id) assigned in OFF but absent in ON

    for (read_name, family_id), off_rows in off.items():
        on_rows = on.get((read_name, family_id))
        if on_rows is None:
            # read+family pair completely absent from ON — assigned read disappeared.
            # Record the KEY ONCE (not once per duplicate row) so the violation count is not inflated.
            missing.append((read_name, family_id))
            continue
        # Compare assigned copies SET-wise: every distinct copy this read+family was assigned to in
        # OFF must still be an ON-assigned copy (dedup so duplicate rows can't inflate flips either).
        off_copies = {r["assigned_copy"] for r in off_rows}
        on_copies  = {r["assigned_copy"] for r in on_rows if r.get("status") == "assigned"}
        for c in sorted(off_copies - on_copies):
            flips.append((read_name, family_id, c,
                          ",".join(sorted(on_copies)) if on_copies else "<absent>"))

    n_off_assigned = sum(len(v) for v in off.values())
    n_on_assigned  = sum(
        1 for rows in on.values()
        for _ in rows  # already filtered to status==assigned
    )
    # Count raw ON assigned (including newly admitted absent-copy reads)
    n_on_total_assigned = 0
    with open(on_path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("status") == "assigned":
                n_on_total_assigned += 1

    print(f"OFF assigned reads : {n_off_assigned}")
    print(f"ON  assigned reads : {n_on_total_assigned}")
    print(f"Net delta (ON-OFF) : {n_on_total_assigned - n_off_assigned:+d}")
    print(f"Missing keys       : {len(missing)}")
    print(f"Copy flips         : {len(flips)}")

    if missing:
        print("\nMISSING (read+family in OFF-assigned but absent in ON):")
        for read_name, family_id in missing[:10]:
            print(f"  {read_name}  {family_id}")
        if len(missing) > 10:
            print(f"  ... and {len(missing) - 10} more")

    if flips:
        print("\nFLIPS (read+family assigned to different copy in ON vs OFF):")
        for read_name, family_id, off_copy, on_copy in flips[:10]:
            print(f"  {read_name}  {family_id}  OFF_copy={off_copy}  ON_copy={on_copy}")
        if len(flips) > 10:
            print(f"  ... and {len(flips) - 10} more")

    n_violations = len(missing) + len(flips)
    if n_violations == 0:
        print("\nSUPERSET CHECK: PASS (0 violations — freeze guarantee holds)")
        return 0
    else:
        print(f"\nSUPERSET CHECK: FAIL ({n_violations} violation(s))")
        return 1


if __name__ == "__main__":
    sys.exit(main())
