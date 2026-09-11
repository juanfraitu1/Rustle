#!/usr/bin/env python3
"""Cross-individual O3 differential (docs/superpowers/specs/2026-09-10-o3-flag-pass-integration-design.md,
"Follow-on"): diffs two --flag-missing-copies family_join.tsv outputs -- same catalog, different --bam --
by flag CATEGORY, not raw read presence (the read-presence approach, ledger section 4l, already showed it
cannot distinguish real copy-number difference from tissue-driven expression difference).

EXPLORATORY, not validated: there is no known-true-positive gorilla case to check any "candidate
differential" against (see the design doc's "What's NOT established" section). Sanity-check the
`no_signal fraction` line before reading anything into the `candidate_differential` count -- if it is not
the overwhelming majority, the shared catalog itself (built from one arm's own reads) is the likely
confound, not individual biology.

usage: o3_cross_individual_diff.py <arm_a>.family_join.tsv <arm_b>.family_join.tsv --label-a testis --label-b fibroblast
"""
import argparse
import collections
import csv


def load(path):
    """Load a `<out>.family_join.tsv` into {(catalog_family_id, catalog_copy_idx): o3_flag}.

    A single catalog copy can appear on more than one row of the same file: `copy_assign --families`
    can independently re-process the same underlying catalog copy under two different LOCALLY-detected
    `family_id` groupings (column 1) when nearby regions/windows both reach it -- the `catalog_family_id`/
    `catalog_copy_idx` columns (the stable catalog identity this script joins on) still agree, but the row
    is duplicated. Silently keeping "whichever row happened to be read last" would be a real bug if the two
    rows ever disagreed on `o3_flag` -- so this is checked and reported, not assumed away.
    """
    rows = {}
    dupe_keys = set()
    conflicts = []
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if "o3_flag" not in r:
                raise SystemExit(f"{path}: no o3_flag column -- was --flag-missing-copies set for this run?")
            k = (r["catalog_family_id"], r["catalog_copy_idx"])
            if k in rows:
                dupe_keys.add(k)
                if rows[k] != r["o3_flag"]:
                    conflicts.append((k, rows[k], r["o3_flag"]))
            rows[k] = r["o3_flag"]
    if conflicts:
        detail = "; ".join(f"{fam}:{cidx} was {v1!r} then {v2!r}" for (fam, cidx), v1, v2 in conflicts)
        raise SystemExit(
            f"{path}: {len(conflicts)} catalog (family_id, copy_idx) key(s) appear on multiple rows with "
            f"DIFFERING o3_flag values -- cannot pick a winner silently: {detail}"
        )
    if dupe_keys:
        print(
            f"NOTE: {path}: {len(dupe_keys)} catalog copy key(s) appear on more than one row "
            f"(same underlying copy processed under >1 local family grouping); every duplicate agreed on "
            f"o3_flag, so the row count for this file ({sum(1 for _ in open(path)) - 1}) is "
            f"{len(dupe_keys)} higher than the distinct-copy count used below ({len(rows)})."
        )
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("arm_a")
    ap.add_argument("arm_b")
    ap.add_argument("--label-a", default="A")
    ap.add_argument("--label-b", default="B")
    a = ap.parse_args()
    if a.label_a == "A" or a.label_b == "B":
        print(
            "WARNING: --label-a/--label-b were not both given explicitly -- output below uses the bare "
            "defaults 'A'/'B', which is ambiguous once this is read back later. Re-run with explicit labels.",
        )
    ra, rb = load(a.arm_a), load(a.arm_b)
    keys = sorted(set(ra) | set(rb))
    cats = collections.Counter()
    candidates = []
    for k in keys:
        fa, fb = ra.get(k, "NA"), rb.get(k, "NA")
        if fa == "missing_copy" and fb == "none":
            cats["candidate_differential"] += 1
            candidates.append((k, a.label_a, a.label_b))
        elif fb == "missing_copy" and fa == "none":
            cats["candidate_differential"] += 1
            candidates.append((k, a.label_b, a.label_a))
        elif fa == "missing_copy" and fb == "missing_copy":
            cats["shared"] += 1
        elif "untestable" in (fa, fb) and "missing_copy" in (fa, fb):
            cats["inconclusive"] += 1
        elif fa == "untestable" and fb == "untestable":
            cats["no_information"] += 1
        elif fa == "none" and fb == "none":
            cats["no_signal"] += 1
        else:
            cats["other"] += 1
    print(f"{len(keys)} copies compared")
    for c, n in sorted(cats.items(), key=lambda kv: -kv[1]):
        print(f"  {c}: {n}")
    na_only = sum(1 for k in keys if ra.get(k, "NA") == "NA" or rb.get(k, "NA") == "NA")
    n_both_testable = len(keys) - na_only
    print(
        f"\nof these, {n_both_testable} were actually scored by O3 in BOTH arms ({na_only} are present in "
        f"only one arm's catalog -- excluded/never run in the other -- and so were never testable there at "
        f"all; do not call the full {len(keys)}-copy union \"testable\")."
    )
    no_signal_frac = cats["no_signal"] / max(1, len(keys))
    print(
        f"\nno_signal fraction: {no_signal_frac:.3f} -- per the design doc, this should be the "
        f"overwhelming majority; if it is not, the shared catalog itself (built from one arm's own "
        f"reads) is the likely confound, not individual biology."
    )
    print(f"\ncandidate differentials ({len(candidates)}):")
    for (fam, cidx), flagged_in, clean_in in candidates:
        print(f"  {fam}:{cidx}  missing_copy in {flagged_in}, none in {clean_in}")


if __name__ == "__main__":
    main()
