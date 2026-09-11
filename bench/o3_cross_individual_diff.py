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
    rows = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if "o3_flag" not in r:
                raise SystemExit(f"{path}: no o3_flag column -- was --flag-missing-copies set for this run?")
            rows[(r["catalog_family_id"], r["catalog_copy_idx"])] = r["o3_flag"]
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
