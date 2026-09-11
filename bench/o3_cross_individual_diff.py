#!/usr/bin/env python3
"""Cross-individual O3 differential (docs/superpowers/specs/2026-09-10-o3-flag-pass-integration-design.md,
"Follow-on"): diffs two --flag-missing-copies family_join.tsv outputs -- same catalog, different --bam,
run on a DIFFERENT ANIMAL (not just a different tissue of the same animal; e.g. the reference run's own
gorilla, OR6737, testis vs a different individual's matched-cell-line fibroblast) -- by flag CATEGORY, not
raw read presence (the read-presence approach, ledger section 4l, already showed it cannot distinguish
real copy-number difference from tissue-driven expression difference, and cannot distinguish either of
those from an individual-level difference).

EXPLORATORY, not validated: there is no known-true-positive gorilla case to check any "candidate
differential" against (see the design doc's "What's NOT established" section). Sanity-check the
`no_signal fraction` line before reading anything into the `candidate_differential` count -- if it is not
the overwhelming majority, the shared catalog itself (built from one arm's own reads) is the likely
confound, not individual biology.

Fix 2 (final whole-branch review): `copy_assign --flag-missing-copies` now writes a distinct `not_tested`
token (instead of overloading `none`) when a copy's `o3_flags` lookup never resolved at all -- fewer than
3 rejected reads reached `detect_missing_copy_pairs`, or (pre-Fix-1) a namespace mismatch. `not_tested` is
NOT a real negative: it must never let a pair count as `candidate_differential`, `shared`, or contribute
to `no_signal` the way a genuinely-tested-and-clean `none` does. It is treated here as equivalent to
`untestable` (both mean "we have no information from this side") for `inconclusive`/`no_information`
purposes; paired with a clean `none` on the other side it falls to `other` (matching the existing
`none`/`untestable` mixed-pair convention), never `no_signal`.

usage: o3_cross_individual_diff.py <arm_a>.family_join.tsv <arm_b>.family_join.tsv --label-a testis --label-b fibroblast [--alpha 0.001]
"""
import argparse
import collections
import csv

# Both `untestable` (a real FlaggedPair result -- p_uncorrected was None, usually an empty/insufficient
# control pool) and `not_tested` (the o3_flags lookup never resolved at all) mean "no information from
# this side" for the purposes of this diff -- neither is a real negative the way `none` is.
NO_INFO = {"untestable", "not_tested"}


def load(path):
    """Load a `<out>.family_join.tsv` into ({(catalog_family_id, catalog_copy_idx): o3_flag}, n_pairs).

    A single catalog copy can appear on more than one row of the same file: `copy_assign --families`
    can independently re-process the same underlying catalog copy under two different LOCALLY-detected
    `family_id` groupings (column 1) when nearby regions/windows both reach it -- the `catalog_family_id`/
    `catalog_copy_idx` columns (the stable catalog identity this script joins on) still agree, but the row
    is duplicated. Silently keeping "whichever row happened to be read last" would be a real bug if the two
    rows ever disagreed on `o3_flag` -- so this is checked and reported, not assumed away.

    `n_pairs` is the count of distinct copies in this file with a REAL computed p-value (flag `missing_copy`
    or `none` -- the only two flags `finalize_flags` ever assigns a p-value to; `untestable` and
    `not_tested` both mean no p-value was computed). This is exactly the denominator
    `finalize_flags`/`copy_assign` used for this arm's own Bonferroni threshold (`alpha / n_pairs`).
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
    n_pairs = sum(1 for v in rows.values() if v in ("missing_copy", "none"))
    return rows, n_pairs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("arm_a")
    ap.add_argument("arm_b")
    ap.add_argument("--label-a", default="A")
    ap.add_argument("--label-b", default="B")
    ap.add_argument(
        "--alpha", type=float, default=0.001,
        help="O3 Bonferroni alpha used for BOTH runs (copy_assign --o3-alpha, default 0.001). This script "
             "cannot read the alpha each run was actually invoked with from family_join.tsv -- pass it "
             "explicitly if either run used a non-default --o3-alpha.",
    )
    a = ap.parse_args()
    if a.label_a == "A" or a.label_b == "B":
        print(
            "WARNING: --label-a/--label-b were not both given explicitly -- output below uses the bare "
            "defaults 'A'/'B', which is ambiguous once this is read back later. Re-run with explicit labels.",
        )
    ra, n_pairs_a = load(a.arm_a)
    rb, n_pairs_b = load(a.arm_b)
    # Per-arm Bonferroni thresholds can differ (different n_pairs -- e.g. one arm's catalog subset excludes
    # more zero-read copies than the other's), so both are reported, not assumed equal.
    thresh = {
        a.label_a: a.alpha / max(1, n_pairs_a),
        a.label_b: a.alpha / max(1, n_pairs_b),
    }
    keys = sorted(set(ra) | set(rb))
    cats = collections.Counter()
    candidates = []
    inconclusives = []
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
        elif (fa in NO_INFO or fb in NO_INFO) and "missing_copy" in (fa, fb):
            cats["inconclusive"] += 1
            flagged_in = a.label_a if fa == "missing_copy" else a.label_b
            no_info_in = a.label_b if fa == "missing_copy" else a.label_a
            no_info_val = fb if fa == "missing_copy" else fa
            inconclusives.append((k, flagged_in, no_info_in, no_info_val))
        elif fa in NO_INFO and fb in NO_INFO:
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
    print(
        f"\nBonferroni thresholds (alpha={a.alpha:g} / n_pairs, n_pairs = copies in that arm with a real "
        f"p-value, i.e. flag missing_copy or none):"
    )
    print(f"  {a.label_a}: {a.alpha:g} / {n_pairs_a} = {thresh[a.label_a]:.3e}")
    print(f"  {a.label_b}: {a.alpha:g} / {n_pairs_b} = {thresh[a.label_b]:.3e}")
    print(f"\ncandidate differentials ({len(candidates)}):")
    for (fam, cidx), flagged_in, clean_in in candidates:
        print(
            f"  {fam}:{cidx}  missing_copy in {flagged_in} (p_threshold={thresh[flagged_in]:.3e}), "
            f"none in {clean_in} (p_threshold={thresh[clean_in]:.3e})"
        )
    print(f"\ninconclusive pairs ({len(inconclusives)}):")
    for (fam, cidx), flagged_in, no_info_in, no_info_val in inconclusives:
        print(
            f"  {fam}:{cidx}  missing_copy in {flagged_in} (p_threshold={thresh[flagged_in]:.3e}), "
            f"{no_info_val} in {no_info_in} (p_threshold={thresh[no_info_in]:.3e})"
        )


if __name__ == "__main__":
    main()
