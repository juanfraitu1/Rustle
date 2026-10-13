#!/usr/bin/env python3
"""Arm (e) of docs/PREREG_popfix_2026-10-12.md: the implementation check of `copy_assign --psv-population` (bar = exact).

    check_e.py RUN_PREFIX POP_TSV            (env YAG_FAMILIES = yag_families.json restricts to the YAG families)
RUN_PREFIX = a `--psv-population FILE --dump-star` run (star_reads.tsv with `positions`, assignments.tsv,
family_join.tsv, population_tier.tsv). A molecule is ELIGIBLE iff it is contested, has a tier (status_fixed != NA), has
>= 1 read-star column, EVERY column has EVERY candidate's genome position present and FIXED in POP_TSV (chrom of a
candidate from family_join.tsv), and its assignments row has origin_rejected == 0 and tie_outside_catalog == 0 (no
later gate changed the call). Where the tier's status_all or status_fixed is `assigned`, BOTH must be `assigned`
with the same copy (copy_all == copy_fixed); anything else is a MISMATCH (a bug). Exact ties among non-assigned
molecules may name a different best copy (ruling R16) and are not compared.
Prints JSON {n_eligible, n_assigned_any, n_mismatch, n_eligible_no_columns, first_mismatches[<= 10]}; exit 1 on any
mismatch. A molecule with no column is counted in n_eligible_no_columns, not n_eligible (it cannot be assigned).
"""
import json, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "yag_o2"))
from sanity import iter_tsv  # noqa: E402
from split_half import load_sites, parse_positions, read_chroms  # noqa: E402


def mismatch_kind(tier):
    s_all, s_fix = tier["status_all"], tier["status_fixed"]
    if s_all == "assigned" and s_fix == "assigned":
        return None if tier["copy_all"] == tier["copy_fixed"] else "both_assigned_different_copy"
    if s_all == "assigned":
        return "only_status_all_assigned"
    if s_fix == "assigned":
        return "only_status_fixed_assigned"
    return None


def main(prefix, pop_tsv):
    sites = load_sites(pop_tsv)
    chroms = read_chroms(prefix)
    only = set(json.load(open(os.environ["YAG_FAMILIES"]))["families"]) if os.environ.get("YAG_FAMILIES") else None
    asg = {(r["read_name"], r["family_id"]): r for r in iter_tsv(prefix + ".assignments.tsv")}
    tiers = {(r["read_name"], r["family_id"]): r for r in iter_tsv(prefix + ".population_tier.tsv")
             if r["status_fixed"] != "NA"}
    n_elig = n_any = n_no_cols = 0
    mism = []
    for r in iter_tsv(prefix + ".star_reads.tsv"):
        fam, key = r["family_id"], (r["read_name"], r["family_id"])
        if (only is not None and fam not in only) or key not in tiers:
            continue
        a = asg.get(key)
        if a is None or a.get("contested") != "1" or a.get("origin_rejected") != "0" or a.get("tie_outside_catalog") != "0":
            continue
        cand = [int(x) for x in r["candidates"].split(",")] if r["candidates"] else []
        gpos = parse_positions(r.get("positions", ""), int(r["n_cols"]), len(cand))
        ch = [chroms.get((fam, ci)) for ci in cand]
        if not all(g is not None and h is not None and (h, g) in sites for col in gpos for h, g in zip(ch, col)):
            continue
        if not gpos:
            n_no_cols += 1
            continue
        n_elig += 1
        tier = tiers[key]
        if tier["status_all"] == "assigned" or tier["status_fixed"] == "assigned":
            n_any += 1
            kind = mismatch_kind(tier)
            if kind:
                mism.append({"read_name": r["read_name"], "family_id": fam, "kind": kind,
                             "status_all": tier["status_all"], "copy_all": tier["copy_all"],
                             "status_fixed": tier["status_fixed"], "copy_fixed": tier["copy_fixed"],
                             "n_dec_fixed": tier["n_dec_fixed"], "n_cols": len(gpos)})
    print(json.dumps({"n_eligible": n_elig, "n_assigned_any": n_any, "n_mismatch": len(mism),
                      "n_eligible_no_columns": n_no_cols, "first_mismatches": mism[:10]}, indent=1))
    return 1 if mism else 0


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit(__doc__)
    sys.exit(main(*sys.argv[1:]))
