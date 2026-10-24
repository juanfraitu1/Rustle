#!/usr/bin/env python3
"""Arm (a) of docs/archive/2026-10/merged/PREREG_popfix_2026-10-12.md: descriptive table sanity on the YAG families, no bar.

    sanity.py POP_TSV YAG_FAMILIES_JSON [RUN_PREFIX]
POP_TSV = psv_population.tsv; YAG_FAMILIES_JSON = prep_yag.py's yag_families.json (its "families" keys are the families).
Writes JSON to stdout:
  status_x_method  per YAG family: n_sites and {status: {method: n}} over the table rows of that family
  column_fractions (only with RUN_PREFIX, a copy_assign --dump-star run with a `positions` column): per YAG family, over
                   the read-star columns of CONTESTED molecules (assignments.tsv contested == 1), the columns at which
                   every candidate that has a genome position is a catalog PSV site of the table (ANY status) and the
                   columns at which every such candidate is FIXED. A column no candidate has a position at is in neither
                   numerator (counted in n_cols_no_position); frac_* divide by all columns, frac_*_of_positioned by the
                   columns with >= 1 position. A candidate's chrom is its copy's in <RUN_PREFIX>.family_join.tsv.

COLUMN RULE, exactly: a column counts when EVERY candidate WITH a position there has its (chrom, pos) in the set;
candidates without a position are ignored (but a column where none has one never counts). This is a per-COLUMN summary,
not the certificate's mask. The prereg (Amendment 1) and the Rust fixed-only certificate mask per PAIR: a column counts
for candidates (a, b) iff BOTH have a position there and both (chrom, pos) are in the set. The two agree only for a
two-candidate column at which both candidates have a position. With >= 3 candidates a column can be in the set for one
pair and not another (here it is all or nothing); a column where a candidate has no position can count here, but never
counts for any pair that includes that candidate.
"""
import collections, json, os, sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "yag_o2"))
from split_half import parse_positions, read_chroms  # noqa: E402


def iter_tsv(path):
    with open(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            yield dict(zip(head, line.rstrip("\n").split("\t")))


def status_method_counts(pop_tsv, families):
    fams = set(families)
    res = {f: {"n_sites": 0, "status_x_method": {}} for f in sorted(fams)}
    for r in iter_tsv(pop_tsv):
        f = r["family_id"]
        if f in fams:
            res[f]["n_sites"] += 1
            by = res[f]["status_x_method"].setdefault(r["status"], {})
            by[r["method"]] = by.get(r["method"], 0) + 1
    return res


def site_sets(pop_tsv):
    """(any_sites, fixed_sites): {(chrom, pos)} of every row / of the FIXED rows."""
    any_s, fixed = set(), set()
    for r in iter_tsv(pop_tsv):
        k = (r["chrom"], int(r["pos"]))
        any_s.add(k)
        if r["status"] == "FIXED":
            fixed.add(k)
    return any_s, fixed


def column_site_fractions(prefix, pop_tsv, families):
    fams = set(families)
    any_s, fixed = site_sets(pop_tsv)
    chroms = read_chroms(prefix)
    contested = {(r["read_name"], r["family_id"]) for r in iter_tsv(prefix + ".assignments.tsv") if r.get("contested") == "1"}
    c = {f: collections.Counter() for f in sorted(fams)}
    for r in iter_tsv(prefix + ".star_reads.tsv"):
        f = r["family_id"]
        if f not in fams or (r["read_name"], f) not in contested:
            continue
        cand = [int(x) for x in r["candidates"].split(",")] if r["candidates"] else []
        n_cols = int(r["n_cols"])
        gpos = parse_positions(r.get("positions", ""), n_cols, len(cand))
        ch = [chroms.get((f, ci)) for ci in cand]
        c[f]["n_molecules"] += 1
        for col in gpos:
            keys = [(h, g) for h, g in zip(ch, col) if g is not None]
            c[f]["n_cols"] += 1
            if not keys:
                c[f]["n_cols_no_position"] += 1
                continue
            c[f]["n_cols_all_any"] += all(k in any_s for k in keys)
            c[f]["n_cols_all_fixed"] += all(k in fixed for k in keys)
    out = {}
    for f, v in c.items():
        n, npos = v["n_cols"], v["n_cols"] - v["n_cols_no_position"]
        out[f] = {k: v[k] for k in ("n_molecules", "n_cols", "n_cols_no_position", "n_cols_all_any", "n_cols_all_fixed")}
        out[f].update({
            "frac_any": round(v["n_cols_all_any"] / n, 4) if n else None,
            "frac_fixed": round(v["n_cols_all_fixed"] / n, 4) if n else None,
            "frac_any_of_positioned": round(v["n_cols_all_any"] / npos, 4) if npos else None,
            "frac_fixed_of_positioned": round(v["n_cols_all_fixed"] / npos, 4) if npos else None})
    return out


def main(pop_tsv, yag_json, prefix=None):
    families = list(json.load(open(yag_json))["families"])
    res = {"status_x_method": status_method_counts(pop_tsv, families)}
    if prefix:
        res["column_fractions"] = column_site_fractions(prefix, pop_tsv, families)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    if len(sys.argv) not in (3, 4):
        raise SystemExit(__doc__)
    main(*sys.argv[1:])
