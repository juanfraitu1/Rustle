#!/usr/bin/env python3
"""Post hoc diagnostics of the O1 hierarchy evaluation (`docs/PREREG_o1_hierarchy_eval_2026-10-13.md`, Results). They
were added after the run and change no verdict. The pinned scorer `hier_match.py` (8f0515ef) is imported, never edited.

    posthoc_checks.py OUT.json --dev HIER_JSON --heldout HIER_JSON
                      [--lift-ends CHAIN SEDEF_BED V1_FA V2_FA] [--window 1000]

Each HIER_JSON is a scored run (`hier_match.py ... OUT.json`). Its `inputs` (members, truth, fold, loci GFF, map) are
reloaded and the hierarchy is re-matched without nulls. The re-match must equal the JSON's `per_truth` placements
(hierarchy and flat) and its hierarchy mean Jaccard, so every number here belongs to the scored run.

1. **distinct-sets.** The one-to-one matching runs over NODES. Identical member sets at several levels (a family and its
   only expansion, ...) are several columns, so two truth families can match one member set through two nodes. The
   re-match uses one column per DISTINCT member set, labelled by its lowest-rank (level, id) (the scorer's column
   order), with the same Hungarian and canonical placement. Reported next to the as-scored and flat matches:
   - the hierarchy - flat gains, and the share of the mean-Jaccard gain that the distinct-set re-match removes;
   - `shared_sets_as_scored`: member sets that carry >= 2 one-to-one matches in the scored run. `extra_matches` counts
     the matches beyond the first on each set, split by whether the truth families on the set have different member
     sets or are twins (duplicate truth families with one set);
   - `duplicate_truth_sets`: evaluable truth families with identical mapped member sets.
2. **size-strata (E2).** The mapped set size is the evaluable truth family's member count in that catalog. Reported:
   - per catalog, the exact rate per size;
   - the dev rate reweighted to the held-out size mix, on strata 2, 3, 4 and >= 5 (and by exact size);
   - the truth families evaluable in both catalogs: their median sizes and 2-member strata.
3. **lift-ends** (`--lift-ends`). SEDEF pair sides are lifted v1.0 -> v2.0 by `v1_to_v2.py`. The sides with equal v1.0
   and v2.0 lengths whose v1.0 interval starts or ends within `--window` bp of a v1.0 chromosome end are compared by
   sequence. The lift's own check sampled 1,000 random same-length intervals out of all sides.
"""
import argparse
import collections
import json
import os
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import hier_match as H  # noqa: E402  (pinned at 8f0515ef: imported read-only)

KEYS = ("exact", "exact_rate", "mean_jaccard", "n_matched", "by_level", "exact_by_level")


def distinct_nodes(nodes):
    """{(level, id): member set} -> one column per DISTINCT member set, labelled by its lowest-rank (level, id) node."""
    out, seen = {}, set()
    for k in sorted(nodes, key=lambda k: (H.RANK[k[0]], k[1])):
        if nodes[k] not in seen:
            seen.add(nodes[k])
            out[k] = nodes[k]
    return out


def duplicate_truth_sets(truth):
    """Groups (labels sorted) of truth families with one identical member set, sorted."""
    by = collections.defaultdict(list)
    for label in sorted(truth):
        by[truth[label]].append(label)
    return sorted(v for v in by.values() if len(v) > 1)


def shared_sets(matches, nodes, truth):
    """Member sets that carry >= 2 one-to-one matches. `extra_matches` = the matches beyond the first on each set,
    split into sets whose truth families have different member sets and sets held by twins of one truth set."""
    by = collections.defaultdict(list)
    for m in matches:
        if m.level:
            by[nodes[(m.level, m.node_id)]].append(m.label)
    groups = [v for v in by.values() if len(v) > 1]
    extra = sum(len(v) - 1 for v in groups)
    twins = sum(len(v) - 1 for v in groups if len({truth[x] for x in v}) == 1)
    return {"shared_sets": len(groups), "matches_on_shared_sets": sum(len(v) for v in groups), "extra_matches": extra,
            "extra_matches_different_truth": extra - twins, "extra_matches_twins": twins}


def stratum(size, top=5):
    """The size stratum: '2', '3', '4', then f'{top}+' for every larger size; the exact size with top=None."""
    return str(size) if top is None or size < top else f"{top}+"


def _order(k):
    return int(k.rstrip("+"))


def rates_by(sizes, exact, key):
    """{label: size}, {label: exact?}, size -> stratum  ->  {stratum: (exact, families)}, strata in size order."""
    c = collections.defaultdict(lambda: [0, 0])
    for label, s in sizes.items():
        c[key(s)][0] += bool(exact[label])
        c[key(s)][1] += 1
    return {k: tuple(c[k]) for k in sorted(c, key=_order)}


def reweight(dev, held_out):
    """The dev exact rate under the held-out stratum mix: sum_k n_held(k) / N_held * exact_dev(k) / n_dev(k).
    dev = {stratum: (exact, families)}, held_out = {stratum: families}. None when a held-out stratum has no dev
    family."""
    if any(dev.get(k, (0, 0))[1] == 0 for k in held_out):
        return None
    n = sum(held_out.values())
    return sum(held_out[k] / n * dev[k][0] / dev[k][1] for k in held_out)


def near_end_sides(regions, lengths, window):
    """Lifted sides ((v1 chrom, s, e), (v2 chrom, s, e)), 0-based half-open, with equal lengths (no indel inside) whose
    v1.0 interval starts within `window` bp of its chromosome's start or ends within `window` bp of its end."""
    return [p for p in regions if p[0][2] - p[0][1] == p[1][2] - p[1][1]
            and (p[0][1] < window or lengths[p[0][0]] - p[0][2] < window)]


def differing_sides(v1_fa, v2_fa, sides):
    """The sides whose v1.0 and v2.0 sequences differ (`samtools faidx`, upper case)."""
    import v1_to_v2 as V
    a, b = V.fetch(v1_fa, [p[0] for p in sides]), V.fetch(v2_fa, [p[1] for p in sides])

    def key(r):
        return f"{r[0]}:{r[1] + 1}-{r[2]}"
    return [p for p in sides if a[key(p[0])] != b[key(p[1])]]


def lift_ends(chain, sedef, v1_fa, v2_fa, window):
    import v1_to_v2 as V
    with open(v1_fa + ".fai") as fh:
        lengths = {f[0]: int(f[1]) for f in (line.split("\t") for line in fh)}
    with open(chain) as fh:
        lifter = V.Lifter(V.read_chain(fh))
    with open(sedef) as fh:
        _, counts, regions = V.lift_rows(lifter, fh, pairs=True)
    near = near_end_sides(regions, lengths, window)
    diff = differing_sides(v1_fa, v2_fa, near)

    def fmt(p):
        return f"{p[0][0]}:{p[0][1]}-{p[0][2]} -> {p[1][0]}:{p[1][1]}-{p[1][2]}"
    return {"window": window, "rows_lifted": counts["rows_lifted"], "sides_lifted": len(regions),
            "same_length_sides": sum(1 for p in regions if p[0][2] - p[0][1] == p[1][2] - p[1][1]),
            "near_end_sides": len(near), "near_end_sides_differing": len(diff),
            "near_end_distinct_intervals": len(set(near)), "near_end_distinct_intervals_differing": len(set(diff)),
            "differing_sides_v1_to_v2_0based": [fmt(p) for p in diff]}


def catalog(path):
    """A scored run's JSON -> (scored JSON, nodes, members, evaluable truth sets), loaded from the run's own inputs."""
    with open(path) as fh:
        scored = json.load(fh)
    a = scored["inputs"]
    nodes, family_of = H.load_nodes(a["members"])
    members = sorted(family_of)
    sets, _, _ = H.truth_sets(H.load_truth(a["truth"]), members, H.load_fold(a["fold"]) if a.get("fold") else None,
                              H.load_loci_gff(a["loci_gff"]) if a.get("loci_gff") else None, a["map"])
    return scored, nodes, members, {k: v for k, v in sets.items() if len(v) >= 2}


def distinct_sets_check(scored, nodes, members, ev):
    """-> (report, as-scored matches). Stops if the re-match is not the scored run."""
    fam = {k: v for k, v in nodes.items() if k[0] == "family"}
    dn = distinct_nodes(nodes)
    runs = {name: H.Scorer(nd, members).run(ev) for name, nd in (("as_scored", nodes), ("distinct_sets", dn),
                                                                  ("flat", fam))}
    pt = {r["label"]: r for r in scored["per_truth"]}
    got = {m.label: (m.level, m.node_id, m.jaccard == 1.0) for m in runs["as_scored"]}
    got_flat = {m.label: m.node_id for m in runs["flat"]}
    if got != {k: (r["level"], r["node"], r["exact"]) for k, r in pt.items()} or \
            got_flat != {k: r["flat_node"] for k, r in pt.items()}:
        raise SystemExit("the re-match is not the scored run: per_truth placements differ")
    out = {"n_nodes": len(nodes), "n_distinct_sets": len(dn)}
    for name, ms in runs.items():
        s = H.summary(ms)
        out[name] = {k: s[k] for k in KEYS}
    if out["as_scored"]["mean_jaccard"] != scored["hierarchy"]["mean_jaccard"]:
        raise SystemExit("the re-match is not the scored run: hierarchy mean Jaccard differs")
    for name in ("as_scored", "distinct_sets"):
        out[f"gain_over_flat_{name}"] = {k: out[name][k] - out["flat"][k] for k in ("exact", "mean_jaccard",
                                                                                    "n_matched")}
    g_as = out["gain_over_flat_as_scored"]["mean_jaccard"]
    g_ds = out["gain_over_flat_distinct_sets"]["mean_jaccard"]
    out["share_of_mean_jaccard_gain_removed"] = (g_as - g_ds) / g_as if g_as else None
    out["shared_sets_as_scored"] = shared_sets(runs["as_scored"], nodes, ev)
    out["duplicate_truth_sets"] = duplicate_truth_sets(ev)
    return out, runs["as_scored"]


def _table(r):
    return {k: {"exact": e, "n": n, "rate": e / n} for k, (e, n) in r.items()}


def size_strata(matches, names):
    """matches {'dev'|'heldout': as-scored matches}; names {'dev'|'heldout': catalog name}."""
    size = {role: {m.label: m.t_size for m in ms} for role, ms in matches.items()}
    exact = {role: {m.label: m.jaccard == 1.0 for m in ms} for role, ms in matches.items()}
    by_size = {role: rates_by(size[role], exact[role], lambda s: stratum(s, None)) for role in matches}
    by_stratum = {role: rates_by(size[role], exact[role], stratum) for role in matches}
    rate = {role: sum(exact[role].values()) / len(exact[role]) for role in matches}
    rw = {"strata_2_3_4_5plus": reweight(by_stratum["dev"], {k: n for k, (_, n) in by_stratum["heldout"].items()}),
          "exact_size": reweight(by_size["dev"], {k: n for k, (_, n) in by_size["heldout"].items()})}
    both = sorted(set(size["dev"]) & set(size["heldout"]))
    two = {role: [x for x in both if size[role][x] == 2] for role in matches}
    return {
        "catalogs": names,
        "exact_rate": rate,
        "observed_gap": rate["heldout"] - rate["dev"],
        "per_size": {role: _table(by_size[role]) for role in matches},
        "per_stratum_2_3_4_5plus": {role: _table(by_stratum[role]) for role in matches},
        "dev_rate_reweighted_to_heldout_sizes": rw,
        "gap_after_reweighting": {k: None if v is None else rate["heldout"] - v for k, v in rw.items()},
        "evaluable_in_both": {
            "n": len(both),
            "median_size": {role: statistics.median(size[role][x] for x in both) for role in matches},
            "two_member_stratum": {role: {"exact": sum(exact[role][x] for x in two[role]), "n": len(two[role])}
                                   for role in matches},
        },
    }


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out")
    ap.add_argument("--dev", required=True, help="the dev catalog's scored hier_match.json")
    ap.add_argument("--heldout", required=True, help="the held-out catalog's scored hier_match.json")
    ap.add_argument("--lift-ends", nargs=4, metavar=("CHAIN", "SEDEF_BED", "V1_FA", "V2_FA"))
    ap.add_argument("--window", type=int, default=1000)
    a = ap.parse_args(argv)
    res, matches, names = {"inputs": vars(a), "distinct_sets": {}}, {}, {}
    for role, path in (("dev", a.dev), ("heldout", a.heldout)):
        names[role] = os.path.basename(os.path.dirname(os.path.abspath(path)))
        res["distinct_sets"][names[role]], matches[role] = distinct_sets_check(*catalog(path))
    res["size_strata"] = size_strata(matches, names)
    if a.lift_ends:
        res["lift_ends"] = lift_ends(*a.lift_ends, a.window)
    with open(a.out, "w") as fh:
        json.dump(res, fh, indent=1)
    print(json.dumps(res, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
