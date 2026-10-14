#!/usr/bin/env python3
"""PREREG_mcl_two_node_2026-10-13 F2: two clusters.tsv compared as member SETS (cluster ids are ranks, so they renumber).

    pair_diff.py BEFORE.tsv AFTER.tsv OUT.json [GRAPH_TSV]

clusters.tsv lists clusters with >= min-size (2) members, one row per member. Allowed: a 2-member cluster in AFTER
whose members are absent from BEFORE (two singletons joined). Every other difference — a member set present on one side
only, or the same set with other size/density/frac_in/corroborated values — is an `other_change`; a 2-member set of
BEFORE whose members are absent from AFTER is a `pair_split`.

With GRAPH_TSV (`mcl_families --dump-graph`) the verdict also requires every joined pair to be an ISOLATED two-node
component of the dump (`not_isolated` 0). The exit status is the verdict: 0 on PASS, 1 on FAIL.
"""
import collections
import json
import sys

VALUES = ("size", "density", "frac_in", "corroborated")


def load(path):
    """-> {frozenset(member 'chrom:start-end'): (size, density, frac_in, corroborated)}."""
    members, vals = collections.defaultdict(list), {}
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ix = {k: hdr.index(k) for k in ("cluster_id", "chrom", "start", "end") + VALUES}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            cid = f[ix["cluster_id"]]
            members[cid].append(f"{f[ix['chrom']]}:{f[ix['start']]}-{f[ix['end']]}")
            vals[cid] = tuple(f[ix[k]] for k in VALUES)
    return {frozenset(m): vals[c] for c, m in members.items()}


def diff(before, after):
    in_before = frozenset().union(*before) if before else frozenset()
    in_after = frozenset().union(*after) if after else frozenset()
    r = collections.Counter(unchanged=0, pairs_joined=0, pairs_split=0, other_changes=0)
    for s, v in after.items():
        if s in before:
            r["unchanged" if before[s] == v else "other_changes"] += 1
        elif len(s) == 2 and not s & in_before:
            r["pairs_joined"] += 1
        else:
            r["other_changes"] += 1
    for s in before:
        if s not in after:
            r["pairs_split" if len(s) == 2 and not s & in_after else "other_changes"] += 1
    r["clusters_before"], r["clusters_after"] = len(before), len(after)
    r["members_before"], r["members_after"] = len(in_before), len(in_after)
    return dict(r)


W_STAR = 0.456885  # inflation 2.8: exact MCL splits a lone pair below this weight; above it a tie that float residue decided


def edge_weights(graph_tsv):
    """`mcl_families --dump-graph` (`u<TAB>v<TAB>weight` per pre-MCL edge) -> {frozenset({u, v}): weight}."""
    out = {}
    with open(graph_tsv) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) == 3:
                out[frozenset(f[:2])] = float(f[2])
    return out


def joined_bands(before, after, weights, w_star=W_STAR):
    """Joined pairs by edge weight: `tie_regime` (w >= w_star: float residue decided them) and `below_w_star` (exact MCL
    split them; the pre-registered rule joins them too). `not_isolated` (always present) counts, ON TOP of the bands, the
    joined pairs with an endpoint that has any other edge in the dump (degree > 1): the fix joins isolated two-node
    components only."""
    in_before = frozenset().union(*before) if before else frozenset()
    degree = collections.Counter(x for e in weights for x in e)
    c = collections.Counter(not_isolated=0)
    for s in after:
        if s not in before and len(s) == 2 and not s & in_before:
            w = weights.get(s)
            c["no_weight" if w is None else ("tie_regime" if w >= w_star else "below_w_star")] += 1
            if any(degree[x] > 1 for x in s):
                c["not_isolated"] += 1
    return dict(sorted(c.items()))


def symmetric_cliques(weights):
    """Components that are complete graphs of >= 3 nodes with one edge weight < 1 (equal to 6 dp): MCL can converge to an
    exact tie there too, which the pair rule does not cover (a uniform clique at weight 1.0 is exactly uniform, so safe).
    -> {"components": n, "members": m}."""
    parent = {}

    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for e in weights:
        a, b = sorted(e)
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[max(ra, rb)] = min(ra, rb)
    nodes, ws, n_edges = collections.defaultdict(set), collections.defaultdict(set), collections.Counter()
    for e, w in weights.items():
        r = find(min(e))
        nodes[r].update(e)
        ws[r].add(w)
        n_edges[r] += 1
    sym = [r for r in nodes if len(nodes[r]) >= 3 and n_edges[r] == len(nodes[r]) * (len(nodes[r]) - 1) // 2
           and len(ws[r]) == 1 and min(ws[r]) < 1.0]
    return {"components": len(sym), "members": sum(len(nodes[r]) for r in sym)}


def main(argv=None):
    a = argv if argv is not None else sys.argv[1:]
    before, after = load(a[0]), load(a[1])
    res = diff(before, after)
    if len(a) > 3:
        w = edge_weights(a[3])
        res["joined_by_weight"] = joined_bands(before, after, w)
        res["uniform_cliques_below_1"] = symmetric_cliques(w)
    isolated = res.get("joined_by_weight", {}).get("not_isolated", 0) == 0
    res["verdict"] = "PASS" if res["pairs_split"] == 0 and res["other_changes"] == 0 and isolated else "FAIL"
    with open(a[2], "w") as fh:
        json.dump(res, fh, indent=1)
    print(json.dumps(res))
    return 0 if res["verdict"] == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
