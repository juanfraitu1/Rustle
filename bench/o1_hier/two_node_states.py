#!/usr/bin/env python3
"""PREREG_mcl_two_node_2026-10-13 F2, post hoc and descriptive: what the two-node fix changed in each catalog.

    two_node_states.py BEFORE_DIR AFTER_DIR GRAPH_DIR OUT.json [SAMPLE ...]

BEFORE_DIR / AFTER_DIR hold `<S>.fam.clusters.tsv` and `<S>.fam.loci.tsv` of the pre-fix and the fixed `mcl_families`
run (`recluster.sh`); GRAPH_DIR holds `<S>.graph.tsv` (`--dump-graph`: the pre-MCL edges). SAMPLE defaults to the four
default catalogs. Per sample, from those files only:

  * the two-node components of the graph (two nodes, one edge) by band, `tie_regime` (w >= w*) and `below_w_star`, each
    with its fates "<state before>-><state after>". A pair is `listed` (a 2-member family in clusters.tsv), `folded`
    (one locus: loci.tsv maps one record onto the other) or `split` (neither);
  * the `loci.tsv` line changes: lines added, lines removed or changed, and where the added ones lie;
  * the uniform cliques below 1 (complete graphs of >= 3 nodes whose edges share one weight < 1, the definition of
    `pair_diff.symmetric_cliques`) with their fates: `listed` (one cluster of exactly those nodes), `folded` (all in one
    locus) or `other`.

`total` sums the samples. The fix is meant to change only two-node components, so a clean run has `still_split_after` 0,
`loci_tsv.removed_or_changed` 0 and `loci_tsv.added_elsewhere` 0.
"""
import collections
import json
import sys

import pair_diff as P

SAMPLES = ("human_A119b", "human_testis", "chimp_PTR", "gorilla_KB3781")
BANDS = ("tie_regime", "below_w_star")


def load_loci(path):
    """`<S>.fam.loci.tsv` (`annotation<TAB>representative` per folded record) -> {annotation: representative}."""
    out = {}
    with open(path) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) == 2:
                out[f[0]] = f[1]
    return out


def components(weights):
    """Connected components of the graph `weights` ({frozenset({u, v}): w}) -> [(nodes, number of edges, set of edge
    weights)], ordered by the smallest node name."""
    parent = {}

    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for e in weights:
        root = find(min(e))
        for x in e:
            rx = find(x)
            if rx != root:
                parent[max(root, rx)] = min(root, rx)
                root = min(root, rx)
    nodes, ws, n_edges = collections.defaultdict(set), collections.defaultdict(set), collections.Counter()
    for e, w in weights.items():
        r = find(min(e))
        nodes[r].update(e)
        ws[r].add(w)
        n_edges[r] += 1
    return [(frozenset(nodes[r]), n_edges[r], ws[r]) for r in sorted(nodes)]


def band(w, w_star):
    return "tie_regime" if w >= w_star else "below_w_star"


def pair_state(pair, clusters, loci):
    """`listed`: a 2-member family; `folded`: loci.tsv puts both records in one locus; else `split`."""
    if pair in clusters:
        return "listed"
    a, b = sorted(pair)
    return "folded" if loci.get(a) == b or loci.get(b) == a else "split"


def clique_state(nodes, clusters, loci):
    """`listed`: one cluster of exactly these nodes; `folded`: all of them in one locus; else `other`."""
    if nodes in clusters:
        return "listed"
    reps = {loci.get(x, x) for x in nodes}
    return "folded" if len(reps) == 1 and next(iter(reps)) in nodes else "other"


def states(before, after, loci_before, loci_after, weights, w_star=P.W_STAR):
    """One sample's accounting. before/after: {frozenset(members): values} (`pair_diff.load`); loci_*: {annotation:
    representative} (`load_loci`); weights: {frozenset({u, v}): w} (`pair_diff.edge_weights`)."""
    comps = components(weights)
    pairs = {nodes: next(iter(ws)) for nodes, n_edges, ws in comps if len(nodes) == 2 and n_edges == 1}
    fates = {b: collections.Counter() for b in BANDS}
    for pr, w in pairs.items():
        fates[band(w, w_star)][f"{pair_state(pr, before, loci_before)}->{pair_state(pr, after, loci_after)}"] += 1
    added = {k: v for k, v in loci_after.items() if loci_before.get(k) != v}
    removed = [k for k, v in loci_before.items() if loci_after.get(k) != v]
    inside, elsewhere = dict.fromkeys(BANDS, 0), 0
    for k, v in added.items():
        pr = frozenset((k, v))
        if pr in pairs:
            inside[band(pairs[pr], w_star)] += 1
        else:
            elsewhere += 1
    cliques = []
    for nodes, n_edges, ws in comps:
        n = len(nodes)
        if n >= 3 and n_edges == n * (n - 1) // 2 and len(ws) == 1 and min(ws) < 1.0:
            cliques.append({"members": sorted(nodes), "weight": next(iter(ws)),
                            "before": clique_state(nodes, before, loci_before),
                            "after": clique_state(nodes, after, loci_after)})
    return {
        "two_node_components": len(pairs),
        "by_band": {b: {"total": sum(fates[b].values()), "fates": dict(sorted(fates[b].items()))}
                    for b in BANDS},
        "still_split_after": sum(n for b in BANDS for f, n in fates[b].items() if f.endswith("->split")),
        "loci_tsv": {"lines_before": len(loci_before), "lines_after": len(loci_after), "added": len(added),
                     "removed_or_changed": len(removed), "added_in_two_node_component": inside,
                     "added_elsewhere": elsewhere},
        "uniform_cliques_below_1": cliques,
    }


def total(per_sample):
    """Sum of the per-sample records."""
    by_band = {b: {"total": 0, "fates": collections.Counter()} for b in BANDS}
    loci, inside = collections.Counter(), collections.Counter()
    fates_before, fates_after = collections.Counter(), collections.Counter()
    n_pairs = still = n_cliques = n_members = 0
    for r in per_sample.values():
        n_pairs += r["two_node_components"]
        still += r["still_split_after"]
        for b in BANDS:
            by_band[b]["total"] += r["by_band"][b]["total"]
            by_band[b]["fates"].update(r["by_band"][b]["fates"])
        for k in ("added", "removed_or_changed", "added_elsewhere"):
            loci[k] += r["loci_tsv"][k]
        inside.update(r["loci_tsv"]["added_in_two_node_component"])
        for c in r["uniform_cliques_below_1"]:
            n_cliques += 1
            n_members += len(c["members"])
            fates_before[c["before"]] += 1
            fates_after[c["after"]] += 1
    return {
        "two_node_components": n_pairs,
        "by_band": {b: {"total": by_band[b]["total"], "fates": dict(sorted(by_band[b]["fates"].items()))}
                    for b in BANDS},
        "still_split_after": still,
        "loci_tsv": {"added": loci["added"], "removed_or_changed": loci["removed_or_changed"],
                     "added_elsewhere": loci["added_elsewhere"],
                     "added_in_two_node_component": {b: inside[b] for b in BANDS}},
        "uniform_cliques_below_1": {"components": n_cliques, "members": n_members,
                                    "fates_before": dict(sorted(fates_before.items())),
                                    "fates_after": dict(sorted(fates_after.items()))},
    }


def main(argv=None):
    a = argv if argv is not None else sys.argv[1:]
    if len(a) < 4:
        sys.exit(__doc__)
    before_dir, after_dir, graph_dir, out = a[:4]
    res = {"samples": {}}
    for s in tuple(a[4:]) or SAMPLES:
        weights = P.edge_weights(f"{graph_dir}/{s}.graph.tsv")
        r = states(P.load(f"{before_dir}/{s}.fam.clusters.tsv"), P.load(f"{after_dir}/{s}.fam.clusters.tsv"),
                   load_loci(f"{before_dir}/{s}.fam.loci.tsv"), load_loci(f"{after_dir}/{s}.fam.loci.tsv"), weights)
        k = r["uniform_cliques_below_1"]
        if P.symmetric_cliques(weights) != {"components": len(k), "members": sum(len(c["members"]) for c in k)}:
            raise RuntimeError(f"{s}: the uniform-clique count differs from pair_diff.symmetric_cliques")
        res["samples"][s] = r
        b = r["by_band"]
        print(f"{s}: {r['two_node_components']} two-node components (tie_regime {b['tie_regime']['total']}, "
              f"below_w_star {b['below_w_star']['total']}); below_w_star fates {b['below_w_star']['fates']}; "
              f"still split after {r['still_split_after']}; loci.tsv +{r['loci_tsv']['added']} "
              f"-{r['loci_tsv']['removed_or_changed']}; uniform cliques below 1: {len(k)}")
    res["total"] = total(res["samples"])
    with open(out, "w") as fh:
        json.dump(res, fh, indent=1)
    t = res["total"]
    tb = t["by_band"]
    print(f"total: {t['two_node_components']} two-node components (tie_regime {tb['tie_regime']['total']}, "
          f"below_w_star {tb['below_w_star']['total']}); below_w_star fates {tb['below_w_star']['fates']}; "
          f"still split after {t['still_split_after']}; loci.tsv +{t['loci_tsv']['added']}; "
          f"uniform cliques below 1: {t['uniform_cliques_below_1']['components']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
