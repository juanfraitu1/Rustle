#!/usr/bin/env python3
"""Post-MCL strongest-edge reassignment of family members (prereg docs/archive/2026-10/merged/PREREG_post_mcl_reassign_2026-10-15.md).

After MCL, a member's `in` is its strongest edge to another member of its own family and its `out` is, per other
family, its strongest edge into that family. A member whose best `out` is strictly greater than its `in` moves to that
family (equal weights between two outside families: the smaller family name). Members are visited in sorted order and
see the moves already made in the pass (asynchronous); passes repeat until one moves nothing (or `max_passes`). Only
members move: an edge to a node in no family is ignored, and no node joins or leaves the catalog. A member that ends
alone in its family is reported as lost (the catalog's minimum family size is 2).
"""
import collections
import csv


def read_dump(path):
    """`mcl_families --dump-graph` (u<TAB>v<TAB>weight per edge) -> {node: [(neighbour, weight)]}."""
    adj = collections.defaultdict(list)
    with open(path) as fh:
        for line in fh:
            u, v, w = line.rstrip("\n").split("\t")
            adj[u].append((v, float(w)))
            adj[v].append((u, float(w)))
    return dict(adj)


def read_clusters(path):
    """clusters.tsv -> {node: family}; nodes are `chrom:start-end` as in the dump and the PAF."""
    with open(path) as fh:
        return {f"{r['chrom']}:{r['start']}-{r['end']}": r["cluster_id"] for r in csv.DictReader(fh, delimiter="\t")}


def _best(v, fam, adj):
    """-> (in, best outside family or None, its out)."""
    f = fam[v]
    w_in, out = 0.0, {}
    for u, w in adj.get(v, ()):
        g = fam.get(u)
        if g is None or u == v:
            continue
        if g == f:
            w_in = max(w_in, w)
        elif w > out.get(g, -1.0):
            out[g] = w
    if not out:
        return w_in, None, 0.0
    g, w_out = min(out.items(), key=lambda kv: (-kv[1], kv[0]))
    return w_in, g, w_out


def reassign(fam_of, adj, max_passes=50):
    """-> (new {node: family}, log): log has moves [(pass, node, from, to, in, out)], passes, converged and lost (the
    members alone in their family at the end)."""
    fam = dict(fam_of)
    moves, passes, converged = [], 0, False
    order = sorted(fam)
    while passes < max_passes:
        passes += 1
        moved = 0
        for v in order:
            w_in, g, w_out = _best(v, fam, adj)
            if g is not None and w_out > w_in:
                moves.append((passes, v, fam[v], g, w_in, w_out))
                fam[v] = g
                moved += 1
        if not moved:
            converged = True
            break
    size = collections.Counter(fam.values())
    lost = sorted(v for v, f in fam.items() if size[f] == 1)
    return fam, dict(moves=moves, passes=passes, converged=converged, lost=lost)


def label_moves(moved, old, new, labels_of):
    """Per moved member with labels (e.g. Soto families): the members sharing one of its labels in its old family vs in
    its new one -> counts of gain (more in the new), loss (fewer), same, unlabelled."""
    out = dict(gain=0, loss=0, same=0, unlabelled=0)
    by_label = collections.defaultdict(set)
    for m, labs in labels_of.items():
        for lab in labs:
            by_label[lab].add(m)
    for m in moved:
        labs = labels_of.get(m)
        if not labs:
            out["unlabelled"] += 1
            continue
        mates = set().union(*(by_label[lab] for lab in labs)) - {m}
        n_old = sum(1 for x in mates if old.get(x) is not None and old.get(x) == old[m])
        n_new = sum(1 for x in mates if new.get(x) is not None and new.get(x) == new[m])
        out["gain" if n_new > n_old else "loss" if n_new < n_old else "same"] += 1
    return out


def pair_changes(before, after):
    """Two `_liftoff.pair_families` results over the same pairs -> joined (split before, shared after), split (the
    reverse), kept, apart (split in both), and the covered counts before and after."""
    out = dict(joined=0, split=0, kept=0, apart=0, covered_before=0, covered_after=0)
    for (_, b0, s0), (_, b1, s1) in zip(before, after):
        out["covered_before"] += b0
        out["covered_after"] += b1
        if not (b0 and b1):
            continue
        key = "kept" if s0 and s1 else "joined" if s1 else "split" if s0 else "apart"
        out[key] += 1
    return out
