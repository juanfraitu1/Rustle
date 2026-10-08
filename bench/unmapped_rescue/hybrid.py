#!/usr/bin/env python3
"""Hybrid net (docs/PREREG_unmapped_rescue_2026-10-08.md Amendment 6): pooled cluster attribution first, per-read attribution for the residual."""


def residual(clusters, cluster_att, reads):
    """reads not rescued by an attributed cluster (abstaining clusters, reads in no cluster)"""
    done = {r for c, rs in clusters.items() if cluster_att.get(c) for r in rs}
    return [r for r in reads if r not in done]


def combine(clusters, cluster_att, read_att, reads):
    """{read: family | None}: the cluster's family when its cluster was attributed, else the read's own attribution"""
    out = {r: read_att.get(r) for r in reads}
    for c, rs in clusters.items():
        if cluster_att.get(c):
            for r in rs:
                out[r] = cluster_att[c]
    return out


def read_metrics(att, truth, d_reads):
    """read level: rescued_correct (D read attributed to its own family), wrong (D read attributed elsewhere, or a background read attributed to any
    family), unattributed_d"""
    ok = wrong = un = 0
    for r, f in att.items():
        t = truth.get(r)
        if r in d_reads:
            if f is None:
                un += 1
            elif f == t:
                ok += 1
            else:
                wrong += 1
        elif t == "bg" and f is not None:
            wrong += 1
    return dict(rescued_correct=ok, wrong=wrong, unattributed_d=un, wrong_fraction_of_joined=wrong / max(1, ok + wrong))
