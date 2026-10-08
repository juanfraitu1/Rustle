#!/usr/bin/env python3
"""Metrics of docs/PREREG_unmapped_rescue_2026-10-08.md section 4 (pure functions)."""
import collections


def majority(reads, truth):
    """most common truth family among the reads with a truth family (not None, not 'bg'); None if there is none"""
    c = collections.Counter(truth[r] for r in reads if truth.get(r) not in (None, "bg"))
    return c.most_common(1)[0][0] if c else None


def cluster_metrics(clusters, truth, d_reads):
    """clusters: {id: [reads]} (>= min size already applied); truth: {read: family | 'bg' | None}; d_reads: the deleted-copy reads (set).
    coverage = D reads in a cluster / D reads; purity = clustered D reads whose cluster majority is their own family / clustered D reads;
    clusters_per_family = clusters whose majority is the family; background_in_clusters = 'bg' reads in any cluster."""
    in_cluster = {}
    for cid, rs in clusters.items():
        for r in rs:
            in_cluster[r] = cid
    maj = {cid: majority(rs, truth) for cid, rs in clusters.items()}
    d_cl = [r for r in d_reads if r in in_cluster]
    pure = sum(1 for r in d_cl if maj[in_cluster[r]] == truth[r])
    per = collections.Counter(m for m in maj.values() if m)
    return dict(n_d=len(d_reads), d_clustered=len(d_cl), coverage=len(d_cl) / len(d_reads) if d_reads else 0.0,
                purity=pure / len(d_cl) if d_cl else None, clusters=len(clusters), clusters_per_family=dict(per),
                background_in_clusters=sum(1 for r in in_cluster if truth.get(r) == "bg"))


def rescue_metrics(clusters, attribution, truth, d_reads, other_roles=frozenset()):
    """clusters: {id: [reads]}; attribution: {cluster id: family | None (abstain)}; truth as in cluster_metrics; d_reads: the deleted-copy reads.
    rescued_correct = D reads of an attributed cluster attributed to their own family; wrong_joins = D reads attributed to another family plus
    background reads in an attributed cluster; unknown_joined = reads without a truth label, or in `other_roles` (e.g. surviving-copy reads), in an attributed
    cluster (reported, not judged);
    cluster_accuracy = attributed clusters whose majority family is the attributed one; copies_reached = families with >= 1 rescued-correct read."""
    ok = wrong = unknown = abstained = 0
    reached = set()
    n_attr = n_abs = n_right = 0
    for cid, rs in clusters.items():
        fam = attribution.get(cid)
        if fam is None:
            n_abs += 1
            abstained += sum(1 for r in rs if r in d_reads)
            continue
        n_attr += 1
        n_right += majority(rs, truth) == fam
        for r in rs:
            t = truth.get(r)
            if r in d_reads:
                if t == fam:
                    ok += 1
                    reached.add(fam)
                else:
                    wrong += 1
            elif t == "bg":
                wrong += 1
            elif t is None or r in other_roles:
                unknown += 1
    return dict(rescued_correct=ok, wrong_joins=wrong, unknown_joined=unknown, abstained_reads=abstained, clusters_attributed=n_attr,
                clusters_abstained=n_abs, cluster_accuracy=(n_right / n_attr) if n_attr else None, copies_reached=len(reached),
                copies_with_unmapped_reads=len({truth[r] for r in d_reads}))
