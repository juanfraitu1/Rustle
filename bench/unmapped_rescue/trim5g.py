#!/usr/bin/env python3
"""Trimming rule T for the untemplated 5' G run of a cluster consensus (docs/PREREG_unmapped_rescue_2026-10-08.md Amendment 13).

best_copy_prefix: how many consensus bases lie before the first dc-megablast HSP on the attributed family's best-covering surviving copy.
decide: remove that prefix iff it is 1 to MAX_G bases, all G."""
import attribute as A

MAX_G = 3


def best_copy_prefix(hsps, family):
    """hsps: [(query, target, qstart, qend)] (1-based, as A.read_blastn) of one cluster. -> prefix length (bases before the first HSP query start on the family's
    copy with the largest union of HSP query spans), or None if the family has no HSP"""
    by = {}
    for q, t, qs, qe in hsps:
        if A.family_of(t) == family:
            by.setdefault(t, []).append((min(qs, qe), max(qs, qe)))
    if not by:
        return None
    best = max(by, key=lambda t: (A.union_length(by[t]), t))
    return min(a for a, _ in by[best]) - 1


def decide(cons, prefix, max_len=MAX_G):
    """-> (consensus, trimmed length, reason) with reason in trimmed | no_hsp | no_prefix | too_long | not_g"""
    if prefix is None:
        return cons, 0, "no_hsp"
    if prefix == 0:
        return cons, 0, "no_prefix"
    if prefix > max_len:
        return cons, 0, "too_long"
    if set(cons[:prefix]) != {"G"}:
        return cons, 0, "not_g"
    return cons[prefix:], prefix, "trimmed"


def prefixes_from_paf(paf_lines, families):
    """Rule T2 (Amendment 14). paf_lines: consensus aligned to the surviving copies' genomic spans (names FAM:copy); families: {cluster key: attributed family or None}.
    -> {cluster key: consensus bases before the best alignment (largest `matches`) on a copy of its family, or None if it has none}; clusters with no family are skipped"""
    best = {}
    for ln in paf_lines:
        f = ln.rstrip("\n").split("\t")
        k = f[0].split("|")[0]
        fam = families.get(k)
        if fam is None or A.family_of(f[5]) != fam:
            continue
        m = int(f[9])
        if k not in best or m > best[k][0]:
            best[k] = (m, int(f[2]))
    return {k: (best[k][1] if k in best else None) for k, fam in families.items() if fam is not None}


def trim_leading_g(cons, max_len=MAX_G):
    """Rule T3 (Amendment 15): remove the maximal leading run of G if it is 1 to max_len long. -> (consensus, trimmed length, reason: trimmed | no_run | too_long)"""
    n = len(cons) - len(cons.lstrip("G"))
    if n == 0:
        return cons, 0, "no_run"
    if n > max_len:
        return cons, 0, "too_long"
    return cons[n:], n, "trimmed"
