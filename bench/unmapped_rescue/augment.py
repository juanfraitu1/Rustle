#!/usr/bin/env python3
"""Augmented-reference realignment (docs/PREREG_unmapped_rescue_2026-10-08.md Amendments 7 and 8). Pure helpers; pipelines in run_augment.py / run_lrpap1.py.

A record is a dict(score, de, qcov, ref, mapq). The genome primary comes from the existing alignment to the masked genome, the consensus record from aligning
the same read with the same minimap2 command to the consensus sequences alone."""
import statistics

QCOV_MIN = 0.80      # the registered hit-coverage floor


def new_primary(genome, cons):
    """(source, record): the primary a concatenated genome + consensus reference would give: the consensus record iff its alignment score is strictly
    higher than the genome primary's; an unmapped read takes a consensus record with query coverage >= 0.80; ties stay on the genome"""
    if genome is None:
        return ("consensus", cons) if cons is not None and cons["qcov"] >= QCOV_MIN else ("none", None)
    if cons is not None and cons["score"] > genome["score"]:
        return "consensus", cons
    return "genome", genome


def move_metrics(rows, cons_family, min_identity=0.98):
    """rows: [(class, true family, genome record | None, consensus record | None, consensus name | None)] with class in D_unm / D_abs / S / bg.
    -> {class: dict(n, moved, moved_to_own_family, move_fraction, median_de_before_moved, median_de_after_moved)}. A move needs the consensus primary
    (and, for unmapped reads, de <= 1 - min_identity); 'own family' = the consensus' cluster majority family equals the read's family."""
    out = {}
    for cls, fam, g, c, cname in rows:
        o = out.setdefault(cls, dict(n=0, moved=0, moved_to_own_family=0, moved_to_a_family_cluster=0, _before=[], _after=[]))
        o["n"] += 1
        src, r = new_primary(g, c)
        if src != "consensus":
            continue
        if g is None and r["de"] > 1 - min_identity:
            continue
        o["moved"] += 1
        o["moved_to_own_family"] += cons_family.get(cname) == fam
        o["moved_to_a_family_cluster"] += cons_family.get(cname) is not None
        if g is not None:
            o["_before"].append(g["de"])
        o["_after"].append(r["de"])
    for o in out.values():
        o["move_fraction"] = o["moved"] / o["n"] if o["n"] else 0.0
        o["median_de_before_moved"] = statistics.median(o.pop("_before")) if o["_before"] else None
        o["median_de_after_moved"] = statistics.median(o.pop("_after")) if o["_after"] else None
    return out
