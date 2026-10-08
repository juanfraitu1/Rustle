#!/usr/bin/env python3
"""Pure helpers of the LRPAP1 worked example (docs/PREREG_unmapped_rescue_2026-10-08.md Amendment 8). Pipeline: run_lrpap1.py."""

PLACED = 0.999     # the registered recovery hit: identity x coverage


def intervals(rows, hap):
    """[(label, accession, start, end)] of the LRPAP1 loci on a haplotype ('mat' | 'pat') from lrpap1_loci.tsv rows. Loci with no site there are skipped;
    loci that lift to the same interval on that haplotype (c01 and c03 on mat) are one site labelled 'c01|c03'."""
    sites = {}
    for r in rows:
        acc = r[f"{hap}_acc"]
        if not acc:
            continue
        sites.setdefault((acc, int(r[f"{hap}_s"]), int(r[f"{hap}_e"])), []).append(r["cid"])
    return [("|".join(c), acc, s, e) for (acc, s, e), c in sites.items()]


def locus_at(iv, acc, start, end, min_frac=0.5):
    """the label of the single locus that holds at least min_frac of the span [start, end) on accession acc, else None"""
    span = end - start
    if span <= 0:
        return None
    hit = [lab for lab, a, s, e in iv if a == acc and max(0, min(end, e) - max(start, s)) >= min_frac * span]
    return hit[0] if len(hit) == 1 else None


def group(cid):
    """scoring group of a copy: p12 and p14 differ by 4 bp in 2,446 on the paternal haplotype (tied), so they are one group"""
    if cid is None:
        return None
    return "p12|p14" if cid in ("p12", "p14") else cid


def placements(records, iv):
    """{locus label: best identity x coverage} of alignment records dict(ref, start, end, de, qcov) that fall on a locus"""
    out = {}
    for r in records:
        lab = locus_at(iv, r["ref"], r["start"], r["end"])
        if lab is None:
            continue
        v = (1 - r["de"]) * r["qcov"]
        out[lab] = max(out.get(lab, 0.0), v)
    return out


def o3_flag(ref_best, other_best, bar=PLACED):
    """the O3 reading of a consensus: it matches a locus of the other haplotype at the bar but nothing on the reference haplotype does"""
    return other_best is not None and other_best >= bar and (ref_best is None or ref_best < bar)
