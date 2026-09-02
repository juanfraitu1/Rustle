#!/usr/bin/env python3
"""Audit: does a family contain two 'copies' that overlap EACH OTHER?

Blind-spot audit item 6 (2026-08-14) recorded, and never re-checked, that
**30 of 194 two-copy families had a copy overlapping its own sister** -- a node-construction
defect, not a definition one. A family whose two copies are the same piece of genome is not
two copies, so the rate bounds how much of every copy count is real.

Pure coordinate truth: no aligner, no annotation, no gene label. Reports the rate over
FAMILIES and over within-family PAIRS with Wilson intervals, plus the overlap as a fraction
of the shorter copy so a 1 bp touch is not scored like a duplicate emission.

Usage: catalog_self_overlap_audit.py <label>=<copies.tsv> [...]
"""
import csv, sys, math, collections, itertools


def wilson(k, n, z=1.96):
    if n == 0:
        return (0.0, 1.0)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def audit(label, path):
    try:
        rows = list(csv.DictReader(open(path), delimiter="\t"))
    except OSError as e:
        print(f"{label}: UNREADABLE ({e})")
        return
    by = collections.defaultdict(list)
    for r in rows:
        by[r["family_id"]].append((r["chrom"], int(r["start"]), int(r["end"])))
    nfam = len(by)
    npair = nov = 0
    fam_ov = set()
    frac = []
    for fid, v in by.items():
        for (c1, s1, e1), (c2, s2, e2) in itertools.combinations(sorted(set(v)), 2):
            npair += 1
            if c1 == c2 and s1 < e2 and s2 < e1:
                nov += 1
                fam_ov.add(fid)
                frac.append((min(e1, e2) - max(s1, s2)) / min(e1 - s1, e2 - s2))
    lo, hi = wilson(len(fam_ov), nfam)
    plo, phi = wilson(nov, npair)
    print(f"\n=== {label} ===")
    print(f"  families={nfam}  copies={len(rows)}  within-family pairs={npair}")
    print(f"  families with a self-overlapping pair: {len(fam_ov)}/{nfam} = "
          f"{len(fam_ov)/max(1,nfam):.4f}  Wilson95 [{lo:.4f},{hi:.4f}]")
    print(f"  overlapping pairs: {nov}/{npair} = {nov/max(1,npair):.4f}  "
          f"Wilson95 [{plo:.4f},{phi:.4f}]")
    if frac:
        frac.sort()
        print(f"  overlap / shorter copy: median {frac[len(frac)//2]:.3f}  "
              f">=50%: {sum(1 for f in frac if f >= 0.5)}  >=90%: {sum(1 for f in frac if f >= 0.9)}")


if __name__ == "__main__":
    for a in sys.argv[1:]:
        label, path = a.split("=", 1)
        audit(label, path)
