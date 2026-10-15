#!/usr/bin/env python3
"""A tolerant reading of O3's "nothing in the genome explains it" (post hoc, descriptive; docs/PREREG_o3_family_best_2026-10-15.md,
Results): an `unexplained` cluster counts only when its reads' median best-alignment divergence (minimap2 `de` on their primary
record, i.e. to the nearest locus the aligner found anywhere in the genome) exceeds delta.

    tolerant.py added-back PREFIX BAM [LABELS]   `copy` clusters (est_div > delta) and their families: as run / tolerant
    tolerant.py removed PREFIX BAM LABELS        eligible masked copies recovered: as run / tolerant (the recovering
                                                 cluster must also be beyond delta), with the copies the tolerant rule loses
LABELS (role D) fixes the family set for added-back; without it every family with a cluster counts.
"""
import collections
import json
import os
import statistics
import sys

import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rd_lib as L  # noqa: E402
import score as S  # noqa: E402


def median_de(names, de):
    v = [de[n] for n in names if n in de]
    return statistics.median(v) if v else None


def beyond(names, de, delta=L.DELTA):
    m = median_de(names, de)
    return m is not None and m > delta


def primary_de(bam_p, want):
    out = {}
    with pysam.AlignmentFile(bam_p) as bam:
        for rec in bam.fetch(until_eof=True):
            if rec.query_name in want and not rec.flag & 0x904 and rec.has_tag("de"):
                out[rec.query_name] = rec.get_tag("de")
    return out


def main(argv):
    cmd, prefix, bam_p, *rest = argv
    placed, clusters, members = S.load(prefix)
    de = primary_de(bam_p, set().union(*members.values()) if members else set())
    mem = lambda c: members[(c["family"], c["host"], c["cluster"])]  # noqa: E731
    if cmd == "added-back":
        fams = set(S.absent_reads(rest[0])) if rest else None
        cp = [c for c in clusters if c["status"] == "unexplained" and float(c["est_div"]) > L.DELTA
              and (fams is None or c["family"] in fams)]
        tol = [c for c in cp if beyond(mem(c), de)]
        res = dict(copy_clusters=len(cp), copy_families=len({c["family"] for c in cp}), tolerant_clusters=len(tol),
                   tolerant_families=sorted({c["family"] for c in tol}))
    else:
        absent = S.absent_reads(rest[0])
        elig, rec, tol, lost = 0, 0, 0, []
        for fam, (copy, d) in sorted(absent.items()):
            p = placed.get(fam, set())
            if len(d & p) < S.MIN_PLACED:
                continue
            elig += 1
            s = L.score_absent(fam, d, p, clusters, members)
            if not s["recovered"]:
                continue
            rec += 1
            bc = max((c for c in clusters if c["family"] == fam and c["status"] == "unexplained"), key=lambda c: len(mem(c) & d & p))
            if beyond(mem(bc), de):
                tol += 1
            else:
                lost.append(dict(family=fam, median_de=round(median_de(mem(bc), de) or 0.0, 4), label=s["label"]))
        res = dict(eligible=elig, recovered=rec, tolerant_recovered=tol, lost=lost)
    print(json.dumps(res, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
