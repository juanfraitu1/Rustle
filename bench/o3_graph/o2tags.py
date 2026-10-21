#!/usr/bin/env python3
"""Which BAM tag or CIGAR statistic picks a tied read's true copy (user 2026-10-20: "how about the analysis of CIGAR and
bam tags such as md, nm, de, ms etc").   o2tags.py BAM COPIES FAMILY CALLS.tsv TRUTH.tsv..."""
import collections
import csv
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def sib(n):
    s = importlib.util.spec_from_file_location(f"o3_graph_{n}_tags", os.path.join(HERE, f"{n}.py"))
    m = importlib.util.module_from_spec(s)
    s.loader.exec_module(m)
    return m


def stats(rec):
    """The statistics of one alignment: AS, ms (intron-free), NM, de, mismatches (NM less indel bases), indel bases,
    clipped read bases, NM + clipped."""
    cig = rec.cigartuples
    ins = sum(n for o, n in cig if o == 1)
    dele = sum(n for o, n in cig if o == 2)
    clip = sum(n for o, n in cig if o in (4, 5))
    nm = rec.get_tag("NM")
    return dict(AS=rec.get_tag("AS"), ms=rec.get_tag("ms"), NM=nm, de=rec.get_tag("de"), mismatch=nm - ins - dele,
                indel=ins + dele, clip=clip, NM_clip=nm + clip)


HIGHER = {"AS", "ms"}                                           # the rest: lower is better


def pick(per, key):
    """The copy a statistic ranks first ({copy: stats}), or None on a tie."""
    v = {c: s[key] for c, s in per.items()}
    best = max(v.values()) if key in HIGHER else min(v.values())
    top = [c for c, x in v.items() if x == best]
    return top[0] if len(top) == 1 else None
def main(bam_p, copies_p, family, calls_p, *truth_ps):
    import pysam
    P, S = sib("o2psv"), sib("summary")
    cps = P.load_copies(copies_p)
    truth = {}
    for p in truth_ps:
        for r in csv.DictReader(open(p), delimiter="\t"):
            if r["class"] in ("unique", "unique_elsewhere"):
                truth[r["read"]] = r["truth"]
    calls = {r["read"]: r["corrector"] for r in csv.DictReader(open(calls_p), delimiter="\t")}
    spans = [(ch, s, e) for ch, v in cps.items() for c, _, _, (s, e) in v if c.rsplit(".", 1)[0] == family]
    best, prim = collections.defaultdict(dict), {}
    with pysam.AlignmentFile(bam_p) as bam:
        for ch, s, e in spans:
            for rec in bam.fetch(ch, s, e):
                if rec.is_unmapped or rec.is_supplementary:
                    continue
                c = S.attribute(rec.get_blocks(), P._strand(rec), cps.get(ch, ()))
                if not c or c.rsplit(".", 1)[0] != family:
                    continue
                st = stats(rec)
                if c not in best[rec.query_name] or st["AS"] > best[rec.query_name][c]["AS"]:
                    best[rec.query_name][c] = st
                if not rec.is_secondary:
                    prim[rec.query_name] = c
    keys = ("AS", "ms", "NM", "de", "mismatch", "indel", "clip", "NM_clip")
    n = collections.defaultdict(collections.Counter)
    for q, c in calls.items():                                  # the tied reads (window 30), both copies aligned
        if q not in truth or q not in prim or len(best[q]) < 2:
            continue
        group = "decided" if c else "abstained"
        t = truth[q]
        for g in (group, "all tied"):
            n[g]["reads"] += 1
            n[g]["primary right"] += prim[q] == t
            for k in keys:
                p = pick(best[q], k)
                n[g][f"{k} decides"] += p is not None
                n[g][f"{k} right"] += p == t
                n[g][f"{k} wrong"] += p is not None and p != t
    for g in ("all tied", "decided", "abstained"):
        c = n[g]
        print(f"== {g}: {c['reads']} reads, minimap2's primary right {c['primary right']}")
        for k in keys:
            print(f"  {k:9s} decides {c[k + ' decides']:5d}  right {c[k + ' right']:5d}  wrong {c[k + ' wrong']:4d}")
    return 0


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:]))
