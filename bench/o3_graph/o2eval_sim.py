#!/usr/bin/env python3
"""Milestone 3 development: score the port's O2 assignments (PREFIX.o2.tsv) on a one-copy simulation run, where every read's
origin is exact (read name <copy>|<allele>_r<i>). Wrong = assigned to a copy that is not its origin (a read of a copy absent from
the copy table is right only on a certified new copy). Baseline: minimap2's primary on the same AS-tied reads.
    o2eval_sim.py RUN_DIR [RUN_DIR ...]
"""
import csv, os, sys, collections, pysam

def copy_at(where, ch, s, e):
    for c, (seg, a, b) in where.items():
        if seg == ch and s < b and a < e:
            return c
    return None

for D in sys.argv[1:]:
    copies = {r["tid"]: (r["chrom"], int(r["start"]), int(r["end"])) for r in csv.DictReader(open(f"{D}/copies.tsv"), delimiter="\t")}
    o2 = list(csv.DictReader(open(f"{D}/port.o2.tsv"), delimiter="\t"))
    prim, mq = {}, {}
    for rd in pysam.AlignmentFile(f"{D}/reads.bam"):
        if not (rd.is_secondary or rd.is_supplementary or rd.is_unmapped):
            prim[rd.query_name] = copy_at(copies, rd.reference_name, rd.reference_start, rd.reference_end)
            mq[rd.query_name] = rd.mapping_quality
    c = collections.Counter()
    for r in o2:
        truth = r["read"].split("|")[0]
        intable = truth in copies
        a = r["assigned"]
        m = prim.get(r["read"])
        mm_ok = m == truth
        if mq.get(r["read"], 0) >= 1:   # baseline 2: minimap2's confident calls (MAPQ >= 1) on the same tied reads
            c["mq1_placed"] += 1
            c["mq1_wrong"] += not mm_ok
        if a in ("ABSTAIN", "HOPELESS"):
            o = "abstain"
            c["hopeless"] += a == "HOPELESS"
        elif a.startswith("new:"):
            o = "right" if not intable else "wrong"
        else:
            o = "right" if a == truth else "wrong"
        c[o] += 1
        c["mm2_" + ("right" if mm_ok else "wrong")] += 1
        if o == "right" and not mm_ok: c["fix"] += 1
        if o == "wrong" and mm_ok: c["break"] += 1
    placed = c["right"] + c["wrong"]
    print(f"{os.path.basename(D)}: tied {len(o2)}; placed {placed} wrong {c['wrong']} ({c['wrong']/max(1,placed):.1%}); abstain {c['abstain']}; "
          f"fixes {c['fix']} breaks {c['break']}; minimap2 wrong {c['mm2_wrong']}/{len(o2)}; hopeless {c['hopeless']}; "
          f"minimap2 MAPQ>=1 placed {c['mq1_placed']} wrong {c['mq1_wrong']}")
