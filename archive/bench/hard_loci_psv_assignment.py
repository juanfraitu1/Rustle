#!/usr/bin/env python3
"""Is DSFAM42 (a near-identical floor case) typical, or do PSVs actually ASSIGN reads in other hard
loci? Characterise every co-located multi-copy family on two axes:
  difficulty  = MAPQ-0 fraction (how badly minimap2 ties the reads)
  PSV success = fraction of reads PSV+junction confidently assigns, and the unique-mapper agreement.
Classify: WIN (hard AND PSV assigns), FLOOR (hard AND can't, like DSFAM42), EASY (minimap2 already ok)."""
import sys, json, time
sys.path.insert(0, "bench")
import copy_assign as CA
import pysam

meta = CA.load_meta(); skel = CA.load_skel(); spliced = CA.load_spliced()
resc = CA.load_rescued(meta, skel, spliced)
fams = CA.colocated_families(meta, rescued_fam=resc)
bam = pysam.AlignmentFile(CA.GGO_BAM, "rb")
rows = []
t0 = time.time()
for i, (fid, dom, copies) in enumerate(fams):
    if i % 10 == 0:
        print(f"  {i}/{len(fams)} ({time.time()-t0:.0f}s)", flush=True)
    try:
        st = CA.assign_family(spliced, bam, copies, meta, skel)
    except Exception as e:
        continue
    n = st["n"]
    if n < 20:
        continue
    mq0 = st["mq0"] / n
    assign_rate = st["assigned_j"] / n
    resolv_rate = st["resolvable_j"] / n
    val = (st["uniq_agree_j"] / st["uniq_j"]) if st["uniq_j"] else float("nan")
    rows.append(dict(family=fid, dom=dom, n_copies=len(copies), reads=n, psv_cols=st["psv_cols"],
                     mq0_frac=round(mq0, 3), assign_rate=round(assign_rate, 3),
                     resolv_rate=round(resolv_rate, 3), uniq_agreement=round(val, 3) if val == val else None,
                     uniq=st["uniq_j"]))

HARD = 0.5   # >=50% MAPQ-0 = a genuinely hard (multimapper-tied) locus
def cls(r):
    if r["mq0_frac"] < HARD: return "EASY (minimap2 ok)"
    return "WIN (hard, PSV assigns)" if r["assign_rate"] >= 0.5 else "FLOOR (hard, PSV can't)"
for r in rows: r["class"] = cls(r)
from collections import Counter
rows.sort(key=lambda r: (-r["mq0_frac"], -r["assign_rate"]))
print(f"\n=== {len(rows)} families (>=20 reads) ===")
print("class tally:", dict(Counter(r["class"] for r in rows)))
hard = [r for r in rows if r["mq0_frac"] >= HARD]
wins = [r for r in hard if r["assign_rate"] >= 0.5]
print(f"\nHARD loci (>={int(HARD*100)}% MAPQ-0): {len(hard)}  |  of those PSV WINS (>=50% assigned): {len(wins)}")
print(f"\n{'family':10s} {'dom':9s} {'cp':>3} {'reads':>6} {'MQ0%':>5} {'PSVc':>5} {'assign%':>7} {'resolv%':>7} {'uniq_agr':>8} class")
for r in hard:
    print(f"{r['family']:10s} {str(r['dom'])[:9]:9s} {r['n_copies']:>3} {r['reads']:>6} "
          f"{100*r['mq0_frac']:>4.0f}% {r['psv_cols']:>5} {100*r['assign_rate']:>6.0f}% "
          f"{100*r['resolv_rate']:>6.0f}% {(r['uniq_agreement'] or 0):>8.2f} {r['class']}")
json.dump(rows, open("bench/hard_loci_psv_assignment.json", "w"), indent=1)
print("\nwrote bench/hard_loci_psv_assignment.json")
