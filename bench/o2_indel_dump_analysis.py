#!/usr/bin/env python3
"""Indel-column diagnostics from --dump-star: for every molecule with >=1 indel column ('0'/'1' alleles),
does the read's indel allele ('0') agree with its SUBSTITUTION-best copy (bk in the flag-off base run)?"""
import csv, sys, statistics
from collections import Counter
base = {r["read_name"]: r for r in csv.DictReader(open(sys.argv[1]), delimiter="\t")}   # ours_odi.assignments.tsv
dump = sys.argv[2]                                                                        # ours_indelA_dump.star_reads.tsv
changed = set(l.strip() for l in open(sys.argv[3])) if len(sys.argv) > 3 else set()
agree = Counter(); per_mol = []; ncols = Counter(); alleles_per_col = Counter()
detail = []
for r in csv.DictReader(open(dump), delimiter="\t"):
    cands = r["candidates"].split(",")
    cols = [c for c in r["columns"].split(",") if c] if r["columns"] else []
    ind = [c for c in cols if c.split(":")[1] == "0" and set(c.split(":")[2]) <= set("01.")]
    if not ind: continue
    ncols[len(ind)] += 1
    b = base.get(r["read_name"])
    if not b: continue
    bk = b["assigned_copy"]                    # substitution-only best (internal unit id == cand index space)
    if bk not in cands: agree["bk_not_cand"] += 1; continue
    ki = cands.index(bk)
    a_ok = a_bad = 0
    for c in ind:
        al = c.split(":")[2]
        alleles_per_col[(al.count("1"), al.count("0"))] += 1
        if al[ki] == "0": a_ok += 1
        elif al[ki] == "1": a_bad += 1
    agree["cols_bk_no_gap"] += a_ok; agree["cols_bk_needs_gap"] += a_bad
    per_mol.append((r["read_name"], b["status"], a_ok, a_bad))
    if r["read_name"] in changed: detail.append((r["read_name"].split("/")[-1], b["status"], b["assigned_copy"], r["status"], r["assigned_copy"], ind))
print("molecules with indel columns:", sum(ncols.values()), "columns/mol:", sorted(ncols.items()))
print("indel columns where the substitution-best copy needs NO gap (consistent):", agree["cols_bk_no_gap"], " needs a gap (the read carries another copy's indel allele):", agree["cols_bk_needs_gap"])
bym = Counter()
for n, st, ok, bad in per_mol: bym[(st, "consistent" if bad == 0 else ("mixed" if ok else "inconsistent"))] += 1
print("per molecule (base status, indel verdict vs substitution-best):", sorted(bym.items()))
print("('1' count, '0' count) per column:", alleles_per_col.most_common(10))
for d in detail[:40]: print("  CHANGED", d[0], "base", d[1], d[2], "-> new", d[3], d[4], "|", " ".join(d[5]))
