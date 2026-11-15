#!/usr/bin/env python3
"""O2 audit 2026-11-13 (measurement only): how far O2's decisions reach into the assembler's read pool, and how often
minimap2's primary is right on the near-tie strata.
  funnel MOLECULES.tsv ASSIGNMENTS.tsv COPIES.tsv   per molecule: good pool (second AS >= 0.98 x best), exact AS tie,
      ms within 10 -> in the assign table -> status -> sole -> assigned copy holds the primary (kept) or not (moved)
  premise ASSIGNMENTS.tsv DUMP [DUMP...]            decidable rows of a strict_truth dump by (AS margin, ms margin)
      stratum: primary right rate (dump al_right, leftmost-position rule), mean 1/n_candidates, O2 right/wrong"""
import collections, csv, sys

def funnel(molecules, assignments, copies):
    span = {}
    for r in csv.DictReader(open(copies), delimiter="\t"): span[(r["family_id"], r["copy_idx"])] = (r["chrom"], int(r["locus_start"]), int(r["locus_end"]))
    rank = {"assigned": 3, "tied": 2, "ambiguous": 1}; a = {}
    for r in csv.DictReader(open(assignments), delimiter="\t"):
        v = (rank[r["status"]], r["sole_candidate"] == "1", int(r["n_candidates"]), (r["family_id"], r["assigned_copy"]))
        if r["read_name"] not in a or v[:2] > a[r["read_name"]][:2] or (v[0] == 3 and not v[1] and a[r["read_name"]][1]): a[r["read_name"]] = v
    c = collections.Counter()
    for l in open(molecules):
        if l[0] == "#": continue
        f = l.rstrip("\n").split("\t"); c["molecules"] += 1
        if int(f[3]) < 2: continue
        b, s, bm, sm = int(f[1]), int(f[2]), int(f[6]), int(f[7])
        cls = [k for k, t in (("pool", s >= 0.98 * b), ("exactAS", s == b), ("ms10", sm >= 0 and bm - sm <= 10)) if t]
        c["multi-placement"] += 1
        if "pool" in cls and int(f[5]) < b: c["pool: primary is not the best AS"] += 1
        v = a.get(f[0])
        for k in cls:
            c[k] += 1; c[f"{k}: records beyond the first"] += int(f[3]) - 1
            if v is None: c[f"{k}: not in the assign table"] += 1; continue
            st = {3: "assigned", 2: "tied", 1: "ambiguous"}[v[0]]
            if st != "assigned": c[f"{k}: {st}"] += 1; continue
            if v[1]: c[f"{k}: assigned, sole candidate (no PSV)"] += 1; continue
            ch, lo, hi = span.get(v[3], ("", 0, 0)); kept = f[4] == ch and f[15] != "NA" and lo - 1000 <= int(f[15]) <= hi
            c[f"{k}: assigned on PSVs, " + ("primary already there" if kept else "MOVED off the primary")] += 1
            c[f"{k}: assigned on PSVs, records beyond the first"] += int(f[3]) - 1
            c[f"{k}: assigned on PSVs, n_candidates >= 3"] += v[2] >= 3
    for k in sorted(c): print(k, c[k], sep="\t")

def premise(assignments, *dumps):
    T = {}
    for d in dumps:
        for l in open(d):
            q, fam, st, v, al, _, t = l.rstrip("\n").split("\t")
            if t != "OUTSIDE" and al != "NA": T[(q, fam)] = (st, v, al)
    c = collections.defaultdict(collections.Counter)
    for r in csv.DictReader(open(assignments), delimiter="\t"):
        k = (r["read_name"], r["family_id"])
        if k not in T or r["as_margin"] in ("", "NA") or r.get("ms_margin", "NA") in ("", "NA"): continue
        a, m = int(r["as_margin"]), int(r["ms_margin"]); st, v, al = T.pop(k)
        s = c[("AS tie" if a == 0 else "AS 1-10" if a <= 10 else "AS > 10") + ", " + ("ms tie" if m == 0 else "ms 1-10")]
        s["n"] += 1; s["prim"] += al == "al_right"; s["chance"] += 1 / max(int(r["n_candidates"]), 1)
        s["sole"] += r["sole_candidate"] == "1"
        if st == "assigned" and r["sole_candidate"] != "1": s["o2_" + v] += 1; s["o2_right_prim_wrong"] += v == "right" and al != "al_right"; s["o2_wrong_prim_right"] += v == "wrong" and al == "al_right"
    print("stratum\tdecidable\tprimary right\tmean 1/candidates\tsole rows\tO2 non-sole right/wrong\tO2 right & primary wrong\tO2 wrong & primary right")
    for k in sorted(c):
        s = c[k]; print(k, s["n"], f"{s['prim'] / s['n']:.3f}", f"{s['chance'] / s['n']:.3f}", s["sole"], f"{s['o2_right']}/{s['o2_wrong']}", s["o2_right_prim_wrong"], s["o2_wrong_prim_right"], sep="\t")

if __name__ == "__main__": {"funnel": funnel, "premise": premise}[sys.argv[1]](*sys.argv[2:])
