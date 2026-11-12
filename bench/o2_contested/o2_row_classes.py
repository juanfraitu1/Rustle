#!/usr/bin/env python3
"""r1656 (docs/PREREG_o2_evidence_free_assignments_2026-11-12.md): decidable assigned rows of copy_assign split by blind
row classes (n_decisive == 0, tie_outside_catalog, sole_candidate) against strict_truth verdicts.
Usage: o2_row_classes.py DUMP[,DUMP2] ASSIGN[,ASSIGN2]"""
import collections, csv, sys

def main(dumps, assigns):
    truth = {}
    for d in dumps.split(","):
        for l in open(d):
            f = l.rstrip("\n").split("\t")
            if f[2] == "assigned" and f[6] != "OUTSIDE": truth[(f[0], f[1])] = f[3]
    c = collections.Counter(); fam = collections.defaultdict(collections.Counter); tot = collections.Counter()
    for a in assigns.split(","):
        for r in csv.DictReader(open(a), delimiter="\t"):
            if r["status"] != "assigned": continue
            nd0 = r["n_decisive"] == "0"; toc = r.get("tie_outside_catalog") == "1"; sole = r["sole_candidate"] == "1"
            tot["assigned rows"] += 1; tot["assigned rows, n_decisive 0"] += nd0
            t = truth.get((r["read_name"], r["family_id"]))
            if not t: continue
            for name, on in (("n_decisive == 0", nd0), ("n_decisive >= 1", not nd0), ("tie_outside_catalog", toc), ("sole_candidate", sole), ("n_decisive >= 1 and not tie_outside", not nd0 and not toc), ("all", True)):
                if on: c[(name, t)] += 1
            if nd0: fam[r["family_id"]][t] += 1
    for name in ("all", "n_decisive == 0", "n_decisive >= 1", "tie_outside_catalog", "sole_candidate", "n_decisive >= 1 and not tie_outside"):
        r, w = c[(name, "right")], c[(name, "wrong")]; print(f"{name}\tright {r}\twrong {w}\twrong rate {w / max(r + w, 1):.3f}")
    print(dict(tot), "families of the n_decisive 0 rows:", [(f, dict(v)) for f, v in sorted(fam.items(), key=lambda x: -sum(x[1].values()))[:6]])

if __name__ == "__main__": main(*sys.argv[1:3])
