#!/usr/bin/env python3
"""Pool the family-level scores of every arm at level M (MCL clusters, registered) and level C (components of the pre-MCL graph), Amendment 2 of docs/PREREG_locus_units_2026-10-06.md.

    pool_levels.py [IDEALDIR WORKDIR]
Per arm and level: reachable (R) and entangled (E) gene-runs found over the four runs (IDEAL-FOUND), E2 / E3 / E4 counts, per-run R found against the bar (NPIP 13, TBC1D3 12), the registered verdict, cluster precision of K*,
Compara F on the windows, and the E4 decomposition at R and E (not in K*_C: no edge path; in K*_C but not in K*_M: partition split). Prints the CANDIDATE_C test of the amendment against D.
"""
import csv
import json
import math
import os
import sys

I = sys.argv[1] if len(sys.argv) > 1 else "/mnt/linuxdisk/tmp/ideal_expression_2026-10-06"
W = sys.argv[2] if len(sys.argv) > 2 else "/mnt/linuxdisk/tmp/entangled_2026-10-06"
ARMS = ["D", "P", "PC", "q", "sd", "sq", "f1", "f1q"]
FAMS = [("NPIP", 14), ("TBC1D3", 13)]


def rd(p):
    return {r["cid"]: r for r in csv.DictReader(open(p), delimiter="\t")}


def paths(arm, fam, rep, level):
    O = f"{W}/{fam}/rep{rep}"
    if level == "C":
        return f"{O}/scoreC_{arm}", f"{O}/famscoreC_{arm}.json"
    if arm == "D":
        return f"{I}/{fam}/rep{rep}/score", f"{I}/{fam}/rep{rep}/famscore.json"
    return f"{O}/score_{arm}", f"{O}/famscore_{arm}.json"


def load(arm, level):
    out = {}
    for fam, _ in FAMS:
        for rep in (1, 2):
            sp, fj = paths(arm, fam, rep, level)
            try:
                rows = rd(f"{sp}.{fam}.copies.tsv")
                summ = json.load(open(sp + ".summary.json"))[fam]
            except Exception:
                return None
            try:
                F = json.load(open(fj))["truths"]["compara"]["f"]
            except Exception:
                F = None
            out[(fam, rep)] = (rows, summ, F)
    return out


def main():
    res = {}
    for arm in ARMS:
        for level in ("M", "C"):
            d = load(arm, level)
            if d:
                res[(arm, level)] = d
    base = {}
    print(f"{'arm':4} {'lvl':3} {'R':>3} {'E':>3} {'R+E':>4} | R E2/E3/E4 | E E2/E3/E4 | NPIP R r1/r2 | TBC1D3 R r1/r2 | verdict N/T | CP min | F NPIP | F TBC1D3")
    for arm in ARMS:
        for level in ("M", "C"):
            d = res.get((arm, level))
            if not d:
                continue
            tR = {k: 0 for k in ("E2", "E3", "E4", "F")}; tE = dict(tR)
            perrun = {}
            cps, Fs = [], {}
            for (fam, rep), (rows, summ, F) in d.items():
                for row in rows.values():
                    if row.get("heldout") == "1":
                        continue
                    grp = "R" if row["in_R"] == "1" else ("E" if row["stratum_in"].startswith("E") else None)
                    if grp is None:
                        continue
                    t = tR if grp == "R" else tE
                    t["E2"] += int(row["E2"] or 0); t["E3"] += int(row["E3"] or 0); t["E4"] += int(row["E4"] or 0); t["F"] += int(row["IDEAL_FOUND"] or 0)
                perrun[(fam, rep)] = summ["R_counts"]["IDEAL_FOUND"]
                cps.append(summ["CP"]); Fs[(fam, rep)] = F
            ver = {}
            for fam, rn in FAMS:
                need = math.ceil(0.9 * rn)
                ver[fam] = "YES" if all(perrun[(fam, r)] >= need and d[(fam, r)][1]["CP"] >= 0.5 for r in (1, 2)) else "NO"
            fm = lambda fam: "/".join("-" if Fs[(fam, r)] is None else f"{Fs[(fam, r)]:.3f}" for r in (1, 2))
            print(f"{arm:4} {level:3} {tR['F']:3} {tE['F']:3} {tR['F'] + tE['F']:4} | {tR['E2']:2}/{tR['E3']:2}/{tR['E4']:2} | {tE['E2']:2}/{tE['E3']:2}/{tE['E4']:2} | {perrun[('NPIP', 1)]:2}/{perrun[('NPIP', 2)]:2} | {perrun[('TBC1D3', 1)]:2}/{perrun[('TBC1D3', 2)]:2} | {ver['NPIP']}/{ver['TBC1D3']} | {min(cps):.2f} | {fm('NPIP')} | {fm('TBC1D3')}")
            if arm == "D":
                base[level] = (tR["F"], tE["F"], Fs, min(cps))
    # E4 decomposition: failures not in K*_C versus in K*_C but not in K*_M
    print("\nE4 failures at level M split by level C (R + E gene-runs):  arm: not in K*_M | of which not in K*_C (no edge path) | in K*_C (partition split)")
    for arm in ARMS:
        m, c = res.get((arm, "M")), res.get((arm, "C"))
        if not m or not c:
            continue
        fm_, nc, sp = 0, 0, 0
        for key in m:
            for cid, row in m[key][0].items():
                if row.get("heldout") == "1" or row["stratum_in"] == "C":
                    continue
                if row["IDEAL_FOUND"] == "1" or row["E4"] == "1" or not row["holder"]:
                    continue
                fm_ += 1
                rc = c[key][0].get(cid)
                if rc and rc["E4"] == "1":
                    sp += 1
                else:
                    nc += 1
        print(f"  {arm:4}: {fm_:3} | {nc:3} | {sp:3}")
    # CANDIDATE_C
    print("\nCANDIDATE_C (pooled IDEAL-FOUND_C above D's; R >= D's; E >= D's - 1; F_C >= D's - .005 in both families; CP_C >= .5 in all runs):")
    for arm in ARMS:
        if arm == "D" or (arm, "C") not in res:
            continue
        d = res[(arm, "C")]
        tR = sum(int(r["IDEAL_FOUND"]) for (k, (rows, s, F)) in d.items() for r in rows.values() if r["in_R"] == "1" and r.get("heldout") != "1")
        tE = sum(int(r["IDEAL_FOUND"]) for (k, (rows, s, F)) in d.items() for r in rows.values() if r["in_R"] != "1" and r["stratum_in"].startswith("E") and r.get("heldout") != "1")
        okF = all(d[k][2] is not None and base["C"][2][k] is not None and d[k][2] >= base["C"][2][k] - 0.005 for k in d)
        okCP = all(d[k][1]["CP"] >= 0.5 for k in d)
        cand = (tR + tE > base["C"][0] + base["C"][1]) and tR >= base["C"][0] and tE >= base["C"][1] - 1 and okF and okCP
        print(f"  {arm:4}: R {tR} (D {base['C'][0]}) E {tE} (D {base['C'][1]}) R+E {tR + tE} (D {base['C'][0] + base['C'][1]}) F-clause {okF} CP-clause {okCP} -> {'CANDIDATE_C' if cand else 'not a candidate'}")


if __name__ == "__main__":
    main()
