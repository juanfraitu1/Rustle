#!/usr/bin/env python3
"""r1657 (docs/PREREG_o2_cross_family_candidates_2026-11-12.md): arm U (--union-certificate) against arm F (BASE with
sole-candidate rows tied) on strict_truth dumps. Usage: o2_union_vs_flag.py NAME=DIR [NAME=DIR ...]
(DIR holds ms10.assignments.tsv, ms10.dump.tsv, u.dump.tsv)"""
import collections, csv, sys

def load(d):
    sole = {(r["read_name"], r["family_id"]) for r in csv.DictReader(open(d + "/ms10.assignments.tsv"), delimiter="\t")
            if r["status"] == "assigned" and r["sole_candidate"] == "1"}
    def dump(p):  # decidable assigned rows -> verdict
        out = {}
        for l in open(p):
            f = l.rstrip("\n").split("\t")
            if f[2] == "assigned" and f[6] != "OUTSIDE": out[(f[0], f[1])] = f[3]
        return out
    return sole, dump(d + "/ms10.dump.tsv"), dump(d + "/u.dump.tsv")

def rw(rows, skip=None):
    c = collections.Counter(v for (q, f), v in rows.items() if f != skip); return c["right"], c["wrong"]

def rate(r, w): return w / max(r + w, 1)

def main(args):
    pooled = {"B": {}, "F": {}, "U": {}}; sole_all = set(); ok12 = True
    for a in args:
        name, d = a.split("="); sole, base, u = load(d)
        f = {k: v for k, v in base.items() if k not in sole}
        for k, rows in (("B", base), ("F", f), ("U", u)): pooled[k].update({(name + q, fam): v for (q, fam), v in rows.items()})
        sole_all |= {(name + q, fam) for (q, fam) in sole if (q, fam) in base}
        (br, bw), (fr, fw), (ur, uw) = rw(base), rw(f), rw(u)
        b1, b2 = ur >= fr, rate(ur, uw) <= rate(fr, fw) + 0.005; ok12 &= b1 and b2
        print(f"{name}\tBASE {br}/{bw} ({rate(br, bw):.3f})\tF {fr}/{fw} ({rate(fr, fw):.3f})\tU {ur}/{uw} ({rate(ur, uw):.3f})\tbar1 {'PASS' if b1 else 'FAIL'}\tbar2 {'PASS' if b2 else 'FAIL'}")
        fate = collections.Counter()  # BASE sole decidable rows under U
        umol = collections.defaultdict(list)
        for (q, fam), v in u.items(): umol[q].append((fam, v))
        for (q, fam) in sole:
            if (q, fam) not in base: continue
            here = u.get((q, fam)); other = [v for g, v in umol[q] if g != fam]
            fate[(base[(q, fam)], "same row " + here if here else "other family " + other[0] if other else "not assigned")] += 1
        for k, n in sorted(fate.items()): print(f"  sole in BASE {k[0]} -> U {k[1]}\t{n}")
    (fr, fw), (ur, uw) = rw(pooled["F"]), rw(pooled["U"])
    print(f"bar0 decidable BASE sole rows pooled {len(sole_all)} ({'PASS' if len(sole_all) >= 20 else 'NO VERDICT'})")
    print(f"bar3 pooled U right - F right = {ur - fr} ({'PASS' if ur - fr >= 10 else 'FAIL'})")
    big = collections.Counter(f for _, f in sole_all).most_common(1)
    if big:
        (fr, fw), (ur, uw) = rw(pooled["F"], big[0][0]), rw(pooled["U"], big[0][0])
        print(f"bar4 minus {big[0][0]}: F {fr}/{fw} ({rate(fr, fw):.3f}) U {ur}/{uw} ({rate(ur, uw):.3f}) ({'PASS' if ur >= fr and rate(ur, uw) <= rate(fr, fw) + 0.005 else 'FAIL'})")
    print("bars 1-2 on every set:", "PASS" if ok12 else "FAIL")

if __name__ == "__main__":
    assert rw({("a", "X"): "right", ("b", "Y"): "wrong"}, "Y") == (1, 0) and rate(3, 1) == 0.25
    main(sys.argv[1:])
