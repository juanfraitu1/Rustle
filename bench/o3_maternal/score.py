#!/usr/bin/env python3
"""R1-R4 and the p12 line of docs/PREREG_o3_maternal_reference_2026-10-08.md section 6, per arm and per run (reference C.REF, env O3_REF).

    O3_REF=mat score.py isocon | inhouse   # reads WR/<dir>/{cands.tsv,cands.ref.paf,cands.other.paf}; writes WR/score/<arm>.json, prints verdicts

Loci passed to evaluate() are dicts: locus, kind (catalog|lrpap1|sex), family, chrom (accession on the TRUTH haplotype C.OTHER), start, end,
n (reads), ident (identity of the nearest reference paralog). Truth is used only here."""
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import common as C  # noqa: E402

MIN_SCORE = 0.999       # registered recovery hit: identity x coverage
ABSENT_KINDS = ("catalog", "lrpap1")      # loci that are absent from the reference; "sex" (control) and "lrpap1_desc" (descriptive, Amendment 2) are not


def recovered_loci(contigs, other_best, loci_, min_score=MIN_SCORE):
    """loci whose interval on the truth haplotype a contig's best hit there (score >= min_score) overlaps"""
    out = set()
    for c in contigs:
        h = other_best.get(c)
        if h and h[0] >= min_score:
            for L in loci_:
                if h[1] == L["chrom"] and h[2] < L["end"] and L["start"] < h[3]:
                    out.add(L["locus"])
    return out


def candidate_class(contigs, ref_best, other_best, recovers):
    """a_recovered (hits a truth locus) | ref (best hit >= .999 on the reference, the reference winning ties) | b_other (>= .999 on the truth
    haplotype, no truth locus) | c_unmatched"""
    if recovers:
        return "a_recovered"
    best = None
    for c in contigs:
        for hap, B in (("ref", ref_best), ("other", other_best)):
            h = B.get(c)
            if h and (best is None or h[0] > best[0]):
                best = (h[0], hap)
    if best is None or best[0] < MIN_SCORE:
        return "c_unmatched"
    return "ref" if best[1] == "ref" else "b_other"


def max_matching(edges, n_left):
    """Kuhn's maximum bipartite matching; edges = [(left, right)]"""
    adj = [[] for _ in range(n_left)]
    for a, b in edges:
        adj[a].append(b)
    match = {}

    def try_(u, seen):
        for v in adj[u]:
            if v in seen:
                continue
            seen.add(v)
            if v not in match or try_(match[v], seen):
                match[v] = u
                return True
        return False
    return sum(1 for u in range(n_left) if try_(u, set()))


def evaluate(cands, loci_, fam_net, ref_best, other_best, floor=2):
    """cands: dicts candidate, family, n_transcripts, contigs[list]; loci_: truth loci (see module doc); fam_net: every family that entered the chain.
    R1: each expressed (>= 3), LARGE (>= 20), beyond-delta catalog locus is recovered by a flagged candidate. R2: no expressed within-delta catalog
    locus is recovered by a flagged candidate. R3: families without an expressed truth locus having a flagged non-recovering candidate / those
    families (bar <= 20%). R4 (reported): sensitivity, precision, bipartite matching. kind 'sex' is never counted; kind 'lrpap1' (p12) is
    reported on the p12 line, not in R1/R2."""
    expressed = [L for L in loci_ if L["kind"] in ABSENT_KINDS and L["n"] >= 3]
    flagged = []
    for c in cands:
        if c["n_transcripts"] < floor:
            continue
        rec = recovered_loci(c["contigs"], other_best, loci_)
        flagged.append(dict(c, recovers=sorted(rec), cls=candidate_class(c["contigs"], ref_best, other_best, rec)))
    got = {l for c in flagged for l in c["recovers"]}
    beyond = 1 - C.DELTA
    t1 = [L["locus"] for L in expressed if L["kind"] == "catalog" and L["n"] >= 20 and L["ident"] < beyond]
    r1 = dict(targets=t1, recovered=[l for l in t1 if l in got], verdict=("PASS" if t1 and all(l in got for l in t1) else "NOT TESTABLE" if not t1 else "FAIL"))
    t2 = [L["locus"] for L in expressed if L["kind"] == "catalog" and L["ident"] >= beyond]
    flagged_within = [l for l in t2 if l in got]
    r2 = dict(targets=t2, flagged=flagged_within, verdict="FAIL" if flagged_within else "PASS")
    fam_with = {L["family"] for L in expressed}
    den = [f for f in fam_net if f not in fam_with]
    classes = {}
    for c in flagged:
        if c["cls"] != "a_recovered" and c["family"] in den:
            classes.setdefault(c["family"], set()).add(c["cls"])
    bad = sorted(classes)
    other_only = [f for f in bad if classes[f] == {"b_other"}]       # post hoc split: sequence found on the truth haplotype at >= .999 but in no truth locus
    frac = len(bad) / len(den) if den else None
    r3 = dict(denominator=len(den), false_families=bad, other_haplotype_only=other_only, unmatched_or_ref=[f for f in bad if f not in other_only], fraction=frac, verdict="NOT TESTABLE" if frac is None else ("PASS" if frac <= 0.20 else "FAIL"))
    idx = {L["locus"]: i for i, L in enumerate(expressed)}
    edges = [(i, idx[l]) for i, c in enumerate(flagged) for l in c["recovers"] if l in idx]
    m = max_matching(edges, len(flagged))
    r4 = dict(expressed=len(expressed), flagged=len(flagged), matched=m,
              sensitivity=(m / len(expressed)) if expressed else None, precision=(m / len(flagged)) if flagged else None)
    return dict(R1=r1, R2=r2, R3=r3, R4=r4, rows=flagged)


def type_hits(mat_best, pat_best, targets, min_score=MIN_SCORE):
    """targets: {name: (hap, acc, s, e)} -> {name: [sequences whose best hit on that haplotype scores >= min_score and overlaps the target]}"""
    out = {}
    for name, (hap, acc, s, e) in targets.items():
        B = mat_best if hap == "mat" else pat_best
        out[name] = sorted(q for q, h in B.items() if h[0] >= min_score and h[1] == acc and h[2] < e and s < h[3])
    return out


def load_loci():
    fate = json.load(open(f"{C.WR}/fate/fate.json"))["loci"]
    out = []
    for r in csv.DictReader(open(f"{C.WR}/truth/loci.tsv"), delimiter="\t"):
        f = fate.get(r["locus"], {})
        out.append(dict(locus=r["locus"], kind=r["kind"], family=r["family"], chrom=r["chrom"], start=int(r["start"]), end=int(r["end"]),
                        n=f.get("n", 0), ident=f.get("paralog_identity") if f.get("paralog_identity") is not None else 0.0))
    return out


def p12_targets():
    rows = {r["name"]: r for r in csv.DictReader(open(f"{C.WR}/truth/lrpap1_loci.tsv"), delimiter="\t")}
    t = {}
    for key, name in (("p12", "LOC134756368"), ("p14", "LOC115932954")):
        r = rows[name]
        if r["pat_acc"]:
            t[f"{key}@pat"] = ("pat", r["pat_acc"], int(r["pat_s"]), int(r["pat_e"]))
        if r["mat_acc"]:
            t[f"{key}@mat"] = ("mat", r["mat_acc"], int(r["mat_s"]), int(r["mat_e"]))
    return t


def main(arm):
    al = C.alias()
    d = f"{C.WR}/{'isoc' if arm == 'isocon' else 'inhouse'}"
    cands = [dict(candidate=r["candidate"], family=r["family"], n_transcripts=int(r["n_transcripts"]), contigs=r["contigs"].split(","))
             for r in csv.DictReader(open(f"{d}/cands.tsv"), delimiter="\t")]
    best = lambda p: C.best_hits(p, al) if os.path.exists(p) and os.path.getsize(p) else {}
    ref_best, other_best = best(f"{d}/cands.ref.paf"), best(f"{d}/cands.other.paf")
    fam_net = [p["fam"] for p in json.load(open(f"{d}/panel.json"))]
    res = evaluate(cands, load_loci(), fam_net, ref_best, other_best)
    res["candidates"] = len(cands)
    if arm == "isocon" and C.REF == "mat":      # the p12 line exists only with the mother as the reference (p12 is mother-absent)
        om, op = best(f"{d}/outputs.pri.paf"), best(f"{d}/outputs.other.paf")
        fam = lambda B: {k: v for k, v in B.items() if k.startswith("LRPAP1|")}
        res["p12"] = type_hits(fam(om), fam(op), p12_targets())
    os.makedirs(f"{C.WR}/score", exist_ok=True)
    json.dump(res, open(f"{C.WR}/score/{arm}.json", "w"), indent=1)
    for k in ("R1", "R2", "R3", "R4"):
        print(k, json.dumps(res[k]))
    print("p12 line (sequences whose best hit is >= .999 on each target):", json.dumps(res.get("p12", "n/a for this arm / reference")))
    print("flagged candidates by class:", {c: sum(1 for r in res["rows"] if r["cls"] == c) for c in ("a_recovered", "ref", "b_other", "c_unmatched")})


if __name__ == "__main__":
    main(sys.argv[1])
