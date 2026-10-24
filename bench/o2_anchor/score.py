#!/usr/bin/env python3
"""Score O2 on anchor-trimmed reads (docs/archive/2026-10/merged/PREREG_o2_anchor_2026-10-13.md: T0, A1, A2 + reported numbers).

    score.py TRUTH.tsv CATALOG_COPIES.tsv TRIMMED.bam RUN_PREFIX OUT.json

Locus-level truth: a catalog copy is right iff its span overlaps the truth copy's span (same chromosome). Any `assigned`
row of a trimmed read in a family other than its truth family is wrong.
"""
import collections, json, math, subprocess, sys


def tsv(path):
    with open(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(head, l.rstrip("\n").split("\t"))) for l in fh]


def overlaps(a, b):
    return a[0] == b[0] and a[1] < b[2] and b[1] < a[2]


def parse_cols(field):
    out = []
    for tok in field.split(",") if field else []:
        pos, obs, al = tok.split(":", 2)
        out.append((int(pos), obs, al))
    return out


def distinguishable(cand, cols, truth_set):
    """YAG Arm A's rule: some truth candidate t such that every rival (not in truth_set) differs from t at >= 2 columns
    where both have a base and the read has a base."""
    ts = [i for i, c in enumerate(cand) if c in truth_set]
    if not ts:
        return False
    rivals = [i for i, c in enumerate(cand) if c not in truth_set]
    for ti in ts:
        ok = True
        for ri in rivals:
            n = sum(1 for _, o, al in cols if o not in ".-" and al[ti] not in ".-" and al[ri] not in ".-" and al[ti] != al[ri])
            if n < 2:
                ok = False
                break
        if ok:
            return True
    return False


def upper95(wrong, n):
    if n == 0:
        return None
    if wrong == 0:
        return round(1 - 0.05 ** (1 / n), 4)
    lo, hi = wrong / n, 1.0
    for _ in range(100):
        mid = (lo + hi) / 2
        cdf = sum(math.comb(n, k) * mid ** k * (1 - mid) ** (n - k) for k in range(wrong + 1))
        lo, hi = (mid, hi) if cdf > 0.05 else (lo, mid)
    return round(hi, 4)


def amended_a1(c, n_as):
    """Amendment 1: d = the truth's own disagreement (1 - T0 rate); PASS iff wrong / assigned <= d + 0.01, >= 30 assigned."""
    d = 1 - c["t0_right"] / c["t0_n"] if c["t0_n"] else None
    if d is None:
        return {"d": None, "verdict": "NO_T0"}
    if d > 0.05:
        return {"d": round(d, 4), "verdict": "STOP (d > 0.05)"}
    v = "UNDERPOWERED" if n_as < 30 else ("PASS" if c["wrong"] / n_as <= d + 0.01 else "FAIL")
    return {"d": round(d, 4), "bar": round(d + 0.01, 4), "assigned": n_as, "wrong": c["wrong"],
            "wrong_own_sibling": c["wrong_own_sibling"], "wrong_own_other": c["wrong_own_other"],
            "wrong_foreign": c["wrong_foreign"], "verdict": v}


def primaries(bam):
    out = {}
    res = subprocess.run(["samtools", "view", "-F", "2308", bam], capture_output=True, text=True, check=True).stdout
    for l in res.splitlines():
        f = l.split("\t")
        import re
        end = int(f[3]) - 1 + sum(int(n) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", f[5]) if op in "MDN=X")
        out[f[0]] = (f[2], int(f[3]) - 1, end, int(f[4]))
    return out


def main(truth_tsv, catalog_tsv, bam, prefix, out):
    span = {(r["family_id"], r["copy_idx"]): (r["chrom"], int(r["start"]), int(r["end"])) for r in tsv(catalog_tsv)}
    truth = {r["read_name"]: r for r in tsv(truth_tsv)}
    prim = primaries(bam)
    join = {(r["family_id"], r["copy_index"]): r["catalog_copy_idx"] for r in tsv(prefix + ".family_join.tsv")}
    rows = collections.defaultdict(list)
    for r in tsv(prefix + ".assignments.tsv"):
        if r["read_name"] in truth:
            rows[r["read_name"]].append(r)
    star = {}
    for r in tsv(prefix + ".star_reads.tsv"):
        if r["read_name"] in truth and r["family_id"] == truth[r["read_name"]]["family_id"]:
            star[r["read_name"]] = r
    c = collections.Counter()
    per_fam = collections.defaultdict(collections.Counter)
    examples = {"wrong": [], "correct": []}
    for q, t in truth.items():
        fam, tspan = t["family_id"], span[(t["family_id"], t["truth_copy"])]
        p = prim.get(q)
        if p is None:
            c["unmapped"] += 1
            continue
        al_right = overlaps(p[:3], tspan)
        if p[3] >= 60:
            c["t0_n"] += 1
            c["t0_right"] += al_right
        own = [r for r in rows.get(q, []) if r["family_id"] == fam]
        if not own:
            c["no_own_row"] += 1
            continue
        if own[0]["contested"] != "1":
            c["own_not_contested"] += 1
            continue
        c["contested"] += 1
        per_fam[fam]["contested"] += 1
        c["aligner_right_contested"] += al_right
        # A1: every assigned row of this molecule
        assigned_right = False
        for r in rows[q]:
            if r["status"] != "assigned":
                continue
            c["assigned_rows"] += 1
            if r["family_id"] != fam:
                c["wrong_foreign"] += 1
                c["wrong"] += 1
                examples["wrong"].append({"read": q, "kind": "foreign", "family": r["family_id"], "copy": r["catalog_copy_idx"],
                                          "truth_family": fam, "truth_copy": t["truth_copy"]})
                continue
            ok = overlaps(span[(fam, r["catalog_copy_idx"])], tspan)
            if ok:
                c["right"] += 1
                assigned_right = True
                c["aligner_right_on_o2_assigned"] += al_right
                per_fam[fam]["right"] += 1
                if len(examples["correct"]) < 20:
                    examples["correct"].append({"read": q, "family": fam, "copy": r["catalog_copy_idx"],
                                                "aligner_primary_right": al_right, "mapq": p[3],
                                                "columns": star.get(q, {}).get("columns", "")[:400]})
            else:
                c["wrong"] += 1
                c["wrong_own"] += 1
                per_fam[fam]["wrong"] += 1
                sib = overlaps(span[(fam, r["catalog_copy_idx"])], span[(fam, t["sibling_copy"])])
                c["wrong_own_sibling" if sib else "wrong_own_other"] += 1
                examples["wrong"].append({"read": q, "kind": "own_sibling" if sib else "own_other", "family": fam, "copy": r["catalog_copy_idx"],
                                          "truth_copy": t["truth_copy"], "aligner_primary_right": al_right,
                                          "columns": star.get(q, {}).get("columns", "")[:400]})
        # A2: distinguishability on the own-family star row
        s = star.get(q)
        if s is None:
            c["no_star_row"] += 1
            continue
        cand = [int(x) for x in s["candidates"].split(",")] if s["candidates"] else []
        tset = {k for k in cand if overlaps(span[(fam, join[(fam, str(k))])], tspan)}
        if distinguishable(cand, parse_cols(s["columns"]), tset):
            c["dist"] += 1
            c["dist_right"] += assigned_right
        else:
            c["notdist"] += 1
            c["notdist_assigned"] += any(r["status"] == "assigned" for r in rows[q])
    n_as = c["assigned_rows"]
    res = {"counts": dict(c),
           "T0_truth_gate": {"n": c["t0_n"], "right": c["t0_right"],
                             "rate": round(c["t0_right"] / c["t0_n"], 4) if c["t0_n"] else None,
                             "pass": bool(c["t0_n"]) and c["t0_right"] / c["t0_n"] >= 0.99},
           "A1_precision": {"assigned": n_as, "wrong": c["wrong"], "upper95_error": upper95(c["wrong"], n_as),
                            "verdict": "UNDERPOWERED" if n_as < 30 else ("PASS" if c["wrong"] / n_as <= 0.01 else "FAIL")},
           "A1_amended": amended_a1(c, n_as),
           "A2_calibration": {"n": c["notdist"], "assigned": c["notdist_assigned"],
                              "rate": round(c["notdist_assigned"] / c["notdist"], 4) if c["notdist"] else None,
                              "verdict": None if not c["notdist"] else ("PASS" if c["notdist_assigned"] / c["notdist"] <= 0.05 else "FAIL")},
           "recall_distinguishable": {"n": c["dist"], "right": c["dist_right"],
                                      "rate": round(c["dist_right"] / c["dist"], 4) if c["dist"] else None},
           "aligner_vs_o2": {"contested": c["contested"],
                             "aligner_right_rate": round(c["aligner_right_contested"] / c["contested"], 4) if c["contested"] else None,
                             "on_o2_assigned_right": {"n": c["right"], "aligner_right": c["aligner_right_on_o2_assigned"]}},
           "per_family": {f: dict(v) for f, v in sorted(per_fam.items(), key=lambda x: -x[1]["contested"])},
           "examples": examples}
    json.dump(res, open(out, "w"), indent=1)
    print(json.dumps({k: v for k, v in res.items() if k not in ("per_family", "examples")}, indent=1))


if __name__ == "__main__":
    main(*sys.argv[1:6])
