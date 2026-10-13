#!/usr/bin/env python3
"""Arm A of docs/PREREG_yag_o2_2026-10-12.md: O2 on simulated chrY reads with truth.

    arm_a.py RUN_PREFIX OUT.json
RUN_PREFIX = a `copy_assign --dump-star` run on reads named `FAM|copy|i` (the truth copy).
distinguishable = the truth copy t is a candidate and every other candidate c has >= 2 columns where t's and c's bases
differ (both defined) — O2's own m >= 2 step.
"""
import collections, json, math, sys


def parse_cols(field):
    out = []
    if not field:
        return out
    for tok in field.split(","):
        pos, obs, al = tok.split(":", 2)
        out.append((int(pos), obs, al))
    return out


def upper95(wrong, n):
    """one-sided 95% upper bound on an error rate (Clopper-Pearson; 'rule of three' at 0 wrong)."""
    if n == 0:
        return None
    if wrong == 0:
        return round(1 - 0.05 ** (1 / n), 4)
    lo, hi = wrong / n, 1.0
    for _ in range(100):  # binomial CDF bisection
        mid = (lo + hi) / 2
        cdf = sum(math.comb(n, k) * mid ** k * (1 - mid) ** (n - k) for k in range(wrong + 1))
        lo, hi = (mid, hi) if cdf > 0.05 else (lo, mid)
    return round(hi, 4)


def main(prefix, out):
    rows = []
    with open(prefix + ".star_reads.tsv") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            rows.append(dict(zip(head, line.rstrip("\n").split("\t"))))
    c = collections.Counter()
    min_dec_hist = collections.Counter()
    import os
    only = set(json.load(open(os.environ["YAG_FAMILIES"]))["families"]) if os.environ.get("YAG_FAMILIES") else None
    for r in rows:
        if only is not None and r["family_id"] not in only:  # YAG families only (user, 2026-10-12)
            continue
        fam, t = r["read_name"].split("|")[0], int(r["read_name"].split("|")[1])
        if fam != r["family_id"]:
            c["truth_family_not_this_family"] += 1
            continue
        cand = [int(x) for x in r["candidates"].split(",")] if r["candidates"] else []
        status = r["status"]
        assigned = status == "assigned"
        correct = assigned and int(r["assigned_copy"]) == t
        if t not in cand:
            c[f"truth_not_candidate|{status}"] += 1
            if assigned:
                c["assigned_wrong"] += 1
            continue
        ti = cand.index(t)
        cols = parse_cols(r.get("columns", ""))
        decs = []
        for ci in range(len(cand)):
            if ci == ti:
                continue
            d = 0
            for _, obs, al in cols:
                a, b = al[ti], al[ci]
                if a not in ".-" and b not in ".-" and a != b and obs not in ".-":
                    d += 1
            decs.append(d)
        md = min(decs) if decs else 0
        min_dec_hist[min(md, 5)] += 1
        dist = bool(decs) and md >= 2
        key = "dist" if dist else "notdist"
        c[f"{key}|n"] += 1
        c[f"{key}|assigned"] += assigned
        c[f"{key}|correct"] += correct
        c[f"{key}|{status}"] += 1
        if assigned and not correct:
            c["assigned_wrong"] += 1
    assigned = c["dist|assigned"] + c["notdist|assigned"] + sum(v for k, v in c.items() if k.startswith("truth_not_candidate|assigned"))
    wrong = c["assigned_wrong"]
    res = {
        "molecules": len(rows),
        "counts": dict(c),
        "min_decisive_columns_vs_closest_rival (5 = >=5)": dict(sorted(min_dec_hist.items())),
        "A_P1_precision": {"assigned": assigned, "wrong": wrong, "upper95_error": upper95(wrong, assigned),
                           "verdict": "UNDERPOWERED" if assigned < 30 else ("PASS" if wrong == 0 else "FAIL")},
        "A_P2_recall_distinguishable": {"n": c["dist|n"], "correct": c["dist|correct"],
                                        "rate": round(c["dist|correct"] / c["dist|n"], 4) if c["dist|n"] else None},
        "A_P3_calibration_not_distinguishable": {"n": c["notdist|n"], "assigned": c["notdist|assigned"],
                                                 "rate": round(c["notdist|assigned"] / c["notdist|n"], 4) if c["notdist|n"] else None},
    }
    r2 = res["A_P2_recall_distinguishable"]
    r2["verdict"] = None if not r2["n"] else ("PASS" if r2["rate"] >= 0.80 else "FAIL")
    r3 = res["A_P3_calibration_not_distinguishable"]
    r3["verdict"] = None if not r3["n"] else ("PASS" if r3["rate"] <= 0.05 else "FAIL")
    json.dump(res, open(out, "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
