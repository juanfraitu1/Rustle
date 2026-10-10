#!/usr/bin/env python3
"""Score docs/PREREG_asj_o2_2026-10-10.md on the per-chromosome JSONL of run_chrom.py.

    score.py dev|heldout|all|X        dev = odd chrmap numbers (1, 3, ..., 23), heldout = even, all = every autosome, X = the chrX control

Pure functions (calls, t1, t2) are tested in test_score.py."""
import collections
import json
import os
import statistics
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import asj_o2 as A  # noqa: E402

OUT = os.environ.get("ASJ_OUT", "/mnt/linuxdisk/tmp/asj_o2")
CHRMAP = "/mnt/linuxdisk/tmp/rna_allele/chrmap.tsv"


def calls(rows, key):
    """BH over every test of the arm in `rows`; -> the tests with q < 0.05 and |dPSI| >= 0.3, each with its gene"""
    tests = [(r, t) for r in rows for t in r[key]]
    qs = A.bh([t["p"] for _, t in tests])
    return [dict(t, q=q, gene=r["gene"], name=r.get("name"), chrom=r["chrom"], anchor=r["anchor"]) for (r, t), q in zip(tests, qs)
            if q < A.Q_MAX and abs(t["dpsi"]) >= A.DPSI_MIN]


def t1(rows):
    """-> (discordant, reads assigned by both halves)"""
    return sum(r["T1"]["discordant"] for r in rows), sum(r["T1"]["both"] for r in rows)


def t2(rows, min_reads=10):
    """-> {class: (anchors with >= min_reads concordance reads, median concordance or None)}"""
    out = {}
    for c in ("ALLELIC", "NON-ALLELIC", "UNRESOLVED"):
        v = [r["anchor"]["conc"] for r in rows if r["anchor"] and r["anchor"]["class"] == c and r["anchor"]["n_conc"] >= min_reads]
        out[c] = (len(v), statistics.median(v) if v else None)
    return out


def chroms(which):
    rows = [ln.rstrip("\n").split("\t") for ln in list(open(CHRMAP))[1:]]
    if which == "X":
        return [r[0] for r in rows if r[1] == "X"]
    auto = [(r[0], int(r[1])) for r in rows if r[1].isdigit()]
    return [c for c, n in auto if which == "all" or (n % 2 == 1) == (which == "dev")]


def main(which):
    rows = []
    for c in chroms(which):
        p = f"{OUT}/{c}.jsonl"
        assert os.path.exists(p), f"missing {p}"
        rows += [json.loads(ln) for ln in open(p)]
    res = dict(half=which, chromosomes=len(chroms(which)), genes=len(rows), genes_with_allelic_alignment=sum(r["allelic_alignment"] for r in rows),
               median_het_sites=statistics.median(r["n_het"] for r in rows))
    a = collections.Counter()
    for r in rows:
        a.update(r["armA"])
    res["armA_reads"] = dict(a)
    res["armA_abstention"] = round(a["abstain"] / sum(a.values()), 4)
    d, b = t1(rows)
    res["T1"] = dict(discordant=d, both=b, rate=round(d / b, 5) if b else None, PASS=(d / b <= 0.005) if b else None)
    t = t2(rows)
    res["T2"] = {k: dict(anchors=n, median_concordance=m) for k, (n, m) in t.items()}
    if t["ALLELIC"][0] >= 10 and t["NON-ALLELIC"][0] >= 10:
        res["T2"]["PASS"] = t["ALLELIC"][1] >= 0.95 and t["NON-ALLELIC"][1] <= 0.80
    else:
        res["T2"]["PASS"] = "UNDERPOWERED"
    ca, cr = calls(rows, "testsA"), calls(rows, "testsR")
    res["armA_tests"] = sum(len(r["testsA"]) for r in rows)
    res["armA_calls"] = dict(calls=len(ca), genes=len({x["gene"] for x in ca}),
                             median_abs_dpsi=round(statistics.median(abs(x["dpsi"]) for x in ca), 3) if ca else None)
    res["armR_tests"] = sum(len(r["testsR"]) for r in rows)
    by = collections.Counter(x["anchor"]["class"] for x in cr)
    res["armR_calls"] = dict(calls=len(cr), genes=len({x["gene"] for x in cr}), by_class=dict(by),
                             genes_by_class={k: len({x["gene"] for x in cr if x["anchor"]["class"] == k}) for k in by})
    resolved = by["ALLELIC"] + by["NON-ALLELIC"]
    res["armR_nonallelic_share_of_resolved_calls"] = round(by["NON-ALLELIC"] / resolved, 4) if resolved else None
    sa = {(x["gene"], tuple(x["intron"])) for x in ca}
    al = [x for x in cr if x["anchor"]["class"] == "ALLELIC"]
    res["armR_allelic_calls_also_armA"] = f"{sum((x['gene'], tuple(x['intron'])) in sa for x in al)}/{len(al)}"
    anchors = [r["anchor"] for r in rows if r["anchor"]]
    res["anchors_by_class"] = dict(collections.Counter(x["class"] for x in anchors))
    for k in ("ALLELIC", "NON-ALLELIC", "UNRESOLVED"):
        xs = [x for x in anchors if x["class"] == k]
        xc = [x["anchor"] for x in cr if x["anchor"]["class"] == k]
        res[f"{k}_anchor_transition_share"] = round(sum(x["transition"] for x in xs) / len(xs), 4) if xs else None
        res[f"{k}_call_transition_share"] = round(sum(x["transition"] for x in xc) / len(xc), 4) if xc else None
        alt = [x["mapq0_alt"] for x in xs if x["mapq0_alt"] is not None]
        ref = [x["mapq0_ref"] for x in xs if x["mapq0_ref"] is not None]
        res[f"{k}_median_mapq0_alt_ref"] = (round(statistics.median(alt), 4) if alt else None, round(statistics.median(ref), 4) if ref else None)
    res["anchors_per_gene"] = round(len(anchors) / len(rows), 4)
    res["armR_calls_per_gene"] = round(len(cr) / len(rows), 4)
    res["top_armA_calls"] = sorted(((round(x["q"], 6), x["name"], x["intron"], round(x["dpsi"], 3)) for x in ca))[:15]
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main(sys.argv[1])
