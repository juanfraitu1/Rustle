#!/usr/bin/env python3
"""PREREG r1797 / r1798 (docs/PREREG_famwise_sim_2026-11-18.md): the statistics the mimic bar is read from, of one truth run.

  mimic_stats.py LABEL READS.tsv S0 S1 S2 [READ_FAMILY.tsv [TIE_SPLIT.tsv]]

READS.tsv = the per-read table of bench/o2_truth/neutral_truth.py; S0 S1 S2 = reads per stratum in the whole population
(the table samples S0 and S1 and must hold every S2 read); READ_FAMILY.tsv = read, family for every read of the
population; TIE_SPLIT.tsv = the tie split. Runs unchanged on a simulated truth run (same file formats).
LABEL.stats.tsv  stat, value, k, n, lo, hi: Wilson 95% for a proportion; for the fam_ rows the 2.5% and 97.5% points of
                 2000 resamplings of the families with >= 50 reads (seed 1).
LABEL.fam.tsv    per family: reads, near ties (S2), their statuses and split classes; most near ties first."""
import collections, csv, random, sys

LOCATED = ("tie", "outside", "competitor unlocatable", "decidable")
NO_HIT = "no near-perfect hit"

def wilson(k, n, z=1.96):
    if not n:
        return 0.0, 1.0
    p, d = k / n, 1 + z * z / n
    c, h = (p + z * z / (2 * n)) / d, z * (p * (1 - p) / n + z * z / (4 * n * n)) ** 0.5 / d
    return max(0.0, c - h), min(1.0, c + h)

def picks(r):
    return dict(a.split("=") for a in r["picks"].split(";"))

def wrong(r):
    """the primary is not on a locus of origin"""
    return picks(r)["primary"] not in r["origin_loci"].split(",")

def famstats(F):
    """F = Counters (n reads, S2 near ties) of the families counted"""
    n = sum(c["n"] for c in F); s2 = sum(c["S2"] for c in F); by = sorted(F, key=lambda c: -c["S2"]); acc = k = 0
    while acc * 2 < s2:
        acc += by[k]["S2"]; k += 1
    sh = sorted(c["S2"] / c["n"] for c in F); q = lambda p: sh[min(len(sh) - 1, int(p * len(sh)))]
    at = lambda x: sum(c["S2"] >= x * c["n"] for c in F) / len(F)
    pooled = lambda X: sum(c["S2"] for c in X) / sum(c["n"] for c in X)
    return dict(pooled=s2 / n, half_in=k, zero=sum(not c["S2"] for c in F) / len(F), ge5=at(.05), ge20=at(.2), ge50=at(.5), p90=q(.9), p95=q(.95),
                pooled_minus_top1=pooled(by[1:]), pooled_minus_top4=pooled(by[4:]),
                top1_of_S2=by[0]["S2"] / s2, top4_of_S2=sum(c["S2"] for c in by[:4]) / s2, mean_family_share=sum(sh) / len(sh))

def main():
    lab, tsv, n0, n1, n2 = sys.argv[1:6]
    famf = sys.argv[6] if len(sys.argv) > 6 else None; splitf = sys.argv[7] if len(sys.argv) > 7 else None
    tot = {"S0": int(n0), "S1": int(n1), "S2": int(n2)}; N = sum(tot.values())
    with open(tsv) as fh:
        R = list(csv.DictReader(fh, delimiter="\t"))
    out = []
    def prop(name, k, n):
        lo, hi = wilson(k, n); out.append((name, f"{k / n:.4f}" if n else "NA", k, n, f"{lo:.4f}", f"{hi:.4f}"))
    def val(name, v):
        out.append((name, v, "", "", "", ""))
    def statuses(name, X, sts):
        for st in sts:
            prop(f"{name}_{st.replace(' ', '_')}", sum(r["status"] == st for r in X), len(X))
    for s in ("S0", "S1", "S2"):
        prop(f"stratum_{s}", tot[s], N)
    S0, S1, S2 = ([r for r in R if r["stratum"] == s] for s in ("S0", "S1", "S2"))
    assert len(S2) == tot["S2"], "S2 capped: set NT_CAP"
    statuses("S2_status", S2, ("tie", NO_HIT, "outside", "competitor unlocatable", "decidable"))
    statuses("S2_located_status", [r for r in S2 if r["status"] != NO_HIT], LOCATED)    # given a locatable read: takes the error tail out
    for s, X in (("S0", S0), ("S1", S1)):
        prop(f"{s}_no_near_perfect_hit", sum(r["status"] == NO_HIT for r in X), len(X))
        statuses(f"{s}_located_status", [r for r in X if r["status"] != NO_HIT], LOCATED)
    for s, X in (("S2", S2), ("S1", S1)):
        d = [r for r in X if r["status"] == "decidable"]
        prop(f"{s}_decidable_primary_wrong", sum(map(wrong, d)), len(d))
    h = collections.Counter(int(r["lead"]) for r in S2)
    for k in (0, 5, 10):
        prop(f"S2_lead_eq_{k}", h[k], len(S2))
    prop("S2_lead_other", len(S2) - h[0] - h[5] - h[10], len(S2))
    l1 = sorted(int(r["lead"]) for r in S1)
    for p in (5, 10, 25, 50, 75, 90, 95):
        val(f"S1_lead_p{p}", l1[min(len(l1) - 1, int(p / 100 * len(l1)))])
    for a, b in ((11, 20), (21, 50), (51, 100)):
        prop(f"S1_lead_{a}_{b}", sum(a <= x <= b for x in l1), len(l1))
    prop("S1_best_AS_locus_not_primary", sum(picks(r)["primary"] != "0" for r in S1), len(S1))
    nl = sorted(len(r["loci"].split(";")) for r in S1 + S2); val("loci_per_read_median_S1S2_sample", nl[len(nl) // 2])
    cls = {}
    if splitf:
        with open(splitf) as fh:
            cls = {r["read"]: r["class"] for r in csv.DictReader(fh, delimiter="\t")}
        c = collections.Counter(cls.get(r["read"], "no row") for r in S2 if r["status"] == "tie")
        prop("S2_tie_several_copies", c["SEVERAL COPIES"], c["SEVERAL COPIES"] + c["ONE PLACE"]); val("S2_tie_without_split_row", c["no row"] + c["UNSPLIT"])
    if famf:
        with open(famf) as fh:
            fam = dict(l.rstrip("\n").split("\t")[:2] for l in fh)
        assert len(fam) == N, (len(fam), N)
        C = collections.defaultdict(collections.Counter)
        for f in fam.values():
            C[f]["n"] += 1
        for r in S2:
            c = C[fam[r["read"]]]; c["S2"] += 1; c[r["status"]] += 1; c["wrong"] += r["status"] == "decidable" and wrong(r)
            if splitf and r["status"] == "tie":
                c[cls.get(r["read"], "no row")] += 1
        cols = ("tie", "decidable", NO_HIT, "outside", "competitor unlocatable", "wrong", "SEVERAL COPIES", "ONE PLACE")
        with open(lab + ".fam.tsv", "w") as w:
            w.write("family\treads\tnear_ties\tshare\ttie\tdecidable\tno near-perfect hit\toutside\tcompetitor unlocatable\tprimary_wrong\tSEVERAL COPIES\tONE PLACE\n")
            for f in sorted(C, key=lambda f: -C[f]["S2"]):
                c = C[f]; w.write(f"{f}\t{c['n']}\t{c['S2']}\t{c['S2'] / c['n']:.4f}\t" + "\t".join(str(c[k]) for k in cols) + "\n")
        F = [c for c in C.values() if c["n"] >= 50]; val("families_with_reads", len(C)); val("families_counted_ge50", len(F))
        base = famstats(F); rng = random.Random(1); B = [famstats([rng.choice(F) for _ in F]) for _ in range(2000)]
        for k, v in base.items():
            b = sorted(x[k] for x in B); fmt = (lambda x: f"{x:.4f}") if isinstance(v, float) else (lambda x: x)
            out.append((f"fam_{k}", fmt(v), "", len(F), fmt(b[50]), fmt(b[1949])))
        val("fam_n_zero", sum(not c["S2"] for c in F))
        for x in (5, 10, 20, 50):
            val(f"fam_n_ge{x}", sum(c["S2"] >= x / 100 * c["n"] for c in F))
        rs = sorted(c["n"] for c in C.values())
        for p in (25, 50, 75, 90):
            val(f"reads_per_family_p{p}", rs[min(len(rs) - 1, int(p / 100 * len(rs)))])
        val("reads_per_family_max", rs[-1]); prop("top5_families_share_of_reads", sum(rs[-5:]), N)
    with open(lab + ".stats.tsv", "w") as w:
        w.write("stat\tvalue\tk\tn\tlo\thi\n")
        for r in out:
            w.write("\t".join(map(str, r)) + "\n")
    print(open(lab + ".stats.tsv").read(), end="")

if __name__ == "__main__":
    main()
