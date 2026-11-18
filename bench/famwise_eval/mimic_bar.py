#!/usr/bin/env python3
"""PREREG r1797 / r1798 (docs/PREREG_famwise_sim_2026-11-18.md), row 1797: the five-item bar a simulated library has to
clear to be called a mimic of HG002. Real against the held-out quarter; the seen quarter is printed beside, no verdict.

  mimic_bar.py --seen PREFIX --held PREFIX --sim PREFIX [--ends REAL_HELD SIM] [--replay TOOL=REAL:SIM ...]

PREFIX.stats.tsv and PREFIX.fam.tsv are of mimic_stats.py; --ends takes two numbers or two read_stats.py --out tables
(both_ends_strict); --replay takes R of the real held-out run and of the simulated replay, once per tool.
F  per registered family the near-tie share (near_ties / reads): Wilson 95% intervals overlap, or larger / smaller
   < 1.5 (one share zero: the intervals alone); a family not in a table fails; >= 8 of 10.
S  the pooled near-tie share of all families but the first four of the list: the same test.
W  primary off its origin among decidable near-tied reads (S2_decidable_primary_wrong, k of n): intervals overlap.
E  both ends complete: |sim - real| <= 0.10.          A  R of the replay, Rustle and StringTie: |sim - real| <= 0.10.
Output: item, detail, seen, held, sim, verdict; then BAR pass, BAR fail (a failed item decides) or, with no failed item
and one whose inputs were not given, BAR incomplete."""
import argparse, csv, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mimic_stats import wilson

FAMILIES = ("CHM13_G0021052", "CHM13_G0044653", "CHM13_G0006774", "CHM13_G0020650", "CHM13_G0049830",
            "CHM13_G0058367", "CHM13_G0020705", "CHM13_G0000527", "CHM13_G0055907", "CHM13_G0002961")
NEED, LEFT_OUT, WITHIN, TOOLS = 8, 4, 0.10, ("Rustle", "StringTie")
WORD = {True: "pass", False: "fail", None: "not given"}

def overlap(k1, n1, k2, n2):
    a, b = wilson(k1, n1), wilson(k2, n2)
    return a[0] <= b[1] and b[0] <= a[1]

def same(k1, n1, k2, n2):
    """two shares k / n -> (pass, why)"""
    if overlap(k1, n1, k2, n2):
        return True, "overlap"
    small, big = sorted((k1 * n2, k2 * n1))
    if not small:
        return False, "disjoint, one share zero"
    ok = 2 * big < 3 * small            # larger / smaller < 1.5 in integers: 0.15 / 0.1 is below 1.5 in floats
    return ok, f"factor {big / small:.2f}" if ok else f"disjoint, factor {big / small:.2f}"

def close(a, b):
    return round(abs(a - b), 9) <= WITHIN    # 0.809 - 0.709 is above 0.1 in floats

def cell(kn):
    if not kn or not kn[1]:
        return "missing"
    lo, hi = wilson(*kn)
    return f"{kn[0]}/{kn[1]} {kn[0] / kn[1]:.4f} [{lo:.4f}, {hi:.4f}]"

def load(prefix):
    """-> (stats rows by name, {first gene id of the family: (near_ties, reads)}; None without a family table)"""
    with open(prefix + ".stats.tsv") as fh:
        stats = {r["stat"]: r for r in csv.DictReader(fh, delimiter="\t")}
    if not os.path.exists(prefix + ".fam.tsv"):
        return stats, None
    with open(prefix + ".fam.tsv") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    fam = {r["family"].split(",")[0].split()[0]: (int(r["near_ties"]), int(r["reads"])) for r in rows}
    assert len(fam) == len(rows), f"{prefix}.fam.tsv: two families with the same first gene id"
    return stats, fam

def pooled(fam):
    rest = [v for g, v in fam.items() if g not in FAMILIES[:LEFT_OUT]]
    return sum(k for k, _ in rest), sum(n for _, n in rest)

def item_f(seen, held, sim):
    rows = []
    for i, g in enumerate(FAMILIES, 1):
        ok, why = same(*sim[g], *held[g]) if g in sim and g in held else (False, "missing")
        rows.append((f"F{i}", f"{g} {why}", cell(seen.get(g)), cell(held.get(g)), cell(sim.get(g)), WORD[ok]))
    n = sum(r[5] == "pass" for r in rows)
    return [("F", f"near-tie share, {n} of {len(FAMILIES)} families (need {NEED})", "", "", "", WORD[n >= NEED])] + rows

def item_s(seen, held, sim):
    ok, why = same(*pooled(sim), *pooled(held))
    return [("S", f"pooled near-tie share without the first {LEFT_OUT} families; {why}", cell(pooled(seen)), cell(pooled(held)), cell(pooled(sim)), WORD[ok])]

def item_w(seen, held, sim):
    def kn(stats):                      # 0 of 0 has the interval 0-1 and would pass against anything
        r = stats.get("S2_decidable_primary_wrong")
        return (int(r["k"]), int(r["n"])) if r and r["n"] not in ("", "0") else None
    s, h, m = kn(seen), kn(held), kn(sim)
    return [("W", "primary off its origin among decidable near-tied reads", cell(s), cell(h), cell(m), WORD[overlap(*m, *h) if h and m else None])]

def both_ends(x):
    """a number, or both_ends_strict of a read_stats.py --out table"""
    try:
        return float(x)
    except ValueError:
        with open(x) as fh:
            return float(dict(l.rstrip("\n").split("\t")[:2] for l in fh)["both_ends_strict"])

def within(item, what, real, sim):
    return (item, f"{what}, |sim - real| {abs(sim - real):.4f} (<= {WITHIN:.2f})", "", f"{real:.4f}", f"{sim:.4f}", WORD[close(real, sim)])

def item_a(replay):
    given = {tool.lower(): (tool, *map(float, v.split(":"))) for tool, v in (x.split("=") for x in replay)}
    rows = [within(f"A{i}", f"{tool} R", real, sim) for i, (tool, real, sim) in enumerate(given.values(), 1)]
    rows += [(f"A{len(rows) + i}", f"{t} R", "", "", "", WORD[None]) for i, t in enumerate((t for t in TOOLS if t.lower() not in given), 1)]
    v = {r[5] for r in rows}            # the registration names both tools: one alone cannot pass
    return [("A", "replay of row 1794, R of every tool", "", "", "", "fail" if "fail" in v else WORD[None] if WORD[None] in v else "pass")] + rows

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for x in ("seen", "held", "sim"):
        ap.add_argument("--" + x, required=True, metavar="PREFIX")
    ap.add_argument("--ends", nargs=2, metavar=("REAL_HELD", "SIM")); ap.add_argument("--replay", action="append", default=[], metavar="TOOL=REAL:SIM")
    a = ap.parse_args()
    (st_seen, seen), (st_held, held), (st_sim, sim) = (load(p) for p in (a.seen, a.held, a.sim))
    rows = []
    for item, what, fn in (("F", "near-tie share of the registered families", item_f), ("S", "pooled near-tie share", item_s)):
        rows += fn(seen or {}, held, sim) if held is not None and sim is not None else [(item, what + ": no family table", "", "", "", WORD[None])]
    rows += item_w(st_seen, st_held, st_sim)
    rows += [within("E", "both ends complete (strict)", *map(both_ends, a.ends))] if a.ends else [("E", "both ends complete (strict)", "", "", "", WORD[None])]
    rows += item_a(a.replay)
    print("item\tdetail\tseen\theld\tsim\tverdict")
    for r in rows:
        print("\t".join(r))
    failed, missing = ([r[0] for r in rows if len(r[0]) == 1 and r[5] == w] for w in ("fail", WORD[None]))
    not_given = "not given " + ", ".join(missing)
    print("BAR fail: " + ", ".join(failed) + ("; " + not_given if missing else "") if failed else "BAR incomplete: " + not_given if missing else "BAR pass")

if __name__ == "__main__":
    main()
