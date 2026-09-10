#!/usr/bin/env python3
"""PREREG hard_locus_bakeoff (5ca5c7e4): compare per-tool derived calls (tool_bakeoff.py --out *.calls.tsv)
on the HARD set (the AS-tied gate's molecules = rows of an assignments.tsv) vs the EASY set (everything else
in the family region), per stratum, per copy, and copy attribution of the O2-assigned molecules.

  python3 bench/hard_locus_bakeoff.py --assign ours_final.assignments.tsv --gtf ours.gtf \
      ours=hard/ours.calls.tsv flair=hard/flair.calls.tsv stringtie=hard/stringtie.calls.tsv isoseq=hard/isoseq.calls.tsv
"""
import argparse, csv, re
from collections import Counter, defaultdict


def load_calls(p):
    out = {}
    with open(p) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[r["molecule"]] = (r["state"], [int(x) for x in r["copies"].split(",") if x])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--assign", required=True)
    ap.add_argument("--gtf", help="our GTF, for copy_index of the transcript carrying each molecule (P6)")
    ap.add_argument("--bam", help="with --min-mult: primaries (-F 2308) give each molecule's intron chain")
    ap.add_argument("--min-mult", type=int, default=0, help="keep only molecules whose exact chain is carried by >= N molecules (support-policy control)")
    ap.add_argument("tools", nargs="+", help="label=calls.tsv; the first is ours")
    a = ap.parse_args()
    tools = [(t.split("=")[0], load_calls(t.split("=")[1])) for t in a.tools]
    ours_l, ours = tools[0]
    assign = {r["read_name"]: r for r in csv.DictReader(open(a.assign), delimiter="\t")}
    allm = set(ours)
    for _, c in tools:
        allm |= set(c)
    if a.min_mult and a.bam:
        import re as _re, subprocess
        def introns(pos, cig):
            o, p = [], pos
            for n, op in _re.findall(r"(\d+)([MIDNSHP=X])", cig):
                n = int(n)
                if op in "M=XD":
                    p += n
                elif op == "N":
                    o.append((p, p + n)); p += n
            return tuple(o)
        chain = {}
        out = subprocess.run(["samtools", "view", "-F", "2308", a.bam], capture_output=True, text=True).stdout
        for ln in out.splitlines():
            f = ln.split("\t", 6)
            chain.setdefault(f[0], (f[2],) + introns(int(f[3]) - 1, f[5]))
        mult = Counter(chain.values())
        keep = {m for m in allm if mult.get(chain.get(m), 0) >= a.min_mult}
        print(f"--min-mult {a.min_mult}: keeping {len(keep)} of {len(allm)} molecules whose exact chain has >= {a.min_mult} molecules")
        allm = keep
    hard = {m for m in allm if m in assign}
    easy = allm - hard
    carried = lambda c, m: c.get(m, ("derived_none", []))[0] != "derived_none"
    strata = {
        "hard (all gate rows)": hard,
        "  contested": {m for m in hard if assign[m]["origin_rejected"] == "0" and int(assign[m]["n_candidates"]) >= 2},
        "    assigned": {m for m in hard if assign[m]["status"] == "assigned" and assign[m]["origin_rejected"] == "0" and int(assign[m]["n_candidates"]) >= 2},
        "    tied": {m for m in hard if assign[m]["status"] == "tied" and assign[m]["origin_rejected"] == "0" and int(assign[m]["n_candidates"]) >= 2},
        "    ambiguous": {m for m in hard if assign[m]["status"] == "ambiguous" and assign[m]["origin_rejected"] == "0" and int(assign[m]["n_candidates"]) >= 2},
        "  tie outside catalog": {m for m in hard if assign[m].get("tie_outside_catalog") == "1"},
        "  origin-rejected": {m for m in hard if assign[m]["origin_rejected"] == "1"},
        "easy (not admitted by the gate)": easy,
    }
    print(f"molecules: {len(allm)} total, hard {len(hard)}, easy {len(easy)}")
    print(f"\n{'stratum':34s} {'n':>6} " + " ".join(f"{l:>10}" for l, _ in tools) + "   (fraction carried = derived_one|multi)")
    rates = {}
    for name, S in strata.items():
        row = []
        for l, c in tools:
            k = sum(1 for m in S if carried(c, m))
            rates[(name, l)] = k / len(S) if S else float("nan")
            row.append(f"{k/len(S):10.3f}" if S else f"{'-':>10}")
        print(f"{name:34s} {len(S):>6} " + " ".join(row))
    # P2 discordance on hard and contested
    for name in ("hard (all gate rows)", "  contested"):
        S = strata[name]
        print(f"\n{name.strip()}: discordance vs ours")
        for l, c in tools[1:]:
            on = sum(1 for m in S if carried(ours, m) and not carried(c, m))
            xn = sum(1 for m in S if carried(c, m) and not carried(ours, m))
            both = sum(1 for m in S if carried(ours, m) and carried(c, m))
            nei = len(S) - on - xn - both
            print(f"   {l:10s} ours-not-{l}: {on:5d}   {l}-not-ours: {xn:5d}   ratio {on/max(1,xn):5.2f}   both {both}  neither {nei}")
    # P3 hard vs easy derived_none
    print("\nP3 derived_none hard / easy:")
    for l, c in tools:
        h = 1 - rates[("hard (all gate rows)", l)]; e = 1 - rates[("easy (not admitted by the gate)", l)]
        print(f"   {l:10s} hard {h:.3f}  easy {e:.3f}  ratio {h/e if e else float('nan'):.2f}")
    # P4 per copy: copies where ours carries >=1 hard molecule and X carries none (by the derived copy)
    print("\nP4 per-copy coverage on the hard set (copies where the tool carries >=1 hard molecule, by derived copy):")
    cov = {}
    for l, c in tools:
        cs = set()
        for m in hard:
            st, cps = c.get(m, ("derived_none", []))
            if st != "derived_none":
                cs.update(cps)
        cov[l] = cs
    for l, _ in tools[1:]:
        print(f"   {l:10s} copies {len(cov[l]):2d}; ours-only {sorted(cov[ours_l]-cov[l])}  {l}-only {sorted(cov[l]-cov[ours_l])}")
    print(f"   {ours_l:10s} copies {len(cov[ours_l])}: {sorted(cov[ours_l])}")
    # P6 copy attribution of the O2-assigned
    asg = strata["    assigned"]
    if a.gtf:
        tx_copy = {}
        for line in open(a.gtf):
            if "\ttranscript\t" not in line:
                continue
            t = re.search(r'transcript_id "([^"]+)"', line); ci = re.search(r'copy_index "([^"]+)"', line)
            if t and ci:
                tx_copy[t.group(1)] = ci.group(1)
    print(f"\nP6 the {len(asg)} O2-assigned molecules (report):")
    for l, c in tools:
        one = [m for m in asg if c.get(m, ("derived_none", []))[0] == "derived_one"]
        agree = sum(1 for m in one if str(c[m][1][0]) == assign[m]["catalog_copy_idx"])
        print(f"   {l:10s} carried {sum(1 for m in asg if carried(c, m)):3d}/{len(asg)}; derived_one {len(one):3d}, of which derived copy == O2 copy: {agree} ({100*agree/max(1,len(one)):.0f}%)")


if __name__ == "__main__":
    main()
