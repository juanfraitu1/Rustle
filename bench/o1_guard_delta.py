#!/usr/bin/env python3
"""TWO-SIDED orientation-guard ledger: what a catalog run GAINS and what it LOSES.

WHY TWO-SIDED. `o1_guard_cost.py` counts only guard-BLOCKED pairs (a qualifying minus record, no
qualifying plus), which measures what a strand fix could RECOVER. That is half the ledger. 1,884/4,778 =
0.3943 of currently-KEPT pairs already involve a single-exon rep, and read-orientation voting is measured
at 0.9650, not 1.0 -- so a wrong vote on a currently-kept pair flips it from '+'-qualifying to
minus-only and LOSES an edge. At ~3.5% error those losses are expected, not hypothetical, and a one-sided
count would read the delta as a pure win.

Usage:  o1_guard_delta.py <off_paf> <off_nodes> <on_paf> <on_nodes> [frozen58.tsv]

Pairs are keyed by `node_key`, never by `idx`: the rep set MOVES under RUSTLE_READ_STRAND (span-aware
locus collapse gates on `a.strand != b.strand` at family_detect.rs:670, inside the skeletons->reps step),
so idx is not a relabelling and any idx-keyed comparison silently compares different objects.
"""
import csv, sys


def load(paf, nodes):
    nex, key = {}, {}
    for r in csv.DictReader(open(nodes), delimiter="\t"):
        i = int(r["idx"]); nex[i] = int(r["n_exon"]); key[i] = r["node_key"]
    plus, minus = set(), set()
    for line in open(paf):
        f = line.rstrip("\n").split("\t")
        if len(f) < 12:
            continue
        q, t = int(f[0]), int(f[5])
        if q == t or q not in key or t not in key:
            continue
        ql, qs, qe = int(f[1]), int(f[2]), int(f[3])
        tl, ts, te = int(f[6]), int(f[7]), int(f[8])
        de = next((float(x[5:]) for x in f[12:] if x.startswith("de:f:")), None)
        ident = (1.0 - de) if de is not None else int(f[9]) / max(int(f[10]), 1)
        cov = (qe - qs) / max(ql, 1) if ql <= tl else (te - ts) / max(tl, 1)   # shipped form; axis follows denominator
        if ident < 0.60 or cov < 0.50:
            continue
        k = tuple(sorted((key[q], key[t])))
        (plus if f[4] == "+" else minus).add(k)
    return plus, minus - plus, nex, key


def main():
    off_p, off_b, _, _ = load(sys.argv[1], sys.argv[2])
    on_p,  on_b,  _, _ = load(sys.argv[3], sys.argv[4])
    print(f"OFF: kept {len(off_p)}  guard-blocked {len(off_b)}")
    print(f"ON : kept {len(on_p)}  guard-blocked {len(on_b)}\n")

    gained = on_p - off_p          # blocked (or absent) in OFF, kept in ON
    lost   = off_p - on_p          # kept in OFF, no longer kept in ON  <-- the arm that was missing
    recovered = off_b & on_p       # specifically: was guard-blocked, now passes
    newblock  = off_p & on_b       # specifically: was kept, now guard-blocked
    print(f"  GAINED (kept in ON, not in OFF)      {len(gained)}")
    print(f"    ...of which were guard-BLOCKED off  {len(recovered)}")
    print(f"  LOST   (kept in OFF, not in ON)      {len(lost)}")
    print(f"    ...of which are now guard-BLOCKED   {len(newblock)}")
    net = len(gained) - len(lost)
    print(f"\n  NET edge-pair change: {net:+d}   (report BOTH arms, never the net alone)")

    if len(sys.argv) > 5:
        frozen = {tuple(sorted((r["node_key_a"], r["node_key_b"])))
                  for r in csv.DictReader(open(sys.argv[5]), delimiter="\t")}
        admitted = frozen & on_p
        print(f"\n  ACCEPTANCE — the {len(frozen)} frozen GENUINE-antisense pairs (both reps junction-determined)")
        print(f"    admitted by ON: {len(admitted)}/{len(frozen)}   "
              f"{'PASS - all still rejected' if not admitted else 'FAIL - the guard let real antisense through'}")
        for k in sorted(admitted)[:5]:
            print(f"      leaked: {k[0]}  x  {k[1]}")


if __name__ == "__main__":
    main()
