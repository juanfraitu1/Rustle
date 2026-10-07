#!/usr/bin/env python3
"""E1 of docs/PREREG_ideal_expression_2026-10-06.md (descriptive): per target copy, the distinct canonical multi-exon chains that appear exactly as a transcript chain in PREFIX.gtf
(any assembled transcript, same chromosome and strand), over all chains and over OBSERVABLE chains (>= 3 reads of the seeding pool P2, STRATAPREFIX.chains.tsv).

    chains.py --truth TRUTHPREFIX --strata STRATAPREFIX --asm ASMPREFIX --out OUT
Writes OUT.E1.tsv and prints one line per family.
"""
import argparse
import collections
import csv
import re


def chain_of(exons):
    ex = sorted(exons)
    return tuple((ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--strata", required=True)
    ap.add_argument("--asm", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rd = lambda p: list(csv.DictReader(open(p), delimiter="\t"))
    targets = {t["cid"]: t for t in rd(a.truth + ".targets.tsv")}
    ex, meta = collections.defaultdict(list), {}
    for ln in open(a.asm + ".gtf"):
        f = ln.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "exon":
            continue
        t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
        meta[t] = (f[0], f[6]); ex[t].append((int(f[3]) - 1, int(f[4])))
    have = set()
    for t, e in ex.items():
        have.add((meta[t][0], meta[t][1], chain_of(e)))
    per = collections.defaultdict(lambda: [0, 0, 0, 0])   # chains, recovered, observable chains, recovered among observable
    for r in rd(a.strata + ".chains.tsv"):
        t = targets[r["cid"]]
        ch = tuple(tuple(map(int, x.split("-"))) for x in r["chain"].split(","))
        got = (t["chrom"], t["strand"], ch) in have
        c = per[r["cid"]]
        c[0] += 1; c[1] += int(got)
        if int(r["observable"]):
            c[2] += 1; c[3] += int(got)
    rows = [dict(cid=cid, name=targets[cid]["name"], family=targets[cid]["family"], chains=c[0], recovered=c[1], observable_chains=c[2], recovered_observable=c[3])
            for cid, c in sorted(per.items())]
    with open(a.out + ".E1.tsv", "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t"); w.writeheader(); w.writerows(rows)
    for fam in sorted({r["family"] for r in rows}):
        rs = [r for r in rows if r["family"] == fam]
        print(f"E1 {fam}: chains {sum(r['chains'] for r in rs)}, recovered exactly {sum(r['recovered'] for r in rs)}; observable chains {sum(r['observable_chains'] for r in rs)}, "
              f"recovered among them {sum(r['recovered_observable'] for r in rs)}; copies with every observable chain recovered {sum(1 for r in rs if r['recovered_observable'] == r['observable_chains'])}/{len(rs)}")


if __name__ == "__main__":
    main()
