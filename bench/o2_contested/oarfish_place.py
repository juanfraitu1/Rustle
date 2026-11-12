#!/usr/bin/env python3
"""r1648 (docs/PREREG_oarfish_placement_2026-11-12.md): picks from an oarfish `.prob` file. A read's posterior is summed
over all sequences of one locus (copies with overlapping spans; OTHER is its own class); assign to the top locus iff
its sum >= THRESH (default 0.99) and it is not OTHER.
Usage: oarfish_place.py OUT.prob COPIES.tsv [NAMES.txt: only these reads] [THRESH] > OUT.picks.tsv"""
import collections, csv, sys

def main(prob, copies, names=None, thresh="0.99"):
    span = {("OTHER", "0"): ("OTHER", 0, 1)}; thresh = float(thresh)
    for r in csv.DictReader(open(copies), delimiter="\t"): span[(r["family_id"], r["copy_idx"])] = (r["chrom"], int(r["start"]), int(r["end"]))
    keep = {l.strip() for l in open(names)} if names else None
    same = lambda a, b: span[a][0] == span[b][0] and span[a][1] < span[b][2] and span[b][1] < span[a][2]
    n = collections.Counter()
    with open(prob) as fh:
        m, _ = map(int, fh.readline().split())
        tx = [("OTHER", "0") if x[0] == "OTHER" else (x[0], x[1]) for x in (fh.readline().rstrip("\n").split("|") for _ in range(m))]
        for line in fh:
            f = line.split("\t"); q = f[0][2:] if f[0].startswith("R|") else f[0]
            if keep is not None and q not in keep: continue
            k = int(f[1]); p = collections.Counter()
            for i, w in zip(f[2:2 + k], f[2 + k:2 + 2 * k]): p[tx[int(i)]] += float(w)
            top = max(p, key=p.get); s = sum(w for c, w in p.items() if c == top or same(c, top)); n["reads"] += 1
            if top[0] == "OTHER": n["abstain_other"] += 1
            elif s < thresh: n["abstain_below"] += 1
            else: n["assigned"] += 1; print(q, top[0], top[1], sep="\t")
    print(dict(n), file=sys.stderr)

if __name__ == "__main__": main(*sys.argv[1:5])
