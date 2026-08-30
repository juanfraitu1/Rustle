#!/usr/bin/env python3
"""merge_precision_arms.py — price the merges a flag CREATES, against the pipeline's own output.

THE GAP THIS FILLS. `family_def_airtight_panel.py` scores edges recomputed from ANNOTATION-derived cDNA
homology and never invokes the pipeline, so it is structurally BLIND to any flag that changes locus
collapse or representative choice: both arms return "identical", which reads as a pass but is a null from
an instrument with zero qualifying candidates (the 2/150 trap). This script instead compares two REAL
catalog runs and asks, of every co-membership pair the flag creates, how many carry DIRECT E_r homology
evidence rather than arriving transitively through a hub.

WHY DIRECT-EDGE SUPPORT IS THE RIGHT STATISTIC. A family is an E_r component, so co-membership is
transitive by construction and a pair being co-family is NOT evidence the two copies resemble each other.
A small number of new edges can fuse large components and manufacture a quadratic number of new pairs.
The BASELINE arm's own direct-edge rate is the comparator: a flag whose new pairs are backed at roughly
the baseline rate is extending families the way the existing ones were built; one whose new pairs fall far
below it is fusing components through hubs.

⚠ Quote the rate WITH its comparator. A bare "15% direct" means nothing without the arm's own ~50%.

Usage: merge_precision_arms.py OFF_ARM_DIR ON_ARM_DIR
  each dir holding cat.copies.tsv and dump/e.edges.tsv
"""
import collections
import csv
import itertools
import sys


def key(c, s, e):
    return (c, int(s), int(e))


def load(arm):
    fam = collections.defaultdict(set)
    with open(f"{arm}/cat.copies.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            fam[r["family_id"]].add(key(r["chrom"], r["start"], r["end"]))
    edges = {}
    with open(f"{arm}/dump/e.edges.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            a = key(r["chrom_i"], r["start_i"], r["end_i"])
            b = key(r["chrom_j"], r["start_j"], r["end_j"])
            edges[tuple(sorted((a, b)))] = (float(r["identity"]), float(r["coverage"]))
    pairs = set()
    for members in fam.values():
        pairs.update(itertools.combinations(sorted(members), 2))
    return fam, edges, pairs


def main(off_dir, on_dir):
    fo, eo, po = load(off_dir)
    fn, en, pn = load(on_dir)
    new, kept = pn - po, pn & po
    rate = lambda ps, E: (sum(1 for p in ps if p in E), len(ps))

    d_off, n_off = rate(po, eo)
    d_keep, n_keep = rate(kept, en)
    d_new, n_new = rate(new, en)

    print(f"  edges         OFF {len(eo):>6}   ON {len(en):>6}   ({len(en)-len(eo):+d})")
    print(f"  co-membership OFF {len(po):>6}   ON {len(pn):>6}   ({len(pn)-len(po):+d})")
    print(f"  created {len(new)}   destroyed {len(po - pn)}")
    print("\n  fraction of co-membership pairs carrying a DIRECT E_r edge:")
    print(f"    OFF baseline      {d_off:>6}/{n_off:<6} = {d_off/max(n_off,1):.2%}   <- the comparator")
    print(f"    kept from OFF     {d_keep:>6}/{n_keep:<6} = {d_keep/max(n_keep,1):.2%}")
    print(f"    NEWLY created     {d_new:>6}/{n_new:<6} = {d_new/max(n_new,1):.2%}")

    famof = {x: f for f, m in fn.items() for x in m}
    conc = collections.Counter(famof.get(a) for a, _ in new)
    print(f"\n  new pairs concentrate in {len(conc)} families; top 3:")
    for f, n in conc.most_common(3):
        print(f"    {f}: {n} new pairs, family size {len(fn[f])}")
    if n_new and d_new / n_new < 0.5 * (d_off / max(n_off, 1)):
        print("\n  VERDICT: new pairs are backed FAR below the baseline rate — hub-driven component"
              "\n  fusion, not homology. Do not ship this flag on.")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
