#!/usr/bin/env python3
"""Separate REAL over-merges from impurity manufactured by overlapping Soto truth intervals.

`union_ab_score.py` credits a predicted copy to EVERY truth member it overlaps. That is the right rule for
recall, but it inflates impurity: 16 pairs of Soto truth members from DIFFERENT families overlap each other
on the genome (AC244669.1/ID_212 sits inside AC244669.2/ID_215 for 52 kb), so a single correctly-placed copy
in such a zone is attributed to two families and its predicted family is scored impure no matter what the
pipeline did. That is a property of the reference catalog, not of our clustering.

This re-scores with each copy assigned to ONE truth member -- the one it overlaps by the most base pairs --
and reports:

  REAL      predicted families still spanning >1 truth family after best-overlap assignment
  ARTIFACT  families impure only under the credit-everything rule

A real over-merge means two copies in the same predicted family whose own best matches disagree. An artifact
means one copy that two truth families both claim. Only the first is our defect.
"""
import argparse, csv
from collections import defaultdict


def truth(bed):
    out = []
    for ln in open(bed):
        f = ln.rstrip("\n").split("\t")
        if len(f) < 4 or "|" not in f[3]:
            continue
        name, fam = f[3].split("|")[0], f[3].split("|")[1]
        out.append((fam, name, f[0], int(f[1]), int(f[2])))
    return out


def copies(path):
    out = []
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if not r.get("chrom"):
                continue
            try:
                out.append((r["family_id"], r["chrom"], int(r["start"]), int(r["end"])))
            except (ValueError, KeyError):
                continue
    return out


def analyse(label, cps, tr):
    by_chrom = defaultdict(list)
    for t in tr:
        by_chrom[t[2]].append(t)

    all_fams = defaultdict(set)   # credit-everything (union_ab_score rule)
    best_fams = defaultdict(set)  # best-overlap only
    contested = []                # copies claimed by >1 truth family
    for fid, c, s, e in cps:
        hits = []
        for tfam, tname, tc, ts, te in by_chrom.get(c, ()):
            ov = min(e, te) - max(s, ts)
            if ov > 0:
                hits.append((ov, tfam, tname))
        if not hits:
            continue
        for _, tfam, _ in hits:
            all_fams[fid].add(tfam)
        ov, tfam, tname = max(hits)
        best_fams[fid].add(tfam)
        if len({h[1] for h in hits}) > 1:
            contested.append((fid, c, s, e, sorted({(h[1], h[2]) for h in hits})))

    impure_all = {k for k, v in all_fams.items() if len(v) > 1}
    impure_best = {k for k, v in best_fams.items() if len(v) > 1}
    artifact = impure_all - impure_best
    return {
        "label": label,
        "impure_all": impure_all,
        "impure_best": impure_best,
        "artifact": artifact,
        "all_fams": all_fams,
        "best_fams": best_fams,
        "contested": contested,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--copies", nargs="+", required=True, help="label=path/to/copies.tsv")
    ap.add_argument("--bed", default="bench/soto/80_fams.chr.bed")
    a = ap.parse_args()
    tr = truth(a.bed)

    print(f"{'run':<16}{'impure (credit-all)':>21}{'impure (best-overlap)':>23}{'ARTIFACT':>10}"
          f"{'worst REAL':>12}")
    print("-" * 82)
    results = []
    for spec in a.copies:
        label, path = spec.split("=", 1)
        r = analyse(label, copies(path), tr)
        worst_real = max((len(r["best_fams"][k]) for k in r["impure_best"]), default=1)
        results.append(r)
        print(f"{label:<16}{len(r['impure_all']):>21}{len(r['impure_best']):>23}"
              f"{len(r['artifact']):>10}{worst_real:>12}")

    for r in results:
        if not r["impure_best"]:
            print(f"\n{r['label']}: NO real over-merge survives best-overlap assignment.")
            for fid in sorted(r["artifact"]):
                print(f"  {fid} scored impure only via contested copies: {sorted(r['all_fams'][fid])}")
            continue
        print(f"\n{r['label']}: {len(r['impure_best'])} REAL over-merge(s)")
        for fid in sorted(r["impure_best"], key=lambda k: -len(r["best_fams"][k]))[:6]:
            print(f"  {fid}: {sorted(r['best_fams'][fid])}")

    if results:
        print(f"\ncontested copies in {results[0]['label']} (one copy, two truth families claim it):")
        for fid, c, s, e, claims in results[0]["contested"][:10]:
            print(f"  {fid} {c}:{s:,}-{e:,} <- {claims}")


if __name__ == "__main__":
    main()
