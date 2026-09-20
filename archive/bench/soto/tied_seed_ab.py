#!/usr/bin/env python3
"""A/B the tied-secondary seeding gate on a single region.

The advisor's objection: for a multimapping read, which placement minimap2 flags PRIMARY is close to
arbitrary (AS ties), so restricting skeleton construction to primary alignments may discard real copies.
`copy_assign --tied-seed` is the lever that answers it: it admits AS-tied SECONDARY alignments that agree
on an intron chain, at loci carrying NO primary skeleton.

This scores both arms against the Soto truth windows and separates the two things the gate can do:

  RECOVERED  a truth member that arm B covers and arm A does not  -> the gate working as intended
  PHANTOM    a copy in arm B overlapping NO truth member          -> the failure mode

Both are needed. Member recall alone cannot distinguish them, which is how a partition-blind evaluation
previously reported a seeding change as a gain when it was an over-merge.
"""
import argparse, csv, sys
from collections import defaultdict


def load_copies(path):
    out = []
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if not row.get("chrom"):
                continue
            try:
                out.append({
                    "family": row.get("family_id", "?"),
                    "chrom": row["chrom"],
                    "start": int(row["start"]),
                    "end": int(row["end"]),
                    "n_exon": int(row.get("n_exon", 0) or 0),
                    "n_reads": int(row.get("n_reads", 0) or 0),
                })
            except (ValueError, KeyError):
                continue
    return out


def load_truth(bed, fam, chrom, lo, hi):
    out = []
    for ln in open(bed):
        f = ln.rstrip("\n").split("\t")
        if len(f) < 4 or not f[3].endswith("|" + fam):
            continue
        c, s, e = f[0], int(f[1]), int(f[2])
        if c == chrom and s >= lo and e <= hi:
            out.append((f[3].split("|")[0], c, s, e))
    return sorted(out, key=lambda x: x[2])


def covers(cp, t):
    _, c, s, e = t
    return cp["chrom"] == c and cp["start"] < e and cp["end"] > s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--a", required=True, help="copies.tsv from the default (primary-only) arm")
    ap.add_argument("--b", required=True, help="copies.tsv from the --tied-seed arm")
    ap.add_argument("--bed", default="bench/soto/80_fams.chr.bed")
    ap.add_argument("--family", default="ID_154")
    ap.add_argument("--region", required=True, help="chr:start-end")
    a = ap.parse_args()

    chrom, rng = a.region.split(":")
    lo, hi = (int(x.replace(",", "")) for x in rng.split("-"))
    truth = load_truth(a.bed, a.family, chrom, lo, hi)
    A, B = load_copies(a.a), load_copies(a.b)

    print(f"region {a.region}   truth members {len(truth)}   "
          f"copies: primary-only {len(A)}, tied-seed {len(B)}\n")

    hdr = f"{'truth member':<12}{'window kb':>10}{'primary-only':>15}{'tied-seed':>12}   verdict"
    print(hdr); print("-" * len(hdr))
    rec = 0
    for t in truth:
        ca = [c for c in A if covers(c, t)]
        cb = [c for c in B if covers(c, t)]
        v = ""
        if cb and not ca:
            v = "RECOVERED by tied-seed"; rec += 1
        elif ca and not cb:
            v = "LOST under tied-seed"
        elif ca and cb:
            v = "both"
        else:
            v = "missed by both"
        print(f"{t[0]:<12}{(t[3]-t[2])/1000:>10.1f}{len(ca):>15}{len(cb):>12}   {v}")

    def phantoms(cps):
        return [c for c in cps if not any(covers(c, t) for t in truth)]

    pa, pb = phantoms(A), phantoms(B)
    print(f"\ncopies overlapping NO truth member:")
    print(f"  primary-only  {len(pa):>3}")
    print(f"  tied-seed     {len(pb):>3}   (net new: {len(pb)-len(pa):+d})")
    new = [c for c in pb if not any(c["start"] < x["end"] and c["end"] > x["start"]
                                    and c["chrom"] == x["chrom"] for x in pa)]
    if new:
        print(f"\n  {len(new)} phantom candidate(s) present ONLY under --tied-seed:")
        for c in sorted(new, key=lambda x: x["start"])[:15]:
            print(f"    {c['chrom']}:{c['start']:,}-{c['end']:,}  "
                  f"{c['n_exon']} exon  {c['n_reads']} reads  {c['family']}")
    print(f"\nSUMMARY: recovered {rec} truth member(s), introduced {len(new)} copy(ies) "
          f"outside every truth window.")
    if rec == 0 and new:
        print("=> the gate paid only costs on this region.")
    elif rec and not new:
        print("=> the gate paid only gains on this region.")


if __name__ == "__main__":
    main()
