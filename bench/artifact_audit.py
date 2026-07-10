#!/usr/bin/env python3
"""Audit a copy_assign catalog for COPY-SET construction artifacts.

Motivation: `CAFAM69` in the 74-family run abstained on 100% of its 1043 reads with `min_p == 1.0`
for every read, despite the family reporting 1204 PSV columns. That only reconciles if two of its
"copies" are the same object -- and they are. Its copies span 164381222-164384848 and
164381237-164384845: the same 3.6 kb interval, offset by 15 bp of boundary wobble. One locus, admitted
twice; the first holds all 1534 reads, the second holds zero. Only the third copy (164442447-164446101,
58 kb away) is a genuinely distinct locus. Those reads are not unassignable -- the copy set is malformed.

Root cause: `family_detect::collapse_loci_groups` unions two de-novo transcripts only when they share an
EXACT intron `(chrom, donor, acceptor)`. It never consults positional overlap, so transcripts lying on top
of each other with wobbled boundaries survive as separate "loci".

Three checks, in decreasing strength:

  1. OVERLAPPING COPIES (structural, NON-CIRCULAR, decisive).
     Two copies of one family whose genomic spans overlap are, by definition, not two loci. Needs no
     annotation and no reference truth -- it is a self-consistency property of the emitted catalog.

  2. DUPLICATE FAMILIES (structural, non-circular).
     The same copy set emitted under two family_ids. Caused by OVERLAPPING INPUT REGIONS in the sweep
     (each region independently detects the same family). Inflates every per-family statistic.

  3. ANNOTATION CROSS-CHECK (corroborating, WEAKER -- read the caveat).
     Two copies landing inside the same annotated gene is *suggestive* of a same-locus artifact. It is
     NOT proof: (a) a genuine COLLAPSED copy legitimately lives inside one annotated gene -- that is
     precisely what O4 discovers; (b) the annotation may merge two real tandem paralogs into a single
     gene record, so "same gene" can mean "annotation under-splits", not "we over-split". Conversely,
     two copies in DISTINCT annotated genes corroborate a real family. Treat check 3 as evidence that
     agrees or disagrees with check 1, never as the verdict.

Usage:
    python bench/artifact_audit.py <prefix>            # reads <prefix>.quant.tsv
    python bench/artifact_audit.py <prefix> --gff GGO_genomic.gff

`<prefix>.quant.tsv` must carry `copy_chrom`, `copy_start`, `copy_end` (emitted by copy_assign).
Exit status is non-zero when any structural (check 1 or 2) artifact is found, so this can gate a sweep.
"""
import argparse
import bisect
import collections
import csv
import re
import sys


def load_copies(quant_path):
    """family_id -> [(copy_tid, chrom, start, end, n_reads_hard)]. Requires span columns."""
    fam = collections.defaultdict(list)
    with open(quant_path) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        need = {"copy_chrom", "copy_start", "copy_end"}
        missing = need - set(rd.fieldnames or [])
        if missing:
            sys.exit(
                f"{quant_path} lacks {sorted(missing)}. Re-run copy_assign with a build that emits copy spans."
            )
        for r in rd:
            fam[r["family_id"]].append(
                (
                    r["copy_tid"],
                    r["copy_chrom"],
                    int(r["copy_start"]),
                    int(r["copy_end"]),
                    int(r.get("n_reads_hard", 0) or 0),
                )
            )
    return fam


def overlapping_copies(copies):
    """Pairs of copies whose spans intersect, with the RECIPROCAL overlap fraction.

    Copies of one family occupy disjoint loci by definition, so any intersection is a defect. The
    reciprocal fraction (overlap / longest span) says which one:
      ~1.0  -> DUPLICATE LOCUS: the same interval admitted twice (boundary wobble). Every read scores
               min_p == 1 against the pair, so the family abstains wholesale and its reads look like K=0.
      <<1.0 -> CONTAINMENT: a long readthrough transcript enclosing a short one. Must NOT be merged --
               merging would let the readthrough absorb genuinely distinct tandem copies.
    """
    bad = []
    for i in range(len(copies)):
        for j in range(i + 1, len(copies)):
            a, b = copies[i], copies[j]
            if a[1] != b[1]:
                continue
            lo, hi = max(a[2], b[2]), min(a[3], b[3])
            if lo < hi:
                longest = max(a[3] - a[2], b[3] - b[2])
                recip = (hi - lo) / longest if longest else 1.0
                kind = "DUPLICATE LOCUS" if recip > 0.9 else "CONTAINMENT (readthrough?)"
                bad.append((a[0], b[0], hi - lo, recip, kind))
    return bad


def duplicate_families(fam):
    """family_ids sharing an identical copy-tid set."""
    sig = collections.defaultdict(list)
    for fid, copies in fam.items():
        sig[tuple(sorted(c[0] for c in copies))].append(fid)
    return {k: v for k, v in sig.items() if len(v) > 1}


def load_genes(gff):
    genes = collections.defaultdict(list)
    with open(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            m = re.search(r"gene=([^;]*)", f[8])
            genes[f[0]].append((int(f[3]), int(f[4]), m.group(1) if m else "?"))
    for c in genes:
        genes[c].sort()
    return genes


def genes_overlapping(genes, chrom, start, end):
    """All annotated genes intersecting [start, end)."""
    lst = genes.get(chrom, [])
    i = bisect.bisect_right(lst, (end, float("inf"), ""))
    return sorted({g for (s, e, g) in lst[:i] if s < end and e > start})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("prefix")
    ap.add_argument("--gff", help="annotation for the corroborating cross-check (optional)")
    args = ap.parse_args()

    fam = load_copies(f"{args.prefix}.quant.tsv")
    print(f"catalog: {len(fam)} families, {sum(len(v) for v in fam.values())} copies\n")

    # --- check 1: overlapping copies (decisive) ---
    print("[1] OVERLAPPING COPIES within a family (structural, non-circular)")
    n_bad = 0
    for fid in sorted(fam):
        for a, b, ov, recip, kind in overlapping_copies(fam[fid]):
            n_bad += 1
            print(f"    {fid:10s} {kind:26s} recip={recip:.2f}  {a} <-> {b}  ({ov} bp)")
    print(f"    => {n_bad} overlapping copy pair(s)\n")

    # --- check 2: duplicate families (structural) ---
    print("[2] DUPLICATE FAMILIES (identical copy sets under >1 family_id)")
    dups = duplicate_families(fam)
    for cs, fids in sorted(dups.items(), key=lambda kv: kv[1]):
        print(f"    {fids}  ({len(cs)} copies)")
    print(f"    => {len(dups)} duplicated family/families spanning {sum(len(v) for v in dups.values())} ids")
    if dups:
        print("       cause: overlapping regions in the input region list; each detects the family once")
    print()

    # --- check 3: annotation cross-check (corroborating only) ---
    if args.gff:
        genes = load_genes(args.gff)
        print("[3] ANNOTATION CROSS-CHECK (corroborating, not decisive -- see module docstring)")
        agree = collections.Counter()
        for fid in sorted(fam):
            hits = [genes_overlapping(genes, c[1], c[2], c[3]) for c in fam[fid]]
            named = [h[0] for h in hits if h]
            if len(named) < 2:
                agree["under-annotated (cannot adjudicate)"] += 1
                continue
            overlaps = overlapping_copies(fam[fid])
            if len(set(named)) < len(named):
                # same gene twice: artifact IF the spans also overlap; otherwise it is either a genuine
                # collapsed copy or an annotation that under-splits a tandem pair.
                if overlaps:
                    agree["same gene AND overlapping spans -> ARTIFACT (both checks agree)"] += 1
                else:
                    agree["same gene, DISJOINT spans -> collapsed copy or under-split annotation"] += 1
            else:
                agree["distinct genes -> family corroborated"] += 1
        for k, v in agree.most_common():
            print(f"    {v:3d}  {k}")
        print()

    if n_bad or dups:
        print("STRUCTURAL ARTIFACTS PRESENT -- these inflate abstention and every per-family statistic.")
        return 1
    print("No structural artifacts.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
