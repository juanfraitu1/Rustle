#!/usr/bin/env python3
"""Arm B input of docs/PREREG_yag_o2_2026-10-12.md: the catalog's families with >= 1 copy on the Y, minus arrays of
> MAX_COPIES copies (the Yq12 satellite), their copies table / FASTA, and the catalog regions overlapping their copies.

    prep_real.py CATALOG_PREFIX YCHROM OUTDIR      (CATALOG_PREFIX.fam.copies.tsv/.fa, CATALOG_PREFIX.regions.txt)
"""
import collections, os, sys

MAX_COPIES = 100


def main(prefix, ychrom, out):
    os.makedirs(out, exist_ok=True)
    rows = [l.rstrip("\n").split("\t") for l in open(prefix + ".fam.copies.tsv")]
    head, rows = rows[0], rows[1:]
    ix = {h: i for i, h in enumerate(head)}
    by = collections.defaultdict(list)
    for r in rows:
        by[r[ix["family_id"]]].append(r)
    # Y-ONLY families (every copy on the Y): the ampliconic families; a family with copies elsewhere would pull whole
    # autosomes' reads into the run (the catalog regions are merged per chromosome)
    keep = {f for f, rs in by.items() if all(r[ix["chrom"]] == ychrom for r in rs) and 2 <= len(rs) <= MAX_COPIES}
    dropped = {f: len(rs) for f, rs in by.items() if any(r[ix["chrom"]] == ychrom for r in rs) and len(rs) > MAX_COPIES}
    with open(f"{out}/copies.tsv", "w") as fh:
        fh.write("\t".join(head) + "\n")
        for r in rows:
            if r[ix["family_id"]] in keep:
                fh.write("\t".join(r) + "\n")
    with open(prefix + ".fam.copies.fa") as fi, open(f"{out}/copies.fa", "w") as fo:
        on = False
        for line in fi:
            if line.startswith(">"):
                on = line[1:].split("|")[0] in keep
            if on:
                fo.write(line)
    spans = collections.defaultdict(list)
    for r in rows:
        if r[ix["family_id"]] in keep:
            spans[r[ix["chrom"]]].append((int(r[ix["start"]]), int(r[ix["end"]])))
    n = 0
    with open(prefix + ".regions.txt") as fi, open(f"{out}/regions.txt", "w") as fo:
        for line in fi:
            chrom, rng = line.strip().rsplit(":", 1)
            a, b = map(int, rng.split("-"))
            if any(s < b and a - 1 < e for s, e in spans.get(chrom, ())):
                fo.write(line)
                n += 1
    with open(f"{out}/ybed.bed", "w") as fh:  # the kept families' Y copies, for read lengths
        for r in rows:
            if r[ix["family_id"]] in keep and r[ix["chrom"]] == ychrom:
                fh.write(f"{ychrom}\t{r[ix['start']]}\t{r[ix['end']]}\n")
    print(f"{len(keep)} families kept ({sum(len(by[f]) for f in keep)} copies), dropped arrays {dropped}, {n} regions")


if __name__ == "__main__":
    main(*sys.argv[1:4])
