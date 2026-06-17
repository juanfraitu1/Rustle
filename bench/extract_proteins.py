#!/usr/bin/env python3
"""Translate the longest CDS per annotated gene -> protein FASTA (DNA/protein tier input).
Proteins stay conserved when nucleotide sequence has diverged, so protein clustering catches the
ancient families the RNA-similarity definition misses. Run with /home/juanfra/miniforge3/bin/python."""
import glob
import os
from collections import defaultdict

import pysam
from Bio.Seq import Seq

GW = "/tmp/gw"
FASTA = "/home/juanfra/winloci_scratch/GGO.fasta"
OUT = "/home/juanfra/winloci_scratch/proteins.fa"


def attrs(f):
    d = {}
    for kv in f.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1); d[k.strip()] = v.strip()
    return d


def main():
    feat_gene = {}; gmeta = {}; cds = defaultdict(list)
    for gff in sorted(glob.glob(os.path.join(GW, "ref_*.gff3"))):
        chrom = os.path.basename(gff)[len("ref_"):-len(".gff3")]
        for line in open(gff):
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            a = attrs(f[8])
            if f[2] == "gene":
                g = a.get("gene") or (a.get("ID", "")[5:] if a.get("ID", "").startswith("gene-") else a.get("ID"))
                if g:
                    gmeta.setdefault(g, (chrom, f[6]))
            elif f[2] in ("mRNA", "transcript"):
                if a.get("ID") and a.get("gene"):
                    feat_gene[(chrom, a["ID"])] = a["gene"]
            elif f[2] == "CDS":
                p = a.get("Parent", "")
                if (chrom, p) in feat_gene:
                    cds[(chrom, p)].append((int(f[3]) - 1, int(f[4])))
    # longest CDS per gene
    best = {}
    for (chrom, rna), ivs in cds.items():
        g = feat_gene[(chrom, rna)]
        tot = sum(e - s for s, e in ivs)
        if g not in best or tot > best[g][2]:
            best[g] = (chrom, sorted(ivs), tot)
    fa = pysam.FastaFile(FASTA)
    n = 0
    with open(OUT, "w") as fh:
        for g in sorted(best):
            chrom, ivs, _ = best[g]
            strand = gmeta.get(g, (chrom, "+"))[1]
            nt = "".join(fa.fetch(chrom, s, e).upper() for s, e in ivs)
            if strand == "-":
                nt = str(Seq(nt).reverse_complement())
            if len(nt) < 90:
                continue
            prot = str(Seq(nt[: len(nt) // 3 * 3]).translate()).rstrip("*").replace("*", "X")
            if len(prot) < 30:
                continue
            fh.write(f">{g}\n")
            for i in range(0, len(prot), 60):
                fh.write(prot[i:i+60] + "\n")
            n += 1
    print(f"wrote {n} proteins -> {OUT}")


if __name__ == "__main__":
    main()
