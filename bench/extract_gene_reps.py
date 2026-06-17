#!/usr/bin/env python3
"""Extract ONE representative spliced transcript sequence per gene, genome-wide, from the RefSeq
annotation + genome FASTA. Feeds the genome-wide cross-chromosome copy-discovery harness
(family_crosschrom_discovery.py). Representative = longest transcript (by summed exon length).

Output:
  /tmp/gene_reps_gw.fa        >GENE  (spliced, +strand-oriented mRNA sequence)
  /tmp/gene_reps_gw.meta.tsv  gene<TAB>chrom<TAB>start<TAB>end<TAB>strand<TAB>len
"""
import glob
import os
import re
from collections import defaultdict

import pysam

GW = "/tmp/gw"
FASTA = "/home/juanfra/winloci_scratch/GGO.fasta"
OUT_FA = "/tmp/gene_reps_gw.fa"
OUT_META = "/tmp/gene_reps_gw.meta.tsv"


def gff_attrs(field):
    d = {}
    for kv in field.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            d[k.strip()] = v.strip()
    return d


def main():
    fa = pysam.FastaFile(FASTA)
    # gene -> chrom,strand ; transcript feature ID -> (gene, exons[])
    gene_meta = {}                         # gene -> (chrom, strand)
    feat_gene = {}                         # (chrom, rna_feature_id) -> gene
    tx_exons = defaultdict(list)           # (chrom, rna_feature_id) -> [(s,e)]
    for gff in sorted(glob.glob(os.path.join(GW, "ref_*.gff3"))):
        chrom = os.path.basename(gff)[len("ref_"):-len(".gff3")]
        for line in open(gff):
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            typ = f[2]
            a = gff_attrs(f[8])
            if typ == "gene":
                g = a.get("gene") or (a.get("ID", "")[len("gene-"):] if a.get("ID", "").startswith("gene-") else a.get("ID"))
                if g:
                    gene_meta.setdefault(g, (chrom, f[6]))
            elif typ in ("mRNA", "transcript", "primary_transcript"):
                fid = a.get("ID", "")
                g = a.get("gene")
                if fid and g:
                    feat_gene[(chrom, fid)] = g
            elif typ == "exon":
                parent = a.get("Parent", "")
                if (chrom, parent) in feat_gene:
                    tx_exons[(chrom, parent)].append((int(f[3]), int(f[4])))

    # per gene: pick the transcript with max summed exon length
    gene_best = {}   # gene -> (chrom, exons, summed_len)
    for (chrom, fid), exons in tx_exons.items():
        g = feat_gene[(chrom, fid)]
        tot = sum(e - s + 1 for s, e in exons)
        if g not in gene_best or tot > gene_best[g][2]:
            gene_best[g] = (chrom, sorted(exons), tot)

    rc = bytes.maketrans(b"ACGTacgtNn", b"TGCAtgcaNn")
    n = 0
    with open(OUT_FA, "w") as fh, open(OUT_META, "w") as mh:
        mh.write("gene\tchrom\tstart\tend\tstrand\tlen\n")
        for g in sorted(gene_best):
            chrom, exons, _ = gene_best[g]
            strand = gene_meta.get(g, (chrom, "+"))[1]
            seq = []
            for s, e in exons:
                try:
                    seq.append(fa.fetch(chrom, s - 1, e).upper())
                except Exception:
                    seq = []
                    break
            if not seq:
                continue
            spliced = "".join(seq)
            if strand == "-":
                spliced = spliced.encode().translate(rc)[::-1].decode()
            if len(spliced) < 100:
                continue
            fh.write(f">{g}\n")
            for i in range(0, len(spliced), 80):
                fh.write(spliced[i:i + 80] + "\n")
            mh.write(f"{g}\t{chrom}\t{exons[0][0]}\t{exons[-1][1]}\t{strand}\t{len(spliced)}\n")
            n += 1
    print(f"wrote {n} gene reps -> {OUT_FA} (+ {OUT_META})")


if __name__ == "__main__":
    main()
