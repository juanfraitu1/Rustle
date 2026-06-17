#!/usr/bin/env python3
"""Extract per-gene INTRON CHAIN (ordered exon-length vector) from RefSeq, genome-wide.
Representative transcript = longest (summed exon length), matching extract_gene_reps.py.
Minimizer-free structural signature for the intron-chain copy-discovery harness.

Output /tmp/gene_chains.tsv:  gene  chrom  strand  n_exon  exon_lens(comma)  intron_lens(comma)
Single-exon genes are emitted (n_exon=1, no introns) so the retrocopy blind spot is explicit.
"""
import glob
import os
import re
from collections import defaultdict

GW = "/tmp/gw"
OUT = "/tmp/gene_chains.tsv"


def attrs(field):
    d = {}
    for kv in field.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            d[k.strip()] = v.strip()
    return d


def main():
    gene_meta = {}
    feat_gene = {}
    tx_exons = defaultdict(list)
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
                    gene_meta.setdefault(g, (chrom, f[6]))
            elif f[2] in ("mRNA", "transcript", "primary_transcript"):
                if a.get("ID") and a.get("gene"):
                    feat_gene[(chrom, a["ID"])] = a["gene"]
            elif f[2] == "exon":
                p = a.get("Parent", "")
                if (chrom, p) in feat_gene:
                    tx_exons[(chrom, p)].append((int(f[3]), int(f[4])))

    gene_best = {}
    for (chrom, fid), exons in tx_exons.items():
        g = feat_gene[(chrom, fid)]
        tot = sum(e - s + 1 for s, e in exons)
        if g not in gene_best or tot > gene_best[g][1]:
            gene_best[g] = (sorted(exons), tot, chrom)

    n = 0
    with open(OUT, "w") as fh:
        fh.write("gene\tchrom\tstrand\tn_exon\texon_lens\tintron_lens\n")
        for g in sorted(gene_best):
            exons, _, chrom = gene_best[g]
            strand = gene_meta.get(g, (chrom, "+"))[1]
            # transcription order: + strand left->right; - strand right->left
            ordered = exons if strand == "+" else exons[::-1]
            exon_lens = [e - s + 1 for s, e in ordered]
            # intron lengths between consecutive exons (genomic order, abs)
            gx = exons
            intron_lens = [gx[i + 1][0] - gx[i][1] - 1 for i in range(len(gx) - 1)]
            if strand == "-":
                intron_lens = intron_lens[::-1]
            fh.write(f"{g}\t{chrom}\t{strand}\t{len(exon_lens)}\t"
                     f"{','.join(map(str, exon_lens))}\t{','.join(map(str, intron_lens))}\n")
            n += 1
    print(f"wrote {n} gene chains -> {OUT}")

    # sizing: exon-count distribution + RABL2 + multi-exon coverage
    counts = defaultdict(int)
    multi = 0
    rabl = {}
    for g in sorted(gene_best):
        exons, _, _ = gene_best[g]
        ne = len(exons)
        counts[ne] += 1
        if ne >= 2:
            multi += 1
        if "RABL2" in g:
            rabl[g] = ne
    print(f"multi-exon (>=2): {multi}/{n} ; single-exon (retrocopy blind spot): {n - multi}")
    print("exon-count buckets (n_exon: n_genes), top by size:")
    for ne, c in sorted(counts.items(), key=lambda kv: -kv[1])[:12]:
        print(f"   n_exon={ne}: {c} genes  (pairs in bucket ~ {c*(c-1)//2})")
    print("RABL2 exon counts:", rabl)
    biggest = max(counts.items(), key=lambda kv: kv[1])
    print(f"biggest bucket: n_exon={biggest[0]} with {biggest[1]} genes")


if __name__ == "__main__":
    main()
