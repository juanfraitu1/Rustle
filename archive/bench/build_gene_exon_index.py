#!/usr/bin/env python3
"""build_gene_exon_index.py — per-gene representative-transcript EXON COORDINATES from the GFF.

The intron index (gene_intron_index.json) stores only intron/exon LENGTHS, not genomic positions,
so it cannot build spliced cDNA sequences for junction-position alignment. This one streaming pass
records, per gene NAME, the representative transcript (most exons; tie -> longest span) as
{chrom, strand, exons:[[start,end],...]} (1-based GFF coords). Used by the genome-wide junction-
POSITION co-threading metric (family_def_cothread_genomewide.py).

Run: /home/juanfra/miniforge3/bin/python bench/build_gene_exon_index.py
Out: /home/juanfra/winloci_scratch/gene_exon_index.json
"""
import collections
import json

GFF = "/mnt/c/Users/jfris/Desktop/GGO_genomic.gff"
OUT = "/home/juanfra/winloci_scratch/gene_exon_index.json"


def attr(s, key):
    for kv in s.split(";"):
        if kv.startswith(key + "="):
            return kv[len(key) + 1:]
    return None


def main():
    gene_name = {}                       # gene ID -> Name (or ID stripped)
    mrna_gene = {}                       # mRNA ID -> gene ID
    mrna_meta = {}                       # mRNA ID -> (chrom, strand)
    exons = collections.defaultdict(list)  # mRNA ID -> [(s,e)]
    n = 0
    with open(GFF) as f:
        for line in f:
            if line.startswith("#"):
                continue
            t = line.rstrip("\n").split("\t")
            if len(t) < 9:
                continue
            typ = t[2]
            if typ == "gene" or typ == "pseudogene":
                gid = attr(t[8], "ID")
                if gid:
                    gene_name[gid] = attr(t[8], "Name") or (gid[5:] if gid.startswith("gene-") else gid)
            elif typ in ("mRNA", "transcript", "lnc_RNA", "ncRNA", "primary_transcript"):
                mid = attr(t[8], "ID")
                if mid:
                    mrna_gene[mid] = attr(t[8], "Parent")
                    mrna_meta[mid] = (t[0], t[6])
            elif typ == "exon":
                par = attr(t[8], "Parent")
                if par:
                    exons[par].append((int(t[3]), int(t[4])))
            n += 1
    # pick representative transcript per gene
    best = {}   # gene ID -> (key, mid)
    for mid, exs in exons.items():
        gid = mrna_gene.get(mid)
        if gid is None:
            continue
        exs_s = sorted(exs)
        key = (len(exs_s), exs_s[-1][1] - exs_s[0][0])
        if gid not in best or key > best[gid][0]:
            best[gid] = (key, mid)
    out = {}
    for gid, (key, mid) in best.items():
        name = gene_name.get(gid, gid[5:] if gid.startswith("gene-") else gid)
        chrom, strand = mrna_meta.get(mid, (None, "+"))
        if chrom is None:
            continue
        out[name] = dict(chrom=chrom, strand=strand, exons=sorted(exons[mid]))
    json.dump(out, open(OUT, "w"))
    print(f"[+] parsed {n} GFF feature lines; {len(out)} genes with exon coords -> {OUT}")
    # spot-check a few panel genes
    for g in ("RABL2A", "RABL2B", "APOBEC3D", "RFPL1", "EEF1A1", "LOC109023808"):
        v = out.get(g)
        if v is None:
            print(f"    {g:16} MISSING")
        else:
            print(f"    {g:16} {v['chrom']} {v['strand']} {len(v['exons'])} exons")


if __name__ == "__main__":
    main()
