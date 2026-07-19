#!/usr/bin/env python3
"""DNA variation-graph ceiling vs Soto SDs (visual demo). Per flagship Soto family, build a base-level
variation graph (abpoa MSA -> GFA) from the member GENOMIC sequences and colour each copy green
(RNA-recovered) / red (DNA-only). The graph REPRESENTS all copies (the DNA ceiling); it does not
"detect" families. See docs/superpowers/specs/2026-07-18-soto-dna-vg-ceiling-design.md."""
import os
import sys
from collections import defaultdict

# ---- paths ----
BED = "bench/soto/80_fams.chr.bed"
MEMFA = "/mnt/linuxdisk/home/juanfraitu/winloci_data/soto_members.fa"
DETECT = "bench/soto/soto_member_detection.tsv"
OUT = "/home/juanfra/winloci_scratch/soto_vg"

# ---- colours (match copy_graph.rs) ----
GREEN = "#1e8e3e"   # RNA-recovered
RED = "#d93025"     # DNA-only (K=0 / silent / coverage)
GREY = "#9aa0a6"    # shared / conserved


def parse_family_members(bed_lines, family_id):
    """(gene, chrom, start, end) for every member of family_id, in file order. BED col4 = GENE|ID_k."""
    out = []
    for ln in bed_lines:
        f = ln.rstrip("\n").split("\t")
        if len(f) < 4:
            continue
        name = f[3]
        if "|" not in name:
            continue
        gene, fam = name.rsplit("|", 1)
        if fam == family_id:
            out.append((gene, f[0], int(f[1]), int(f[2])))
    return out


def read_fasta(path):
    seqs, cur, buf = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if cur is not None:
                seqs[cur] = "".join(buf)
            cur = line[1:].strip().split()[0]
            buf = []
        else:
            buf.append(line.strip())
    if cur is not None:
        seqs[cur] = "".join(buf)
    return seqs


def member_seq(fa, chrom, start, end):
    """Genomic sequence for a BED member. soto_members.fa headers are 1-based (start+1)."""
    return fa[f"{chrom}:{start + 1}-{end}"]
