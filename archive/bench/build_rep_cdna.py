#!/usr/bin/env python3
"""Pick the longest cDNA per annotated gene (the genes.bed vertex set) -> rep_cdna.fa,
header = gene symbol. This is the DNA representative for all-vs-all paralogy.
"""
import sys

GENES_BED = "/home/juanfra/winloci_scratch/unmapped_poc/genes.bed"
GFF = "/home/juanfra/winloci_scratch/GGO_tx.gff"
CDNA = "/tmp/test_cdna.fa"
OUT = "/home/juanfra/winloci_scratch/rep_cdna.fa"

# 1) vertex universe
genes = set()
with open(GENES_BED) as f:
    for line in f:
        genes.add(line.rstrip("\n").split("\t")[3])
print(f"vertex genes: {len(genes):,}", flush=True)

# 2) feature ID -> gene symbol  (ID=rna-XXX / gene-XXX  ;  gene=SYMBOL)
id2gene = {}
with open(GFF) as f:
    for line in f:
        if line.startswith("#"):
            continue
        c = line.split("\t")
        if len(c) < 9:
            continue
        attr = c[8]
        fid = None
        gene = None
        for kv in attr.rstrip().split(";"):
            if kv.startswith("ID="):
                fid = kv[3:]
            elif kv.startswith("gene="):
                gene = kv[5:]
        if fid and gene and gene in genes:
            id2gene[fid] = gene
print(f"feature IDs mapped to vertex genes: {len(id2gene):,}", flush=True)

# 3) stream cDNA, keep longest per gene
best = {}     # gene -> (length, header_id)
seqs = {}     # header_id -> sequence (only for current best; rebuild at end)
cur_id, cur_len, cur_seq = None, 0, []
records = {}  # header_id -> seq  (store all; 386MB fits in 17GB)

def flush(hid, seq_parts):
    if hid is None:
        return
    seq = "".join(seq_parts)
    gene = id2gene.get(hid)
    if gene is None:
        return
    L = len(seq)
    if gene not in best or L > best[gene][0]:
        best[gene] = (L, hid)
        records[gene] = seq

with open(CDNA) as f:
    for line in f:
        if line.startswith(">"):
            flush(cur_id, cur_seq)
            cur_id = line[1:].split()[0]
            cur_seq = []
        else:
            cur_seq.append(line.strip())
    flush(cur_id, cur_seq)

print(f"genes with a cDNA representative: {len(records):,}", flush=True)
with open(OUT, "w") as o:
    for gene, seq in records.items():
        o.write(f">{gene}\n")
        for i in range(0, len(seq), 70):
            o.write(seq[i:i+70] + "\n")
print(f"wrote {OUT}", flush=True)
