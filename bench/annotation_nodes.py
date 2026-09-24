#!/usr/bin/env python3
"""Prereg Addenda AI/AJ/AK: gene-level node tables from an annotation, for `bench/node_graph_mcl.py prep`.

usage: annotation_nodes.py {refseq|cat|ensembl} GFF CONTIGS OUT_NODES_TSV
  refseq:  gene+pseudogene records, exon union from exon lines' `gene=` (`guided_min.load_genes`).
  cat:     CAT/Liftoff GENCODE `gene` records, exon union from the exons of the gene's transcripts (exon -> Parent
           transcript -> Parent gene); a gene with no exon keeps its span.
  ensembl: Ensembl gene / ncRNA_gene / pseudogene records (contigs named N, renamed chrN), exon union from the exons of
           every feature whose Parent is the gene.
Writes the node table (idx order = sorted by chrom, start), `<OUT>.names.tsv` (idx, name, biotype) and `<OUT>.cds.tsv`
(idx, strand, CDS segments `start-end:phase` of the gene's transcript with the longest CDS; genes without CDS omitted).
"""
import collections
import re
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import guided_pipeline as gp  # noqa: E402


def load_genes(gff, contigs):
    genes, exons = {}, collections.defaultdict(list)
    for line in open(gff):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[0] not in contigs:
            continue
        if f[2] in ("gene", "pseudogene"):
            n = re.search(r"(?:^|;)Name=([^;]+)", f[8])
            if n:
                genes[n.group(1)] = (f[0], int(f[3]) - 1, int(f[4]), f[6])
        elif f[2] == "exon":
            g = re.search(r"(?:^|;)gene=([^;]+)", f[8])
            if g:
                exons[g.group(1)].append((int(f[3]) - 1, int(f[4])))
    return genes, {n: gp.merge(exons.get(n) or [(g[1], g[2])]) for n, g in genes.items()}


def attrs(col):
    return dict(kv.split("=", 1) for kv in col.strip().split(";") if "=" in kv)


def longest(cds_by_tx):
    """transcript -> [(start0, end, phase)] -> the longest CDS's segments (ties: first transcript id)."""
    best = None
    for t in sorted(cds_by_tx):
        segs = cds_by_tx[t]
        L = sum(e - s for s, e, _ in segs)
        if best is None or L > best[0]:
            best = (L, sorted(segs))
    return best[1] if best else None


def refseq(gff, contigs):
    genes, exons = load_genes(gff, contigs)
    biotype, n_records = {}, collections.Counter()
    cds = collections.defaultdict(lambda: collections.defaultdict(list))
    for line in open(gff):
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[0] not in contigs:
            continue
        if f[2] in ("gene", "pseudogene"):
            a = attrs(f[8])
            if "Name" in a:
                n_records[a["Name"]] += 1
                biotype[a["Name"]] = a.get("gene_biotype", f[2])
        elif f[2] == "CDS":
            a = attrs(f[8])
            if "gene" in a:
                cds[a["gene"]][a.get("Parent", "?")].append((int(f[3]) - 1, int(f[4]), int(f[7]) if f[7] in "012" else 0))
    dup = sum(1 for v in n_records.values() if v > 1)
    print(f"refseq: {len(genes)} named genes ({dup} names on > 1 record; last record kept, as load_genes)")
    return [(c, s, e, st, exons[n], n, biotype.get(n, "?"), longest(cds[n]) if n in cds else None)
            for n, (c, s, e, st) in genes.items()]


def gencode_like(gff, contigs, gene_types, rename=None):
    genes, parent, exons = {}, {}, collections.defaultdict(list)
    cds = collections.defaultdict(list)
    for line in open(gff):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        chrom = rename(f[0]) if rename else f[0]
        if chrom not in contigs:
            continue
        a = attrs(f[8])
        if f[2] in gene_types:
            genes[a["ID"]] = (chrom, int(f[3]) - 1, int(f[4]), f[6], a.get("Name", a.get("gene_id", a["ID"])),
                              a.get("gene_biotype", a.get("gene_type", a.get("biotype", "?"))))
        elif f[2] == "exon":
            exons[a.get("Parent")].append((int(f[3]) - 1, int(f[4])))
        elif f[2] == "CDS":
            cds[a.get("Parent")].append((int(f[3]) - 1, int(f[4]), int(f[7]) if f[7] in "012" else 0))
        elif "ID" in a and "Parent" in a:
            parent[a["ID"]] = a["Parent"]
    gex = collections.defaultdict(list)
    gcds = collections.defaultdict(dict)
    for t, ex in exons.items():
        if parent.get(t) in genes:
            gex[parent[t]] += ex
    for t, segs in cds.items():
        if parent.get(t) in genes:
            gcds[parent[t]][t] = segs
    print(f"{len(genes)} genes, {sum(1 for g in genes if gex.get(g))} with exons, {len(gcds)} with CDS")
    return [(c, s, e, st, gp.merge(gex.get(g) or [(s, e)]), name, bt, longest(gcds[g]) if g in gcds else None)
            for g, (c, s, e, st, name, bt) in genes.items()]


def cat(gff, contigs):
    return gencode_like(gff, contigs, ("gene",))


def ensembl(gff, contigs):
    return gencode_like(gff, contigs, ("gene", "ncRNA_gene", "pseudogene"), rename=lambda c: c if c.startswith("chr") else "chr" + c)


def main():
    kind, gff, contigs, out = sys.argv[1:5]
    rows = sorted({"refseq": refseq, "cat": cat, "ensembl": ensembl}[kind](gff, set(contigs.split(","))), key=lambda r: r[:7])
    with open(out, "w") as fh, open(out + ".names.tsv", "w") as fn, open(out + ".cds.tsv", "w") as fc:
        fh.write("idx\tchrom\tstart\tend\tstrand\tn_exon\tn_reads\texons\n")
        fn.write("idx\tname\tbiotype\n")
        fc.write("idx\tstrand\tcds\n")
        for i, (c, s, e, st, ex, name, bt, cd) in enumerate(rows):
            fh.write(f"{i}\t{c}\t{s}\t{e}\t{st}\t{len(ex)}\t0\t{','.join(f'{x}-{y}' for x, y in ex)}\n")
            fn.write(f"{i}\t{name}\t{bt}\n")
            if cd:
                fc.write(f"{i}\t{st}\t{','.join(f'{x}-{y}:{p}' for x, y, p in cd)}\n")
    print(f"wrote {out}: {len(rows)} nodes")


if __name__ == "__main__":
    main()
