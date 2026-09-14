#!/usr/bin/env python3
"""Prereg Addendum AI/AJ: gene-level node tables from an annotation, for `bench/node_graph_mcl.py prep`.

usage: annotation_nodes.py {refseq|cat} GFF CONTIGS OUT_NODES_TSV
  refseq: gene+pseudogene records, exon union from exon lines' `gene=` (`guided_min.load_genes`).
  cat:    CAT/Liftoff GENCODE `gene` records, exon union from the exons of the gene's transcripts (exon -> Parent
          transcript -> Parent gene); a gene with no exon keeps its span.
Writes the node table (idx order = sorted by chrom, start) and `<OUT>.names.tsv` (idx, name, biotype).
"""
import collections
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import guided_min  # noqa: E402
import guided_pipeline as gp  # noqa: E402


def attrs(col):
    return dict(kv.split("=", 1) for kv in col.strip().split(";") if "=" in kv)


def refseq(gff, contigs):
    genes, exons = guided_min.load_genes(gff, contigs)
    biotype, n_records = {}, collections.Counter()
    for line in open(gff):
        f = line.rstrip("\n").split("\t")
        if len(f) > 8 and f[0] in contigs and f[2] in ("gene", "pseudogene"):
            a = attrs(f[8])
            if "Name" in a:
                n_records[a["Name"]] += 1
                biotype[a["Name"]] = a.get("gene_biotype", f[2])
    dup = sum(1 for v in n_records.values() if v > 1)
    print(f"refseq: {len(genes)} named genes ({dup} names on > 1 record; last record kept, as load_genes)")
    return [(c, s, e, st, exons[n], n, biotype.get(n, "?")) for n, (c, s, e, st) in genes.items()]


def cat(gff, contigs):
    genes, tx_gene, exons = {}, {}, collections.defaultdict(list)
    for line in open(gff):
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[0] not in contigs:
            continue
        a = attrs(f[8])
        if f[2] == "gene":
            genes[a["ID"]] = (f[0], int(f[3]) - 1, int(f[4]), f[6], a.get("Name", a["ID"]),
                              a.get("gene_biotype", a.get("gene_type", "?")))
        elif f[2] == "transcript":
            tx_gene[a["ID"]] = a.get("Parent")
        elif f[2] == "exon":
            exons[a.get("Parent")].append((int(f[3]) - 1, int(f[4])))
    gex = collections.defaultdict(list)
    for t, ex in exons.items():
        if tx_gene.get(t) in genes:
            gex[tx_gene[t]] += ex
    print(f"cat: {len(genes)} genes, {sum(1 for g in genes if gex.get(g))} with exons")
    return [(c, s, e, st, gp.merge(gex.get(g) or [(s, e)]), name, bt) for g, (c, s, e, st, name, bt) in genes.items()]


def main():
    kind, gff, contigs, out = sys.argv[1:5]
    rows = sorted({"refseq": refseq, "cat": cat}[kind](gff, set(contigs.split(","))))
    with open(out, "w") as fh, open(out + ".names.tsv", "w") as fn:
        fh.write("idx\tchrom\tstart\tend\tstrand\tn_exon\tn_reads\texons\n")
        fn.write("idx\tname\tbiotype\n")
        for i, (c, s, e, st, ex, name, bt) in enumerate(rows):
            fh.write(f"{i}\t{c}\t{s}\t{e}\t{st}\t{len(ex)}\t0\t{','.join(f'{x}-{y}' for x, y in ex)}\n")
            fn.write(f"{i}\t{name}\t{bt}\n")
    print(f"wrote {out}: {len(rows)} nodes")


if __name__ == "__main__":
    main()
