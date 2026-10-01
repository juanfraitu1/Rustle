#!/usr/bin/env python3
"""Gene table + exon-sum FASTA for docs/PREREG_rna_allele_haplotype_count_2026-10-01.md (gene sets, step 2).

Every RefSeq gene/pseudogene record of `GGO_genomic.gff` on a `_pri` chromosome listed in chrmap.tsv with >= 1 exon: exon union over
its transcripts (0-based half-open). Plus the S_fam copies (gorilla NPIP/TBC1D3 copy truth of 2026-09-29, gene_id = cid, exons = union of
its truth transcripts). The exon-sum sequence (exon union concatenated, transcript strand) is written for the paralog and PAR tests.

    python3 genes.py --gff GGO_genomic.gff --chrmap chrmap.tsv --fam-gtf truth.ggo.gtf --fam-copies copies.ggo.tsv --fasta GGO.fasta \
        --out genes.tsv --fa exonsum.fa
"""
import argparse
import collections
import csv
import re
import subprocess

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def merge(iv):
    out = []
    for a, b in sorted(iv):
        if out and a <= out[-1][1]:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for k in ("gff", "chrmap", "fam_gtf", "fam_copies", "fasta", "out", "fa"):
        ap.add_argument("--" + k.replace("_", "-"), required=True)
    a = ap.parse_args(argv)
    keep = {r["pri"]: r["chrom"] for r in csv.DictReader(open(a.chrmap), delimiter="\t")}
    genes, tx_gene, ex = {}, {}, collections.defaultdict(list)
    for ln in open(a.gff):
        if ln[0] == "#":
            continue
        f = ln.rstrip("\n").split("\t")
        if len(f) < 9 or f[0] not in keep:
            continue
        at = dict(x.split("=", 1) for x in f[8].split(";") if "=" in x)
        i = at.get("ID", "")
        if f[2] in ("gene", "pseudogene"):
            genes[i] = dict(gene_id=i, name=at.get("Name", ""), biotype=at.get("gene_biotype", ""), chrom=f[0], strand=f[6],
                            desc=at.get("description", ""))
        elif f[2] == "exon":
            ex[at.get("Parent", "")].append((int(f[3]) - 1, int(f[4])))
        elif "Parent" in at and at["Parent"].startswith("gene-") and i:
            tx_gene[i] = at["Parent"]
    gex = collections.defaultdict(list)
    for p, v in ex.items():
        g = p if p in genes else tx_gene.get(p)
        if g in genes:
            gex[g].extend(v)
    rows = []
    for g, r in genes.items():
        if gex.get(g):
            rows.append(dict(r, exons=merge(gex[g]), set="refseq"))
    # S_fam copies: gene_id = cid, exons from the truth GTF
    fx = collections.defaultdict(list)
    for ln in open(a.fam_gtf):
        f = ln.rstrip("\n").split("\t")
        if len(f) > 8 and f[2] == "exon":
            fx[re.search(r'gene_id "([^"]+)"', f[8]).group(1)].append((int(f[3]) - 1, int(f[4])))
    for c in csv.DictReader(open(a.fam_copies), delimiter="\t"):
        if c["chrom"] in keep and fx.get(c["cid"]):
            rows.append(dict(gene_id=c["cid"], name=c["name"], biotype="S_fam:" + c["family"], chrom=c["chrom"], strand=c["strand"],
                             desc="", exons=merge(fx[c["cid"]]), set="S_fam"))
    # exon-sum sequences, one samtools call per chromosome region list
    by_chrom = collections.defaultdict(list)
    for r in rows:
        by_chrom[r["chrom"]].append(r)
    with open(a.out, "w") as out, open(a.fa, "w") as fa:
        out.write("gene_id\tname\tbiotype\tset\tchrom\tchrom_label\tstrand\texons\texonic_bp\n")
        for chrom, rr in by_chrom.items():
            regs = [f"{chrom}:{x + 1}-{y}" for r in rr for x, y in r["exons"]]
            seqs = {}
            for i in range(0, len(regs), 20000):
                res = subprocess.run(["samtools", "faidx", a.fasta] + regs[i:i + 20000], capture_output=True, text=True, check=True).stdout
                cur = None
                for line in res.splitlines():
                    if line.startswith(">"):
                        cur = line[1:]; seqs[cur] = []
                    else:
                        seqs[cur].append(line.strip())
            for r in rr:
                s = "".join("".join(seqs[f"{chrom}:{x + 1}-{y}"]) for x, y in r["exons"])
                if r["strand"] == "-":
                    s = s.translate(COMP)[::-1]
                bl = ",".join(f"{x}-{y}" for x, y in r["exons"])
                out.write(f"{r['gene_id']}\t{r['name']}\t{r['biotype']}\t{r['set']}\t{chrom}\t{keep[chrom]}\t{r['strand']}\t{bl}\t{len(s)}\n")
                fa.write(f">{r['gene_id']}\n")
                for j in range(0, len(s), 80):
                    fa.write(s[j:j + 80] + "\n")
    c = collections.Counter((r["set"], keep[r["chrom"]] == "X") for r in rows)
    print(len(rows), "records;", dict(c))


if __name__ == "__main__":
    main()
