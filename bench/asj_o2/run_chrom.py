#!/usr/bin/env python3
"""Per-chromosome driver for docs/PREREG_asj_o2_2026-10-10.md (resumable; exit 75 = run again; prints DONE).

    run_chrom.py <primary chromosome>          e.g. NC_073224.2; chrX (no other haplotype) runs arm R only, genes >= 5 Mb from either end

Writes OUT/<chrom>.genes.tsv (expressed multi-exon genes), OUT/<chrom>.paf (gene spans on the other haplotype), OUT/<chrom>.jsonl (one gene per line)."""
import json
import os
import subprocess
import sys
import time

import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import asj_o2 as A  # noqa: E402

OUT = os.environ.get("ASJ_OUT", "/mnt/linuxdisk/tmp/asj_o2")
BAM = "/mnt/linuxdisk/home/juanfraitu/fibroblasts/GCA_029281585.2_flnc_mm.bam"
GFF = "/mnt/linuxdisk/home/juanfraitu/winloci_data/GGO_genomic.gff"
PRIMARY_FA = "/mnt/linuxdisk/tmp/rustle_figures/liftoff/gorilla/genome.fa"
HAP_FA = "/mnt/linuxdisk/home/juanfraitu/gorilla_haps/{}.fa"
CHRMAP = "/mnt/linuxdisk/tmp/rna_allele/chrmap.tsv"
FLANK, MIN_READS, X_END = 1000, 20, 5_000_000
BUDGET = int(os.environ.get("ASJ_BUDGET", "540"))
T0 = time.time()


def primaries(bam, c, s, e):
    for r in bam.fetch(c, s, e):
        if not (r.is_secondary or r.is_supplementary or r.is_unmapped) and r.query_sequence:
            yield r


def genes_of(chrom, bam, clen):
    p = f"{OUT}/{chrom}.genes.tsv"
    if os.path.exists(p):
        return [tuple(x if i in (0, 3, 4, 5) else int(x) for i, x in enumerate(ln.rstrip("\n").split("\t"))) for ln in open(p)]
    tp = f"{OUT}/genes.tsv"
    if not os.path.exists(tp):
        with open(GFF) as fh, open(tp + ".tmp", "w") as o:
            for g in A.gene_table(fh):
                o.write("\t".join(map(str, g)) + "\n")
        os.replace(tp + ".tmp", tp)
    keep = []
    for ln in open(tp):
        f = ln.rstrip("\n").split("\t")
        g = (f[0], int(f[1]), int(f[2]), f[3], f[4], f[5], int(f[6]))
        if g[0] != chrom or g[6] < 2:
            continue
        if chrom_is_x(chrom) and (g[1] < X_END or g[2] > clen - X_END):
            continue
        n = 0
        for _ in primaries(bam, chrom, g[1], g[2]):
            n += 1
            if n >= MIN_READS:
                keep.append(g)
                break
    with open(p, "w") as o:
        for g in keep:
            o.write("\t".join(map(str, g)) + "\n")
    return keep


def chrom_is_x(chrom):
    return CHR[chrom]["chrom"] == "X"


def het_of(chrom, genes, clen):
    """{gene id: (het, snvs, runs, indels)} from the gene spans' alignment to the other haplotype's chromosome"""
    b = CHR[chrom]["B_name"]
    if not b:
        return {}
    paf = f"{OUT}/{chrom}.paf"
    if not os.path.exists(paf + ".done"):
        pri = pysam.FastaFile(PRIMARY_FA)
        with open(f"{OUT}/{chrom}.q.fa", "w") as o:
            for g in genes:
                s, e = max(0, g[1] - FLANK), min(clen, g[2] + FLANK)
                o.write(f">{g[3]}|{s}\n{pri.fetch(chrom, s, e)}\n")
        subprocess.run(f"samtools faidx {HAP_FA.format(CHR[chrom]['B_hap'])} {b} > {OUT}/{chrom}.B.fa", shell=True, check=True)
        subprocess.run(f"minimap2 -x asm5 -c --cs -t 4 {OUT}/{chrom}.B.fa {OUT}/{chrom}.q.fa > {paf}", shell=True, check=True, stderr=subprocess.DEVNULL)
        open(paf + ".done", "w").write("ok")
        for x in (f"{OUT}/{chrom}.q.fa", f"{OUT}/{chrom}.B.fa"):
            os.remove(x)
    out = {}
    for q, f in A.allelic_alignments(open(paf)).items():
        gid, off = q.rsplit("|", 1)
        cs = next(t[5:] for t in f[12:] if t.startswith("cs:Z:"))
        snvs, runs, indels = A.cs_walk(cs, f[4], int(f[2]), int(f[3]), int(off))
        out[gid] = (A.het_sites(snvs, indels), snvs, runs, indels)
    return out


def main(chrom):
    os.makedirs(OUT, exist_ok=True)
    bam = pysam.AlignmentFile(BAM, "rb")
    clen = bam.get_reference_length(chrom)
    genes = genes_of(chrom, bam, clen)
    het = het_of(chrom, genes, clen)
    jp = f"{OUT}/{chrom}.jsonl"
    done = {json.loads(ln)["gene"] for ln in open(jp)} if os.path.exists(jp) else set()
    pri = pysam.FastaFile(PRIMARY_FA)
    with open(jp, "a") as o:
        for g in genes:
            if g[3] in done:
                continue
            if time.time() - T0 > BUDGET:
                print(f"{chrom}: paused at {len(done)} of {len(genes)} genes; run again")
                sys.exit(75)
            h, snvs, runs, indels = het.get(g[3], ({}, {}, [], []))
            reads = [(r.cigartuples, r.reference_start, r.query_sequence, r.mapping_quality) for r in primaries(bam, chrom, g[1], g[2])]
            lo = min([g[1]] + [r[1] for r in reads])
            hi = max([g[2]] + [r[1] + sum(n for op, n in r[0] if op in (0, 2, 3, 7, 8)) for r in reads])
            seq = pri.fetch(chrom, lo, hi).upper()
            res = A.process_gene(reads, h, snvs, runs, indels, lambda p: seq[p - lo] if lo <= p < hi else "N")
            res.update(gene=g[3], name=g[4], chrom=chrom, start=g[1], end=g[2], allelic_alignment=g[3] in het)
            o.write(json.dumps(res) + "\n")
            o.flush()
            done.add(g[3])
    print(f"{chrom}: {len(genes)} genes, {sum(1 for g in genes if g[3] in het)} with an allelic alignment; DONE")


CHR = {}
for _ln in list(open(CHRMAP))[1:]:
    _f = _ln.rstrip("\n").split("\t") + [""] * 7
    CHR[_f[0]] = dict(chrom=_f[1], B_hap=_f[5], B_name=_f[6])

if __name__ == "__main__":
    main(sys.argv[1])
