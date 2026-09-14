#!/usr/bin/env python3
"""Run the guided catalog's own construction (all-vs-all asm20 of loci spans -> mcl_families graph + MCL) on a node set
from `bench/denovo_shared_def.py` (read-derived nodes), so de novo and guided differ only in where nodes come from.

prep:  node_graph_mcl.py prep --nodes ac.nodes.tsv --genome genome.fa --outdir DIR [--chunk-bp 8000000]
align: node_graph_mcl.py align --outdir DIR --chunk K [--threads 4]
mcl:   node_graph_mcl.py mcl --outdir DIR --mcl-families BIN [mcl_families args ...]  -> DIR/nodes_mcl.clusters.tsv + copies.tsv
"""
import argparse
import csv
import glob
import os
import subprocess

import pysam


def cmd_prep(a):
    os.makedirs(a.outdir, exist_ok=True)
    genome = pysam.FastaFile(a.genome)
    nodes = list(csv.DictReader(open(a.nodes), delimiter="\t"))
    with open(f"{a.outdir}/nodes.gff", "w") as g:
        g.write("##gff-version 3\n")
        for r in nodes:
            s, e = int(r["start"]), int(r["end"])
            name = f"SD{r['idx']}"
            g.write(f"{r['chrom']}\tsd\tgene\t{s + 1}\t{e}\t.\t{r['strand'] if r['strand'] in '+-' else '.'}\t.\tID=gene-{name};Name={name}\n")
            for b in r["exons"].split(","):
                x, y = map(int, b.split("-"))
                g.write(f"{r['chrom']}\tsd\texon\t{x + 1}\t{y}\t.\t.\t.\tgene={name}\n")
    order = sorted(nodes, key=lambda r: int(r["end"]) - int(r["start"]))
    with open(f"{a.outdir}/spans.fa", "w") as fh:
        for r in nodes:
            s, e = int(r["start"]), int(r["end"])
            fh.write(f">{r['chrom']}:{s + 1}-{e}\n{genome.fetch(r['chrom'], s, e).upper()}\n")
    chunks, cur, bp = [], [], 0
    for r in order:
        L = int(r["end"]) - int(r["start"])
        if cur and bp + L > a.chunk_bp:
            chunks.append(cur)
            cur, bp = [], 0
        cur.append(r)
        bp += L
    if cur:
        chunks.append(cur)
    for i, ch in enumerate(chunks):
        with open(f"{a.outdir}/q{i:03d}.fa", "w") as fh:
            for r in ch:
                s, e = int(r["start"]), int(r["end"])
                fh.write(f">{r['chrom']}:{s + 1}-{e}\n{genome.fetch(r['chrom'], s, e).upper()}\n")
    print(f"{len(nodes)} nodes -> nodes.gff, spans.fa, {len(chunks)} query chunks")


def cmd_align(a):
    q = f"{a.outdir}/q{a.chunk:03d}.fa"
    out = f"{a.outdir}/q{a.chunk:03d}.paf"
    if os.path.exists(out):
        return
    mmi = f"{a.outdir}/spans.mmi"
    if not os.path.exists(mmi):
        subprocess.run(["minimap2", "-x", "asm20", "-t", str(a.threads), "-d", mmi, f"{a.outdir}/spans.fa"], check=True,
                       stderr=subprocess.DEVNULL)
    with open(out + ".tmp", "w") as fh:
        subprocess.run(["minimap2", "-x", "asm20", "-c", "-X", "-N", "50", "-p", "0.1", "-t", str(a.threads), mmi, q],
                       stdout=fh, stderr=subprocess.DEVNULL, check=True)
    os.replace(out + ".tmp", out)
    print(f"wrote {out}")


def cmd_mcl(a, extra):
    paf = f"{a.outdir}/all.paf"
    with open(paf, "w") as fh:
        for f in sorted(glob.glob(f"{a.outdir}/q*.paf")):
            fh.write(open(f).read())
    out = f"{a.outdir}/{a.tag}"
    subprocess.run([a.mcl_families, "--paf", paf, "--gff", f"{a.outdir}/nodes.gff", "--dump-graph", f"{out}.graph.tsv",
                    "--out", out, *extra], check=True, stdout=subprocess.DEVNULL, stderr=open(f"{out}.err", "w"))
    rows = list(csv.DictReader(open(f"{out}.clusters.tsv"), delimiter="\t"))
    with open(f"{out}.copies.tsv", "w") as fh:
        fh.write("family_id\tchrom\tstart\tend\n")
        for r in rows:
            fh.write(f"{r['cluster_id']}\t{r['chrom']}\t{int(r['start']) - 1}\t{r['end']}\n")
    print(f"{out}: {len({r['cluster_id'] for r in rows})} clusters, {len(rows)} members")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("prep")
    p.add_argument("--nodes", required=True)
    p.add_argument("--genome", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--chunk-bp", type=int, default=8_000_000)
    p = sub.add_parser("align")
    p.add_argument("--outdir", required=True)
    p.add_argument("--chunk", type=int, required=True)
    p.add_argument("--threads", type=int, default=4)
    p = sub.add_parser("mcl")
    p.add_argument("--outdir", required=True)
    p.add_argument("--mcl-families", required=True)
    p.add_argument("--tag", default="nodes_mcl")
    a, extra = ap.parse_known_args()
    if a.cmd == "mcl":
        cmd_mcl(a, extra)
    else:
        {"prep": cmd_prep, "align": cmd_align}[a.cmd](a)


if __name__ == "__main__":
    main()
