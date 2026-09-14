#!/usr/bin/env python3
"""Guided mode with minimal annotation, genome scale: keep a random fraction of the annotated genes as seeds, find the
hidden copies by gene-body homology (`guided_pipeline.gene_body_chains`: asm20 chains, identity >= 0.80, aligned >= 0.50
of min(query, extrapolated span); chains overlapping a seed dropped; leader clustering), give each candidate the seed's
exons projected through its chain, and run the guided catalog's construction (`mcl_families` on the node graph) on
seeds + candidates.

prep:  guided_min.py prep --gff genes.gff --genome genome.fa --contigs c1,c2 --frac 0.5 --seed 1 --outdir DIR
align: guided_min.py align --outdir DIR --chunk K [--threads 4]          (seed spans -> contigs, asm20)
nodes: guided_min.py nodes --outdir DIR                                   (writes DIR/min.nodes.tsv)
graph: guided_min.py graph --outdir DIR --genes-paf allgenes.paf [--threads 4] --mcl-families BIN [mcl args]
  The node graph reuses the annotation all-vs-all PAF for seed-seed pairs and aligns candidate spans against all node
  spans; `mcl_families` then clusters exactly as the guided catalog does. Writes DIR/min.copies.tsv.
"""
import argparse
import collections
import glob
import os
import random
import re
import subprocess
import sys

import pysam

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


def cmd_prep(a):
    os.makedirs(a.outdir, exist_ok=True)
    contigs = a.contigs.split(",")
    genes, exons = load_genes(a.gff, set(contigs))
    names = sorted(genes)
    rng = random.Random(a.seed)
    seeds = set(rng.sample(names, round(a.frac * len(names))))
    genome = pysam.FastaFile(a.genome)
    with open(f"{a.outdir}/seeds.tsv", "w") as fh:
        fh.write("name\tchrom\tstart\tend\tstrand\tseed\texons\n")
        for n in names:
            c, s, e, st = genes[n]
            fh.write(f"{n}\t{c}\t{s}\t{e}\t{st}\t{int(n in seeds)}\t{','.join(f'{x}-{y}' for x, y in exons[n])}\n")
    with open(f"{a.outdir}/target.fa", "w") as fh:
        for c in contigs:
            fh.write(f">{c}\n{genome.fetch(c)}\n")
    order = sorted((n for n in names if n in seeds), key=lambda n: genes[n][2] - genes[n][1])
    chunks, cur, bp = [], [], 0
    for n in order:
        L = genes[n][2] - genes[n][1]
        if cur and bp + L > a.chunk_bp:
            chunks.append(cur)
            cur, bp = [], 0
        cur.append(n)
        bp += L
    if cur:
        chunks.append(cur)
    for i, ch in enumerate(chunks):
        with open(f"{a.outdir}/s{i:03d}.fa", "w") as fh:
            for n in ch:
                c, s, e, _ = genes[n]
                fh.write(f">{n}\n{genome.fetch(c, s, e).upper()}\n")
    print(f"{len(names)} genes, {len(seeds)} seeds ({a.frac}), {len(chunks)} seed chunks")


def cmd_align(a):
    q, out = f"{a.outdir}/s{a.chunk:03d}.fa", f"{a.outdir}/s{a.chunk:03d}.paf"
    if os.path.exists(out):
        return
    mmi = f"{a.outdir}/target.asm20.mmi"
    if not os.path.exists(mmi):
        subprocess.run(["minimap2", "-x", "asm20", "-t", str(a.threads), "-d", mmi, f"{a.outdir}/target.fa"], check=True,
                       stderr=subprocess.DEVNULL)
    with open(out + ".tmp", "w") as fh:
        subprocess.run(["minimap2", "-c", "-x", "asm20", "-N", "50", "-p", "0.1", "-t", str(a.threads), mmi, q], stdout=fh,
                       stderr=subprocess.DEVNULL, check=True)
    os.replace(out + ".tmp", out)


def read_seeds(outdir):
    rows = {}
    for line in open(f"{outdir}/seeds.tsv"):
        if line.startswith("name\t"):
            continue
        n, c, s, e, st, seed, ex = line.rstrip("\n").split("\t")
        rows[n] = {"chrom": c, "start0": int(s), "end": int(e), "strand": st, "seed": seed == "1",
                   "exons": [tuple(map(int, b.split("-"))) for b in ex.split(",")]}
    return rows


def candidates_from(chains, query_nodes, occupied):
    """Candidate loci from chains of `query_nodes` (name -> dict with chrom/start0/end/strand/exons) that overlap
    no `occupied` interval; leader clustering; exons projected through the chain."""
    free = [c for c in chains if c["q"] in query_nodes and not any(s < c["e"] and c["s"] < e for s, e in occupied[c["chrom"]])]
    out = []
    for c in gp.leaders(free):
        r = query_nodes[c["q"]]
        qex = [(s - r["start0"], e - r["start0"]) for s, e in r["exons"]]
        ex = gp.project_exons(c, qex)
        ex = [(max(s, c["xs"]), min(e, c["xe"])) for s, e in ex if min(e, c["xe"]) > max(s, c["xs"])]
        if not ex:
            continue
        strand = r["strand"] if c["strand"] == "+" else {"+": "-", "-": "+"}.get(r["strand"], ".")
        out.append((c["chrom"], min(ex[0][0], c["xs"]), max(ex[-1][1], c["xe"]), strand, ex, "candidate", c["q"]))
    return out, len(free)


def cmd_nodes(a):
    gp.MIN_ID, gp.MIN_COV = a.min_id, a.min_cov
    rec = read_seeds(a.outdir)
    seeds = {n: r for n, r in rec.items() if r["seed"]}
    chains = [c for f in sorted(glob.glob(f"{a.outdir}/s*.paf")) for c in gp.gene_body_chains(f)]
    occupied = collections.defaultdict(list)
    for r in seeds.values():
        occupied[r["chrom"]].append((r["start0"], r["end"]))
    nodes = [(r["chrom"], r["start0"], r["end"], r["strand"], r["exons"], "seed", n) for n, r in sorted(seeds.items())]
    new, n_free = candidates_from(chains, seeds, occupied)
    print(f"round 1: chains {len(chains)}, free {n_free}, candidates {len(new)}")
    rnd = 1
    while new and rnd < a.rounds:
        nodes.extend(new)
        for c, s0, e0, *_ in new:
            occupied[c].append((s0, e0))
        rnd += 1
        qn = {f"cand{rnd}_{i}": {"chrom": c, "start0": s0, "end": e0, "strand": st, "exons": ex}
              for i, (c, s0, e0, st, ex, _k, _src) in enumerate(new)}
        genome = pysam.FastaFile(f"{a.outdir}/target.fa")
        qfa = f"{a.outdir}/round{rnd}.fa"
        with open(qfa, "w") as fh:
            for k, r in qn.items():
                fh.write(f">{k}\n{genome.fetch(r['chrom'], r['start0'], r['end']).upper()}\n")
        paf = f"{a.outdir}/round{rnd}.paf"
        with open(paf, "w") as fh:
            subprocess.run(["minimap2", "-c", "-x", "asm20", "-N", "50", "-p", "0.1", "-t", str(a.threads),
                            f"{a.outdir}/target.asm20.mmi", qfa], stdout=fh, stderr=subprocess.DEVNULL, check=True)
        new, n_free = candidates_from(gp.gene_body_chains(paf), qn, occupied)
        print(f"round {rnd}: free chains {n_free}, new candidates {len(new)}")
    if new and rnd >= a.rounds:
        nodes.extend(new)
    n_cand = sum(1 for x in nodes if x[5] == "candidate")
    nodes.sort(key=lambda x: (x[0], x[1], x[2]))
    with open(f"{a.outdir}/{a.tag}.nodes.tsv", "w") as fh:
        fh.write("idx\tchrom\tstart\tend\tstrand\tn_exon\tn_reads\texons\tkind\tsource\n")
        for i, (c, s, e, st, ex, kind, src) in enumerate(nodes):
            fh.write(f"{i}\t{c}\t{s}\t{e}\t{st}\t{len(ex)}\t0\t{','.join(f'{x}-{y}' for x, y in ex)}\t{kind}\t{src}\n")
    print(f"seeds {len(seeds)}; candidate loci {n_cand}; nodes {len(nodes)}")


def cmd_graph(a, extra):
    import csv
    nodes = list(csv.DictReader(open(f"{a.outdir}/{a.tag}.nodes.tsv"), delimiter="\t"))
    genome_fa = f"{a.outdir}/target.fa"
    genome = pysam.FastaFile(genome_fa) if os.path.exists(genome_fa + ".fai") else None
    if genome is None:
        pysam.faidx(genome_fa)
        genome = pysam.FastaFile(genome_fa)
    key = lambda r: f"{r['chrom']}:{int(r['start']) + 1}-{r['end']}"
    with open(f"{a.outdir}/nodes.gff", "w") as g:
        g.write("##gff-version 3\n")
        for r in nodes:
            name = f"N{r['idx']}"
            g.write(f"{r['chrom']}\tmin\tgene\t{int(r['start']) + 1}\t{r['end']}\t.\t{r['strand'] if r['strand'] in '+-' else '.'}\t.\tID=gene-{name};Name={name}\n")
            for b in r["exons"].split(","):
                x, y = map(int, b.split("-"))
                g.write(f"{r['chrom']}\tmin\texon\t{x + 1}\t{y}\t.\t.\t.\tgene={name}\n")
    seed_keys = {key(r) for r in nodes if r["kind"] == "seed"}
    with open(f"{a.outdir}/allspans.fa", "w") as fh:
        for r in nodes:
            fh.write(f">{key(r)}\n{genome.fetch(r['chrom'], int(r['start']), int(r['end'])).upper()}\n")
    with open(f"{a.outdir}/candspans.fa", "w") as fh:
        for r in nodes:
            if r["kind"] == "candidate":
                fh.write(f">{key(r)}\n{genome.fetch(r['chrom'], int(r['start']), int(r['end'])).upper()}\n")
    paf = f"{a.outdir}/{a.tag}.graph.paf"
    if not os.path.exists(paf):
        with open(paf + ".tmp", "w") as out:
            for line in open(a.genes_paf):
                f = line.split("\t", 6)
                if f[0] in seed_keys and f[5] in seed_keys:
                    out.write(line)
            res = subprocess.run(["minimap2", "-x", "asm20", "-c", "-N", "50", "-p", "0.1", "-t", str(a.threads),
                                  f"{a.outdir}/allspans.fa", f"{a.outdir}/candspans.fa"], capture_output=True, text=True, check=True)
            out.write(res.stdout)
        os.replace(paf + ".tmp", paf)
    outp = f"{a.outdir}/{a.tag}"
    subprocess.run([a.mcl_families, "--paf", paf, "--gff", f"{a.outdir}/nodes.gff", "--dump-graph", f"{outp}.graph.tsv", "--out", outp, *extra],
                   check=True, stdout=subprocess.DEVNULL, stderr=open(f"{outp}.err", "w"))
    rows = list(csv.DictReader(open(f"{outp}.clusters.tsv"), delimiter="\t"))
    with open(f"{outp}.copies.tsv", "w") as fh:
        fh.write("family_id\tchrom\tstart\tend\n")
        for r in rows:
            fh.write(f"{r['cluster_id']}\t{r['chrom']}\t{int(r['start']) - 1}\t{r['end']}\n")
    print(f"{outp}: {len({r['cluster_id'] for r in rows})} clusters, {len(rows)} members")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("prep")
    for k in ("--gff", "--genome", "--contigs", "--outdir"):
        p.add_argument(k, required=True)
    p.add_argument("--frac", type=float, default=0.5)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--chunk-bp", type=int, default=4_000_000)
    p = sub.add_parser("align")
    p.add_argument("--outdir", required=True)
    p.add_argument("--chunk", type=int, required=True)
    p.add_argument("--threads", type=int, default=4)
    p = sub.add_parser("nodes")
    p.add_argument("--outdir", required=True)
    p.add_argument("--min-id", type=float, default=0.80)
    p.add_argument("--min-cov", type=float, default=0.50)
    p.add_argument("--rounds", type=int, default=1)
    p.add_argument("--tag", default="min")
    p.add_argument("--threads", type=int, default=4)
    p = sub.add_parser("graph")
    p.add_argument("--outdir", required=True)
    p.add_argument("--genes-paf", required=True)
    p.add_argument("--tag", default="min")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--mcl-families", required=True)
    a, extra = ap.parse_known_args()
    if a.cmd == "graph":
        cmd_graph(a, extra)
    else:
        {"prep": cmd_prep, "align": cmd_align, "nodes": cmd_nodes}[a.cmd](a)


if __name__ == "__main__":
    main()
