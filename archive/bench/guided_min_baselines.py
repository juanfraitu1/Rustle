#!/usr/bin/env python3
"""Comparators for `guided_min.py` (Addendum AH): all annotated genes and the seed genes alone, each through the guided
catalog's construction on the annotation PAF restricted to the substrate.
usage: guided_min_baselines.py --outdir DIR --genes-paf PAF --gff GFF --mcl-families BIN [mcl args]"""
import argparse, csv, os, subprocess
ap = argparse.ArgumentParser()
for k in ("--outdir", "--genes-paf", "--gff", "--mcl-families"):
    ap.add_argument(k, required=True)
a, extra = ap.parse_known_args()
seeds = {}
for line in open(f"{a.outdir}/seeds.tsv"):
    if line.startswith("name\t"):
        continue
    n, c, s, e, st, sd, ex = line.rstrip("\n").split("\t")
    if sd == "1":
        seeds[f"{c}:{int(s) + 1}-{e}"] = (c, s, e, st, ex, n)
contigs = {v[0] for v in seeds.values()}
with open(f"{a.outdir}/base_all.paf", "w") as fa, open(f"{a.outdir}/base_seeds.paf", "w") as fs:
    for line in open(a.genes_paf):
        f = line.split("\t", 6)
        if f[0].rsplit(":", 1)[0] in contigs and f[5].rsplit(":", 1)[0] in contigs:
            fa.write(line)
            if f[0] in seeds and f[5] in seeds:
                fs.write(line)
with open(f"{a.outdir}/base_seeds.gff", "w") as g:
    for k, (c, s, e, st, ex, n) in seeds.items():
        g.write(f"{c}\tseed\tgene\t{int(s) + 1}\t{e}\t.\t{st}\t.\tID=gene-{n};Name={n}\n")
        for b in ex.split(","):
            x, y = map(int, b.split("-"))
            g.write(f"{c}\tseed\texon\t{x + 1}\t{y}\t.\t.\t.\tgene={n}\n")
for tag, paf, gff in (("base_all", f"{a.outdir}/base_all.paf", a.gff), ("base_seeds", f"{a.outdir}/base_seeds.paf", f"{a.outdir}/base_seeds.gff")):
    out = f"{a.outdir}/{tag}"
    subprocess.run([a.mcl_families, "--paf", paf, "--gff", gff, "--out", out, *extra], check=True,
                   stdout=subprocess.DEVNULL, stderr=open(f"{out}.err", "w"))
    rows = list(csv.DictReader(open(f"{out}.clusters.tsv"), delimiter="\t"))
    with open(f"{out}.copies.tsv", "w") as fh:
        fh.write("family_id\tchrom\tstart\tend\n")
        for r in rows:
            fh.write(f"{r['cluster_id']}\t{r['chrom']}\t{int(r['start']) - 1}\t{r['end']}\n")
    print(tag, len(rows))
