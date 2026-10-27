#!/usr/bin/env python3
"""Attribute polish losses of observed reference transcripts (missed.py class `polish`) to the polish step, via
polish_port.py on an `--assembly-polish none` run (gorilla KB3781; paths fixed, run in /mnt/linuxdisk/tmp/stgaps/miss)."""
import sys
sys.path.insert(0, "/mnt/linuxdisk/home/juanfraitu/rustle_stgaps/bench/stgaps")
import polish_port as P
from collections import Counter, defaultdict
data = P.load("ggo_nopolish.gtf")
exons, key, reads, gene, chain, span = data
sites = P.polish(data)
def cs(t): return (key[t][0], ",".join(f"{a}-{b-1}" for a, b in chain[t]))  # (exon end 1-based, next start-1) form
# port check: kept multi-exon chains vs base.gtf
base = P.load("/mnt/linuxdisk/tmp/stgaps/fwl/ggo/base.gtf")
bc = {(base[1][t][0], ",".join(f"{a}-{b-1}" for a, b in base[4][t])) for t in base[4] if base[4][t]}
kept = {cs(t) for t in chain if chain[t] and t not in sites}
print("port kept multi chains", len(kept), "base", len(bc), "kept-base", len(kept - bc), "base-kept", len(bc - kept))
by = defaultdict(list)
for t in chain:
    if chain[t]:
        by[cs(t)].append(t)
rows = [l.rstrip("\n").split("\t") for l in open("ggo.tsv")][1:]
import re
ref = P.load("/mnt/linuxdisk/tmp/stgaps/fw/ggo_ref.gtf")
rk = {t: (ref[1][t][0], ",".join(f"{a}-{b-1}" for a, b in ref[4][t])) for t in ref[4] if ref[4][t]}
cnt = Counter()
out = open("ggo_polish_steps.tsv", "w")
for r in rows:
    if r[6] != "polish":
        continue
    ts = by.get(rk[r[0]], [])
    st = sorted({sites.get(t, "kept") for t in ts})
    cnt[(r[1], "+".join(st))] += 1
    out.write("\t".join(r[:6] + ["+".join(st), ";".join(f"{t}:{reads.get(t,0)}:{gene.get(t)}" for t in ts)]) + "\n")
for k, v in sorted(cnt.items()):
    print(k, v)
# frac losses: is the gene's top transcript overlapping the dropped one, or elsewhere (another copy in the same gene_id)?
gbest = {}
for t in gene:
    if t in sites and sites[t] in ("ism", "ism_mono", "floor", "shadow"):
        continue
    g = gene[t]
    if g not in gbest or reads.get(t, 0) > reads.get(gbest[g], 0):
        gbest[g] = t
c2 = Counter()
for r in rows:
    if r[6] != "polish":
        continue
    for t in by.get(rk[r[0]], []):
        if sites.get(t) != "frac":
            continue
        b = gbest[gene[t]]
        ov = key[b] == key[t] and span[b][0] < span[t][1] and span[b][1] > span[t][0]
        gs = [u for u in gene if gene[u] == gene[t]] if False else None
        c2[(r[1], "top overlaps" if ov else "top elsewhere", "ratio<0.005" if reads.get(t,0) < 0.005*reads.get(b,0) else "ratio 0.005-0.02")] += 1
for k, v in sorted(c2.items()):
    print(k, v)
