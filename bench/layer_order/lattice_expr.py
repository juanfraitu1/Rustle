#!/usr/bin/env python3
"""Nested edge-test lattice, step 2: testis read counts for every node of V (lattice/nodes.tsv).

Rule = bench/layer_order/lo_expr_recount.py (itself = heavy/scripts/expr_counts.py): human_testis.t2t.bam, primary reads
only (samtools -F 2308), read blocks split at N and D, a read counts for gene G if >= 1 block overlaps >= 1 bp of an exon
of G (strand ignored); 'unique' = the read's blocks hit exons of exactly one RefSeq record genome-wide. Exon-less records
count on their gene body. Only the gene set differs (all of V instead of U). Check: the 68 U genes must reproduce
integrate_slim/expr_recount.tsv (any and unique).

usage: lattice_expr.py   -> lattice/expr_counts.tsv, lattice/expr_counts.out
"""
import bisect
import csv
import gzip
import re
import subprocess
import sys
import time
from collections import defaultdict

sys.path.insert(0, "/mnt/c/Users/jfris/Desktop/Rustle/bench/layer_order")
from lattice_common import INT, OUT, tsv  # noqa: E402

BAM = "/mnt/linuxdisk/home/juanfraitu/_from_wsl/human_val/human_testis.t2t.bam"
GFF = "/mnt/linuxdisk/home/juanfraitu/winloci_data/Reference/chm13v2.0_RefSeq_full.gff.gz"
CIG = re.compile(r"(\d+)([MIDNSHP=X])")
T0 = time.time()
LOG = []


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    LOG.append(s)


want = {r["gene_id"] for r in tsv(f"{OUT}/nodes.tsv")}
genes, parent, exraw = {}, {}, []
with gzip.open(GFF, "rt") as fh:
    for ln in fh:
        if ln[0] == "#":
            continue
        f = ln.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        t = f[2]
        if t not in ("gene", "pseudogene", "exon") and "Parent=" not in f[8]:
            continue
        a = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
        fid, par = a.get("ID"), a.get("Parent")
        if t in ("gene", "pseudogene"):
            s, e = int(f[3]) - 1, int(f[4])
            if fid in genes:
                c0, s0, e0 = genes[fid]
                if c0 == f[0]:
                    genes[fid] = (c0, min(s0, s), max(e0, e))
            else:
                genes[fid] = (f[0], s, e)
            continue
        if t == "exon" and par:
            exraw.append((par.split(",")[0], f[0], int(f[3]) - 1, int(f[4])))
        if fid and par:
            parent.setdefault(fid, par.split(",")[0])


def up(p):
    k = 0
    while p not in genes:
        p = parent.get(p)
        k += 1
        if p is None or k > 10:
            return None
    return p


exons = defaultdict(set)
for p, c, s, e in exraw:
    g = up(p)
    if g:
        exons[g].add((c, s, e))
for g, (c, s, e) in genes.items():
    if not exons.get(g):
        exons[g] = {(c, s, e)}
say(f"[gff] gene/pseudogene ids {len(genes)}; V genes {len(want)} (missing from GFF {len(want - set(genes))}); "
    f"{time.time() - T0:.0f}s")
per = defaultdict(list)
for g, lst in exons.items():
    for c, s, e in lst:
        per[c].append((s, e, g))
idx = {}
for c, lst in per.items():
    lst.sort()
    idx[c] = ([x[0] for x in lst], lst, max(x[1] - x[0] for x in lst))


def hits(c, bs, be):
    if c not in idx:
        return set()
    st, lst, ml = idx[c]
    i, j = bisect.bisect_left(st, bs - ml), bisect.bisect_left(st, be)
    return {lst[k][2] for k in range(i, j) if lst[k][1] > bs and lst[k][0] < be}


spans = defaultdict(list)
for g in want:
    c, s, e = genes[g]
    spans[c].append((s, e))
bed = f"{OUT}/expr_counts.windows.bed"
nwin = nbp = 0
with open(bed, "w") as out:
    for c in sorted(spans):
        iv = sorted(spans[c])
        cs, ce = iv[0]
        for s, e in iv[1:] + [(None, None)]:
            if s is not None and s <= ce:
                ce = max(ce, e)
                continue
            out.write(f"{c}\t{cs}\t{ce}\n")
            nwin += 1
            nbp += ce - cs
            if s is not None:
                cs, ce = s, e
p = subprocess.Popen(["samtools", "view", "-M", "-F", "2308", "-L", bed, BAM], stdout=subprocess.PIPE, text=True)
n_any, n_uni = defaultdict(int), defaultdict(int)
seen = set()
n_rec = n_hit = 0
for ln in p.stdout:
    f = ln.split("\t", 6)
    key = (f[0], f[1], f[2], f[3])
    if key in seen:
        continue
    seen.add(key)
    n_rec += 1
    c, pos = f[2], int(f[3]) - 1
    hit = set()
    for n, op in CIG.findall(f[5]):
        n = int(n)
        if op in "M=X":
            hit |= hits(c, pos, pos + n)
            pos += n
        elif op in "DN":
            pos += n
    if not hit:
        continue
    n_hit += 1
    for g in hit & want:
        n_any[g] += 1
    if len(hit) == 1:
        (g,) = hit
        if g in want:
            n_uni[g] += 1
assert p.wait() == 0
say(f"[bam] windows {nwin} ({nbp} bp); primary records {n_rec}; on >= 1 exon {n_hit}; {time.time() - T0:.0f}s")
prev = {"gene-" + r["gene_id"]: r for r in tsv(f"{INT}/expr_recount.tsv")}
diff = [g for g in prev if g in want and (n_any[g], n_uni[g]) != (int(prev[g]["n_reads_any"]), int(prev[g]["n_reads_unique"]))]
say(f"[check] genes also in integrate_slim/expr_recount.tsv: {sum(1 for g in prev if g in want)}; differing any/unique: "
    f"{len(diff)} {diff[:10]}")
with open(f"{OUT}/expr_counts.tsv", "w") as out:
    out.write("gene_id\tn_reads_any\tn_reads_unique\n")
    for g in sorted(want):
        out.write(f"{g}\t{n_any[g]}\t{n_uni[g]}\n")
for t in (1, 3):
    say(f"[summary] V genes with any-overlap reads >= {t}: {sum(1 for g in want if n_any[g] >= t)}")
with open(f"{OUT}/expr_counts.out", "w") as fh:
    fh.write("\n".join(LOG) + "\n")
