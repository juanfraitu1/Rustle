#!/usr/bin/env python3
"""Assemble every number the artifact renders into one JSON blob."""
import json, math, statistics, subprocess

SRC = "/mnt/linuxdisk/home/juanfraitu/ecband/edge_ec2.tsv"
COLS = ["reads_ov00", "reads_ov25", "reads_ov50", "reads_ov80"]

rows = []
with open(SRC) as fh:
    h = fh.readline().rstrip("\n").split("\t"); ci = {c: i for i, c in enumerate(h)}
    for l in fh:
        f = l.rstrip("\n").split("\t")
        rows.append((float(f[ci["identity"]]), float(f[ci["coverage"]]),
                     [int(f[ci[c]]) for c in COLS]))
n = len(rows)
PRIM = 0                                  # any-overlap: most generous to E_c

def wilson(k, m, z=1.96):
    if m == 0: return [0.0, 0.0]
    p = k/m; d = 1 + z*z/m
    c = (p + z*z/(2*m))/d
    hw = z*math.sqrt(p*(1-p)/m + z*z/(4*m*m))/d
    return [max(0.0, c-hw), min(1.0, c+hw)]

def band(lo, hi):
    sub = [r for r in rows if lo <= r[0] < hi]
    k = sum(1 for r in sub if r[2][PRIM] >= 1)
    d = [r[2][PRIM] for r in sub if r[2][PRIM] >= 1]
    return {"lo": lo, "hi": min(hi, 1.0), "n": len(sub), "k": k,
            "p": (k/len(sub) if sub else 0), "ci": wilson(k, len(sub)),
            "medDepth": (statistics.median(d) if d else 0),
            "maxDepth": (max(d) if d else 0)}

W, LO, NB = 0.02, 0.60, 20
buckets = [[] for _ in range(NB)]
for r in rows:
    b = min(NB-1, max(0, int((r[0]-LO)/W)))
    buckets[b].append(r)
assert sum(len(b) for b in buckets) == n
fine = []
for i, sub in enumerate(buckets):
    k = sum(1 for r in sub if r[2][PRIM] >= 1)
    fine.append({"lo": round(LO+i*W, 4), "hi": round(LO+(i+1)*W, 4), "n": len(sub),
                 "k": k, "p": (k/len(sub) if sub else None), "ci": wilson(k, len(sub))})

BANDS = ((0.60, 0.80), (0.80, 0.90), (0.90, 0.95), (0.95, 1.01))
coarse = [band(a, b) for a, b in BANDS]

def band_at(lo, hi, j):
    sub = [r for r in rows if lo <= r[0] < hi]
    k = sum(1 for r in sub if r[2][j] >= 1)
    return {"n": len(sub), "k": k, "p": (k/len(sub) if sub else 0)}
robust = {"rules": [0.0, 0.25, 0.50, 0.80],
          "bands": [{"lo": a, "hi": min(b, 1.0),
                     "cells": [band_at(a, b, j) for j in range(4)]} for a, b in BANDS]}
top    = [band(a, b) for a, b in ((0.95, 0.97), (0.97, 0.99), (0.99, 1.01))]

ids = sorted(r[0] for r in rows)
withs   = [r[0] for r in rows if r[2][PRIM] >= 1]
without = [r[0] for r in rows if r[2][PRIM] == 0]

# overlap-rule robustness, read back from the run log
sweep = []
for line in open("/mnt/linuxdisk/home/juanfraitu/ecband/run2.err"):
    if line.startswith("OV=") and "E_r supported" in line:
        ov = float(line.split("=")[1].split(":")[0])
        sup = int(line.split("supported ")[1].split("/")[0])
        extra = int(line.rsplit(" ", 1)[1])
        sweep.append({"ov": ov, "k": sup, "p": sup/n, "extra": extra})

theory = [{"p": round(0.60 + i*0.005, 4),
           "score": max(0.0, 5*(0.60+i*0.005) - 4),
           "seed": (0.60+i*0.005)**15} for i in range(81)]

out = {
  "n": n,
  "fine": fine, "coarse": coarse, "robust": robust, "top": top, "theory": theory, "sweep": sweep,
  "depth": [{"m": m, "k": sum(1 for r in rows if r[2][PRIM] >= m)}
            for m in (1, 2, 5, 10, 25, 50)],
  "shared": sum(1 for r in rows if r[2][PRIM] >= 1),
  "extra": sweep[0]["extra"] if sweep else None,
  "medianIdentity": statistics.median(ids),
  "medianWith": statistics.median(withs),
  "medianWithout": statistics.median(without),
  "topBandShare": sum(1 for r in rows if r[0] >= 0.95)/n,
  "cmds": (
    "# reads: the substrate's own alignment, unmodified\n"
    "minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes \\\n"
    "         GGO.fasta  GCA_029281585.2_flnc.fastq.gz\n\n"
    "# E_r: the shipped edge dump for this catalog\n"
    "arm_f2/dump/e.edges.tsv   3,141 edges over 3,598 loci\n"
    "  identity >= 0.60 and coverage >= 0.50 of the shorter, one record, forward only\n\n"
    "# E_c: recomputed here on the identical loci\n"
    "samtools view npip3.bam | ecband2.py     # tie = AS >= 0.95 * the read's best AS\n\n"
    "# scoring constants, measured not assumed\n"
    "minimap2 -ax splice:hq --eqx ref.fa reads_with_k_mismatches.fa\n"
    "  2000 matches, 0 mm -> AS 2000    |  1950 matches, 50 mm -> AS 1750\n"
    "  => A = 1 per match, B = 4 per mismatch\n"),
}
json.dump(out, open("/mnt/linuxdisk/home/juanfraitu/ecband/artifact.json", "w"))
def show(b): return f"{b['lo']:.2f}-{b['hi']:.2f} n={b['n']:5d} k={b['k']:4d} p={b['p']:.4f}"
print("coarse:"); [print("  ", show(b), f"medDepth={b['medDepth']:.0f} max={b['maxDepth']}") for b in coarse]
print("top:");    [print("  ", show(b)) for b in top]
print("sweep:");  [print("  ", s) for s in sweep]
print("depth:");  [print("  ", d, f"{d['k']/n:.4f}") for d in out["depth"]]
print("robust:")
for b in robust["bands"]:
    print(f"   {b['lo']:.2f}-{b['hi']:.2f} " + " ".join(f"{c['p']:7.4f}" for c in b["cells"]))
for k in ("shared","extra","medianIdentity","medianWith","medianWithout","topBandShare"):
    print(f"  {k} = {out[k]}")
