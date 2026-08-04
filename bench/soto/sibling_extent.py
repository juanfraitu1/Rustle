#!/usr/bin/env python3
"""Can a locus boundary be fixed by keeping only the part its FAMILY SIBLINGS also have?

Usage: sibling_extent.py <copies.tsv> <bam> <genome.fa> [chr1,chr15]

THE PROBLEM. Two boundary estimators, both biased, measured on the same 114 detected truth members:
    representative span      median 0.55x of truth   (47% below 0.5x)
    union extent of reads    median 11.27x of truth  (95% reach >= 0.8x)
The reads DO cover the locus; the pipeline inherits the boundaries of one chosen chain. But read extent
overshoots badly because reads spill into neighbouring transcription.

THE HYPOTHESIS. The spill is into UNIQUE sequence. A locus is a copy, so the real locus is duplicated and
its siblings carry it, while whatever the reads ran into next is not duplicated and no sibling carries it.
So: take the generous read extent, align it against the read extents of its predicted family siblings, and
keep only the sub-interval covered by at least one sibling alignment. The family supplies the bound; no
annotation is consulted.

THE RISK THIS TESTS FOR. These loci sit inside segmental duplications, where homology extends far past the
gene -- a ±10 kb padding test earlier found homology continuing for 73% of pairs at the 5' side and 76% at
the 3' side, with the 90th percentile saturating the pad. If that dominates, sibling homology will NOT trim
and the output will look like the read extent. That outcome refutes the idea; it does not need rescuing.
"""
import os
import subprocess
import sys
from collections import defaultdict

COP, BAM, GENOME = sys.argv[1], sys.argv[2], sys.argv[3]
CHROMS = set(sys.argv[4].split(",")) if len(sys.argv) > 4 else None
BED = "bench/soto/80_fams.gene_preferred.bed"
CACHE = "/tmp/claude-1000/-mnt-c-Users-jfris-Desktop/161b66ca-e17e-4e70-b1cf-8b479d6c328d/scratchpad"
MIN_ID = 0.90

copies = []
for i, ln in enumerate(open(COP)):
    if i == 0:
        continue
    f = ln.rstrip("\n").split("\t")
    if len(f) < 7 or (CHROMS and f[3] not in CHROMS):
        continue
    copies.append({"fam": f[0], "chrom": f[3], "s": int(f[4]), "e": int(f[5]), "nexon": int(f[6])})

# ---- generous read extent per copy (the upper-biased estimator) -----------------------------------
for c in copies:
    out = subprocess.run(["samtools", "view", "-F", "2308", BAM, f"{c['chrom']}:{c['s']}-{c['e']}"],
                         capture_output=True, text=True).stdout
    lo, hi = c["s"], c["e"]
    starts, ends, spans = defaultdict(int), defaultdict(int), []
    for ln in out.splitlines():
        f = ln.split("\t")
        p = int(f[3])
        q, n = p, 0
        for ch in f[5]:
            if ch.isdigit():
                n = n * 10 + ord(ch) - 48
            else:
                if ch in "MDN=X":
                    q += n
                n = 0
        spans.append((p, q))
        starts[p // 100] += 1
        ends[q // 100] += 1
    keep = [(a, b) for (a, b) in spans if starts[a // 100] >= 3 or ends[b // 100] >= 3]
    if keep:
        lo, hi = min(a for a, _ in keep), max(b for _, b in keep)
    c["rlo"], c["rhi"] = lo, hi

fa = f"{CACHE}/sibext.fa"
name = {}
with open(f"{CACHE}/sibext.txt", "w") as fh:
    for i, c in enumerate(copies):
        fh.write(f"{c['chrom']}:{c['rlo']+1}-{c['rhi']}\n")
        name[f"{c['chrom']}:{c['rlo']}-{c['rhi']}"] = i
if os.path.exists(fa):
    os.remove(fa)
subprocess.run(f"samtools faidx {GENOME} -r {CACHE}/sibext.txt > {fa} 2>/dev/null", shell=True, check=True)

paf = f"{CACHE}/sibext.paf"
if os.path.exists(paf):
    os.remove(paf)
with open(paf, "w") as fh:
    fh.write(subprocess.run(
        ["minimap2", "-x", "asm20", "-c", "--eqx", "-X", "--no-long-join", "-t", "4",
         "-N", "50", "-p", "0.1", fa, fa], capture_output=True, text=True, check=True).stdout)

# ---- keep only the query sub-interval covered by a SIBLING (same predicted family) -----------------
cov_iv = defaultdict(list)
for ln in open(paf):
    f = ln.split("\t")
    def back(n):
        c, rng = n.rsplit(":", 1)
        a, b = rng.split("-")
        return f"{c}:{int(a)-1}-{b}"
    q, t = back(f[0]), back(f[5])
    if q not in name or t not in name or q == t:
        continue
    qi, ti = name[q], name[t]
    if copies[qi]["fam"] != copies[ti]["fam"]:
        continue
    ident = 1.0 - float(next((x[5:] for x in f[12:] if x.startswith("de:f:")), "1"))
    if ident < MIN_ID:
        continue
    cov_iv[qi].append((int(f[2]), int(f[3])))

for i, c in enumerate(copies):
    iv = sorted(cov_iv.get(i, []))
    if not iv:
        c["slo"], c["shi"] = c["s"], c["e"]      # no sibling support: fall back to the rep span
        continue
    merged = [list(iv[0])]
    for a, b in iv[1:]:
        if a <= merged[-1][1] + 500:
            merged[-1][1] = max(merged[-1][1], b)
        else:
            merged.append([a, b])
    best = max(merged, key=lambda m: m[1] - m[0])
    c["slo"], c["shi"] = c["rlo"] + best[0], c["rlo"] + best[1]

# ---- score all four estimators against truth -------------------------------------------------------
members = []
for ln in open(BED):
    ch, s, e, nm, *_ = ln.rstrip("\n").split("\t")
    if CHROMS and ch not in CHROMS:
        continue
    members.append((ch, int(s), int(e), nm.split("|")[0]))

rows = []
for (ch, ms, me, g) in members:
    ov = [(min(me, c["e"]) - max(ms, c["s"]), c) for c in copies
          if c["chrom"] == ch and ms < c["e"] and me > c["s"]]
    if not ov:
        continue
    c = max(ov, key=lambda x: x[0])[1]
    T = me - ms
    rows.append((g, c["nexon"], (c["e"] - c["s"]) / T, (c["rhi"] - c["rlo"]) / T,
                 (c["shi"] - c["slo"]) / T))


def med(v):
    v = sorted(v)
    return v[len(v) // 2] if v else 0.0


print(f"detected truth members: {len(rows)}   (identity floor for sibling support: {MIN_ID})\n")
print(f"{'estimator':<28}{'median':>9}{'>=0.8x':>9}{'<0.5x':>9}{'>2x':>7}")
for lab, i in (("representative (shipped)", 2), ("read extent, >=3 support", 3),
               ("SIBLING-BOUNDED", 4)):
    v = [r[i] for r in rows]
    print(f"{lab:<28}{med(v):>9.2f}{100*sum(1 for x in v if x>=0.8)/len(v):>8.0f}%"
          f"{100*sum(1 for x in v if x<0.5)/len(v):>8.0f}%{100*sum(1 for x in v if x>2)/len(v):>6.0f}%")

print(f"\n{'gene':<14}{'nexon':>6}{'rep':>7}{'reads':>8}{'sibling':>9}")
for g, nx, a, b, c in sorted(rows, key=lambda r: r[2])[:10]:
    print(f"{g:<14}{nx:>6}{a:>7.2f}{b:>8.2f}{c:>9.2f}")
