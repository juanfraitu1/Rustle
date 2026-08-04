#!/usr/bin/env python3
"""Replay the E_r tail offline from a rep-audit dump: edges -> components -> >=2-loci gate -> score.

Usage: replay_er.py <catalog.log with RUSTLE_LOCUS_AUDIT=1> <genome.fa> [chr1,chr15]

WHY THIS IS NOT THE PROXY THAT FAILED BEFORE. The 2026-08-02 identity-floor prediction (offline 0.553,
pipeline 0.433) was computed by regrouping `copies.fa` -- the set of copies ALREADY EMITTED under the old
rule, i.e. a node set pre-filtered by the very gate the change perturbs. That is an upper bound by
construction. This script instead consumes the `[rep-audit]` dump, which is the complete pre-gate node set
that `distinct_locus_reps` is about to run on, so the gate can be applied here exactly as the pipeline
applies it.

It still is not the pipeline: the partition here is plain connected components, while the shipped path runs
a gamma-quasi-clique decomposition. So VALIDATE against a real arm before trusting a delta -- the script
prints the emitted copy/family counts precisely so they can be compared to the run's own totals.

The expensive part (all-vs-all minimap2) is cached to disk keyed by substrate, so sweeping identity and
coverage after the first call costs milliseconds instead of 18 minutes per point.
"""
import os
import subprocess
import sys
from itertools import combinations

LOG, GENOME = sys.argv[1], sys.argv[2]
CHROMS = set(sys.argv[3].split(",")) if len(sys.argv) > 3 else None
BED = "bench/soto/80_fams.gene_preferred.bed"
CACHE = "/tmp/claude-1000/-mnt-c-Users-jfris-Desktop/161b66ca-e17e-4e70-b1cf-8b479d6c328d/scratchpad"

# ---- the complete pre-gate node set -----------------------------------------------------------------
reps = []
for ln in open(LOG, errors="replace"):
    if not ln.startswith("[rep-audit]"):
        continue
    f = ln.rstrip("\n").split("\t")
    chrom, s, e = f[2], int(f[3]), int(f[4])
    if CHROMS and chrom not in CHROMS:
        continue
    reps.append({"chrom": chrom, "s": s, "e": e, "strand": f[5],
                 "n_reads": int(f[6]), "nexon": int(f[8])})
idx = {f"{r['chrom']}:{r['s']}-{r['e']}": i for i, r in enumerate(reps)}
print(f"reps in dump: {len(reps)}   (single-exon {sum(1 for r in reps if r['nexon']==1)})")

members = []
for ln in open(BED):
    c, ms, me, name, *_ = ln.rstrip("\n").split("\t")
    if CHROMS and c not in CHROMS:
        continue
    members.append((c, int(ms), int(me), name.split("|")[0], name.split("|")[1]))

# ---- genomic-span substrate, aligned once and cached -------------------------------------------------
span_fa = f"{CACHE}/replay_spans.fa"
if not os.path.exists(span_fa):
    open(f"{CACHE}/replay_regions.txt", "w").write(
        "\n".join(f"{r['chrom']}:{r['s']+1}-{r['e']}" for r in reps) + "\n")
    subprocess.run(f"samtools faidx {GENOME} -r {CACHE}/replay_regions.txt > {span_fa}",
                   shell=True, check=True)

paf_path = f"{CACHE}/replay_spans.paf"
if not os.path.exists(paf_path):
    with open(paf_path, "w") as fh:
        for tier in (["-x", "asm20"], ["-k", "11", "-w", "5"]):
            fh.write(subprocess.run(
                ["minimap2", *tier, "-c", "--eqx", "-X", "--no-long-join", "-t", "4",
                 "-N", "50", "-p", "0.1", span_fa, span_fa],
                capture_output=True, text=True, check=True).stdout)

best = {}
for ln in open(paf_path):
    f = ln.split("\t")
    def norm(n):
        c, rng = n.rsplit(":", 1)
        a, b = rng.split("-")
        return f"{c}:{int(a)-1}-{b}"
    q, t = norm(f[0]), norm(f[5])
    if q not in idx or t not in idx or q == t:
        continue
    cov = (int(f[3]) - int(f[2])) / min(int(f[1]), int(f[6]))
    ident = 1.0 - float(next((x[5:] for x in f[12:] if x.startswith("de:f:")), "1"))
    k = (idx[q], idx[t]) if idx[q] < idx[t] else (idx[t], idx[q])
    if k not in best or cov > best[k][0]:
        best[k] = (cov, ident)
print(f"cached alignment records: {len(best)} rep pairs\n")


def components(edges, n):
    p = list(range(n))
    def find(x):
        while p[x] != x:
            p[x] = p[p[x]]
            x = p[x]
        return x
    for a, b in edges:
        ra, rb = find(a), find(b)
        if ra != rb:
            p[ra] = rb
    out = {}
    for i in range(n):
        out.setdefault(find(i), []).append(i)
    return list(out.values())


def distinct_loci(group):
    """>= 2 SPATIALLY DISTINCT loci: same-chromosome overlapping reps count once."""
    seen = []
    for i in group:
        r = reps[i]
        if not any(r["chrom"] == c and r["s"] < e and r["e"] > s for (c, s, e) in seen):
            seen.append((r["chrom"], r["s"], r["e"]))
    return len(seen)


def score(min_id, min_cov):
    E = [k for k, (cov, ident) in best.items() if cov >= min_cov and ident >= min_id]
    fams = [g for g in components(E, len(reps)) if len(g) >= 2 and distinct_loci(g) >= 2]
    lab = {}
    for fi, g in enumerate(fams):
        for i in g:
            lab[i] = fi
    # map each truth member to the family of the best-overlapping emitted rep
    pred, det = {}, 0
    for mi, (c, ms, me, g, fam) in enumerate(members):
        ov = [(min(me, reps[i]["e"]) - max(ms, reps[i]["s"]), lab[i])
              for i in lab if reps[i]["chrom"] == c and ms < reps[i]["e"] and me > reps[i]["s"]]
        if ov:
            pred[g] = f"F{max(ov)[1]}"
            det += 1
        else:
            pred[g] = f"__unplaced_{mi}"
    tp = fp = fn = 0
    for (a, b) in combinations(members, 2):
        st, sp = a[4] == b[4], pred[a[3]] == pred[b[3]]
        tp += st and sp
        fp += sp and not st
        fn += st and not sp
    P = tp / (tp + fp) if tp + fp else 0.0
    R = tp / (tp + fn) if tp + fn else 0.0
    F1 = 2 * P * R / (P + R) if P + R else 0.0
    return len(fams), sum(len(g) for g in fams), det, P, R, F1


print(f"{'id':>5} {'cov':>5} {'fams':>6} {'copies':>7} {'detect':>7} {'P':>7} {'R':>7} {'F1':>7}")
for min_id in (0.70, 0.80, 0.90, 0.93, 0.95, 0.98):
    for min_cov in (0.30, 0.50, 0.70):
        nf, nc, det, P, R, F1 = score(min_id, min_cov)
        print(f"{min_id:5.2f} {min_cov:5.2f} {nf:6d} {nc:7d} {det:7d} {P:7.3f} {R:7.3f} {F1:7.3f}")
