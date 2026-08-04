#!/usr/bin/env python3
"""Consensus locus representative as a max-weight path on a READ-WITNESSED splice graph.

Usage: cothread_chain.py <copies.tsv> <bam> [chr1,chr15]

THE DEFECT THIS FIXES. Selecting the maximum-weight set of pairwise NON-OVERLAPPING junctions has no
locality constraint: two junctions from different genes do not overlap, so the DP will happily chain them.
That is why weakly-supported chains blew up (DNM1P50 22x, DNM1P51 12x, AMY2A 6x true size) and why raising
the support floor to 20 reads appeared to help -- it did not fix the chaining, it just discarded two thirds
of the loci. 3 reads per junction is the floor used everywhere else in the pipeline and should be kept.

THE CONSTRAINT. Junction B may follow junction A only if at least one READ carries both. Nodes are
junctions, edges are read-witnessed successions, and the representative is the maximum-weight path. This is
a splice graph, and the DP is a longest path on a DAG (junctions sorted by coordinate), so it stays exact
and threshold-light. Locality is now enforced by the reads rather than by a support threshold.

Three arms are compared on the same loci: the shipped single-exon stub, the unconstrained max-weight set at
a given floor, and the co-observed path at the same floor. Any difference is the constraint alone.
"""
import subprocess
import sys
from bisect import bisect_right
from collections import defaultdict

COP, BAM = sys.argv[1], sys.argv[2]
CHROMS = set(sys.argv[3].split(",")) if len(sys.argv) > 3 else None
BED = "bench/soto/80_fams.gene_preferred.bed"

copies = []
for i, ln in enumerate(open(COP)):
    if i == 0:
        continue
    f = ln.rstrip("\n").split("\t")
    if len(f) < 7 or (CHROMS and f[3] not in CHROMS):
        continue
    copies.append((f[3], int(f[4]), int(f[5]), int(f[6])))

members = []
for ln in open(BED):
    c, s, e, nm, *_ = ln.rstrip("\n").split("\t")
    if CHROMS and c not in CHROMS:
        continue
    members.append((c, int(s), int(e), nm.split("|")[0]))


def read_junctions(chrom, lo, hi):
    """Per-read junction lists (in order), so both support counts and successions are available."""
    out = subprocess.run(["samtools", "view", "-F", "2308", BAM, f"{chrom}:{lo}-{hi}"],
                         capture_output=True, text=True).stdout
    per_read = []
    for ln in out.splitlines():
        f = ln.split("\t")
        p, n, js = int(f[3]), 0, []
        for ch in f[5]:
            if ch.isdigit():
                n = n * 10 + ord(ch) - 48
            else:
                if ch == "N":
                    js.append((p, p + n))
                    p += n
                elif ch in "MDX=":
                    p += n
                n = 0
        if js:
            per_read.append(js)
    return per_read


def weights_and_succ(per_read):
    w = defaultdict(int)
    succ = defaultdict(set)
    for js in per_read:
        for j in js:
            w[j] += 1
        for a, b in zip(js, js[1:]):
            succ[a].add(b)
    return w, succ


def unconstrained(w, succ, floor):
    """Max-weight pairwise non-overlapping set. No locality: this is the version that blows up."""
    iv = sorted(((a, b, c) for (a, b), c in w.items() if c >= floor), key=lambda x: x[1])
    if not iv:
        return []
    ends = [b for _, b, _ in iv]
    best = [0] * (len(iv) + 1)
    take = [None] * len(iv)
    for i, (a, b, c) in enumerate(iv):
        j = bisect_right(ends, a) - 1
        cand = c + (best[j + 1] if j >= 0 else 0)
        if cand > best[i]:
            best[i + 1], take[i] = cand, j
        else:
            best[i + 1], take[i] = best[i], None
    ch, i = [], len(iv) - 1
    while i >= 0:
        if best[i + 1] != best[i]:
            ch.append(iv[i])
            i = take[i] if take[i] is not None else -1
        else:
            i -= 1
    return sorted(ch)


def co_observed(w, succ, floor):
    """Max-weight path where every step A->B is witnessed by a read carrying both."""
    nodes = sorted([j for j, c in w.items() if c >= floor])
    if not nodes:
        return []
    idx = {j: i for i, j in enumerate(nodes)}
    best = [w[j] for j in nodes]
    prev = [None] * len(nodes)
    for i, j in enumerate(nodes):                    # nodes sorted by start; successors come later
        for nxt in succ.get(j, ()):
            k = idx.get(nxt)
            if k is None or k <= i or nxt[0] < j[1]:  # must exist, be later, and not overlap
                continue
            if best[i] + w[nxt] > best[k]:
                best[k] = best[i] + w[nxt]
                prev[k] = i
    end = max(range(len(nodes)), key=lambda i: best[i])
    ch, i = [], end
    while i is not None:
        ch.append((nodes[i][0], nodes[i][1], w[nodes[i]]))
        i = prev[i]
    return sorted(ch)


data = []
for (c, ms, me, g) in members:
    ov = [(min(me, pe) - max(ms, ps), ps, pe, nx) for (pc, ps, pe, nx) in copies
          if pc == c and ms < pe and me > ps]
    if not ov:
        continue
    _, ps, pe, nx = max(ov)
    if nx > 1:
        continue
    pr = read_junctions(c, ps, pe)
    w, succ = weights_and_succ(pr)
    data.append((g, me - ms, (pe - ps) / (me - ms), w, succ))

N = len(data)


def report(lab, fn, floor):
    v, used = [], 0
    for (g, T, stub, w, succ) in data:
        ch = fn(w, succ, floor) if w else []
        if ch:
            v.append((ch[-1][1] - ch[0][0]) / T)
            used += 1
        else:
            v.append(stub)
    v.sort()
    n = len(v)
    ib = sum(1 for x in v if 0.5 <= x <= 2)
    print(f"{lab:<34}{used:>6}{v[n//2]:>9.2f}{100*ib/n:>10.0f}%"
          f"{100*sum(1 for x in v if x>2)/n:>6.0f}%{ib:>10}")


print(f"stub-represented truth members: {N}   (hybrid: chain where found, else keep the stub)\n")
print(f"{'arm':<34}{'used':>6}{'median':>9}{'in 0.5-2x':>11}{'>2x':>6}{'#correct':>10}")
base = sorted(d[2] for d in data)
ibb = sum(1 for x in base if 0.5 <= x <= 2)
print(f"{'shipped (all stubs)':<34}{'-':>6}{base[N//2]:>9.2f}{100*ibb/N:>10.0f}%"
      f"{100*sum(1 for x in base if x>2)/N:>6.0f}%{ibb:>10}")
for fl in (3, 5, 10, 20):
    report(f"unconstrained set, >= {fl}", unconstrained, fl)
for fl in (3, 5, 10, 20):
    report(f"CO-OBSERVED path, >= {fl}", co_observed, fl)
