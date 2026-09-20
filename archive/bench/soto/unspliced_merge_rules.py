#!/usr/bin/env python3
"""Is there an EMPIRICAL, ground-truth-testable rule for merging unspliced read-clusters?

Section 12: unspliced (single-exon) copies cannot span a locus -- their extent is bounded by read length, so
GUSBP1's 231 kb gene is covered by 11 separate ~9 kb clusters. Merging them by proximity would recover the
size but reintroduce readthrough over-merge. So: is there a rule that separates "these two clusters are the
same locus" from "these two clusters are different loci", judged against the Soto members?

CRITICAL constraint on candidate rules: `cluster_unspliced` is single-linkage span-overlap clustering, so if
ANY read overlapped both clusters they would already BE one cluster. No unspliced read can bridge two
clusters -- that witness is unavailable by construction. The evidence has to come from reads the unspliced
path never used:

  SPLICED-BRIDGE  n spliced primary reads with aligned blocks in BOTH clusters  <- a certificate, not a
                  threshold: one such read is direct evidence of a single transcription unit
  GAP-COVERAGE    is the intervening region covered, or is it a real gap?
  GAP-BP          raw distance (the naive proximity rule, as a baseline to beat)
  SAME-STRAND     cheap necessary condition

Labels: POSITIVE = both clusters inside the SAME Soto member; NEGATIVE = inside DIFFERENT members. Negatives
are restricted to the same chromosome within MAX_GAP so the decision is non-trivial (distant pairs are easy).
"""
import subprocess, re, sys, itertools, random
from collections import defaultdict
import numpy as np

SAM  = "/home/juanfra/miniforge3/bin/samtools"
BAM  = "/home/juanfra/winloci_scratch/soto_cache/soto_regions.bam"
CAT  = "/home/juanfra/winloci_scratch/soto_cache/definitive.copies.tsv"
BED  = "bench/soto/80_fams.gene_preferred.bed"
MAX_GAP   = 300_000
MAX_PAIRS = 400
random.seed(0)

truth = []
for ln in open(BED):
    c, s, e, name, *_ = ln.rstrip("\n").split("\t")
    truth.append((name.split("|")[0], c, int(s), int(e)))

# unspliced copies only, attributed to the truth member that contains their midpoint
clusters = defaultdict(list)      # truth-member key -> [(chrom,start,end,strand)]
for i, ln in enumerate(open(CAT)):
    if i == 0: continue
    f = ln.rstrip("\n").split("\t")
    if len(f) < 9 or int(f[6]) != 1: continue         # n_exon == 1
    c, s, e, strand = f[3], int(f[4]), int(f[5]), f[7]
    mid = (s + e) // 2
    for k, (gn, tc, ts, te) in enumerate(truth):
        if tc == c and ts <= mid <= te:
            clusters[k].append((c, s, e, strand)); break

pos, neg = [], []
for k, v in clusters.items():
    for a, b in itertools.combinations(sorted(v, key=lambda x: x[1]), 2):
        pos.append((a, b, truth[k][0], truth[k][0]))
keys = sorted(clusters)
for k1, k2 in itertools.combinations(keys, 2):
    if truth[k1][1] != truth[k2][1]: continue
    for a in clusters[k1]:
        for b in clusters[k2]:
            lo, hi = sorted([a, b], key=lambda x: x[1])
            if 0 < hi[1] - lo[2] < MAX_GAP:
                neg.append((lo, hi, truth[k1][0], truth[k2][0]))
if len(pos) > MAX_PAIRS: pos = random.sample(pos, MAX_PAIRS)
if len(neg) > MAX_PAIRS: neg = random.sample(neg, MAX_PAIRS)
print(f"pairs: POSITIVE (same Soto member) {len(pos)}   NEGATIVE (different members, gap<{MAX_GAP//1000}kb) {len(neg)}")

def blocks(pos0, cig):
    """reference blocks of an alignment, splitting at N (introns)"""
    out, ref, cur = [], pos0, pos0
    for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cig):
        n = int(n)
        if op in "M=X": ref += n
        elif op == "D": ref += n
        elif op == "N":
            out.append((cur, ref)); ref += n; cur = ref
    out.append((cur, ref))
    return out

def features(a, b):
    c = a[0]; lo, hi = min(a[1], b[1]), max(a[2], b[2])
    gap = max(0, max(a[1], b[1]) - min(a[2], b[2]))
    out = subprocess.run([SAM, "view", "-F", "2308", BAM, f"{c}:{lo+1}-{hi}"],
                         capture_output=True, text=True).stdout
    bridge = 0; gapcov = 0; adepth = 0; bdepth = 0
    gs, ge = min(a[2], b[2]), max(a[1], b[1])
    for ln in out.splitlines():
        f = ln.split("\t")
        if len(f) < 6: continue
        p0 = int(f[3]) - 1; cig = f[5]
        blks = blocks(p0, cig)
        ina = any(x < a[2] and a[1] < y for x, y in blks)
        inb = any(x < b[2] and b[1] < y for x, y in blks)
        if ina: adepth += 1
        if inb: bdepth += 1
        if ina and inb and "N" in cig:
            bridge += 1                                  # SPLICED read with blocks in both clusters
        if ge > gs and any(x < ge and gs < y for x, y in blks):
            gapcov += 1
    denom = max(min(adepth, bdepth), 1)
    return dict(gap=gap, bridge=bridge, bridge_frac=bridge/denom,
                gapcov=gapcov, gapcov_frac=gapcov/denom,
                same_strand=int(a[3] == b[3]))

rows = []
for lab, pairs in ((1, pos), (0, neg)):
    for a, b, ga, gb in pairs:
        rows.append((lab, features(a, b), ga, gb))

def report(name, key, thresholds):
    print(f"\n  rule: {name}")
    for t in thresholds:
        tp = sum(1 for l, f, *_ in rows if l == 1 and f[key] >= t)
        fp = sum(1 for l, f, *_ in rows if l == 0 and f[key] >= t)
        fn = sum(1 for l, f, *_ in rows if l == 1 and f[key] < t)
        prec = tp / max(tp + fp, 1); rec = tp / max(tp + fn, 1)
        f1 = 2 * prec * rec / max(prec + rec, 1e-9)
        print(f"    >= {t:<8} TP={tp:<4} FP={fp:<4} FN={fn:<4}  precision={prec:.3f} recall={rec:.3f} F1={f1:.3f}")

def report_le(name, key, thresholds):
    print(f"\n  rule: {name}")
    for t in thresholds:
        tp = sum(1 for l, f, *_ in rows if l == 1 and f[key] <= t)
        fp = sum(1 for l, f, *_ in rows if l == 0 and f[key] <= t)
        fn = sum(1 for l, f, *_ in rows if l == 1 and f[key] > t)
        prec = tp / max(tp + fp, 1); rec = tp / max(tp + fn, 1)
        f1 = 2 * prec * rec / max(prec + rec, 1e-9)
        print(f"    <= {t:<8} TP={tp:<4} FP={fp:<4} FN={fn:<4}  precision={prec:.3f} recall={rec:.3f} F1={f1:.3f}")

report("SPLICED-BRIDGE (n spliced reads with blocks in both)", "bridge", [1, 2, 3, 5, 10])
report("GAP-COVERAGE (n reads covering the intervening gap)", "gapcov", [1, 3, 10, 30])
report_le("GAP-BP (naive proximity baseline)", "gap", [1000, 10_000, 50_000, 150_000])
print("\n  baseline: always-merge  precision=%.3f" % (len(pos)/max(len(pos)+len(neg),1)))

# --- GAP-MATCHED CONTROL -------------------------------------------------------------------------------
# The gap rule above may be an artifact: positives are pairs INSIDE one member (gap <= member length) while
# negatives are BETWEEN members, so "gap" partly encodes "member length vs inter-member distance" rather than
# same-locus evidence. Re-evaluate inside gap bins, where that confound is held constant. A rule that only
# works across bins is measuring the sampling, not the biology.

BINS = [(0, 10_000), (10_000, 50_000), (50_000, 150_000), (150_000, 300_000)]
print("\n=== GAP-MATCHED CONTROL: discrimination WITHIN gap bins ===")
print(f"{'gap bin':>18} {'n+':>5} {'n-':>5} {'base':>6} {'bridge>=1 prec':>15} {'rec':>6} {'gapcov>=3 prec':>15}")
for lo, hi in BINS:
    P = [r for r in rows if r[0] == 1 and lo <= r[1]['gap'] < hi]
    N = [r for r in rows if r[0] == 0 and lo <= r[1]['gap'] < hi]
    if not P and not N: continue
    base = len(P) / max(len(P) + len(N), 1)
    tp = sum(1 for r in P if r[1]['bridge'] >= 1); fp = sum(1 for r in N if r[1]['bridge'] >= 1)
    prec = tp / max(tp + fp, 1); rec = tp / max(len(P), 1)
    tp2 = sum(1 for r in P if r[1]['gapcov'] >= 3); fp2 = sum(1 for r in N if r[1]['gapcov'] >= 3)
    prec2 = tp2 / max(tp2 + fp2, 1)
    print(f"{f'{lo//1000}-{hi//1000}kb':>18} {len(P):>5} {len(N):>5} {base:>6.3f} "
          f"{prec:>15.3f} {rec:>6.3f} {prec2:>15.3f}")

print("\n=== COMBINED RULES (vs always-merge baseline %.3f) ===" % (len(pos)/max(len(pos)+len(neg),1)))
def ev(name, fn):
    tp = sum(1 for l, f, *_ in rows if l == 1 and fn(f)); fp = sum(1 for l, f, *_ in rows if l == 0 and fn(f))
    fn_ = sum(1 for l, f, *_ in rows if l == 1 and not fn(f))
    prec = tp / max(tp + fp, 1); rec = tp / max(tp + fn_, 1)
    print(f"  {name:<44} prec={prec:.3f} rec={rec:.3f} F1={2*prec*rec/max(prec+rec,1e-9):.3f}")
ev("bridge>=1", lambda f: f['bridge'] >= 1)
ev("gap<=50kb", lambda f: f['gap'] <= 50_000)
ev("bridge>=1 OR gap<=50kb", lambda f: f['bridge'] >= 1 or f['gap'] <= 50_000)
ev("bridge>=1 AND same_strand", lambda f: f['bridge'] >= 1 and f['same_strand'])
ev("bridge>=1 OR gap<=10kb", lambda f: f['bridge'] >= 1 or f['gap'] <= 10_000)
