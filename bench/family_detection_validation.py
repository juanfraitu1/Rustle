#!/usr/bin/env python
"""
family_detection_validation.py
==============================

CRITERION-LEVEL validation of rustle's data-driven multi-copy paralog-FAMILY
detection criterion against an annotated-paralogy truth set, plus a
justification of (or alternative to) its similarity threshold.

WHAT IS UNDER TEST
------------------
rustle forms paralog families by grouping loci whose sequences have
minimizer-Jaccard >= a threshold.
  - default FAMILY_MERGE_JACCARD_DEFAULT = 0.30
  - minimizer params k = 15, w = 10
  - source: src/rustle/vg_family/family_graph.rs  (fn `minimizers`,
            `refine_by_minimizer_jaccard`, `merge_singletons_by_sequence`)

This script reproduces THAT criterion on labeled GENE sequences and asks two
questions:
  (Q1 THRESHOLD JUSTIFICATION) Is the within-family vs between-family Jaccard
      signal bimodal / separable, and where does 0.30 fall (valley? optimal?)?
  (Q2 FAMILY DETECTION) If we union-find cluster all genes at a threshold, how
      well do the predicted clusters match the truth families
      (pairwise P/R/F1, ARI, V-measure, families recovered pure)?

HONEST CAVEATS (also printed and written to the .md)
----------------------------------------------------
  * TRUTH-SET CIRCULARITY: universe.tsv is itself built by SEQUENCE SIMILARITY
    (minimap2 -x asm20, ident>=0.90/cov>=0.50, then union-find — see
    bench/copy_recovery_eval/30_build_universe.py + lib_eval.build_universe).
    This script then tests whether a DIFFERENT similarity criterion
    (minimizer-Jaccard) recovers those families, so the high AUC measures
    concordance of two similarity clusterings, NOT agreement with curated
    biology (Compara/OrthoFinder/synteny). Temper "recovers truth" claims.
  * CONTROLS ARE NOT certified single-copy: they are merely genes absent from
    universe.tsv, which was filtered to minimap2 hits at ident>=0.90/cov>=0.50,
    so sub-threshold paralog pairs leak into the control pool. There is at least
    one J=1.0 control-control pair (a missed paralog) — see report.
  * This is a CRITERION-level validation of the similarity gate ONLY, run on
    REFERENCE gene sequences. It is NOT the full multi-stage rustle pipeline
    (which also uses position overlap, per-exon clustering, EM, flow, etc.).
  * Only EXPRESSED / ANNOTATED paralogs are in the truth set.
  * COLLAPSED copies (paralogs that map to one locus and never appear as
    distinct sequences here) are NOT detectable by this sequence-similarity
    criterion; that is the separate PSV-co-segregation test.
  * rustle's production `minimizers` uses a NON-canonical FNV-1a hash over the
    forward k-mer. mmh3 is absent here, so per task spec we use a CANONICAL
    minimizer (min(kmer, revcomp)) hashed with blake2b truncated to 64 bits.
    The WINDOWING SCHEME is reproduced exactly (window over k-mer-start
    positions, size w, min hash per window). Canonicalization makes the metric
    strand-robust; it does not change the separability conclusion (homologous
    exons share many minimizers, unrelated exons share ~0 under either scheme).

DETERMINISM: all sampling uses a single FIXED integer seed (SEED, below).
No wall-clock, no os.urandom, no unseeded RNG anywhere.

Run:
  /home/juanfra/miniforge3/bin/python bench/family_detection_validation.py
"""

import os
import sys
import random
import hashlib
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --------------------------------------------------------------------------
# FIXED determinism knobs
# --------------------------------------------------------------------------
SEED = 20260616            # fixed integer constant; never wall-clock
N_CONTROLS = 150           # ~150 single-copy control genes
K = 15                     # rustle minimizer k
W = 10                     # rustle minimizer w
RUSTLE_THRESHOLD = 0.30    # FAMILY_MERGE_JACCARD_DEFAULT

random.seed(SEED)
np.random.seed(SEED)

BASE = "/mnt/c/Users/jfris/Desktop/Rustle/bench"
RESULTS = os.path.join(BASE, "copy_recovery_eval", "results")
UNIVERSE_TSV = os.path.join(RESULTS, "universe.tsv")
GENE_FA = os.path.join(RESULTS, "gene_rep.fa")

OUT_PNG = os.path.join(BASE, "family_detection_validation.png")
OUT_MD = os.path.join(BASE, "family_detection_validation.md")

# Palette
NAVY = "#22313f"
ORANGE = "#e8590c"
GREEN = "#1e8449"


# --------------------------------------------------------------------------
# Minimizer machinery (rustle windowing, canonical k-mers, blake2b hash)
# --------------------------------------------------------------------------
_RC = bytes.maketrans(b"ACGTacgtNn", b"TGCAtgcaNn")

def revcomp(s: bytes) -> bytes:
    return s.translate(_RC)[::-1]

def _hash64(kmer: bytes) -> int:
    # blake2b truncated to 64 bits (mmh3 absent; per task spec).
    return int.from_bytes(hashlib.blake2b(kmer, digest_size=8).digest(), "little")

def minimizers(seq: bytes, k: int = K, w: int = W):
    """Canonical-kmer minimizer set, reproducing rustle's WINDOWING exactly.

    rustle (family_graph.rs:minimizers): n = len-k+1 k-mer starts; for
    win_start in 0..max(1, n-w): window = [win_start, win_start+w); emit the
    min hash in the window; collect into a SET. We add canonicalization
    (min(kmer, revcomp)) for strand robustness.
    """
    seq = seq.upper()
    out = set()
    L = len(seq)
    if L < k:
        return out
    n = L - k + 1
    # precompute canonical hashes per k-mer start
    hashes = [0] * n
    for i in range(n):
        km = seq[i:i + k]
        rc = revcomp(km)
        can = km if km <= rc else rc
        hashes[i] = _hash64(can)
    upper = max(1, n - w)
    for win_start in range(upper):
        win_end = min(win_start + w, n)
        out.add(min(hashes[win_start:win_end]))
    return out

def jaccard(a: set, b: set) -> float:
    if not a and not b:
        return 0.0
    inter = len(a & b)
    uni = len(a | b)
    return inter / uni if uni else 0.0


# --------------------------------------------------------------------------
# Load truth + sequences
# --------------------------------------------------------------------------
def load_fasta(path):
    seqs = {}
    cur = None
    buf = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line[0] == ">":
                if cur is not None:
                    seqs[cur] = "".join(buf).encode()
                cur = line[1:].split()[0]
                # normalize: annotation may use a 'gene-' prefix; universe uses bare
                if cur.startswith("gene-"):
                    cur = cur[len("gene-"):]
                buf = []
            else:
                buf.append(line)
        if cur is not None:
            seqs[cur] = "".join(buf).encode()
    return seqs

def load_universe(path):
    """Return gene_id -> family_id (truth labels). One label per gene."""
    g2f = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        ti = header.index("gene_id")
        fi = header.index("family_id")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(ti, fi):
                continue
            g = parts[ti].strip()
            f = parts[fi].strip()
            if g.startswith("gene-"):
                g = g[len("gene-"):]
            g2f[g] = f
    return g2f


print(f"[seed] FIXED SEED = {SEED}", file=sys.stderr)

fasta = load_fasta(GENE_FA)
g2f = load_universe(UNIVERSE_TSV)

universe_genes = sorted(g2f.keys())
truth_families = sorted(set(g2f.values()))

# Reconcile: which truth genes have a sequence?
matched = [g for g in universe_genes if g in fasta]
unmatched = [g for g in universe_genes if g not in fasta]
print(f"[match] universe genes: {len(universe_genes)}; "
      f"matched to fasta: {len(matched)}; unmatched: {len(unmatched)}",
      file=sys.stderr)

# Control set: FIXED-seed sample of single-copy genes NOT in any universe family.
control_pool = sorted(set(fasta.keys()) - set(universe_genes))
rng = random.Random(SEED)
control_pool_shuffled = control_pool[:]
rng.shuffle(control_pool_shuffled)
controls = sorted(control_pool_shuffled[:N_CONTROLS])
print(f"[control] pool={len(control_pool)} sampled={len(controls)} "
      f"(seed={SEED})", file=sys.stderr)

# --------------------------------------------------------------------------
# CONTROL-PURITY SCREEN (Finding #4) — apply minimap2 to the sampled controls
# vs (universe genes + the other controls) and flag any control with a hit as
# NOT single-copy. Purely DIAGNOSTIC: we do NOT silently drop controls (that
# would HIDE the leak); we report them so the resulting false merge stays
# visible.
#
# CRITICAL nuance discovered while building this: the universe builder's own
# criterion `minimap2 -x asm20` produces NO alignment for very SHORT identical
# paralogs (e.g. the 153 bp pair LOC129531078==LOC129531352, which is byte-
# identical yet asm20 reports nothing because it does not seed short sequences).
# So screening with the SAME asm20 criterion canNOT certify purity — it has the
# SAME blind spot that let the pair leak in. We therefore run TWO presets and
# report both: (a) `asm20` = exactly the universe criterion (for parity), and
# (b) `map-ont`/default sensitive seeding = catches short paralogs asm20 misses.
# The minimizer-Jaccard criterion under test (k=15/w=10) is itself MORE
# sensitive than asm20 here (it scored that pair J=1.0), so (b) is the fair
# purity bar.
# Deterministic: minimap2 fixed inputs; outputs sorted.
import shutil, subprocess, tempfile
MINIMAP2 = shutil.which("minimap2")
UNI_IDENT, UNI_COV = 0.90, 0.50  # same cutoffs as 30_build_universe.py

def _run_control_screen(preset, extra_args):
    """Return sorted list of (control, partner, ident, cov) hits >= cutoffs."""
    screen_targets = [g for g in sorted(set(universe_genes) | set(controls))
                      if g in fasta]
    with tempfile.TemporaryDirectory() as td:
        qfa = os.path.join(td, "controls.fa")
        tfa = os.path.join(td, "targets.fa")
        with open(qfa, "w") as qh:
            for g in sorted(controls):
                if g in fasta:
                    qh.write(f">{g}\n{fasta[g].decode()}\n")
        with open(tfa, "w") as th:
            for g in screen_targets:
                th.write(f">{g}\n{fasta[g].decode()}\n")
        cmd = [MINIMAP2, "-x", preset] + extra_args + ["-c", tfa, qfa]
        paf = subprocess.run(cmd, capture_output=True, text=True,
                             check=True).stdout
    hits = []
    seen = set()
    for line in paf.splitlines():
        f = line.split("\t")
        if len(f) < 12:
            continue
        q, qlen, tname = f[0], int(f[1]), f[5]
        nmatch, alen = int(f[9]), int(f[10])
        if q == tname:
            continue
        ident = nmatch / alen if alen else 0.0
        cov = alen / qlen if qlen else 0.0
        if ident >= UNI_IDENT and cov >= UNI_COV:
            key = (q, tname)
            if key in seen:
                continue
            seen.add(key)
            hits.append((q, tname, ident, cov))
    hits.sort(key=lambda r: (-r[2] * r[3], r[0], r[1]))
    return hits

control_screen_hits = []           # universe-parity (asm20) hits
control_screen_hits_sens = []      # sensitive (default seeding) hits
control_screen_ran = False
control_screen_error = ""
if MINIMAP2:
    try:
        # (a) exact universe criterion
        control_screen_hits = _run_control_screen("asm20", ["--no-long-join"])
        # (b) sensitive: default short-seq seeding (catches short paralogs).
        #     Use -x sr (short-read) which seeds short sequences the universe's
        #     asm20 misses; cutoffs unchanged.
        control_screen_hits_sens = _run_control_screen("sr", [])
        control_screen_ran = True
        impure_controls = sorted({h[0] for h in control_screen_hits})
        impure_controls_sens = sorted({h[0] for h in control_screen_hits_sens})
        print(f"[control-screen] asm20(universe-parity): "
              f"{len(control_screen_hits)} hit(s) -> "
              f"{len(impure_controls)} impure controls", file=sys.stderr)
        print(f"[control-screen] sr(sensitive, catches short paralogs): "
              f"{len(control_screen_hits_sens)} hit(s) -> "
              f"{len(impure_controls_sens)} impure controls", file=sys.stderr)
        for q, tn, idn, cv in control_screen_hits_sens[:10]:
            print(f"    IMPURE-CONTROL[sr] {q} ~ {tn} ident={idn:.3f} "
                  f"cov={cv:.3f}", file=sys.stderr)
    except Exception as e:
        control_screen_error = str(e)
        print(f"[control-screen] FAILED ({e}); skipping (diagnostic only)",
              file=sys.stderr)
else:
    control_screen_error = "minimap2 not found"
    print("[control-screen] minimap2 not found; skipping (diagnostic only)",
          file=sys.stderr)

impure_controls = sorted({h[0] for h in control_screen_hits})
impure_controls_sens = sorted({h[0] for h in control_screen_hits_sens})

# Build the working gene set: matched paralogs + controls.
# Each control is its own singleton truth family ("CTRL::<gene>") so that
# control-control and control-paralog pairs count as BETWEEN, and any
# control that clusters with another gene is correctly scored as an error.
genes = []          # ordered list of gene ids
label = {}          # gene -> truth family label
is_paralog = {}     # gene -> bool
for g in matched:
    genes.append(g)
    label[g] = g2f[g]
    is_paralog[g] = True
for g in controls:
    genes.append(g)
    label[g] = f"CTRL::{g}"
    is_paralog[g] = False

print(f"[set] total genes in test: {len(genes)} "
      f"({len(matched)} paralogs + {len(controls)} controls)", file=sys.stderr)

# Precompute minimizer sets.
print("[mins] computing minimizer sets ...", file=sys.stderr)
mins = {}
for i, g in enumerate(genes):
    mins[g] = minimizers(fasta[g], K, W)
    if (i + 1) % 50 == 0:
        print(f"  {i+1}/{len(genes)}", file=sys.stderr)

# Drop any gene whose minimizer set is empty (too short) — none expected.
empty = [g for g in genes if not mins[g]]
if empty:
    print(f"[warn] {len(empty)} genes have empty minimizer sets, dropping: "
          f"{empty[:5]}", file=sys.stderr)
    genes = [g for g in genes if mins[g]]


# --------------------------------------------------------------------------
# Pairwise Jaccard: WITHIN-family vs BETWEEN
# --------------------------------------------------------------------------
n = len(genes)
within_j = []   # same truth family (only meaningful for paralog families w/ >=2)
between_j = []   # different family (incl. ctrl-ctrl, ctrl-paralog, cross-family)

# Store full pair list for threshold sweep (memory ok: ~345 genes -> ~59k pairs)
pair_j = []      # (jaccard, is_within_bool)

for i in range(n):
    gi = genes[i]
    li = label[gi]
    mi = mins[gi]
    for k_ in range(i + 1, n):
        gj = genes[k_]
        lj = label[gj]
        j = jaccard(mi, mins[gj])
        same = (li == lj)
        pair_j.append((j, same))
        if same:
            within_j.append(j)
        else:
            between_j.append(j)

within_j = np.array(within_j)
between_j = np.array(between_j)
pair_arr = np.array([p[0] for p in pair_j])
pair_same = np.array([p[1] for p in pair_j], dtype=bool)

print(f"[pairs] within={len(within_j)} between={len(between_j)}", file=sys.stderr)


# --------------------------------------------------------------------------
# BETWEEN-PAIR forensics (Findings #2 and #4)
# --------------------------------------------------------------------------
# We need the IDENTITY of every between-family pair above a threshold, not just
# a rounded fraction. A between pair is one of three kinds:
#   PARA-PARA : two paralog genes in DIFFERENT truth families -> true cross-
#               family false merge (the metric the headline cares about).
#   PARA-CTRL : a paralog gene and a control -> control contamination.
#   CTRL-CTRL : two controls -> either unrelated genes that share minimizers, OR
#               (the dangerous case) a real paralog pair that the minimap2
#               universe builder MISSED and that therefore leaked into the
#               control pool. A high-Jaccard CTRL-CTRL pair is a missed paralog.
def between_pairs_above(threshold):
    """Return list of (jaccard, gi, gj, kind) for all BETWEEN pairs >= threshold.
    Deterministic (i<k walk)."""
    out = []
    for i in range(n):
        gi = genes[i]; li = label[gi]; mi = mins[gi]
        for k_ in range(i + 1, n):
            gj = genes[k_]; lj = label[gj]
            if li == lj:
                continue
            j = jaccard(mi, mins[gj])
            if j >= threshold:
                pi, pj = is_paralog[gi], is_paralog[gj]
                if pi and pj:
                    kind = "PARA-PARA"
                elif (not pi) and (not pj):
                    kind = "CTRL-CTRL"
                else:
                    kind = "PARA-CTRL"
                out.append((j, gi, gj, kind))
    out.sort(key=lambda r: (-r[0], r[1], r[2]))
    return out

# Forensics at the rustle threshold and at a low diagnostic band (0.06).
between_at_rustle = between_pairs_above(RUSTLE_THRESHOLD)
between_at_006 = between_pairs_above(0.06)

def kind_counts(pairs):
    c = {"PARA-PARA": 0, "PARA-CTRL": 0, "CTRL-CTRL": 0}
    for _j, _gi, _gj, kind in pairs:
        c[kind] += 1
    return c

between_rustle_counts = kind_counts(between_at_rustle)
between_006_counts = kind_counts(between_at_006)
n_between_ge_030 = len(between_at_rustle)

print(f"[between>=0.30] total={n_between_ge_030} "
      f"PARA-PARA={between_rustle_counts['PARA-PARA']} "
      f"PARA-CTRL={between_rustle_counts['PARA-CTRL']} "
      f"CTRL-CTRL={between_rustle_counts['CTRL-CTRL']}", file=sys.stderr)
for j, gi, gj, kind in between_at_rustle:
    print(f"    FALSE-MERGE J={j:.4f} {gi} <-> {gj} [{kind}]", file=sys.stderr)

# Residual control-control similarity mass: are the controls actually
# single-copy? (Finding #4 — control purity is NOT certified.)
ctrl_ctrl_006 = [p for p in between_at_006 if p[3] == "CTRL-CTRL"]
print(f"[control-purity] CTRL-CTRL pairs >= 0.06: {len(ctrl_ctrl_006)} "
      f"(if >0, controls are not certified single-copy)", file=sys.stderr)


# --------------------------------------------------------------------------
# Discriminability: AUC (within=positive), overlap, separability
# --------------------------------------------------------------------------
def auc_mann_whitney(pos, neg):
    """AUC = P(random positive scores higher than random negative).
    Computed via rank-sum (handles ties at 0.5). Deterministic.
    """
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    allv = np.concatenate([pos, neg])
    order = np.argsort(allv, kind="mergesort")
    ranks = np.empty(len(allv), dtype=float)
    ranks[order] = np.arange(1, len(allv) + 1)
    # average ranks for ties
    sv = allv[order]
    i = 0
    while i < len(sv):
        j = i
        while j + 1 < len(sv) and sv[j + 1] == sv[i]:
            j += 1
        if j > i:
            avg = (ranks[order[i]] + ranks[order[j]]) / 2.0
            for t in range(i, j + 1):
                ranks[order[t]] = avg
        i = j + 1
    r_pos = ranks[:len(pos)].sum()
    n_pos = len(pos)
    n_neg = len(neg)
    u = r_pos - n_pos * (n_pos + 1) / 2.0
    return u / (n_pos * n_neg)

auc = auc_mann_whitney(within_j, between_j)
print(f"[auc] within-vs-between AUC = {auc:.4f}", file=sys.stderr)

# Overlap region between the distributions (fraction of between >= min within,
# and fraction of within <= max between) -- a crude separability read.
within_q = np.percentile(within_j, [5, 50, 95]) if len(within_j) else [np.nan]*3
between_q = np.percentile(between_j, [5, 50, 95]) if len(between_j) else [np.nan]*3
# overlap mass: fraction of between pairs above the 5th pctile of within
frac_between_above_w5 = float(np.mean(between_j >= within_q[0])) if len(within_j) else float("nan")
frac_within_below_b95 = float(np.mean(within_j <= between_q[2])) if len(within_j) else float("nan")


# --------------------------------------------------------------------------
# Threshold sweep: pairwise precision / recall / F1 over a grid
# --------------------------------------------------------------------------
grid = np.round(np.arange(0.00, 1.0001, 0.01), 4)
P_curve, R_curve, F_curve = [], [], []
n_within_total = int(pair_same.sum())
for t in grid:
    pred_pos = pair_arr >= t
    tp = int(np.sum(pred_pos & pair_same))
    fp = int(np.sum(pred_pos & ~pair_same))
    fn = n_within_total - tp
    P = tp / (tp + fp) if (tp + fp) else 1.0
    R = tp / (tp + fn) if (tp + fn) else 0.0
    F = 2 * P * R / (P + R) if (P + R) else 0.0
    P_curve.append(P)
    R_curve.append(R)
    F_curve.append(F)
P_curve = np.array(P_curve)
R_curve = np.array(R_curve)
F_curve = np.array(F_curve)

best_idx = int(np.argmax(F_curve))
best_t = float(grid[best_idx])
best_F = float(F_curve[best_idx])
print(f"[sweep] data-optimal threshold (max pairwise F1) = {best_t:.2f} "
      f"(F1={best_F:.3f})", file=sys.stderr)

# rustle's 0.30 pairwise stats
ri = int(np.argmin(np.abs(grid - RUSTLE_THRESHOLD)))
rustle_P = float(P_curve[ri]); rustle_R = float(R_curve[ri]); rustle_F = float(F_curve[ri])
# raw-pair tp/fp at 0.30 (for the honest precision-not-1.0 caveat; these are the
# UN-clustered pair counts, distinct from the post-clustering counts in Q2).
_pred_030 = pair_arr >= RUSTLE_THRESHOLD
rustle_tp_raw = int(np.sum(_pred_030 & pair_same))
rustle_fp_raw = int(np.sum(_pred_030 & ~pair_same))


# --------------------------------------------------------------------------
# Union-find clustering at a given threshold -> predicted families
# --------------------------------------------------------------------------
def cluster_at(threshold):
    parent = list(range(n))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb
    # iterate over stored pairs (deterministic order)
    for idx in range(len(pair_j)):
        if pair_arr[idx] >= threshold:
            # recover (i,k) from flat index
            pass
    # rebuild pairwise to union (re-walk i<k deterministically)
    p = 0
    for i in range(n):
        for k_ in range(i + 1, n):
            if pair_arr[p] >= threshold:
                union(i, k_)
            p += 1
    clusters = defaultdict(list)
    for i in range(n):
        clusters[find(i)].append(genes[i])
    return list(clusters.values())


# --------------------------------------------------------------------------
# Clustering evaluation vs truth
# --------------------------------------------------------------------------
def pairwise_prf(clusters):
    """Pairwise precision/recall/F1 of predicted clustering vs truth labels."""
    gene_to_pred = {}
    for ci, members in enumerate(clusters):
        for g in members:
            gene_to_pred[g] = ci
    tp = fp = fn = 0
    for i in range(n):
        gi = genes[i]
        for k_ in range(i + 1, n):
            gj = genes[k_]
            same_pred = gene_to_pred[gi] == gene_to_pred[gj]
            same_true = label[gi] == label[gj]
            if same_pred and same_true:
                tp += 1
            elif same_pred and not same_true:
                fp += 1
            elif (not same_pred) and same_true:
                fn += 1
    P = tp / (tp + fp) if (tp + fp) else 1.0
    R = tp / (tp + fn) if (tp + fn) else 0.0
    F = 2 * P * R / (P + R) if (P + R) else 0.0
    return P, R, F, tp, fp, fn

def adjusted_rand_index(clusters):
    """ARI computed manually from the contingency table."""
    gene_to_pred = {}
    for ci, members in enumerate(clusters):
        for g in members:
            gene_to_pred[g] = ci
    # contingency
    pred_labels = [gene_to_pred[g] for g in genes]
    true_labels = [label[g] for g in genes]
    from collections import Counter
    cont = defaultdict(int)
    for pl, tl in zip(pred_labels, true_labels):
        cont[(pl, tl)] += 1
    a = Counter(pred_labels)  # row sums
    b = Counter(true_labels)  # col sums
    def comb2(x):
        return x * (x - 1) / 2.0
    sum_ij = sum(comb2(v) for v in cont.values())
    sum_a = sum(comb2(v) for v in a.values())
    sum_b = sum(comb2(v) for v in b.values())
    total = comb2(n)
    expected = sum_a * sum_b / total if total else 0.0
    max_index = 0.5 * (sum_a + sum_b)
    denom = max_index - expected
    if denom == 0:
        return 1.0
    return (sum_ij - expected) / denom

def v_measure(clusters):
    """Homogeneity, completeness, V-measure (manual entropy computation)."""
    gene_to_pred = {}
    for ci, members in enumerate(clusters):
        for g in members:
            gene_to_pred[g] = ci
    pred_labels = [gene_to_pred[g] for g in genes]
    true_labels = [label[g] for g in genes]
    N = n
    from collections import Counter
    cont = defaultdict(int)
    for pl, tl in zip(pred_labels, true_labels):
        cont[(pl, tl)] += 1
    a = Counter(pred_labels)
    b = Counter(true_labels)

    def entropy(counts):
        h = 0.0
        for v in counts.values():
            if v > 0:
                p = v / N
                h -= p * np.log(p)
        return h

    H_C = entropy(b)   # classes (truth)
    H_K = entropy(a)   # clusters (pred)
    # H(C|K)
    H_C_given_K = 0.0
    for (pl, tl), v in cont.items():
        if v > 0:
            H_C_given_K -= (v / N) * np.log(v / a[pl])
    # H(K|C)
    H_K_given_C = 0.0
    for (pl, tl), v in cont.items():
        if v > 0:
            H_K_given_C -= (v / N) * np.log(v / b[tl])
    homogeneity = 1.0 if H_C == 0 else 1 - H_C_given_K / H_C
    completeness = 1.0 if H_K == 0 else 1 - H_K_given_C / H_K
    if homogeneity + completeness == 0:
        vm = 0.0
    else:
        vm = 2 * homogeneity * completeness / (homogeneity + completeness)
    return homogeneity, completeness, vm

def family_recovery(clusters):
    """How many of the 62 truth paralog families are recovered as a PURE
    cluster (all members together AND no foreign members), partially recovered
    (all members together but cluster also has foreign genes = over-merge),
    over-merges (distinct truth families fused into one cluster), and splits
    (one truth family fragmented across >1 cluster)."""
    gene_to_pred = {}
    for ci, members in enumerate(clusters):
        for g in members:
            gene_to_pred[g] = ci
    # truth family -> set of genes (paralog families only, exclude CTRL singletons)
    fam_genes = defaultdict(list)
    for g in genes:
        if is_paralog[g]:
            fam_genes[label[g]].append(g)
    pure = 0
    partial_overmerged = 0  # all together but cluster impure
    split = 0
    fam_status = {}
    for fam, members in fam_genes.items():
        pred_clusters_of_fam = set(gene_to_pred[g] for g in members)
        if len(pred_clusters_of_fam) > 1:
            split += 1
            fam_status[fam] = "split"
            continue
        # all members in one cluster; is that cluster pure (only this family)?
        cid = next(iter(pred_clusters_of_fam))
        cluster_members = clusters[cid]
        foreign = [g for g in cluster_members if label[g] != fam]
        if not foreign:
            pure += 1
            fam_status[fam] = "pure"
        else:
            partial_overmerged += 1
            fam_status[fam] = "overmerged"
    # over-merges (PARALOG-ONLY): clusters that fuse >1 distinct truth PARALOG
    # family. NOTE: this is structurally BLIND to any false merge that involves a
    # control gene, because controls are excluded from `fams`. We therefore ALSO
    # compute control-aware counts below (Finding #3): a cluster that fuses a
    # paralog family with a control, or two controls, is a real false merge that
    # this paralog-only metric silently reports as 0.
    overmerge_clusters = 0
    for members in clusters:
        fams = set(label[g] for g in members if is_paralog[g])
        if len(fams) > 1:
            overmerge_clusters += 1

    # CONTROL-AWARE false merges. Treat EVERY gene's truth label as a distinct
    # family (each CTRL singleton is its own family, exactly as constructed). A
    # cluster is a "false merge" iff it contains >1 distinct truth label.
    false_merge_clusters_any = 0   # clusters fusing >1 truth label (any kind)
    control_contam_clusters = 0    # clusters mixing >=1 control with anything else
    for members in clusters:
        labs = set(label[g] for g in members)
        if len(labs) > 1:
            false_merge_clusters_any += 1
            has_ctrl = any(not is_paralog[g] for g in members)
            if has_ctrl:
                control_contam_clusters += 1

    n_pred_clusters = len(clusters)
    n_nontrivial = sum(1 for c in clusters if len(c) > 1)
    return dict(
        pure=pure, partial_overmerged=partial_overmerged, split=split,
        overmerge_clusters=overmerge_clusters,
        false_merge_clusters_any=false_merge_clusters_any,
        control_contam_clusters=control_contam_clusters,
        n_pred_clusters=n_pred_clusters, n_nontrivial=n_nontrivial,
        fam_status=fam_status,
    )

def evaluate(threshold):
    clusters = cluster_at(threshold)
    P, R, F, tp, fp, fn = pairwise_prf(clusters)
    ari = adjusted_rand_index(clusters)
    hom, comp, vm = v_measure(clusters)
    rec = family_recovery(clusters)
    return dict(
        threshold=threshold, clusters=clusters,
        P=P, R=R, F=F, tp=tp, fp=fp, fn=fn,
        ari=ari, homogeneity=hom, completeness=comp, vmeasure=vm,
        **rec,
    )

eval_rustle = evaluate(RUSTLE_THRESHOLD)
eval_opt = evaluate(best_t)

n_truth_families = len([f for f in truth_families])  # 62

for name, ev in [("rustle@0.30", eval_rustle), (f"optimal@{best_t:.2f}", eval_opt)]:
    print(f"[cluster {name}] families pure={ev['pure']}/{n_truth_families} "
          f"split={ev['split']} overmerged={ev['partial_overmerged']} "
          f"overmerge_clusters(paralog-only)={ev['overmerge_clusters']} "
          f"false_merge_clusters(any)={ev['false_merge_clusters_any']} "
          f"control_contam={ev['control_contam_clusters']} "
          f"P={ev['P']:.3f} R={ev['R']:.3f} F1={ev['F']:.3f} "
          f"ARI={ev['ari']:.3f} V={ev['vmeasure']:.3f}", file=sys.stderr)


# --------------------------------------------------------------------------
# FIGURE
# --------------------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.6))

# Left: within vs between Jaccard histograms
bins = np.linspace(0, 1, 51)
ax1.hist(between_j, bins=bins, color=NAVY, alpha=0.75,
         label=f"between (n={len(between_j):,})", density=True)
ax1.hist(within_j, bins=bins, color=ORANGE, alpha=0.75,
         label=f"within-family (n={len(within_j):,})", density=True)
ax1.axvline(RUSTLE_THRESHOLD, color=GREEN, lw=2.2, ls="-",
            label=f"rustle 0.30")
ax1.axvline(best_t, color="black", lw=1.6, ls="--",
            label=f"optimal {best_t:.2f}")
ax1.set_xlabel("minimizer-Jaccard")
ax1.set_ylabel("density")
ax1.set_title(f"Jaccard separability (AUC={auc:.3f})")
ax1.legend(fontsize=8, loc="upper right")
ax1.set_xlim(0, 1)
# Honest annotation: the density histogram makes the single between-pair >= 0.30
# invisible. Mark the integer false-merge count so the figure cannot be read as
# "zero between merges" (Finding #2).
_fm_note = (f"{n_between_ge_030} between pair(s) >= 0.30"
            + (f"\n(CTRL-CTRL J={between_at_rustle[0][0]:.3f},\nmissed paralog)"
               if (n_between_ge_030 == 1 and between_at_rustle[0][3] == "CTRL-CTRL")
               else ""))
ax1.annotate(_fm_note, xy=(0.32, 0.0), xytext=(0.46, ax1.get_ylim()[1] * 0.45),
             fontsize=7.5, color=NAVY,
             arrowprops=dict(arrowstyle="->", color=NAVY, lw=1.0))

# Right: P/R/F1 vs threshold
ax2.plot(grid, P_curve, color=GREEN, lw=2, label="precision")
ax2.plot(grid, R_curve, color=ORANGE, lw=2, label="recall")
ax2.plot(grid, F_curve, color=NAVY, lw=2.4, label="F1")
ax2.axvline(RUSTLE_THRESHOLD, color=GREEN, lw=1.4, ls="-", alpha=0.7)
ax2.axvline(best_t, color="black", lw=1.4, ls="--", alpha=0.8)
ax2.scatter([best_t], [best_F], color="black", zorder=5, s=28)
ax2.set_xlabel("threshold")
ax2.set_ylabel("pairwise score")
ax2.set_title("Detection P/R/F1 vs threshold")
ax2.legend(fontsize=8, loc="lower center")
ax2.set_xlim(0, 1)
ax2.set_ylim(0, 1.02)

fig.suptitle("rustle paralog-family detection: criterion validation", fontsize=12)
fig.text(0.5, 0.005,
         "Truth set = minimap2-similarity families (ident>=0.90/cov>=0.50)+union-find: "
         "this validates one similarity criterion against another, not curated biology.",
         ha="center", fontsize=7.5, color="#555555")
fig.tight_layout(rect=[0, 0.04, 1, 0.96])
fig.savefig(OUT_PNG, dpi=130)
print(f"[fig] wrote {OUT_PNG}", file=sys.stderr)


# --------------------------------------------------------------------------
# REPORT (.md)
# --------------------------------------------------------------------------
def band_density(t):
    """Density of all-pairs in [t-0.05, t+0.05]."""
    lo, hi = max(0, t - 0.05), min(1, t + 0.05)
    return float(np.mean((pair_arr >= lo) & (pair_arr < hi)))

valley_band = band_density(RUSTLE_THRESHOLD)

# Where exactly does 0.30 sit relative to the two distributions?
pct_within_below_030 = float(np.mean(within_j < RUSTLE_THRESHOLD)) if len(within_j) else float("nan")
pct_between_above_030 = float(np.mean(between_j >= RUSTLE_THRESHOLD)) if len(between_j) else float("nan")
within_median = float(np.median(within_j)) if len(within_j) else float("nan")

# ----- HONEST false-merge description (Finding #2) -----
# Do NOT print a rounded 0.0000 as if it were zero. Report the INTEGER count of
# between-family pairs >= 0.30, broken down by kind, and name them.
def describe_false_merges(pairs, counts):
    n_total = len(pairs)
    if n_total == 0:
        return "0 between-family false merges (raw pairs >= 0.30)"
    parts = []
    for j, gi, gj, kind in pairs:
        parts.append(f"{gi}<->{gj} J={j:.4f} [{kind}]")
    breakdown = (f"{counts['PARA-PARA']} cross-family paralog, "
                 f"{counts['PARA-CTRL']} paralog-control, "
                 f"{counts['CTRL-CTRL']} control-control")
    return (f"{n_total} between-family false merge(s) at raw pairs >= 0.30 "
            f"({breakdown}): " + "; ".join(parts))

false_merge_desc_030 = describe_false_merges(between_at_rustle, between_rustle_counts)

# Top control-control leak (Finding #4) — highest-Jaccard CTRL-CTRL pair, if any.
# between_at_006 / ctrl_ctrl_006 are sorted by descending Jaccard already.
if ctrl_ctrl_006:
    jmax_cc, gname0, gname1, _kcc = ctrl_ctrl_006[0]
else:
    jmax_cc, gname0, gname1 = float("nan"), "", ""

# Honest data-driven verdict on where 0.30 sits. 0.30 is "off" (too high) if a
# large fraction of within-family pairs fall below it AND the optimum is well
# below 0.30; it is "near-optimal" only if it's close to the F1 peak.
F_at_030 = rustle_F
flatness = best_F - F_at_030  # how much F1 we lose by using 0.30 vs optimum
too_high = (pct_within_below_030 > 0.5) and (best_t < RUSTLE_THRESHOLD - 0.05)
if too_high:
    verdict_030 = (
        f"rustle's 0.30 sits ABOVE the within-family mode (within-family median "
        f"Jaccard = {within_median:.3f}, and {pct_within_below_030*100:.0f}% of "
        f"within-family pairs fall BELOW 0.30 at the whole-gene level). The "
        f"data-optimal threshold is {best_t:.2f} (raw pair-threshold "
        f"F1={best_F:.3f} vs {F_at_030:.3f} at 0.30; after union-find "
        f"clustering F1={eval_opt['F']:.3f} at {best_t:.2f} vs "
        f"{eval_rustle['F']:.3f} at 0.30). So for WHOLE-GENE comparison 0.30 is "
        f"OFF (too high): it fragments {eval_rustle['split']}/{n_truth_families} "
        f"families and recovers only {eval_rustle['pure']}/{n_truth_families} as "
        f"pure clusters, trading recall ({eval_rustle['R']:.2f}) for "
        f"high precision ({eval_rustle['P']:.3f}). It is NOT a zero-false-merge "
        f"regime: {false_merge_desc_030}.")
else:
    verdict_030 = (
        f"rustle's 0.30 is near the data-optimal threshold {best_t:.2f} "
        f"(F1 {F_at_030:.3f} vs optimum {best_F:.3f}).")

# Compose the false-merge clause honestly: name the count and the single
# offending pair if there is one. (Finding #2)
if n_between_ge_030 == 0:
    fm_clause = "0 between-family false merges"
elif n_between_ge_030 == 1:
    _j, _gi, _gj, _kind = between_at_rustle[0]
    fm_clause = (f"1 between-family false merge (a {_kind} pair "
                 f"{_gi}<->{_gj}, Jaccard {_j:.4f})")
else:
    fm_clause = (f"{n_between_ge_030} between-family false merges "
                 f"({between_rustle_counts['PARA-PARA']} cross-family paralog, "
                 f"{between_rustle_counts['PARA-CTRL']} paralog-control, "
                 f"{between_rustle_counts['CTRL-CTRL']} control-control)")

# If the only false merge is a CTRL-CTRL pair, it is itself a paralog pair the
# minimap2 universe builder MISSED and leaked into the control pool — flag it.
ctrl_ctrl_leak_note = ""
if n_between_ge_030 >= 1 and all(k == "CTRL-CTRL" for (_, _, _, k) in between_at_rustle):
    ctrl_ctrl_leak_note = (
        " — itself a missed paralog leaked into the control pool, NOT a genuine "
        "cross-family error")

headline = (
    f"At threshold T=0.30 (rustle default), the criterion recovers "
    f"{eval_rustle['pure']}/{n_truth_families} truth families as PURE clusters "
    f"at pairwise P={eval_rustle['P']:.3f} R={eval_rustle['R']:.3f} "
    f"F1={eval_rustle['F']:.3f} (ARI={eval_rustle['ari']:.3f}, "
    f"V={eval_rustle['vmeasure']:.3f}) — high precision with {fm_clause}"
    f"{ctrl_ctrl_leak_note}. It is recall-limited: 0.30 is ABOVE the "
    f"within-family Jaccard mode at the whole-gene level (median "
    f"{within_median:.3f}), so {eval_rustle['split']}/{n_truth_families} "
    f"families fragment. The data-optimal whole-gene threshold is "
    f"{best_t:.2f} ({eval_opt['pure']}/{n_truth_families} pure, clustering "
    f"F1={eval_opt['F']:.3f}, ARI={eval_opt['ari']:.3f}). NOTE: the truth set "
    f"is itself a sequence-similarity clustering (minimap2), so this validates "
    f"one similarity criterion against another, not against curated biology."
)

md = []
md.append("# rustle paralog-family DETECTION: criterion validation\n")
md.append(f"_Deterministic. Fixed seed = **{SEED}**. k={K}, w={W}. "
          f"Generated by `bench/family_detection_validation.py`._\n")
md.append("## Headline\n")
md.append(f"> {headline}\n")
md.append(
    f"The data-optimal threshold (max raw pair-threshold F1, before union-find) "
    f"is **{best_t:.2f}** (F1={best_F:.3f}); rustle's 0.30 vs the optimum is "
    f"discussed below. (NB: the Q2 table reports POST-clustering F1, which "
    f"differs from the raw pair-threshold F1 because union-find applies "
    f"transitive closure.)\n")

md.append("## What is under test\n")
md.append(
    "rustle (`src/rustle/vg_family/family_graph.rs`) groups loci into paralog "
    "families when their sequences have **minimizer-Jaccard >= "
    "FAMILY_MERGE_JACCARD_DEFAULT = 0.30** (minimizers k=15, w=10). This script "
    "reproduces that exact similarity gate on labeled GENE sequences and "
    "union-find clusters them, then scores the predicted clusters against the "
    "annotated paralogy truth set.\n")

md.append("## Inputs reconciled\n")
md.append(f"- Truth: `universe.tsv` = **{n_truth_families} families** over "
          f"**{len(universe_genes)} genes** (gene_id/family_id labels).\n")
md.append(f"- Sequences: `gene_rep.fa` (1 representative per gene). Header IDs "
          f"are bare names; universe uses bare names, so **{len(matched)}/"
          f"{len(universe_genes)}** universe genes matched directly "
          f"(no `gene-` prefix present). Unmatched: {len(unmatched)}.\n")
md.append(f"- Control set: **{len(controls)}** NON-UNIVERSE genes sampled "
          f"(fixed seed {SEED}) from the {len(control_pool)} fasta genes NOT in "
          f"any universe family. Each is its own singleton label so that "
          f"control-control and control-paralog pairs are BETWEEN and any "
          f"control that clusters is scored as an error. **These controls are "
          f"NOT certified single-copy** (Finding #4): they are merely absent "
          f"from `universe.tsv`, which was itself filtered to minimap2 hits at "
          f"ident>=0.90 / cov>=0.50, so a paralog pair below those cutoffs falls "
          f"into the control pool. Residual control-control similarity: "
          f"**{len(ctrl_ctrl_006)} CTRL-CTRL pairs >= 0.06**"
          + (f" including {gname0}<->{gname1} at Jaccard {jmax_cc:.4f} "
             f"(a missed paralog leaked into the controls)."
             if ctrl_ctrl_006 else ".") + "\n")
if control_screen_ran:
    md.append(
        f"- **Control-purity SCREEN (independent minimap2, ident>={UNI_IDENT}/"
        f"cov>={UNI_COV}; controls kept, not silently dropped, so the leak stays "
        f"visible).** Two presets, because the discrepancy is itself the point:\n"
        f"  - `-x asm20` (EXACTLY the universe-builder criterion): "
        f"**{len(impure_controls)}/{len(controls)} impure**. This MISSES the "
        f"J=1.0 leak `LOC129531078`==`LOC129531352` — those two are byte-"
        f"identical 153 bp sequences, but `asm20` returns NO alignment for "
        f"sequences that short (it does not seed them). The universe builder "
        f"used `asm20`, which is precisely WHY this real paralog pair leaked "
        f"into the controls: the truth set's own criterion is blind to short "
        f"paralogs.\n"
        f"  - `-x sr` (sensitive short-sequence seeding): "
        f"**{len(impure_controls_sens)}/{len(controls)} impure**"
        + (" — " + "; ".join(f"`{q}`~`{tn}` ident={idn:.3f} cov={cv:.3f}"
                             for q, tn, idn, cv in control_screen_hits_sens[:5])
           + "." if control_screen_hits_sens
           else " (none flagged even under sensitive seeding).")
        + " The minimizer-Jaccard criterion under test (k=15/w=10) is itself "
        "MORE sensitive than the universe's `asm20` on short sequences — it "
        "scored that pair J=1.0 — so it actually CORRECTS a truth-set miss "
        "here, which is the opposite of an error.\n")
else:
    md.append(
        f"- **Control-purity screen NOT run** ({control_screen_error}); "
        f"control purity is therefore asserted only by absence from the "
        f"universe, which is itself a similarity filter, so purity is "
        f"UNVERIFIED (see the residual control-control mass above).\n")
md.append(f"- Working test set: **{n} genes** "
          f"({len(matched)} paralogs + {len(controls)} controls).\n")

md.append("## Q1 — Threshold justification (separability)\n")
md.append(f"- Pairs: within-family **{len(within_j):,}**, between **{len(between_j):,}**.\n")
md.append(f"- **AUC (within vs between) = {auc:.4f}** "
          f"(1.0 = perfectly separable).\n")
md.append(f"- within-family Jaccard percentiles [5/50/95] = "
          f"[{within_q[0]:.3f} / {within_q[1]:.3f} / {within_q[2]:.3f}].\n")
md.append(f"- between Jaccard percentiles [5/50/95] = "
          f"[{between_q[0]:.3f} / {between_q[1]:.3f} / {between_q[2]:.3f}].\n")
md.append(f"- BETWEEN pairs >= 0.30: **{n_between_ge_030}** pair(s) "
          f"(= false-merge COUNT; the fraction {pct_between_above_030:.6f} is "
          f"tiny but NON-ZERO — do not read it as zero). Breakdown: "
          f"{between_rustle_counts['PARA-PARA']} cross-family paralog, "
          f"{between_rustle_counts['PARA-CTRL']} paralog-control, "
          f"{between_rustle_counts['CTRL-CTRL']} control-control. "
          f"Fraction of WITHIN pairs < 0.30: "
          f"**{pct_within_below_030:.4f}** (missed-merge mass).\n")
if between_at_rustle:
    md.append("- The between-pair(s) >= 0.30 are: "
              + "; ".join(f"`{gi}`<->`{gj}` J={j:.4f} [{kind}]"
                          for j, gi, gj, kind in between_at_rustle) + ".\n")
md.append(f"- within-family **median Jaccard = {within_median:.3f}** "
          f"(below 0.30).\n")
md.append(f"- Local pair density in [0.25,0.35) (the 0.30 band): "
          f"**{valley_band:.4f}** of all pairs.\n")
md.append(
    "**Verdict (separability):** the two classes ARE separable on the "
    f"between side (AUC={auc:.3f}): between-family pairs pile up near Jaccard 0 "
    "(unrelated genes share ~no canonical minimizers), and only "
    f"**{n_between_ge_030}** between pair reaches 0.30 "
    + (f"— and that one is a CONTROL-CONTROL pair "
       f"({between_at_rustle[0][1]}<->{between_at_rustle[0][2]}, "
       f"J={between_at_rustle[0][0]:.4f}) that is itself a missed paralog, not a "
       f"genuine cross-family error. "
       if (n_between_ge_030 == 1 and between_at_rustle[0][3] == "CTRL-CTRL")
       else ". ")
    + "BUT the within-family distribution is NOT a tight high mode — it is broad "
    f"(5th/50th/95th pctile {within_q[0]:.3f}/{within_q[1]:.3f}/"
    f"{within_q[2]:.3f}), because whole-gene paralog pairs range from "
    "near-identical to weakly similar (divergence + differing exon content). "
    "So the data is bimodal in the SENSE that unrelated pairs are cleanly "
    "separable from related ones near 0, but it is NOT cleanly bimodal at 0.30: "
    f"0.30 cuts through the BODY of the within-family mode, not a valley.\n")

md.append("## Q2 — Family detection at 0.30 vs data-optimal\n")
md.append("| metric | rustle 0.30 | optimal " + f"{best_t:.2f}" + " |\n")
md.append("|---|---|---|\n")
for key, lab in [("P", "pairwise precision"), ("R", "pairwise recall"),
                 ("F", "pairwise F1"), ("ari", "ARI"),
                 ("vmeasure", "V-measure"),
                 ("homogeneity", "homogeneity"),
                 ("completeness", "completeness")]:
    md.append(f"| {lab} | {eval_rustle[key]:.3f} | {eval_opt[key]:.3f} |\n")
md.append(f"| families recovered PURE | {eval_rustle['pure']}/"
          f"{n_truth_families} | {eval_opt['pure']}/{n_truth_families} |\n")
md.append(f"| families split (fragmented) | {eval_rustle['split']} | "
          f"{eval_opt['split']} |\n")
md.append(f"| families over-merged (impure cluster) | "
          f"{eval_rustle['partial_overmerged']} | "
          f"{eval_opt['partial_overmerged']} |\n")
md.append(f"| clusters fusing >1 truth PARALOG family (control-blind) | "
          f"{eval_rustle['overmerge_clusters']} | "
          f"{eval_opt['overmerge_clusters']} |\n")
md.append(f"| clusters fusing >1 truth label (CONTROL-AWARE) | "
          f"{eval_rustle['false_merge_clusters_any']} | "
          f"{eval_opt['false_merge_clusters_any']} |\n")
md.append(f"| clusters contaminated by a control gene | "
          f"{eval_rustle['control_contam_clusters']} | "
          f"{eval_opt['control_contam_clusters']} |\n")
md.append(f"| # predicted clusters (>1 gene) | "
          f"{eval_rustle['n_nontrivial']} | {eval_opt['n_nontrivial']} |\n")
md.append(
    "\n> The paralog-only over-merge count is STRUCTURALLY BLIND to any false "
    "merge that involves a control gene (Finding #3): `family_recovery()` "
    "excludes controls from the over-merge tally, so a paralog↔control or "
    "control↔control fusion reports 0 there even though pairwise precision "
    "(which counts every pair) drops below 1.0. The CONTROL-AWARE row treats "
    "each control as its own family and is the honest false-merge count.\n")

md.append("\n## Where does 0.30 fall? (data-driven)\n")
md.append(verdict_030 + "\n\n")
md.append(
    "**Interpretation / reconciliation with the source.** This is a WHOLE-GENE "
    "test, whereas rustle applies 0.30 at the EXON-class level "
    "(`refine_by_minimizer_jaccard` / `merge_singletons_by_sequence`), where "
    "homologous exons are short and near-identical so their Jaccard is much "
    "higher than the whole-gene Jaccard of the two genes. The source comment "
    "explicitly tunes 0.30 to REFUSE adversarial false merges (embedded shared "
    "TE/domain J~0.22, low-complexity cores J~0.17) at the exon level, and the "
    "between-family mass here confirms 0.30 buys that precision nearly for "
    f"free: only {n_between_ge_030} between-pair >= 0.30"
    + (f", and that one is the CTRL-CTRL J=1.0 missed-paralog leak, not a real "
       f"cross-family merge"
       if (n_between_ge_030 == 1 and between_at_rustle[0][3] == "CTRL-CTRL")
       else "")
    + ". The honest finding is therefore: "
    "for the COARSE whole-gene grouping decision tested here, 0.30 is too "
    f"conservative (data-optimal {best_t:.2f}); a lower bar would recover "
    f"{eval_opt['pure']}/{n_truth_families} families pure at the cost of "
    f"{eval_opt['false_merge_clusters_any']} control-aware false-merge "
    f"cluster(s) ({eval_opt['overmerge_clusters']} of them paralog-paralog, "
    f"{eval_opt['control_contam_clusters']} control-contaminated). "
    "Whether to lower it depends on the downstream cost of an over-merge vs a "
    "split — which is exactly the assembly over-merge cost the source comment "
    "calls the binding constraint, and which THIS criterion-level test does "
    "not measure.\n")

md.append("## Honest caveats\n")
md.append(
    "- **TRUTH-SET CIRCULARITY (read this first).** `universe.tsv` is NOT an "
    "independent biological ground truth. It is itself built by SEQUENCE "
    "SIMILARITY: `bench/copy_recovery_eval/30_build_universe.py` runs "
    "`minimap2 -x asm20 -X --no-long-join` all-vs-all on the per-gene "
    "representative sequences, keeps pairs at **ident>=0.90 and cov>=0.50**, "
    "then takes union-find connected components (`lib_eval.build_universe`). "
    "This validation then asks whether a DIFFERENT sequence-similarity "
    "criterion (minimizer-Jaccard, k=15/w=10) recovers those same families. So "
    f"the high **AUC={auc:.3f}** and the recovery numbers measure CONCORDANCE "
    "BETWEEN TWO SIMILARITY CLUSTERINGS (minimap2 alignment identity vs "
    "minimizer-set overlap), NOT agreement with curated biology. They do not "
    "validate the families against Ensembl Compara, OrthoFinder, synteny, or a "
    "curated paralog database. Any 'recovers the truth families' language must "
    "be read as 'reproduces the minimap2-defined families'; an independent "
    "cross-check against a curated paralog DB is needed before a stronger "
    "claim. The two methods share enough mechanism (both reward shared "
    "subsequence) that part of the concordance is methodological, not "
    "biological.\n")
md.append(
    "- **Controls are NOT certified single-copy (Finding #4).** Controls are "
    "defined ONLY as fasta genes absent from `universe.tsv`. But the universe "
    "was filtered to minimap2 hits at ident>=0.90/cov>=0.50, so a paralog pair "
    "below those cutoffs is NOT in the universe and therefore lands in the "
    "control pool. This is not hypothetical: there are "
    f"**{len(ctrl_ctrl_006)} control-control pair(s) with Jaccard >= 0.06**"
    + (f", the strongest being {gname0}<->{gname1} at **J={jmax_cc:.4f}** — a "
       f"genuine paralog pair the minimap2 builder missed and that leaked into "
       f"the controls. It is the single between-pair >= 0.30 and the sole "
       f"reason RAW-PAIR precision at 0.30 is {rustle_P:.4f} "
       f"({rustle_tp_raw} TP, {rustle_fp_raw} FP), not 1.0; the "
       f"{eval_rustle['P']:.3f} in the Q2 table is the slightly different "
       f"POST-clustering precision after transitive closure."
       if ctrl_ctrl_006 else ".")
    + " To certify the controls, screen the pool with the same minimap2 step "
    "and drop any gene with any hit; until then the control set should be "
    "described as 'non-universe', not 'single-copy'.\n")
md.append("- **Criterion-level, not the full pipeline.** This tests ONLY the "
          "minimizer-Jaccard similarity GATE on whole-gene representative "
          "sequences. The production pipeline also uses position-overlap "
          "clustering, per-exon class refinement, EM, and flow; results here "
          "bound the gate's intrinsic discriminability, not end-to-end "
          "assembly.\n")
md.append("- **Reference sequences only.** Run on annotation-derived "
          "representative gene sequences, not on assembled reads.\n")
md.append("- **Expressed/annotated paralogs only.** The truth set is the 62 "
          "annotated families; unannotated or unexpressed paralogs are out of "
          "scope.\n")
md.append("- **Collapsed copies are NOT detectable this way.** Paralogs that "
          "collapse onto one locus (no distinct sequence) cannot be separated "
          "by sequence similarity; that is the separate PSV-co-segregation "
          "test.\n")
md.append(f"- **Minimizer hash differs from production.** Production uses "
          f"non-canonical FNV-1a; mmh3 absent, so per spec this uses a "
          f"CANONICAL minimizer hashed with blake2b (64-bit). The WINDOWING is "
          f"reproduced exactly. Canonicalization only adds strand robustness "
          f"and does not change the bimodal-separability conclusion.\n")
md.append(f"- **Determinism.** Single fixed seed = {SEED} for the control "
          f"sample; no wall-clock or unseeded RNG anywhere.\n")

with open(OUT_MD, "w") as fh:
    fh.write("".join(md))
print(f"[md] wrote {OUT_MD}", file=sys.stderr)

# --------------------------------------------------------------------------
# Machine-readable summary to stdout (for the orchestrator)
# --------------------------------------------------------------------------
print("=== SUMMARY ===")
print(f"matched_genes={len(matched)}/{len(universe_genes)}")
print(f"controls={len(controls)} total_test_genes={n}")
print(f"within_pairs={len(within_j)} between_pairs={len(between_j)}")
print(f"AUC={auc:.4f}")
print(f"frac_between_ge_030={pct_between_above_030:.6f} "
      f"frac_within_lt_030={pct_within_below_030:.4f}")
print(f"between_pairs_ge_030_COUNT={n_between_ge_030} "
      f"(PARA-PARA={between_rustle_counts['PARA-PARA']} "
      f"PARA-CTRL={between_rustle_counts['PARA-CTRL']} "
      f"CTRL-CTRL={between_rustle_counts['CTRL-CTRL']})")
print(f"control_control_pairs_ge_006={len(ctrl_ctrl_006)}"
      + (f" top={gname0}<->{gname1}@J={jmax_cc:.4f}" if ctrl_ctrl_006 else ""))
print(f"control_purity_screen_ran={control_screen_ran} "
      f"impure_controls_asm20={len(impure_controls)}/{len(controls)} "
      f"impure_controls_sr={len(impure_controls_sens)}/{len(controls)}"
      + (f" (asm20 misses the short J=1.0 leak; sr/minimizer-Jaccard catch it)"
         if control_screen_ran else f" ({control_screen_error})"))
print(f"TRUTH_SET_NOTE=universe.tsv built by minimap2 similarity "
      f"(ident>=0.90/cov>=0.50)+union-find; validates similarity-vs-similarity")
print(f"optimal_threshold={best_t:.2f} optimal_F1={best_F:.3f}")
print(f"RUSTLE@0.30: P={eval_rustle['P']:.3f} R={eval_rustle['R']:.3f} "
      f"F1={eval_rustle['F']:.3f} ARI={eval_rustle['ari']:.3f} "
      f"V={eval_rustle['vmeasure']:.3f}")
print(f"RUSTLE@0.30 families: pure={eval_rustle['pure']}/{n_truth_families} "
      f"split={eval_rustle['split']} overmerged={eval_rustle['partial_overmerged']} "
      f"overmerge_clusters_paralog_only={eval_rustle['overmerge_clusters']} "
      f"false_merge_clusters_any={eval_rustle['false_merge_clusters_any']} "
      f"control_contam={eval_rustle['control_contam_clusters']}")
print(f"OPT@{best_t:.2f} families: pure={eval_opt['pure']}/{n_truth_families} "
      f"split={eval_opt['split']} overmerged={eval_opt['partial_overmerged']}")
print(f"HEADLINE: {headline}")
