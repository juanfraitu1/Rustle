#!/usr/bin/env python3
"""Score a catalog against Soto 2025's families, per the pre-registered rule of ledger §6bw.

Truth: `bench/soto/80_fams.chr.bed` — 83 families / 362 member genes, CHM13 v2.0, the same build
as `A119b.t2t.bam`. Prediction: `cat.copies.tsv` from one or two arms.

⚠ **This is RE-IMPLEMENTATION CONCORDANCE, never "replication" and never "validation"** — three of
the four legs of any Soto comparison consume Soto's own files (2026-08-02 audit).
⚠ **Precision is UNDERSTATED by construction**: Soto's set is CAT-bounded, so a real copy CAT missed
scores as a false positive.
⚠ **Soto = SD98, i.e. >=98% identity.** It adjudicates the HIGH-identity stratum only; it cannot
speak to the 0.80-0.90 band that is 86.31% of this catalog. Output is stratified for that reason.

Three units, all reported, because quoting one is how this metric misleads (§6bw):
  DETECTION  — is a Soto member emitted as a copy at all (any-overlap AND >=50% of the gene)?
  PAIR       — over gene pairs: Soto puts A,B together; do we? The unit the edge rule asserts.
  FAMILY     — exact set match; reported for completeness, expected low.

Usage: soto_adjudicate.py <bed> <arm_dir> [<arm_dir2>] [--reads <bam>]
"""
import csv
import sys
import os
import math
import collections
import itertools


def wilson(k, n, z=1.96):
    if n == 0:
        return (0.0, 1.0)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def load_truth(bed):
    genes, fam = {}, {}
    for line in open(bed):
        f = line.rstrip("\n").split("\t")
        if len(f) < 4:
            continue
        name = f[3]
        g, fid = name.split("|") if "|" in name else (name, "NA")
        genes[name] = (f[0], int(f[1]), int(f[2]))
        fam[name] = fid
    return genes, fam


def assign(copies, genes, min_frac):
    """copy -> the Soto gene it best overlaps, requiring >= min_frac of the GENE."""
    out = {}
    for i, r in enumerate(copies):
        c, s, e = r["chrom"], int(r["start"]), int(r["end"])
        best = None
        for name, (gc, gs, ge) in genes.items():
            if gc != c:
                continue
            ov = min(ge, e) - max(gs, s)
            if ov <= 0:
                continue
            frac = ov / max(1, ge - gs)
            if frac >= min_frac and (best is None or frac > best[0]):
                best = (frac, name)
        if best:
            out[i] = best[1]
    return out


def score(arm_dir, genes, fam, min_frac, label):
    copies = list(csv.DictReader(open(os.path.join(arm_dir, "cat.copies.tsv")), delimiter="\t"))
    a = assign(copies, genes, min_frac)
    detected = set(a.values())
    print(f"\n--- {label}   (gene-overlap floor {min_frac:.0%})")
    print(f"  copies emitted: {len(copies)}   mapped to a Soto gene: {len(a)}   "
          f"UNLABELLED: {len(copies)-len(a)}")
    lo, hi = wilson(len(detected), len(genes))
    print(f"  DETECTION: {len(detected)}/{len(genes)} = {len(detected)/len(genes):.4f} "
          f"Wilson95 [{lo:.4f},{hi:.4f}]")

    # ---- PAIR unit
    pred = collections.defaultdict(set)          # our family -> set of Soto genes
    for i, r in enumerate(copies):
        if i in a:
            pred[r["family_id"]].add(a[i])
    our_pairs = set()
    for fid, gs in pred.items():
        for x, y in itertools.combinations(sorted(gs), 2):
            our_pairs.add((x, y))
    truth_pairs = set()
    byfam = collections.defaultdict(list)
    for g, f in fam.items():
        byfam[f].append(g)
    for f, gs in byfam.items():
        for x, y in itertools.combinations(sorted(gs), 2):
            truth_pairs.add((x, y))
    truth_det = {(x, y) for (x, y) in truth_pairs if x in detected and y in detected}
    tp = len(our_pairs & truth_pairs)
    fp = len(our_pairs - truth_pairs)
    fn_det = len(truth_det - our_pairs)
    fn_all = len(truth_pairs - our_pairs)
    P = tp / (tp + fp) if tp + fp else float("nan")
    R_det = tp / (tp + fn_det) if tp + fn_det else float("nan")
    R_all = tp / (tp + fn_all) if tp + fn_all else float("nan")
    F = 2 * P * R_det / (P + R_det) if P + R_det else float("nan")
    print(f"  PAIR  asserted={len(our_pairs)}  TP={tp}  FP={fp}")
    print(f"        precision              = {P:.4f}  ⚠ UNDERSTATED (Soto is CAT-bounded)")
    print(f"        recall | both detected = {R_det:.4f}   ({tp}/{tp+fn_det})")
    print(f"        recall | ALL Soto pairs= {R_all:.4f}   ({tp}/{tp+fn_all})  ⚠ bounded by EXPRESSION")
    print(f"        F1 (detected denom)    = {F:.4f}")

    # ---- FAMILY unit
    truth_sets = {frozenset(g for g in gs if g in detected) for gs in byfam.values()}
    truth_sets = {t for t in truth_sets if len(t) >= 2}
    pred_sets = {frozenset(v) for v in pred.values() if len(v) >= 2}
    exact = len(truth_sets & pred_sets)
    print(f"  FAMILY exact set match: {exact}/{len(truth_sets)} "
          f"(multi-member truth families that are fully detected)")
    ident = edge_identity(arm_dir, copies, a)
    if ident:
        band_table(label, our_pairs, truth_pairs, ident)
    return dict(copies=copies, assign=a, detected=detected, our_pairs=our_pairs,
                truth_pairs=truth_pairs, tp=tp, fp=fp, P=P, R_det=R_det, ident=ident)


BANDS = [(0.60, 0.70), (0.70, 0.80), (0.80, 0.90), (0.90, 1.001)]


def edge_identity(arm_dir, copies, a):
    """Map each ASSERTED gene pair to the identity of the E_r edge that links the two copies.

    ⭐ Required by §6bw: Soto is SD98 (>=98% identity), so it can only adjudicate the HIGH band.
    Reporting P/R without this stratification would let a high-band result be read as covering the
    0.80-0.90 band, which is 86.31% of the catalog and which Soto cannot see.
    """
    ef = os.path.join(arm_dir, "dump", "e.edges.tsv")
    if not os.path.exists(ef):
        return {}
    # copy interval -> Soto gene
    g_of = {}
    for i, r in enumerate(copies):
        if i in a:
            g_of[(r["chrom"], int(r["start"]), int(r["end"]))] = a[i]
    out = {}
    for r in csv.DictReader(open(ef), delimiter="\t"):
        ki = (r["chrom_i"], int(r["start_i"]), int(r["end_i"]))
        kj = (r["chrom_j"], int(r["start_j"]), int(r["end_j"]))
        gi, gj = g_of.get(ki), g_of.get(kj)
        if not gi or not gj or gi == gj:
            continue
        pair = tuple(sorted((gi, gj)))
        ident = float(r["identity"])
        if pair not in out or ident > out[pair]:
            out[pair] = ident
    return out


def band_table(label, our_pairs, truth_pairs, ident):
    print(f"\n  === {label}: asserted pairs by E_r identity band (⚠ Soto = SD98, >=98% id) ===")
    print(f"  {'band':>14s} {'asserted':>9s} {'TP':>5s} {'FP':>5s} {'precision':>10s}")
    seen = 0
    for lo, hi in BANDS:
        sub = [p for p in our_pairs if p in ident and lo <= ident[p] < hi]
        seen += len(sub)
        tp = sum(1 for p in sub if p in truth_pairs)
        fp = len(sub) - tp
        pr = tp / len(sub) if sub else float("nan")
        print(f"  [{lo:.2f},{hi:.2f}) {len(sub):>9} {tp:>5} {fp:>5} {pr:>10.4f}")
    nomap = len(our_pairs) - seen
    print(f"  {'no edge found':>14s} {nomap:>9}   ⚠ pair asserted transitively, not by a direct edge")


def main():
    bed, arms = sys.argv[1], [x for x in sys.argv[2:] if not x.startswith("--")]
    genes, fam = load_truth(bed)
    nfam = len(set(fam.values()))
    single = sum(1 for f, n in collections.Counter(fam.values()).items() if n == 1)
    print(f"TRUTH: {len(genes)} member genes, {nfam} families "
          f"({single} singletons — excluded from the pair unit by construction)")

    res = {}
    for frac, tag in ((0.0001, "any-overlap"), (0.50, ">=50% of the gene")):
        print(f"\n================ {tag} ================")
        for d in arms:
            res[(d, frac)] = score(d, genes, fam, frac, os.path.basename(d.rstrip("/")))

    if len(arms) == 2:
        A, B = res[(arms[0], 0.50)], res[(arms[1], 0.50)]
        print(f"\n=== ⚖️ THE ADJUDICATION (>=50% floor): what the ON arm did to LABELLED copies ===")
        ka = {(r["chrom"], r["start"], r["end"]) for i, r in enumerate(A["copies"]) if i in A["assign"]}
        kb = {(r["chrom"], r["start"], r["end"]) for i, r in enumerate(B["copies"]) if i in B["assign"]}
        lost_lab = ka - kb
        print(f"  Soto-LABELLED copies: {len(ka)} -> {len(kb)}   LOST: {len(lost_lab)}")
        print(f"  Soto members detected: {len(A['detected'])} -> {len(B['detected'])}   "
              f"LOST members: {sorted(A['detected'] - B['detected'])}")
        print(f"  pair precision {A['P']:.4f} -> {B['P']:.4f}   "
              f"pair recall(det) {A['R_det']:.4f} -> {B['R_det']:.4f}")
        print(f"  asserted pairs {len(A['our_pairs'])} -> {len(B['our_pairs'])}   "
              f"TP {A['tp']} -> {B['tp']}   FP {A['fp']} -> {B['fp']}")


if __name__ == "__main__":
    main()
