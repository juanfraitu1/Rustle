#!/usr/bin/env python3
"""flank_homology_differentiator.py — test whether upstream/downstream flanking
homology separates real paralog edges from domain-share over-merges.

HYPOTHESIS.  Real paralog families arise from whole-locus duplications, so the
sequence *outside* the transcribed body (upstream / downstream) should also be
homologous.  Domain-share over-merges share only one exon/domain; their flanks
are unrelated.  This is orthogonal to the five settled exclusion axes
(nucleotide/protein/TE/topology/architecture) because it looks at sequence
*adjacent* to, rather than inside, the transcript.

METHOD.  On the same 5,571 refined-E_r direct edges used by rna_levers_explore,
extract 2 kb upstream and 2 kb downstream of each de-novo locus (oriented
5'->3' relative to transcription), align same-orientation flanks with minimap2,
and compute per-edge flank-homology features.  Evaluate vs the DNA label
in_dna_loose (real edge) and the genuine-over-merge residual, using the same
component folds stored in rna_levers_explore.tsv.

OUTPUTS.  bench/flank_homology_edge.tsv  (per-edge features)
          bench/flank_homology_edge.json (summary metrics)

Deterministic: PYTHONHASHSEED=0, fixed seeds, sorted iteration.
"""
import os
import sys

if os.environ.get("PYTHONHASHSEED") != "0":
    os.environ["PYTHONHASHSEED"] = "0"
    os.execv(sys.executable, [sys.executable] + sys.argv)

import csv
import json
import math
import subprocess
import tempfile
from collections import defaultdict

import numpy as np
import pysam
import scipy.stats as ss

BENCH = os.path.dirname(os.path.abspath(__file__))
SCRATCH = "/home/juanfra/winloci_scratch"
sys.path.insert(0, BENCH)
from rna_only_edge_oracle import auc_mw, standardize, logistic_fit  # noqa: E402

META_TSV = os.path.join(SCRATCH, "denovo_transcripts.meta.tsv")
FASTA = os.path.join(SCRATCH, "GGO.fasta")
LEVER_TSV = os.path.join(BENCH, "rna_levers_explore.tsv")
UNIV_TSV = os.path.join(BENCH, "ri_sharedlen_universal.tsv")

OUT_TSV = os.path.join(BENCH, "flank_homology_edge.tsv")
OUT_JSON = os.path.join(BENCH, "flank_homology_edge.json")

MINIMAP_PRESET = "asm20"
NA = "NA"
FLANK = 2000
MASK_REPEATS = False

COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G",
              "a": "t", "t": "a", "g": "c", "c": "g", "N": "N", "n": "n"}


def revcomp(s):
    return "".join(COMPLEMENT.get(b, "N") for b in reversed(s))


def load_meta():
    d = {}
    with open(META_TSV) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            d[row["id"]] = dict(
                chrom=row["chrom"],
                start=int(row["start"]),
                end=int(row["end"]),
                strand=row["strand"],
                n_exon=int(row["n_exon"]),
                n_reads=int(row["n_reads"]),
            )
    return d


def load_dn_map():
    """gene-pair frozenset -> (dn_a, dn_b), preserving file order."""
    d = {}
    with open(UNIV_TSV) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            k = frozenset((row["gene_a"], row["gene_b"]))
            d[k] = (row["dn_a"], row["dn_b"])
    return d


def load_lever():
    rows = []
    with open(LEVER_TSV) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            rows.append(row)
    return rows


def mask_lowercase(seq):
    """Convert soft-masked (lowercase) bases to N; preserves length."""
    return "".join("N" if c.islower() else c for c in seq)


def flank_seq(fa, meta, dn_id, which):
    """Return oriented 5'->3' flank sequence for DN locus.

    which = 'up' (upstream) or 'down' (downstream).
    Handles contig boundaries by clipping; returns empty string if flank is 0 bp.
    """
    m = meta.get(dn_id)
    if m is None:
        return ""
    chrom, start, end, strand = m["chrom"], m["start"], m["end"], m["strand"]
    # reference length for clipping
    reflen = fa.get_reference_length(chrom)
    if strand == "+":
        if which == "up":
            s = max(0, start - FLANK)
            e = start
            seq = fa.fetch(chrom, s, e)
        else:  # down
            s = end
            e = min(reflen, end + FLANK)
            seq = fa.fetch(chrom, s, e)
    else:  # '-' strand; extract the genomic complement and reverse it
        if which == "up":
            # upstream in transcription orientation = genomic [end .. end+FLANK), RC
            s = end
            e = min(reflen, end + FLANK)
            seq = revcomp(fa.fetch(chrom, s, e))
        else:  # down
            # downstream in transcription orientation = genomic [start-FLANK .. start), RC
            s = max(0, start - FLANK)
            e = start
            seq = revcomp(fa.fetch(chrom, s, e))
    if MASK_REPEATS:
        seq = mask_lowercase(seq)
    return seq


def write_fasta(path, seqs):
    with open(path, "w") as fh:
        for name, seq in sorted(seqs.items()):
            fh.write(f">{name}\n{seq}\n")


def run_minimap2(target_fa, query_fa, out_paf):
    cmd = [
        "/home/juanfra/miniforge3/bin/minimap2",
        "-x", MINIMAP_PRESET,
        "-c",          # output CIGAR / identity tags
        "-X",          # all-vs-all
        "-t", "4",
        target_fa, query_fa,
    ]
    with open(out_paf, "w") as fh:
        subprocess.run(cmd, stdout=fh, check=True)


def parse_paf(paf_path):
    """Parse all-vs-all PAF; return dict (qname,tname) -> best record.

    For each unordered pair we keep the alignment maximizing nmatch.
    Returns a dict keyed by frozenset({qname,tname}) mapping to a dict with:
        qlen, tlen, qspan, tspan, nmatch, alen, dv, id
    """
    best = {}
    with open(paf_path) as fh:
        for line in fh:
            if line.startswith("@") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            qn, qlen, qs, qe = f[0], int(f[1]), int(f[2]), int(f[3])
            tn, tlen, ts, te = f[5], int(f[6]), int(f[7]), int(f[8])
            nmatch, alen = int(f[9]), int(f[10])
            if qn == tn:
                continue
            # parse tags for dv (sequence divergence) if present
            dv = None
            for tag in f[12:]:
                if tag.startswith("dv:f:"):
                    dv = float(tag.split(":")[-1])
                    break
            # approximate identity if dv missing
            if dv is None:
                id_approx = nmatch / alen if alen > 0 else 0.0
            else:
                id_approx = max(0.0, 1.0 - dv)
            k = frozenset((qn, tn))
            # keep the hit with the most matching bases
            cur = best.get(k)
            if cur is None or nmatch > cur["nmatch"]:
                best[k] = dict(
                    qlen=qlen, tlen=tlen,
                    qspan=qe - qs, tspan=te - ts,
                    nmatch=nmatch, alen=alen,
                    dv=dv, id=id_approx,
                )
    return best


def flank_features(up_best, down_best, dn_a, dn_b):
    """Compute per-edge flank features from upstream/downstream PAF dicts."""
    def get(paf, a, b):
        rec = paf.get(frozenset((a, b)))
        if rec is None:
            return None, None
        # reciprocal alignment fraction: min of oriented spans / lengths
        af = min(rec["qspan"] / rec["qlen"], rec["tspan"] / rec["tlen"]) if rec["qlen"] and rec["tlen"] else 0.0
        af = min(1.0, af)
        return af, rec["id"]

    up_af, up_id = get(up_best, f"up:{dn_a}", f"up:{dn_b}")
    down_af, down_id = get(down_best, f"down:{dn_a}", f"down:{dn_b}")

    afs = [x for x in (up_af, down_af) if x is not None]
    ids = [x for x in (up_id, down_id) if x is not None]

    return dict(
        up_aln_frac=fmt(up_af), up_id=fmt(up_id),
        down_aln_frac=fmt(down_af), down_id=fmt(down_id),
        max_flank_aln_frac=fmt(max(afs) if afs else None),
        min_flank_aln_frac=fmt(min(afs) if afs else None),
        mean_flank_aln_frac=fmt(sum(afs) / len(afs) if afs else None),
        best_flank_id=fmt(ids[afs.index(max(afs))] if afs and ids else None),
    )


def fmt(x):
    if x is None or (isinstance(x, float) and x != x):
        return NA
    if isinstance(x, float):
        return f"{x:.6f}"
    return str(x)


def parse_float(x):
    if x is None or x == NA or x == "":
        return float("nan")
    return float(x)


def oof_logistic_auc(feat_cols, y, folds):
    """Out-of-fold logistic AUC.  Standardize and fit inside each training fold."""
    feat_cols = np.asarray(feat_cols, dtype=float)
    y = np.asarray(y, dtype=int)
    folds = np.asarray(folds, dtype=int)
    pred = np.full(len(y), np.nan)
    for k in sorted(set(folds)):
        train = folds != k
        test = folds == k
        Xtr, mu, sd = standardize(feat_cols[train])
        w = logistic_fit(Xtr, y[train])
        Xte = (feat_cols[test] - mu) / sd
        Xte = np.hstack([np.ones((Xte.shape[0], 1)), Xte])
        pred[test] = 1.0 / (1.0 + np.exp(-np.clip(Xte @ w, -30, 30)))
    pos = pred[y == 1]
    neg = pred[y == 0]
    a, _, _ = auc_mw(list(pos), list(neg))
    return a


def evaluate(rows):
    """Compute univariate and combined AUCs.  Returns metrics dict."""
    y = np.array([1 if r["in_dna_loose"] == "1" else 0 for r in rows])
    folds = np.array([int(r["fold"]) for r in rows])

    feats = ["up_aln_frac", "down_aln_frac", "max_flank_aln_frac",
             "min_flank_aln_frac", "mean_flank_aln_frac", "best_flank_id"]

    def col(name):
        return np.array([parse_float(r[name]) for r in rows])

    # univariate pooled AUCs
    univariate = {}
    for name in feats:
        x = col(name)
        pos = x[y == 1]
        neg = x[y == 0]
        # drop nan
        pos = pos[~np.isnan(pos)]
        neg = neg[~np.isnan(neg)]
        a, np_, nn = auc_mw(list(pos), list(neg))
        univariate[name] = dict(auc=round(a, 6), n_pos=int(np_), n_neg=int(nn))

    # OOF logistic combinations.  Missing flank homology = no detectable alignment,
    # which we model as 0 for the additive combination (the TSV keeps NA).
    core = col("core_recip")
    aln = col("aln_frac")
    baseline_X = np.column_stack([core, aln])
    baseline_auc = oof_logistic_auc(baseline_X, y, folds)

    combos = {}
    for name in ["max_flank_aln_frac", "mean_flank_aln_frac", "min_flank_aln_frac"]:
        x = col(name)
        x = np.where(np.isnan(x), 0.0, x)
        X = np.column_stack([core, aln, x])
        a = oof_logistic_auc(X, y, folds)
        combos[name] = dict(auc=round(a, 6), lift=round(a - baseline_auc, 6))

    # residual-crack: edges touching the irreducible roster, kept by aln_frac
    residual = []
    for r in rows:
        if r["residual"] == "1":
            residual.append(r)
    res_metrics = {}
    if residual:
        ry = [1 if r["in_dna_loose"] == "1" else 0 for r in residual]
        for name in ["max_flank_aln_frac", "mean_flank_aln_frac"]:
            vals = [parse_float(r[name]) for r in residual]
            pos = [v for v, lab in zip(vals, ry) if lab == 1 and v == v]
            neg = [v for v, lab in zip(vals, ry) if lab == 0 and v == v]
            a, np_, nn = auc_mw(pos, neg)
            res_metrics[name] = dict(auc=round(a, 6), n_pos=int(np_), n_neg=int(nn))

    return dict(
        n_edges=len(rows),
        n_pos=int(y.sum()),
        n_neg=int((1 - y).sum()),
        baseline_core_aln_oof_auc=round(baseline_auc, 6),
        univariate=univariate,
        oof_combinations=combos,
        residual_crack=res_metrics,
    )


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--flank", type=int, default=2000,
                        help="bp of upstream/downstream sequence to extract")
    parser.add_argument("--mask-repeats", action="store_true",
                        help="mask soft-masked (lowercase) bases to N before alignment")
    parser.add_argument("--out-suffix", default="",
                        help="suffix appended to output filenames before extension")
    args = parser.parse_args()

    global FLANK, MASK_REPEATS, OUT_TSV, OUT_JSON
    FLANK = args.flank
    MASK_REPEATS = args.mask_repeats
    if args.out_suffix:
        base, ext = os.path.splitext(OUT_TSV)
        OUT_TSV = f"{base}_{args.out_suffix}{ext}"
        base, ext = os.path.splitext(OUT_JSON)
        OUT_JSON = f"{base}_{args.out_suffix}{ext}"

    print(f"[*] flank size = {FLANK} bp; mask_repeats={MASK_REPEATS}; "
          f"outputs: {OUT_TSV}, {OUT_JSON}", flush=True)
    print("[*] loading meta / edge tables ...", flush=True)
    meta = load_meta()
    dn_map = load_dn_map()
    lever_rows = load_lever()

    # restrict to edges present in both tables
    edges = []
    for row in lever_rows:
        k = frozenset((row["gene_a"], row["gene_b"]))
        if k not in dn_map:
            continue
        dn_a, dn_b = dn_map[k]
        row = dict(row)
        row["dn_a"] = dn_a
        row["dn_b"] = dn_b
        edges.append(row)
    print(f"    edges with DN mapping: {len(edges)}", flush=True)

    # collect unique DN loci
    needed_dn = set()
    for r in edges:
        needed_dn.add(r["dn_a"])
        needed_dn.add(r["dn_b"])
    print(f"    unique DN loci: {len(needed_dn)}", flush=True)

    print("[*] extracting flanking sequences ...", flush=True)
    fa = pysam.FastaFile(FASTA)
    up_seqs, down_seqs = {}, {}
    short_count = 0
    for dn in sorted(needed_dn):
        up = flank_seq(fa, meta, dn, "up")
        down = flank_seq(fa, meta, dn, "down")
        if len(up) < FLANK // 2 or len(down) < FLANK // 2:
            short_count += 1
        up_seqs[f"up:{dn}"] = up if up else "N"
        down_seqs[f"down:{dn}"] = down if down else "N"
    print(f"    loci with short (<{FLANK//2} bp) flank: {short_count}", flush=True)

    with tempfile.TemporaryDirectory(prefix="flank_hom_") as td:
        up_fa = os.path.join(td, "up.fa")
        down_fa = os.path.join(td, "down.fa")
        up_paf = os.path.join(td, "up.paf")
        down_paf = os.path.join(td, "down.paf")
        write_fasta(up_fa, up_seqs)
        write_fasta(down_fa, down_seqs)

        print("[*] running minimap2 upstream ...", flush=True)
        run_minimap2(up_fa, up_fa, up_paf)
        print("[*] running minimap2 downstream ...", flush=True)
        run_minimap2(down_fa, down_fa, down_paf)

        print("[*] parsing PAF ...", flush=True)
        up_best = parse_paf(up_paf)
        down_best = parse_paf(down_paf)

    print("[*] computing edge features ...", flush=True)
    out_rows = []
    for r in edges:
        feats = flank_features(up_best, down_best, r["dn_a"], r["dn_b"])
        out_row = dict(r)
        out_row.update(feats)
        # drop heavy intermediate columns
        out_row.pop("dn_a", None)
        out_row.pop("dn_b", None)
        out_rows.append(out_row)

    # write TSV
    cols = ["gene_a", "gene_b", "in_dna_loose", "cls", "core_recip", "aln_frac",
            "fold", "residual",
            "up_aln_frac", "up_id", "down_aln_frac", "down_id",
            "max_flank_aln_frac", "min_flank_aln_frac", "mean_flank_aln_frac",
            "best_flank_id"]
    with open(OUT_TSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        for row in out_rows:
            w.writerow({c: row.get(c, NA) for c in cols})
    print(f"[+] wrote {OUT_TSV}", flush=True)

    print("[*] evaluating ...", flush=True)
    metrics = evaluate(out_rows)
    metrics["params"] = dict(flank_bp=FLANK, minimap_preset=MINIMAP_PRESET,
                             mask_repeats=MASK_REPEATS)
    with open(OUT_JSON, "w") as fh:
        json.dump(metrics, fh, indent=2, sort_keys=True)
    print(f"[+] wrote {OUT_JSON}", flush=True)

    print("\n=== SUMMARY ===")
    print(f"n_edges = {metrics['n_edges']}  pos={metrics['n_pos']} neg={metrics['n_neg']}")
    print(f"baseline core+aln OOF AUC = {metrics['baseline_core_aln_oof_auc']}")
    print("univariate pooled AUC:")
    for name, v in sorted(metrics["univariate"].items()):
        print(f"    {name:24s} {v['auc']:.4f}  (n_pos={v['n_pos']} n_neg={v['n_neg']})")
    print("OOF logistic lift over core+aln:")
    for name, v in sorted(metrics["oof_combinations"].items()):
        print(f"    {name:24s} {v['auc']:.4f}  lift={v['lift']:+.4f}")
    if metrics["residual_crack"]:
        print("residual-crack AUC:")
        for name, v in sorted(metrics["residual_crack"].items()):
            print(f"    {name:24s} {v['auc']:.4f}  (n_pos={v['n_pos']} n_neg={v['n_neg']})")


if __name__ == "__main__":
    main()
