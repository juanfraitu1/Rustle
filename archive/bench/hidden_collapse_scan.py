#!/usr/bin/env python3
"""Task (c): scan single-copy ANNOTATED recall loci for HIDDEN collapsed-copy PSV signal.

The headroom probe's main caveat: it detects COLLAPSED paralog loci only via ANNOTATED overlap.
Cross-mapping can pile a divergent / unannotated paralog copy's reads onto a locus the annotation
calls single-copy, creating copy-discriminating PSV signal the annotation hides. If that is common,
the read-coherence recall isoforms at "single-copy" loci might actually be PSV-resolvable -> the
headroom would be > 0.

This detector calls PSVs DE NOVO from the BAM pileup (annotation-free) and tests CO-SEGREGATION
(the copy_split identifiability signal). The hard confound is a HETEROZYGOUS gene: linked het SNPs
also co-segregate into 2 read groups. Discriminators (baked in, pressure-tested by the workflow):
  - LINKED COLUMN COUNT  : diploid het of a transcript carries ~few linked SNPs (genome het ~1/1kb);
                           a collapsed paralog at >1% divergence carries many fixed PSVs. -> n_coseg.
  - GROUP COUNT          : >=3 co-segregating allele groups cannot be one diploid locus -> >=3 copies.
  - MULTIMAPPING         : collapsed paralog reads multimap (secondary alns / MAPQ 0); het maps unique.
We emit per-locus features + a tiered verdict; the honest headroom = COLLAPSED_LIKE count.

Reuses the VALIDATED co-segregation logic from copy_split_sim (faithful port of copy_split.rs).
Deterministic (no RNG). CLI: --chrom to parallelize one contig.
"""
import argparse
import os
import sys
from collections import defaultdict

import numpy as np
import pysam

sys.path.insert(0, os.path.dirname(__file__))
from copy_split_sim import (ReadObs, split_readchain_by_psv, distinct_columns,  # noqa: E402
                            read_consistent_with)

BASES = set("ACGT")
DUMMY_CHAIN = ((0, 0),)  # all reads at one locus share a placeholder chain (we test allele linkage)


_NOT_UNMAPPED = lambda r: not r.is_unmapped  # include secondary/supplementary (multimap signal)
MAX_CAND = 400        # cap candidate columns processed per locus (a real collapsed block is < this)


def call_psv_columns(bam, fasta, chrom, start, end, d_min, c_min, f_min, max_depth):
    """Two-stage. Stage 1 (FAST, count_coverage = C-level): per-base A/C/G/T counts over the
    whole locus -> candidate PSV columns (minor allele >= c_min, freq >= f_min, depth in range)
    with NO per-read Python loop. The 99% of loci with < min_k candidate columns are rejected
    here in milliseconds. Stage 2 (only for survivors): a single-column pileup at each candidate
    column to recover per-read alleles for the co-segregation test, plus the multimapping
    signature (mapq / secondary) from those reads.
    Returns (cols, read_cols, (frac_sec, frac_mq0, n_reads_at_cols))."""
    L = end - start
    if L <= 0:
        return [], {}, (0.0, 0.0, 0)
    cov = bam.count_coverage(chrom, start, end, quality_threshold=0, read_callback=_NOT_UNMAPPED)
    stack = np.vstack([np.asarray(cov[i], dtype=np.int64) for i in range(4)])  # 4 x L (A,C,G,T)
    depth = stack.sum(axis=0)
    srt = np.sort(stack, axis=0)
    minor = srt[2]   # second-largest base count per column
    ok = (depth >= d_min) & (depth <= max_depth) & (minor >= c_min) & (minor >= f_min * depth)
    cand_local = np.nonzero(ok)[0]
    if len(cand_local) == 0:
        return [], {}, (0.0, 0.0, 0)
    if len(cand_local) > MAX_CAND:
        cand_local = cand_local[:MAX_CAND]
    cols = [start + int(i) for i in cand_local]

    # stage 2: ONE pileup over the candidate span. Recovers per-read alleles from PRIMARY
    # (SEQ-bearing) reads, the MAPQ-0 fraction among them (minimap2 MAPQ 0 = read maps equally
    # well elsewhere => a paralog copy exists), and the secondary-record density (cross-mapping)
    # — all in a single pass. Secondary alignments here carry no SEQ, so alleles come only from
    # primaries. The span is capped so scattered (non-block) candidate sets stay cheap.
    # Phase the DENSEST WINDOW_CAP-wide cluster of candidates (a collapsed copy's divergent block
    # is contiguous and transcript-sized; this finds it wherever it sits, not just at cols[0]).
    WINDOW_CAP = 30000
    if cols[-1] - cols[0] < WINDOW_CAP:
        span_lo, span_hi = cols[0], cols[-1] + 1
    else:
        best_i, best_cnt, j = 0, 0, 0
        for i in range(len(cols)):
            if j < i:
                j = i
            while j + 1 < len(cols) and cols[j + 1] - cols[i] < WINDOW_CAP:
                j += 1
            if j - i + 1 > best_cnt:
                best_cnt, best_i = j - i + 1, i
        span_lo = cols[best_i]
        span_hi = span_lo + WINDOW_CAP
    candset = {pos: ci for ci, pos in enumerate(cols) if span_lo <= pos < span_hi}
    read_cols = defaultdict(dict)
    seen = {}            # primary read_name -> mapq
    sec_names = set()    # secondary/supp read names touching candidate columns
    for pc in bam.pileup(chrom, span_lo, span_hi, truncate=True, max_depth=max_depth,
                         stepper="nofilter", ignore_overlaps=False, min_base_quality=0):
        ci = candset.get(pc.reference_pos)
        if ci is None:
            continue
        for pr in pc.pileups:
            aln = pr.alignment
            if aln.is_secondary or aln.is_supplementary:
                sec_names.add(aln.query_name)
                continue
            if pr.is_del or pr.is_refskip or pr.query_position is None:
                continue
            seq = aln.query_sequence
            if seq is None:
                continue
            base = seq[pr.query_position]
            if base not in BASES:
                continue
            read_cols[aln.query_name][ci] = base
            if aln.query_name not in seen:
                seen[aln.query_name] = aln.mapping_quality
    n = len(seen)
    frac_mq0 = (sum(1 for q in seen.values() if q == 0) / n) if n else 0.0
    denom = n + len(sec_names - set(seen))
    frac_sec = (len(sec_names - set(seen)) / denom) if denom else 0.0
    return cols, read_cols, (frac_sec, frac_mq0, n)


def scan_locus(bam, fasta, chrom, start, end, args):
    cols, read_cols, sig = call_psv_columns(
        bam, fasta, chrom, start, end, args.d_min, args.c_min, args.f_min, args.max_depth)
    frac_sec0, frac_mq0_0, n_aln0 = sig
    m = len(cols)
    if m < args.min_k:
        return dict(verdict="NO_SIGNAL", n_psv=m, n_coseg=0, n_copies=0,
                    frac_sec=0.0, frac_mq0=0.0, n_aln=0, minor_frac=0.0)

    # build ReadObs over the m candidate columns; only reads covering >=2 columns inform linkage.
    base2int = {"A": 65, "C": 67, "G": 71, "T": 84}
    reads = []
    for rname, colmap in read_cols.items():
        if len(colmap) < 2:
            continue
        vec = tuple(base2int.get(colmap.get(ci)) if ci in colmap else None for ci in range(m))
        reads.append(ReadObs(DUMMY_CHAIN, vec))
    if len(reads) < args.min_reads * 2:
        return dict(verdict="NO_SIGNAL", n_psv=m, n_coseg=0, n_copies=0,
                    frac_sec=frac_sec0, frac_mq0=frac_mq0_0, n_aln=n_aln0, minor_frac=0.0)

    out = split_readchain_by_psv(reads, m, args.min_k, args.min_reads)
    copies = [c for c in out if c.identifiable]
    n_copies = len(copies)

    frac_sec, frac_mq0, n_aln = frac_sec0, frac_mq0_0, n_aln0

    if n_copies < 2:
        return dict(verdict="NO_SIGNAL", n_psv=m, n_coseg=0, n_copies=n_copies,
                    frac_sec=frac_sec, frac_mq0=frac_mq0, n_aln=n_aln, minor_frac=0.0)

    # n_coseg = columns that actually discriminate the two largest copies (the linked block size).
    copies_sorted = sorted(copies, key=lambda c: -c.read_count)
    n_coseg = distinct_columns(copies_sorted[0].allele_vector, copies_sorted[1].allele_vector)
    tot = sum(c.read_count for c in copies)
    minor_frac = copies_sorted[1].read_count / tot if tot else 0.0

    # ---- tiered verdict (het vs collapsed) ----
    multimap = frac_sec >= args.sec_thr or frac_mq0 >= args.mq0_thr
    if n_copies >= 3 or n_coseg >= args.high_col or multimap:
        verdict = "COLLAPSED_LIKE"
    elif n_copies == 2 and n_coseg < args.high_col and not multimap:
        verdict = "HET_LIKE"
    else:
        verdict = "AMBIGUOUS"
    return dict(verdict=verdict, n_psv=m, n_coseg=n_coseg, n_copies=n_copies,
                frac_sec=frac_sec, frac_mq0=frac_mq0, n_aln=n_aln, minor_frac=minor_frac)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", default="/home/juanfra/winloci_scratch/GGO.bam")
    ap.add_argument("--fasta", default="/home/juanfra/winloci_scratch/GGO.fasta")
    ap.add_argument("--loci", default=os.path.join(os.path.dirname(__file__),
                                                   "single_copy_recall_loci.tsv"))
    ap.add_argument("--chrom", default=None, help="restrict to one contig (parallelism)")
    ap.add_argument("--region", default=None, help="chrom:start-end single-locus control mode")
    ap.add_argument("--out", default=None)
    # PSV calling thresholds
    ap.add_argument("--d-min", dest="d_min", type=int, default=6)     # min depth to call a column
    ap.add_argument("--c-min", dest="c_min", type=int, default=3)     # min minor-allele read count
    ap.add_argument("--f-min", dest="f_min", type=float, default=0.20)  # min minor-allele freq
    ap.add_argument("--max-depth", dest="max_depth", type=int, default=1500)
    # co-segregation thresholds
    ap.add_argument("--min-k", dest="min_k", type=int, default=3)     # min co-seg columns to split
    ap.add_argument("--min-reads", dest="min_reads", type=int, default=3)
    # het-vs-collapsed discriminators
    ap.add_argument("--high-col", dest="high_col", type=int, default=8)  # >=this linked cols -> collapsed
    ap.add_argument("--sec-thr", dest="sec_thr", type=float, default=0.10)
    ap.add_argument("--mq0-thr", dest="mq0_thr", type=float, default=0.30)
    args = ap.parse_args()

    bam = pysam.AlignmentFile(args.bam, "rb")
    fasta = pysam.FastaFile(args.fasta)

    if args.region:
        c, rng = args.region.split(":")
        s, e = (int(x.replace(",", "")) for x in rng.split("-"))
        r = scan_locus(bam, fasta, c, s, e, args)
        print(f"{args.region}\t" + "\t".join(f"{k}={v}" for k, v in r.items()))
        return

    loci = []
    with open(args.loci) as fh:
        fh.readline()
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 5:
                continue
            if args.chrom and p[0] != args.chrom:
                continue
            loci.append((p[0], int(p[1]), int(p[2]), int(p[3]), int(p[4])))

    counts = defaultdict(int)
    rows = []
    for chrom, s, e, niso, nfsm in loci:
        r = scan_locus(bam, fasta, chrom, s, e, args)
        counts[r["verdict"]] += 1
        # dump every locus with a linked PSV block (n_copies>=2) so the het/collapse boundary
        # can be re-thresholded without re-scanning; NO_SIGNAL loci kept as counts only.
        if r["verdict"] in ("COLLAPSED_LIKE", "AMBIGUOUS", "HET_LIKE"):
            rows.append((chrom, s, e, niso, nfsm, r))

    out = args.out or (f"/tmp/gw/hidden_collapse_{args.chrom}.tsv" if args.chrom
                       else "/tmp/gw/hidden_collapse_all.tsv")
    with open(out, "w") as fh:
        fh.write("chrom\tstart\tend\tn_isoforms\tn_fsm\tverdict\tn_psv\tn_coseg\t"
                 "n_copies\tfrac_sec\tfrac_mq0\tn_aln\tminor_frac\n")
        for chrom, s, e, niso, nfsm, r in rows:
            fh.write(f"{chrom}\t{s}\t{e}\t{niso}\t{nfsm}\t{r['verdict']}\t{r['n_psv']}\t"
                     f"{r['n_coseg']}\t{r['n_copies']}\t{r['frac_sec']:.3f}\t{r['frac_mq0']:.3f}\t"
                     f"{r['n_aln']}\t{r['minor_frac']:.3f}\n")
    total = sum(counts.values())
    summary = " ".join(f"{k}={v}" for k, v in sorted(counts.items()))
    print(f"[{args.chrom or 'ALL'}] loci={total} {summary} -> {out}")


if __name__ == "__main__":
    main()
