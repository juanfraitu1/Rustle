#!/usr/bin/env python3
"""Per-locus EVIDENCE for the hidden-collapse verification (task c).

Given a region (a COLLAPSED_LIKE hit), dump the hard evidence needed to judge whether it is a
REAL hidden collapsed/cross-mapped paralog copy or a false-positive mode:

  base-change spectrum  -> RNA editing (A>G on +, T>C on -; clustered, transition-only) vs paralog
                           divergence (all 12 substitution types, ~2:1 transition bias but mixed).
  copy haplotypes       -> the >=2 co-segregating allele groups, read counts, minor-copy fraction
                           (clean balanced copies vs noisy single-column artifact).
  mapping structure     -> MAPQ histogram, MAPQ-0 fraction (this locus's reads map elsewhere too),
                           secondary-record fraction (reads from elsewhere spill here) = cross-mapping.
  reference context     -> low-complexity / homopolymer / STR fraction around the block (alignment
                           artifacts pile mismatches at repeats).

Reuses the detector's column caller + the validated co-segregation port. Deterministic.
"""
import argparse
import os
import sys
from collections import defaultdict, Counter

import pysam

sys.path.insert(0, os.path.dirname(__file__))
from hidden_collapse_scan import call_psv_columns, BASES  # noqa: E402
from copy_split_sim import (ReadObs, split_readchain_by_psv, distinct_columns)  # noqa: E402

DUMMY_CHAIN = ((0, 0),)


def lowcomplexity_frac(seq):
    """Fraction of positions in homopolymer runs (>=4) or simple di/tri repeats — a proxy for the
    low-complexity context where aligners pile spurious mismatches."""
    seq = seq.upper()
    if not seq:
        return 0.0
    flagged = [False] * len(seq)
    # homopolymer runs >= 4
    i = 0
    while i < len(seq):
        j = i
        while j + 1 < len(seq) and seq[j + 1] == seq[i]:
            j += 1
        if j - i + 1 >= 4:
            for k in range(i, j + 1):
                flagged[k] = True
        i = j + 1
    return sum(flagged) / len(seq)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", default="/home/juanfra/winloci_scratch/GGO.bam")
    ap.add_argument("--fasta", default="/home/juanfra/winloci_scratch/GGO.fasta")
    ap.add_argument("--region", required=True)
    ap.add_argument("--d-min", dest="d_min", type=int, default=6)
    ap.add_argument("--c-min", dest="c_min", type=int, default=3)
    ap.add_argument("--f-min", dest="f_min", type=float, default=0.20)
    ap.add_argument("--max-depth", dest="max_depth", type=int, default=1500)
    ap.add_argument("--min-k", dest="min_k", type=int, default=3)
    ap.add_argument("--min-reads", dest="min_reads", type=int, default=3)
    args = ap.parse_args()

    bam = pysam.AlignmentFile(args.bam, "rb")
    fasta = pysam.FastaFile(args.fasta)
    c, rng = args.region.split(":")
    s, e = (int(x.replace(",", "")) for x in rng.split("-"))

    cols, read_cols, sig = call_psv_columns(bam, fasta, c, s, e,
                                            args.d_min, args.c_min, args.f_min, args.max_depth)
    frac_sec, frac_mq0, n_reads = sig
    m = len(cols)
    print(f"# region {args.region}  len={e-s}  candidate_PSV_columns={m}  "
          f"reads_at_cols={n_reads}  frac_mq0={frac_mq0:.3f}  frac_sec={frac_sec:.3f}")
    if m == 0:
        return

    # base-change spectrum (ref -> minor) over candidate columns
    spectrum = Counter()
    transitions = {"AG", "GA", "CT", "TC"}
    n_ti = n_tv = 0
    col_index = {pos: ci for ci, pos in enumerate(cols)}
    for pos in cols:
        ref = fasta.fetch(c, pos, pos + 1).upper()
        ci = col_index[pos]
        tally = Counter()
        for cm in read_cols.values():
            if ci in cm:
                tally[cm[ci]] += 1
        if not tally:
            continue
        order = tally.most_common()
        major = order[0][0]
        minor = order[1][0] if len(order) > 1 else None
        if minor and ref in BASES:
            change = f"{ref}>{minor}" if minor != ref else f"{ref}>{major}"
            sub = (ref + minor) if minor != ref else (ref + major)
            spectrum[change] += 1
            if sub in transitions:
                n_ti += 1
            else:
                n_tv += 1
    print("# base-change spectrum (ref>minor):",
          ", ".join(f"{k}:{v}" for k, v in spectrum.most_common()))
    edit_like = spectrum.get("A>G", 0) + spectrum.get("T>C", 0)
    tot_spec = sum(spectrum.values())
    print(f"#   transitions={n_ti} transversions={n_tv}  "
          f"A>G+T>C fraction={edit_like/max(tot_spec,1):.2f}  "
          f"(>0.8 => RNA-editing-like; mixed 12-type + Ti:Tv~2 => paralog divergence)")

    # co-segregation haplotypes
    reads = []
    for cm in read_cols.values():
        if len(cm) < 2:
            continue
        vec = tuple(ord(cm[ci]) if ci in cm else None for ci in range(m))
        reads.append(ReadObs(DUMMY_CHAIN, vec))
    out = split_readchain_by_psv(reads, m, args.min_k, args.min_reads)
    copies = sorted([x for x in out if x.identifiable], key=lambda z: -z.read_count)
    print(f"# co-segregating copies={len(copies)} (informative reads={len(reads)})")
    for i, cp in enumerate(copies[:4]):
        hap = "".join(chr(b) if b else "." for b in cp.allele_vector)
        print(f"#   copy{i} reads={cp.read_count} hap={hap}")
    if len(copies) >= 2:
        print(f"#   linked cols (copy0 vs copy1) = "
              f"{distinct_columns(copies[0].allele_vector, copies[1].allele_vector)}")
        tot = sum(cp.read_count for cp in copies)
        print(f"#   minor-copy fraction = {copies[1].read_count/tot:.2f}")

    # MAPQ histogram (primary reads at the locus)
    mapq = Counter()
    nprim = 0
    for aln in bam.fetch(c, s, e):
        if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
            continue
        nprim += 1
        mapq[min(aln.mapping_quality, 60)] += 1
        if nprim >= 5000:
            break
    print("# primary MAPQ hist (val:count):",
          ", ".join(f"{k}:{v}" for k, v in sorted(mapq.items())))

    # reference low-complexity around the candidate block
    blk_lo, blk_hi = cols[0], cols[-1] + 1
    refseq = fasta.fetch(c, max(0, blk_lo - 20), blk_hi + 20)
    print(f"# block span={blk_hi-blk_lo}bp  ref low-complexity(homopolymer>=4) frac="
          f"{lowcomplexity_frac(refseq):.2f}")


if __name__ == "__main__":
    main()
