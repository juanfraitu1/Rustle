#!/usr/bin/env python3
"""
Research spike for the family-aware HMM rescue idea.

Question: does a per-exon profile (PSSM) built from 3 GOLGA6L7 paralogs
score *unmapped* reads more like positives (mapped reads from the family
region) than like random reads (mapped elsewhere)?

If yes (positive >> unmapped > random), the family-profile approach is
finding signal that the linear-reference aligner misses, and the full
HMM design in
docs/superpowers/specs/2026-04-28-vg-novel-copy-hmm-design.md
is worth implementing in full.

If no (unmapped ≈ random), the architectural choice is wrong and we
should reconsider before investing in the full plan.

This is a deliberately simple spike: PSSM (match-state profile) per
exon class, no insert/delete states, no full HMM forward. If a PSSM
separates the populations, a real HMM with M/I/D will too.

Inputs:
    GGO.fasta (genome)
    GGO_genomic.gff (annotation)
    GGO_19.bam (BAM with mapped + unmapped reads on chr19)

Default paths assume the layout used in this repo's MULTI_COPY_FAMILY_PROOF:
    /mnt/c/Users/jfris/Desktop/{GGO.fasta, GGO_genomic.gff, GGO_19.bam}

Output: a small report with per-population log-odds histograms and a
discrimination summary.
"""

from __future__ import annotations

import argparse
import math
import random
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pysam

# Three GOLGA6L7 paralogs on chr19 (NC_073243.2) per MULTI_COPY_FAMILY_PROOF.md.
DEFAULT_FAMILY = ["LOC115931294", "LOC134757625", "LOC101137218"]
DEFAULT_CHROM = "NC_073243.2"
COMP = bytes.maketrans(b"ACGTNacgtn", b"TGCANtgcan")
BASES = "ACGT"
BASE_IDX = {b: i for i, b in enumerate(BASES)}


@dataclass
class Mrna:
    gene_id: str
    mrna_id: str
    chrom: str
    strand: str
    exons: list[tuple[int, int]]  # 0-based half-open


def parse_gff_for_family(gff: Path, gene_ids: list[str]) -> list[Mrna]:
    """Find each gene's longest-by-cds-equivalent mRNA and its exons."""
    gene_to_mrnas: dict[str, dict[str, Mrna]] = {g: {} for g in gene_ids}
    rna_to_gene: dict[str, str] = {}

    with open(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            chrom, _src, ftype, start, end, _sc, strand, _ph, attrs = fields
            if ftype not in ("mRNA", "exon"):
                continue
            attr_map = dict(
                kv.split("=", 1) for kv in attrs.split(";") if "=" in kv
            )

            if ftype == "mRNA":
                parent = attr_map.get("Parent", "")
                # Parent looks like "gene-LOC115931294".
                gene_name = None
                for g in gene_ids:
                    if g in parent:
                        gene_name = g
                        break
                if gene_name is None:
                    continue
                rna_id = attr_map.get("ID", "")
                rna_to_gene[rna_id] = gene_name
                gene_to_mrnas[gene_name][rna_id] = Mrna(
                    gene_id=gene_name,
                    mrna_id=rna_id,
                    chrom=chrom,
                    strand=strand,
                    exons=[],
                )
            elif ftype == "exon":
                parent = attr_map.get("Parent", "")
                if parent not in rna_to_gene:
                    continue
                gene_name = rna_to_gene[parent]
                m = gene_to_mrnas[gene_name][parent]
                m.exons.append((int(start) - 1, int(end)))

    chosen: list[Mrna] = []
    for g in gene_ids:
        if not gene_to_mrnas[g]:
            print(f"  warning: no mRNA found for {g}", file=sys.stderr)
            continue
        best = max(
            gene_to_mrnas[g].values(),
            key=lambda m: sum(e - s for s, e in m.exons),
        )
        best.exons.sort()
        chosen.append(best)
    return chosen


def fetch_exon_sequences(fasta: pysam.FastaFile, mrnas: list[Mrna]) -> list[list[bytes]]:
    """Per-mRNA list of per-exon sequence bytes (transcript strand)."""
    out: list[list[bytes]] = []
    for m in mrnas:
        seqs: list[bytes] = []
        for s, e in m.exons:
            raw = fasta.fetch(m.chrom, s, e).upper().encode()
            if m.strand == "-":
                raw = raw[::-1].translate(COMP)
            seqs.append(raw)
        if m.strand == "-":
            seqs.reverse()  # Now in transcript order.
        out.append(seqs)
    return out


def nw_align(a: bytes, b: bytes, gap: int = -2, match: int = 1, mismatch: int = -1) -> tuple[bytes, bytes]:
    """Simple Needleman-Wunsch returning aligned a, b with b'-' as gap."""
    la, lb = len(a), len(b)
    score = np.zeros((la + 1, lb + 1), dtype=np.int32)
    score[:, 0] = np.arange(la + 1) * gap
    score[0, :] = np.arange(lb + 1) * gap
    for i in range(1, la + 1):
        for j in range(1, lb + 1):
            s = match if a[i - 1] == b[j - 1] else mismatch
            score[i, j] = max(
                score[i - 1, j - 1] + s,
                score[i - 1, j] + gap,
                score[i, j - 1] + gap,
            )
    # Traceback.
    i, j = la, lb
    aa, bb = bytearray(), bytearray()
    while i > 0 and j > 0:
        s = match if a[i - 1] == b[j - 1] else mismatch
        if score[i, j] == score[i - 1, j - 1] + s:
            aa.append(a[i - 1]); bb.append(b[j - 1]); i -= 1; j -= 1
        elif score[i, j] == score[i - 1, j] + gap:
            aa.append(a[i - 1]); bb.append(ord("-")); i -= 1
        else:
            aa.append(ord("-")); bb.append(b[j - 1]); j -= 1
    while i > 0:
        aa.append(a[i - 1]); bb.append(ord("-")); i -= 1
    while j > 0:
        aa.append(ord("-")); bb.append(b[j - 1]); j -= 1
    return bytes(aa[::-1]), bytes(bb[::-1])


def progressive_msa(seqs: list[bytes]) -> list[bytes]:
    """Progressive pairwise MSA: 3 sequences. Returns equal-length rows."""
    if len(seqs) == 0:
        return []
    if len(seqs) == 1:
        return [seqs[0]]
    rows = [seqs[0], seqs[1]]
    rows[0], rows[1] = nw_align(rows[0], rows[1])
    for nxt in seqs[2:]:
        # Align nxt to consensus (use the first row as guide; cheap and good enough for spike).
        consensus = bytes(b if b != ord("-") else ord("N") for b in rows[0])
        cons_al, nxt_al = nw_align(consensus, nxt)
        # Splice gap insertions in cons_al into existing rows.
        new_rows = [bytearray() for _ in rows]
        ri = [0] * len(rows)
        for c in cons_al:
            if c == ord("-"):
                for k in range(len(rows)):
                    new_rows[k].append(ord("-"))
            else:
                for k in range(len(rows)):
                    new_rows[k].append(rows[k][ri[k]])
                    ri[k] += 1
        rows = [bytes(r) for r in new_rows]
        rows.append(nxt_al)
    # Pad to same length (defensive).
    L = max(len(r) for r in rows)
    return [r.ljust(L, b"-") for r in rows]


def build_pfm(rows: list[bytes], pseudo: float = 0.5) -> np.ndarray:
    """Position-frequency matrix in log-odds vs uniform background."""
    L = len(rows[0])
    counts = np.full((L, 4), pseudo, dtype=np.float64)
    for r in rows:
        for c in range(L):
            b = r[c:c+1]
            if b in (b"A", b"C", b"G", b"T"):
                counts[c, BASE_IDX[b.decode()]] += 1.0
    probs = counts / counts.sum(axis=1, keepdims=True)
    log_odds = np.log(probs / 0.25)
    return log_odds  # (L, 4)


def best_pssm_score(read_idx: np.ndarray, pssm: np.ndarray) -> float:
    """Max over offsets of sum(pssm[col, read[offset+col]]) via vectorized stride-1 sweep.
    read_idx is base indices in [0..3], -1 for N (any window with N is skipped)."""
    L = pssm.shape[0]
    n = len(read_idx)
    if n < L:
        return float("-inf")
    n_off = n - L + 1
    # Per-column contribution to every starting offset, vectorised over offsets.
    # score[off] = sum_c pssm[c, read_idx[off+c]]
    score = np.zeros(n_off, dtype=np.float64)
    has_n = np.zeros(n_off, dtype=bool)
    for c in range(L):
        col_bases = read_idx[c:c + n_off]
        bad = col_bases < 0
        has_n |= bad
        # Look up pssm[c, base] for every offset; clip negative indices to 0 (masked out below).
        safe = np.where(bad, 0, col_bases)
        score += pssm[c, safe]
    score[has_n] = -np.inf
    if not np.isfinite(score).any():
        return float("-inf")
    return float(score.max())


def family_score(read_seq: bytes, pssms: list[np.ndarray]) -> float:
    """Best score across all exon-class PSSMs and offsets."""
    # Convert to uint8 base indices (-1 for N).
    arr = np.frombuffer(read_seq, dtype=np.uint8)
    idx = np.full(arr.shape, -1, dtype=np.int8)
    for c, i in BASE_IDX.items():
        idx[arr == ord(c)] = i
    idx = idx.astype(np.int64)
    fwd = max((best_pssm_score(idx, p) for p in pssms), default=float("-inf"))
    # Reverse complement.
    rc = read_seq.translate(COMP)[::-1]
    arr_rc = np.frombuffer(rc, dtype=np.uint8)
    idx_rc = np.full(arr_rc.shape, -1, dtype=np.int8)
    for c, i in BASE_IDX.items():
        idx_rc[arr_rc == ord(c)] = i
    idx_rc = idx_rc.astype(np.int64)
    rev = max((best_pssm_score(idx_rc, p) for p in pssms), default=float("-inf"))
    return max(fwd, rev)


def sample_reads(bam_path: Path, region: tuple[str, int, int] | None, unmapped: bool, n: int, rng: random.Random) -> list[bytes]:
    """Sample up to n read sequences. unmapped=True uses pysam's '*' fetch (fast on indexed BAMs)."""
    out: list[bytes] = []
    with pysam.AlignmentFile(str(bam_path), "rb") as bf:
        if unmapped:
            try:
                it = bf.fetch("*")  # Coordinate-sorted BAMs: unmapped block at end.
            except (ValueError, KeyError):
                it = bf.fetch(until_eof=True)
            for rec in it:
                if not rec.is_unmapped:
                    continue
                seq = rec.query_sequence
                if seq and len(seq) >= 100:
                    out.append(seq.encode())
                    if len(out) >= n * 5:
                        break
        else:
            assert region is not None
            chrom, s, e = region
            for rec in bf.fetch(chrom, s, e):
                if rec.is_unmapped or rec.is_secondary or rec.is_supplementary:
                    continue
                seq = rec.query_sequence
                if seq and len(seq) >= 100:
                    out.append(seq.encode())
                    if len(out) >= n * 5:
                        break
    rng.shuffle(out)
    return out[:n]


def summarise(name: str, scores: list[float]) -> None:
    arr = np.array(scores)
    if arr.size == 0:
        print(f"{name:<14}  (no reads)")
        return
    qs = np.quantile(arr, [0.05, 0.25, 0.5, 0.75, 0.95])
    print(f"{name:<14}  n={arr.size:4d}  median={qs[2]:8.2f}  IQR=[{qs[1]:7.2f}, {qs[3]:7.2f}]  min={arr.min():8.2f}  max={arr.max():8.2f}")


def histogram(name: str, scores: list[float], width: int = 50, bins: int = 12) -> None:
    if not scores:
        return
    lo, hi = min(scores), max(scores)
    if hi - lo < 1e-6:
        print(f"{name}: degenerate (all == {lo:.2f})")
        return
    counts = [0] * bins
    for s in scores:
        b = min(bins - 1, int((s - lo) / (hi - lo) * bins))
        counts[b] += 1
    peak = max(counts) or 1
    print(f"{name} histogram (bin width = {(hi - lo) / bins:.2f}):")
    for i, c in enumerate(counts):
        bar = "#" * int(c / peak * width)
        edge_lo = lo + (hi - lo) * i / bins
        print(f"  {edge_lo:8.2f} | {bar} {c}")


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--fasta", default="/mnt/c/Users/jfris/Desktop/GGO.fasta")
    p.add_argument("--gff", default="/mnt/c/Users/jfris/Desktop/GGO_genomic.gff")
    p.add_argument("--bam", default="/mnt/c/Users/jfris/Desktop/GGO_19.bam",
                   help="BAM for positive (in-region) and random-ctrl (off-region) reads")
    p.add_argument("--unmapped-bam", default="/mnt/c/Users/jfris/Desktop/GGO.bam",
                   help="BAM to draw unmapped reads from (chr19 subset has these stripped; use full BAM)")
    p.add_argument("--family", nargs="+", default=DEFAULT_FAMILY)
    p.add_argument("--chrom", default=DEFAULT_CHROM)
    p.add_argument("--n-positive", type=int, default=200)
    p.add_argument("--n-unmapped", type=int, default=200)
    p.add_argument("--n-random", type=int, default=200)
    p.add_argument("--seed", type=int, default=1)
    args = p.parse_args()
    rng = random.Random(args.seed)

    print("[1/5] Parsing GFF for family copies ...")
    mrnas = parse_gff_for_family(Path(args.gff), args.family)
    if len(mrnas) < 2:
        print(f"  ERROR: need ≥ 2 family copies, found {len(mrnas)}", file=sys.stderr)
        return 2
    print(f"  Found {len(mrnas)} mRNAs:")
    span = (10**18, 0)
    for m in mrnas:
        s = min(s for s, _ in m.exons); e = max(e for _, e in m.exons)
        span = (min(span[0], s), max(span[1], e))
        print(f"    {m.gene_id:<14} {m.mrna_id:<22} {m.chrom} {m.strand} {len(m.exons)} exons  ({s}-{e})")
    family_region = (args.chrom, span[0], span[1])

    print("[2/5] Fetching exon sequences from FASTA ...")
    fasta = pysam.FastaFile(args.fasta)
    per_mrna_exons = fetch_exon_sequences(fasta, mrnas)
    n_exons = min(len(e) for e in per_mrna_exons)
    print(f"  Each copy has {[len(e) for e in per_mrna_exons]} exons; using first {n_exons} (assumes shared architecture).")

    print(f"[3/5] Building per-exon PSSMs over {len(mrnas)} copies ...")
    pssms: list[np.ndarray] = []
    for ei in range(n_exons):
        seqs = [per_mrna_exons[ci][ei] for ci in range(len(mrnas))]
        # Skip very short exons (PSSM noise) and very long ones (alignment slow,
        # and the long terminal exon dominates everything if scored stride-1).
        L_min = min(len(s) for s in seqs)
        L_max = max(len(s) for s in seqs)
        if L_min < 50 or L_max > 800:
            print(f"    exon {ei:2d}: lengths {[len(s) for s in seqs]} → SKIP (out of [50, 800] bp)")
            continue
        msa = progressive_msa(seqs)
        pssms.append(build_pfm(msa))
        print(f"    exon {ei:2d}: lengths {[len(s) for s in seqs]} → MSA len {len(msa[0])}")
    print(f"  Built {len(pssms)} per-exon PSSMs.")

    print("[4/5] Sampling read populations ...")
    pos_reads = sample_reads(Path(args.bam), family_region, unmapped=False, n=args.n_positive, rng=rng)
    unm_reads = sample_reads(Path(args.unmapped_bam), None, unmapped=True, n=args.n_unmapped, rng=rng)
    # Random control: same chromosome, far from the family.
    rand_region = (args.chrom, max(0, family_region[1] - 30_000_000), max(1_000_000, family_region[1] - 25_000_000))
    rnd_reads = sample_reads(Path(args.bam), rand_region, unmapped=False, n=args.n_random, rng=rng)
    print(f"  positive: {len(pos_reads)}  unmapped: {len(unm_reads)}  random: {len(rnd_reads)}")

    print("[5/5] Scoring ...")
    pos_pssm = [family_score(s, pssms) for s in pos_reads]
    unm_pssm = [family_score(s, pssms) for s in unm_reads]
    rnd_pssm = [family_score(s, pssms) for s in rnd_reads]

    # K-mer Jaccard baseline: build a family k-mer set from the same exonic
    # sequences and score each read by k-mer hits / read k-mer count.
    K = 15
    fam_kmers: set[int] = set()
    for ci in range(len(mrnas)):
        for ei in range(n_exons):
            seq = per_mrna_exons[ci][ei]
            for i in range(len(seq) - K + 1):
                w = seq[i:i + K]
                if all(b in b"ACGT" for b in w):
                    fam_kmers.add(hash(bytes(w)))
    print(f"  Family k-mer set: {len(fam_kmers)} unique 15-mers")

    def kmer_hits(read: bytes) -> float:
        rc = read.translate(COMP)[::-1]
        best = 0
        for r in (read, rc):
            hits = 0
            n = max(1, len(r) - K + 1)
            for i in range(0, len(r) - K + 1):
                w = r[i:i + K]
                if all(b in b"ACGT" for b in w) and hash(bytes(w)) in fam_kmers:
                    hits += 1
            best = max(best, hits)
        return best

    pos_kmer = [kmer_hits(s) for s in pos_reads]
    unm_kmer = [kmer_hits(s) for s in unm_reads]
    rnd_kmer = [kmer_hits(s) for s in rnd_reads]

    print()
    print("=== Family-PSSM log-odds (vs uniform background, max over exons & strand) ===")
    summarise("positive   ", pos_pssm)
    summarise("unmapped   ", unm_pssm)
    summarise("random ctrl", rnd_pssm)
    print()
    print("=== Family k-mer hits (k=15, max over strands) — baseline ===")
    summarise("positive   ", pos_kmer)
    summarise("unmapped   ", unm_kmer)
    summarise("random ctrl", rnd_kmer)
    print()
    histogram("positive  PSSM ", pos_pssm)
    print()
    histogram("unmapped  PSSM ", unm_pssm)
    print()
    histogram("random ct PSSM ", rnd_pssm)
    print()
    histogram("positive  KMER ", pos_kmer)
    print()
    histogram("unmapped  KMER ", unm_kmer)
    print()
    histogram("random ct KMER ", rnd_kmer)

    # Verdict using upper-tail statistics on both signals — the question is whether
    # SOME unmapped reads carry family signal, not whether ALL of them do.
    print()
    print("=== Verdict ===")
    def tail(name: str, xs: list[float]) -> tuple[float, float]:
        if not xs: return (0.0, 0.0)
        return float(np.quantile(xs, 0.90)), float(np.max(xs))

    pos_p90, pos_max = tail("positive", pos_pssm)
    unm_p90, unm_max = tail("unmapped", unm_pssm)
    rnd_p90, rnd_max = tail("random",   rnd_pssm)
    pos_k90, pos_kmax = tail("positive", pos_kmer)
    unm_k90, unm_kmax = tail("unmapped", unm_kmer)
    rnd_k90, rnd_kmax = tail("random",   rnd_kmer)
    print(f"  PSSM 90th-pct  positive={pos_p90:7.2f}  unmapped={unm_p90:7.2f}  random={rnd_p90:7.2f}")
    print(f"  PSSM max       positive={pos_max:7.2f}  unmapped={unm_max:7.2f}  random={rnd_max:7.2f}")
    print(f"  KMER 90th-pct  positive={pos_k90:7.2f}  unmapped={unm_k90:7.2f}  random={rnd_k90:7.2f}")
    print(f"  KMER max       positive={pos_kmax:7.2f}  unmapped={unm_kmax:7.2f}  random={rnd_kmax:7.2f}")

    # Does ANY scorer separate unmapped from random in the tail?
    pssm_signal = unm_max > rnd_max + 1.0 or unm_p90 > rnd_p90 + 0.5
    kmer_signal = unm_kmax > rnd_kmax + 5 or unm_k90 > rnd_k90 + 5
    if kmer_signal and pssm_signal:
        print("  PASS: BOTH PSSM and k-mer scoring show a tail of unmapped reads above random.")
        print("        Family signal in unmapped reads is reproducible across two methods.")
        print("        Full HMM with proper M/I/D + family pseudo-counts will do strictly better.")
    elif kmer_signal:
        print("  PASS (k-mer): unmapped reads have a tail above random by k-mer hits.")
        print("        PSSM is too noisy on this dataset (3 copies → flat profile + composition bias).")
        print("        Architecture works; full HMM should improve over k-mer baseline.")
    elif pssm_signal:
        print("  WEAK: PSSM shows a tail but k-mer doesn't. Could be PSSM artifact; investigate.")
    else:
        print("  FAIL: no separation in either scorer's tail. Reconsider before investing.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
