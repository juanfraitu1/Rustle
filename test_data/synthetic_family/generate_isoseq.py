#!/usr/bin/env python3
"""Generate realistic IsoSeq-like synthetic reads for two paralogous gene copies.

Creates a BAM file that simulates PacBio IsoSeq CCS reads mapping to two gene copies
(paralogs) with ~98% sequence identity, different isoform expression profiles, and
multi-mapping reads between copies.

Gene family layout on chr_test (20kb reference):
  Copy A: exons at 1000-1200, 2000-2300, 3000-3150, 4000-4500
  Copy B: exons at 10000-10200, 11000-11300, 12000-12150, 13000-13500

Copy B differs from A by ~2% (diagnostic SNPs planted in exons).

Transcripts:
  A1: full 4-exon (major at copy A)
  A2: exon-3-skip (minor at copy A)
  B1: full 4-exon (only at copy B)

IsoSeq characteristics simulated:
  - CCS-quality sequences (Q20, ~1% residual error)
  - Poly-A tails (20-40bp, soft-clipped)
  - Variable TSS/TES (±5bp jitter at terminal exons)
  - ts:A tag for splice strand
  - Proper NH tags for multi-mappers
  - MD tags for mismatch tracking
"""

import random
import subprocess
import sys
import os

random.seed(42)

# ── Reference ──────────────────────────────────────────────────────────────────
REF_PATH = os.path.join(os.path.dirname(__file__), "ref.fa")
with open(REF_PATH) as f:
    lines = f.readlines()
REF_SEQ = "".join(l.strip() for l in lines[1:])

# ── Gene coordinates (0-based, half-open) ──────────────────────────────────────
GENE_A_EXONS = [(1000, 1200), (2000, 2300), (3000, 3150), (4000, 4500)]
GENE_A_SKIP  = [(1000, 1200), (2000, 2300), (4000, 4500)]
GENE_B_EXONS = [(10000, 10200), (11000, 11300), (12000, 12150), (13000, 13500)]

# ── Planted diagnostic SNPs ────────────────────────────────────────────────────
# Copy B has these substitutions relative to copy A (at corresponding positions).
# These are at fixed offsets within exons to help test --vg-snp mode.
# Format: (copy_A_ref_pos, copy_B_ref_pos, ref_base, alt_base)
DIAGNOSTIC_SNPS = [
    # Exon 1: pos 50 within exon
    (1050, 10050, None, None),
    # Exon 2: pos 100, 200 within exon
    (2100, 11100, None, None),
    (2200, 11200, None, None),
    # Exon 3: pos 75 within exon
    (3075, 12075, None, None),
    # Exon 4: pos 100, 250, 400 within exon
    (4100, 13100, None, None),
    (4250, 13250, None, None),
    (4400, 13400, None, None),
]

def complement(base):
    return {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}[base]

def rev_comp(seq):
    return "".join(complement(b) for b in reversed(seq))

def extract_spliced_seq(exons):
    """Extract concatenated exon sequences from reference."""
    return "".join(REF_SEQ[s:e] for s, e in exons)

def apply_snps_to_seq(seq, exons, snp_positions_in_ref):
    """Apply diagnostic SNP at reference positions that fall within these exons."""
    seq = list(seq)
    offset = 0
    for s, e in exons:
        for snp_ref_pos in snp_positions_in_ref:
            if s <= snp_ref_pos < e:
                idx = offset + (snp_ref_pos - s)
                orig = seq[idx]
                # Deterministic substitution
                alts = [b for b in "ACGT" if b != orig]
                seq[idx] = alts[snp_ref_pos % 3]
        offset += (e - s)
    return "".join(seq)

def add_isoseq_errors(seq, error_rate=0.01):
    """Add IsoSeq CCS-level errors (~1% residual after CCS correction)."""
    seq = list(seq)
    for i in range(len(seq)):
        if random.random() < error_rate:
            error_type = random.random()
            if error_type < 0.5:
                # Substitution
                alts = [b for b in "ACGT" if b != seq[i]]
                seq[i] = random.choice(alts)
            elif error_type < 0.8:
                # Insertion (we skip for CIGAR simplicity — just substitute)
                alts = [b for b in "ACGT" if b != seq[i]]
                seq[i] = random.choice(alts)
            else:
                # Deletion represented as substitution for simplicity
                alts = [b for b in "ACGT" if b != seq[i]]
                seq[i] = random.choice(alts)
    return "".join(seq)

def make_polya_tail(length=None):
    """Generate a poly-A tail (20-40bp) with occasional non-A bases."""
    if length is None:
        length = random.randint(20, 40)
    tail = []
    for _ in range(length):
        if random.random() < 0.05:
            tail.append(random.choice("CGT"))
        else:
            tail.append("A")
    return "".join(tail)

def make_polyt_prefix(length=None):
    """Generate a poly-T prefix (reverse complement of poly-A)."""
    if length is None:
        length = random.randint(20, 40)
    prefix = []
    for _ in range(length):
        if random.random() < 0.05:
            prefix.append(random.choice("ACG"))
        else:
            prefix.append("T")
    return "".join(prefix)

def jitter_terminal(pos, max_jitter=5, direction="start"):
    """Add small jitter to terminal exon boundaries (TSS/TES variation)."""
    delta = random.randint(-max_jitter, max_jitter)
    return max(0, pos + delta)

def build_cigar(exons, clip_left=0, clip_right=0):
    """Build CIGAR string from exon coordinates."""
    parts = []
    if clip_left > 0:
        parts.append(f"{clip_left}S")
    for i, (s, e) in enumerate(exons):
        parts.append(f"{e - s}M")
        if i < len(exons) - 1:
            intron = exons[i + 1][0] - e
            parts.append(f"{intron}N")
    if clip_right > 0:
        parts.append(f"{clip_right}S")
    return "".join(parts)

def compute_md(seq_aligned, ref_seq_at_exons):
    """Compute MD tag from aligned sequence vs reference."""
    md_parts = []
    match_count = 0
    for i in range(min(len(seq_aligned), len(ref_seq_at_exons))):
        if seq_aligned[i] == ref_seq_at_exons[i]:
            match_count += 1
        else:
            md_parts.append(str(match_count))
            md_parts.append(ref_seq_at_exons[i])
            match_count = 0
    md_parts.append(str(match_count))
    return "".join(md_parts)

def count_mismatches(seq, ref):
    return sum(1 for a, b in zip(seq, ref) if a != b)

# ── Generate reads ─────────────────────────────────────────────────────────────
sam_lines = []

def add_read(name, flag, pos_1based, exons, seq, clip_left_seq="", clip_right_seq="",
             nh=1, ts="+", mate_name="*", mate_pos=0):
    """Add a SAM record."""
    cigar = build_cigar(exons, len(clip_left_seq), len(clip_right_seq))
    full_seq = clip_left_seq + seq + clip_right_seq
    qual = "I" * len(full_seq)
    mapq = 60 if nh == 1 else 3

    # Compute MD and NM from the aligned portion
    ref_at_exons = extract_spliced_seq(exons)
    md = compute_md(seq, ref_at_exons)
    nm = count_mismatches(seq, ref_at_exons[:len(seq)])

    tags = f"NH:i:{nh}\tts:A:{ts}\tNM:i:{nm}\tMD:Z:{md}"
    sam_lines.append(
        f"{name}\t{flag}\tchr_test\t{pos_1based}\t{mapq}\t{cigar}\t{mate_name}\t{mate_pos}\t0\t{full_seq}\t{qual}\t{tags}"
    )

# Copy A SNP positions (for modifying copy-B sequences to look like copy-A alignments)
snp_a_positions = [s[0] for s in DIAGNOSTIC_SNPS]
snp_b_positions = [s[1] for s in DIAGNOSTIC_SNPS]

# ── Unique reads at Copy A ─────────────────────────────────────────────────────

# A1: full 4-exon transcript, 30 reads
for i in range(30):
    exons = list(GENE_A_EXONS)
    # TSS/TES jitter on terminal exons
    exons[0] = (jitter_terminal(exons[0][0], 3), exons[0][1])
    exons[-1] = (exons[-1][0], jitter_terminal(exons[-1][1], 3))

    ref_seq_for_read = extract_spliced_seq(exons)
    seq = add_isoseq_errors(ref_seq_for_read, error_rate=0.005)

    # Add poly-A tail (soft-clipped)
    polya = make_polya_tail()

    # Flag 16 = reverse strand (IsoSeq maps to - strand for this gene)
    add_read(f"uniq_A1_{i}", 16, exons[0][0] + 1, exons, seq,
             clip_right_seq=polya, nh=1, ts="+")

# A2: exon-skip (skip exon 3), 20 reads
for i in range(20):
    exons = list(GENE_A_SKIP)
    exons[0] = (jitter_terminal(exons[0][0], 3), exons[0][1])
    exons[-1] = (exons[-1][0], jitter_terminal(exons[-1][1], 3))

    ref_seq_for_read = extract_spliced_seq(exons)
    seq = add_isoseq_errors(ref_seq_for_read, error_rate=0.005)

    polya = make_polya_tail()
    add_read(f"uniq_A2_{i}", 16, exons[0][0] + 1, exons, seq,
             clip_right_seq=polya, nh=1, ts="+")

# ── Unique reads at Copy B ─────────────────────────────────────────────────────

# B1: full 4-exon, 25 reads (with diagnostic SNPs)
for i in range(25):
    exons = list(GENE_B_EXONS)
    exons[0] = (jitter_terminal(exons[0][0], 3), exons[0][1])
    exons[-1] = (exons[-1][0], jitter_terminal(exons[-1][1], 3))

    ref_seq_for_read = extract_spliced_seq(exons)
    # Copy B has diagnostic SNPs already baked into the reference
    # (the reference itself encodes copy B differently)
    seq = add_isoseq_errors(ref_seq_for_read, error_rate=0.005)

    polya = make_polya_tail()
    add_read(f"uniq_B1_{i}", 16, exons[0][0] + 1, exons, seq,
             clip_right_seq=polya, nh=1, ts="+")

# ── Multi-mapping reads (primary at A, supplementary at B) ────────────────────
# These represent reads from a highly similar paralog that the aligner can't
# confidently assign. They map equally well to both copies.

for i in range(6):
    # Primary at copy A
    exons_a = list(GENE_A_EXONS)
    exons_a[0] = (jitter_terminal(exons_a[0][0], 2), exons_a[0][1])
    exons_a[-1] = (exons_a[-1][0], jitter_terminal(exons_a[-1][1], 2))
    seq_a = add_isoseq_errors(extract_spliced_seq(exons_a), error_rate=0.005)
    polya = make_polya_tail()

    add_read(f"multi_{i}", 16, exons_a[0][0] + 1, exons_a, seq_a,
             clip_right_seq=polya, nh=2, ts="+")

    # Supplementary at copy B (same read, different mapping)
    exons_b = list(GENE_B_EXONS)
    exons_b[0] = (exons_a[0][0] + 9000, exons_b[0][1])
    exons_b[-1] = (exons_b[-1][0], exons_a[-1][1] + 9000)
    seq_b = add_isoseq_errors(extract_spliced_seq(exons_b), error_rate=0.005)

    add_read(f"multi_{i}", 2064, exons_b[0][0] + 1, exons_b, seq_b,
             clip_right_seq=polya, nh=2, ts="+")

# ── Multi-mapping reads (primary at B, supplementary at A) ────────────────────

for i in range(4):
    exons_b = list(GENE_B_EXONS)
    exons_b[0] = (jitter_terminal(exons_b[0][0], 2), exons_b[0][1])
    exons_b[-1] = (exons_b[-1][0], jitter_terminal(exons_b[-1][1], 2))
    seq_b = add_isoseq_errors(extract_spliced_seq(exons_b), error_rate=0.005)
    polya = make_polya_tail()

    add_read(f"multi_r_{i}", 16, exons_b[0][0] + 1, exons_b, seq_b,
             clip_right_seq=polya, nh=2, ts="+")

    exons_a = list(GENE_A_EXONS)
    exons_a[0] = (exons_b[0][0] - 9000, exons_a[0][1])
    exons_a[-1] = (exons_a[-1][0], exons_b[-1][1] - 9000)
    seq_a = add_isoseq_errors(extract_spliced_seq(exons_a), error_rate=0.005)

    add_read(f"multi_r_{i}", 2064, exons_a[0][0] + 1, exons_a, seq_a,
             clip_right_seq=polya, nh=2, ts="+")

# ── Write SAM and convert to BAM ──────────────────────────────────────────────
output_dir = os.path.dirname(__file__)
sam_path = os.path.join(output_dir, "reads.sam")
bam_path = os.path.join(output_dir, "reads_sorted.bam")

header = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr_test\tLN:20000\n"
with open(sam_path, "w") as f:
    f.write(header)
    for line in sam_lines:
        f.write(line + "\n")

print(f"Wrote {len(sam_lines)} alignments to {sam_path}")

# Sort and index
subprocess.run(["samtools", "sort", sam_path, "-o", bam_path], check=True)
subprocess.run(["samtools", "index", bam_path], check=True)
print(f"Sorted BAM: {bam_path}")

# Summary stats
primary = sum(1 for l in sam_lines if "\t16\t" in l or "\t0\t" in l)
supp = sum(1 for l in sam_lines if "\t2064\t" in l or "\t2048\t" in l)
print(f"Primary alignments: {primary}, Supplementary: {supp}")
print(f"Copy A unique: 30 (A1) + 20 (A2) = 50")
print(f"Copy B unique: 25 (B1)")
print(f"Multi-mapping: 6 (A→B) + 4 (B→A) = 10 reads, 20 alignments")
