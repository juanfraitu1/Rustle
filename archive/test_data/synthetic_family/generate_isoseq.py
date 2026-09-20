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
REF_HEADER = lines[0].strip()  # ">chr_test"
REF_SEQ = "".join(l.strip() for l in lines[1:])

# ── Gene coordinates (0-based, half-open) ──────────────────────────────────────
GENE_A_EXONS = [(1000, 1200), (2000, 2300), (3000, 3150), (4000, 4500)]
GENE_A_SKIP  = [(1000, 1200), (2000, 2300), (4000, 4500)]
GENE_B_EXONS = [(10000, 10200), (11000, 11300), (12000, 12150), (13000, 13500)]

# ── Planted diagnostic SNPs ────────────────────────────────────────────────────
# Copy A and copy B carry DIFFERENT bases at these corresponding exon-relative
# offsets. This is what makes the two paralog copies distinguishable by the
# fingerprint-EM: copy A genome base == base_A, copy B genome base == base_B,
# with base_A != base_B. A read's allele at these sites is therefore diagnostic
# of its copy of origin.
#
# base_A = the existing ref.fa base at copy_A_ref_pos (kept as-is, the "X" allele).
# base_B = a deterministically-chosen base != base_A (the "Y" allele) that we
#          plant into the copy-B region of the reference.
#
# Format: (copy_A_ref_pos, copy_B_ref_pos, base_A, base_B)  [bases filled below]
_SNP_POSITIONS = [
    # Exon 1: pos 50 within exon
    (1050, 10050),
    # Exon 2: pos 100, 200 within exon
    (2100, 11100),
    (2200, 11200),
    # Exon 3: pos 75 within exon
    (3075, 12075),
    # Exon 4: pos 100, 250, 400 within exon
    (4100, 13100),
    (4250, 13250),
    (4400, 13400),
]

def _pick_alt_base(base_a):
    """Deterministically pick copy B's base: the first ACGT base that differs
    from copy A's base. Guarantees base_A != base_B (a real substitution)."""
    return next(b for b in "ACGT" if b != base_a)

# Derive base_A (the genome's current base at posA) and base_B (the planted
# divergent base) for every diagnostic site.
DIAGNOSTIC_SNPS = []
for pos_a, pos_b in _SNP_POSITIONS:
    base_a = REF_SEQ[pos_a].upper()
    base_b = _pick_alt_base(base_a)
    DIAGNOSTIC_SNPS.append((pos_a, pos_b, base_a, base_b))

# Plant copy B's divergent allele into the per-copy reference. Copy A keeps its
# existing base (base_A); copy B's genome position posB now carries base_B.
# rustle reads ref.fa to build per-copy exon sequences, so this is what makes
# the fingerprint sites actually diagnostic.
_ref_list = list(REF_SEQ)
for pos_a, pos_b, base_a, base_b in DIAGNOSTIC_SNPS:
    _ref_list[pos_b] = base_b
REF_SEQ = "".join(_ref_list)

# Per-placement allele maps {ref_pos: base}. A read always carries its TRUE
# copy's allele at the diagnostic sites, regardless of which copy it is *placed*
# on. These are keyed by the *placement's* ref positions so apply_snps_to_seq
# can overwrite the spliced read sequence at the right spliced offset.
#
#   A_PLACEMENT_TRUE_A : read placed at copy A, true copy A → carries X at posA
#   B_PLACEMENT_TRUE_A : read placed at copy B, true copy A → carries X at posB
#                        (vs copy B's genome Y → a recorded mismatch = signal)
#   B_PLACEMENT_TRUE_B : read placed at copy B, true copy B → carries Y at posB
#   A_PLACEMENT_TRUE_B : read placed at copy A, true copy B → carries Y at posA
#                        (vs copy A's genome X → a recorded mismatch = signal)
A_PLACEMENT_TRUE_A = {pos_a: base_a for pos_a, _, base_a, _ in DIAGNOSTIC_SNPS}
B_PLACEMENT_TRUE_A = {pos_b: base_a for _, pos_b, base_a, _ in DIAGNOSTIC_SNPS}
B_PLACEMENT_TRUE_B = {pos_b: base_b for _, pos_b, _, base_b in DIAGNOSTIC_SNPS}
A_PLACEMENT_TRUE_B = {pos_a: base_b for pos_a, _, _, base_b in DIAGNOSTIC_SNPS}

def complement(base):
    return {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}[base]

def rev_comp(seq):
    return "".join(complement(b) for b in reversed(seq))

def extract_spliced_seq(exons):
    """Extract concatenated exon sequences from reference."""
    return "".join(REF_SEQ[s:e] for s, e in exons)

def apply_snps_to_seq(seq, exons, allele_map):
    """Overwrite the read's spliced sequence at diagnostic reference positions
    with a specific allele.

    `allele_map` is {ref_pos: base}. For every diagnostic ref_pos that falls
    inside one of `exons`, the corresponding base in the spliced read sequence
    is set to `base`. This is used to make a multi-mapper carry its TRUE copy's
    allele at the OTHER copy's placement: e.g. a true-copy-A read placed at copy
    B carries copy A's base there, which (vs copy B's divergent genome base)
    becomes a recorded mismatch — exactly the signal the fingerprint scores.
    """
    seq = list(seq)
    offset = 0
    for s, e in exons:
        for ref_pos, base in allele_map.items():
            if s <= ref_pos < e:
                seq[offset + (ref_pos - s)] = base
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
    # True copy A, placed at copy A: lock copy A's allele (X) at the diagnostic
    # sites so CCS error injection can't corrupt the signal.
    seq = apply_snps_to_seq(seq, exons, A_PLACEMENT_TRUE_A)

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
    # True copy A, placed at copy A: lock copy A's allele (X).
    seq = apply_snps_to_seq(seq, exons, A_PLACEMENT_TRUE_A)

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
    # Copy B's divergent alleles are baked into REF_SEQ (and ref.fa) above, so
    # the slice already carries Y. Add CCS errors, then re-lock Y at the
    # diagnostic sites so errors can't corrupt the signal.
    seq = add_isoseq_errors(ref_seq_for_read, error_rate=0.005)
    seq = apply_snps_to_seq(seq, exons, B_PLACEMENT_TRUE_B)

    polya = make_polya_tail()
    add_read(f"uniq_B1_{i}", 16, exons[0][0] + 1, exons, seq,
             clip_right_seq=polya, nh=1, ts="+")

# ── Multi-mapping reads ─────────────────────────────────────────────────────
# Use DETERMINISTIC, EVENLY-SPACED TSS positions to avoid collapsing with:
#   (a) unique reads (which have TSS at 1000±3 = 997-1003 for copy A, 10000±3 for B)
#   (b) each other
# Group A (primary at copy A): TSS at 1005, 1010, ..., 1070 (14 reads, 5-bp steps)
# Group B (primary at copy B): TSS at 10075, 10080, ..., 10140 (14 reads)
#   → supplementary at copy A: 1075, 1080, ..., 1140
# All TES fixed to avoid inter-read exonmatch collapse.
N_MULTI = 14

for i in range(N_MULTI):
    # Fixed TSS offset: 5 + i*5 from exon 1 start (1000), safely inside exon (1000-1200)
    tss_offset_a = 5 + i * 5  # 5, 10, 15, ..., 70
    tes_jitter_a = 0           # fixed TES to avoid exonmatch collapse

    exons_a = list(GENE_A_EXONS)
    exons_a[0] = (exons_a[0][0] + tss_offset_a, exons_a[0][1])
    # TES fixed (no jitter) so no two multi-mappers have the same exon structure
    seq_a = add_isoseq_errors(extract_spliced_seq(exons_a), error_rate=0.005)
    # TRUE copy A, primary placement at copy A: carry copy A's allele (X).
    seq_a = apply_snps_to_seq(seq_a, exons_a, A_PLACEMENT_TRUE_A)
    polya = make_polya_tail()

    add_read(f"multi_{i}", 16, exons_a[0][0] + 1, exons_a, seq_a,
             clip_right_seq=polya, nh=2, ts="+")

    # Supplementary at copy B: offset by 9000
    exons_b = list(GENE_B_EXONS)
    exons_b[0] = (exons_b[0][0] + tss_offset_a, exons_b[0][1])
    seq_b = add_isoseq_errors(extract_spliced_seq(exons_b), error_rate=0.005)
    # TRUE copy A, but placed at copy B: the read STILL carries copy A's allele
    # (X) at copy B's diagnostic positions. vs copy B's genome base (Y) this is a
    # recorded mismatch, so the fingerprint scores X -> matches A, mismatches B.
    seq_b = apply_snps_to_seq(seq_b, exons_b, B_PLACEMENT_TRUE_A)

    add_read(f"multi_{i}", 2064, exons_b[0][0] + 1, exons_b, seq_b,
             clip_right_seq=polya, nh=2, ts="+")

for i in range(N_MULTI):
    # TSS offset for group B: 75 + i*5 from exon 1 start, no overlap with group A (5-70)
    tss_offset_b = 75 + i * 5  # 75, 80, 85, ..., 140

    exons_b = list(GENE_B_EXONS)
    exons_b[0] = (exons_b[0][0] + tss_offset_b, exons_b[0][1])
    seq_b = add_isoseq_errors(extract_spliced_seq(exons_b), error_rate=0.005)
    # TRUE copy B, primary placement at copy B: carry copy B's allele (Y).
    seq_b = apply_snps_to_seq(seq_b, exons_b, B_PLACEMENT_TRUE_B)
    polya = make_polya_tail()

    add_read(f"multi_r_{i}", 16, exons_b[0][0] + 1, exons_b, seq_b,
             clip_right_seq=polya, nh=2, ts="+")

    exons_a = list(GENE_A_EXONS)
    exons_a[0] = (exons_a[0][0] + tss_offset_b, exons_a[0][1])
    seq_a = add_isoseq_errors(extract_spliced_seq(exons_a), error_rate=0.005)
    # TRUE copy B, but placed at copy A: the read STILL carries copy B's allele
    # (Y) at copy A's diagnostic positions. vs copy A's genome base (X) this is a
    # recorded mismatch, so the fingerprint scores Y -> matches B, mismatches A.
    seq_a = apply_snps_to_seq(seq_a, exons_a, A_PLACEMENT_TRUE_B)

    add_read(f"multi_r_{i}", 2064, exons_a[0][0] + 1, exons_a, seq_a,
             clip_right_seq=polya, nh=2, ts="+")

# ── Rewrite reference with the planted copy-B divergence ──────────────────────
# REF_SEQ now carries copy B's divergent bases (Y) at the posB diagnostic sites
# while copy A keeps its bases (X). Persist it to ref.fa so rustle's per-copy
# exon sequences (and therefore the fingerprint sites) actually differ.
output_dir = os.path.dirname(__file__)
with open(REF_PATH, "w") as f:
    f.write(REF_HEADER + "\n")
    for i in range(0, len(REF_SEQ), 80):
        f.write(REF_SEQ[i:i+80] + "\n")
# Rebuild the .fai so coordinates stay authoritative.
subprocess.run(["samtools", "faidx", REF_PATH], check=True)
print(f"Rewrote reference with {len(DIAGNOSTIC_SNPS)} planted copy-B SNPs: {REF_PATH}")

# ── Write SAM and convert to BAM ──────────────────────────────────────────────
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
print(f"Multi-mapping: {N_MULTI} (A→B) + {N_MULTI} (B→A) = {2*N_MULTI} reads, {4*N_MULTI} alignments")
