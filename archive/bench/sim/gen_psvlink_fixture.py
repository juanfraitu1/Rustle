#!/usr/bin/env python3
"""
gen_psvlink_fixture.py — Deterministic PSV->junction LINKAGE fixture (Layer-2 "C").

Design
------
chrSIM_PSVLINK (~50 kb) contains two near-identical 3-exon paralog copies:
  Copy A : exons at [10000,10300), [10800,11100), [11700,12000) — 3×300 bp
  Copy B : exons at [30000,30300), [30800,31100), [31700,32000) — 3×300 bp

Both copies are on the '+' strand. Copy B's ENTIRE gene span (exons AND introns,
[gene_start, gene_end)) is copied verbatim from copy A, so minimap2 aligns reads from
both copies with BYTE-IDENTICAL CIGARs (identical splice boundaries) — this is what
makes the family graph's per-copy exon nodes EQUAL LENGTH, which the v1 PSV-column
extractor requires (it does no indel re-alignment). Copying introns too (not just
exons) is essential: divergent intron sequence lets minimap2 slide the splice site a
few bp differently per copy, jittering the exon node lengths and silently zeroing out
the PSV columns (the (E) gate then fails). With the whole span copied, the only
differences are the PSVs.

We then sprinkle PSVs (paralog-sequence variants) into copy B's exon1 and exon3 ONLY
(the two exons a middle-exon-skip isoform RETAINS) — 3 PSVs per exon, 6 total — kept
well INSIDE the exons (away from the splice boundaries) so they never perturb the
alignment. These 6 positions are the ONLY thing distinguishing copy A from copy B over
the skipped isoform's footprint. exon2 is left identical between the copies (it does
not matter for the skip isoform and keeps the demonstration focused on the retained
exons).

  - copy A allele at each PSV = the original (reference / pre-mutation) base.
  - copy B allele at each PSV = a deterministically mutated base.
So a read carrying COPY-A's sequence reads the copy-A allele at every PSV; a read
carrying copy-B's sequence reads the copy-B allele.

The exon-skip isoform (the recovery target)
-------------------------------------------
  COPY A expresses an exon-SKIP isoform (exon1 + exon3, SKIPPING exon2) in addition
  to its full 3-exon isoform.  COPY B expresses ONLY its full 3-exon isoform — it has
  no skip.  The skip reads use COPY A's sequence, so they carry copy-A's PSV alleles
  at the 6 retained-exon PSVs AND span the skip junction (exon1.end -> exon3.start).

Why this is resolvable ONLY via the linked PSV
----------------------------------------------
The two copies are near-identical, so copy-A reads cross-map to copy B and vice-versa
(minimap2 emits secondary alignments at the sibling locus) → a Layer-2 family forms
with the 6 PSVs as diagnostic columns (the (E) identifiability gate wants >= 3).

The skip junction, BY COORDINATES ALONE, is ambiguous: a 2-exon (exon1+exon3) chain
at copy A's coordinates looks structurally just like the same chain mirrored onto
copy B (the copies share intron/exon structure). Nothing in the splice graph says
which copy the skip belongs to. The ONLY positive evidence is WITHIN the molecule:
each skip read spans >= 3 copy-A PSV alleles AND the skip junction, so PSV->junction
linkage assigns that junction (and the assembled 2-exon isoform) to copy A. Suppress
the PSV channel (or set the per-read PSV floor absurdly high) and the skip vanishes —
no other path recovers it (verified with RUSTLE_VG_LAYER2_NO_MULTI_ISOFORM=1, which
isolates Part A's multi-isoform discovery so only the (C) PSV channel can emit it).

Read mix (deterministic; exact spliced transcript sequences, 0 errors, HiFi-like)
--------------------------------------------------------------------------------
  copy A full : 30 reads, copy-A sequence, full 3-exon            → primary at A
  copy A skip : 10 reads, copy-A sequence, exon1+exon3 (skip ex2) → primary at A
                (carry copy-A PSV alleles AND span the skip junction → link to A)
  copy B full : 30 reads, copy-B sequence, full 3-exon            → primary at B
Each copy's near-identical sibling produces SECONDARY alignments (cross-mapping) so
the Layer-2 family forms. >= min_linked (2) confidently-A-assigned skip reads link
the skip junction to copy A (10 skip reads ≫ 2).

Outputs (written to bench/fixtures/)
--------------------------------------
  sim_psvlink.fa        — synthetic genome FASTA (chrSIM_PSVLINK, 50 kb)
  sim_psvlink.reads.fa  — read FASTA (for provenance / re-alignment)
  sim_psvlink.bam       — minimap2 splice:hq alignment (sorted, indexed)
  sim_psvlink.bam.bai
  sim_psvlink.fa.fai    — produced by samtools faidx

Usage
-----
  python3 bench/sim/gen_psvlink_fixture.py --out-dir bench/fixtures
"""

import argparse
import os
import random
import subprocess
import sys

# ---------------------------------------------------------------------------
# Parameters (all deterministic / fixed seeds)
# ---------------------------------------------------------------------------
SEED          = 42
GENOME_LEN    = 50_000
CHROM         = "chrSIM_PSVLINK"

# Gene copy A
A_EXONS = [(10000, 10300), (10800, 11100), (11700, 12000)]  # [start, end) 0-based

# Gene copy B — same exon/intron structure, shifted +20 kb
B_OFFSET = 20_000
B_EXONS  = [(s + B_OFFSET, e + B_OFFSET) for s, e in A_EXONS]

# PSVs: 3 per RETAINED exon (exon1 and exon3 — the exons the skip isoform keeps),
# placed deterministically at 1/4, 1/2, 3/4 through each retained exon. exon2 (the
# SKIPPED exon, index 1) gets NO PSVs — the demonstration lives entirely on the
# retained-exon footprint that a skip read spans. 6 PSVs total → the (E) gate (>= 3
# diagnostic columns) passes and a single skip read spans all 6.
RETAINED_EXON_IDX = [0, 2]
PSVS_PER_EXON     = 3
PSV_SEED          = SEED + 99   # separate seed for PSV mutation

BASES = "ACGT"

# Read mix producing the PSV-linkage condition
A_FULL_READS = 30   # copy A full 3-exon isoform — primary at A
A_SKIP_READS = 10   # copy A exon-skip (exon1+exon3) — primary at A; THE recovery target
B_FULL_READS = 30   # copy B full 3-exon isoform — primary at B
B_SKIP_READS = 0    # copy B expresses NO skip isoform


def make_genome(seed: int) -> str:
    """Return a genome string of GENOME_LEN random bases (fixed seed)."""
    rng = random.Random(seed)
    return "".join(rng.choice(BASES) for _ in range(GENOME_LEN))


def build_psv_positions():
    """
    Deterministically pick PSV positions inside copy B's RETAINED exons (exon1, exon3)
    at 1/4, 1/2, 3/4 through each. Returns a list of 0-based genome positions in B.
    """
    positions = []
    for ei in RETAINED_EXON_IDX:
        bs, be = B_EXONS[ei]
        length = be - bs
        for k in range(1, PSVS_PER_EXON + 1):
            positions.append(bs + length * k // (PSVS_PER_EXON + 1))
    return positions


def apply_psvs(genome_str: str, psv_positions: list, seed: int) -> str:
    """Mutate specific positions to a different base (deterministic)."""
    psv_rng = random.Random(seed)
    g = list(genome_str)
    for pos in psv_positions:
        orig = g[pos]
        alts = [b for b in BASES if b != orig]
        g[pos] = psv_rng.choice(alts)
    return "".join(g)


def splice_transcript(genome: str, exons) -> str:
    """Concatenate exon sequences from genome (0-based half-open intervals)."""
    return "".join(genome[s:e] for s, e in exons)


def write_fasta(path: str, records: list):
    """Write [(name, seq), ...] to a FASTA file (60 bp/line)."""
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i+60] + "\n")


def run(cmd: str, check: bool = True) -> int:
    print(f"  $ {cmd}", file=sys.stderr)
    result = subprocess.run(cmd, shell=True, capture_output=False)
    if check and result.returncode != 0:
        raise RuntimeError(f"Command failed (exit {result.returncode}): {cmd}")
    return result.returncode


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default="bench/fixtures",
                    help="Directory to write output files into")
    args = ap.parse_args()

    out = args.out_dir
    os.makedirs(out, exist_ok=True)

    genome_fa = os.path.join(out, "sim_psvlink.fa")
    reads_fa  = os.path.join(out, "sim_psvlink.reads.fa")
    bam_path  = os.path.join(out, "sim_psvlink.bam")

    # -----------------------------------------------------------------------
    # 1.  Build genome
    # -----------------------------------------------------------------------
    print("Building genome...", file=sys.stderr)
    base_genome = make_genome(SEED)

    # Copy COPY A's ENTIRE gene span (exons AND introns) verbatim into copy B's span.
    # Copying the introns too — not just the exons — is what makes minimap2 produce
    # byte-identical CIGARs for both copies, so the family graph's per-copy exon nodes
    # are EQUAL LENGTH (the v1 PSV-column extractor requires equal length; divergent
    # introns would let the aligner slide the splice site a few bp per copy, jittering
    # node lengths and zeroing out the PSV columns).
    a_gene_start = A_EXONS[0][0]
    a_gene_end   = A_EXONS[-1][1]
    gene_len     = a_gene_end - a_gene_start
    g = list(base_genome)
    for i in range(gene_len):
        g[B_OFFSET + a_gene_start + i] = base_genome[a_gene_start + i]
    genome_pre_psv = "".join(g)  # copy A and copy B gene spans are byte-identical here

    # Apply PSVs ONLY to copy B's retained exons (exon1, exon3).
    psv_positions = build_psv_positions()
    final_genome  = apply_psvs(genome_pre_psv, psv_positions, PSV_SEED)

    # copy A uses the pre-PSV (reference) sequence; copy B uses the PSV-applied one.
    genome_a = genome_pre_psv
    genome_b = final_genome  # only differs from genome_a inside copy B's retained exons

    # Verify retained-exon PSV count over the skip footprint (exon1 + exon3).
    a_skip_seq = splice_transcript(genome_a, [A_EXONS[i] for i in RETAINED_EXON_IDX])
    b_skip_seq = splice_transcript(genome_b, [B_EXONS[i] for i in RETAINED_EXON_IDX])
    n_psv_skip = sum(a != b for a, b in zip(a_skip_seq, b_skip_seq))
    skip_ident = (len(a_skip_seq) - n_psv_skip) / len(a_skip_seq) * 100

    # Full-isoform identity (all 3 exons; exon2 is identical so this equals n_psv_skip).
    a_full_exonseq = splice_transcript(genome_a, A_EXONS)
    b_full_exonseq = splice_transcript(genome_b, B_EXONS)
    n_psv_full = sum(a != b for a, b in zip(a_full_exonseq, b_full_exonseq))
    full_ident = (len(a_full_exonseq) - n_psv_full) / len(a_full_exonseq) * 100

    write_fasta(genome_fa, [(CHROM, final_genome)])
    print(f"  Genome written: {genome_fa}  ({GENOME_LEN} bp)", file=sys.stderr)
    print(f"  Copy A exons (0-based): {A_EXONS}", file=sys.stderr)
    print(f"  Copy B exons (0-based): {B_EXONS}", file=sys.stderr)
    print(f"  Skip junction (copy A, 0-based): {A_EXONS[0][1]} -> {A_EXONS[2][0]}", file=sys.stderr)
    print(f"  PSV positions in B (0-based, retained exons only): {psv_positions}", file=sys.stderr)
    print(f"  PSVs over skip footprint (exon1+exon3): {n_psv_skip}  identity: {skip_ident:.2f}%",
          file=sys.stderr)
    print(f"  PSVs over full isoform (exon2 identical): {n_psv_full}  identity: {full_ident:.2f}%",
          file=sys.stderr)

    # -----------------------------------------------------------------------
    # 2.  Generate reads
    # -----------------------------------------------------------------------
    print("Generating reads...", file=sys.stderr)

    seq_a_full = splice_transcript(genome_a, A_EXONS)
    seq_a_skip = splice_transcript(genome_a, [A_EXONS[0], A_EXONS[2]])  # exon1 + exon3
    seq_b_full = splice_transcript(genome_b, B_EXONS)

    records = []
    # Copy A: 30 full + 10 skip (the skip is the recovery target).
    for i in range(A_FULL_READS):
        records.append((f"copyA_full_{i:03d}", seq_a_full))
    for i in range(A_SKIP_READS):
        records.append((f"copyA_skip_{i:03d}", seq_a_skip))
    # Copy B: 30 full, NO skip.
    for i in range(B_FULL_READS):
        records.append((f"copyB_full_{i:03d}", seq_b_full))

    write_fasta(reads_fa, records)
    total = A_FULL_READS + A_SKIP_READS + B_FULL_READS + B_SKIP_READS
    print(f"  Reads written: {reads_fa}  ({total} reads total)", file=sys.stderr)
    print(f"    copyA_full: {A_FULL_READS}  copyA_skip: {A_SKIP_READS}  copyB_full: {B_FULL_READS}",
          file=sys.stderr)

    # -----------------------------------------------------------------------
    # 3.  Index genome
    # -----------------------------------------------------------------------
    print("Indexing genome (samtools faidx)...", file=sys.stderr)
    run(f"samtools faidx {genome_fa}")

    # -----------------------------------------------------------------------
    # 4.  Align with minimap2 (-N 5 ≥ copy count → secondaries at the sibling locus)
    # -----------------------------------------------------------------------
    print("Aligning reads (minimap2 splice:hq)...", file=sys.stderr)
    mm2_cmd = (
        f"minimap2 -ax splice:hq -uf --secondary=yes -N 5 "
        f"--MD -Y "
        f"{genome_fa} {reads_fa} "
        f"| samtools sort -o {bam_path}"
    )
    run(mm2_cmd)
    run(f"samtools index {bam_path}")

    # -----------------------------------------------------------------------
    # 5.  Verify PSV-linkage conditions
    # -----------------------------------------------------------------------
    print("\n=== Fixture verification ===", file=sys.stderr)
    run(f"samtools idxstats {bam_path}")

    # The skip junction in CIGAR is an N gap from exon1.end to exon3.start. Its CIGAR
    # intron length = exon3.start - exon1.end. We detect the copy-A skip by an N whose
    # length matches AND whose alignment starts at copy A's exon1 (POS ~= 10001 1-based).
    a_skip_intron = A_EXONS[2][0] - A_EXONS[0][1]   # copy A skip junction intron length
    b_skip_intron = B_EXONS[2][0] - B_EXONS[0][1]   # would-be copy B skip intron length
    a_full_pos    = A_EXONS[0][0] + 1               # 1-based POS of a copy-A alignment
    b_full_pos    = B_EXONS[0][0] + 1               # 1-based POS of a copy-B alignment

    def view(args):
        r = subprocess.run(f"samtools view {args} {bam_path}",
                           shell=True, capture_output=True, text=True)
        return r.stdout

    # Count copy-A skip spliced alignments: POS at copy A's exon1 AND a single N of the
    # skip-junction length (one intron => 2-exon alignment). Primary alignments only.
    sam_primary = view("-F 0x100")
    a_skip_aln = 0
    a_full_aln = 0
    b_skip_aln = 0
    for line in sam_primary.strip().split("\n"):
        if not line:
            continue
        f = line.split("\t")
        if f[2] != CHROM:
            continue
        pos = int(f[3])
        cigar = f[5]
        n_introns = cigar.count("N")
        has_a_skip_N = f"{a_skip_intron}N" in cigar
        has_b_skip_N = f"{b_skip_intron}N" in cigar
        # copy A region: pos near a_full_pos (exon1 start of copy A)
        near_a = abs(pos - a_full_pos) <= 5
        near_b = abs(pos - b_full_pos) <= 5
        if near_a and n_introns == 1 and has_a_skip_N:
            a_skip_aln += 1
        elif near_a and n_introns == 2:
            a_full_aln += 1
        if near_b and n_introns == 1 and has_b_skip_N:
            b_skip_aln += 1

    print(f"  copy-A SKIP spliced alignments (POS~{a_full_pos}, 1 intron of {a_skip_intron}bp): {a_skip_aln}",
          file=sys.stderr)
    print(f"  copy-A FULL spliced alignments (POS~{a_full_pos}, 2 introns):                    {a_full_aln}",
          file=sys.stderr)
    print(f"  copy-B SKIP spliced alignments (POS~{b_full_pos}, 1 intron of {b_skip_intron}bp): {b_skip_aln}",
          file=sys.stderr)

    # PSV positions must actually differ between the copies (sanity on the genome).
    psv_diff_ok = all(genome_a[p - B_OFFSET] != genome_b[p] for p in psv_positions)

    ok = (a_skip_aln >= 2 and a_full_aln >= 1 and b_skip_aln == 0
          and n_psv_skip >= 3 and psv_diff_ok)
    if ok:
        print("\nPSV-LINKAGE STRUCTURE CONFIRMED:", file=sys.stderr)
        print(f"  copy A has {a_skip_aln} skip alignment(s) + {a_full_aln} full → both isoforms present at A",
              file=sys.stderr)
        print(f"  copy B has {b_skip_aln} skip alignment(s) → copy B expresses NO skip", file=sys.stderr)
        print(f"  {n_psv_skip} PSVs over the skip footprint distinguish the copies (alleles differ)",
              file=sys.stderr)
        print("  The copy-A skip junction is coordinate-ambiguous between the copies; only the", file=sys.stderr)
        print("  copy-A PSV alleles carried within the skip reads resolve it to copy A.", file=sys.stderr)
        print("  PASS", file=sys.stderr)
    else:
        print("\nWARNING: PSV-linkage structure may not be correct.", file=sys.stderr)
        print(f"  Expected: a_skip>=2, a_full>=1, b_skip==0, n_psv_skip>=3, psv_diff_ok=True", file=sys.stderr)
        print(f"  Got: a_skip={a_skip_aln} a_full={a_full_aln} b_skip={b_skip_aln} "
              f"n_psv_skip={n_psv_skip} psv_diff_ok={psv_diff_ok}", file=sys.stderr)
        print("  WARN", file=sys.stderr)

    print(f"\nFiles produced:", file=sys.stderr)
    for f in [genome_fa, genome_fa + ".fai", reads_fa, bam_path, bam_path + ".bai"]:
        size = os.path.getsize(f) if os.path.exists(f) else -1
        print(f"  {f}  ({size:,} bytes)", file=sys.stderr)


if __name__ == "__main__":
    main()
