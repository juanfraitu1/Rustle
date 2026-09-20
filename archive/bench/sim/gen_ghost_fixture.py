#!/usr/bin/env python3
"""
gen_ghost_fixture.py — Deterministic GHOST-copy paralog fixture for M5 recovery.

Design
------
chrSIM_GHOST (~50 kb) contains two homologous 3-exon genes:
  Copy A ("dominant")  : exons at [10000,10300), [10800,11100), [11700,12000) — 3×300 bp
  Copy B ("ghost")     : exons at [30000,30300), [30800,31100), [31700,32000) — 3×300 bp

Both copies are on the '+' strand.  The copies share ~99.3% exon identity — only 6
PSVs (two per exon, placed deterministically) distinguish copy B from copy A.  This
high similarity means every copy-B-derived read aligns with its PRIMARY to copy A
(the same splice graph, lower genomic coordinates win minimap2's tie-break) and only
a SECONDARY alignment lands at copy B.

GHOST mechanism (distinct from the starved fixture)
---------------------------------------------------
  Copy A: 40 full-isoform reads + 10 exon-skip reads → 50 primary reads at A.
  Copy B: ZERO primary reads.  There is NO read carrying copy B's exact sequence,
          so nothing ever aligns its PRIMARY to copy B.
          8 cross-map reads use copy-A's sequence → their PRIMARY lands at copy A
          (lower coords win the tie), producing a SECONDARY alignment at copy B.

Result
------
  copy A region [10000-12000]: 50 primary reads → well-assembled.
  copy B region [30000-32000]:  0 primary reads → Layer 1 builds NO bundle here.
                                8 secondary reads → these get locus=None (there is no
                                  bundle to attach them to) and are dropped by Layer 1.

Why this needs M5 (all-secondary new-copy discovery), NOT homology-transfer
---------------------------------------------------------------------------
In the `sim_starved` fixture copy B has 1 primary read, so a bundle forms, the
side-index attaches the secondaries to that locus, and the existing Layer-2
homology coordinate-transfer path reconstructs copy B.  Here copy B has ZERO
primaries: no bundle, no locus, the secondaries are `locus=None`.  The ONLY way to
recover copy B is M5 (--vg-layer2-new-copies): cluster the locus=None secondaries
(8 alignments ≥ 3 distinct reads, each with 1 placement ≤ max_placements=3, tracing
B's 3-exon chain whose terminal exons ~300bp ≫ the 25bp min-exon gate) and emit the
cluster's consensus chain as a recovered new copy.

  - default --vg / --vg-layer2 (M5 OFF): copy B ABSENT (no bundle ever formed).
  - --vg-layer2-new-copies (M5 ON): copy B recovered as a multi-exon transcript.

Outputs (written to bench/fixtures/)
--------------------------------------
  sim_ghost.fa        — synthetic genome FASTA (chrSIM_GHOST, 50 kb)
  sim_ghost.reads.fa  — read FASTA (for provenance / re-alignment)
  sim_ghost.bam       — minimap2 splice:hq alignment (sorted, indexed)
  sim_ghost.bam.bai
  sim_ghost.fa.fai    — produced by samtools faidx

Usage
-----
  python3 bench/sim/gen_ghost_fixture.py --out-dir bench/fixtures
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
CHROM         = "chrSIM_GHOST"

# Gene copy A
A_EXONS = [(10000, 10300), (10800, 11100), (11700, 12000)]  # [start, end) 0-based

# Gene copy B — same exon/intron structure, shifted +20 kb
B_OFFSET = 20_000
B_EXONS  = [(s + B_OFFSET, e + B_OFFSET) for s, e in A_EXONS]

# PSVs: 2 per exon, placed at 1/3 and 1/2 of each exon
# These 6 positions distinguish copy B from copy A in minimap2 alignments.
# They are ONLY applied to copy B's exons.
PSV_SEED = SEED + 99   # separate seed for PSV mutation

BASES = "ACGT"

# Read mix producing the GHOST condition
A_FULL_READS    = 40   # copy A full 3-exon isoform — primary at A
A_SKIP_READS    = 10   # copy A exon-skip (exon1+exon3) — primary at A
B_PRIMARY_READS = 0    # copy B has ZERO primary reads — NO bundle forms at B
B_XMAP_READS    = 8    # copy-B-expression reads using copy-A sequence — primary at A
                       #   (secondary at B, locus=None: M5 clusters these to recover B)


def make_genome(seed: int) -> str:
    """Return a genome string of GENOME_LEN random bases (fixed seed)."""
    rng = random.Random(seed)
    return "".join(rng.choice(BASES) for _ in range(GENOME_LEN))


def build_psv_positions():
    """
    Deterministically pick 2 PSV positions per exon (1/3 and 1/2 through each exon).
    Returns a list of 0-based genome positions within B_EXONS.
    """
    positions = []
    for bs, be in B_EXONS:
        length = be - bs
        positions.append(bs + length // 3)
        positions.append(bs + length // 2)
    return positions


def apply_psvs(genome_str: str, psv_positions: list, seed: int) -> str:
    """
    Mutate specific positions to a different base (deterministic).
    Returns the modified genome string.
    """
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

    genome_fa = os.path.join(out, "sim_ghost.fa")
    reads_fa  = os.path.join(out, "sim_ghost.reads.fa")
    bam_path  = os.path.join(out, "sim_ghost.bam")

    # -----------------------------------------------------------------------
    # 1.  Build genome
    # -----------------------------------------------------------------------
    print("Building genome...", file=sys.stderr)
    base_genome = make_genome(SEED)

    # Make copy B's exon positions identical to copy A (pre-PSV)
    g = list(base_genome)
    for (as_, ae), (bs, be) in zip(A_EXONS, B_EXONS):
        a_seq = base_genome[as_:ae]
        for i, ch in enumerate(a_seq):
            g[bs + i] = ch
    genome_pre_psv = "".join(g)  # A and B exons are identical here

    # Apply PSVs only to copy B's exons (6 positions, 2 per exon)
    psv_positions = build_psv_positions()
    final_genome  = apply_psvs(genome_pre_psv, psv_positions, PSV_SEED)

    # Verify: A exon seq vs B exon seq
    a_exon_seq = splice_transcript(genome_pre_psv, A_EXONS)  # pre-PSV = A
    b_exon_seq  = splice_transcript(final_genome, B_EXONS)   # post-PSV = B
    n_psv_in_tx = sum(a != b for a, b in zip(a_exon_seq, b_exon_seq))
    identity    = (len(a_exon_seq) - n_psv_in_tx) / len(a_exon_seq) * 100

    write_fasta(genome_fa, [(CHROM, final_genome)])
    print(f"  Genome written: {genome_fa}  ({GENOME_LEN} bp)", file=sys.stderr)
    print(f"  Copy A exons (0-based): {A_EXONS}", file=sys.stderr)
    print(f"  Copy B exons (0-based): {B_EXONS}", file=sys.stderr)
    print(f"  PSV positions in B (0-based): {psv_positions}", file=sys.stderr)
    print(f"  PSVs in transcript: {n_psv_in_tx}  identity: {identity:.2f}%", file=sys.stderr)

    # -----------------------------------------------------------------------
    # 2.  Generate reads
    # -----------------------------------------------------------------------
    print("Generating reads...", file=sys.stderr)

    seq_a_full = splice_transcript(final_genome, A_EXONS)
    seq_a_skip = splice_transcript(final_genome, [A_EXONS[0], A_EXONS[2]])
    seq_b_full = splice_transcript(final_genome, B_EXONS)  # has PSVs

    a_b_ident_tx = sum(a == b for a, b in zip(seq_a_full, seq_b_full)) / len(seq_a_full) * 100
    print(f"  A_full vs B_full transcript identity: {a_b_ident_tx:.2f}%", file=sys.stderr)

    records = []

    # Copy A: 40 full + 10 skip
    for i in range(A_FULL_READS):
        records.append((f"copyA_full_{i:03d}", seq_a_full))
    for i in range(A_SKIP_READS):
        records.append((f"copyA_skip_{i:03d}", seq_a_skip))

    # Copy B primary: ZERO reads — no read carries copy B's exact sequence, so
    # nothing ever aligns its primary to copy B → Layer 1 builds NO bundle at B.
    for i in range(B_PRIMARY_READS):
        records.append((f"copyB_primary_{i:03d}", seq_b_full))

    # Copy B cross-map: 8 reads using A sequence → primary at A, secondary at B.
    # These represent copy-B expression that cross-maps because the copies are
    # near-identical. Their secondaries at copy B carry locus=None (no bundle to
    # attach to); M5 clusters them to recover the ghost copy.
    for i in range(B_XMAP_READS):
        records.append((f"copyB_xmap_{i:03d}", seq_a_full))

    write_fasta(reads_fa, records)
    total = A_FULL_READS + A_SKIP_READS + B_PRIMARY_READS + B_XMAP_READS
    print(f"  Reads written: {reads_fa}  ({total} reads total)", file=sys.stderr)
    print(f"    copyA_full: {A_FULL_READS}  copyA_skip: {A_SKIP_READS}", file=sys.stderr)
    print(f"    copyB_primary: {B_PRIMARY_READS}  copyB_xmap: {B_XMAP_READS}", file=sys.stderr)

    # -----------------------------------------------------------------------
    # 3.  Index genome
    # -----------------------------------------------------------------------
    print("Indexing genome (samtools faidx)...", file=sys.stderr)
    run(f"samtools faidx {genome_fa}")

    # -----------------------------------------------------------------------
    # 4.  Align with minimap2
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
    # 5.  Verify GHOST conditions
    # -----------------------------------------------------------------------
    print("\n=== Fixture verification ===", file=sys.stderr)
    run(f"samtools idxstats {bam_path}")

    def count_reads(flag_filter, region):
        r = subprocess.run(
            f"samtools view -c {flag_filter} {bam_path} '{CHROM}:{region}'",
            shell=True, capture_output=True, text=True
        )
        return r.stdout.strip()

    prim_a  = count_reads("-F 0x100", "10000-12000")
    prim_b  = count_reads("-F 0x100", "30000-32000")
    sec_a   = count_reads("-f 0x100", "10000-12000")
    sec_b   = count_reads("-f 0x100", "30000-32000")
    tot_sec = subprocess.run(
        f"samtools view -c -f 0x100 {bam_path}",
        shell=True, capture_output=True, text=True
    ).stdout.strip()

    print(f"Primary alignments  — copy A [10000-12000]: {prim_a}", file=sys.stderr)
    print(f"Primary alignments  — copy B [30000-32000]: {prim_b}", file=sys.stderr)
    print(f"Secondary alignments — copy A region:       {sec_a}", file=sys.stderr)
    print(f"Secondary alignments — copy B region:       {sec_b}", file=sys.stderr)
    print(f"Total secondary alignments:                  {tot_sec}", file=sys.stderr)

    # GHOST criteria check: ZERO primaries at copy B (no bundle), >= 3 secondaries
    # at copy B (M5's distinct-read floor is 3).
    prim_b_int = int(prim_b) if prim_b.isdigit() else -1
    sec_b_int  = int(sec_b) if sec_b.isdigit() else -1

    if prim_b_int == 0 and sec_b_int >= 3:
        print("\nGHOST STRUCTURE CONFIRMED:", file=sys.stderr)
        print(f"  copy B has {prim_b_int} primary read(s) → NO bundle forms (Layer 1 builds nothing at B)", file=sys.stderr)
        print(f"  copy B has {sec_b_int} secondary reads → locus=None; M5 clusters them to recover the copy", file=sys.stderr)
        print("  Run rustle --vg-layer2 (M5 OFF) → copy B ABSENT; add --vg-layer2-new-copies → recovered.", file=sys.stderr)
    else:
        print("\nWARNING: ghost structure may not be correct.", file=sys.stderr)
        print(f"  Expected: prim_b == 0, sec_b >= 3; got prim_b={prim_b}, sec_b={sec_b}", file=sys.stderr)

    # Spliced alignment check
    splice_check = subprocess.run(
        f"samtools view {bam_path} | awk '$6 ~ /N/' | head -3",
        shell=True, capture_output=True, text=True
    )
    if splice_check.stdout.strip():
        print("\nSpliced alignments confirmed (CIGAR with N):", file=sys.stderr)
        for line in splice_check.stdout.strip().split("\n"):
            fields = line.split("\t")
            print(f"  {fields[0]}  FLAG={fields[1]}  POS={fields[3]}  CIGAR={fields[5]}", file=sys.stderr)
    else:
        print("\nWARNING: no spliced alignments detected!", file=sys.stderr)

    print(f"\nFiles produced:", file=sys.stderr)
    for f in [genome_fa, genome_fa + ".fai", reads_fa, bam_path, bam_path + ".bai"]:
        size = os.path.getsize(f) if os.path.exists(f) else -1
        print(f"  {f}  ({size:,} bytes)", file=sys.stderr)


if __name__ == "__main__":
    main()
