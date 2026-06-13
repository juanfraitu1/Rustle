#!/usr/bin/env python3
"""
gen_starved_fixture.py — Deterministic genuine-starvation paralog fixture.

Design
------
chrSIM (~50 kb) contains two homologous 3-exon genes:
  Copy A ("dominant")  : exons at [10000,10300), [10800,11100), [11700,12000) — 3×300 bp
  Copy B ("starved")   : exons at [30000,30300), [30800,31100), [31700,32000) — 3×300 bp

Both copies are on the '+' strand.  The copies share ~99.3% exon identity — only 6
PSVs (two per exon, placed deterministically) distinguish copy B from copy A.  This
high similarity means most copy-B-derived reads align with their PRIMARY to copy A
(the same splice graph, lower genomic coordinates win minimap2 tie-break).

Starvation mechanism
--------------------
  Copy A: 40 full-isoform reads + 10 exon-skip reads → 50 primary reads at A.
  Copy B: 1 read with the EXACT copy-B transcript sequence → primaries at copy B
           (the 6 PSVs give it a slight alignment advantage at B over A, so it
           stays at B as primary).
           20 reads with the copy-A sequence → their primary lands at copy A
           (these are "copy-B-expression" reads that cross-map, stealing B's
           signal to A).  They produce SECONDARY alignments at copy B.

Result
------
  copy A region [10000-12000]: 70 primary reads → well-assembled.
  copy B region [30000-32000]: 1 primary read  → bundle forms, NO transcript emitted
                                70 secondary reads → side-index signals copy B's
                                                     coordinates for Layer-2 recovery.

This is the controlled case for Layer-2 homology coordinate-transfer:
  - Bundle at copy B materialises (1 primary read is enough for the bundle to form
    and for the side-index to attach secondaries to that locus).
  - No transcript emitted from copy B in default --vg (1 read falls below assembly
    threshold — readthr/isofrac cut).
  - 70 secondary reads carry copy B's splice coordinates into the B region, which
    Layer-2 can use to reconstruct the full 3-exon isoform.

Outputs (written to bench/fixtures/)
--------------------------------------
  sim_starved.fa        — synthetic genome FASTA (chrSIM_STARVED, 50 kb)
  sim_starved.reads.fa  — read FASTA (for provenance / re-alignment)
  sim_starved.bam       — minimap2 splice:hq alignment (sorted, indexed)
  sim_starved.bam.bai
  sim_starved.fa.fai    — produced by samtools faidx

Usage
-----
  python3 bench/sim/gen_starved_fixture.py --out-dir bench/fixtures
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
CHROM         = "chrSIM_STARVED"

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

# Read mix producing genuine starvation
A_FULL_READS  = 40   # copy A full 3-exon isoform — primary at A
A_SKIP_READS  = 10   # copy A exon-skip (exon1+exon3) — primary at A
B_PRIMARY_READS = 1  # copy B full 3-exon with exact B sequence — primary at B
                     #   (PSVs push this read to prefer B over A)
B_XMAP_READS  = 20   # copy-B-expression reads using copy-A sequence — primary at A
                     #   (secondary at B: carries B's splice coords to the side-index)


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

    genome_fa = os.path.join(out, "sim_starved.fa")
    reads_fa  = os.path.join(out, "sim_starved.reads.fa")
    bam_path  = os.path.join(out, "sim_starved.bam")

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

    # Copy B primary: 1 read with exact B sequence → primary at B
    for i in range(B_PRIMARY_READS):
        records.append((f"copyB_primary_{i:03d}", seq_b_full))

    # Copy B cross-map: 20 reads using A sequence → primary at A, secondary at B
    # These represent copy-B expression that cross-maps because copies are near-identical
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
    # 5.  Verify starvation conditions
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

    # Starvation criteria check
    prim_b_int = int(prim_b) if prim_b.isdigit() else -1
    sec_b_int  = int(sec_b) if sec_b.isdigit() else -1

    if prim_b_int >= 1 and sec_b_int >= 10:
        print("\nSTARVATION STRUCTURE CONFIRMED:", file=sys.stderr)
        print(f"  copy B has {prim_b_int} primary read(s) → bundle materialises", file=sys.stderr)
        print(f"  copy B has {sec_b_int} secondary reads → side-index can attach B coords", file=sys.stderr)
        print("  Run rustle --vg to confirm copy B is ABSENT from assembled output.", file=sys.stderr)
    else:
        print("\nWARNING: starvation structure may not be correct.", file=sys.stderr)
        print(f"  Expected: prim_b >= 1, sec_b >= 10; got prim_b={prim_b}, sec_b={sec_b}", file=sys.stderr)

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
