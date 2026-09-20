#!/usr/bin/env python3
"""
gen_paralog_fixture.py — Deterministic 2-copy paralog fixture generator.

Design
------
chrSIM (~50 kb) contains two homologous 3-exon genes:
  Copy A ("expressed")  : exons at [10000,10300), [10800,11100), [11700,12000) — 3×300 bp
  Copy B ("starved")    : exons at [30000,30300), [30800,31100), [31700,32000) — 3×300 bp

Both copies are on the '+' strand.  Copy B is derived from copy A by applying
~1.5% base substitutions at deterministic positions (PSVs), making the copies
~98.5% identical.  A handful of PSV positions (listed in output) are sprinkled
across all three exons to give the Layer-2 PSV machinery something to anchor.

Expression
----------
  Copy A: 40 reads covering the full 3-exon isoform (A_full)
           10 reads covering a 2-exon isoform (skip exon 2) (A_skip)
  Copy B:  5 reads covering its full 3-exon isoform (B_full)

Reads are exact spliced transcript sequences (0 errors, HiFi-like).
They are forward-oriented (IsoSeq-refine convention).

Outputs (written to bench/fixtures/)
--------------------------------------
  sim_paralog.fa        — synthetic genome FASTA
  sim_paralog.reads.fa  — read sequences (FASTA, no quality, for minimap2 -a)
  sim_paralog.bam       — minimap2 splice:hq alignment (sorted, indexed)
  sim_paralog.bam.bai
  sim_paralog.fa.fai    — produced by samtools faidx

Usage
-----
  python3 bench/sim/gen_paralog_fixture.py --out-dir bench/fixtures
"""

import argparse
import os
import random
import subprocess
import sys

# ---------------------------------------------------------------------------
# Parameters (all deterministic / fixed)
# ---------------------------------------------------------------------------
SEED          = 42
GENOME_LEN    = 50_000
CHROM         = "chrSIM"

# Gene copy A
A_EXONS = [(10000, 10300), (10800, 11100), (11700, 12000)]  # [start, end) 0-based

# Gene copy B  — same intron lengths, just shifted by 20 kb
B_OFFSET = 20_000
B_EXONS = [(s + B_OFFSET, e + B_OFFSET) for s, e in A_EXONS]

# Divergence: ~1.5% of exonic bases mutated in copy B
DIVERGENCE_RATE = 0.015

# Expression
A_FULL_READS = 40
A_SKIP_READS = 10  # exon-skip: skip middle exon of copy A
B_FULL_READS = 5

BASES = "ACGT"


def make_genome(seed: int) -> str:
    """Return a genome string of GENOME_LEN random bases with fixed seed."""
    rng = random.Random(seed)
    return "".join(rng.choice(BASES) for _ in range(GENOME_LEN))


def apply_psv(genome: str, exons, rate: float, seed: int):
    """
    Mutate positions within exon intervals at the given rate.
    Returns (mutated_genome_list, psv_positions_0based).
    """
    rng = random.Random(seed)
    g = list(genome)
    psvs = []
    for start, end in exons:
        for pos in range(start, end):
            if rng.random() < rate:
                orig = g[pos]
                alts = [b for b in BASES if b != orig]
                g[pos] = rng.choice(alts)
                psvs.append(pos)
    return "".join(g), psvs


def splice_transcript(genome: str, exons) -> str:
    """Concatenate exon sequences from genome (0-based half-open intervals)."""
    return "".join(genome[s:e] for s, e in exons)


def write_fasta(path: str, records: list):
    """Write [(name, seq), ...] to a FASTA file."""
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            # wrap at 60
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i+60] + "\n")


def generate_reads(reads_path: str, genome_a: str, genome_b: str):
    """
    Emit read FASTA records.
    Returns a list of (read_name, origin_copy, isoform) for provenance.
    """
    records = []
    provenance = []

    # Copy A full 3-exon isoform
    seq_a_full = splice_transcript(genome_a, A_EXONS)
    for i in range(A_FULL_READS):
        name = f"copyA_full_{i:03d}"
        records.append((name, seq_a_full))
        provenance.append((name, "A", "full"))

    # Copy A 2-exon skip isoform (exon 1 + exon 3, skip exon 2)
    seq_a_skip = splice_transcript(genome_a, [A_EXONS[0], A_EXONS[2]])
    for i in range(A_SKIP_READS):
        name = f"copyA_skip_{i:03d}"
        records.append((name, seq_a_skip))
        provenance.append((name, "A", "skip"))

    # Copy B full 3-exon isoform
    seq_b_full = splice_transcript(genome_b, B_EXONS)
    for i in range(B_FULL_READS):
        name = f"copyB_full_{i:03d}"
        records.append((name, seq_b_full))
        provenance.append((name, "B", "full"))

    write_fasta(reads_path, records)
    return provenance, seq_a_full, seq_a_skip, seq_b_full


def run(cmd: str, check=True):
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

    genome_fa   = os.path.join(out, "sim_paralog.fa")
    reads_fa    = os.path.join(out, "sim_paralog.reads.fa")
    bam_path    = os.path.join(out, "sim_paralog.bam")
    sam_tmp     = os.path.join(out, "sim_paralog.tmp.sam")

    # -----------------------------------------------------------------------
    # 1.  Build genome
    # -----------------------------------------------------------------------
    print("Building genome...", file=sys.stderr)
    base_genome = make_genome(SEED)

    # Copy A occupies [10000,12000); copy B is the mutated version at [30000,32000)
    # We embed copy A's sequence directly into the genome, then build genome_b
    # by applying PSVs only to copy B's region.
    # But: the genome is one contiguous sequence.  Both copies share the SAME
    # underlying random genome; their exons are at different positions and thus
    # have different base content (they are both random).  We then MAKE them
    # homologous by copying copy A's exon sequence into copy B's exon positions,
    # then applying divergence to copy B.

    g = list(base_genome)

    # Overwrite copy B's exon regions with copy A's sequence (pre-divergence)
    for (as_, ae), (bs, be) in zip(A_EXONS, B_EXONS):
        a_seq = base_genome[as_:ae]
        for i, ch in enumerate(a_seq):
            g[bs + i] = ch

    genome_pre_psv = "".join(g)

    # Now apply PSVs only inside copy B's exons
    genome_b_str, psvs = apply_psv(genome_pre_psv, B_EXONS,
                                   DIVERGENCE_RATE, seed=SEED + 1)
    # genome_a uses the original (pre-PSV) sequence in its exons
    genome_a_str = genome_pre_psv

    # The final genome for chrSIM: copy A region from genome_a_str,
    # copy B region from genome_b_str, filler elsewhere from base_genome.
    final_g = list(genome_b_str)  # genome_b_str has copy B PSVs; copy A is identical there
    # (genome_b_str only differs from genome_a_str in B_EXONS positions)

    write_fasta(genome_fa, [(CHROM, "".join(final_g))])
    print(f"  Genome written: {genome_fa}  ({GENOME_LEN} bp)", file=sys.stderr)
    print(f"  Copy A exons (0-based): {A_EXONS}", file=sys.stderr)
    print(f"  Copy B exons (0-based): {B_EXONS}", file=sys.stderr)
    print(f"  PSV count in copy B exons: {len(psvs)}", file=sys.stderr)
    print(f"  PSV positions (0-based, first 20): {psvs[:20]}", file=sys.stderr)

    # Compute actual identity
    a_exon_seq = splice_transcript(genome_a_str, A_EXONS)
    b_exon_seq  = splice_transcript("".join(final_g), B_EXONS)
    matches = sum(a == b for a, b in zip(a_exon_seq, b_exon_seq))
    identity = matches / len(a_exon_seq) * 100
    print(f"  Copy A/B exon identity: {identity:.2f}%  ({len(psvs)} PSVs / {len(a_exon_seq)} exonic bp)",
          file=sys.stderr)

    # -----------------------------------------------------------------------
    # 2.  Generate reads
    # -----------------------------------------------------------------------
    print("Generating reads...", file=sys.stderr)
    genome_for_a = genome_a_str   # copy A uses the pre-psv sequence
    genome_for_b = "".join(final_g)  # copy B uses the psv-applied sequence

    provenance, seq_a_full, seq_a_skip, seq_b_full = generate_reads(
        reads_fa, genome_for_a, genome_for_b
    )
    total_reads = A_FULL_READS + A_SKIP_READS + B_FULL_READS
    print(f"  Reads written: {reads_fa}  ({total_reads} reads)", file=sys.stderr)
    print(f"    copyA_full: {A_FULL_READS}  copyA_skip: {A_SKIP_READS}  copyB_full: {B_FULL_READS}",
          file=sys.stderr)
    print(f"  seq_a_full length: {len(seq_a_full)} bp  (3 exons × 300 bp)",
          file=sys.stderr)
    print(f"  seq_b_full length: {len(seq_b_full)} bp", file=sys.stderr)
    a_b_ident = sum(a == b for a, b in zip(seq_a_full, seq_b_full)) / len(seq_a_full) * 100
    print(f"  Copy A_full vs B_full transcript identity: {a_b_ident:.2f}%", file=sys.stderr)

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
    # 5.  Report
    # -----------------------------------------------------------------------
    print("\n=== Fixture verification ===", file=sys.stderr)
    run(f"samtools idxstats {bam_path}")
    sec_count_rc = subprocess.run(
        f"samtools view -c -f 0x100 {bam_path}",
        shell=True, capture_output=True, text=True
    )
    sec_count = sec_count_rc.stdout.strip()
    print(f"Secondary alignments (0x100): {sec_count}", file=sys.stderr)

    splice_check = subprocess.run(
        f"samtools view {bam_path} | awk '$6 ~ /N/' | head -5",
        shell=True, capture_output=True, text=True
    )
    if splice_check.stdout.strip():
        print("Spliced alignments (CIGAR with N) found:", file=sys.stderr)
        for line in splice_check.stdout.strip().split("\n"):
            fields = line.split("\t")
            print(f"  {fields[0]}  FLAG={fields[1]}  POS={fields[3]}  CIGAR={fields[5]}", file=sys.stderr)
    else:
        print("WARNING: no spliced alignments detected!", file=sys.stderr)

    # Primary vs secondary counts per region
    prim_a_rc = subprocess.run(
        f"samtools view -c -F 0x100 {bam_path} '{CHROM}:10000-12000'",
        shell=True, capture_output=True, text=True
    )
    prim_b_rc = subprocess.run(
        f"samtools view -c -F 0x100 {bam_path} '{CHROM}:30000-32000'",
        shell=True, capture_output=True, text=True
    )
    sec_a_rc = subprocess.run(
        f"samtools view -c -f 0x100 {bam_path} '{CHROM}:10000-12000'",
        shell=True, capture_output=True, text=True
    )
    sec_b_rc = subprocess.run(
        f"samtools view -c -f 0x100 {bam_path} '{CHROM}:30000-32000'",
        shell=True, capture_output=True, text=True
    )
    print(f"Primary alignments in copy A region [10000-12000]: {prim_a_rc.stdout.strip()}", file=sys.stderr)
    print(f"Primary alignments in copy B region [30000-32000]: {prim_b_rc.stdout.strip()}", file=sys.stderr)
    print(f"Secondary alignments in copy A region:             {sec_a_rc.stdout.strip()}", file=sys.stderr)
    print(f"Secondary alignments in copy B region:             {sec_b_rc.stdout.strip()}", file=sys.stderr)

    int_sec = int(sec_count) if sec_count.isdigit() else 0
    if int_sec > 0:
        print(f"\nSUCCESS: {int_sec} secondary alignments present — multimapping structure confirmed.", file=sys.stderr)
        print("Layer-2 PSV recovery target: copy B's 3-exon isoform via secondary side-index.", file=sys.stderr)
    else:
        print("\nWARNING: 0 secondary alignments. Layer-2 recovery test will not trigger.", file=sys.stderr)
        print("Consider increasing homology or -N parameter.", file=sys.stderr)

    # Summary for README / logging
    print(f"\nFiles produced:", file=sys.stderr)
    for f in [genome_fa, genome_fa + ".fai", reads_fa, bam_path, bam_path + ".bai"]:
        size = os.path.getsize(f) if os.path.exists(f) else -1
        print(f"  {f}  ({size:,} bytes)", file=sys.stderr)


if __name__ == "__main__":
    main()
