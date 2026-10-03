#!/usr/bin/env python3
"""Generate the fixture of `tests/o3_candidates.rs` (plan docs/superpowers/plans/2026-10-02-o3-candidates.md, task 8): a two-copy family
whose second copy is ABSENT from the reference, so `o3_candidates` must flag exactly one candidate copy.

Layout (contig `chrT`, 60 kb random sequence; the gene is on `+`):

    copy A   chrT:10000-13100   exons 10000-10300 (300 bp), 11300-11500 (200 bp), 12700-13100 (400 bp); GT-AG introns of 1000 and
                                1200 bp; spliced length 900 bp. In the reference.
    copy B   (was at ~40 kb)    the same gene, its exons with exactly 3% substitutions (27 of 900 bases, each to another base). NOT in
                                the reference: genome.fa holds plain random sequence there, so B's reads can only align to A.

Reads: 60 per copy, each the spliced transcript from a 5' start drawn uniformly in 0..150 (the 3' end is fixed) with 0.2% random
substitutions per base (A_00..A_59, B_00..B_59: the names carry the truth). They are aligned with the pipeline's flags
`minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes`, sorted and indexed (samtools).

copies.tsv: one family `MCL0` with one copy (A) in `mcl_families`' `P.fam.copies.tsv` layout (COPIES_HEADER, src/bin/mcl_families.rs);
`locus_start`/`locus_end` = A's span. copies.fa: A's spliced exon sum, header `>MCL0|0|chrT:10000-13100|+|nexon=3`.

Deterministic: one seeded Mersenne Twister, drawn only through `random()` (the call Python keeps stable across versions). Run from
anywhere with minimap2 and samtools on PATH: `python3 make_fixture.py` writes genome.fa(.fai), reads.bam(.bai), copies.tsv, copies.fa and
README next to itself (relative paths in every command, so the BAM header names no local directory).
"""
import os
import random
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SEED = 20261002
CHROM, GENOME_LEN = "chrT", 60_000
EXONS = [(10_000, 10_300), (11_300, 11_500), (12_700, 13_100)]  # copy A, 0-based half-open
COPY_B_DIVERGENCE = 0.03
READS_PER_COPY, MAX_5P_TRIM, READ_SUB_RATE = 60, 150, 0.002
BASES = "ACGT"

rng = random.Random(SEED)


def below(n):
    return int(rng.random() * n)


def rand_seq(n):
    return "".join(BASES[below(4)] for _ in range(n))


def other_base(b):
    return [x for x in BASES if x != b][below(3)]


def substitute(seq, rate):
    """Each base independently replaced by another base with probability `rate`."""
    return "".join(other_base(b) if rng.random() < rate else b for b in seq)


def main():
    genome = list(rand_seq(GENOME_LEN))
    exon_seqs = [rand_seq(e - s) for s, e in EXONS]
    for (s, e), x in zip(EXONS, exon_seqs):
        genome[s:e] = x
    for (_, donor), (acceptor, _) in zip(EXONS, EXONS[1:]):  # canonical GT..AG introns
        genome[donor:donor + 2] = "GT"
        genome[acceptor - 2:acceptor] = "AG"
    genome = "".join(genome)
    tx_a = "".join(exon_seqs)
    assert tx_a == "".join(genome[s:e] for s, e in EXONS) and len(tx_a) == 900

    # copy B: exactly 3% of the transcript's bases substituted (positions drawn without replacement)
    n_sub = round(COPY_B_DIVERGENCE * len(tx_a))
    positions = set()
    while len(positions) < n_sub:
        positions.add(below(len(tx_a)))
    tx_b = "".join(other_base(b) if i in positions else b for i, b in enumerate(tx_a))
    assert sum(a != b for a, b in zip(tx_a, tx_b)) == n_sub == 27

    with open(os.path.join(HERE, "genome.fa"), "w") as f:
        f.write(f">{CHROM}\n")
        for i in range(0, GENOME_LEN, 60):
            f.write(genome[i:i + 60] + "\n")
    with open(os.path.join(HERE, "reads.fq"), "w") as f:
        for tag, tx in (("A", tx_a), ("B", tx_b)):
            for k in range(READS_PER_COPY):
                read = substitute(tx[below(MAX_5P_TRIM + 1):], READ_SUB_RATE)
                f.write(f"@{tag}_{k:02d}\n{read}\n+\n{'I' * len(read)}\n")

    start, end = EXONS[0][0], EXONS[-1][1]
    exons = ",".join(f"{s}-{e}" for s, e in EXONS)
    header = ("family_id\tcopy_idx\ttid\tchrom\tstart\tend\tn_exon\tstrand\tn_reads\texons\tmax_family_identity\tsource\tgene_id\t"
              "core_hull\tsd_depth\tcore_bp\trep_frac\tmember_status\tlocus_start\tlocus_end")
    row = [ "MCL0", "0", f"DN_{CHROM}_{start}_A", CHROM, str(start), str(end), str(len(EXONS)), "+", str(2 * READS_PER_COPY), exons,
            "1.000", "fixture", "geneA", "NA", "1", str(len(tx_a)), "0.000", "kept", str(start), str(end)]
    assert len(row) == len(header.split("\t"))
    with open(os.path.join(HERE, "copies.tsv"), "w") as f:
        f.write(header + "\n" + "\t".join(row) + "\n")
    with open(os.path.join(HERE, "copies.fa"), "w") as f:
        f.write(f">MCL0|0|{CHROM}:{start}-{end}|+|nexon={len(EXONS)}\n{tx_a}\n")

    run = lambda cmd: subprocess.run(cmd, shell=True, check=True, cwd=HERE)
    run("minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes -t 1 genome.fa reads.fq 2> /dev/null"
        " | samtools sort -o reads.bam - && samtools index reads.bam && samtools faidx genome.fa")
    os.remove(os.path.join(HERE, "reads.fq"))

    version = lambda cmd: subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True).stdout.splitlines()[0].strip()
    with open(os.path.join(HERE, "README"), "w") as f:
        f.write(f"Generated by make_fixture.py (seed {SEED}) with Python {sys.version.split()[0]}, "
                f"minimap2 {version('minimap2 --version')}, {version('samtools --version')}.\n")


if __name__ == "__main__":
    main()
