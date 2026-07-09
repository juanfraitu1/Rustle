#!/usr/bin/env python3
"""PERFECT-COVERAGE K=0 EXPERIMENT.

Two questions, one planted genome:

  Q1. With PERFECT coverage (every read full-length, spanning every PSV column), does the
      abstained ("tied") mass collapse to EXACTLY the exonically-identical copies?
      -> i.e. is coverage the only other cause of abstention, and K=0 the only hard wall?

  Q2. For those K=0 copies, is provenance really absent from the BAM — or does a read that
      carries FLANKING sequence (readthrough past the transcript into the divergent flank)
      still identify its copy?

PLANT (one chrom, 4 co-located copies of one gene):
  A, B : exons diverged (d=0.01)  -> exonic PSVs exist       -> assignable
  C, D : exons IDENTICAL (d=0)    -> NO exonic PSV at all    -> K=0 (the hard wall)
         but their INTRONS and 3' FLANKS are diverged        -> provenance exists OFF-exon

READS (perfect coverage: full-length, no truncation):
  <copy>_exon_i   : the spliced transcript only            -> C/D must be TIED (provably)
  <copy>_flank_i  : spliced transcript + that copy's 3' flank (readthrough) -> carries the
                    copy-specific flank bases

Q1 is answered by the assignment status per copy. Q2 is answered two ways:
  (a) does the pipeline assign the flank-bearing C/D reads?
  (b) DIRECT test: align each flank-bearing read's flank portion to flank_C vs flank_D and
      see whether it recovers the true copy (this is the information-content question,
      independent of whether the current pipeline exploits it).

Outputs FASTQ + FASTA + BAM under $OUT, then prints the two answers.
Run FOREGROUND (small; no sweeps).
"""
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sim_genome import RNG, Chrom, make_gene, mutate, mutate_intron, copy_layout_seqs, randseq
from sim_reads import simulate_reads

OUT = "/home/juanfra/winloci_scratch/k0flank"
MM2 = "/home/juanfra/miniforge3/bin/minimap2"
SAM = "/home/juanfra/miniforge3/bin/samtools"
READS = 60          # per copy per class -> perfect coverage, all full-length
FLANK_LEN = 300     # 3' flank appended for the readthrough class
FLANK_DIV = 0.05    # flank divergence between copies (flanks evolve fast)
EXON_DIV = 0.01     # A,B exonic divergence (planted PSVs)
INTRON_DIV = 0.02   # C,D intron divergence (identical exons, divergent introns)
ERR = 0.003         # HiFi error

os.makedirs(OUT, exist_ok=True)


def build():
    exons, introns = make_gene(5)
    flank_base = randseq(FLANK_LEN)
    copies = {}      # name -> dict(genomic, spliced, flank)

    # A, B : diverged exons AND introns -> exonic PSVs exist
    for name in ("A", "B"):
        ex = [mutate(e, EXON_DIV) for e in exons]
        it = [mutate_intron(i, EXON_DIV) for i in introns]
        g, sp, _ = copy_layout_seqs(ex, it)
        copies[name] = dict(genomic=g, spliced=sp, flank=mutate(flank_base, FLANK_DIV))

    # C, D : EXONS IDENTICAL (the K=0 wall) but introns + flanks diverged
    for name in ("C", "D"):
        ex = list(exons)                                   # <-- identical exons, no mutation
        it = [mutate_intron(i, INTRON_DIV) for i in introns]
        g, sp, _ = copy_layout_seqs(ex, it)
        copies[name] = dict(genomic=g, spliced=sp, flank=mutate(flank_base, FLANK_DIV))

    assert copies["C"]["spliced"] == copies["D"]["spliced"], "C/D must be exonically identical"
    assert copies["A"]["spliced"] != copies["B"]["spliced"]

    # lay them out co-located on one chrom, each followed by ITS OWN flank
    chrom = Chrom("k0chr")
    chrom.add_bg(3000)
    starts = {}
    for name in ("A", "B", "C", "D"):
        starts[name] = chrom.add_copy(copies[name]["genomic"] + copies[name]["flank"])
        chrom.add_bg(2000)
    chrom.add_bg(3000)

    with open(f"{OUT}/k0.fasta", "w") as f:
        f.write(f">{chrom.name}\n{chrom.seq()}\n")
    subprocess.run([SAM, "faidx", f"{OUT}/k0.fasta"], check=True)
    return copies, starts, chrom


def make_reads(copies):
    """perfect coverage: full-length reads, two classes per copy."""
    recs = []
    for name, c in copies.items():
        # exon-only: the spliced transcript (spans every exonic column)
        for i, (r, q) in enumerate(simulate_reads(c["spliced"], READS, err=ERR, seed=hash(name) % 9999)):
            recs.append((f"{name}_exon_{i}", r, q))
        # flank-bearing: transcript + this copy's own 3' flank (readthrough)
        seq = c["spliced"] + c["flank"]
        for i, (r, q) in enumerate(simulate_reads(seq, READS, err=ERR, seed=(hash(name) + 7) % 9999)):
            recs.append((f"{name}_flank_{i}", r, q))
    with open(f"{OUT}/k0.fq", "w") as f:
        for nm, r, q in recs:
            f.write(f"@{nm}\n{r}\n+\n{q}\n")
    return len(recs)


def align():
    with open(f"{OUT}/k0.sam", "w") as sam:
        subprocess.run([MM2, "-ax", "splice", "-uf", "-C5", "--secondary=no",
                        f"{OUT}/k0.fasta", f"{OUT}/k0.fq"], stdout=sam,
                       stderr=subprocess.DEVNULL, check=True)
    subprocess.run(f"{SAM} sort -o {OUT}/k0.bam {OUT}/k0.sam && {SAM} index {OUT}/k0.bam",
                   shell=True, check=True)


def flank_discriminability(copies):
    """Q2(b): does the FLANK portion alone recover the true copy? (information content, pipeline-independent)"""
    def dist(a, b):
        # simple Levenshtein (short strings; sim only)
        prev = list(range(len(b) + 1))
        for i, ca in enumerate(a, 1):
            cur = [i]
            for j, cb in enumerate(b, 1):
                cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb)))
            prev = cur
        return prev[-1]

    correct = total = 0
    for name in ("C", "D"):
        seq = copies[name]["spliced"] + copies[name]["flank"]
        for r, _ in simulate_reads(seq, 20, err=ERR, seed=(hash(name) + 7) % 9999):
            tail = r[-FLANK_LEN:]                       # the flank portion of the read
            dC = dist(tail, copies["C"]["flank"])
            dD = dist(tail, copies["D"]["flank"])
            call = "C" if dC < dD else ("D" if dD < dC else "?")
            correct += (call == name)
            total += 1
    return correct, total


def main():
    copies, starts, chrom = build()
    n = make_reads(copies)
    align()
    lo = min(starts.values()) - 500
    hi = max(starts.values()) + len(copies["D"]["genomic"]) + FLANK_LEN + 500
    region = f"{chrom.name}:{max(0,lo)}-{hi}"
    print(f"planted 4 copies (A,B exonic-PSVs | C,D EXONICALLY IDENTICAL, divergent introns+flanks)")
    print(f"  C/D spliced identical: {copies['C']['spliced'] == copies['D']['spliced']}")
    print(f"  reads: {n} (perfect coverage, full-length; classes: exon-only, flank-bearing)")
    print(f"  region: {region}")
    print(f"\nRUN NEXT (foreground):")
    print(f"  copy_assign --bam {OUT}/k0.bam --fasta {OUT}/k0.fasta --region {region} "
          f"--min-copies 2 --skip-poa-diagnostic --out {OUT}/k0")
    c, t = flank_discriminability(copies)
    print(f"\nQ2(b) DIRECT flank information content (pipeline-independent):")
    print(f"  flank portion recovers the true copy: {c}/{t} ({100*c/max(t,1):.1f}%)")
    print(f"  -> if high, the K=0 read's provenance IS in the BAM (in the flank), just unused.")
    with open(f"{OUT}/region.txt", "w") as f:
        f.write(region + "\n")


if __name__ == "__main__":
    main()
