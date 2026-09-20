#!/usr/bin/env python
"""STARVED-COPY rescue exhibit: does using the multimapping (secondary) reads recover a copy that minimap2
left below the support gate? Plant a 2-copy family where the MINOR copy has only 2 primary reads (below
GATE_MIN_READS=3 -> dropped) but, being near-identical to the dominant sibling, collects many AS-tied SECONDARY
reads. `copy_assign --recover-copies` feeds those secondaries in so the starved copy can clear the gate.

Run: /home/juanfra/miniforge3/bin/python bench/sim_starved.py   (then bench/sim_starved_run.sh)
"""
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(__file__))
from sim_reads import simulate_reads, write_fastq
import sim_genome as sg   # reuse randseq / make_intron / copy_layout_seqs / RNG

OUT = "/home/juanfra/winloci_scratch"
SAM = "/home/juanfra/miniforge3/bin/samtools"


def gene_with_psvs(n_exons, n_psv):
    """one base gene + a set of PSV (exon,offset) sites; returns (exons, introns, psv_sites)."""
    exons = [sg.randseq(sg.RNG.randint(500, 800)) for _ in range(n_exons)]
    introns = [sg.make_intron(sg.RNG.randint(400, 1200)) for _ in range(n_exons - 1)]
    flat = [(ei, off) for ei, e in enumerate(exons) for off in range(len(e))]
    psv_sites = sg.RNG.sample(flat, n_psv)
    return exons, introns, psv_sites


def make_copy(exons, introns, psv_sites, alleles):
    ex = [list(e) for e in exons]
    for (ei, off), a in zip(psv_sites, alleles):
        ex[ei][off] = a
    return sg.copy_layout_seqs(["".join(e) for e in ex], introns)


def main():
    truth, reads = [], []
    chrom = sg.Chrom("sc")
    chrom.add_bg(30000)

    # 3 near-identical copies, ~3.5 kb, 6 PSVs. copies 0,1 healthy (40 reads each, so the family FORMS);
    # copy2 STARVED = 2 primary reads (< gate 3), so without the flag it is dropped and the family looks like
    # 2 copies. --recover-copies feeds copy2's AS-tied secondaries (the 80 reads of its near-identical siblings
    # that minimap2 placed on copy2's locus) into the rescue, so copy2 clears the support gate.
    exons, introns, psv = gene_with_psvs(5, 6)
    a0 = [exons[ei][off] for (ei, off) in psv]                                 # copy0 = reference alleles
    a1 = ["ACGT"[("ACGT".index(b) + 1) % 4] for b in a0]                       # copy1 = +1 per PSV
    a2 = ["ACGT"[("ACGT".index(b) + 2) % 4] for b in a0]                       # copy2 = +2 per PSV
    for ci, (alleles, nr) in enumerate([(a0, 40), (a1, 40), (a2, 1)]):
        g, spliced, intr = make_copy(exons, introns, psv, alleles)
        start = chrom.add_copy(g)
        chrom.add_bg(sg.RNG.randint(2500, 5000))
        role = "STARVED(2 reads)" if nr <= 2 else "healthy"
        truth.append(("starved", ci, "sc", start, start + len(g), nr, role))
        for i, rq in enumerate(simulate_reads(spliced, nr, err=0.003, indel=0.0008, seed=hash(("sv", ci)) & 0xffff)):
            reads.append((f"SC|starved|{ci}|{i}", rq))
    chrom.add_bg(20000)

    with open(f"{OUT}/simsc.fasta", "w") as fh:
        fh.write(f">{chrom.name}\n{chrom.seq()}\n")
    subprocess.run([SAM, "faidx", f"{OUT}/simsc.fasta"], check=True)
    with open(f"{OUT}/simsc_truth.tsv", "w") as fh:
        fh.write("family\tcopy\tchrom\tstart\tend\tn_reads\trole\n")
        for t in truth:
            fh.write("\t".join(str(x) for x in t) + "\n")
    with open(f"{OUT}/simsc_reads.fastq", "w") as fh:
        for name, rq in reads:
            write_fastq(fh, name, rq)
    print(f"genome sc={chrom.pos:,}bp | planted dominant(60 reads) + STARVED(2 reads) near-identical pair (5 PSVs)")
    print(f"reads: {len(reads)} | wrote simsc.fasta, simsc_truth.tsv, simsc_reads.fastq")


if __name__ == "__main__":
    main()
