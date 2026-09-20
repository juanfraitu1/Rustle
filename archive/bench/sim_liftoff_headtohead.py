#!/usr/bin/env python3
"""HEAD-TO-HEAD: exon-sum+minimap2 vs Liftoff -copies at finding gene copies across divergence.

Plant ONE base gene, then N copies at increasing divergence (0..30%), spaced 50 kb apart on one chrom.
Emit the genome FASTA + a GFF annotating the 0%-divergence copy as the reference gene. Then:
  - Liftoff -copies -sc <t> lifts that gene onto the genome -> which copies does it report?
  - exon-sum: minimap2 the 0%-copy's spliced transcript vs the genome -> which loci at identity >= <t>?
Ground truth = all N are copies of the gene. Report, per copy (by divergence), whether each method found it.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sim_genome import make_gene, copy_layout, Chrom

OUT = "/home/juanfra/winloci_scratch"
DIVS = [0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30]
SPACER = 50000

chrom = Chrom("hh1")
chrom.add_bg(20000)
exons, introns = make_gene(5)
copies = []   # (div, genomic_start, spliced_seq, intron_offsets)
for d in DIVS:
    g, spliced, intr = copy_layout(exons, introns, d)
    start = chrom.add_copy(g)
    chrom.add_bg(SPACER)
    copies.append((d, start, len(g), spliced, intr))

with open(f"{OUT}/hh.fasta", "w") as f:
    f.write(f">hh1\n{chrom.seq()}\n")

# GFF: annotate ONLY the d=0.0 copy (the "reference gene" Liftoff will lift + search for copies of)
d0, s0, glen0, spliced0, intr0 = copies[0]
# exon intervals within the copy (0-based) = complement of introns
exon_iv = []
prev = 0
for (a, b) in intr0:
    exon_iv.append((prev, a)); prev = b
exon_iv.append((prev, glen0))
with open(f"{OUT}/hh_ref.gff", "w") as f:
    f.write("##gff-version 3\n")
    ge, gs = s0 + glen0, s0 + 1
    f.write(f"hh1\tsim\tgene\t{gs}\t{ge}\t.\t+\t.\tID=HHGENE\n")
    f.write(f"hh1\tsim\tmRNA\t{gs}\t{ge}\t.\t+\t.\tID=HHMRNA;Parent=HHGENE\n")
    for k, (es, ee) in enumerate(exon_iv):
        f.write(f"hh1\tsim\texon\t{s0+es+1}\t{s0+ee}\t.\t+\t.\tID=HHEX{k};Parent=HHMRNA\n")

# the d=0 copy's spliced transcript (the exon-sum query for our method)
with open(f"{OUT}/hh_query.fa", "w") as f:
    f.write(f">HHMRNA\n{spliced0}\n")

# truth table of copy loci by divergence
with open(f"{OUT}/hh_truth.tsv", "w") as f:
    f.write("div\tstart\tend\n")
    for d, s, gl, _, _ in copies:
        f.write(f"{d}\t{s}\t{s+gl}\n")
print(f"planted {len(copies)} copies of ONE gene at divergences {DIVS} on hh1 ({SPACER//1000}kb apart)")
for d, s, gl, _, _ in copies:
    print(f"  div={d:.2f}  hh1:{s}-{s+gl}")
