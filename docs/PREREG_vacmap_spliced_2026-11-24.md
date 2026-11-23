# PREREG: VACmap on Iso-Seq reads, against the genome and against isoform sequences (row 1936, 2026-11-24)

## Current state
Registered, not run.

User (2026-11-24): the advisor asked for VACmap as well; "we both know it does not work with spliced alignment but what if
we can force it, maybe an option in vacmap or what if we extracted the loci with all isoforms". Written before any VACmap
alignment of these reads exists. Known at this point (setup only): VACmap (github micahvista/VACmap, cloned 2026-11-24) has
no spliced option (`-mode H|L|S|R|asm`; penalties `-globalpenalty`, `-localpenalty`); it reports a read as a primary record
plus supplementary records, so an intron can only appear as a deletion or as a split.

## Setup
- Sample and truth: HG002 Kinnex on CHM13, the 42,770 selected reads and the repaired truth table of row 1785
  (`r1785_truth_repair/hg/{reads.fa,reads.tsv}`); reused, as in rows 1776-1791, 1844.
- Arm G ("forced"): `vacmap -ref chm13v2.0.fa -read reads.fa -mode L` (HiFi mode), other options default.
- Arm GR: the same with `-mode R` (fixed variation penalty, described as more sensitive to gene conversion).
- Arm T (the user's idea): reference = the spliced sequences of every transcript of `annot_related/chm13.related.gtf`
  that overlaps a candidate locus of the table on the same strand (gffread); `vacmap -mode L` of the reads against them.
- Control TM: minimap2 `-ax map-hifi` (no spliced model) against the same transcript sequences, so that the effect of
  the reference is not credited to the aligner.
- Pick of an arm = the candidate locus of the read's table row that holds the most aligned bases of the read's primary
  and supplementary records on the locus strand (arms T, TM: the transcript's exons in genome coordinates); none =
  "elsewhere", counted wrong when the read is decidable. Baseline = the table's `primary` pick (minimap2 `splice:hq`).

## Step 0: what the alignments look like (reported, no bar)
Per arm: reads aligned; records per read; share of reads whose records together cover >= 90% of the bases the minimap2
primary aligned; arm G / GR: share of the minimap2 introns that reappear as a deletion or a split with both ends within
10 bases; MAPQ distribution on the near ties.

## Step 1: placement (the decision)
Decidable near ties (stratum S2) of the repaired table; baseline primary 553 right / 47 wrong of 600.
- GAIN: arm wrong <= 31 (a third of the room recovered, the bar of rows 1777-1783) AND, on the decidable clear-lead reads
  (S1), arm wrong <= primary wrong + 3.
- NO GAIN: arm wrong >= 47.
- PARTIAL: between.
Reported beside: fixes and breaks against the primary; the 36 decidable TBC1D3 reads of row 1844 (primary 16 / 20).

## Predictions (mine)
G and GR: most reads come back as several records, short exons are lost; placement NO GAIN, more wrong than the primary
(the candidates of a near tie are equally good for any scorer of the reference sequence). T: NO GAIN, close to TM; worse
where the annotation lacks the read's isoform. I do not expect VACmap's MAPQ to separate right from wrong picks.

## Limits
One sample, one individual, the seen quarter of the reads; the truth is the catalogue-free test with its known limits.
The transcript reference depends on the CHM13 annotation (Liftoff / CAT models). A NO GAIN says VACmap as released does
not place these reads better; it says nothing about its use on DNA.
