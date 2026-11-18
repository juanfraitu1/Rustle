# PREREG r1777: assignment of near-tied reads on locus graphs built from the produced GTF and the sample's own alleles (2026-11-17)

> Naming corrected 2026-11-18: "Eichler's rule" here always meant our adaptation to one reference of the Eichler lab's 10-point filter (Guitart 2024); the paper's own procedure is different and was not tested. Definition: `docs/copy_assignment_definition.md`. Rules, bars and predictions are unchanged.

Written before any count of the arm. Workspace `r1777_allele_graph`, script `bench/o2_truth/allele_graph.py`.
User go-ahead 2026-11-17 after r1776 (room 71 of 710) and after being told the lever is r1374's.

## What is tested
User and advisor: graphs from the produced GTF; reads with near-tied secondary alignments are assigned on them.
The rows r1643 / r1649 / r1696 built such graphs with one sequence per copy and equalled minimap2's primary. Here each locus
graph also carries the alleles the sample's own reads show at that locus. Same lever as r1374 (sample's own variant calls,
HG002 showcase PASS, held-out FAIL 6.3% wrong); new here: catalogue-free truth on both HG002 haplotypes (r1776), all family
genes, bar on the measured room.

## Sample and reads
HG002 Kinnex on CHM13, the reads and truth table of r1776 (`r1776_hg002_room/hg/reads.tsv`, `mm.sel.bam`), unchanged.
Scored reads: all 7,770 near-tie reads (S2) and, as a guard, the first 5,000 clear-lead reads (S1) of the table (hash order).

## The graph
- Regions = the placements of the scored reads (the table's loci with AS >= 0.8 x best). All records of
  `HG002_kinnex.t2t.bam` on those regions are extracted; the assembler (`rustle_pipeline.sh assemble`, current defaults) is run on them.
- Nodes = merged exons of the produced GTF. A locus without a produced model has no nodes and so no alleles.
- Anchor reads of a locus = reads with a PRIMARY alignment of mapping quality 60 there (from the full-depth BAM).
- Allele = a substitution inside a node shown by >= 3 anchor reads and >= 0.2 of the anchor reads covering the position.
  Substitutions only (indels in homopolymers are not trusted).
- A read's score at a placement on the graph: AS' = AS + 5 x (its mismatches there that equal an allele of the graph);
  5 = match 1 + mismatch 4 of `splice:hq`.

## Arms
- primary: minimap2's primary (the figure to beat: 639 right / 71 wrong on the 710 decidable S2 reads).
- graph: the adapted Eichler filter on the graph scores: the placement with the highest AS' when it leads every other by more than 10;
  otherwise abstain.
- graph-or-primary: graph's call, else the primary.
- graph5 (reported, not judged): the same with a lead of >= 5.

## Readings, fixed now (decidable S2 reads)
- PASS: graph-or-primary wrong <= 47 (a third of the room recovered, net) AND the graph's own calls wrong <= 5% AND, on the
  decidable guard S1 reads, graph-or-primary wrong <= the primary's wrong + 3.
- FAIL: graph-or-primary wrong >= 71 (no net gain), or the graph's own calls wrong > 10% (the primary's level).
- Anything else: PARTIAL.
- r1374's failure mode, reported next to the verdict: the number of 'outside' S2 reads (origin at none of the placements,
  659) the graph calls. If it calls more than 20% of them, the verdict carries the note "calls reads of copies the reference lacks".
- Reported: calls on tie and no-hit reads; reads moved off the primary and their right / wrong; alleles and nodes counted.

## Predictions
PARTIAL. The graph calls fewer than 150 decidable S2 reads, more than 90% of them right, net gain 5-20 reads (wrong 51-66).
Guard holds. Fewer than 20% of outside reads called.

## Limits stated beforehand
- One sample, and HG002 was r1374's development donor (different families, rule and truth).
- The GTF is assembled from reads at the candidate regions only, not from a whole-sample run.
- Anchors can themselves be misplaced reads; nothing here certifies them.
- The S1 guard reads are among their own locus's anchors (one read of the >= 3 needed).
- Decidable reads only; ties (4,248) stay unplaceable by sequence.
