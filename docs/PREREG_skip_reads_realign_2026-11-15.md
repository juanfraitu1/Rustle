# PREREG row 1746 (2026-11-15): do the reads that skip a very short exon carry its bases? (the one graph candidate of row 1744)

Row 1744, post hoc: annotated internal exons of <= 15 bp are skipped by at least as many reads as spell them in 36%
(OR6737), 59% (A119b), 30% (KB3781) of chain x exon units, against under 1% for exons over 60 bp, and the default emits
only the skip chain for a third to half of them. Not shown there: that the skip reads are misalignments and not real
skipping. Written and committed before any read was examined.

## Measure (`bench/famcmp/skip_realign.py`)

Units as row 1744 (annotated chain x internal exon no annotated chain skips), bins <= 15, 16-30, 31-60 bp, at most
1,500 units per bin (contig order shuffled with a fixed seed). Per unit up to 20 full-chain primary reads and 20 skip
reads. For each read, the read bases between 25 bp before the donor and 25 bp after the acceptor (bases inserted at the
junction included) are compared by edit distance with the spliced reference with the exon (A) and without it (B):
`has_exon` when A is better by >= 2 edits, `skips` when B is better by >= 2, else `unclear`.

## Gate

Full-chain reads of the <= 15 and 16-30 bp bins pooled: `has_exon` >= 0.90 on each sample, else INVALID.

## Reading rules (skip reads; OR6737 and A119b decide, KB3781 reported)

- **Misalignment:** `has_exon` share of the skip reads at <= 15 bp exons >= 0.50 on both samples, and the 31-60 bp
  control <= 0.10 on both. Reading: the aligner drops the exon although the read has it; a local realignment of such
  reads against the junctions other reads of the locus show (the graph use) is worth registering, with these units
  as its ceiling.
- **Real skipping:** share <= 0.10 on both samples. Reading: unannotated skipping, no alignment loss; closed, vg not
  pursued for splicing.
- **Otherwise:** mixed; report the shares and the unit counts, no follow-up registered from this row.

## Limits

Annotated exons only, so this is a ceiling measurement, not a de novo rule. A read with the exon partly present or
with errors in a <= 5 bp exon can fall in `unclear` (the 2-edit margin). The cap of 20 reads per unit weights units,
not reads. No default changes.
