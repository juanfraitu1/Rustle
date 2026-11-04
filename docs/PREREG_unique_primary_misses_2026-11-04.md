# PREREG r1534: what the 177 unique-primary no-model members are (step 0 only; 2026-11-04)

After r1526-r1533 the r1521 panel's 212 expressed no-model members split into 35 multimapper misses (closed) and 177 whose
primaries are uniquely placed (primary_molecules >= 2, mapq0_share < 0.5, base = 0 in ~/r1526/panel.genes.tsv). r1533 called
them "unspliced / single-exon emission" without checking the reads. This row checks, before any arm is designed.

Panel fact fixed before the count: family_panel.py hit = >= 1 SPLICED model (>= 2 exons) on the gene's strand overlapping its
exons. A truth single-exon member can therefore only be hit by a spliced model of a neighbour.

## Step 0 rule (fixed before the numbers; ~/r1534/step0_classes.py)
Per member, primary records (not secondary / supplementary) overlapping the member's exon span in the A119b BAM, deduped by
name. Per read: own = >= 50% of aligned bases inside the member's exons; foreign = < 50% own and >= 50% inside exons of a
different same-strand CAT gene; spliced = >= 1 N op >= 30 bp in the CIGAR. Base models overlapping the member's exons on its
strand are counted by exon number (single-exon vs spliced) from fusedgw/R1/x.gtf.
Member classes, assigned in this order:
- A `single_exon_truth`: the CAT gene has 1 exon (panel hit impossible by definition); sub-count: a base single-exon model
  overlaps / none.
- B1 `overlapping_gene`: foreign reads >= 50% of the member's primaries.
- B2 `unspliced_reads`: own reads >= 2 and spliced share among own reads < 0.5.
- B3 `spliced_reads_no_model`: own reads >= 2 and spliced share >= 0.5 (reads are spliced; the assembler emitted no spliced
  model) ; sub-count: a base single-exon model overlaps / none.
- B4 `few_own_reads`: own reads < 2.
Decision rule: the largest of A, B2, B3 names the next row; A -> a scorer-definition row (count a single-exon model as a hit
for a single-exon gene), registered separately before it is applied; B3 -> an assembly row on those members; B2 -> the
r1518/r1520 trade, no new row (the user's call). B1 and B4 are recorded, not acted on. No arm in this row; no bars.
