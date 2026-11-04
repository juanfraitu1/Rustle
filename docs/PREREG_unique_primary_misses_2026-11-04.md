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

## Step 0 RESULT (~/r1534/step0_classes.py -> step0.out, step0_members.tsv; 17 s)
| class (registered) | members | with a base single-exon model on the exons |
|---|---|---|
| A single_exon_truth | 51 | 0 |
| B1 overlapping_gene | 1 | 0 |
| B2 unspliced_reads | 20 | 0 |
| B3 spliced_reads_no_model | 52 | 0 |
| B4 few_own_reads (own < 2 under the 50% rule; the panel's P counts any exon overlap) | 53 | 0 |
Decision rule: B3 (52) > A (51) > B2 (20) -> the next row is an assembly row on the B3 members. B3 examples: SPANXA2 58 own
reads, 53 sharing one intron chain (MAPQ 13-22, unique primaries); SIGLEC10 32; ZNF679 31 (22 distinct chains); FAM86B1 27.

### After the numbers (labelled)
1. Context of the 177 in the base GTF (~/r1534/check_context.py -> check_context.out): no model of any kind on the span at
   103 (A 15, B2 12, B3 36, B4 39); only an opposite-strand model at 48; inside a same-strand spliced model's span at 24.
   Own reads agree with the CAT strand at 6,799 of 7,185 (SIGLEC10 9/32, FAM86B1 2/27 do not).
2. Region-slice assembly (A119b records in +-200 kb windows, `assemble --contig`, same binary, no other change) EMITS models the
   genome-wide base lacks at 4 of 7 tried loci: CTAGE15 (1,783 own reads; base: no model touching the exon; slice: 2-3-exon
   models with 39 / 15 / 8 reads), PPP1R2B (slice 18-read model; base none), ZNF679 (slice adds + strand models of 7, 9, 10
   reads; base has only - strand ones), GK3P (3 vs 1). Same at TTC30B (none in both), SIGLEC10, IGFL3 (identical). Files
   /anvme/workspace/v138da13-rustle/r1534/{s,six}/x.gtf, regions.txt. The cause is NOT identified here; a genome-wide run
   loses models that a local run emits at high-read loci. That is the lever to diagnose next, before any emission rule.
3. Confounds of item 2, cleared (advisor): settings match (both runs strict junctions, same mono-primed rule, binary of
   05:01 today); the genome-wide AS table (`--as-table a1.molecules.tsv`) changes nothing (s_tab, six_tab identical counts);
   but with SECONDARIES STRIPPED (`samtools view -F 0x900`) every window-only model vanishes: CTAGE15 6 -> 0, SPANXA2 1 -> 0,
   PPP1R2B 1 -> 0, ZNF679 20 -> 4 (= R1), GK3P 3 -> 1 (= R1). Item 2 is therefore the window form of the r1518 secondary-seeded
   blob, not a loss of the member's own primaries. Withdrawn as a lever.
4. The single-exon truth members (class A) by default form NO model from their own primaries: CTAGE15 1,787 unspliced
   primaries, 4% internally primed, strand vote on or off -> 0 models; no single-exon model survives to the output (the MONO PRIMED line counts 5 on the
   whole chr7 window). With the r1489 opt-in `RUSTLE_MONO_MODEL=cluster` (part of the layered PREFIX.final.gtf, which the panel's
   base x.gtf does not carry: R1 ran RUSTLE_ASSEMBLY_LAYERS=off) the same primaries-only windows emit single-exon models:
   CTAGE15 25 models / 350 reads, PPP1R2B 16 / 265, GK3P 8 / 209, TTC30B 8 / 33 (fragmented by the 50 bp start-end window).
   Two things keep them off the panel: the panel judges the unlayered x.gtf, and its hit rule needs a spliced model.
   Files /anvme/workspace/v138da13-rustle/r1534/{s_prim,s_tab,s_nostrand,s_mono,six_prim,six_tab,six_mono}/x.gtf.
5. B4 (53): under the 50% rule IGFL3 has 0 of 315 primaries at its span with half their bases in its exons, ADH1A 1 of 270;
   the panel's P rule (any exon overlap on the strand) counts them as expressed. Up to 53 of the 177 may be neighbours' reads,
   not misses. A scorer row for A or B4 needs a span guard (r1518-style model span <= 2x gene), or the 100 kb blobs hit for free.
   The mono layer was verified in windows only; no genome-wide A119b final.gtf exists to confirm CTAGE15 appears there, and
   the families stage does not consume final.gtf by the SPLIT design; whether the panel should judge final.gtf is the user's call.
