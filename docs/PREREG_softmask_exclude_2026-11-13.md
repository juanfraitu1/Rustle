# PREREG r1690: the reference soft mask as the exclusion rule for satellite / rRNA families (2026-11-13)

User goal (2026-11-13): "look for ways in which the soft masking could be used to avoid rRNA and other low complexity
regions while still being able to call gene families". Written after the mask-coverage table below (BEDs and FASTA only),
before any per-copy or per-family mask value was computed.

Today's route is row 1559: `RUSTLE_EXCLUDE_BED=censat.bed` on the assign stage (families with most copies in a T2T censat
element are not assigned). It needs an annotation only CHM13 has. The soft mask (lowercase bases) ships with every genome.

## Step 0a, the instrument (job 4470776, `r1690_softmask/mask.*.tsv`, `mask_classes.py`): lowercase share of each class
| | CHM13 v2.0 (T2T) | gorilla GGO.fa (NCBI) | HG002 paternal |
|---|---|---|---|
| genome | 0.617 | 0.473 | 0.402 |
| RepeatMasker Satellite | 1.000 | 0.976 | no BED |
| Simple_repeat / Low_complexity | 1.000 / 1.000 | 0.959 / 0.853 | |
| rRNA | 1.000 | 0.205 | |
| tRNA / snRNA | 1.000 / 1.000 | 0.029 / 0.289 | |
| LINE / SINE / LTR | 1.000 each | 0.513 / 0.820 / 0.374 | |
| outside every RepeatMasker interval | 0.159 | 0.298 | |
| censat: hor, hsat1A/1B/2/3, mon, bsat, gsat, dhor | 0.978-1.000 | | |
| censat: generic `censat` class (in the satellite label) | 0.847 | | |
| censat: rDNA arrays / ct (transition) | 0.624 / 0.629 | | |
The masks are not the same instrument: CHM13's contains its RepeatMasker run whole (so satellite recall there is partly
circular with the censat label), the gorilla mask is a different masker that covers satellites and simple repeats and
not rRNA or most interspersed repeats. Consequence stated now: **no soft-mask rule can reach rRNA on the gorilla FASTA**;
satellites are covered in both; a threshold tuned to interspersed-repeat content will not transfer.

## Ways (per copy from `PREFIX.fam.copies.tsv` exons, case read from the FASTA; a family falls under a way when strictly
more than half of its copies do, the row-1559 majority rule, judged on the base copies table)
- W1: exonic lowercase fraction >= 0.9.
- W1b: exonic lowercase fraction >= 0.5.
- W2: lowercase fraction of the copy's span (introns included) >= 0.9.
- W3: fewer than 100 uppercase exonic bp.
Edge level (a family link must rest on aligned uppercase sequence) changes the default families and is a separate row.

## Labels
- A119b: per family, the class holding most copies: satellite (censat hor / hsat* / mon / bsat / gsat / dhor / censat),
  rDNA (censat rDNA or >= 50% of exon bp in RepeatMasker rRNA), ct, none. ct is reported in its own column and is in no bar.
- Gorilla OR6737 and KB3781: RepeatMasker class with >= 50% of the exon bp (Satellite, rRNA, other).
- HG002: the rDNA scaffold loci of row 1673 (descriptive; no satellite label).

## Bars (A119b decides which way is frozen; then one frozen way transfers)
- Reach: >= 0.90 of the satellite-label families fall under the way (rDNA and ct reported, no bar).
- Cost 1: pair recall of the default families (pairscore, default level) with the way's families removed drops by <= 0.005,
  precision does not drop; every removed truth family and member is listed by name with the repeat class of its exons.
- Cost 2: removed families that carry no satellite / rDNA / ct label hold <= 1% of the default run's assigned rows.
- Unexplained: removed families with no label and no RepeatMasker-majority exon class are counted, not barred.
- Freeze the way with the highest reach among those passing both costs on A119b (none passes: STOP, report the table).
- Transfer (frozen way, no change): gorilla OR6737 and KB3781, reach >= 0.90 on Satellite-label families and Cost 1.
- Only a transfer PASS proposes an opt-in switch (soft mask in place of the BED on the same skip-list hatch) plus one
  confirming assign run on A119b (non-satellite assigned reads within 1%, wall reported).

## Outcome on A119b (2026-11-13, run `audit/g2/r76e`: 1,380 families, 4,815 assigned rows; `r1690_softmask/a76.*`, `ps.a76.txt`): no way passes, STOP as registered
Labels: satellite 251 families, rDNA 15, ct 365, none 749. The 628 families the censat BED excludes today are 365 ct,
241 satellite, 15 rDNA, 7 none: by family count more than half of what row 1559 removes is transition-zone sequence
(row 1559 listed ct in the BED and reported it apart: 746 of the 9,776 assigned reads it dropped).
| way | satellite reach (bar 0.90) | rDNA | ct | none-label families removed / their assigned rows (bar <= 48) | pair precision / recall (full 0.912 / 0.782) | protein-coding truth pairs lost |
|---|---|---|---|---|---|---|
| W1 exon lowercase >= 0.9 | 223 of 251 (0.888) | 6 of 15 | 34 | 60 / 125 | 0.912 / 0.782 | 0 |
| W1b exon lowercase >= 0.5 | 243 (0.968) | 12 | 186 | 263 / 2,935 | 0.902 / 0.646 | 297 |
| W2 span lowercase >= 0.9 | 225 (0.896) | 6 | 24 | 43 / 55 | 0.912 / 0.779 | 6 (SPANXN2-5) |
| W3 < 100 uppercase exonic bp | 216 (0.861) | 6 | 22 | 35 / 121 | 0.912 / 0.782 | 0 |
- W1 and W3 remove no protein-coding truth pair. The ignored non-coding pairs fall 49,940 -> 27,595; those pairs are
  unjudged here (see the corrections section for what they are). Both miss the
  reach bar (W1 by 3 families) and the assigned-row bar.
- chr13 13.9-15.6 Mb (the assign hotspot): 761 copies in 157 families, 758 copies hsat1A, median exon lowercase 1.00;
  155 of the 157 families fall under W1.
- What W1 removes without a label: 36 LINE-majority and 13 LTR-majority families (transposon copies clustered as
  families; MCL13 is one LINE family of 150 copies), 5 Unknown, 1 low-complexity, 5 with no repeat class. 120 of the 125
  assigned rows are three LINE families (MCL13 51, MCL216 49, MCL445 20).
- Satellite families W1 misses (28): mostly the generic `censat` class (131 of their 165 copies), median exon lowercase 0.75.
- rDNA: 6 of 15 families; rDNA-family copies have median exon lowercase 0.80 (the arrays are 0.62 masked).
- W1b is not usable: half of ordinary exons are >= 50% lowercase on CHM13 (none-label copies: median 0.43).

## Amendment 1 (2026-11-13, written before any gorilla or HG002 copy table was read): descriptive transfer of W1
No way was frozen, so nothing is proposed for adoption by this row. Because the question is whether the mask works where
no censat BED exists, W1 (the only way with no pair lost and reach within 3 families of the bar) is run unchanged on
OR6737 (`r1662_edgemiss/e2e_or`), KB3781 (`r1662_edgemiss/kb_all`, with its assigned rows) and HG002 (`e2e_hg`):
families under W1 by RepeatMasker label (gorilla), pair precision / recall with them removed, true pairs lost by name.
Descriptive, no bar; whether a soft-mask switch is built is the user's call.

## Descriptive transfer of W1 (2026-11-13, `r1690_softmask/{or,kb,hg}.*`, `ps.{or,kb,hg}.txt`)
| sample | families | under W1 (copies) | what they are | pair precision / recall, full -> W1 removed | protein-coding truth pairs lost | rDNA-label family |
|---|---|---|---|---|---|---|
| OR6737 | 379 | 5 (88) | 2 LINE, 1 simple repeat, 2 no class | 0.359 / 0.744 -> same | 0 | 1 (16 copies), not under W1 |
| KB3781 | 225 | 2 (10) | no RepeatMasker class | 0.219 / 0.450 -> same | 0 | 1 (3 copies), not under W1 |
| HG002 | 537 | 30 (233) | no repeat annotation to label them | 0.703 / 0.355 -> same | 0 | no label |
- The gorilla family tables hold no satellite-majority family at all (0 of 379, 0 of 225): the A119b satellite load does not
  exist there, so reach cannot be measured on gorilla. The rule is harmless there (no pair lost, 5 and 2 families).
- The gorilla rDNA families are not reached, as Step 0a predicted (rRNA is 0.21 lowercase in that FASTA).
- (Gorilla pair precision here is the protein-hit truth with same-human-gene pairs counted false; see r1655.)

## Reading
1. Exon lowercase >= 0.9 on the majority of a family's copies sets aside families of satellite, rRNA and
   transposon-class exons without losing one protein-coding truth pair in any of four samples (non-coding pairs have no truth). It is not the censat BED: it keeps the 365 transition-zone
   families the BED removes (34 fall under it) and adds 49 LINE / LTR families the BED keeps.
2. It covers the assign hotspot (155 of 157 chr13 hsat1A families), which is what the BED was adopted for.
3. rRNA is reached only where the mask is RepeatMasker-made and the exons are the rRNA itself (CHM13: 5 of 6 such
   families); spacer transcripts of the arrays and the gorilla / HG002 rDNA loci are not. Mask provenance differs per genome, so
   any lower threshold (W1b) is unsafe: half of ordinary CHM13 exons are >= 50% lowercase.
4. Span-based masking (W2) is worse than exon-based: it removes SPANXN, a real family inside repeats.
Scope: only the family-level skip list was tested. "Covers the hotspot" is a family count, not an assign wall.
Not run at the time of this verdict: the confirming assign run with W1 in place of the BED (since run, see the switch section: PASS, 601 s vs 346 s); the edge-level way (still not run: it changes the default families and needs its own row);
`pairscore --noncoding` (sequence acceptance cannot tell shared LINE sequence from duplication);
a sequence-based rRNA rule (alignment to an rDNA unit), which is what rRNA would need.

## Corrections after review (2026-11-13, `r1690_softmask/{genes_w1,ae}.txt`, `ps` reruns; all post hoc)
- **What the W1 families are, by annotated gene (A119b, copies mapped with famtruth.map_copy):** satellite-label 223
  families / 1,384 copies, 1,323 on no gene; rDNA-label 6 / 410, 383 on rRNA genes (33,548 gene pairs: this is the fall in
  ignored non-coding pairs); ct-label 34 / 124, 86 on no gene; none-label 60 / 333, 296 on no gene (MCL13: 145 of 150
  copies unannotated, the rest lncRNA; MCL216 and MCL445 wholly unannotated). "Transposon family" is the exon repeat
  class, not a verified reading of each family.
- **Protein-coding genes with a copy inside a W1 family (24 copies, 21 genes), none in a hit pair:** TBC1D3 / B / D / E /
  G / H / I / K, ANKRD30BL, FAM95C, AL391987.1 (satellite-label families); ANKRD20A1, ANKRD20A3P, PGM5, TM2D3 (ct);
  AC007326.4, ANKRD36, ANKRD36B, BLOC1S6, FAM228B, FAM90A1 (none). Their assigned reads would go with the family.
- **Default path (famreps, `r1662_edgemiss/e2e`; rule judged on the base table, reps table filtered):** 0.907 / 0.809 with
  and without the W1 families, 0 protein-coding pairs lost; satellite reach 224 of 250 (0.896). The table above is the
  pre-famreps run r76e, the one with assign rows.
- **rDNA split (CHM13):** copies whose exons are RepeatMasker rRNA: 383 of 603 at exon lowercase >= 0.9 (median 1.00), 5
  of 6 such families under W1 (the sixth, 220 copies, has exons running into spacer). Copies in a censat rDNA array with
  no rRNA exon class (spacer): 27 of 473 (median 0.29), 1 of 9 families.
- **HG002 rDNA label (row 1673's table `r1673_sechap/softmask.out`, not rerun):** the 103 real rDNA-array loci on unplaced
  scaffolds have median exon lowercase 0.38 and none reaches 0.9.
- Three HG002 copies had no exon bases to score (NaN) and count as not under any way.

## Amendment 2 (2026-11-13; user: "ok, lets do that"): the opt-in switch and its confirming assign run. Written before the run.
`RUSTLE_EXCLUDE_SOFTMASK=FRAC` on the assign stage (`tools/softmask_families.py`, stdlib only, reads the FASTA through its
.fai): families with most copies at >= FRAC lowercase exon bases join the skip list, judged on the base copies table like
`RUSTLE_EXCLUDE_BED`; both may be set. Treated like the BED in the driver's guards and parallel defaults. Its list on r76e
at 0.9 equals the bench W1 list (323 families).
**Run (A119b, default path products of `r1662_edgemiss/asg_def`, same binaries and driver for both arms, one job, 2 h limit):**
arm BED (`RUSTLE_EXCLUDE_BED`, today's route) and arm MASK (`RUSTLE_EXCLUDE_SOFTMASK=0.9`, no BED).
- C1: on families skipped by neither arm, MASK's assigned rows are within 1% of BED's.
- C2: MASK assign wall <= 2 x BED's and the job ends inside its limit (peak memory reported).
- Reported, no bar: assigned rows MASK gains in the transition-zone families the BED skips, and loses in the families only
  MASK skips.
PASS both: the switch stays as the documented opt-in for genomes without a censat BED. Otherwise it stays in the tree
with the failure written next to it, or is removed if C2 fails (an exclusion that does not bound the wall has no use).

## Outcome, Amendment 2 (2026-11-13, job 4470790, `r1690_softmask/asg_{bed,mask}`, `compare.txt`): PASS both bars
| | BED arm | MASK arm (`RUSTLE_EXCLUDE_SOFTMASK=0.9`) |
|---|---|---|
| families skipped | 628 | 323 (260 in both) |
| assign wall | 346 s | 601 s (1.74 x; bar 2 x) |
| rows / assigned rows | 207,377 / 5,399 | 296,370 / 5,788 |
| assigned molecules | 3,797 | 4,113 |
| C1: assigned rows on families skipped by neither | 5,251 | 5,227 (-0.46%; bar 1%; 12 families differ) |
- MASK assigns 561 rows in the 368 families only the BED skips (transition zone; MCL192 106, MCL228 62), and gives up 148
  rows in the 63 families only it skips (MCL13 50, MCL216 49, MCL445 20: the LINE-class families).
- `check_run_tables.py` clean on both arms. Peak memory was not captured (no GNU time on the node).
- The switch stays as the opt-in for genomes without a censat BED. Defaults unchanged; rRNA is still not covered.
