# PREREG: multi-copy mode, the copies with own reads and no transcript (row 1865, 2026-11-21)

User (2026-11-21): "lets do 2 and then 1" (row 1863: 2 = copies with own reads and no Rustle transcript; 1 = tie groups).

## Diagnosis (job 4491345, `r1861_multicopy_mode/why.{a1,or}.txt`; done before this text, not a test)
Of the 14 copies with >= 2 own spliced molecules and no Rustle transcript, 13 have every own read on a DIFFERENT intron
chain (PLA2G4B: 22 primary reads at MAPQ 60, 22 chains, 30 of 43 introns in >= 2 reads; SPTBN5 5 reads, 5 chains; TSPY2 12
reads, 10 chains), and Rustle emits nothing at the locus at any stage, with either seed pool: its floor is 2 reads on the
SAME chain. StringTie, which assembles from the splice graph, has transcripts at 6 of the 14 (PRY2, TSPY2, XAGE5, PLA2G4B,
SPTBN5, LOC129534646). The 14th (LOC129534646) has Rustle transcripts that the copy mapping gives to an overlapping gene.

## Test: the existing one-read tiers on these loci
Stage `oneread` (rows 1766 / 1806) on row 1858's `good` assemblies of A119b and OR6737. Arms for `toolcompare.py`:
`good` (as it is), `good+strict` (+ `x.oneread.gtf`), `good+ranked20` (+ the best 20% of `x.oneread.ranked.gtf`),
`good+floor1` (the floor-1 assembly `x.oneread.d/f1.gtf`), `stringtie`.
- An arm QUALIFIES if it models >= 7 of the 14 copies AND its junction precision over the named families is not below
  StringTie's on either sample. The qualifying arm with the highest junction precision is the mode's candidate.
- None qualifies: the one-read tiers do not solve it; the next step is a graph-path rule, said, not built here.

## Predictions (mine)
`good+floor1` models most of the 14 at junction precision below StringTie; `good+strict` adds few (its ORF >= 900 nt test
and "no longer chain contains it" test reject fragments); `good+ranked20` in between and the likeliest to qualify.

## Limits
A119b and OR6737 are development samples for the mode; fourteen copies; junction precision is against annotated exon
boundaries. A candidate is judged afterwards on testis samples not used here.

## Result (2026-11-21): no registered arm qualifies; a post hoc arm does
Jobs 4491349 / 4491350 (`r1865_multicopy_oneread/cmp.{a1,or}.tsv`, `why.{a1,or}.txt`). Copies / junction precision /
spliced transcripts over the named families:

| arm | A119b (127 expressed) | OR6737 (36) | of the 14 copies |
|---|---|---|---|
| `good` | 123 / 0.49 / 1,079 | 27 / 0.58 / 132 | 0 |
| `good+strict` | 123 / 0.49 / 1,087 | 28 / 0.60 / 134 | 1 |
| `good+ranked20` | 124 / 0.49 / 1,092 | 29 / 0.59 / 150 | 2 |
| `good+floor1` | 127 / 0.39 / 1,801 | 36 / 0.61 / 303 | 11 |
| `stringtie` | 107 / 0.47 / 958 | 28 / 0.59 / 155 | 6 |

`good+strict` and `good+ranked20` miss the 7-copy bar; `good+floor1` misses junction precision on A119b (0.39 against
0.47). NONE QUALIFIES. Predictions: floor 1 as predicted; `ranked20` likeliest to qualify was wrong.

POST HOC (job 4491364, `rescue.py`, `cmp2.{a1,or}.tsv`; made after seeing the table, not a registered arm): `good+rescue` =
`good` plus the spliced floor-1 transcripts at places where `good` has no transcript on that strand. A119b 127 / 0.49 /
1,093, OR6737 35 / 0.66 / 174, 11 of the 14 copies (PRY, PRY2, TSPY2, XAGE5, PLA2G4B, six SPTBN5). It meets the bar's
numbers on the samples it was designed on. Genome-wide it adds 17,997 (A119b) and 11,798 (OR6737) one-read transcripts whose
precision is not measured. To be judged on unused testis samples before it is called the mode.

## Item 1, tie groups (2026-11-21; a report format, not a test; job 4491370, `r1863_copy_evidence/{ev2,tiegroups}.py`, `groups.{a1,or}.tsv`)
Copies of a family joined when >= 2 molecules tie exactly (AS) between them; per group the scores allow at least
max(1, copies with own reads) and at most all its copies. On the `good+rescue` assembly:
A119b 123 copies with reads: 85 with own reads, 38 tied only; expressed copies between 88 and 123. TSPY 6 to 32 (one
group of 29 holding 3 copies with own reads), CT47A 2 to 7, XAGE 2 to 4, HSFY 1 to 2, BPY2 2 to 3; NPIP, GOLGA8, GOLGA6L,
TBC1D3, SPATA31, RBMY, DAZ, PRY exact. OR6737: 36 copies, all with own reads. Borrowed copies (reads score better
elsewhere): 9 on A119b, 13 on OR6737 (8 MAPKBP1, 4 SPTBN5).

# PREREG: the candidate mode on four unused testis samples (row 1867, 2026-11-21)

User (2026-11-21): "yes lets do that" (judge `good+rescue` on the unused testis samples, with a genome-wide check of the added
one-read transcripts). Written before any of these assemblies exists. Nothing of the candidate is changed.

## Samples and families
Testis of chimp `PTR_mm`, bonobo `PPA_mm`, Sumatran `PAB_mm` and Bornean `PPY_mm` orangutan (workspace `reads/TAG/`), never
used for seeding pools or the mode (PPY / PAB were used for the one-read ranking, rows 1766-1808: said). Families = the 15
human families + the 4 gorilla ones, found in each RefSeq GFF by gene symbol or description
(`r1867_multicopy_judge/targets.py`, written and run on the annotation only: 230 / 262 / 309 / 263 genes).

## Arms
`primary` (`assemble --seed-pool primary`), `good`, `good+rescue` (`rescue.py` of row 1865 unchanged: `good` + the spliced
transcripts of the floor-1 assembly where `good` has none on that strand), and the de novo assemblies of
`apetools/TAG_mm`: `stringtie` (`st_dn.gtf`), `flair`, `isoquant`. `toolcompare.py --targets`.

## Rule
- LOCI PASS: `good+rescue` models more copies than StringTie in the total over the four samples, and at least as many on
  >= 3 of the 4.
- JUNCTIONS PASS: its junction precision over the four samples together (annotated / all, summed) is not below StringTie's.
- The mode is CONFIRMED on the named families only if both pass. FLAIR and IsoQuant are reported beside, no bar.
- GENOME-WIDE CHECK: share of multi-exon transcripts whose intron chain is a RefSeq chain, for `good`, for the transcripts
  the rescue adds, and for StringTie. If the added transcripts are below HALF the precision of `good` on >= 3 of 4 samples,
  the rescue must be confined to family loci before it becomes a mode (said as the consequence, not built here).

## Predictions (mine)
Loci pass. Junctions close to StringTie, either side. Added transcripts well below half of `good`'s chain precision.

## Limits
RefSeq copy annotation of the ape T2T assemblies; a copy counts with any spliced transcript; tie-only copies are counted as
copies (the tie-group range is reported for the mode's arm, not for the tools).

## Result of row 1867 (2026-11-21): CONFIRMED on the named families; the rescue must be confined
Jobs 4491373-4491376 (assemblies), 4491379 (comparison rerun for PPA / PAB), 4491381 (pooled counts);
`r1867_multicopy_judge/{cmp,chain,raw}.TAG.tsv`, commit in `TAG/commit.txt`.

Copies modeled / junction precision / spliced transcripts:

| arm | PTR (157 expressed) | PPA (187) | PPY (183) | PAB (197) | four together |
|---|---|---|---|---|---|
| `primary` | 94 / 0.81 / 223 | 86 / 0.79 / 220 | 165 / 0.77 / 1,186 | 135 / 0.79 / 676 | 480 / 0.785 / 2,305 |
| `good` | 129 / 0.81 / 327 | 139 / 0.77 / 475 | 176 / 0.74 / 1,453 | 174 / 0.75 / 1,154 | 618 / 0.760 / 3,409 |
| `good+rescue` | 157 / 0.76 / 439 | 178 / 0.72 / 635 | 183 / 0.73 / 1,496 | 190 / 0.73 / 1,231 | 708 / 0.733 / 3,801 |
| `stringtie` | 121 / 0.73 / 424 | 126 / 0.74 / 378 | 173 / 0.66 / 1,801 | 167 / 0.70 / 910 | 587 / 0.699 / 3,513 |
| `flair` | 89 / 0.70 / 261 | 86 / 0.70 / 270 | 165 / 0.50 / 2,640 | 152 / 0.55 / 1,383 | 492 / 0.561 / 4,554 |
| `isoquant` | 81 / 0.80 / 201 | 95 / 0.79 / 208 | 158 / 0.81 / 553 | 129 / 0.80 / 341 | 463 / 0.801 / 1,303 |

- LOCI PASS: 708 against 587 in total, more on 4 of 4. JUNCTIONS PASS: 0.733 against 0.699 (below StringTie on PPA alone,
  0.72 against 0.74). CONFIRMED on the named families. Largest gains over StringTie: TSPY (PTR 21 / 6, PPA 25 / 3, PAB
  22 / 10), RBMY (PPA 28 / 13, PAB 6 / 1), NPIP (PTR 44 / 35, PPA 48 / 38). Behind StringTie: TBC1D3 on PPA (9 / 10) and
  PAB (24 / 25), GOLGA8 on PPA (10 / 12).
- GENOME-WIDE CHECK: chain precision against RefSeq, `good` / added by the rescue / StringTie: PTR 0.365 / 0.096 / 0.350,
  PPA 0.324 / 0.059 / 0.281, PPY 0.178 / 0.026 / 0.167, PAB 0.203 / 0.065 / 0.234. The added transcripts are below half of
  `good` on 4 of 4: the rescue MUST BE CONFINED to family loci before it is a mode. `good` itself keeps the chain
  precision of `primary` within 0.006 on every sample.
- Predictions: all three held.
- Deviations: 1 (PPA) and 4 (PAB) target genes the loader does not hold were dropped so `toolcompare.py` could run; the
  pooled junction precision was computed with a copy of the script printing raw counts. Tie-group ranges not run here.

# PREREG: the mode as a driver stage, rescue confined to family loci (row 1869, 2026-11-21)

User (2026-11-21): "lets wire it in as a new piece of code". Built before this text: stage `multicopy` in
`tools/rustle_pipeline.sh` + `tools/multicopy_rescue.py` (self-test in `tools/test_selftests.py`). The stage keeps a floor-1
candidate only when the families stage, run on the assembly + candidates, makes its locus a copy of a family of >= 2 loci.
Written before the stage has run on any sample.

## Run
Stage `multicopy` on the `good` assemblies of rows 1865 / 1867 (A119b, OR6737, PTR, PPA, PPY, PAB), with
`RUSTLE_MULTICOPY_F1` = the floor-1 GTF those rows already made (same BAM, same pool); on OR6737 also once without it, to
exercise the whole stage (the two outputs are compared). SD damping BEDs for A119b / OR6737 as in row 1855, none for the apes.

## Rule (the four apes; they were used by row 1867, so this is a check that the confinement keeps that result, not a new test)
- KEPT: `multicopy.gtf` models more named copies than StringTie in total and at least as many on >= 3 of 4, at pooled
  junction precision not below StringTie's (row 1867's rule), AND the added transcripts are <= 25% of the unconfined
  rescue's on every sample.
- Otherwise NOT KEPT, said with which half failed. Reported beside: copies retained of the unconfined gain (708 - 618 = 90),
  chain precision of the added transcripts against RefSeq, A119b and OR6737 tables.

## Predictions (mine)
Added transcripts fall by more than 90%. Most of the gain is retained where the family already has copies in the assembly;
copies of families with no assembled member (SPTBN5-like) may be lost. KEPT, narrowly on the copy total.

## Result of row 1869 (2026-11-21): KEPT
Jobs 4491407-4491413, 4491416; `r1869_multicopy_stage/{cmp,chain}.TAG.tsv`.

| sample | `good` | `multicopy` (stage) | unconfined rescue | StringTie | added transcripts, stage / unconfined |
|---|---|---|---|---|---|
| PTR | 129 / 0.81 | 155 / 0.76 | 157 / 0.76 | 121 / 0.73 | 1,136 / 11,992 (9.5%) |
| PPA | 139 / 0.77 | 177 / 0.72 | 178 / 0.72 | 126 / 0.74 | 1,472 / 12,867 (11.4%) |
| PPY | 176 / 0.74 | 182 / 0.73 | 183 / 0.73 | 173 / 0.66 | 1,779 / 17,964 (9.9%) |
| PAB | 174 / 0.75 | 187 / 0.73 | 190 / 0.73 | 167 / 0.70 | 1,339 / 14,897 (9.0%) |
| four together | 618 / 0.760 | 701 / 0.732 (6,201 / 8,471) | 708 / 0.733 | 587 / 0.699 | |
| A119b | 123 / 0.49 | 126 / 0.49 | 127 / 0.49 | 107 / 0.47 | 1,166 / 17,997 (6.5%) |
| OR6737 | 27 / 0.58 | 34 / 0.66 | 35 / 0.66 | 28 / 0.59 | 1,262 / 11,798 (10.7%) |

Row 1867's rule holds (701 against 587, more on 4 of 4, junction precision 0.732 against 0.699) and the added transcripts
are 9.0-11.4% of the unconfined rescue's on every ape (bar 25%): KEPT. 83 of the 90 copies of the unconfined gain are
retained. Predictions held, except that copies of families with no assembled member were NOT lost (OR6737 34 of 35).
Not improved by the confinement: the chain precision of the added transcripts against RefSeq, 0.068 / 0.035 / 0.031 / 0.041
(unconfined 0.096 / 0.059 / 0.026 / 0.065; `good` 0.365 / 0.324 / 0.178 / 0.203). The confinement cuts their number by ten,
not their error rate; about a third of the rescued loci sit in families with no copy from the base assembly.

Found by the end-to-end run: the stage first inherited stage `oneread`'s pooling-guard bypass, which another row (1868)
made the default the same day; on OR6737 that gives 27,916 candidates and 1,688 rescued loci instead of 11,798 and 445. The
judged runs had it off, so the stage now has its own switch, `RUSTLE_MULTICOPY_GATE` (unset = 0). With it the end-to-end
run on OR6737 is byte-identical to the run that reuses the floor-1 GTF.

# PREREG: counts that do not multiply over identical copies, tie groups in the output (row 1871, 2026-11-21)

User (2026-11-21): implement (1) read weighting for tied copies and (2) identical copies reported as a group; on the
aligner's cap of 50 secondaries the user holds that 50 is plenty since only near ties are used. Built before this text:
stage `multicopy` now runs `tools/eq_classes.py` on its GTF (ties on minimap2's ms, exact band) and
`multicopy_rescue.py annotate` (self-test passes): transcript attributes `reads_placed`, `reads_em`, `reads_max`,
`tie_group`, `tie_group_loci`, `reads_tieadj`; `PREFIX.multicopy.groups.tsv`. Seen before writing: three identical TSPY10
copies of A119b each carry the same five transcripts with reads 2, 2, 8, 72, 10. No output of the new step has been seen.

## Checks (A119b and OR6737, stage rerun on row 1869's inputs)
- a. STRUCTURE: the transcripts and exons of `multicopy.gtf` are those of row 1869's file (attributes added, nothing else).
- b. NO MULTIPLICATION: over the loci of the TSPY copies of A119b, the sum of `reads_em` is within 10% of the number of
  distinct molecules with a spliced record on those copies (`r1863_copy_evidence/mode.a1.tsv.reads`), while the sum of
  `reads_max` exceeds it at least fivefold.
- c. GROUPS: the 14 TSPY copies row 1870 called IDENTICAL (by AS and bases) fall in at most 2 tie groups (ms, exact band).
All three hold = DONE. Otherwise said which failed; nothing is tuned to make them hold.

## The cap of 50 secondaries (a measurement, for the user's point)
A119b: 367,525 of 21,391,966 molecules (1.7%) have 51 records, the cap. One pass over the BAM: of those, how many have
ALL their records within 2% of their best AS (only then can the cap have cut a near tie), and how many of these lie on a
named-family copy. The user's position holds if that is under 5% of the capped molecules and none of the named families
is affected.

## Predictions (mine)
a holds. b holds. c: one or two groups. Cap: the user is right (most capped molecules are repeats whose 51st record is
far below the best).

## Result (2026-11-21, jobs 4491484-4491486, `r1871_multicopy_counts`, `check.py`): a HOLDS, c HOLDS, b FAILS as written
- a. STRUCTURE holds: md5 of transcript + exon coordinates equals row 1869's file on A119b and on OR6737.
- c. GROUPS holds: all 14 IDENTICAL TSPY copies are in one tie group (G4: 39 chrY loci; 29 of the 32 annotated TSPY copies;
  TSPY2 and two TSPY3 copies with own reads stand alone).
- b. NO MULTIPLICATION fails both halves: over the 32 TSPY loci sum `reads_em` 127.6 against 234 molecules (0.55, bar
  within 10%); sum `reads_max` 1,058 (4.5 x, bar 5 x). Cause, found after: my reference count was the wrong population.
  Of the 234 molecules 79 score better elsewhere (not counted here by design), and group G4 holds 7 loci that are not
  annotated TSPY copies and take a share. Inside the group the counts are conserved: 31 placed + 134 tied = 165
  molecules, sum `reads_em` 164.7, sum `reads_max` 1,443 (8.7 x). POST HOC, not the registered check.
- Seen and not predicted: the EM gives `reads_em` 0.0 to the three identical TSPY10 copies (72-read transcript ->
  `reads_tieadj` 0.0): with an abundance prior the tied reads go to the copies that also have placed reads. For a copy
  with no placed read, `reads_tieadj` 0 is one end of the range 0..`reads_max`, not evidence of absence; the group's
  `min_expressed`..`max_expressed` (G4: 8..39) is the statement to quote.
- Totals: A119b 679 tie groups, 2,234 loci, 1,798..2,234 expressed, 12,849 of 220,943 transcripts in a group; OR6737 213
  groups, 540 loci, 374..540, 2,049 of 70,754 transcripts. Cost: the step took 47 min on A119b (stage 7 -> 54 min), 2 min
  on OR6737.
- Cap of 50: of 367,525 capped molecules 8,896 (2.4%) have all 51 records within 2% of the best AS; 299,272 (81%) have
  at most 5 such records; none of the 8,896 lies on a named-family copy. Under the bar: the user's position holds. My
  prediction a, c and cap right; b wrong.
