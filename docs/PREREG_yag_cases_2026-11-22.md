# PREREG: the cases one by one, recently diverged copies, tandem copies, the Y-ampliconic genes, and IsoCon (row 1906, 2026-11-22)

User (2026-11-22): "look carefully at the present cases: copies that have recently diverged so their sequences are still
similar, tandem copies that could be covered in readthroughs that wrongly merge the exons of one another or contain other
mapping artifacts; the Y ampliconic genes are a pretty good example, and IsoCon demonstrated its usability by
disambiguating YAGs".
Written before any of the new measures is computed. Seen: rows 1903-1905 (copy-level scores of seven tools and two
filters on three planted sets, seed 1903, with silent copies counted as false whether or not they are identical to an
expressed copy); row 1881 (joins); nothing per divergence, nothing at the level of distinct sequences, no IsoCon run on
planted reads.

## What is new here
1. UNITS. Copies whose planted isoforms are identical in sequence form one unit: no read can tell them apart, for any
   method. A unit is expressed when any of its copies is. Unit recall = expressed units with a transcript on any member;
   unit precision = called units that are expressed. This removes the calls nobody could avoid.
2. DIVERGENCE. For every silent copy, the edits between its first isoform and the nearest expressed copy of its family
   (minimap2 asm20, NM + unaligned bases): 0, 1-5, 6-20, > 20. Silent copies called, per class and tool: the "recently
   diverged" case measured directly. Same for expressed copies missed.
3. A description of each family: copies, units, pairwise edits, layout (chromosome span, strands, neighbours within 200 kb).
4. IsoCon as an arm: per family the planted reads pooled (as its targeted design does), `IsoCon pipeline`, its
   `final_candidates.fa` aligned to CHM13 (`minimap2 -ax splice:hq`, primaries) and scored like the others; and at
   sequence level: candidates equal to an isoform of an expressed copy / of silent copies only / of none.
5. A second planted set, seed 1906, same settings: every arm of rows 1903-1905 (the EM test included: its validation on
   data it was not looked at on) + IsoCon. Seed 1903 is reported beside it as development data.
Joined transcripts per family and tool are reported from the scorer as before.

## Claims (seed 1906; half-silent scenario unless said)
- Y1 (near ties used with evidence beat both simple strategies): mode + EM test has unit precision >= Rustle with
  primaries only and >= the mode alone, with unit recall within one unit of the mode's, in E, D and S.
- Y2 (against the other tools): on the Y-ampliconic families (BPY2, CDY, DAZ, HSFY, PRY, RBMY, TSPY, VCY) mode + EM test
  has unit recall >= every other tool's (IsoCon included) and unit precision >= StringTie's.
- Y3 (recently diverged copies): among silent copies 1-5 edits from an expressed copy, mode + EM test calls no more than
  StringTie does.
Each is reported as it falls; none tunes anything.

## Predictions (mine)
Most silent copies called by every tool are identical to an expressed copy (class 0), so unit precision is above 0.9 for
all tools and the differences shrink. Y1: precision holds, recall fails in D (the low-count copies of row 1905).
Y2: holds on recall, unsure on precision. Y3: holds. IsoCon: separates copies that differ by several bases, merges or
drops those 1-2 edits apart at 20-60 reads, and writes no transcript for some low-count copies; its exact-sequence rate
is the highest of all.

## Result (2026-11-22, job 4495087, `r1906_yag_cases/{s1906,s1903}.x/*/{score,units,units.yag}.tsv`): Y1 HOLDS, Y2 FAILS on recall, Y3 HOLDS
The families (`s1906/plant/families.tsv`): 133 copies = 108 distinct units. TSPY 35 copies / 20 units (largest 9
identical), CT47A 7 / 3, BPY2 3 / 1, PRY 2 / 1, RBMY 6 / 4, XAGE 6 / 5; GOLGA8, TBC1D3, GOLGA6L all distinct with
nearest copies 2, 4, 5 edits apart. Edits are on the first isoform; units compare every planted isoform.

Seed 1906 (not looked at before), half silent: 62 expressed units, 46 silent; 63 silent copies = 18 identical to an
expressed copy, 14 at 1-5 edits, 5 at 6-20, 26 further.
| tool | unit recall | false units | unit precision | silent copies called: identical / 1-5 / 6-20 / > 20 |
|---|---|---|---|---|
| Rustle, primaries only (ignore) | 58 of 62 | 12 | 0.829 | 18 / 6 / 1 / 3 |
| Rustle mode (take all near ties) | 59 | 30 | 0.663 | 18 / 14 / 4 / 11 |
| mode + EM test | 59 | 4 | 0.937 | 17 / 1 / 0 / 1 |
| primaries + EM test | 58 | 4 | 0.935 | 18 / 1 / 0 / 0 |
| StringTie | 59 | 12 | 0.831 | 17 / 7 / 1 / 3 |
| FLAIR | 58 | 12 | 0.829 | 10 / 6 / 1 / 4 |
| IsoQuant | 58 | 10 | 0.853 | 17 / 6 / 1 / 2 |
| bambu | 59 | 12 | 0.831 | 18 / 6 / 1 / 3 |
| IsoCon | 62 | 9 | 0.873 | 16 / 4 / 0 / 4 |
All expressed (108 units), units found: mode 104, mode + EM 103, StringTie 105, FLAIR 104, IsoQuant 101, bambu 104,
IsoCon 107. One dominant: mode 103, mode + EM 102, StringTie 104, FLAIR 106, IsoCon 107.
- Y1 holds: mode + EM unit precision 0.937 against 0.829 (primaries) and 0.663 (mode), unit recall 59 = 59 (S),
  103 against 104 (E), 102 against 103 (D).
- Y2 fails on recall: Y-ampliconic families, S, units found mode + EM 24 of 26, StringTie 25, IsoCon 26; unit precision
  0.889 against StringTie 0.758, IsoCon 0.812 (false units 3 against 8 and 6). E / D: 36 of 40 against StringTie 38 / 39,
  IsoCon 39.
- Y3 holds: silent copies 1-5 edits from an expressed one, called by mode + EM 1 of 14, StringTie 7, the mode alone 14.
- IsoCon: every expressed sequence of S is among its candidates exactly, in every family (E: all but one in CDY and
  DAZ), and no candidate equals a sequence only silent copies have. It also writes 1.5-3 candidates per expressed
  sequence that are 5'-truncated versions (`inside_expressed`), and its false units come from aligning those to the
  genome, where a truncated candidate fits a silent copy as well.
- Joined transcripts: 1-6 per tool and scenario, none stands out. Transcripts off the planted copies: mode 45-66,
  mode + EM 11-14, the other tools 6-29.
- Seed 1903 (development) shows the same: mode + EM false units 5 (precision 0.921), others 6-9, the mode 26.
Predictions: "most called silent copies are identical" right for the tools that ignore near ties (17-18 of 28), "unit
precision above 0.9 for all" wrong (0.83-0.87); Y1 "recall fails in D" wrong (one unit); Y2 as feared, on recall;
Y3 right; IsoCon "merges copies 1-2 edits apart" wrong (it recovers them all), "exact rate highest" right on recall.
Reading: the copy-level failure of row 1905 (10 expressed copies lost in D) was mostly identical copies of one unit; at
unit level the EM test costs one unit. Its remaining weakness is 2 Y units that StringTie and IsoCon find.
Limits: as rows 1903-1905; IsoCon was given reads pooled by family from the truth (its targeted design), the other tools
all reads together.

# PREREG: replication of the EM test before it becomes anything by default (row 1907, 2026-11-22)

Seen: rows 1903-1906 (seeds 1903 and 1906 on the human named families; gorilla OR6737's own families scored at copy
level without the EM test). Since row 1906 the EM test is an arm of `run_tools.sh` (`mode+em`, `primary+em`) and the
second step of the driver's `o2` under `RUSTLE_O2_EXTEND=1` (-> `PREFIX.o2.gtf`); unit scores are written per scenario.
Not seen: any of the runs below. Workspace `r1907_o2_replication`.

## Runs
- Human named families, five new seeds 1907-1911, scenarios E, D, S, all tools (+ IsoCon in S).
- Gorilla OR6737, its own 379 families (1,328 copies, the planted set of the driver's `o2bench` check), E and S, the EM
  arms and unit scores added to the existing alignments.
- One real run through the driver: `o2` with `RUSTLE_O2_EXTEND=1` on gorilla OR6737 (row 1888's `usual` assembly inputs),
  to see that the wired extension completes and what it drops (no truth there).

## Claims, per human seed
- R1 = Y1: S, unit precision of mode + EM above primaries-only and above the mode; units found within one of the mode's
  in E, D and S.
- R2 = Y3: S, silent copies 1-5 edits from an expressed copy called by mode + EM <= by StringTie.
- R3: S, false units of mode + EM below every other genome-based tool's (StringTie, FLAIR, IsoQuant, bambu).
- R4 (gorilla, S): unit precision of mode + EM >= StringTie's, and units found by mode + EM >= 0.98 x the mode's.
VALIDATED ON PLANTED READS = R1, R2, R3 each hold on at least 4 of the 5 seeds, and R4 holds. Then, on the user's word of
2026-11-22 ("once validated it will just be used to extend the assembly to multi-copy gene families and use the
near-tied multimappers"), the extension with the test becomes the default of the `o2` stage (`RUSTLE_O2_EXTEND=0` to turn
it off). Otherwise it stays opt-in and the failing claim is the next thing to work on.

## Predictions (mine)
R1, R2, R3 hold on all five seeds. R4: precision holds; the recall half is the risk (gorilla copies have 30 reads and
many are single-exon; the EM may starve more of them).

# PREREG: does the toss cost loci when reads are few? A depth sweep (row 1908, 2026-11-22)

User (2026-11-22): with near ties the primary is a coin toss, so seeding on primaries alone should miss loci and also
make false ones; and a read that fits several identical loci may be real at all of them, so an unassignable read should
count for all of them.
Seen: rows 1903-1907 so far, all at 60 reads per copy (20 in the minority copies of D), where every tool reaches nearly
every expressed copy; in D (seed 1906) units found: primaries only 98 of 108, the mode 103, mode + EM 102, StringTie 104.
Not seen: any run below 20 reads per copy.

## Run (workspace `r1908_depth_sweep`)
Human named families, seed 1908, reads per copy 3, 5, 10, 20; scenarios E and S; all tools and both EM arms.
Measured the way the user's second point asks: an expressed copy counts as found when it has a transcript, identical or
not (copy recall); a call is false only when no member of the unit is expressed (false units, unit precision).

## Claims
- D1 (the toss costs loci at low depth): at 3 and at 5 reads per copy, in E, Rustle with primaries only and StringTie
  each find at least 10 fewer expressed copies than the mode.
- D2 (near ties with the test keep them): at 3 and 5 reads per copy, mode + EM finds at least 10 more expressed copies
  than StringTie in E, and in S has no more false units than StringTie.
- D3: the gap of D1 shrinks with depth: smaller at 20 than at 5 reads per copy.

## Predictions (mine)
D1 holds for primaries-only Rustle (it needs two reads on a chain at a locus) and probably for StringTie. D3 holds.
D2 is the open one: the EM test drops a locus under 2 reads, which may cost at 3 reads per copy what the near ties gain.

## Result of row 1908 (2026-11-22, job 4495090, `r1908_depth_sweep/r{3,5,10,20}/{E,S}/{score,units}.tsv`): D1 HOLDS, D3 HOLDS, D2 FAILS by one copy at 5 reads
All copies expressed, expressed copies with a transcript, of 133:
| reads per copy | Rustle, primaries | StringTie | FLAIR | IsoQuant | bambu | Rustle mode | mode + EM test |
|---|---|---|---|---|---|---|---|
| 3 | 44 | 91 | 52 | 76 | 77 | 121 | 115 |
| 5 | 63 | 117 | 106 | 112 | 113 | 128 | 126 |
| 10 | 97 | 124 | 114 | 125 | 126 | 126 | 124 |
| 20 | 117 | 126 | 118 | 127 | 129 | 129 | 127 |
Half silent, expressed copies found of 70 / false units: 3 reads: primaries 13 / 0, StringTie 40 / 2, mode 57 / 27,
mode + EM 41 / 0; 5 reads: 29 / 1, 59 / 2, 65 / 30, 60 / 1; 10 reads: 44 / 2, 64 / 5, 66 / 31, 66 / 2; 20 reads: 61 / 5,
66 / 7, 66 / 30, 66 / 2.
- D1 holds: at 3 and 5 reads per copy primaries-only Rustle finds 77 and 65 fewer expressed copies than the mode,
  StringTie 30 and 11 fewer.
- D2: mode + EM against StringTie in E, +24 at 3 reads (holds), +9 at 5 reads (bar 10: fails by one); false units in S
  0 against 2 and 1 against 2 (holds).
- D3 holds: the mode's lead over StringTie is 30, 11, 2, 3 copies at 3, 5, 10, 20 reads; over primaries-only Rustle 77,
  65, 29, 12.
- With half the copies silent and 3 reads per copy the EM test costs 16 of the mode's 57 expressed copies and leaves it
  level with StringTie (41 against 40); from 5 reads up it keeps the lead or ties, always with the fewest false units.
- The mode alone has 27-31 false units at every depth.
Predictions: D1 right, D3 right, D2 open and it fell one copy short at 5 reads; the EM cost at 3 reads in S as feared.
Reading: the user's argument holds. At low coverage the aligner's toss leaves copies without enough primaries, and the
near-tied secondaries seed them; at 20-60 reads per copy every tool gets there by chance, which is why rows 1903-1907
showed no recall lead. Primaries-only Rustle is the weakest arm at low depth (its floor of two reads on one chain).

## Result of row 1907 (2026-11-22, job 4495089, `r1907_o2_replication`): R2, R3, R4 HOLD; R1 holds on 2 of 5 seeds. NOT VALIDATED by the registered rule; the extension stays opt-in
Half silent, per seed (units found of expressed / false units / unit precision / silent copies 1-5 edits away called):
| seed | Rustle primaries | mode | mode + EM | StringTie | FLAIR | IsoQuant | bambu | IsoCon |
|---|---|---|---|---|---|---|---|---|
| 1907 | 60 of 63 / 10 / 0.857 / 4 | 61 / 26 / 0.701 / 9 | 60 / 2 / 0.968 / 0 | 60 / 11 / 0.845 / 5 | 59 / 14 / 0.808 / 5 | 59 / 9 / 0.868 / 3 | 60 / 10 / 0.857 / 4 | 61 / 7 / 0.897 / 3 |
| 1908 | 57 of 60 / 10 / 0.851 / 4 | 57 / 36 / 0.613 / 13 | 57 / 2 / 0.966 / 0 | 58 / 14 / 0.806 / 5 | 59 / 12 / 0.831 / 5 | 56 / 7 / 0.889 / 2 | 57 / 11 / 0.838 / 4 | 60 / 11 / 0.845 / 5 |
| 1909 | 56 of 60 / 12 / 0.824 / 2 | 57 / 36 / 0.613 / 9 | 55 / 4 / 0.932 / 0 | 59 / 13 / 0.819 / 3 | 58 / 12 / 0.829 / 4 | 55 / 12 / 0.821 / 3 | 57 / 14 / 0.803 / 3 | 59 / 12 / 0.831 / 3 |
| 1910 | 57 of 59 / 10 / 0.851 / 3 | 57 / 26 / 0.687 / 9 | 56 / 2 / 0.966 / 1 | 59 / 13 / 0.819 / 6 | 57 / 12 / 0.826 / 5 | 55 / 11 / 0.833 / 6 | 57 / 12 / 0.826 / 5 | 58 / 7 / 0.892 / 3 |
| 1911 | 60 of 62 / 8 / 0.882 / 3 | 60 / 26 / 0.698 / 11 | 59 / 2 / 0.967 / 0 | 60 / 9 / 0.870 / 5 | 61 / 6 / 0.910 / 2 | 58 / 6 / 0.906 / 3 | 61 / 9 / 0.871 / 4 | 61 / 7 / 0.897 / 3 |
Units found of 108, E / D: mode 104 / 104 (1910: 101 in D), mode + EM 102 / 103, 102 / 103, 103 / 103, 103 / 100, 103 / 103;
StringTie 106 / 100, 103 / 101, 105 / 104, 103 / 104, 105 / 102.
- R1: precision above both simple strategies on all five seeds; "units within one of the mode's" fails on three
  (two units short in E on seeds 1907 and 1908, in S on 1909). Holds on 2 of 5: FAILS the 4-of-5 rule.
- R2 holds on 5 of 5 (0, 0, 0, 1, 0 against StringTie's 5, 5, 3, 6, 5).
- R3 holds on 5 of 5 (false units 2, 2, 4, 2, 2 against the best other genome-based tool's 9, 7, 12, 11, 6).
- R4 (gorilla OR6737's own families, 1,227 units, S): mode + EM 607 of 684 units, 6 false, precision 0.990; the mode 608,
  178 false, 0.774; StringTie 609, 13, 0.979; FLAIR 682, 40, 0.945. Holds (my prediction of a recall problem was wrong).
- The wired extension on the real gorilla run completed: 445 rescued loci; 20,120 loci labelled expressed 18,868, shared
  46, homology_only 1,206.
Reading: the test's benefit is stable (false units a third to a seventh of the other tools', on every seed and on
gorilla). Its cost is one to two units per scenario, which the registered rule did not allow. Since the extension now
labels instead of dropping (user and advisor, 2026-11-22), that cost is one to two expressed units carrying the label
`homology_only`, not missing transcripts. The default is not changed.
Predictions: R1 on all five wrong (2 of 5); R2, R3 right; R4 right on precision, wrong about the recall risk.

# PREREG: validation of what is now the default, the labels (row 1911, 2026-11-22)

Row 1907's rule was written for a filter that DROPS loci, and failed on "units within one of the mode's". On the user's
word the default of `o2` is now the extension that drops nothing and LABELS every locus (`o2_status` expressed / shared /
homology_only). That form has not been scored. Seen: rows 1905-1908 (the dropping form), one labelled planted set by
hand (seed 1906, S: of 47 silent members with a transcript 28 homology_only; of 67 expressed members 1 homology_only).
Not seen: seeds 1912-1916, any label on gorilla.

## Measures (`bench/o2_families/label_score.py`; units = copies with identical planted isoforms)
A unit's label = the best label among the loci of its copies (expressed > shared > homology_only); no transcript = absent.
- AVAILABLE: units, expressed or silent, with a transcript (the advisor's measure), against every other tool.
- Expressed units labelled expressed or shared (right), homology_only (wrong), absent.
- Silent units with a transcript labelled homology_only (right) or expressed / shared (wrong).

## Claims, human named families, seeds 1912-1916, half silent (and E for L2)
- L1 (available): the labelled extension has a transcript at more units (expressed + silent) than every other
  genome-based tool, on every seed.
- L2 (expressed units keep their label): at most 2 expressed units are labelled homology_only, in S and in E, on every
  seed.
- L3 (silent units are told apart): at least 80% of the silent units that have a transcript are labelled homology_only,
  on every seed.
- L4 (gorilla OR6737's own families, S): L2 as a share, at most 2% of expressed units with a transcript labelled
  homology_only; L3 at least 80%.
VALIDATED = L1, L2, L3 each on at least 4 of 5 seeds, and L4. Whatever the outcome the default stays the user's call
(it is on, on their word); a failure is reported with what it would take to fix.

## Predictions (mine)
L1 holds on all seeds. L2 holds on 4 or 5 (the 1-2 units of row 1907). L3 holds (about 85-95%). L4: L3 holds; the 2%
half is the risk on gorilla (30 reads per copy).

## Result of row 1911 (2026-11-22, job 4495102, `r1911_label_validation`): L1 5 of 5, L2 5 of 5, L3 4 of 5, L4 HOLDS. VALIDATED by the registered rule
Half silent, per seed:
| seed | units available, labelled extension (expressed + silent) | best other genome-based tool | expressed units labelled homology_only (S / E) | silent units with a transcript labelled homology_only |
|---|---|---|---|---|
| 1912 | 61 + 28 = 89 | StringTie 69 | 2 / 1 | 22 of 28 (0.786) |
| 1913 | 59 + 25 = 84 | StringTie 72 | 1 / 1 | 20 of 25 (0.800) |
| 1914 | 59 + 27 = 86 | bambu 74 | 0 / 1 | 23 of 27 (0.852) |
| 1915 | 57 + 29 = 86 | StringTie 73 | 0 / 1 | 25 of 29 (0.862) |
| 1916 | 58 + 27 = 85 | StringTie, FLAIR 70 | 0 / 1 | 25 of 27 (0.926) |
- L1 holds on 5 of 5: 84-89 units available against 68-74.
- L2 holds on 5 of 5: 0-2 expressed units carry `homology_only` (of 59-62 in S, 1 of 108 in E).
- L3 holds on 4 of 5 (seed 1912: 0.786, bar 0.80).
- L4, gorilla OR6737's own families: 1 of 608 expressed units with a transcript labelled `homology_only` (0.2%; E: 3 of
  1,105); 172 of 178 silent units with a transcript labelled `homology_only` (0.966). Units available 786 against FLAIR
  722, bambu 624, StringTie 622. Holds.
- With all copies expressed the labelled extension has a transcript at 104 of 108 units on every seed (StringTie 103-105,
  FLAIR 103-106): no lead there, as in rows 1903-1907.
Predictions: L1, L2, L3 right (L3 lower than my 85-95% on two seeds); L4 right, the 2% risk did not show.
So the form that is the default (label, drop nothing) passes the rule written for it before its data were seen.
