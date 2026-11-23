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

## Where IsoCon's lead comes from (unregistered diagnosis, 2026-11-22; seed 1906, E + D + S pooled; user: what does IsoCon's consensus step use that ours lacks)
The 133 planted copies by the transcripts they were planted from: 106 with canonical introns only, 9 single-exon, 18
with a non-canonical splice site or a zero-length "intron" (lift-over models: TSPY2 has GA-AG, AT-AG and a 0-base intron).
Expressed copies with a transcript:
| class | Rustle primaries | Rustle mode | StringTie | FLAIR | IsoQuant | bambu | IsoCon |
|---|---|---|---|---|---|---|---|
| canonical introns (265) | 259 | 264 | 259 | 228 | 263 | 265 | 256 |
| single exon (23) | 17 | 18 | 22 | 23 | 12 | 18 | 23 |
| non-canonical or 0-length intron (48) | 40 | 42 | 45 | 45 | 43 | 42 | 48 |
- IsoCon's lead is entirely in the last two classes. On copies with ordinary introns it finds FEWER than the mode
  (256 against 264) and merges copies one edit apart (TSPY4 into TSPY10).
- The units the mode misses are not near-identical copies: TSPY2 (741 edits from its nearest copy; 60 uniquely placed
  reads, 31 on one chain; Rustle writes no transcript, StringTie two), GOLGA6L2 (4,159 edits; 60 unique reads, 23
  chains; none), HSFY2's single-exon copy (60 reads; none). Rustle drops them for their splice sites or for being
  single-exon, not for any tie.
- So IsoCon's consensus step has no evidence our provenance test lacks for telling copies apart on these reads. What it
  has is no model of introns: it reports sequences, so splice-site and single-exon rules cannot remove a copy. Its
  other property, finding copies the reference lacks, is not exercised here (sample = reference).
- For this benchmark: 18 of 133 planted copies carry lift-over artefacts no real read would show as introns; scores
  should be given for the canonical class beside the total.

# PREREG: which step of the assembly loses the minor isoforms (row 1916, 2026-11-22)

User (2026-11-22): do the three follow-ups in order. Steps 1 and 2 are answered below from data in hand; this row is
step 3, the low isoform recall on gorilla's own families.
Seen (unregistered diagnosis on the gorilla own-family planted set, E): single-exon copies 147 of 1,328, found by FLAIR
138, the mode 35, StringTie 2, IsoQuant and bambu 0 (Rustle's single-exon floor is the 0.82 quantile of the multi-exon
read counts per contig: 18-58 reads on most contigs, 566 and 597 on two, in a BAM where tied reads count at every
copy); multi-exon isoforms: main isoform found by the mode 0.90-0.94, StringTie 0.99, bambu 1.00; second isoform (10 of
a copy's 30 reads) 0.64-0.74 against 0.82-0.86 and 0.89-0.94; the mode's loss is the same (about a third) whatever the
difference from the main isoform, and the floor-1 assembly of the same reads holds 89-100% of them.
Not seen: any arm below.

## Arms (workspace `r1916_isoform_loss`; `assemble --seed-pool good` on the planted gorilla E BAM, one switch each)
default; RUSTLE_MIN_CHAIN_READS=1; RUSTLE_ASSEMBLY_LAYERS=off; RUSTLE_POLISH_JUNCTION_SNAP=off;
RUSTLE_POLISH_ISM_RATIO=0; RUSTLE_POLISH_ISOFORM_FRACTION=0; RUSTLE_POLISH_RETAINED_RATIO=0; RUSTLE_POLISH_TSS=off;
RUSTLE_BRIDGE_REGROUP=off; and `--seed-pool primary`.
Measures: main and second multi-exon isoforms reproduced exactly; isoform precision (score.py).
A switch EXPLAINS the loss when it lifts the second-isoform recall by at least 0.10 over default. Descriptive otherwise;
nothing becomes a default from this row (one planted set whose isoforms are Rustle's own models).

## Prediction (mine)
RUSTLE_POLISH_ISM_RATIO=0 explains it: 30% of the planted reads are 5'-truncated, a second isoform that shares its 3'
part with the main one is read as its fragment and dropped. MIN_CHAIN_READS=1 lifts recall too, by adding one-read
chains, at a large precision cost.

## Result of row 1916 (2026-11-23, job 4496483 + checks on the same data): the loss is a fault of my planted reads
One switch each on the gorilla own-family planted set (E): main / second multi-exon isoform found, isoform precision.
default 0.910 / 0.653, 0.884; RUSTLE_MIN_CHAIN_READS=1 1.000 / 0.980, 0.894; layers off, junction snap off, isoform
fraction 0, retained-intron ratio 0, TSS off, bridge regroup off: identical to default; ISM ratio 0: 0.910 / 0.662 at
precision 0.582; primaries only 0.875 / 0.608. StringTie 0.992 / 0.793, FLAIR 0.936 / 0.804, bambu 0.997 / 0.940.
Only the read floor explains it (prediction wrong: I named the ISM filter). Why a floor of 2 loses isoforms that have
7-20 reads (found after): 67% of the planted primaries reproduce their chain exactly and 1,795 of 1,818 isoforms have
>= 2 such reads, but Rustle's floor counts DISTINCT records, coordinate duplicates once (row 1850), and my full-length
planted reads of one transcript all start and end on the same bases. Isoforms whose exact-chain reads share one
(start, end): 475, found 32% (99% at floor 1); with two or more distinct pairs: 1,343, found 99.6%.
So the multi-exon deficit of Rustle on planted reads is produced by the simulator, not by the assembler: real reads of
a transcript do not share both ends unless they are PCR duplicates. It has depressed every Rustle number of rows
1903-1911 (the other tools do not collapse duplicates); the comparisons between Rustle arms (primaries, mode, labels)
are affected alike, the comparisons against other tools are biased against Rustle.
Steps 1 and 2 of the same follow-up, from data in hand: (1) TSPY2 and GOLGA6L2 are dropped by `--assembly-junctions
strict`; of the 18 planted copies with non-canonical or zero-length annotated introns, real A119b reads carry those
introns in 3 (GOLGA6L2 228 of 413 reads, GOLGA8Q 28 of 65, GOLGA6L22 4 of 116) and in none of the other 15 (TSPY2: 8
real reads, all canonical): lift-over artefacts, the rule is not changed. (2) the single-exon misses are CDY1B (overlaps
the annotated CDY1, the transcript is credited there) and a 285-base HSFY2 fragment; single-exon models fall under the
floor = 0.82 quantile of the multi-exon read counts of the contig, which in a family-only BAM where tied reads count
at every copy is 18-58 reads on most contigs and 537-597 on some.
`plant.py` now takes 0..20 bases off each end of every read by default (`--ragged 0` = the old reads).

# PREREG: the planted comparison with ragged read ends (row 1917, 2026-11-23)

Seen: row 1916 above (gorilla, the old reads). Not seen: any planted set with ragged ends.
Runs (workspace `r1917_ragged_reads`): human named families, fresh seed 1917, E / D / S, all tools, labels; gorilla
OR6737's own families re-planted (E, S), all tools. Same scorers.
- G1 (the artefact is the cause): on gorilla E, Rustle with near-tie secondaries reproduces >= 0.90 of the second
  multi-exon isoforms (old reads 0.653) and >= 0.97 of the main ones.
- G2 (other tools unaffected): StringTie's and bambu's second-isoform recall changes by under 0.03 from the old reads.
- G3 (human, the comparison redone): on seed 1917 the mode's isoform recall is >= every other genome-based tool's in E,
  D and S, and its copy recall in E >= StringTie's.
Whatever falls, rows 1903-1911 are to be read with this row beside them.
Predictions (mine): G1 holds, G2 holds, G3 holds on isoform recall except against bambu (it was 0.977 on gorilla),
copy recall holds.

## Result of row 1917 (2026-11-23, job 4496486, `r1917_ragged_reads`): G1 HOLDS, G2 HOLDS, G3 holds on isoform recall and FAILS on copy recall
Gorilla OR6737's own families, ragged ends (copy recall / isoform recall / isoform precision; second multi-exon isoform found):
| tool | E | S | second isoform, E (old reads) |
|---|---|---|---|
| Rustle, primaries | 0.978 / 0.982 / 0.950 | 0.973 / 0.979 / 0.874 | 0.994 |
| Rustle, secondaries | 0.959 / 0.964 / 0.909 | 0.966 / 0.973 / 0.717 | 0.978 (0.653) |
| Rustle mode | 0.959 / 0.964 / 0.909 | 0.966 / 0.973 / 0.713 | 0.978 |
| mode + EM test | 0.956 / 0.960 / 0.912 | 0.963 / 0.970 / 0.871 | - |
| StringTie | 0.883 / 0.848 / 0.973 | 0.886 / 0.844 / 0.879 | 0.790 (0.793) |
| FLAIR | 0.999 / 0.942 / 0.958 | 1.000 / 0.949 / 0.859 | 0.835 |
| IsoQuant | 0.841 / 0.814 / 0.984 | 0.850 / 0.817 / 0.905 | 0.790 |
| bambu | 0.888 / 0.902 / 0.728 | 0.892 / 0.901 / 0.678 | 0.940 (0.940) |
- G1 holds: second isoform 0.978, main 0.999. G2 holds: StringTie 0.790 (0.793), bambu 0.940 (0.940).
- On this set Rustle now has the highest isoform recall of all tools (0.96-0.98; FLAIR 0.94-0.95) and the second
  highest copy recall after FLAIR; with the old reads it had the lowest isoform recall. S, silent copies called of 616:
  the mode 270, mode + EM 76, Rustle primaries 85, StringTie 86, FLAIR 110, IsoQuant 69, bambu 84.
Human named families, seed 1917, ragged ends (copy recall / isoform recall / isoform precision):
| tool | E | D | S |
|---|---|---|---|
| Rustle, primaries | 0.970 / 0.820 / 0.800 | 0.970 / 0.812 / 0.812 | 0.957 / 0.806 / 0.592 |
| Rustle mode | 0.970 / 0.824 / 0.629 | 0.970 / 0.824 / 0.654 | 0.957 / 0.815 / 0.449 |
| mode + EM test | 0.962 / 0.820 / 0.641 | 0.955 / 0.812 / 0.662 | 0.943 / 0.798 / 0.576 |
| StringTie | 1.000 / 0.669 / 0.773 | 0.992 / 0.661 / 0.790 | 1.000 / 0.661 / 0.582 |
| FLAIR | 0.992 / 0.816 / 0.741 | 0.985 / 0.791 / 0.781 | 1.000 / 0.734 / 0.532 |
| IsoQuant | 0.940 / 0.590 / 0.815 | 0.947 / 0.598 / 0.841 | 0.914 / 0.565 / 0.569 |
| bambu | 0.970 / 0.745 / 0.622 | 0.970 / 0.703 / 0.712 | 0.957 / 0.710 / 0.458 |
- G3: isoform recall of the mode is the highest of the genome-based tools in E, D and S (holds). Copy recall in E: the
  mode 129 of 133, StringTie 133, FLAIR 132 (fails): the four are TSPY2 and GOLGA6L2 (non-canonical annotated introns),
  the HSFY2 fragment and CDY1B (overlaps CDY1), the artefact copies of row 1916's steps 1 and 2.
- S, units: false units the mode 24, mode + EM 4, primaries + EM 3, StringTie 10, FLAIR 10, IsoQuant 8, bambu 8. Labels:
  84 units available against 73 (StringTie, FLAIR); 20 of 24 silent units with a transcript labelled homology_only; no
  expressed unit mislabelled.
Predictions: G1, G2 right; G3 isoform recall right (and above bambu, which I had excepted); copy recall wrong.
Rows 1903-1911 are superseded for every Rustle-against-other-tool number; their findings about false copies, the EM
test and the labels repeat here.

# PREREG: the three planted results run again with ragged read ends (row 1918, 2026-11-23)

User (2026-11-23): "lets document this, then run all three in order". Rows 1908 (depth sweep), 1906 (units, divergence,
Y-ampliconic families, IsoCon) and 1911 (the labels) were made with reads whose ends coincide. Seen: row 1917 (seed
1917 at 60 reads per copy and gorilla). Not seen: anything below. Workspace `r1918_rerun_ragged`. Claims and bars are
those of the original rows, unchanged; fresh seeds.

## 1. Depth sweep (seed 1918; 3, 5, 10, 20 reads per copy; E and S; all tools)
- D1: at 3 and 5 reads per copy, in E, Rustle with primaries only and StringTie each find >= 10 fewer expressed copies
  than the mode. D2: mode + EM >= 10 more than StringTie there, with no more false units in S. D3: the mode's lead over
  StringTie is smaller at 20 reads than at 5.
## 2. Units, divergence, Y-ampliconic families, IsoCon (seed 1919; 60 reads; E, D, S; IsoCon in all three)
- Y1: mode + EM unit precision above primaries-only and above the mode, units found within one of the mode's (E, D, S).
- Y2: on the Y-ampliconic families unit recall of mode + EM >= every other tool's, IsoCon included, and unit precision
  >= StringTie's. Y3: silent copies 1-5 edits from an expressed one called by mode + EM <= by StringTie.
- Beside it, by class of the planted transcripts (canonical / single exon / non-canonical): expressed copies found per
  tool, the table of the unregistered IsoCon diagnosis.
## 3. Labels (seeds 1920-1924, E and S; and the gorilla set of row 1917)
- L1 more units available than every other genome-based tool (S); L2 at most 2 expressed units labelled homology_only
  (S and E); L3 >= 80% of silent units with a transcript labelled homology_only; each on >= 4 of 5 seeds; L4 gorilla:
  <= 2% and >= 80%.

## Predictions (mine)
1: D1 fails for StringTie at 5 reads and may fail at 3 (with distinct ends Rustle's primaries-only arm no longer loses
copies to the duplicate rule, and the lead of the mode over it shrinks); D3 holds. This is the result most likely to
change. 2: Y1 holds, Y3 holds, Y2 fails on recall by the artefact copies. 3: L1-L4 hold.

## Result of row 1918 (2026-11-23, jobs 4496514-4496516, `r1918_rerun_ragged`)
### 1. Depth sweep (seed 1918): D1 holds at 3 reads and FAILS at 5 for StringTie; D2 holds at 3 and FAILS at 5; D3 HOLDS
All copies expressed, expressed copies with a transcript of 133 (isoform recall):
| reads per copy | Rustle primaries | Rustle mode | mode + EM | StringTie | FLAIR | IsoQuant | bambu |
|---|---|---|---|---|---|---|---|
| 3 | 81 (0.34) | 125 (0.67) | 120 (0.65) | 95 (0.31) | 45 (0.18) | 79 (0.31) | 80 (0.32) |
| 5 | 115 (0.60) | 128 (0.74) | 117 (0.67) | 121 (0.52) | 107 (0.43) | 110 (0.48) | 114 (0.51) |
| 10 | 129 (0.73) | 129 (0.80) | 127 (0.79) | 132 (0.64) | 126 (0.64) | 125 (0.57) | 128 (0.69) |
| 20 | 129 (0.80) | 129 (0.82) | 127 (0.81) | 132 (0.66) | 131 (0.78) | 126 (0.60) | 129 (0.70) |
- D1: behind the mode at 3 / 5 reads: Rustle primaries 44 / 13, StringTie 30 / 7 (bar 10). D2: mode + EM against
  StringTie +25 / -4; false units in S 0 against 1 at both. D3: the lead is 30, 7, -3, -3 copies at 3, 5, 10, 20 reads.
- Half silent at 3 reads: expressed copies found the mode 67 of 70 (22 false units), mode + EM 52 (0), StringTie 39 (1).
- Against row 1908 (old reads): Rustle's primaries arm 44 -> 81 copies at 3 reads; the copy lead of near ties is now
  confined to about 3 reads per copy. The isoform-recall lead is at every depth (0.67 against 0.31 at 3 reads, 0.80
  against 0.64 at 10).
### 2. Units, divergence, Y-ampliconic families, IsoCon (seed 1919): Y1 HOLDS, Y3 HOLDS, Y2 FAILS on recall
Half silent (62 expressed units, 46 silent): units found / false units / unit precision / silent copies 1-5 edits away
called: Rustle primaries 60 / 11 / 0.845 / 5; the mode 60 / 24 / 0.714 / 17; mode + EM 60 / 3 / 0.952 / 1; primaries + EM
60 / 3 / 0.952 / 1; StringTie 62 / 13 / 0.827 / 9; FLAIR 62 / 16 / 0.795 / 13; IsoQuant 57 / 9 / 0.864 / 5; bambu 59 / 11 /
0.843 / 5; IsoCon 62 / 12 / 0.838 / 8. E / D units found of 108: mode 104 / 104, mode + EM 103 / 103, StringTie 108 / 106,
FLAIR 108 / 107, IsoCon 107 / 97.
- Y1 holds (0.952 against 0.845 and 0.714; units 103 / 103 / 60 against the mode's 104 / 104 / 60). Y3 holds (1 against
  StringTie 9). Y2 fails on recall: Y-ampliconic units, S, mode + EM 23 of 25 against StringTie, FLAIR, IsoCon 25;
  precision 0.885 against 0.758, 0.735, 0.735.
- Expressed copies found by class of the planted transcripts (E + D + S): canonical introns (267): every Rustle arm 267,
  StringTie 267, bambu 267, FLAIR 264, IsoQuant 262, IsoCon 245; single exon (22): IsoCon 22, FLAIR 21, StringTie 20,
  Rustle 17, IsoQuant 12; non-canonical or 0-length intron (47): StringTie 47, FLAIR 47, IsoCon 46, Rustle 40-42.
  On copies with ordinary introns Rustle misses none; all its misses are in the two artefact classes.
- IsoCon in D finds 115 of 133 copies (97 of 108 units), with 20 joined transcripts: it is the weakest tool when
  expression is uneven. Its candidates are now substrings of the expressed sequences (the reads are ragged), 0-1 equal.
- Isoform recall E / D / S: the mode 0.824 / 0.816 / 0.815, FLAIR 0.820 / 0.795 / 0.782, bambu 0.749 / 0.711 / 0.726,
  IsoCon 0.678 / 0.561 / 0.685, StringTie 0.674 / 0.661 / 0.637, IsoQuant 0.598 / 0.611 / 0.605.
### 3. Labels (seeds 1920-1924 + gorilla): L1 5 of 5, L2 5 of 5, L3 4 of 5, L4 HOLDS. VALIDATED again
Half silent, units available labelled extension against the best other genome-based tool: 92 / 75, 90 / 76, 88 / 75,
91 / 75, 86 / 75. Expressed units labelled homology_only: S 1, 1, 1, 1, 0; E 2, 1, 2, 1, 1. Silent units with a
transcript labelled homology_only: 0.909, 0.806, 0.786, 0.906, 0.815. Gorilla (1,227 units): 860 units available against
FLAIR 725; 1 of 664 expressed units mislabelled; 193 of 196 silent ones labelled homology_only.
With all copies expressed the labelled extension has 104 of 108 units, StringTie and FLAIR 108.
Predictions: 1 as predicted (fails at 5 reads, holds at 3; D3 holds). 2 as predicted. 3 as predicted.
