# PREREG: frozen assembly rules on two unseen samples, chimp PTR and bonobo PPA (2026-11-15)

User (2026-11-15): "can we use both of these samples to test the assembly? Maybe for determining if there are more
empirical rules that work no matter what sample we use". Written before any assembly of either sample exists.

PTR (`reads/PTR/PTR_mm.bam`, 1,509,097 reads on mPanTro3) and PPA (`reads/PPA/PPA_mm.bam`, 1,575,695 reads on
mPanPan1) are the only samples no rule has been fitted or tested on. They are therefore used as frozen tests only:
nothing is tuned on them here. New rules are mined on the four spent samples and scored on these two once.

## Data
Default `assemble` (driver, binary `~/rustle_target/release` at d7404c03), multi-exon exact intron chains vs RefSeq
(RS_2026_05, `refseq_gff2gtf.py`; the converter reproduces the gorilla truth file model for model). Same library
protocol as gorilla OR6737 (ccs, lima, refine; Eichler lab). Workspace `ptr_ppa_rules`.

## Tests (each sample on its own; a rule "works on every sample" only if both pass)
1. **r1677 stable classes** `cont5 >= full`, `cont3 >= full`, retained intron (`chain_rules.py` Part 1 table, the
   r1677 criterion unchanged): ratio of share correct, class vs rest, below 1 in every read stratum with >= 50
   chains on each side, and stratum-pooled ratio <= 0.67. At least two strata must qualify, else "not measurable".
2. **`trunc5` stage**: precision of the kept set minus the read-count curve at the same number of correct chains
   >= +0.02 (seen: +0.025 / +0.024 / +0.013 / +0.015).
3. **r1677 logistic score**, fitted on OR6737 + A119b pooled, frozen: score curve minus read-count curve at 0.90
   of the default correct chains >= +0.02 (KB3781: +0.006 FAIL). The swap half (B) is not run.
4. **`shiftdrop` indel rule**: dropped chains in RefSeq <= 0.05; their moved intron annotated in RefSeq or CAT
   <= 0.15; correct chains lost <= 0.2%; precision of the filtered GTF not below the equal-size read cut. Sanity,
   not a bar: >= 0.80 of coding RefSeq introns are in CAT. The event-strata bar of row 1739 is not run.

## Limits stated now
- Depth is about 1/10 of A119b: high read strata may not reach 50 per side.
- PPA RefSeq is thinner than PTR (81k vs 109k mRNA): only curve-against-curve differences are compared across
  samples, never absolute precision.
- Tissue: both testis (user, 2026-11-15, stated after submission and before any result): with OR6737 and A119b these are testis samples three and four, weaker evidence than
  the fibroblast test.
- Retained-intron ratio (row 1745) is another session's job and is not repeated here.

## RESULT (register row 1747; job 4472256, workspace `ptr_ppa_rules`: `job.4472256.out`, `rules.ptr_ppa.txt`, `frozen.txt`)
Default assembly, multi-exon chains / correct / precision: PTR 53,970 / 20,982 / 0.389 (64 s), PPA 56,149 / 18,096 /
0.322 (172 s). CAT sanity 0.821 / 0.835 (met).
| test | PTR | PPA | verdict |
|---|---|---|---|
| 1. `cont5 >= full` pooled ratio (4 strata each, all below 1) | 0.10 | 0.24 | PASS both |
| 1. `cont3 >= full` | 0.22 | 0.37 | PASS both |
| 1. retained intron (3 strata each) | 0.32 | 0.37 | PASS both |
| 2. `trunc5` vs read-count curve (bar +0.02) | +0.034 (4,940 dropped, 0.028 correct, keeps 0.993) | +0.013 (2,854 dropped, 0.046 correct, keeps 0.993) | PASS / FAIL |
| 3. frozen logistic score at 0.90 (bar +0.02) | +0.085 (0.507 vs 0.422) | +0.039 (0.401 vs 0.362) | PASS both |
| 4. `shiftdrop` indel: dropped in RefSeq <= 0.05 | 0.079 (353 chains) | 0.029 (485) | FAIL / PASS |
| 4. moved intron annotated <= 0.15 | 0.139 | 0.095 | PASS both |
| 4. correct lost <= 0.2% | 0.13% | 0.08% | PASS both |
| 4. precision vs equal-size read cut | 0.3908 vs 0.3898 | 0.3248 vs 0.3237 | PASS both |
**VERDICT.** The three chain classes hold on six samples now. The frozen score passes on both (it failed KB3781
fibroblast, +0.006): it transfers between testis samples, not shown across tissues. `trunc5` as a filter fails
PPA: over six samples +0.025 / +0.024 / +0.013 / +0.015 / +0.034 / +0.013, the class is always 0.03-0.05 correct
but the gain over a read cut clears the bar on three of six. `shiftdrop` fails one bar on PTR (28 of 353
dropped chains are in RefSeq). Nothing changes in the defaults; both stages stay opt-in. PTR and PPA are now spent
for these four rules.

## Addendum: truncation classes on the two orangutans PPY and PAB (pre-registered 2026-11-15 before the feature tables exist)
User: "lets run the truncation class table". Seen already on these samples (descriptive, no class feature): default
assembly precision 0.192 / 0.228, exact share by read stratum, RefSeq + CAT lift, read length 1,755 / 1,120 bp
against 2,454 / 2,763 for PTR / PPA. Workspace `pongo_asm`, both testis-protocol libraries of the same lab.
- **Bar (test 1 unchanged):** `cont5 >= full`, `cont3 >= full`, retained intron: ratio below 1 in every stratum with
  >= 50 chains per side, pooled ratio <= 0.67, at least two strata, each sample.
- **Fragment explanation of the low precision (descriptive, prediction stated now):** if shorter molecules cause it,
  the share of default chains in `cont5 >= full` or `cont3 >= full` is higher on both orangutans than on PTR and
  PPA in the 21+ read stratum, and the exact share of chains outside both classes in that stratum comes within
  0.05 of PTR / PPA. If neither holds, fragments are not the explanation.
- Reported beside, not a bar: `trunc5` rule gain over the read-count curve. The frozen score and `shiftdrop` are not
  run; PPY and PAB stay unseen for them.

## RESULT addendum (register row 1748; job 4472265, workspace `pongo_asm`: `feat.{ppy,pab}.tsv`, `rules.txt`, `feat.4472265.out`)
| | PPY | PAB | PTR | PPA |
|---|---|---|---|---|
| `cont5 >= full` pooled ratio | 0.09 PASS | 0.08 PASS | 0.10 | 0.24 |
| `cont3 >= full` | 0.20 PASS | 0.22 PASS | 0.22 | 0.37 |
| retained intron | 0.35 PASS | 0.42 PASS | 0.32 | 0.37 |
| chains in either truncation class, 2 / 3-5 / 6-20 / 21+ reads | 0.22 / 0.24 / 0.22 / 0.10 | 0.23 / 0.28 / 0.24 / 0.14 | 0.16 / 0.17 / 0.12 / 0.05 | 0.13 / 0.12 / 0.09 / 0.04 |
| exact share outside both classes, same strata | 0.11 / 0.19 / 0.34 / 0.63 | 0.16 / 0.25 / 0.40 / 0.62 | 0.26 / 0.40 / 0.60 / 0.81 | 0.17 / 0.29 / 0.48 / 0.74 |
| `trunc5` rule vs read-count curve (not a bar) | +0.027 (14,987 dropped, 0.016 correct) | +0.043 (12,317, 0.016) | +0.034 | +0.013 |
**VERDICT.** Classes PASS on both: eight samples. Fragment explanation: first half holds (class share two to three
times PTR / PPA in every stratum), second half does not (outside the classes the 21+ stratum is still 0.63 / 0.62
against 0.81 / 0.74, beyond the 0.05 margin). Truncated chains are more common on the orangutans but explain only
part of the lower precision. PPY and PAB are spent for the three classes and `trunc5`; unseen for the frozen score
and `shiftdrop`.

## Addendum 2: are the unannotated orangutan introns shared between the two species (pre-registered 2026-11-15 before any lifted table)
User: "yes lets do that". Seen already: the gap classes of PPY / PAB at 21+ reads (one new intron sharing a site 0.081
/ 0.070 of chains; no annotated intron 0.036 / 0.064), and the map's self-check (annotated PPY introns landing on an
annotated PAB intron 0.844). Not seen: any new intron lifted.
Map: `genomes/xmap/PPY_to_PAB.paf` + `lift.sh` (direction PPY -> PAB only). Script `bench/famcmp/xspecies_introns.py`.
Unit = intron of a PPY default chain, lifted to PAB. **Shared** = >= 2 primary reads of the PAB library splice exactly
at the lifted position (reads, not PAB's assembly, so PAB's assembler filters do not enter).
Classes, by the PPY chain's reads (6-20, 21+):
- Y, yardstick: introns of chains exact in RefSeq or CAT;
- N1: the new intron of chains with one new intron sharing a site;
- N0: introns of chains with no annotated intron;
- reported beside: the same for chains with 2+ new introns, and the share of N1 / N0 lifted introns that PAB's own
  RefSeq or CAT annotates.
Rates are over introns that lift. **Verdict per class in the 21+ stratum:** real when its shared rate is >= 0.5 x
the Y rate; noise when <= 0.2 x; mixed between. Not informative if the Y rate is below 0.5 or the class has < 100
lifted introns.
Limits stated now: sharing needs the gene expressed in the other animal, so a testis-specific real intron counts and
a real intron of a gene PAB does not express does not; an artefact that recurs at the same sequence in both species
(repeat-mediated) would count as shared; PAB reads are shorter (1,120 bp), which lowers every rate alike.

## RESULT addendum 2 (register row 1750; workspace `pongo_asm/xspecies`)
First run stopped on paftools' partial-lift names (`_t5` / `_t3`); `lift.sh` now drops them (none among these
introns), rerun unchanged otherwise. 85,209 PPY introns in, 83,272 lifted once.
| PPY chain reads | class | introns | lifted | shared in PAB reads | x yardstick | annotated in PAB (RefSeq or CAT) |
|---|---|---|---|---|---|---|
| 21+ | Y annotated chains | 39,036 | 38,360 | 0.989 | 1.00 | 0.995 |
| 21+ | N1 one new intron | 695 | 673 | 0.914 | 0.92 | 0.189 |
| 21+ | N0 no annotated intron | 428 | 402 | 0.898 | 0.91 | 0.303 |
| 21+ | 2+ new introns | 269 | 251 | 0.900 | 0.91 | 0.195 |
| 6-20 | Y | 35,742 | 34,972 | 0.917 | 1.00 | 0.976 |
| 6-20 | N1 | 3,699 | 3,592 | 0.693 | 0.76 | 0.110 |
| 6-20 | N0 | 3,389 | 3,177 | 0.655 | 0.71 | 0.211 |
| 6-20 | 2+ new | 1,951 | 1,845 | 0.718 | 0.78 | 0.145 |
**VERDICT: N1 and N0 are real by the registered bar (0.92 / 0.91 of the yardstick at 21+ reads; bar 0.5).** Nine in
ten unannotated introns of well-supported Bornean chains are spliced at the same position by Sumatran reads, and
only 19-30% of them are annotated in the Sumatran RefSeq or CAT. This part of the orangutan precision gap is
missing annotation, not assembly error. Not shown: that the whole chain is right (unit = intron), and the 5'
shortened class is untouched by this test. A recurrent artefact at the same sequence would also count as shared.

## Addendum 3: sister-species reads as a chain-level truth, and the rules rescored with it (pre-registered 2026-11-15 before any chain table)
User: "ok, lets do all of that in order" (chain-level check both directions and chimp-bonobo; rescore assemblies and
rules; then the 5' shortened class). Seen: addendum 2 (intron level, PPY -> PAB only). Pairs: PTR <-> PPA,
PPY <-> PAB, four maps `genomes/xmap/A_to_B.paf`. Scripts `bench/famcmp/xspecies_chains.py`, `sister_score.py`;
workspace `sister`. Feature tables rebuilt with the real 1-read-floor assemblies and a chain key.
**Sister-confirmed** = every intron of the chain lifts once, to one contig, and >= 1 primary read of the sister
library has exactly the lifted intron chain. Rates are over chains that lift whole; the lift rate is reported per class.
An exact sister read shows a chain is reproducible, not that it is full length: both libraries of a pair can
fragment the same gene the same way. Hence:
1. **Validity gate (STOP for end questions).** Default chains, confirmed rate of the class `cont5 >= full` (full >= 1;
   0.02-0.05 correct on eight samples) against chains exact in RefSeq, standardised to the class's own distribution
   over intron count (1, 2-3, 4-7, 8+) x reads (2, 3-5, 6-20, 21+). If the ratio is >= 0.5 in a sample, sister truth
   is INVALID for end questions there (`trunc5`, 5' class) and is read only for intron-structure questions.
2. **Every rate is reported on two populations:** all chains, and chains that are not a contiguous sub-chain of a
   longer RefSeq or assembled chain. Class comparisons are standardised over intron count x reads as in 1.
3. **Assembly precision:** RefSeq, and RefSeq or sister-confirmed (descriptive).
4. **`trunc5`** (read only if the gate passes in that sample): dropped chains sister-confirmed vs kept; gain over the
   read-count curve with the wide label, bar +0.02.
5. **1-read class** (chains of the 1-read-floor run that the default lacks, one read): wide-correct share against the
   default's 2-read chains, standardised over intron count, population "not a sub-chain". The floor is discarding
   real isoforms when the 1-read share is >= 0.8 x the 2-read share in all four samples; noise when <= 0.5 x.
6. **One-shot tests on PPY and PAB, verdict label RefSeq** (comparable with the six earlier samples), wide label
   beside: frozen r1677 logistic score, bar +0.02 at 0.90 of the default correct chains; `shiftdrop` indel rule,
   the four bars of test 4 above. After this PPY and PAB are spent for every existing rule.
Limits: sister truth needs the isoform expressed and sequenced in the other animal, so it undercounts rare and
long chains (PAB reads average 1,120 bp); a pair shares protocol and tissue, so shared artefacts count as confirmed.

## RESULT addendum 3 (register row 1751; workspace `sister`: `score.txt`, `feat.*.tsv`, `sister.*.tsv`, `floor1/`)
Maps: annotated introns landing on an annotated intron of the sister, PPY->PAB 0.844, PAB->PPY 0.822, PTR->PPA
0.871, PPA->PTR 0.785. Default chains are in the 1-read-floor runs (99.7-99.8%).
| | PTR | PPA | PPY | PAB |
|---|---|---|---|---|
| default chains lifting whole | 0.901 | 0.887 | 0.931 | 0.939 |
| RefSeq-exact chains sister-confirmed | 0.853 | 0.843 | 0.723 | 0.860 |
| RefSeq-absent chains sister-confirmed | 0.592 | 0.583 | 0.408 | 0.614 |
| **1 gate**, `cont5` class vs RefSeq-exact (standardised), all chains | 0.594 / 0.774 = 0.77 | 0.706 / 0.742 = 0.95 | 0.455 / 0.581 = 0.78 | 0.706 / 0.797 = 0.89 |
| gate, not a sub-chain (165 / 153 / 1,285 / 588 chains) | 0.49 | 0.83 | 0.42 | 0.44 |
| 3 precision RefSeq -> RefSeq or sister | 0.389 -> 0.715 | 0.322 -> 0.671 | 0.192 -> 0.499 | 0.228 -> 0.673 |
| 5 1-read class vs 2-read chains, wide label, not a sub-chain | 0.183 / 0.528 = 0.35 | 0.147 / 0.439 = 0.33 | 0.085 / 0.298 = 0.28 | 0.166 / 0.467 = 0.35 |
| 6 frozen score, RefSeq label (bar +0.02) | +0.085 | +0.041 | +0.072 PASS | +0.097 PASS |
| 6 `shiftdrop`: dropped in RefSeq / moved intron annotated / correct lost / precision vs read cut | (row 1747) | (row 1747) | 0.031 / 0.129 / 0.05% / 0.1930 vs 0.1928 PASS | 0.046 / 0.139 / 0.05% / 0.2285 vs 0.2283 PASS |
| `shiftdrop` dropped chains sister-confirmed (all default chains) | | | 0.300 (0.469) | 0.497 (0.671) |
**VERDICT.**
1. **Gate: sister truth is INVALID for end questions on all four samples** (ratio 0.77-0.95). A chain is its introns,
   so a sister read cut anywhere in the same exon reproduces a 5' truncated chain. Test 4 (`trunc5`) and the wide
   column of the frozen score are therefore not read (they come out negative, as a truth that confirms truncations must).
   Row 3's wide precision is an upper bound for the same reason.
2. **1-read class: noise by the registered bar on all four** (0.28-0.35 x the 2-read chains, bar 0.5; the same with
   RefSeq alone, 0.26-0.38). The 2-read floor is not discarding isoforms as good as the ones it keeps.
3. **Frozen score PASS on both orangutans** with the RefSeq label: it passes on PTR, PPA, PPY, PAB (orangutan tissue
   not confirmed), the fibroblast failed.
4. **`shiftdrop` PASS on both orangutans**, all four bars. Its dropped chains are sister-confirmed less often than
   the average chain, but 30-50% of them are: the indel beside the junction recurs in the sister reads, which either
   means part of the class is real or the alignment artefact is sequence-driven. Not resolved here.
PPY and PAB are now spent for every existing rule.

## Addendum 4: the 5' shortened class: alternative starts or cut molecules (pre-registered 2026-11-15 before any end table)
Step 3 of the user's list. Seen: annotated chains missing 5' introns are 0.09 / 0.14 of 21+ read chains on PPY / PAB
against 0.05 / 0.04 on PTR / PPA; chain-level sister confirmation cannot judge ends (addendum 3 gate). Not seen:
any 5' end position of these chains, in either species.
Idea: a start site puts read 5' ends at one position, and the sister species starts at the homologous position;
a cut molecule ends anywhere in the exon, and the sister's cuts fall elsewhere. So the unit is the 5' END POSITION.
Script `bench/famcmp/xspecies_starts.py`; workspace `sister/starts`; four samples, each with its sister.
Classes of default chains (RefSeq, same strand):
- Y: chain exact in RefSeq (yardstick);
- S3: not in RefSeq, a contiguous sub-chain of a RefSeq chain missing introns at the 3' side only (control: its 5'
  end is an annotated start);
- S5: the same missing introns at the 5' side only (the class in question).
Measured per chain: own concentration = share of the reads with exactly the chain whose 5' end is within 10 bp of
the commonest one; sister concentration = share of sister reads with exactly the lifted chain whose 5' end is
within 10 bp of the lifted commonest position. A chain qualifies with >= 6 own and >= 3 sister reads and is
**start-like** when both concentrations are >= 0.5.
**Verdict per sample:** S5 = alternative starts when its start-like share is >= 0.5 x the Y share; cut molecules
when <= 0.2 x; mixed between. Not informative if the Y share is < 0.3, S5 has < 100 qualifying chains, or the S3
control is below 0.5 x Y (then the measure does not recognise annotated starts and nothing is read).
Limits: both libraries of a pair share protocol, so a sequence-driven stop of reverse transcription could recur at
the homologous base and count as a start; no CAGE exists for these species; only chains expressed in both animals
are judged; RefSeq starts themselves are often predicted (XM models).

## RESULT addendum 4 (register row 1752; workspace `sister/starts`: `{ptr,ppa,ppy,pab}.{txt,tsv}`)
Qualifying chains (>= 6 own, >= 3 sister reads) and the share start-like (both concentrations >= 0.5):
| class | PTR | PPA | PPY | PAB |
|---|---|---|---|---|
| Y exact in RefSeq | 9,425: 0.270 | 8,544: 0.265 | 8,299: 0.263 | 7,521: 0.263 |
| S3 control (3' introns missing) | 859: 0.229 (0.85 x) | 647: 0.233 (0.88 x) | 811: 0.202 (0.77 x) | 870: 0.220 (0.83 x) |
| S5 (5' introns missing) | 1,439: 0.163 (0.61 x) | 1,064: 0.207 (0.78 x) | 2,627: 0.138 (0.53 x) | 3,488: 0.120 (0.46 x) |
| Y own concentration >= 0.5 / sister >= 0.5 | 0.423 / 0.344 | 0.413 / 0.342 | 0.401 / 0.355 | 0.440 / 0.328 |
**VERDICT: NOT INFORMATIVE as registered** (Y share 0.26-0.27, below the 0.3 floor on all four samples). Read 5' ends
of annotated chains are within 10 bp of one position in only 40-44% of chains even in their own library, so a
10 bp window does not recognise annotated starts and cannot separate a start from a cut. Descriptive only: S5 is
start-like 0.46-0.78 x as often as annotated chains, the control 0.77-0.88 x, and the two orangutans are the lowest.
Missing for a verdict: a negative control (a class known to be cut molecules) and a start measure that tolerates
broad start regions. The 5' shortened class stays unexplained.
