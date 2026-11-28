# PREREG: which read end should a transcript take? (row 1992, 2026-11-28)

User goal: improve the transcript starts and ends (TSS / TES); look at what gffcompare and the annotation miss, and at what other tools do.

## What is known (not to be rediscovered)
- gffcompare matches a spliced transcript by its intron chain and never looks at its two ends; an end can be 5 kb off and
  the transcript is still "=". Nothing in the headline numbers sees ends.
- Today's rule: both ends = the most extended end among the reads with exactly the chain; one model per chain.
- Row 1718 (CAGE peaks, human): the most extended 5' end lies in a CAGE peak 0.417 (HG002) / 0.288 (A119b); the median of
  the reads 0.864 / 0.711; the commonest start 0.891 / 0.737. The median layer exists (`--polish-tss-end median`, row 1487).
- It is NOT the default because: against RefSeq starts it hurts every non-human ape library by 15-25 points (row 1834;
  those starts are predicted models at the most upstream read, and there is no CAGE for them), simulated sets -1.2 to
  -1.5, mouse heart -2; and moved ends changed the families on A119b (rows 1490 / 1491).
- 3' ends: a tail-placed end helps mouse ES and hurts WTC11 (row 1860); the 3' extension from unspliced reads is opt-in (1488).
- Step 0 of today (`bench/o2_truth/end_accuracy.py`, nearest annotated end of the same chain): on simulated fly / mouse
  Rustle's ends are the best of the four tools (0.97 / 0.98 within 10 bp at both ends). On mouse ES its 5' end is like
  StringTie's and IsoQuant's (0.49 within 50 bp; past the annotated start by > 50 bp in 0.45) and behind FLAIR's (0.68);
  its 3' end is the worst of the four (within 10 bp 0.46 against 0.52 StringTie, 0.55 FLAIR, 0.59 IsoQuant; past the
  annotated end by > 50 bp in 0.37 against 0.15-0.22).

## What is tested: an end that is extended AND supported
The median cuts real full-length starts when most reads are cut short (simulated sets, heart); the most extended end
follows one stray read. Candidates per chain, from the reads with exactly the chain (`bench/o2_truth/end_choice.py`):
- `ext` most extended (today); `med` lower median (row 1487); `mode` centre of the densest 20 bp window;
- `sup2` the most extended read end that has another read end within 10 bp (none: `mode`);
- `sup10` the most extended read end with at least max(2, 10% of the chain's reads) ends within 10 bp (none: `mode`).
Judges: human WTC11 r1: 5' in a refTSS CAGE peak, 3' within 25 bp of a polyA site of the atlas (both lifted to CHM13),
and the nearest annotated end of the chain; mouse ES, heart: annotation; simulated fly, mouse: the truth; gorilla OR6737: RefSeq
(reported, expected to favour `ext` for the reason above).
- **D1 (5'):** `sup2` or `sup10` is in a CAGE peak at least 0.20 more often than `ext` on WTC11 and within 0.05 of `mode`.
- **D2 (no cost on complete reads):** on both simulated sets that candidate is within 0.005 of `ext` at 10 bp, at both ends.
- **D3 (3'):** that candidate is at a polyA site more often than `ext` on WTC11 and within 10 bp of the annotated end at
  least 0.05 more often than `ext` on mouse ES.
- **D4:** on mouse heart it is not below `ext` at 50 bp (the median was, by 2 points).
Prediction: D1 holds for `sup10`, D2 holds for both, D3 holds weakly, D4 holds; on gorilla RefSeq both lose to `ext`.

## Result on the six libraries (job 4517236)
Signal = within 25 bp of a CAGE peak (5') / a polyA site (3'); background 500 bp away 0.08 / 0.01. Annotated chains: share
within 10 bp / 50 bp of the nearest annotated end of the chain.
| library | end | `ext` (today) | `med` | `mode` | `sup2` | `sup10` |
|---|---|---|---|---|---|---|
| WTC11 r1 | 5' signal | 0.755 | 0.855 | 0.855 | 0.818 | 0.847 |
| WTC11 r1 | 5' annotated | 0.221 / 0.611 | 0.401 / 0.824 | 0.418 / 0.817 | 0.289 / 0.716 | 0.365 / 0.795 |
| WTC11 r1 | 3' signal | 0.717 | 0.776 | 0.787 | 0.767 | 0.788 |
| WTC11 r1 | 3' annotated | 0.605 / 0.677 | 0.569 / 0.613 | 0.587 / 0.628 | 0.635 / 0.698 | 0.638 / 0.689 |
| mouse ES | 5' annotated | 0.142 / 0.495 | 0.222 / 0.665 | 0.234 / 0.655 | 0.191 / 0.599 | 0.267 / 0.682 |
| mouse ES | 3' annotated | 0.463 / 0.560 | 0.534 / 0.611 | 0.552 / 0.620 | 0.533 / 0.630 | 0.587 / 0.659 |
| mouse heart | 5' annotated | 0.249 / 0.653 | 0.207 / 0.630 | 0.216 / 0.625 | 0.271 / 0.690 | 0.277 / 0.688 |
| mouse heart | 3' annotated | 0.552 / 0.641 | 0.517 / 0.618 | 0.528 / 0.627 | 0.577 / 0.666 | 0.582 / 0.663 |
| fly simulated | 5' truth | 0.967 / 0.994 | 0.643 / 0.986 | 0.924 / 0.986 | 0.942 / 0.989 | 0.943 / 0.989 |
| fly simulated | 3' truth | 0.979 / 0.999 | 0.855 / 0.999 | 0.845 / 0.999 | 0.971 / 0.999 | 0.971 / 1.000 |
| mouse simulated | 5' truth | 0.971 / 0.993 | 0.637 / 0.988 | 0.937 / 0.987 | 0.951 / 0.989 | 0.952 / 0.989 |
| mouse simulated | 3' truth | 0.979 / 0.995 | 0.906 / 0.998 | 0.905 / 0.998 | 0.976 / 0.996 | 0.976 / 0.997 |
| gorilla OR6737, RefSeq | 5' | 0.342 / 0.599 | 0.081 / 0.359 | 0.097 / 0.369 | 0.224 / 0.523 | 0.168 / 0.471 |
| gorilla OR6737, RefSeq | 3' | 0.518 / 0.593 | 0.636 / 0.710 | 0.655 / 0.726 | 0.583 / 0.659 | 0.658 / 0.729 |
- **D1 FAILS by its bar:** `sup10` is near a CAGE peak 0.847 against 0.755 (+0.09, bar +0.20; within 0.05 of `mode`). Row 1718's
  0.417 for `ext` was inside the peak itself; with 25 bp around it today's rule already does better than I assumed.
- **D2 FAILS:** on simulated reads `sup10` is 0.019-0.024 under `ext` at the 5' end and 0.003-0.008 at the 3' end (bar 0.005).
- **D3 HOLDS:** polyA site 0.717 -> 0.788 on WTC11; mouse ES 3' within 10 bp 0.463 -> 0.587.
- **D4 HOLDS:** mouse heart 5' 0.653 -> 0.688 at 50 bp (the median 0.630).
- Gorilla RefSeq: 5' loses as predicted (0.599 -> 0.471; the median 0.359); **3' GAINS** (0.518 -> 0.658 at 10 bp), not predicted.
Reading: `sup10` is the best or within 0.02 of the best on every REAL library at both ends, and unlike the median it does
not fall under today's rule on mouse heart. The 3' end is the clear case: better on all four real libraries including the
ape, for 0.003-0.008 on simulated reads. The 5' end gains on human and mouse and loses on ape RefSeq and 0.02 on simulated.

## Confirmation on unused libraries, rule fixed (`sup10`, window 10 bp)
WTC11 r2 (CAGE, polyA), mouse adrenal, mouse cortex, gorilla KB3781, zebrafish simulated.
- **X1 (3'):** within 10 bp of the annotated end above `ext` by >= 0.02 on WTC11 r2, adrenal, cortex and KB3781; at a polyA
  site more often on WTC11 r2; zebrafish simulated not more than 0.01 below `ext`.
- **X2 (5'):** WTC11 r2 near a CAGE peak >= +0.05 and within 50 bp of the annotated start >= +0.10; adrenal and cortex not
  below `ext` at 50 bp; KB3781 (RefSeq) reported, expected lower; zebrafish simulated reported.

## Confirmation result (job 4517322), `ext` -> `sup10`
| library | 5' signal | 5' annotated 10 / 50 bp | 3' signal | 3' annotated 10 / 50 bp |
|---|---|---|---|---|
| WTC11 r2 | 0.762 -> 0.836 | 0.252 / 0.657 -> 0.378 / 0.809 | 0.746 -> 0.795 | 0.618 / 0.686 -> 0.614 / 0.666 |
| mouse adrenal | | 0.235 / 0.645 -> 0.279 / 0.694 | | 0.541 / 0.634 -> 0.567 / 0.648 |
| mouse cortex | | 0.197 / 0.564 -> 0.254 / 0.662 | | 0.523 / 0.605 -> 0.588 / 0.666 |
| gorilla KB3781, RefSeq | | 0.352 / 0.613 -> 0.167 / 0.473 | | 0.490 / 0.564 -> 0.683 / 0.749 |
| zebrafish simulated | | 0.966 / 0.995 -> 0.938 / 0.987 | | 0.979 / 0.996 -> 0.974 / 0.997 |
- **X1 HOLDS on 3 of 4 annotated judges** (adrenal +0.026, cortex +0.065, KB3781 +0.193; WTC11 r2 -0.004 by the annotation
  but +0.049 at polyA sites); zebrafish -0.005.
- **X2 HOLDS:** WTC11 r2 near a CAGE peak +0.074, within 50 bp of the annotated start +0.152; adrenal +0.049, cortex +0.098.
  KB3781 RefSeq starts 0.613 -> 0.473 (expected); zebrafish simulated -0.028 at 10 bp.

## Built: `tools/end_choice.py`, opt-in `RUSTLE_SUPPORTED_ENDS=3|5|both` in `all` (job 4517340)
The tool reproduces the bench numbers; chains, transcript and line counts identical, exons valid (ES, fly, cortex).
Mouse ES against the other tools (nearest annotated end of the chain; 5' within 50 bp, 3' within 10 bp):
Rustle today 0.494, 0.462; with both ends moved 0.682, 0.586; FLAIR 0.677, 0.547; IsoQuant 0.492, 0.593; StringTie 0.482, 0.516.
Fly simulated: Rustle today 0.967 / 0.979 at 10 bp, with both ends 0.943 / 0.971; IsoQuant 0.971 / 0.913, StringTie 0.963 / 0.950, FLAIR 0.616 / 0.850.

## Families with the rule on at both ends (job 4517354)
Planted families, seeds 1998 / 1999, 3 and 60 reads per copy, E and S: 20-135 chains moved per set; copies found, false
copies, isoforms found, transcripts, exact transcripts and the three label counts are IDENTICAL with and without the rule
in all eight sets (main assembly and the mode). Gorilla KB3781: 35,025 of 60,287 chains moved (5' 29,515; 3' 17,454);
families with >= 2 copies 260 -> 260, the same sets of transcripts, 724 copies. cargo test 1,591 passed.
(Rows 1490 / 1491 had found changed families on A119b with the median 5' end and the 3' EXTENSION; this rule only moves
ends inward. A119b itself was not rerun.)
State: opt-in `RUSTLE_SUPPORTED_ENDS=3|5|both`; default unchanged (user's call).

## Decision (user, 2026-11-28): the 3' end uses the rule by default
`RUSTLE_SUPPORTED_ENDS` unset = 3; `5` / `both` opt-in; 0 = off. Check (job 4517502, fly through the whole `all`): 196 of
7,221 chains moved at the 3' end, none at the 5' end; every end equals the bench tool's output (7,221 of 7,221); with 0 no
transcript is tagged; a second run redoes nothing; no warning, `PREFIX.final.gtf` made; cargo test 1,591 passed.
Seen in the same table: the one-read chains of row 1990 (401 on fly) take their ends from a single read, so the starts of
`PREFIX.gtf` as a whole are 0.943 within 10 bp of the truth against 0.967 before those chains existed; with the rule off
the 3' end is 0.972, with it 0.964 on this simulated set (the known simulated cost).
Every `PREFIX.gtf` made with `all` before this date has the assembler's 3' ends.

## Addendum before running (2026-11-28, user: how do we compare on transcripts that do NOT begin or end at the expected TSS / TES?)
A transcript end away from the annotated one is either a real alternative start / polyA site or an artifact (a cut-short
or over-extended read). The annotation cannot tell; CAGE peaks and the polyA-site atlas can (human, CHM13).
HG002 (`reads/HSA/HG002_kinnex.t2t.bam`; StringTie, FLAIR, IsoQuant de novo assemblies of the same reads exist), Rustle
with ends off, with the 3' rule (today's default) and with both. `bench/o2_truth/unexpected_ends.py`: transcripts with an
annotated chain, each end classed against the nearest annotated end of the chain: expected (within 50 bp), short (> 50 bp
inside), past (> 50 bp beyond); per class the share within 25 bp of a CAGE peak (5') / polyA site (3'). Also the
transcripts with an unannotated chain.
- **U1:** Rustle's ends that are `past` the annotated one have the signal less often than its `expected` ends by >= 0.30
  (they are mostly over-extension), at both ends, with the rule off.
- **U2:** the rule lowers the NUMBER of `past` ends by >= 30% at the end it is applied to and raises the signal share of
  the ends that remain away from the expected one (short + past together) by >= 0.05.
- **U3:** with the rule at both ends, Rustle's away-from-expected ends carry the signal at least as often as the best of
  the three other tools minus 0.05, at both ends.
Prediction: U1 holds at the 5' end and fails at the 3' end (ends past the annotated one are often real distal polyA sites);
U2 holds; U3 holds at the 3' end, fails at the 5' end (FLAIR's supported windows).

## Result of the addendum (job 4517844, HG002): ends away from the annotated ones
Transcripts with an annotated chain; per class: share of the transcripts, and (signal) the share within 25 bp of a CAGE
peak / polyA site. Background 0.08 / 0.01.
| tool | 5' expected | 5' short | 5' past | 5' away: signal | 3' expected | 3' short | 3' past | 3' away: signal |
|---|---|---|---|---|---|---|---|---|
| Rustle, ends off | 0.504 (0.906) | 0.062 (0.565) | 0.434 (0.440) | 0.455 | 0.545 (0.845) | 0.251 (0.439) | 0.204 (0.050) | 0.265 |
| Rustle, default (3') | as above | | | 0.455 | 0.579 (0.883) | 0.385 (0.575) | 0.036 (0.096) | 0.534 |
| Rustle, both ends | 0.774 (0.963) | 0.104 (0.656) | 0.122 (0.728) | 0.695 | 0.579 (0.883) | 0.385 (0.575) | 0.036 (0.096) | 0.534 |
| StringTie | 0.619 (0.938) | 0.059 (0.667) | 0.323 (0.391) | 0.433 | 0.508 (0.892) | 0.419 (0.565) | 0.073 (0.052) | 0.488 |
| FLAIR | 0.799 (0.975) | 0.156 (0.740) | 0.045 (0.703) | 0.732 | 0.500 (0.900) | 0.486 (0.554) | 0.014 (0.107) | 0.541 |
| IsoQuant | 0.426 (0.862) | 0.050 (0.535) | 0.524 (0.397) | 0.409 | 0.493 (0.809) | 0.153 (0.338) | 0.355 (0.041) | 0.130 |
Annotated chains 29,648 (Rustle), 27,456, 38,678, 16,088. Unannotated chains, 5' / 3' signal: Rustle off 0.447 / 0.345,
default 0.447 / 0.364, both 0.487 / 0.364 (92,032 chains); StringTie 0.385 / 0.271 (88,123); FLAIR 0.629 / 0.404
(690,596); IsoQuant 0.415 / 0.178 (51,622).
- **U1 HOLDS at both ends** (my 3' prediction was wrong): a 3' end PAST the annotated one is at a polyA site 5% of the
  time, in every tool; it is over-extension, not a distal polyA site. A 5' end past the annotated one has a CAGE peak 0.44.
- **U2 HOLDS:** the rule removes 82% of the past 3' ends (6,049 -> 1,082) and 72% of the past 5' ends (12,879 -> 3,615);
  the ends that stay away from the expected one carry the signal 0.265 -> 0.534 (3') and 0.455 -> 0.695 (5').
- **U3 HOLDS at both ends** (my 5' prediction was wrong): 0.534 against FLAIR's 0.541, 0.695 against FLAIR's 0.732.
- A 3' end SHORT of the annotated one is at a polyA site in 0.55-0.58 of the cases for Rustle, StringTie and FLAIR: more
  than half are real proximal polyA sites. Rustle reports one model per chain, so a chain with two real ends keeps one;
  FLAIR reports up to two. Not addressed here (`altends` exists as an opt-in stage).
All 3' ends of annotated-chain transcripts at a polyA site: Rustle off 0.581, default 0.736.
