# PREREG: single-exon models admitted on family evidence (row 1919, 2026-11-23)

User (2026-11-23): "can we try to see if we can cleverly admit more single exon transcripts to reduce the edge
stringtie and flair have over rustle?"
Written before any run of the new switch. Seen: rows 1831, 1832 (lowering the single-exon floor genome-wide buys
1.7-3.7 x the single-exon truth at half to a quarter of the precision on every real library: no value qualifies);
rows 1897, 1901 (single-exon models classed against multi-exon ones; `--polish-mono-shadow-spare` opt-in: +11 to +33%
annotated single-exon genes, lower single-exon precision); row 1918 (human named families, expressed single-exon copies
found: Rustle 17 of 22, StringTie 20, FLAIR 21); and today, on the gorilla own-family planted set of row 1917:
expressed single-exon copies found Rustle primaries 119 of 147, the mode 92, FLAIR 146, StringTie / IsoQuant / bambu 0;
half silent 58 / 53 / 77 of 77 and silent single-exon copies called 1 / 38 (9 with the EM test) / 12 of 70.
So the edge to close is FLAIR's (and StringTie's on the human set); on gorilla StringTie writes no single-exon model.

## The rule (`RUSTLE_MULTICOPY_MONO=1`, opt-in; `tools/multicopy_rescue.py candidates --mono`)
The read-count floor for single-exon models (0.82 quantile of the multi-exon read counts of the contig) asks "is this
well expressed"; a family asks "is this a copy of something else we assembled". Inside stage `multicopy`: the floor-1
run is made without the single-exon floor (`RUSTLE_POLISH_MONO_QUANTILE=0`; the shadow and internal-priming rules and
the built-in 2-read floor stay); its single-exon transcripts at places the assembly leaves empty are candidates; a
candidate is kept only when the families stage puts its locus in a family of >= 2 loci (`multicopy_rescue "mono"`).
Then the provenance labels as for every locus. Nothing changes without the switch.

## Runs (workspace `r1919_mono_family_rescue`)
- Planted, fresh seed 1925, ragged ends: gorilla OR6737's own families (E, S) and the human named families (E, S);
  arms mode, mode + mono rescue, each with the EM test; StringTie, FLAIR, IsoQuant, bambu on the same BAMs.
- Real: OR6737 and A119b, stage `multicopy` with the switch on the assemblies of rows 1888 / 1869. For every set of
  single-exon models (the default's, the rescued ones, StringTie's, FLAIR's): how many, and the share with >= 50% of
  their bases on an exon of an annotated gene of that strand (RefSeq gorilla, CAT/Liftoff human).

## Claims
- M1 (the edge closes): gorilla planted, expressed single-exon copies found by mode + mono rescue + EM >= 0.9 x FLAIR's,
  in E and in S.
- M2 (not by calling everything): gorilla S, silent single-exon copies called by mode + mono rescue + EM <= FLAIR's.
- M3 (nothing else moves): multi-exon copies found and isoform recall of mode + mono rescue within 0.005 of the mode's.
- M4 (real data, the precision of what is added): on OR6737 and on A119b the rescued single-exon models' annotated
  share is >= the default's single-exon models' share minus 0.05.
M1-M3 = it works on planted reads. M4 decides whether it is worth proposing beyond an opt-in; whatever falls it stays
off by default until the user says otherwise.

## Predictions (mine)
M1 holds (the planted single-exon copies are family members by construction and the families stage will group them).
M2 holds with the EM test. M3 holds. M4 fails on at least one sample (row 1897: unpaired single-exon models are
annotated 6-43% of the time; family membership should lift that, I expect 30-60% against 60-80% for the default's).

## Result (2026-11-23, jobs 4496561-4496563): M1 FAILS, M2 HOLDS, M3 holds for multi-exon copies, M4 FAILS on gorilla. Stays an opt-in; not proposed
Planted, gorilla OR6737's own families, seed 1925 (expressed single-exon copies found / silent single-exon copies called):
| tool | E (147 expressed) | S (76 expressed, 71 silent) | isoform recall / precision, E |
|---|---|---|---|
| Rustle, primaries | 124 | 55 / 0 | 0.982 / 0.954 |
| mode | 101 | 52 / 36 | 0.966 / 0.914 |
| mode + EM | 100 | 52 / 8 | 0.961 / 0.917 |
| mode + single-exon rescue | 119 | 58 / 46 | 0.975 / 0.915 |
| mode + single-exon rescue + EM | 118 | 58 / 8 | 0.970 / 0.918 |
| StringTie | 0 | 0 / 0 | 0.848 / 0.967 |
| FLAIR | 146 | 76 / 11 | 0.943 / 0.962 |
| IsoQuant, bambu | 0 | 0 / 0 | 0.814 / 0.984, 0.901 / 0.740 |
- M1 fails: 118 of 147 and 58 of 76 against FLAIR's 146 and 76 (bar 0.9 x: 131 and 68). The rescue lifts the mode from
  101 to 119, to about where Rustle with primaries only already is (124).
- M2 holds: 8 silent single-exon copies called against FLAIR's 11 (46 without the EM test).
- M3: multi-exon copies found unchanged (1,181 / 635); isoform recall + 0.009 (E) and + 0.005 (S) from the added models,
  precision unchanged.
- Human named families: only 2 (E) and 1 (S) copies are single-exon in every planted isoform; the rescue finds none
  (StringTie 2 and 0, FLAIR 2 and 1).
Real samples, single-exon models and the share with >= 50% of their bases on an annotated exon of their strand:
| set | OR6737 | A119b |
|---|---|---|
| the default's single-exon models | 534, 56.4% | 10,844, 10.8% |
| rescued on family evidence | 728, 13.3% | 4,317, 6.9% |
| every single-exon model of the floor-free run | 3,617, 22.8% | 29,687, 7.1% |
| StringTie | 758, 50.8% | 38,635, 3.8% |
| FLAIR | 14,850, 22.0% | 468,993, 5.3% |
- M4 fails on OR6737 (13.3% against 56.4%) and holds by the letter on A119b (6.9% against 10.8%, within 0.05) where
  every set is poorly annotated. The family filter does not select annotated single-exon models: on gorilla the rescued
  ones are annotated LESS often than the unfiltered floor-free ones (13.3% against 22.8%).
- Named copies (toolcompare) unchanged by the rescue: OR6737 34 / 0.66, A119b 126 / 0.49.
Predictions: M1 wrong (I expected it to hold), M2 right, M3 right, M4 as feared and worse (13%, I said 30-60%).
Reading: membership in a de novo family is not evidence that a single-exon model is a gene; families of single-exon
loci are largely repeats and processed fragments. FLAIR's planted lead comes with 22% (gorilla) and 5% (human)
annotated single-exon models on real data, 15 k and 469 k of them: it is not a standard to match. On real gorilla
Rustle's default single-exon models are the most often annotated of the three tools.

# PREREG: the single-exon models primaries alone would pass, restored in the mode (row 1920, 2026-11-23)

Seen: row 1919 (planted gorilla, seed 1925: expressed single-exon copies found Rustle primaries 124, the mode 101;
real OR6737: the mode's default single-exon models 534 at 56.4% on an annotated exon). Idea told to the user before any
run: the single-exon floor is a quantile of the multi-exon read counts, near-tie secondaries inflate those counts, so
the mode drops single-exon models the same reads pass when only primaries are counted. Restoring them admits nothing a
default Rustle run would not write.
First as an emulation, no assembler change: `multicopy_rescue.py mono_primary MODE.gtf PRIMARY.gtf OUT.gtf` adds the
single-exon transcripts of the primaries-only assembly at places the mode's assembly has no transcript on that strand.
Not seen: any output of it.

## Runs (workspace `r1920_mono_from_primaries`)
- Planted gorilla OR6737's own families, fresh seed 1926 (E, S), all tools; arms mode, mode + primaries' single-exon
  models, each with the EM test.
- Real OR6737 and A119b: both pools assembled today, stage `multicopy`, then the addition. The added models: how many,
  share with >= 50% of bases on an annotated exon of their strand, against the mode's own single-exon models.

## Claims
- Q1: planted, E and S: mode + primaries' single-exon models + EM finds at least (the primaries arm's count - 2)
  expressed single-exon copies.
- Q2: S: silent single-exon copies called (with EM) <= FLAIR's.
- Q3: real, both samples: the added models' annotated share >= the mode's default single-exon models' share - 0.05.
- Q4: multi-exon copies found and isoform recall unchanged (within 0.005; isoform recall may only rise).
All four = worth building properly (the floor taken from primaries inside the assembler, or the extra assembly inside
the stage) and proposing; the default stays the user's call.

## Predictions (mine)
All four hold. The risk is Q3 on A119b, where the default's single-exon models are annotated only 10.8% of the time
and small shares move easily.

## Result (row 1920, 2026-11-23; jobs 4496575 planted, 4496576 OR6737, 4496577 A119b)
Planted gorilla, seed 1926 (expressed single-exon copies found / silent called):

| arm | E (147 expressed) | S (75 expressed, 72 silent) | isoform recall E / S |
|---|---|---|---|
| Rustle primaries | 120 / 0 | 59 / 1 | 0.983 / 0.977 |
| mode + EM | 98 / 0 | 51 / 5 | 0.962 / 0.970 |
| mode + primaries' single-exon + EM | 119 / 0 | 58 / 5 | 0.973 / 0.977 |
| FLAIR | 147 / 0 | 75 / 10 | 0.943 / 0.938 |
| StringTie, IsoQuant, bambu | 0 | 0 | |

Real samples, models added and their share with >= 50% on an annotated exon:

| sample | added | annotated | the mode's own single-exon models | primaries assembly's |
|---|---|---|---|---|
| OR6737 | 1 | 1 (1.000) | 534 at 0.564 | 311 at 0.598 |
| A119b | 224 | 6 (0.027) | 10,844 at 0.108 | 9,369 at 0.044 |

- Q1 HOLDS (119 against 120 - 2; 58 against 59 - 2). Q2 HOLDS (5 against FLAIR's 10). Q4 HOLDS (multi-exon found 1,173 = 1,173
  and 637 = 637; isoform recall rises 0.962 -> 0.973, 0.970 -> 0.977).
- Q3 FAILS on A119b (0.027 against the bar 0.108 - 0.05 = 0.058); on OR6737 it holds on ONE model.
- **Verdict: three of four, not built, nothing proposed.** My prediction (all four) was wrong at the place I named as the risk.
- What it shows beyond the rule: the loss the planted benchmark measures (98 against 120) hardly exists in the real
  samples. On OR6737 the mode drops one single-exon model the primaries would pass; on A119b 224, of which 6 are
  annotated. The planted families are dense in single-exon copies with near-tied relatives, real libraries are not.
  The single-exon gap to FLAIR on planted reads should be quoted as a benchmark figure, not as lost genes.
- Side fact: on A119b the mode's single-exon models are MORE often annotated (0.108) than the primaries assembly's
  (0.044), so seeding near ties does not degrade that class there.
- The emulation `multicopy_rescue.py mono_primary` stays in the tool (self-test covered), called by no stage.

### A119b run in full (row 1920, job 4496577: every step exit 0, 8,597 s; primaries assembly 2,368 s, with near-tie secondaries 2,594 s, families + multi-copy stage 3,530 s)

| set | single-exon models | >= 50% on an annotated exon | share | distinct annotated genes |
|---|---|---|---|---|
| primaries-only assembly | 9,369 | 412 | 0.044 | 408 |
| multi-copy mode | 10,844 | 1,176 | 0.108 | 522 |
| added back from primaries | 224 | 6 | 0.027 | 6 |
| StringTie | 38,635 | 1,473 | 0.038 | 1,375 |
| FLAIR | 468,993 | 24,897 | 0.053 | 10,697 |

The mode's single-exon models have the highest annotated share of any set on this sample; StringTie and FLAIR reach
more annotated genes in absolute numbers by reporting 3.6 x and 43 x as many models. OR6737 (job 4496576, 2,316 s):
primaries 311 at 0.598, mode 534 at 0.564, StringTie 758 at 0.508, FLAIR 14,850 at 0.220.
