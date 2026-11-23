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
