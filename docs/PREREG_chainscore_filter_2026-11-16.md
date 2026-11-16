# PREREG: the frozen r1677 chain score as a filter with a fixed cut, opt-in stage `chainscore` (2026-11-16, register 1761)

User (2026-11-16, after row 1760): "do #2" = build the score as a filter stage, cut fixed on spent samples, tested as its own row.
Written before the cut is computed and before the score has been applied to HG002.

## Rule (nothing below is chosen after seeing a test sample)
- Features: `chain_rules_feat.py`, all columns, standardised within the sample (`chain_rules.design`), default multi-exon chains only.
- Weights: the row 1677 fit, logistic, pooled over the default chains of OR6737 and A119b
  (`r1677_chainrules/feat.{or,a}.tsv`), the same fit rows 1747 / 1751 / 1760 scored as "FROZEN pooled".
- Cut: a chain is dropped when its fitted probability is below TAU. TAU = the largest value (grid 0.001) at which OR6737 and A119b
  each keep >= 0.95 of their correct default chains. Computed once from those two tables and stored with the weights in
  `tools/data/chainscore_r1677.json`.
- Single-exon transcripts are never dropped.

## Test, one sample, once
HG002 Kinnex on CHM13 (`r1692_trunc5_hg/feat.hg.tsv`: 153,033 default chains, 35,745 correct by CHM13 RefSeq). The feature table exists
since row 1692 (`trunc5`); the score has never been computed on it. Bars:
1. correct chains kept >= 0.93 of the default's;
2. precision of the kept set minus the read-count curve at the same number of correct chains >= +0.02.
Both = PASS. Reported beside it, not a test (all seen by the score curve): the same two numbers on OR6737, A119b, KB3781, the four
testis apes and the four cell-line samples.

## What a verdict means
The stage `chainscore` (driver, opt-in, `PREFIX.chainscore.gtf`; `PREFIX.gtf` untouched) is built either way, as `trunc5` and
`shiftdrop` were. PASS = the filter is put to the user as a default candidate for `PREFIX.final.gtf`. FAIL = it stays opt-in.

## Limits stated now
- HG002 shares genome, annotation and species with A119b, one of the two fit samples; it is a different individual, tissue (cell line vs
  testis) and library (Kinnex).
- Every table here is an assembly made before transfrag extension became the default (2026-11-16). The stage will run on extended
  assemblies; one confirming run on a new-default assembly is reported with the result, without a bar.
- A probability cut assumes the share of wrong chains is comparable across samples; the reported table shows how far kept recall moves.

## RESULT (register 1761; `r1761_chainscore/eval.txt`, `confirm.4473789.out`)
TAU = 0.061. Fit samples (in-sample): OR6737 keeps 0.951 of correct chains, 0.352 -> 0.413 (+0.040 over the read cut); A119b 0.978,
0.176 -> 0.214 (+0.031).
| sample | default chains / correct / precision | kept | correct kept (share) | precision kept | read cut there | difference |
|---|---|---|---|---|---|---|
| **HG002 (the test)** | 153,033 / 35,745 / 0.234 | 124,669 | 34,727 (**0.972**) | 0.279 | 0.253 | **+0.026** |
| KB3781 fibroblast | 72,102 / 25,845 / 0.358 | 60,424 | 24,110 (0.933) | 0.399 | 0.396 | +0.003 |
| PTR testis | 53,970 / 20,982 / 0.389 | 42,880 | 20,033 (0.955) | 0.467 | 0.403 | +0.064 |
| PPA testis | 56,149 / 18,096 / 0.322 | 46,895 | 17,142 (0.947) | 0.366 | 0.341 | +0.025 |
| PPY testis | 104,858 / 20,182 / 0.192 | 81,781 | 19,630 (0.973) | 0.240 | 0.198 | +0.042 |
| PAB testis | 70,461 / 16,067 / 0.228 | 54,180 | 15,631 (0.973) | 0.289 | 0.233 | +0.055 |
| PPY fibroblast | 77,528 / 23,491 / 0.303 | 62,047 | 22,664 (0.965) | 0.365 | 0.320 | +0.046 |
| PAB fibroblast | 70,941 / 24,673 / 0.348 | 58,562 | 23,259 (0.943) | 0.397 | 0.379 | +0.018 |
| PPA fibroblast | 97,744 / 20,800 / 0.213 | 76,669 | 20,380 (0.980) | 0.266 | 0.218 | +0.047 |
| PTR lymphoblastoid | 89,021 / 30,351 / 0.341 | 75,167 | 28,758 (0.948) | 0.383 | 0.362 | +0.021 |
**Verdict: PASS as registered** (HG002 0.972 >= 0.93, +0.026 >= +0.02). The filter removes 15-23% of the chains for 2-7% of the
correct ones. On the seen samples it clears both bars on 8 of 10; KB3781 (+0.003) and PAB fibroblast (+0.018) are the two where it is
little better than raising the read floor. Stage `chainscore` built opt-in (driver, test, README); on PAB fibroblast it reproduces the
table row exactly, and on that sample's extended (new default) assembly it keeps 0.942 of the correct chains, 0.346 -> 0.387, +0.013
over the read cut. Default candidate for `PREFIX.final.gtf`; not switched.

**User decision 2026-11-16: default.** `assemble` computes the score on `PREFIX.gtf` and takes the same chain list out of
`PREFIX.final.gtf` (`RUSTLE_CHAINSCORE=0` = off; a failing filter leaves the layered assembly unfiltered and says so). `PREFIX.gtf` and
the families input are byte-identical with and without it. PTR testis, both new defaults, job 4473791: `final.gtf` 51,913 chains /
21,753 correct / 0.419 unfiltered, 45,418 / 20,993 / 0.462 filtered. Not measured: the filter on `final.gtf` of any other sample, and the
weights were fitted on assemblies made before transfrag extension.
