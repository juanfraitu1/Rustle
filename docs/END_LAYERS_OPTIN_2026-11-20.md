# Row 1846: the three end layers of `final.gtf` become opt-in (user decision, 2026-11-20)

User (2026-11-20): "lets turn all those defaults you mentioned into opt ins" (the single-exon model layer, the 3' UTR
extension and the median 5' end). Basis: row 1831 (`docs/PREREG_mono_end_rules_2026-11-20.md`): mono-model pays on the two ape
libraries and fails on zebrafish, mouse ES and heart; the 3' extension has an effect on ape testis only; the median 5' end
helps A119b and mouse ES and hurts the other libraries; row 1834: no cap-signal gate predicts where it helps.
Change: `tools/rustle_pipeline.sh`, new `RUSTLE_FINAL_END_LAYERS` (`tssend`, `tesext`, `monomodel`, `all`; unset = none).
`RUSTLE_FINAL_END_LAYERS=all` = the output before. The sub-chain layer of `final.gtf` (`RUSTLE_FINAL_SUBCHAIN`), the chain
score filter and the 5' truncation rule are untouched. `PREFIX.gtf` is not affected by construction. Workspace
`r1846_end_layers_optin`.

## Check written before the runs (not a bar for the decision, which is made)
Row 1831 switched one layer off at a time, on seven datasets, before the junction replacement was a default. Runs: the new
default and `RUSTLE_FINAL_END_LAYERS=all` on the twelve datasets of rows 1837-1845. Reported for `x.final.gtf` with row 1831's
scorer: multi-exon chains and correct chains, single-exon models / precision / truth found, 5' and 3' ends within 50 bp among
correct chains, exonic base sensitivity and precision. Parity: the `all` arm equals row 1845's `jr20` byte for byte; `x.gtf`
is identical between the arms.
Flag to the user: a dataset where the new default loses more than 0.5% of the correct chains or more than 1 point of chain
precision in `final.gtf`; and the size of what the ape libraries give up (ends, bases, single-exon), which row 1831 predicts.
Known interaction: the chain score filter and the truncation rule read end features; row 1831 saw 3-49 MORE correct chains
with the median 5' end off.

## Result (2026-11-20): chains unaffected; ends and single-exon models move as row 1831 said; A119b gives up most

Jobs 4490847-4490858, 24 runs exit 0; `r1846_end_layers_optin/score.NAME.txt` (row 1831's scorer). Parity: the `all` arm's
`x.final.gtf` equals row 1845's `jr20` byte for byte on all twelve; `x.gtf` is identical between the arms. Test suite 1,585
passed, 0 failed (driver test adapted).

`x.final.gtf`, `all` (before) -> new default (no end layer):
| dataset | correct chains / precision | 5' within 50 bp | 3' within 50 bp | base sens. | base precision | single-exon: truth found, model precision |
|---|---|---|---|---|---|---|
| fly simulated | 7,016 / 0.967 -> 7,020 / 0.967 | 0.979 -> 0.994 | 0.999 -> 0.999 | 0.684 -> 0.696 | 0.998 -> 0.996 | 342, 0.989 -> 343, 0.921 |
| mouse simulated | 27,041 / 0.959 -> 27,078 / 0.959 | 0.982 -> 0.994 | 0.995 -> 0.996 | 0.739 -> 0.737 | 0.998 -> 0.998 | 301, 0.969 -> 85, 0.966 |
| zebrafish simulated | 12,360 / 0.895 -> 12,363 / 0.896 | 0.980 -> 0.995 | 0.998 -> 0.998 | 0.758 -> 0.766 | 0.994 -> 0.994 | 80, 1.000 -> 144, 0.861 |
| mouse ES | 17,679 / 0.367 -> 17,703 / 0.366 | 0.673 -> 0.499 | 0.559 -> 0.559 | 0.356 -> 0.360 | 0.766 -> 0.760 | 270, 0.103 -> 275, 0.121 |
| mouse heart | 13,069 / 0.590 -> 13,086 / 0.589 | 0.634 -> 0.654 | 0.641 -> 0.640 | 0.255 -> 0.257 | 0.907 -> 0.899 | 128, 0.351 -> 59, 0.557 |
| mouse adrenal | 14,362 / 0.591 -> 14,380 / 0.590 | 0.641 -> 0.647 | 0.634 -> 0.633 | 0.259 -> 0.261 | 0.917 -> 0.906 | 152, 0.466 -> 65, 0.591 |
| mouse cortex | 15,656 / 0.500 -> 15,676 / 0.499 | 0.625 -> 0.567 | 0.606 -> 0.606 | 0.337 -> 0.340 | 0.842 -> 0.828 | 179, 0.509 -> 91, 0.596 |
| OR6737 | 25,565 / 0.384 -> 25,586 / 0.382 | 0.368 -> 0.602 | 0.607 -> 0.601 | 0.444 -> 0.453 | 0.692 -> 0.691 | 134, 0.222 -> 74, 0.216 |
| A119b | 37,936 / 0.221 -> 37,990 / 0.220 | 0.712 -> 0.460 | 0.596 -> 0.557 | 0.473 -> 0.467 | 0.410 -> 0.265 | 315, 0.169 -> 234, 0.021 |
| chimpanzee | 20,807 / 0.420 -> 20,818 / 0.418 | 0.427 -> 0.578 | 0.672 -> 0.656 | 0.333 -> 0.340 | 0.708 -> 0.687 | 145, 0.313 -> 104, 0.166 |
| bonobo | 17,840 / 0.340 -> 17,848 / 0.339 | 0.311 -> 0.463 | 0.604 -> 0.583 | 0.263 -> 0.268 | 0.691 -> 0.682 | 96, 0.231 -> 71, 0.217 |
| KB3781 | 25,667 / 0.380 -> 25,686 / 0.379 | 0.367 -> 0.615 | 0.583 -> 0.579 | 0.405 -> 0.411 | 0.675 -> 0.673 | 85, 0.232 -> 56, 0.247 |
(mouse rows here leave out the SIRV and ERCC contigs; row 1845's table kept the `chr` contigs only, hence 17,679 against 17,678.)

- **Chains: nothing to flag.** Every dataset gains correct chains (+3 to +54, at most +0.14%); chain precision moves by at most
  -0.002.
- **5' ends:** better without the layer on nine datasets (the four ape libraries other than A119b by 15-25 points, the
  simulated sets by 1.2-1.5, heart and adrenal slightly), worse on A119b (-25 points), mouse ES (-17) and cortex (-6), as rows
  1831 / 1834 found.
- **3' ends and bases:** A119b gives up 4 points of 3' ends within 50 bp and **15 points of base precision (0.410 -> 0.265)**;
  chimpanzee and bonobo 2 points of each; elsewhere under 1.5 points. A119b's base precision falls because the 3' layer and the
  single-exon model layer together had removed and rebuilt its single-exon models (now 10,848 models, 231 matching, against
  9,212 with 1,554).
- **Single-exon:** without the model layer fewer annotated single-exon transcripts are found on eight datasets (mouse
  simulated 301 -> 85, heart, adrenal, cortex about half, the apes by a quarter to a half), more on zebrafish (80 -> 144);
  model precision rises on the four mouse libraries and falls on fly, zebrafish, A119b and chimpanzee.
- SIRV contigs: unchanged (67 reported, 48 correct, 19 wrong).

**To flag for the user:** the decision costs A119b most (the library the layers were tuned on): 5' ends, base precision and
single-exon precision; `RUSTLE_FINAL_END_LAYERS=all` restores it for that library, `tesext monomodel` for ape testis in
general. The single-exon model layer is the least one-sided of the three: it finds more single-exon truth on eight of twelve
datasets. Keeping `monomodel` on by default while the other two stay opt-in is an alternative the user may prefer.
