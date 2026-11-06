# PREREG: are the unannotated families signal or noise? (row 1556, 2026-11-06 23:15, user-requested)

## Question
212 of A119b's 1,565 families have no copy on any GENCODE gene (all biotypes, lncRNA and pseudogene included);
they hold 8,994 of the 12,294 assigned reads and the assign wall. Are they transcribed gene families or
alignment artefacts of a repeat block?

## Checks and thresholds (fixed before looking)
1. Splice motif: for every copy's junction (two-exon models: one intron each; multi-exon: all introns), the
   dinucleotides at the intron ends on the copy's strand. Canonical = GT..AG (or GC..AG). Per group, the fraction
   of canonical introns. Noise mark if skipped < 0.50 while kept >= 0.90; signal mark if skipped >= 0.80.
2. Repeat context: fraction of copies (by count) overlapping the T2T cenSat v2.0 intervals (satellite arrays,
   rDNA, pericentromeric), downloaded to /anvme/workspace/v138da13-rustle/r1556/dl/censat.bed. Noise mark if
   skipped >= 0.50 while kept <= 0.10.
3. Decisiveness: median n_decisive and median margin of assigned reads, skipped vs kept (gw_off assignments).
   Noise mark if skipped median n_decisive <= 2 or skipped median margin < half of kept's.
Verdict: NOISE if >= 2 noise marks and no signal mark; SIGNAL if check 1 gives the signal mark and check 3 gives
no noise mark; otherwise UNDETERMINED, and check 4 (families on a second sample) is the next row.
## Prediction
Noise: motif ~0.35 (row 1519's hotspot figure), cenSat overlap high for the chr13 block, decisiveness low (the
hot9 yield of 0.1% suggests one-or-two-site decisions). Expected NOISE verdict; if SIGNAL, the annotated-only
restriction (row 1551) must not be adopted.
