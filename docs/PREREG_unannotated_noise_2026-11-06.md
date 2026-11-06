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

## Result (2026-11-06 23:35): SIGNAL; prediction (NOISE) FAILED; row 1551 must not be adopted
Correction before the numbers were read: the first motif pass used 1-based acceptor coordinates on 0-based exon
starts (0/9943 canonical in BOTH groups, top motif "GT-CA" = AG shifted by one); fixed and rerun; check 1 below is
the corrected pass.
1. Splice motif: kept 9,861/9,943 canonical (0.99); skipped 831/951 (0.87) >= 0.80 = SIGNAL mark. The skipped
   remainder is 120 AT-AC introns (12.6%; the U12 minor-spliceosome motif, itself a splice signal; genome-wide
   AT-AC is ~0.3%), so 951/951 of the skipped introns carry a known splice motif. The nine hotspot families:
   75 of 85 introns GT/GC-AG.
2. cenSat: skipped 1,228/1,388 copies (0.88) in cenSat intervals; kept 2,840/6,599 (0.43). Kept > 0.10, so no
   noise mark by the registered rule (the v2.0 track's broad ct/p_arm intervals cover many ordinary genes; a
   satellite-only split was not computed: post hoc).
3. Decisiveness: skipped assignments median n_decisive 9, median margin 34.5; kept 3 and 6.9. No noise mark;
   the unannotated assignments rest on MORE sites and wider margins than the annotated ones.
Verdict by the registered rule: SIGNAL (check 1 signal mark, check 3 no noise mark). The 212 unannotated
families are spliced at canonical and U12 motifs in a pericentromeric block and their 8,994 assignments are
the most decisive in the run. The annotated-only restriction (row 1551) trades away the pipeline's strongest
results and stays opt-in, not recommended. The assign cost is where the results are.
