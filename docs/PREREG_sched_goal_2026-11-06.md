# PREREG: targeted scheduling goal (row 1550, 2026-11-06 20:50, user-set goal)

## Goal (user's words)
assign wall on singlenode (72 cores, 250 GB) <= chr13 critical path + 10% (fetch + detect, ~4,550 s), six assign
tables byte-identical to gw_off, RSS < 250 GB; end to end from a clean directory <= 2 h, measured n = 3.

## Bars, fixed before any number
- A. assign wall <= 1.10 x (chr13 reads_in_region + chr13 detect_and_assign) of the SAME run; the ~4,550 s figure
  is the gw12/gw13 expectation (202-208 + 3,922-3,935), not the bar.
- B. six tables identical to /anvme/workspace/v138da13-rustle/asgsweep/gw_off/a1.assign.*; RSS < 250 GB.
- C. end to end (tools/rustle_pipeline.sh from a clean directory, staged BAM, the gw7 recipe + the adopted
  schedule) <= 7,200 s on each of three runs; assign tables of each run identical to gw_off.
## Mechanism under test
Row 1548 Amendment 2 (`RUSTLE_DISPATCH_DEFER=1`, bin13): chr13 holds the empty pool first, deferred regions cycle a
try-acquire queue. gw15 = bar A/B. If A/B pass, the flag becomes the default (unset path keeps the gw7 schedule
under RUSTLE_DISPATCH_DEFER=0) and three e2e runs follow for bar C.
## Prediction
gw15: chr13 held at t ~ 0, wall 4,400-4,900 s; A passes if the other regions' tail after chr13 releases stays
under ~400 s; identical. e2e: assemble ~2,700 + families ~320 + assign ~4,600 = ~7,600 s > 7,200: bar C is
predicted to FAIL on the assemble stage's share unless assign lands near 4,200 s. Stated now so the e2e verdict
cannot be reinterpreted.

## Amendment 1 (2026-11-06 21:00, before any e2e run): bar C identity clause withdrawn
A clean-directory run re-assembles, and assemble is not byte-reproducible (row 1543 Result 7: 7,691 vs 7,988
copies), so "assign tables identical to gw_off" cannot hold for an e2e run. Bar C = time only: each of three
clean runs <= 7,200 s; each run's assign table set is kept and its copy count reported. Identity stays bar B.

## Result gw15 (job 4450704, bin14, RUSTLE_DISPATCH_DEFER=1, staged BAM, singlenode 72 cores)
| bar | value | verdict |
|---|---|---|
| A wall <= 1.1 x chr13 (241.5 + 4,040.2 = 4,281.7 s; bar 4,710 s) | 6,638 s | FAIL |
| B six tables identical to gw_off | all six `same` | PASS |
| B RSS < 250 GB | 215.1 GB | PASS |
The mechanism did what it was built for: chr13 held the pool at log line 745 (t ~ 0), fetched in 241 s and ran
4,040 s. The wall is 2,356 s longer than chr13 because the tail is NOT small: while chr13 held 3.3M of the 4M
read budget (its index ESTIMATE; the actual fetch is 1.42M records), the xchrom families (MCL3/5/8/10 at 500-580k
records each, detect 460-495 s each) and the chromosome hulls chr1/chr2/chr7/chr19 (23k-123k reads) could not
acquire permits and all ran after chr13 released (log lines 195k-266k of 269k). 932 deferrals. Prediction was
wrong on the tail size (~400 s predicted, ~2,350 s observed). Bar C not run (A failed, flag stays opt-in).
Next registered lever (row 1558): release the estimate-minus-actual permits once the fetch count is known.
