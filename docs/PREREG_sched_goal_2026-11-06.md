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
