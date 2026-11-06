# PREREG: span-first dispatch + budget deferral (row 1548, 2026-11-06 13:40)

## Evidence that reopens scheduling (gw12 log, docs/PREREG_bam_stage_2026-11-06.md)
chr13's `budget estimate` printed at log line 1,076 of 17,438 (dispatched early, flat order) but `budget held`
at line 17,145, after the last other region's line (17,144). Seven threads kept taking small permits, so the
pool never had 3.3M of 4M free; chr13 starved for ~2,000 s, then ran alone 202 + 3,922 s. Wall 6,228 s is the
two phases serialised. Result 4's closing analysis ("they can only overlap by sharing those cores") was
arithmetic; gw12 shows the other 935 regions are fetch-bound (fetch sum 7,799 s of their ~2,100 s x 8 threads).
gw8 (largest-first) failed because the other big regions grabbed the pool first; gw9 (turnstile) because
waiting big regions blocked dispatch threads. Both are addressed here, not retried.

## Rule (opt-in `RUSTLE_DISPATCH_DEFER=1`; unset = current loop, byte-identical and same schedule)
1. Dispatch order within a chunk = hull span (hi - lo) descending, index-stable; output slots keep flat order.
   chr13 (105.5 Mb) is position 0 and takes the empty pool at t = 0. No BAM pre-pass (gw8's 480 s).
2. A region whose budget estimate cannot be taken immediately (`try_acquire`) is deferred, not blocked on;
   the thread moves to the next region. Position 0 always blocks (pool empty, holds at once).
3. After the main list is exhausted, threads drain the deferred list with the blocking acquire.
Identity by construction: the per-region computation and the merge order are unchanged.

## Bars (gw13: singlenode, gw7 recipe, bin11 flag set, staged BAM)
- assign <= 5,400 s; six tables identical to gw_off; RSS < 250 GB.
- Functional pre-check first (<= 15 min): ~40 small regions, budget forced tiny so deferral triggers, flag
  set vs unset: tables identical, every region computed, log shows `deferred` lines.
## Prediction
chr13 holds at t ~ 0; 911 regions <= 0.66M (6,685 s work) run beside it on 7 threads (~1,000 s); the 24
regions > 0.66M (3,046 s work, mostly fetch) are deferred and run after chr13 releases at ~4,150 s; wall
4,500-5,000 s (PASS), chr13 detect 3,900-4,300 s (sharing cores), RSS <= 200 GB, identical.
