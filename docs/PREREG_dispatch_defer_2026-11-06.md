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

## Pre-check result (2026-11-06 14:10): PASS
pre13 (job 4450386; 14 regions, 585 families, budget forced to 50k so almost every region defers): unset 1,471 s,
defer 1,388 s, 475 region records each, 93 deferrals, five tables identical, no hang. gw13 submitted.
Clock note: the "13:00/13:10/13:40" stamps in the gw12 Result and this file's header were estimates; actual
wall-clock was ~12:00-12:25.

## Amendment 1 (2026-11-06 18:10, gw13 running, before gw14): sort key was wrong
gw13 at 30 min: 933 of 936 regions done, 20 deferrals, but chr13's `budget held` is at log line 17,065, after the
rest. Hull span put chr3 (200 Mb hull) and chr10 ahead of chr13 (105 Mb hull, 3.3M records), so chr13 was not
position 0 and starved exactly as in gw12. gw13's wall outcome is decided (~1,800 + 4,124 s, FAIL); it is left
running only as the genome-wide identity check of the deferral code. gw14: bin12 = bin11 with the dispatch key
changed to the BAM-index record estimate (in-memory, no decoding; the quantity the budget is charged in), rule
otherwise unchanged. Bars and prediction as registered: assign <= 5,400 s, identical, RSS < 250 GB; wall
4,500-5,000 s with chr13 held at t ~ 0.

## Result gw14 (2026-11-06 20:05): cancelled at 58 min, outcome decided (FAIL); two loop defects
gw14 (job 4450579): 931 of 936 regions deferred; chr13 `budget held` at log line 8,131 (threads 1-7 raced it for
the pool with try_acquire before it could hold); the main list was exhausted in seconds, and the drain phase then
took deferred regions largest-first with the blocking acquire, so seven threads blocked behind chr13 (597 regions
done at 58 min; gw9's serialisation). Cancelled per the no-waiting rule. gw13 left running as the identity check.

## Amendment 2 (2026-11-06 20:10, before gw15): head gate + rolling retry
1. Head gate: threads other than position 0 wait until position 0 holds its budget (or fails) before their
   first acquire, so chr13 takes the empty pool. Cost: one acquire's latency (~0 s, pool empty).
2. Deferred regions are a rolling queue: a thread pops the front, try_acquire; on failure pushes it back and
   sleeps 200 ms when it has cycled the whole queue. No blocking acquire anywhere after position 0. A region
   that needs the whole pool fits once everything else has finished.
Bars and prediction unchanged (assign <= 5,400 s, identical, RSS < 250 GB; wall 4,500-5,000 s). Pre-check
(pre13 recipe) must pass first: identical tables, all regions computed, no hang.

## Result gw13 (2026-11-06 20:30): identity PASS genome-wide; wall 6,667 s FAIL as decided at 30 min
gw13 (job 4450483, bin11 span key): 6,667 s, RSS 180 GB, six tables identical to gw_off, 20 deferrals, chr13
held at line 17,065 and detect 3,935 s; the deferred regions ran after chr13 (last detect line 17,478), so the
tail added ~440 s over gw12. The deferral path is genome-wide identical; the schedule was the defect (Amendment 1-2).
