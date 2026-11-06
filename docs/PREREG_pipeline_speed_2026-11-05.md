# PREREG PIPELINE SPEED, GENOME-WIDE END TO END (register 1543, 2026-11-05, before any new number)

Goal (user 2026-11-05 /goal): "dramatically speed up all the pipeline so genome-wide analysis can be run end to end quickly".
Baseline (A119b, CHM13, one Fritz node, 72 cores; measured): assemble ~47 min per arm (two layer arms share the node), families
~5 min, assign 19,358 s = 5.4 h (r1507 OFF arm, 2026-11-03, `/anvme/workspace/v138da13-rustle/asgsweep/gw_off`). End to end
~6.3 h. Candidates and flag are opt-in stages that take days genome-wide: OUT of this row's end-to-end definition.
Where the assign wall goes (from that run's timing lines): the pipeline's assign runs regions SEQUENTIALLY (`--region-threads 1`
is the driver default; docs/PERFORMANCE.md wrongly says 8, corrected in this row); window fetch 11,285 s of which 9,839 s is the
740 cross-chromosome family keys (small families decoded window by window, e.g. MCL40: 1,371 reads in 351 s) and 1,445 s plain
regions; detect_and_assign 7,067 s of which the chr13 hotspot region is 4,564 s (89 families; per family read-star minimap2 and
"parallel per-read assign" over ~1M records are equal at ~28k CPU-s each summed).
Targets (wall, same node, same inputs): assign <= 5,400 s (1.5 h; stretch 3,600 s); end to end assemble + families + assign
<= 2.5 h. Identity bars, every lever: the six assign tables (`assignments, families, quant, famcn_readonly, family_join,
copy_conversions` where written) byte-identical to gw_off's; assemble and families products untouched (no change to those
stages unless a later amendment registers one with its own identity bar). Peak RSS recorded (`/usr/bin/time -v`), must stay
under the 250 GB node.
Levers, in order, each developed on the slim A119b panel (`~/msamp/ws/slim`, 391 s) or the chr13 9-family slice (1,664 s) with
the identity cmp, then confirmed genome-wide once:
L1 region threads in the pipeline driver (already cmp-proven identical): default to the node's cores / star threads, memory
   permitting; the xchrom keys' fetch parallelises 8x.
L2 finish r1507's window store for the cross-chromosome keys without the 250 GB peak: serve xchrom windows from a per-chromosome
   sweep that releases records once every key needing that chromosome is served (the ON arm cut fetch to 116 s; the OOM was on
   ~xchrom~MCL3's 579k reads held with everything else).
L3 "parallel per-read assign" over records that produce no row (hotspot families see ~1M records for ~30k molecules): profile on
   the 9-family slice; if a record with no record on any copy of the family has no side effect, exit early. Identical by
   construction only if the six tables cmp equal; otherwise not a lever here.
Out of scope (product-changing, registered rows with KB3781 truth if ever): own-family windows instead of region-union windows;
read-star reuse across families; slim BAM for assign.
Stop rule: a lever whose cmp differs is reverted the same day; the row PASSES when the genome-wide confirmation run meets the
targets with identical tables; otherwise the measured wall is recorded and the row FAILS with the remaining gap named.

## Amendment 1 (2026-11-05, labelled; written after the slim/hotspot identity arms, before the genome-wide confirmation)
L3 WITHDRAWN: the "parallel per-read assign" timer (family.rs t_assign) encloses the read-star timer (t_star), so the two
~28k CPU-s figures are the same cost counted twice, not a second lever. The hotspot floor is read-star minimap2 CPU.
Added lever L4, scheduler: StarSlots (family.rs) gives every family a permit from a pool of cores / RUSTLE_STAR_THREADS (8),
so families fan out even when several regions run (the ae73fc57 auto rule set family threads to 1 whenever ACTIVE_REGIONS>1,
which would have serialised the chr13 hotspot under L1). RUSTLE_FAMILY_SLOTS=0 restores the old rule.
L2 as built: RUSTLE_WINDOW_SWEEP=1 builds the per-chromosome sweep store for RUSTLE_SWEEP_BATCH keys at a time (default 64)
and drops it after the batch; the whole-sweep store of r1507 is gone.
Identity results so far (frozen binaries ~/r1543/bin): slim panel batch 4, region threads 1: six tables same, 1,576 s, RSS 7.4 GB.
Hotspot 9 families, slots on: five tables same, 1,626 s (ref 1,664 s), RSS 91 GB.

## Amendment 2 (2026-11-05, labelled; written after the killed genome-wide run, before the relaunch)
Genome-wide arm (sweep on, batch 64, 8 region threads, slots on; job 4442328): OOM-killed (exit 137) at 239.8 GB RSS after
573 of 936 regions, 10,095 s at the kill. The six DIFFER lines are the incomplete tables of a killed run, not a logic
difference; the stop rule's revert applies to L2 only. L1 and L4 identity still rest on slim b4r8 and hot9 (six/five tables same).
L2 WITHDRAWN on two measured grounds: (a) the nine per-batch sweeps are serialised before each batch's compute and cost
4,600 s of the 10,095 s (slim: 38 sweeps 2,616 s on both arms against a 391 s baseline); the direct xchrom fetch they replace
(9,839 s sequential) parallelises under L1 instead. (b) `flat.chunks(batch)` makes batch N+1 wait for batch N, so the chr13
hotspot (4,523 s, same as alone) held the other seven region threads idle. The code stays gated behind RUSTLE_WINDOW_SWEEP
(default off); it is closed, not adopted. Saved log: /anvme/workspace/v138da13-rustle/r1543/gw_killed/.
Relaunch arm (registered before its numbers): sweep OFF, RUSTLE_REGION_THREADS=8, slots on (RUSTLE_FAMILY_SLOTS=1), same
six-table identity bars against gw_off, wall <= 5,400 s, RSS < 250 GB. Prediction: 4,600-5,200 s, chr13 the floor, provided
the eight direct xchrom fetches plus chr13's 91 GB fit (no prior measurement of that combination). Same environment run on the
slim panel as the identity arm for this exact combination (b4r8 had the sweep on; hot9 ran one region).

## Amendment 3 (2026-11-05, labelled; after the second kill, before the third genome-wide run)
Relaunch arm (job 4442736, sweep off, 8 region threads, slots): OOM-killed at 241.4 GB, 3,433 s, 911 of 936 regions done,
the chr13 hotspot still running. Slim under the same environment: six tables same, 526 s, 8.9 GB (identity of L1+L4 holds).
Cause, from gw_off's reads_in_region counts: the eight heaviest regions (chr13 1.42M records; ~xchrom~MCL3/5/8/10 0.5-0.58M;
MCL18/44/70 0.19-0.29M) sum to 4.4M records; chr13 alone peaks at 91 GB (~64 KB per record), so eight at once need ~280 GB.
Neither window span (chr13 is 2.6 Mb) nor the copies table's n_reads predicts this before the fetch.
Lever L5, read budget: after the fetch (cheap: the records Vec), a region acquires its record count from a pool of
RUSTLE_REGION_READ_BUDGET records (default 2,500,000 ~ 160 GB) before detect_and_assign and releases it after; a region above
the budget takes the whole pool and runs alone. Active only with region threads > 1; star slots and the budget share one
`Permits` type (family.rs). Identity by construction (ordering only); still checked on the slim panel in the same job.
Third genome-wide arm: sweep off, 8 region threads, slots, budget 2.5M; bars unchanged (six tables, <= 5,400 s, < 250 GB).
Prediction: 4,800-5,600 s (chr13 floor plus some queueing of the four 0.5M-record keys behind it), RSS 160-200 GB.

## Amendment 4 (2026-11-05, labelled; after the third kill, before the fourth genome-wide run)
Third arm (job 4442751, budget on the fetched count before detect_and_assign): OOM-killed at 242.6 GB, 2,345 s. Slim same
environment: six tables same, 597 s, 9.9 GB. At the kill 912 regions were fetched and 911 finished: only chr13 was in
detect_and_assign; the other seven threads were mid-FETCH of the heavy cross-chromosome keys (MCL3/5/8/10/18/44/80, none
had logged reads_in_region). The memory is in the fetch itself, so a budget taken after the fetch cannot bound it.
L5 as rebuilt: the budget is acquired BEFORE the fetch from the BAM index alone (`BamIndexCache::estimate_records`: the
windows' index-chunk compressed bytes / the BAM's mean compressed bytes per record, the `deep_families_estimate` method of
the max-window-records guard; seconds, no decoding, good to a small factor) and held until the region is done. Default
2,500,000 estimated records. The log's `budget estimate` lines are compared with the real `reads_in_region` counts in the
result. Bars unchanged. Prediction: 4,800-5,800 s, RSS < 200 GB.

## Amendment 5 (2026-11-05 16:45, before gw5): gw4 deadlock, region threads moved off rayon

gw4 (job 4442782, pre-fetch index budget) did not OOM: the assign log stops at 08:56 after 334/936 regions,
last line `budget estimate ~xchrom~MCL378: 682239 records`, and the job idled to its 8 h wall. Cause: regions
ran as a rayon `par_iter` on a pool shared with the inner per-copy work; a worker holding budget permits
enters an inner join, steals another region's job, and blocks in `Permits::acquire` on permits it holds
itself. Fix: regions now run on `--region-threads` plain OS threads over a shared index (`std::thread::scope`);
inner `par_iter`s use the global pool bounded by `--threads`. Serial path unchanged. No scoring code touched.
Relaunch gw5 + slim5 with the same env (sweep 0, region threads 8, family slots 1, budget 2.5M).
Prediction unchanged from Amendment 4: wall 4,800–5,800 s, RSS < 200 GB, six tables identical.

## Amendment 6 (2026-11-05 19:05, before gw6): gw5 OOM with the budget enforced; allocator A/B

gw5 (job 4446337, OS-thread regions, budget 2.5M): no deadlock; OOM-killed at 242.6 GB after 6,609 s with
755/936 regions done. The `budget estimate` line prints BEFORE the permit wait, so the eight regions with an
estimate and no completion were mostly queued, not fetching: at most ~2.05M estimated records (MCL3 fetched,
579k real, in detect; MCL61/62/621/623 could fit) were held, ~130 GB at 64 KB/record. gw3 reached the same
242.6 GB with only chr13 live. So >100 GB sits OUTSIDE the budgeted regions, and lowering the budget cannot
fix it. Wall: 6,609 s at 755/936 with chr19 (3.2M est), MCL53 (2.57M), MCL4 (1.8M) still queued, so this
configuration misses the 5,400 s bar regardless; chr13 holds 90% of the pool for 4,577 s and the heavy
regions serialise behind it.
Hypotheses: H1 glibc per-thread arena retention (no custom allocator; freed small allocations from one
thread's arena are not reusable by others and not trimmed); H2 retained `RegionWork` output (estimated a few
GB, unlikely). Arms, same binary and env as gw5 plus `[timing] budget held` and a 30 s RSS series:
gw6a exports MALLOC_ARENA_MAX=2; gw6b does not. Predictions: H1 -> gw6a peak < 200 GB and finishes, gw6b
ratchets monotonically to the kill; H2 -> both ratchet and RSS does not step down when chr13 finishes.
gw6a may run slower from malloc lock contention. Wall bar is NOT expected to pass in either arm; the arm
result decides whether the budget can be raised (memory not in records) or the row FAILS on wall.

## Amendment 7 (2026-11-05 21:10, before gw7): chr13 detect is the memory, not retention

gw6b (default glibc): chr13's index estimate is 3.34M > the 2.5M cap, so it ran ALONE after all 935 other
regions (done at minute 39, 48 GB). RSS then rose from 48 GB to > 200 GB during chr13's detect alone: 94
families on 9 family threads, the eight largest ~30k molecules x 16-31 copies, read-star 900-1,400 s each.
The serial baseline ran chr13 the same way on a clean heap and fit; after 48 GB of retained memory it does
not. gw6a (MALLOC_ARENA_MAX=2): RSS flat at 88 GB but 110 regions in 92 min (~8x slower); cancelled, not a
usable switch. H1 (allocator retention) is a ~48 GB contributor, not the cause.
gw7 (launched before gw6b finished, to save a wait): same binary, RUSTLE_FAMILY_THREADS=4 RUSTLE_STAR_THREADS=16
(documented identity-preserving switches, docs/PERFORMANCE.md), RUSTLE_REGION_READ_BUDGET=4000000 so chr13
(3.34M) no longer takes the whole pool and other regions overlap it. Predictions: chr13 peak ~100 GB, node
peak 180-220 GB, wall 4,800-5,400 s, six tables identical. If OOM: family threads 3 and budget 3.5M as
Amendment 8. If wall > 5,400 s with identity: FAIL, gap named.

## Result 1 (2026-11-05 21:50): gw6b completed, genome-wide identity holds

gw6b (job 4446733; sweep off, region threads 8, family slots 1, budget 2.5M, default glibc): assign exit 0,
wall 8,697 s (baseline 19,358 s; 2.2x), peak RSS 222.2 GB, six tables byte-identical to gw_off. Timeline:
935 regions done by minute 39 (48 GB), then chr13 alone 4,576 s (its estimate 3.34M exceeds the 2.5M cap, so
it takes the whole pool) with RSS rising to 222 GB. Wall bar (<= 5,400 s) NOT met by this configuration;
identity and RSS bars met. L1 (region threads) + L4 (permits) + L5 (budget) are confirmed genome-wide.
Remaining gap to the bar: 3,300 s, all of it chr13 running serially after the rest. gw7 (Amendment 7) is
the arm that overlaps chr13 with the other regions.

## Result 2 (2026-11-05 23:35): gw7 completed; wall bar FAILS by 1,189 s, identity and RSS hold

gw7 (job 4447078; region threads 8, family slots 1, family threads 4, star threads 16, budget 4M): assign
exit 0, wall 6,589 s (baseline 19,358 s; 2.9x; gw6b 8,697 s), peak RSS 198.7 GB (gw6b 222 GB), six tables
byte-identical to gw_off. chr13 overlapped the other regions (budget held at minute 29) and its detect took
3,957 s (gw6b 4,576 s alone). Amendment 7 predictions: node peak 180-220 GB met; wall 4,800-5,400 s NOT
met. Verdict per the pre-committed rule: assign bar (<= 5,400 s) FAIL, gap 1,189 s = the chr13 region's own
detect (eight ~30k-molecule x 16-31-copy families, read-star 900-1,400 s each at 16 star threads); nothing
outside that region can shorten the run further. Identity and RSS bars PASS. Shipped: gw7 settings are the
driver defaults (tools/rustle_pipeline.sh: region threads 8, family threads 4, star threads 16, budget 4M;
LOW_MEM=1 keeps regions sequential), since they beat gw6b on every bar with identical outputs. Amendment 8
(family threads 3, budget 3.5M) not run: memory passed, so it would only slow chr13. Stretch 3,600 s and
end-to-end 2.5 h not reached; closing the chr13 gap needs a cheaper read-star for 30k-molecule families, a
new row.

## Amendment 9 (2026-11-05 23:55, after Result 2, before gw8): the gap is chr13's late start

Correction to Result 2: "nothing outside that region can shorten the run" was wrong. gw7 timeline: chr13's
budget was not held until minute 29 (the pool had to drain below 0.66M while small regions ran), its fetch
was then contended (868 s), and its detect (3,957 s) finished at 6,560 s. chr13 idled ~1,700 s. gw8: same
binary logic and env as gw7 plus largest-estimate-first dispatch (regions sorted by BAM-index estimate
descending before the OS-thread loop; output order unchanged). Predictions: chr13 starts at t=0 with the
whole pool, wall ~ fetch + detect = 4,500-4,900 s (PASS); peak RSS <= 220 GB (gw7 overlapped chr13 with the
rest at 199 GB; more overlap early); six tables identical. End to end then ~ 47 + 5 + 80 min ~ 2.2 h, under
the 2.5 h bar; at gw7's 6,589 s end to end is ~2.7 h (FAIL by arithmetic). If gw8 OOMs: row stays FAIL on
wall at gw7's figure; if wall > 5,400 s: FAIL, gap named. L3 (per-read early exit) was never attempted.

## Result 3 (2026-11-06 01:45): gw8 largest-first dispatch did not help; row closes FAIL on wall at gw7's figure

gw8 (job 4447200): exit 0, wall 7,071 s (gw7 6,589 s), peak RSS 194.3 GB, six tables identical. chr13's
detect was 3,962 s (gw7 3,957 s), but its budget was held only after 236 regions had finished: the eight
threads race for the pool, the seven other big regions (1.1-3.2M each) took it first, and chr13 (3.3M of 4M)
had to wait for the pool to drain again. Amendment 9's prediction (4,500-4,900 s) FAILED. Ordering the
dispatch list is not enough; chr13 needs to hold the pool before any other region acquires. Not run tonight.
Final verdict for row 1543: assign wall bar FAIL (best 6,589 s, gw7, 2.9x), identity PASS, RSS PASS, end to
end ~2.7 h FAIL by arithmetic. Shipped defaults = gw7 settings; largest-first sort stays (harmless, identical).
Next row: budget acquisition in dispatch order (a turnstile so the largest region holds first), predicted
wall ~ 870 + 3,960 = 4,800-5,000 s.

## Amendment 10 (2026-11-06 02:00, after Result 3, before gw9): budget acquired in dispatch order

gw9: gw8 binary plus a turnstile: a region waits for its dispatch position before acquiring the read budget
and releases the turn once it holds it, so the largest region (chr13, position 0) holds the pool before any
other region can take permits; everything else unchanged (env as gw7). Predictions: chr13 budget held at
t ~ 0, wall ~ 870 + 3,960 = 4,800-5,100 s (PASS if <= 5,400), peak RSS <= 220 GB, six tables identical.
If wall > 5,400 s: FAIL stands at gw7's figure.

## Result 4 (2026-11-06 02:40): gw9 turnstile cancelled at 30 min; scheduling cannot meet the bar

gw9 (job 4447403): chr13 held the pool at t = 0 as predicted, but only 1 other region finished in 30 min
(gw7: 467). The turnstile makes position k+1 wait until position k HOLDS its permits; positions 1-7 are
1.1-3.2M-record regions that do not fit in the 0.66M left beside chr13, so every region queued behind them.
Projected wall >= 4,830 (chr13) + 2,300 (the rest at 8 threads) > 7,000 s: outcome decided, cancelled per
the no-waiting rule. Turnstile code reverted (largest-first sort kept). Amendment 10 prediction FAILED.
Closing analysis: chr13 alone needs fetch + detect ~ 4,800 s on 64 cores (family 4 x star 16); the other 935
regions need ~39 min of 8-thread wall. They can only overlap by sharing those cores, which is what gw7 did
(6,589 s). Meeting 5,400 s needs chr13's own detect cut (cheaper read-star for the eight 30k-molecule
families), not a different schedule. Row 1543 verdict unchanged: FAIL on wall at 6,589 s (2.9x), PASS
identity + RSS, end to end ~2.7 h FAIL. gw7 settings remain the shipped defaults.

## Amendment 11 (2026-11-06 02:55, before gw10): more cores, since the floor is CPU-bound

chr13's eight heavy families cost ~153k core-seconds (eight x ~1,200 s x 16 threads) and scale near
linearly with star threads (MCL1 1,452 s at 8 vs 248 s at 64, docs/PERFORMANCE.md). gw7 at 72 cores is near
the CPU floor, so the remaining lever is a bigger node. gw10: spr1tb (104 cores, 1 TB), gw8 binary (sort kept,
no turnstile), env as gw7 except RUSTLE_STAR_THREADS=24 (family 4 x star 24 = 96), --threads 96,
RUSTLE_REGION_READ_BUDGET=6M. Predictions: wall 4,600-5,200 s (PASS if <= 5,400), peak RSS 220-300 GB (fine on
1 TB; the 250 GB bar is a singlenode bar and does not apply to this node class, stated now), six tables
identical. If wall > 5,400 s: FAIL stands; the bar is then out of reach without changing read-star itself.

## Result 5 (2026-11-06 04:25): gw10 on 104 cores, 5,507 s; wall bar missed by 107 s

gw10 (job 4447416, spr1tb 104 cores, family 4 x star 24, budget 6M, largest-first sort): exit 0, wall 5,507 s,
peak RSS 203.6 GB, six tables identical. chr13 held the pool at t = 0 and its detect took 3,011 s (gw7
3,957 s), but only 11 other regions finished in the first 64 min: the sort put the seven other big regions
(1.1-3.2M records) on the free threads with contended fetches, so the 925 small regions ran late. Amendment 11
prediction (4,600-5,200 s) FAILED by 307 s; the bar by 107 s. End to end by arithmetic 47 + 5 + 92 min = 2.4 h,
under 2.5 h on this node class only. Verdict for row 1543 stays FAIL on the assign bar (best 5,507 s on 104
cores, 6,589 s on the 250 GB singlenode class the bar was written for), PASS identity and RSS. Scheduling
and cores are exhausted; the next row is a per-family pilot with early exit: the eight chr13 hotspot families
(all de novo copies at chr13 13.9-15.1 Mb, no CAT gene) align ~200k reads for ~120 unique assignments.

## Amendment 12 (2026-11-06 04:40, before e2e): time the pipeline end to end with the shipped defaults

The end-to-end bar (<= 2.5 h) was never measured. e2e job: `rustle_pipeline.sh all --to assign` from a clean
directory on singlenode (72 cores, 250 GB), shipped defaults, no RUSTLE_* speed variables exported, --no-cache.
The flag stage is excluded (--to assign), as in the baseline arithmetic. Prediction: 47 + 5 + 110 min ~ 2.7 h
(FAIL); the six assign tables identical to gw_off. This measures, it does not change the verdict.

## Correction (2026-11-06 05:30, before gw11 and e2e rerun): gw10 and the first e2e ran the TURNSTILE binary

`~/r1543/bin8/copy_assign` was rebuilt for gw9 (turnstile) at 01:44 and never rebuilt after the source
revert; gw10 and e2e job 4447459 used it. So Result 5 measured turnstile + largest-first on 104 cores, not
"gw8 binary, sort kept, no turnstile" as Amendment 11 states (gw10's log fits: chr13 held at 0 regions,
11 regions in 64 min). Result 5's wall (5,507 s) and identity stand as measured; its attribution is wrong
and a clean 104-core arm is untested. Also: the largest-first sort is NOT harmless (Result 3 claimed so):
gw8 vs gw7 on the same node and env differ only by the sort, 7,071 vs 6,589 s. Sort reverted in source;
the shipped source now equals the gw7 logic. New binary `bin9` built from HEAD.
gw11: spr1tb, gw7 env with RUSTLE_STAR_THREADS=24, budget 6M, --threads 96, bin9 (no ordering trick).
Prediction: with 6M, chr13 (3.3M) fits beside the others from the start; wall 3,900-4,800 s (PASS),
RSS above gw10's 204 GB (fine on 1 TB), six tables identical.
e2e rerun with bin9 (Amendment 12 otherwise unchanged); e2e job 4447459 cancelled, no numbers taken from it.

## Result 6 (2026-11-06 07:55): gw11, clean 104-core arm, 4,514 s; assign bar PASSES on that node class

gw11 (job 4447489, spr1tb 104 cores, bin9 = shipped gw7 logic, family 4 x star 24, --threads 96, budget 6M):
exit 0, wall 4,514 s, peak RSS 217.1 GB, six tables identical to gw_off. chr13 detect 3,269 s (gw7 3,957 s
at star 16). Correction prediction (3,900-4,800 s, RSS > 204 GB, identical) held. Against the registered bars:
assign <= 5,400 s PASS on the 104-core class (stretch 3,600 missed); the 250 GB bar was written for singlenode
and does not apply here (Amendment 11). On the 250 GB singlenode class the best figure stays gw7 6,589 s, FAIL.
gw10 (5,507 s) vs gw11 (4,514 s) on the same node and env isolates turnstile + sort at ~1,000 s, confirming
the Correction. Settings for the 104-core class are not shipped as defaults (they assume 1 TB and 96 cores);
`RUSTLE_STAR_THREADS=24 RUSTLE_REGION_READ_BUDGET=6000000 --threads 96` is the recorded recipe.

## Result 7 (2026-11-06 08:05): e2e measured 8,905 s = 2.47 h on singlenode; end-to-end bar PASSES by 95 s

e2e (job 4447488, singlenode 72 cores / 250 GB, `all --to assign`, clean directory, shipped defaults, bin9):
exit 0, wall 8,905 s, peak RSS 214.0 GB. Stages: assemble 2,700 s, families 324 s, assign 5,865 s. Bar
<= 2.5 h (9,000 s) PASS with 95 s to spare; Amendment 12 predicted 2.7 h FAIL and was wrong because assign
ran 5,865 s, not 6,589 s. The six assign tables are NOT identical to gw_off and the comparison is not
interpretable: the fresh assembly differs from gw_off's upstream (a1.gtf 2,193,334 vs 2,178,364 lines;
a1.fam.copies.tsv 7,691 vs 7,988 copies; a1.families.tsv identical), so assign saw a different copy set, which
also explains the shorter assign wall. Assign identity on a fixed input is established by gw7 and gw11, not by
this run. Where the upstream difference comes from (gw_off's assembly predates the 2026-11-02 layers split and
the r1525 strand vote) is not part of this row. Margin is thin: on the same class gw7 assign was 6,589 s, which
would put the end to end at 2.67 h FAIL; the bar passes on this measurement but not with margin.

## Correction to Result 7 (2026-11-06 08:20, after advisor review; appended, Result 7 text left as written)

Two claims in Result 7 do not hold. (a) "the different copy set also explains the shorter assign wall" is
wrong: e2e chr13 detect was 3,913 s (gw7 3,957 s) and the chr13 13.9-15.1 Mb hotspot has 438 copies in e2e vs
435 in gw_off. The env is the same as gw7 (sweep off and slots on are the binary defaults). The 724 s gap
between e2e assign (5,865 s) and gw7 (6,589 s) is therefore run-to-run variation in when chr13 obtained its
budget, n = 1 each, not a property of the input. e2e's second and third regions were slower (1,316 / 1,025 s
vs 716 / 709 s), consistent with different overlap, not with less work. (b) "PASS by 95 s" overclaims: one run
at a 1% margin on a shared node. The honest verdict line is "8,905 s measured, at the bar, n = 1". With gw7's
assign the same arithmetic gives 2.67 h. (c) Amendment 12's identity prediction was unsatisfiable by
construction: a clean directory regenerates the copy set under current defaults (layers split, r1525 strand
vote shipped after gw_off's upstream), so e2e was never a valid identity test; that is a design error in the
amendment, not a result. Identity rests on gw7 and gw11, which ran on gw_off's own copy set.

## Goal state (2026-11-06 08:20)

Assign: 2.9x with byte-identical outputs on the 250 GB class (6,589 s, misses the 5,400 s bar); 4,514 s on
104 cores (bar met on that class only). End to end from a clean directory: 8,905 s (~2.5 h) vs ~6.4 h baseline
arithmetic, at the bar with no margin. Further speed-ups (per-family pilot + early exit, dropping the genomic
read-star form where units explain the reads, per-bubble graph vote) change outputs and need the user's call.

## Goal state update (2026-11-06 13:10)
Identity-preserving levers are exhausted: r1545 pilot early exit (changes outputs, no time), r1546 unit form
(different answer), r1547 BAM staged on NVMe (gw12 6,228 s identical, FAIL; fetch -31% but chr13 detect 3,922 s
is the critical path). 250 GB class stands at 6,228-6,589 s identical; bar met only on 104 cores (gw11 4,514 s);
e2e at the bar n = 1. Next step is the user's: accept the 104-core recipe, or authorise an output-changing lever
(yield-rate gate, per-bubble graph vote) under a new PREREG.
