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
