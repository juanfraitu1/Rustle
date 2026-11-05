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
