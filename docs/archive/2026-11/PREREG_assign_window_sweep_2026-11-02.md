# PREREG ASSIGN WINDOW SWEEP (register 1507, 2026-11-02)

Why: r1505 left the window hotspot standing: on the slim A119b assign the 73 cross-chromosome families spend 1,427 s of 1,952 s
re-reading the same chr17/chr19 BGZF blocks window by window (MCL0: 870k records decoded twice for 3,079 kept reads), while one
multithreaded sweep of that BAM is ~240 s.

Change (`copy_assign`, the `--tied-prefilter` path with an indexed BAM, which the pipeline always takes): before the region sweep,
(1) one indexed query per chromosome over the hull of every region's copy windows, with the 4-worker BGZF reader, bucketing each
record into every window it overlaps (noodles' own overlap rule: 1-based start <= hi and end >= lo + 1) as the light pass does
today; (2) per region, the tied names from those light records (the same dedup and `as_tied` rule the window loop applies);
(3) a second per-chromosome query keeping the full records whose name is tied in the window's region. The window loop then takes
its light and full results from the store instead of querying. Records reach a window in file (coordinate) order, as a query
returns them, so every product is the same by construction. Unset / no index / `--keep-coordinate-duplicates` / streaming:
unchanged path. Opt-out `RUSTLE_WINDOW_SWEEP=0`.

Decision rule (pre-declared), slim A119b run (`~/msamp/ws/A`, flags of r1503, 64 threads):
1. all five `x.assign.*` tables byte-identical to r1503 (`assignments`, `families`, `quant`, `famcn_readonly`, `family_join`);
2. assign wall <= 976 s (>= 2x against 1,952 s).
Then a genome-wide identity check on A119b assign is required before the default flips (its wall is reported, no bar: the chr13
hotspot is read-star bound). If 1 fails anywhere: fix or drop. If 2 fails: keep opt-in (`RUSTLE_WINDOW_SWEEP=1`), record the gain.

## RESULT r1507, panel (2026-11-03): identical, 5.0x -> opt-in `RUSTLE_WINDOW_SWEEP=1` pending the genome-wide check
Job 4420759, slim A119b chr17+chr19 through `rustle_pipeline.sh assign` (64 threads, binary of commit below):
assign wall 1,952 s (r1503) -> **391 s**; window sweep 786 windows of 76 regions over 2 chromosomes: light 72.1 s, tie + full
46.3 s; the window loop's own `reads_in_region` brackets sum to 29.6 s (MCL0 178 -> 2.5 s). All six assign tables
(`assignments`, `families`, `quant`, `famcn_readonly`, `family_join`, `copy_conversions`) byte-identical to r1503; `params.tsv`
differs only in the two input paths. Bars 1 and 2 met. Genome-wide A119b identity (sweep on vs off, same binary) pending.

## RESULT genome-wide (2026-11-03, jobs 4420941 on / 4420942 off, slim A119b, singlenode 250 GB, no memory cap)
Off: exit 0, 19,358 s, reads_in_region sum 11,285 s. On: OOM-KILLED at 7,359 s (exit 137, one oom_kill event) while
assigning `~xchrom~MCL3` (579,016 reads, 220 copies); its window sweep itself was done (4,427 windows of 936 regions over
25 chromosomes: light 510 s, tie + full 357 s; reads_in_region sum 116 s = 97x less fetch time). No assign table written,
so the identity bar is NOT met: `RUSTLE_WINDOW_SWEEP` stays opt-in. Cause to fix before a rerun: the window store keeps
every chromosome's fetched records alive while the largest family runs; release a chromosome's store once its families are
done (or cap the store), then rerun with `RUSTLE_TIMING=1` and the peak-RSS line.
