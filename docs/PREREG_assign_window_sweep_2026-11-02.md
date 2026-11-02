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
