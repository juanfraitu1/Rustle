# PREREG ASSIGN WINDOW FETCH (register 1505, 2026-11-02)

Why: on the A119b chr17+chr19 slim assign (r1503, `~/msamp/ws/A/x.assign.log`, 64 threads) the `reads_in_region` bracket of
`copy_assign` took 1,737 s of 1,952 s wall (89%); read-star minimap2 216 s, per-read assign 222 s. The 73 cross-chromosome
families (`~xchrom~` keys) cost 1,427 s: `~xchrom~MCL0` fetched 3,079 reads from 86 copy windows in 181.7 s (~2 s per window),
while the plain chr19 region fetched 115,192 reads in 214 s. The cost scales with windows, not reads. Each window opens its own
BAM reader (`BamIndexCache::reads_in_region_mode`, family.rs) and the bracket also holds the readthrough flagging.

Plan: (a) profile one family (`--only-families` MCL0, RUSTLE_TIMING=1, perf or finer timers) to attribute the per-window cost;
(b) remove it (one reader per family, windows fetched in coordinate order in one pass, or whatever (a) names), keeping the record
set and order the window loop produces today.

Decision rule (pre-declared): the change ships as default only if, on the slim A119b run (`~/msamp/ws/A`, same binary flags):
1. every `x.assign.*` product is byte-identical to the r1503 run (assignments, families, quant, famcn_readonly, params minus run
   stamps);
2. assign wall falls by >= 2x on the panel (1,952 s -> <= 976 s).
If 1 fails the change is wrong, not a trade-off: fix or drop. If only 2 fails, record the measured gain and keep it (identical
products cost nothing) but do not claim the hotspot fixed.

## RESULT r1505 (2026-11-02): identical products, bar 2 not met -> kept, hotspot not fixed
`~/asg_prof/r1505`, MCL0 alone (`--only-families`, 16 threads, same flags as r1503). perf (job 4419069): 35% of cycles in
single-threaded BGZF inflate, 14% RecordBuf conversion, 21% minimap2; MCL0's 86 windows hold 870,215 BAM records for 3,079 kept
reads, decoded twice (light pass, then the tied names in full).
| arm | chr17 region (20.7k reads) | chr19 region (119.8k reads) | ~xchrom~MCL0 (86 windows) | wall |
|---|---|---|---|---|
| base | 140.7 s | 142.4 s | 178.2 s | 494 s |
| multithreaded reader on every query | 90.9 s | 98.3 s | 622.2 s | 844 s |
| shipped: multithreaded reader for regions >= 2 Mb, name prefilter on the full pass | 93.8 s | 97.7 s | 166.9 s | 390 s |
All five assign tables byte-identical in both arms. The multithreaded reader restarts its pipeline on every index seek, so a
window (mostly seeks) loses 3.5x while a contiguous region gains 1.5x. Decision rule: bar 1 met, bar 2 (>= 2x) not met on the
window families -> kept as default (identical, free), the window hotspot stands. Next lever (not done): one sweep per chromosome
serving every family's light pass (the 73 cross-chromosome families' windows re-read the same chr17/chr19 blocks; one sweep of the
slim BAM is ~240 s against 1,427 s of window fetches).
