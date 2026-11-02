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
