# Pipeline performance: where the time goes and the switches that change it

Living document (started 2026-11-02). Every number here was measured on Fritz (Intel Xeon Platinum 8360Y, 72 cores, `singlenode`)
on A119b (human testis, CHM13) unless stated. Switches marked **identical** were checked to leave every product byte-identical.

## Where the time goes

| stage | A119b wall | what dominates |
|---|---|---|
| assemble (genome-wide, good seeding pool) | ~47 min per arm (32 threads, two arms on one node) | BAM decoding (register 983) |
| families | ~5 min | all-vs-all minimap2 over locus spans |
| assign | **10.75 h** | read-star minimap2 inside copy_assign: 92% of CPU (perf, 2026-11-02) |
| candidates (`candidate_copies`) | days genome-wide (pass A never finished on a 57-family panel) | one indexed BAM query per copy, ~0.85 GB read per query; phase 1 (clusters) over a family pool since r1506 (2.0x on the r1504 panel, identical) |
| assign on a sample without the chr13 hotspot (slim A119b chr17+19, r1503) | 33 min | the BAM window fetch: 89% of wall; a cross-chromosome family's copy windows decode ~1M records twice for ~3k kept reads (r1505, perf: 35% inflate, 14% record conversion) |

The assign stage's time is one chr13 region (5.2-110.7 Mb): 89 families with every copy in chr13:13.9-15.6 Mb hold 86% of the
unresolved rows, from ~30k reads that sit in a median 62 of those families each. Per family the read-star step aligns the
family's molecules against its copies' genomic loci (`-x splice`); on MCL1 (28,788 molecules x 23 copies) it is ~1,450 s at 8
threads, and the region's families ran one after another (`~/asg_prof/`, perf: 67% of CPU in minimap2 chaining `mg_lchain_dp`).

## Switches

| variable / flag | effect | measured | products |
|---|---|---|---|
| `RUSTLE_TIMING=1` | per-region and read-star timing lines in the copy_assign log | — | identical |
| `RUSTLE_STAR_THREADS=N` (default 8) | minimap2 threads of the read-star step | MCL1 read-star 1,452 s (8) / 423 s (32) / 248 s (64); wall 1,754 / 723 / 852 s | **identical** (4 tables) |
| `RUSTLE_FAMILY_THREADS=N` (unset = auto: a region running ALONE gets cores / star threads, capped at its families; otherwise 1; `1` forces sequential) | a region's families in parallel, merged in family order (commit ddcd2e2d; auto 2026-11-02) | 9 chr13 hotspot families: 4,227 s sequential vs 1,664 s auto (2.5x). 57-family panel: auto fans out one region only; walls 555 / 649 / 885 s across reruns are read I/O noise (`reads_in_region` 83 s for 2 reads), so the panel cannot resolve +-20% | **identical** (hotspot 5 tables, panel 6 tables, every rerun) |
| `--region-threads N` / `RUSTLE_REGION_THREADS` | regions in parallel on plain OS threads (not a rayon pool: a worker blocked on a permit while holding others deadlocked, r1543 gw4) | r1543 genome-wide A119b at 8: 6,589 s vs 19,358 s sequential (2.9x) with family 4 x star 16 and budget 4M (gw7), peak 199 GB; these are the pipeline defaults since r1543 (`LOW_MEM=1` keeps 1) | **identical** (6 tables, genome-wide gw6b) |
| `RUSTLE_REGION_READ_BUDGET=N` (binary default 2,500,000, pipeline default 4,000,000; active only with region threads > 1) | a region holds its BAM-index record estimate out of N permits from before its fetch until done; a region above N takes the pool alone | estimate within 0.6-0.85x of real counts; chr13 (3.34M est, 1.42M real) exceeds 2.5M so it runs alone at the tail: 935 regions done at 39 min / 48 GB, then chr13 4,576 s alone rising to 222 GB (its 94 families on 9 family threads are the memory, not the fetch); at 4M with family 4 x star 16 chr13 overlaps the rest: 3,957 s, node peak 199 GB, wall 6,589 s (gw7) | identical |
| `MALLOC_ARENA_MAX=2` | glibc arena cap | RSS flat at 88 GB but ~8x slower (110 regions in 92 min vs 935). **Not a lever.** | — |
| slim BAM (`slimbam`, ratio 0.8) | drops secondaries no step reads | chr20: 12% fewer records; assemble 46 -> 36-44 s | assemble **identical**; assign **changes** (6 of 36 rows lost) |
| mm2-fast (`RUSTLE_MINIMAP2`, `~/tools_home/mm2fast/src/minimap2`, v2.24 + AVX-512) | accelerated chaining | MCL1 at 64 threads: read-star 299 s vs **256 s for our minimap2 2.30 default**; 2.30 with the `--max-chain-skip=1000000` mm2-fast requires takes 3,031 s, so mm2-fast only accelerates a setting we do not use. **Not a lever.** | **changes** (15,946 of ~31,800 assignment lines differ vs 2.30; statuses nearly equal) |

| `RUSTLE_WINDOW_SWEEP=1` (r1507) | one BGZF sweep per chromosome serves every family's copy windows (window store) | slim A119b panel: assign 1,952 -> 391 s (5.0x; MCL0 window fetch 178 -> 2.5 s), 6 tables identical. Genome-wide A119b: fetch 11,285 -> 116 s but OOM-killed (250 GB) on `~xchrom~MCL3` (579k reads); OFF arm 19,358 s | **identical** on the panel; genome-wide not completed: **opt-in** until the store releases finished chromosomes |
| multithreaded BAM reader (r1505) | regions >= 2 Mb only | 1.5x on multi-Mb regions, 3.5x slower on copy windows (seek restarts) | **identical** (5 tables) |
| `candidate_copies` phase 1 over a family pool (r1506) | clusters step parallel across families | clusters 385 -> 189 s, total 846 -> 646 s on the r1504 panel | **identical** (9 products) |
| `candidate_copies --molecules --min-ms-ratio 0.5` (r1513) | nets keep informative secondaries only | 447 vs 416 s on 40 families: nets are primaries-bound, no gain | **changes** (6 vs 3 reference-absent contigs) |

Threads multiply: region threads x family threads x star threads should not exceed the node's cores by much (the panel run at
8 x 8 x 8 was slower than 8 x 1 x 8). Memory genome-wide (250 GB node): eight heavy cross-chromosome fetches at once ~280 GB
(~64 KB per record live); the chr13 hotspot's detect alone climbs ~150 GB.

## Not levers (measured)

Caching alignments by copy sequence (all 590 hotspot copies are distinct); `-G 100k` and query dedup (§6w9: slower / 5.5%
duplicates); Roaring bitmaps and smarter comparisons in assemble (register 983: BAM iteration is the wall); NumPy / bitsets / LRU
caches in copy_assign (our own code is ~5 s of a ~1,450 s family; the rest is minimap2); mm2-fast (slower than minimap2 2.30's
default chaining heuristic on our read-star job, see the table).

## Open levers, by expected gain

1. Done: auto family threads (the hotspot runs alone at the tail and fans out); region threads + read budget + family 4 x star 16 (r1543, 2.9x genome-wide; remaining gap is chr13 starting late: budget held only at minute 29, then 3,957 s of detect; largest-first dispatch (gw8) did not help, 7,071 s: the threads race for the budget; a budget turnstile (gw9) serialised the big regions; scheduling is closed, the next lever is a cheaper read-star for chr13's eight 30k-molecule families).
   Remaining gap to the 5,400 s bar is chr13 serialised after the rest (4,576 s); overlapping it needs its detect peak halved
   (family threads 4 x star threads 16, gw7 pending).
2. One read-star alignment per read per region, reused across families (the hotspot families' padded loci overlap heavily); may
   change products -> needs an identity check, a registration if not identical.
3. `candidate_copies` pass A from the sequential sweep (needed before the candidates lever can be tested at all).
4. One BGZF sweep per chromosome serving every family's light window pass in `copy_assign` (r1505: the 73 cross-chromosome families
   of the slim A119b run spend 1,427 s re-reading the same chr17/chr19 blocks; one multithreaded sweep of that BAM is ~240 s).
   Identical by construction if each window's records are bucketed in coordinate order. The multithreaded reader itself is used
   only on regions >= 2 Mb (1.5x there; 3.5x slower on windows, it restarts on every index seek).
