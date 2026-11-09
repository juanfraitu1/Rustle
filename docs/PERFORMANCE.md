# Pipeline performance: where the time goes and the switches that change it

Living document (started 2026-11-02). Every number here was measured on Fritz (Intel Xeon Platinum 8360Y, 72 cores, `singlenode`)
on A119b (human testis, CHM13) unless stated. Switches marked **identical** were checked to leave every product byte-identical.

## Where the time goes

Current (`all` on A119b, 72 threads, `audit/e2e0`, 2026-11-09; censat exclusion and SD damping on):

| stage | A119b wall | what dominates |
|---|---|---|
| assemble | 595 s | BAM decoding; the sweep is I/O-bound on /anvme (~610 MB/s one stream) |
| families | 465 s -> **385 s** (audit 22: minimap2 `-K 4G`, one query batch, identical) | all-vs-all minimap2 over locus spans: the index build (~90 s on ~3 threads) and the last batch's tail leave the node idle |
| assign | 549 s | the last region, `~xchrom~MCL13` (150 copies): read-star 123 s of its 150 s |
| index | 124 s, once per genome | minimap2 `-d`; shared next to the FASTA |
| flag scan | 283-314 s, peak 16.6 GB | BAM sweeps. Audit rounds 17-20: before them the old scan ran over 5 h on chr13 alone and hit the job limit |
| flag align + verdict | 41 s, peak 41 GB | the 12.3 GB splice index in memory |
| **`all` end to end, from scratch** | **1,328 s (22 min)**, `audit/e2e2` (1,572 s in e2e1, before the background flag scan of audit 39). The scan runs under families and assign only when assemble runs in the same call: a resumed run that skips assemble, or `--low-mem`, scans after assign (about +290 s) | assemble 397 (as_table cached; +~190 s once per new BAM), families 372, assign 475, index 0 (shared), flag 327. Every product identical to e2e0. Gorilla KB3781 with no extra env: 756 s |
| `all --low-mem` (gorilla KB3781) | 1,549 s, node peak 43.7 GB in families (audit 40) | the split index changes 7 of 16,467 flag rows, no verdict class; every other product identical. Families' minimap2 memory scales with threads (gorilla 38 GB at 72, 18 GB at 16, 14 GB at 8; `-K` 4G vs 500M changes it by 2 GB; PAF identical in all, audit 43), so since audit 43 `--low-mem` runs families on at most 16 threads (127 -> 219 s on gorilla) |

History (2026-11-02, before the region threads, the censat exclusion and the window fixes):

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
| `--region-threads N` / `RUSTLE_REGION_THREADS` | regions in parallel on plain OS threads (not a rayon pool: a worker blocked on a permit while holding others deadlocked, r1543 gw4) | r1543 genome-wide A119b at 8: 6,589 s vs 19,358 s sequential (2.9x) with family 4 x star 16 and budget 4M (gw7), peak 199 GB; these are the pipeline defaults since r1543 (`LOW_MEM=1` keeps 1); on 104 cores / 1 TB with star 24, budget 6M, `--threads 96`: 4,514 s, 217 GB (gw11); BAM staged on NVMe instead of NFS home (r1547 gw12): 6,228 s, fetch sum 7,799 vs 11,308 s, chr13 fetch stable at ~200 s instead of 189-868 s, wall gain only 361 s because chr13 detect (3,922 s) is the critical path | **identical** (6 tables, genome-wide gw6b) |
| `RUSTLE_REGION_READ_BUDGET=N` (binary default 2,500,000, pipeline default 4,000,000; active only with region threads > 1) | a region holds its BAM-index record estimate out of N permits from before its fetch until done; a region above N takes the pool alone | estimate within 0.6-0.85x of real counts; chr13 (3.34M est, 1.42M real) exceeds 2.5M so it runs alone at the tail: 935 regions done at 39 min / 48 GB, then chr13 4,576 s alone rising to 222 GB (its 94 families on 9 family threads are the memory, not the fetch); at 4M with family 4 x star 16 chr13 overlaps the rest: 3,957 s, node peak 199 GB, wall 6,589 s (gw7) | identical |
| `MALLOC_ARENA_MAX=2` | glibc arena cap | RSS flat at 88 GB but ~8x slower (110 regions in 92 min vs 935). **Not a lever.** | — |
| slim BAM (`slimbam`, ratio 0.8) | drops secondaries no step reads | chr20: 12% fewer records; assemble 46 -> 36-44 s | assemble **identical**; assign **changes** (6 of 36 rows lost) |
| mm2-fast (`RUSTLE_MINIMAP2`, `~/tools_home/mm2fast/src/minimap2`, v2.24 + AVX-512) | accelerated chaining | MCL1 at 64 threads: read-star 299 s vs **256 s for our minimap2 2.30 default**; 2.30 with the `--max-chain-skip=1000000` mm2-fast requires takes 3,031 s, so mm2-fast only accelerates a setting we do not use. **Not a lever.** | **changes** (15,946 of ~31,800 assignment lines differ vs 2.30; statuses nearly equal) |

| `RUSTLE_WINDOW_SWEEP=1` (r1507) | one BGZF sweep per chromosome serves every family's copy windows (window store) | slim A119b panel: assign 1,952 -> 391 s (5.0x; MCL0 window fetch 178 -> 2.5 s), 6 tables identical. Genome-wide A119b: fetch 11,285 -> 116 s but OOM-killed (250 GB) on `~xchrom~MCL3` (579k reads); OFF arm 19,358 s | **identical** on the panel. Audit round 12 (2026-11-09, current defaults, genome-wide A119b): no OOM (38.8 GB) and 7/7 tables identical, but **2,914 s vs 546 s** for the default fetch. **Opt-in**, not a lever at current defaults |
| multithreaded BAM reader (r1505) | regions >= 2 Mb only | 1.5x on multi-Mb regions, 3.5x slower on copy windows (seek restarts) | **identical** (5 tables) |
| `candidate_copies` phase 1 over a family pool (r1506) | clusters step parallel across families | clusters 385 -> 189 s, total 846 -> 646 s on the r1504 panel | **identical** (9 products) |
| `candidate_copies --molecules --min-ms-ratio 0.5` (r1513) | nets keep informative secondaries only | 447 vs 416 s on 40 families: nets are primaries-bound, no gain | **changes** (6 vs 3 reference-absent contigs) |
| flag scan, pass 2 one query over the fired loci (audit 17-19) | chr13 > 5 h -> 370 s; A119b chr20-22+X 823 -> 225 s | one BAM query per mixture-firing locus re-read the nested chr13 hotspot (4.4 TB read); per-span queries re-read the BAI coarse bins (84 GB for chr13's 9 GB) | **identical** (gorilla 8 contigs, A119b chr20-22+X, chr13 hotspot window) |
| flag scan, 30 Mb locus groups in parallel (audit 20) | whole-genome A119b 509 -> 283 s, peak 16.6 GB | one thread per contig left chr13 decoding alone; 10 Mb groups read 511 GB for 294 s | **identical** |

Threads multiply: region threads x family threads x star threads should not exceed the node's cores by much (the panel run at
8 x 8 x 8 was slower than 8 x 1 x 8). Memory genome-wide (250 GB node): eight heavy cross-chromosome fetches at once ~280 GB
(~64 KB per record live); the chr13 hotspot's detect alone climbs ~150 GB.

## Not levers (measured)

Caching alignments by copy sequence (all 590 hotspot copies are distinct); `-G 100k` and query dedup (§6w9: slower / 5.5%
duplicates); Roaring bitmaps and smarter comparisons in assemble (register 983: BAM iteration is the wall); NumPy / bitsets / LRU
caches in copy_assign (our own code is ~5 s of a ~1,450 s family; the rest is minimap2); mm2-fast (slower than minimap2 2.30's
default chaining heuristic on our read-star job, see the table).

## Open levers, by expected gain

1. Done: auto family threads (the hotspot runs alone at the tail and fans out); region threads + read budget + family 4 x star 16 (r1543, 2.9x genome-wide; remaining gap is chr13 starting late: budget held only at minute 29, then 3,957 s of detect; largest-first dispatch (gw8) did not help, 7,071 s: the threads race for the budget; a budget turnstile (gw9) serialised the big regions; scheduling is closed, 104 cores (gw10, star 24) 5,507 s identical; the next lever is a per-family read pilot with early exit: chr13's eight hotspot families align ~200k reads for ~120 assignments).
   Remaining gap to the 5,400 s bar is chr13 serialised after the rest (4,576 s); overlapping it needs its detect peak halved
   (family threads 4 x star threads 16, gw7 pending).
2. One read-star alignment per read per region, reused across families (the hotspot families' padded loci overlap heavily); may
   change products -> needs an identity check, a registration if not identical.
3. `candidate_copies` pass A from the sequential sweep (needed before the candidates lever can be tested at all).
4. One BGZF sweep per chromosome serving every family's light window pass in `copy_assign` (r1505: the 73 cross-chromosome families
   of the slim A119b run spend 1,427 s re-reading the same chr17/chr19 blocks; one multithreaded sweep of that BAM is ~240 s).
   Identical by construction if each window's records are bucketed in coordinate order. The multithreaded reader itself is used
   only on regions >= 2 Mb (1.5x there; 3.5x slower on windows, it restarts on every index seek).
