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
| `--region-threads N` / `RUSTLE_REGION_THREADS` | regions in parallel | the pipeline's assign uses 8 | identical (merge in region order) |
| slim BAM (`slimbam`, ratio 0.8) | drops secondaries no step reads | chr20: 12% fewer records; assemble 46 -> 36-44 s | assemble **identical**; assign **changes** (6 of 36 rows lost) |
| mm2-fast (`RUSTLE_MINIMAP2`, `~/tools_home/mm2fast/src/minimap2`, v2.24 + AVX-512) | accelerated chaining | MCL1 at 64 threads: read-star 299 s vs **256 s for our minimap2 2.30 default**; 2.30 with the `--max-chain-skip=1000000` mm2-fast requires takes 3,031 s, so mm2-fast only accelerates a setting we do not use. **Not a lever.** | **changes** (15,946 of ~31,800 assignment lines differ vs 2.30; statuses nearly equal) |

Threads multiply: region threads x family threads x star threads should not exceed the node's cores by much (the panel run at
8 x 8 x 8 was slower than 8 x 1 x 8).

## Not levers (measured)

Caching alignments by copy sequence (all 590 hotspot copies are distinct); `-G 100k` and query dedup (§6w9: slower / 5.5%
duplicates); Roaring bitmaps and smarter comparisons in assemble (register 983: BAM iteration is the wall); NumPy / bitsets / LRU
caches in copy_assign (our own code is ~5 s of a ~1,450 s family; the rest is minimap2); mm2-fast (slower than minimap2 2.30's
default chaining heuristic on our read-star job, see the table).

## Open levers, by expected gain

1. Done: auto family threads (the hotspot runs alone at the tail and fans out).
2. One read-star alignment per read per region, reused across families (the hotspot families' padded loci overlap heavily); may
   change products -> needs an identity check, a registration if not identical.
3. `candidate_copies` pass A from the sequential sweep (needed before the candidates lever can be tested at all).
4. One BGZF sweep per chromosome serving every family's light window pass in `copy_assign` (r1505: the 73 cross-chromosome families
   of the slim A119b run spend 1,427 s re-reading the same chr17/chr19 blocks; one multithreaded sweep of that BAM is ~240 s).
   Identical by construction if each window's records are bucketed in coordinate order. The multithreaded reader itself is used
   only on regions >= 2 Mb (1.5x there; 3.5x slower on windows, it restarts on every index seek).
