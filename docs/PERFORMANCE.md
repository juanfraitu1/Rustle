# Pipeline performance: where the time goes and the switches that change it

Living document (started 2026-11-02). Every number here was measured on Fritz (Intel Xeon Platinum 8360Y, 72 cores, `singlenode`)
on A119b (human testis, CHM13) unless stated. Switches marked **identical** were checked to leave every product byte-identical.

## Where the time goes

Current (`all` on A119b, 72 threads, `audit/e2e0`, 2026-11-09; censat exclusion and SD damping on):

| stage | A119b wall | what dominates |
|---|---|---|
| assemble | 595 s | BAM inflate and record decoding (audit round 48: /anvme gives ~1 GB/s in total, a warm cache saves only ~25 s of a 200 s pass) |
| families | 465 s -> **385 s** (audit 22: minimap2 `-K 4G`, one query batch, identical) | all-vs-all minimap2 over locus spans: the index build (~90 s on ~3 threads) and the last batch's tail leave the node idle |
| assign | 549 s | the last region, `~xchrom~MCL13` (150 copies): read-star 123 s of its 150 s |
| index | 124 s, once per genome | minimap2 `-d`; shared next to the FASTA |
| flag scan | 283-314 s, peak 16.6 GB | BAM sweeps. Audit rounds 17-20: before them the old scan ran over 5 h on chr13 alone and hit the job limit |
| flag align + verdict | 41 s, peak 41 GB | the 12.3 GB splice index in memory |
| **`all` end to end, from scratch** | **1,328 s (22 min)**, `audit/e2e2` (1,572 s in e2e1, before the background flag scan of audit 39). The scan runs under families and assign only when assemble runs in the same call: a resumed run that skips assemble, or `--low-mem`, scans after assign (about +290 s) | assemble 397 (as_table cached; +~190 s once per new BAM), families 372, assign 475, index 0 (shared), flag 327. Every product identical to e2e0. Gorilla KB3781 with no extra env: 756 s |
| `all --low-mem` (gorilla KB3781) | 1,549 s, node peak 43.7 GB in families (audit 40); with the caps of audits 43 and 45 families ~25 GB node, flag align 24.9 GB RSS | the split index changes 7 of 16,467 flag rows, no verdict class; every other product identical. Families' minimap2 memory scales with threads (gorilla 38 GB at 72, 18 GB at 16, 14 GB at 8; `-K` 4G vs 500M changes it by 2 GB; PAF identical in all, audit 43), so since audit 43 `--low-mem` runs families on at most 16 threads (127 -> 219 s on gorilla), and since audit 45 flag's align too (43.6 -> 25.4 GB, 171 -> 203 s) |

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
| `RUSTLE_EXCLUDE_SOFTMASK=0.9` (opt-in, register 1690) | assign skips families whose copies are mostly soft-masked (>= that share of exon bases lowercase in `--fasta`; `tools/softmask_families.py`). The BED-free route for genomes without a censat annotation: it takes satellite and transposon-class families, keeps the transition-zone families the censat BED removes, and does not take rRNA (for that: `tools/sequence_families.py PREFIX.fam.base.copies.fa RDNA.fa`, register 1691; since wired as `RUSTLE_EXCLUDE_FASTA`, next row) | A119b assign alone: 601 s against 346 s with the censat BED; 323 families skipped against 628; 5,788 assigned rows against 5,399; rows on families both keep 5,227 against 5,251 |
| `RUSTLE_EXCLUDE_FASTA=REF.fa` (opt-in, register 1691) | assign skips families whose copies are mostly covered (>= 50% of the copy) by minimap2 asm20 alignments to REF.fa (`tools/sequence_families.py`). For rRNA / rDNA families, which the soft mask does not reach: `tools/data/rdna_human.fa` (human 45S unit + spacer + 5S; also finds gorilla rDNA families) or `tools/data/rdna_human.masked.fa`. Joined with the BED and soft-mask lists; needs python3 and minimap2. | HG002 on its paternal assembly (row 1741, 2026-11-15, one run per arm at the driver's thread defaults on a 250 GB node): with the switch `assign` finished in 1,631 s at 31.5 GB (20 of 537 families skipped); without it the run was OOM-killed at 1,128 s with a 239.8 GB process. An older default run (2026-10-31, other binaries / settings) did complete, so this is one kill, not a limit; whether the memory comes from the rDNA families or from large families in parallel is not known. On that older run the listed families hold 18% of assignment rows and 0 assigned rows. Still opt-in. |
| `RUSTLE_SOLE_CANDIDATE=1` (opt-in, the default before 2026-11-13) | assign reports a contested molecule with exactly one candidate copy as `assigned` (`sole_candidate = 1`). Default now: `tied` (user ruling; registers 787 / 1656 / 1695). A119b: 5,784 assigned rows with the switch, 3,555 without; no other row changes. | |
| stage `repeatloci` (opt-in, rows 1698 / 1702) | `rustle_pipeline.sh ... repeatloci` writes `PREFIX.repeat_loci.tsv`: loci whose exons are >= 0.9 soft-masked or >= 50% covered by the rDNA unit (`RUSTLE_REPEAT_FASTA`). Marks only; no other product changes. | KB3781: 46 of 14,860 loci marked |
| assign defaults by machine (audit round 54) | `rustle_pipeline.sh` picks 24 regions and a 16M-record budget when the node has >= 48 cores and >= 128 GB and `RUSTLE_EXCLUDE_BED` or `RUSTLE_EXCLUDE_SOFTMASK` is set; otherwise 8 / 4M. Both variables still override. The rule tests only that an exclusion BED is set, not what is in it: it was measured with the T2T censat BED, and a BED that leaves the satellite families in should come with `RUSTLE_REGION_THREADS=8` | A119b assign alone, tables identical in every arm: 8 / 4M 484 s (34 GB peak), 16 / 8M 349 s (40 GB), 24 / 16M 286 s (40 GB), 36 / 32M 281 s (43 GB). Raising regions without the budget does little (16 / 4M 444 s, 24 / 4M 435 s) |
| assemble defaults by machine (audit rounds 50, 58) | 24 contigs swept at once under an 80M-record budget on a node with >= 48 cores and >= 128 GB, else 8 / 20M; inflate workers per reader scale with the cores (`RUSTLE_BAM_THREADS`, a quarter of the cores, 4 to 16) | A119b assemble alone, GTFs identical: workers 4 / 8 / 16 = 397 / 309 / 275 s; then contigs 8 / 16 / 25 = 291 / 237 / 212 s (9 / 11 / 13 GB per process) |
| flag in the background (audit round 55) | `all` runs the whole flag stage (scan, align, verdict) beside families and assign, unless `--candidates`, `--low-mem` or a first-run index build | A119b: 42 s no longer after assign |
| `MALLOC_ARENA_MAX=2` | glibc arena cap | RSS flat at 88 GB but ~8x slower (110 regions in 92 min vs 935). **Not a lever.** | — |
| slim BAM (`slimbam`, ratio 0.8) | drops secondaries no step reads | chr20: 12% fewer records; assemble 46 -> 36-44 s | assemble **identical**; assign **changes** (6 of 36 rows lost) |
| mm2-fast (`RUSTLE_MINIMAP2`, `~/tools_home/mm2fast/src/minimap2`, v2.24 + AVX-512) | accelerated chaining | MCL1 at 64 threads: read-star 299 s vs **256 s for our minimap2 2.30 default**; 2.30 with the `--max-chain-skip=1000000` mm2-fast requires takes 3,031 s, so mm2-fast only accelerates a setting we do not use. **Not a lever.** | **changes** (15,946 of ~31,800 assignment lines differ vs 2.30; statuses nearly equal) |

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

Assign router (r1600, 2026-11-09): reading the BAM once and scheduling from exact counts was slower and bigger than the default (583 s vs 553 s; 58.5 vs 36.4 GB), tables identical. Dropped.
- **Where the BAM lives matters more than any switch (row 1741, 2026-11-15):** `assign` reading `hg002.pat.bam` from `~/work` stalled on three regions of 8-10 M records (under 1 MB/s, one thread in disk wait, 88 minutes without finishing). The same regions finished in the 27-minute run once the BAM was copied to the workspace (`/anvme/.../stage/`, copy 1,646 s for 131 GB). Stage the BAM on `/anvme` before a whole-sample `assign`.

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
5. Data-structure survey (2026-11-25; grep counts and a clippy pass with the performance and ownership lints, 110 hits in 134k
   lines; nothing measured, nothing changed). None of these is expected to move wall time (inflate and minimap2 are the bound);
   the gains are memory per record and type safety. In order of payoff:
   - `src/copy_assign_bin.rs:8861`: `primary.clone()` deep-copies the region's whole primary read pool under `--assemble-only`
     and the original is never used again (clippy `redundant_clone`). A move is a one-word change, identical by construction.
   - `chrom: String` in every read and locus struct (35 structs, 181 `chrom.clone()` / `to_string()` sites; pass 1 clones it per
     junction to probe a `BTreeMap<(String, (u64, u64)), u32>`) -> a `u32` contig id. Ids assigned in lexicographic name order
     keep every `BTreeMap` iteration order, so products can stay identical. Estimate: ~50 bytes and one allocation per
     `PrimaryRead`.
   - `lazy_reads` (`src/family.rs:50917`): the op list is cloned into a `Cigar` for `exons_from_cigar`, and the read name is
     allocated before the `Only(names)` filter. Two allocations per record; about ten lines; the test
     `lazy_and_buffered_readers_agree_on_clipped_reverse_tagged_records` covers it.
   - `cigar: Vec<(char, u64)>` (`AlignedRead`, `PileRead`): 16 bytes per op -> a `#[repr(u8)]` op enum + `u32` length, 8 bytes,
     with exhaustive `match` instead of `_ => {}`.
   - `Vec<Option<u8>>` for PSV alleles and observations (99 sites): 2 bytes per cell -> `Option<NonZeroU8>` or a base enum, 1 byte.
   - Strand held as `char`, `String`, `Option<String>`, `Vec<char>` and `u8` (68 sites) -> one enum. `(u64, u64)` and
     `(String, u64, u64)` tuples (about 440 and 121 sites) -> half-open interval newtypes with one overlap rule (the helpers mix
     closed and half-open today). Correctness, not speed.
   - Maps and sets keyed by `String` (322 sites; not classified into read names and other ids) -> interned ids or `read_name_hash`.
   - Rest of the clippy list (45 `push_str(&format!())`, 13 `a = b.clone()`, 13 by-value parameters only read, 9 `&mut Vec`
     where a slice does, 6 needless `collect()`): cosmetic. Rerun: `cargo clippy --lib --bins -- -A clippy::all -W clippy::perf
     -W clippy::redundant_clone -W clippy::needless_collect -W clippy::needless_pass_by_value -W clippy::ptr_arg`.
   - Not zero-cost, so not in this list: packing `seq` to 2 or 4 bits and keeping qualities only at PSV columns (the two vectors
     are the memory of a full read). Interval trees and rust-lapper: not a lever, the lookups that matter already use a sorted
     vector with `partition_point`.
