# PREREG r1600: assign reads the BAM once (router), schedules from exact counts (2026-11-09)

Design agreed with the user on 2026-11-09 (approach A of three; B = a slim BAM per run, C = repairing the r1507 window
store). This file is the spec and the decision rule. It is written and committed before step 0 and before any code.

## Why

On A119b (`audit/e2e2`, current defaults), assign takes 502 s of wall time at 8 region threads.

- Its 416 regions read 1,451,051 records through one indexed query per copy window. 231 of the regions are
  cross-chromosome families, and they hold 665,732 of those records.
- Summed over threads, reading costs 2,229 s (`[timing] reads_in_region`). Detection and assignment cost 1,051 s.
- Overlapping windows, and windows that fall in the same BGZF block, decode the same blocks again and again.
- Regions are admitted by a record budget charged from BAM-index estimates (0.6-0.85x of the real counts, r1543), so
  which regions run together depends on estimates and timing.

The window store (`RUSTLE_WINDOW_SWEEP=1`, r1507) sweeps once per chromosome per batch. At current defaults it was 5.3x
slower (audit round 12: 2,914 s vs 546 s), so it stays opt-in.

## Design

### 1. Window plan
Before any region runs, list every region's windows: the `region_windows` map, or the region itself when it has no windows.
Each window carries its region key and its position in that region's window list.

### 2. Router (new)
- From the cached BAI, take the chunks of every window and merge them into sorted, non-overlapping BGZF ranges.
- Split the ranges into N parts by compressed file offset. A record belongs to the part in which it starts.
- Each worker opens its own reader, seeks to its part and decodes it. For each record it finds every window that overlaps
  the record, using the noodles query test: same reference id, and the closed interval [alignment_start, alignment_end]
  intersects the window's [lo+1, hi]. The record is appended to those windows' lists.
- Raw record bytes are stored once; a window's list holds indices into the store.
- Parts are joined in part order. Every window's list is therefore in file order, which is the order its own indexed query
  returns, at any worker count.

### 3. Conversion split out of the query
`BamIndexCache::reads_in_region_from` (`src/family.rs`) becomes two parts:

- the query, which yields records;
- `convert_window_records(records, chrom, lo, hi, mode)`, which holds every per-window rule: the
  `RUSTLE_GTF_SECONDARY_AS_RATIO` buffered filter, the read-home table and its report, and the Full / Light / Only(names)
  modes.

Under the router, a region calls the conversion on its routed lists. Nothing from `pool_windows` onward changes.

### 4. Scheduling
- Region cost is its exact routed record count, summed over its windows.
- Start order is by count, descending; ties go by position in `flat`.
- The `RUSTLE_REGION_READ_BUDGET` permit pool and the "above the budget runs alone" rule stay, charged with the exact count.
- `estimate_records` is removed.
- Wall time still varies with compute time. The start order, the charges and every region's input are fixed by the input.

### 5. Memory
- The store keeps raw BAM records.
- A region's references are released when the region finishes, and a record is freed when no window holds it any more.
- If step 0 predicts a peak above today's assign peak, the spec gains a file-backed store (written to the run directory and
  memory-mapped) before implementation, recorded as an addendum here. Otherwise the store stays in memory.

### 6. Failure
- A read or decode error fails assign before any table is written. The error names the file offset and the window.
- There is no fallback to per-window queries.
- The path with no `.bai` (`reads_in_region`, decoding the file sequentially) is unchanged.

### Not changed
- the assemble stage's reading
- the streaming path
- `--tied-prefilter`, which uses the same conversion with the Light / Only modes
- the readthrough buffer (`rt_buf`), which is fed from the converted records as today

## Step 0 (measure before code; no bars)
On A119b, with e2e2's families products:

1. today's assign peak RSS at defaults (`audit/rss.py`), with its wall;
2. the compressed bytes covered by the merged chunk ranges of every window, from the index only;
3. the routed record count, with duplicates and distinct, and the raw record bytes: one decode pass over the merged ranges
   by a throwaway probe;
4. the router's predicted peak, store bytes plus today's peak minus today's read buffers. If it is above (1), add the
   file-backed store as an addendum here before any implementation.

## Tests (code)
1. **Router equals queries.** On a fixture BAM, every window's routed list equals its own noodles query, record for
   record, byte for byte. The windows cover:
   - a record ending one base before the window, and one ending on its first base;
   - a record starting on the window's last base, and one starting one base after it;
   - one record in two windows;
   - windows on two references;
   - an unmapped record that carries a position;
   - part splits at 1, 2 and 7 workers.
2. **Thread independence.** Router output is byte-identical at 1 and 8 workers.
3. **Full suite.** `cargo test --profile dev-opt` passes. The streaming, `--tied-prefilter` and no-`.bai` tests stay.

## Acceptance (all must pass; registered before any arm runs)
| # | run | bar |
|---|---|---|
| A1 | A119b assign, `RUSTLE_ASSIGN_ROUTER=1`, defaults, on e2e2's families products | all 7 assign tables byte-identical to e2e2 |
| A2 | same run | assign wall <= 250 s (e2e2: 502 s) |
| A3 | same run | peak RSS <= step 0 (1) |
| A4 | A119b assign, router, `--threads 16` | all 7 tables byte-identical to e2e2 |
| A5 | gorilla KB3781 assign, router, the 2x2 grid region-threads {1, 8} x `--threads` {16, 72} (4 runs) | every table byte-identical to e2e_kb |

- Tables (A119b): `assign.assignments.tsv`, `assign.copy_conversions.tsv`, `assign.famcn_readonly.tsv`,
  `assign.families.tsv`, `assign.family_join.tsv`, `assign.quant.tsv`, `assign.params.tsv`. `params.tsv` records
  settings, not the environment, and the router adds no row, so it must match outright.
- Any table that differs is a FAIL, whatever the speed; there is no "explained difference" exception.
- Identical tables with A2 or A3 missed is a FAIL as registered, recorded with its numbers.
- Gorilla stands in for A119b's thread matrix because A119b at region-threads 1 takes hours.

## Rollout
- **During development:** the router is opt-in, `RUSTLE_ASSIGN_ROUTER=1`, so every comparison uses one binary.
- **PASS:**
  - the router becomes the default and the variable is removed;
  - deleted: the window store (`RUSTLE_WINDOW_SWEEP`, `RUSTLE_SWEEP_BATCH`), `estimate_records`, and the r1548
    scheduler extras (`RUSTLE_DISPATCH_DEFER`, the head gate, `RUSTLE_BUDGET_SHRINK`);
  - PERFORMANCE.md and the O2 ledger are updated.
- **FAIL:**
  - a negative-results register row, with the numbers;
  - the router code is removed, not kept as an option.

## Not covered
- HG002: a 131 GB BAM, and assign there needs > 250 GB.
- HG002 and OR6737 are not in the bars.

## Outcome
(empty until the acceptance runs)
