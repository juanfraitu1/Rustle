# Assign BAM router (r1600) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** `copy_assign --families` reads every copy window's records in one parallel pass over the BAM, and then
schedules regions by their exact record counts. The output stays byte-identical.

**Architecture:**
- **Router.** A new library module, `src/bam_router.rs`:
  - merges every window's BAI chunks into one set of byte ranges;
  - decodes those ranges in N parts by file offset;
  - gives each window its overlapping records, in file order.
- **Conversion.** The per-window record-to-read conversion comes out of
  `BamIndexCache::reads_in_region_from` (`src/family.rs`). Today's query path and the router then share it.
- **Wiring.** `copy_assign` uses the router under `RUSTLE_ASSIGN_ROUTER=1` until the acceptance runs decide.

**Tech Stack:** Rust; noodles-bam 0.70, noodles-bgzf 0.33, noodles-csi 0.40, noodles-core 0.15; bash and Slurm for the
acceptance runs.

**Spec:** `docs/PREREG_assign_router_2026-11-09.md` (r1600). It holds the decision rule, and the bars are copied below. **Addendum 2 (two passes) governs Tasks 2-4; their text below is revised to match (2026-11-10).**

## Global Constraints

- **Output.** Every assign table must be byte-identical to the query path. No "explained difference" is accepted.
- **Bars (A119b assign):**
  - wall <= 250 s (e2e2: 502 s);
  - peak RSS <= today's peak, measured in Task 1;
  - all 7 tables byte-identical to `audit/e2e2`.
  - The 7 tables: `assign.assignments.tsv`, `assign.copy_conversions.tsv`, `assign.famcn_readonly.tsv`,
    `assign.families.tsv`, `assign.family_join.tsv`, `assign.quant.tsv`, `assign.params.tsv`.
- **Overlap test.** Exactly noodles' `intersects`:
  - the same reference id;
  - `(alignment_start..=alignment_end)` intersects `Region::interval()`, the region parsed from `"{chrom}:{lo+1}-{hi}"`;
  - a record without a start or an end never intersects.
- **No fallback.** A router error fails assign. The error names the compressed offset and the window.
- **Environment:** the opt-in variable is `RUSTLE_ASSIGN_ROUTER=1` during development. On a pass it is removed.
- **Housekeeping (user rules):**
  - Big outputs go under `/anvme/workspace/v138da13-rustle/audit/r1600/`, never `$HOME`.
  - Jobs set `PATH=$HOME/.pixi/detached-envs/rustle-7434162673922408322/envs/default/bin:$PATH`.
  - Peak RSS comes from `audit/rss.py`, since the nodes have no `/usr/bin/time`.
  - Before every commit, check `git diff --cached --stat` and stage files by name; the index is shared with other
    sessions.
  - Commit messages end with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
  - Push nothing unless the user asks.
- **Tests:** iterate with `cargo test --profile dev-opt`. Timing numbers come only from `--release` binaries.

## Review Focus

1. **A window on a contig missing from the BAM header.** Today's query fails with noodles' "invalid reference sequence
   name". The router must fail on the same window, before decoding anything. Test in Task 3.
2. **Empty or inverted windows** (`hi <= lo`). They must give exactly what the query gives, whether that is an empty list
   or an error. Test in Task 3.
3. **The same window listed twice** (one region, or two regions sharing a window). Each listing gets the full list. Test in
   Task 3.
4. **A record whose BGZF block starts in one part and continues into the next.** It is decoded once, by the part where it
   starts. The 7-worker split in Task 3 forces part boundaries inside a block.
5. **A region with no window entry** (`region_windows` has no key, so the region itself is its one window). Under the
   router it must get the same records as its query. Covered by Task 4's fixture runs with `--region`.

---

### Task 1: Step 0, measure before code (throwaway probe)

**Files:**
- Create, not committed: `/anvme/workspace/v138da13-rustle/audit/r1600/step0/` (job script, window dump, probe)
- Modify, not committed (a local patch reverted after this task): `src/copy_assign_bin.rs` near line 7362. When
  `RUSTLE_DUMP_WINDOWS=PATH` is set, write `region\twindow_chrom\tlo\thi` for every window of every region in `flat`
  order, then exit 0 before any region computes.
- Modify: `docs/PREREG_assign_router_2026-11-09.md`, adding a "Step 0 RESULT" section.

**Interfaces:**
- Produces: the four step-0 numbers in the prereg, and the decision whether a file-backed store addendum is needed.

- [ ] **Step 1: Measure today's peak.** Run a job on A119b: assign alone, from e2e2's families products, release binary,
  under `audit/rss.py`. Expected: wall within 10% of 502 s, and a "Maximum resident set size" line. Record both.
- [ ] **Step 2: Dump the windows.** Apply the local dump patch, build `--release`, and run with
  `RUSTLE_DUMP_WINDOWS=…/windows.tsv`. Expected: 416 distinct region keys.
- [ ] **Step 3: Run the probe.** It is a throwaway `examples/route_probe.rs`, not committed. It opens the BAI with
  noodles, collects `index.query(ref_id, interval)` chunks for every window and merges them. It prints:
  - the total compressed bytes;
  - one decode pass over the merged ranges, counting per window the records that pass the overlap test;
  - the sum of counts (with duplicates), the distinct records, and the sum of raw record bytes (`record.as_ref().len()`).

  Expected: the summed counts equal 1,451,051 minus any records the conversion drops. Note the difference rather than
  forcing it to match.
- [ ] **Step 4: Decide the store.** Predicted router peak = Step 1 peak - R + D, where:
  - D is the distinct raw bytes from Step 3;
  - R is the raw bytes of the 8 regions with the most routed raw bytes, which today hold their records at once at region
    threads 8.

  If the predicted peak is above the Step 1 peak, write the file-backed store addendum into the prereg, and stop the plan
  for the user's review. Otherwise write "in-memory store" and the numbers.
- [ ] **Step 5: Clean up and commit.** Revert the dump patch and delete the example (`git diff --stat` shows only the
  prereg). Commit the prereg with message `r1600 step 0: today's assign peak, window bytes, routed records`.

### Task 2: Split the per-record conversion out of the query

**Files:**
- Modify: `src/family.rs`, the `ratio <= 0.0` loop of `BamIndexCache::reads_in_region_from` (near 51080).
- Test: `tests/bam_window_convert.rs` (new).

**Interfaces:**
- Produces, on `impl BamIndexCache` in `rustle::denovo_assemble`:
  - `pub fn convert_record(&self, record: &noodles_bam::Record, chrom: &str, mode: ReadMode<'_>) -> Result<(Option<PrimaryRead>, Option<BamRead>)>`.
    It holds the body of that loop when no read-home table is set: the `Only` name skip, `RecordBuf` conversion,
    `alignment_read_from_record`, `bam_read_for`. The loop calls it; the home-table and buffered-ratio paths are untouched.
  - `pub fn header(&self) -> &noodles_sam::Header`
  - `pub fn index(&self) -> &noodles_bam::bai::Index`

- [ ] **Step 1: Write the failing test.** `convert_record_equals_reads_in_region`. For each fixture BAM
  (`tests/fixtures/{xfam_conflict,o2_origin_resolve,candidate_copies,same_chrom_supplement}/reads.bam`) and each `@SQ`
  contig, use 3 windows: the whole contig, `[100, 2000)` and `[0, 1)`. Run a noodles query, map every record through
  `convert_record`, and collect the `Some`s.
  - The result must equal `cache.reads_in_region_mode(bam, chrom, lo, hi, mode)`, compared with `format!("{:?}")`.
  - Check this for `Full` and for `Light`, and for `Only(names)` with the names of the first two records.
- [ ] **Step 2: Run it and see it fail.** `cargo test --profile dev-opt --test bam_window_convert` → compile error: no
  `convert_record`.
- [ ] **Step 3: Implement.** Move the loop body into `convert_record`; the logic is unchanged.
- [ ] **Step 4: Run the tests.** The new test passes, and the full `cargo test --profile dev-opt` passes with the same
  count as before plus 1.
- [ ] **Step 5: Commit** `family.rs: per-record conversion split from the BAM query (r1600 task 2, no output change)`.

### Task 3: The router

**Files:**
- Create: `src/bam_router.rs`, registered in `src/lib.rs` as `pub mod bam_router;`.
- Test: `tests/bam_router.rs` (new). It writes its own small BAM with noodles plus `noodles_bam::fs::index` into
  `CARGO_TARGET_TMPDIR`.

**Interfaces:**
- Consumes: `BamIndexCache::{header, index}` from Task 2.
- Produces, in `rustle::bam_router`:
  - `pub fn route<T: Clone + Send>(bam_path: &str, cache: &BamIndexCache, windows: &[(String, u64, u64)], workers: usize, convert: &(dyn Fn(&noodles_bam::Record, &[usize]) -> anyhow::Result<Option<T>> + Sync)) -> anyhow::Result<Vec<Vec<T>>>`.
    - The result has one list per input window, in input order, each in file order.
    - `convert` is called once per record that overlaps at least one window. Its second argument is the sorted list of
      overlapping window indices.
    - A `Some(t)` is cloned into every one of those windows; `None` adds the record to none of them.
    - A caller that needs per-window selection (pass 2) returns `None` when no listed window wants the record, and filters
      per window afterwards with the indices it already knows. Simplest form: return `Some(t)` and let pass 2's caller
      drop the copies whose region's tied set lacks the name.
- Algorithm (the part a signature does not fix):
  1. Per window, parse `"{chrom}:{lo+1}-{hi}"`, resolve the reference id through the header (an error names the window),
     and collect `index.query(id, interval)` chunks.
  2. Merge every chunk into sorted, non-overlapping `(start_vpos, end_vpos)` ranges.
  3. Split the ranges into `workers` contiguous groups of about equal compressed bytes. Cut only at range boundaries.
  4. Each worker opens its own `noodles_bam::io::Reader` (single-threaded BGZF), seeks to each of its ranges, and reads
     while `virtual_position < end`. Overlapping windows come from per-reference windows sorted by `lo`, with a prefix
     max of `hi`. The overlap test is noodles' `intersects`.
  5. Each worker returns its `(window, T)` pushes in its own order. Groups are concatenated in group order, which is file
     order.
  - Errors carry `format!("routing window {chrom}:{lo}-{hi} at compressed offset {}", vpos.compressed())`, or the offset
    alone for a decode error.

- [ ] **Step 1: Write the failing tests.**
  - `routed_equals_query_on_edge_windows`. Build a BAM with these records on `c1`:
    - `[1000, 1999]`;
    - `[2000, 2999]`;
    - a spliced record `[1500, 9000]` with `N`;
    - one record on `c2`;
    - an unmapped record placed at `c1:1500`.

    Windows: `[999, 1000)`, `[1000, 1001)`, `[1999, 2000)`, `[2000, 2001)`, `[3000, 3001)`, `[0, 20000)` on `c1`, one
    window on `c2`, and `[1500, 1500)` (inverted or empty). For every window, the routed list (record bytes) must equal
    that window's noodles query, both as a record sequence and as a `Vec<Vec<u8>>`.
  - `routed_is_identical_at_1_2_7_and_8_workers`: on the same BAM plus 5,000 random records across 40 BGZF blocks (seed
    fixed), the `route` output is the same for every worker count.
  - `convert_is_called_once_per_record`: count the calls with an atomic counter; the count must equal the number of
    distinct records overlapping any window.
  - In every test, `T` is the record's `RecordBuf` rendered with `{:?}`. The query side is rendered the same way.
  - `duplicate_windows_each_get_the_full_list`: the same window given twice → both lists equal.
  - `unknown_contig_fails_naming_the_window`: a window on `c9` → `route` errors, and the message contains `c9`.
- [ ] **Step 2: Run them and see them fail.** `cargo test --profile dev-opt --test bam_router` → unresolved
  `rustle::bam_router`.
- [ ] **Step 3: Implement `route`** as specified above.
- [ ] **Step 4: Run the tests.** All 4 pass. Then run `routed_is_identical_at_1_2_7_and_8_workers` 3 times; it passes
  each time.
- [ ] **Step 5: Commit** `bam_router: one pass over every window's BAI chunks, file order at any worker count (r1600 task 3)`.

### Task 4: Wire the router into copy_assign behind `RUSTLE_ASSIGN_ROUTER=1`

**Files:**
- Modify: `src/copy_assign_bin.rs`:
  - near 7338 (`bam_cache`): build the window plan and call `route`;
  - near 7576-7615: the region budget charged with the exact count;
  - near 7762-7790: the window fetch takes routed lists;
  - near 8741-8770: the start order.
- Test: `tests/copy_assign_router.rs` (new). It follows `tests/copy_assign_xfam.rs::determinism_across_repeats_and_region_threads`.

**Interfaces:**
- Consumes: `route` (Task 3); `convert_record` (Task 2).
- Behaviour when `RUSTLE_ASSIGN_ROUTER=1` and the Addendum 2 scope holds (`--tied-prefilter`, a `.bai` present, not
  streaming, `gtf_secondary_as_ratio() <= 0.0`, `read_home_table().is_none()`, `RUSTLE_WINDOW_SWEEP` unset):
  - **Plan.** For each region in `flat` order, its `region_windows` entry, or `vec![region]`, concatenated, with
    per-region offsets.
  - **Pass 1, before the batch loop.** Call `route` with `convert = cache.convert_record(rec, chrom_of(rec), Light)`.
    Then, per region and in parallel over regions, run exactly today's code from `pool_windows` through the tied gate.
    The code that today reads the windows takes the routed lists instead, and the readthrough buffer is fed from them in
    window order as today. Keep per region `(pr, light_spans, tied)` plus anything else the existing code computes
    there.
  - **Pass 2.** Take the union of the tied sets. `convert` returns `None` when the record's name is not in the union,
    and otherwise `convert_record(rec, chrom, Full)`. Afterwards, per window, keep only the items whose name is in that
    window's region's tied set. Then apply today's `fseen` de-duplication per region.
  - **compute.** It takes its region's prepared `(pr, full, light_spans)` from a `Mutex<Vec<Option<...>>>` with `take()`,
    in place of reading. Everything after is unchanged.
  - **Budget and order.** The budget is charged the region's pass-1 routed count. The start order is that count,
    descending, with ties by `flat` index.
  - **Timing lines:** `[timing] router pass1`, `[timing] router gate`, `[timing] router pass2`.
- Outside the scope, every line behaves as today.

- [ ] **Step 1: Write the failing test.** `router_arm_is_byte_identical` in `tests/copy_assign_router.rs`.
  - For both fixtures (`xfam_conflict` with its `--region` arguments copied from `copy_assign_xfam.rs`, and
    `o2_origin_resolve` copied from `tests/o2_origin_resolve.rs`), run `copy_assign` with `RUSTLE_ASSIGN_ROUTER` unset and
    with it set to `1`.
  - Each run uses `--region-threads 1` and `--region-threads 4`.
  - Assert that every output file (`OUT.*`) is byte-identical across all 4 runs of a fixture.
  - Also assert that the router arm's stderr contains `[timing] router pass2` when `RUSTLE_TIMING=1`. This proves the
    arm took the router.
  - Both arms pass `--tied-prefilter`, the pipeline's setting.
- [ ] **Step 2: Run it and see it fail.** `cargo test --profile dev-opt --test copy_assign_router` → FAIL on the
  `[timing] router` assertion.
- [ ] **Step 3: Implement the wiring** as above.
- [ ] **Step 4: Run the tests.** The new test passes. Then run the full `cargo test --profile dev-opt` and
  `cargo build --release`.
- [ ] **Step 5: Commit** `copy_assign: RUSTLE_ASSIGN_ROUTER=1 reads every window in one routed pass, budget and order from exact counts (r1600 task 4)`.

### Task 5: Acceptance runs (A1-A5, as registered)

**Files:**
- Create: `/anvme/workspace/v138da13-rustle/audit/r1600/acc/{a1,a4,a5}/run.sh`
- Modify: `docs/PREREG_assign_router_2026-11-09.md`, its `## Outcome` section.

**Interfaces:**
- Consumes: the release binaries from Task 4, frozen by copying `target/release` to
  `/anvme/workspace/v138da13-rustle/target_r1600`; e2e2's families products; e2e_kb's gorilla products.

- [ ] **Step 1: Submit A1-A3.** A119b assign with `RUSTLE_ASSIGN_ROUTER=1 RUSTLE_TIMING=1`, driver defaults, under
  `rss.py`. Compare the 7 tables with `cmp` against `audit/e2e2/a1.assign.*`, and print wall and peak.
- [ ] **Step 2: Submit A4.** The same with `--threads 16`, comparing the 7 tables.
- [ ] **Step 3: Submit A5.** Gorilla KB3781 assign over the 2x2 grid, region-threads {1, 8} x `--threads` {16, 72}.
  Each run compares every `kb.assign.*` table with `audit/e2e_kb`.
- [ ] **Step 4: Apply the registered rule.**
  - PASS only if A1, A4 and A5 are all identical, A2's wall is <= 250 s, and A3's peak is <= the step-0 peak.
  - Write the outcome with every number into the prereg, and add an O2 ledger line.
- [ ] **Step 5: Commit** `r1600 outcome: <PASS|FAIL as registered> (A1..A5 numbers)`.

### Task 6: Rollout, by outcome

**Files (PASS):**
- `src/copy_assign_bin.rs`:
  - the router becomes unconditional when a `.bai` is present and the run is not streaming;
  - `RUSTLE_ASSIGN_ROUTER` is removed;
  - the window store is deleted (`WindowStore`, `build_store`, `sweep_on`, `RUSTLE_WINDOW_SWEEP`, `RUSTLE_SWEEP_BATCH`);
  - so are `MAY_DEFER`, `IS_HEAD`, `RETRY`, the head gate, the `deferred` queue, `RUSTLE_DISPATCH_DEFER` and
    `RUSTLE_BUDGET_SHRINK`.
- `src/family.rs`: `estimate_records`, if no other caller remains (grep first).
- Tests touching deleted variables: grep `tests/` for each name and update or delete those tests.
- `docs/PERFORMANCE.md`: replace the window-store and region-budget rows with one router row.
- `docs/MODULE_STATUS.md`: update, then run its enforcing test.

**Files (FAIL):**
- `docs/NEGATIVE_RESULTS_REGISTER.md`: a row 1600 with the numbers.
- Remove `src/bam_router.rs`, its test and the Task 4 wiring. Keep Task 2's split if it stayed byte-identical (it is
  smaller code either way).

- [ ] **Step 1: Grep.** Grep every name to delete across `src/ tests/ tools/ docs/ bench/`, and list the hits.
- [ ] **Step 2: Apply the outcome's file list.**
- [ ] **Step 3: Run the tests.** `cargo test --profile dev-opt` passes, then the release build.
- [ ] **Step 4: Recheck on PASS.** Rerun A1 with the default binary (no variable). The 7 tables must be byte-identical
  to e2e2.
- [ ] **Step 5: Commit** `assign: router is the default; window store, index estimates and r1548 scheduler extras removed (r1600 PASS)`,
  or the FAIL equivalent.
