# Audit rounds, November 2026 (ponytail-audit, goal: 100 rounds)

Goal (user, 2026-11-09): run 100 rounds of `/ponytail:ponytail-audit`. Apply the fixes that bring the pipeline closer to
its objectives: fast, and self-contained (one run directory, stage hatches, Rust hot paths, Python for orchestration and
scoring). The suggestions of the last rounds need not be applied.

Load assumed: one analyst per run, on one Slurm node (72 cores, ~250 GB). Inputs are a 20-100 GB Iso-Seq BAM and a
3 Gbp genome. Correctness means byte-identical products unless a round says otherwise.

Each round lists its findings and what was done. "Applied" means changed and verified; the check is named.

## Round 1: `tools/rustle_pipeline.sh` and the binaries of the default path

Default path: `all` = assemble -> families -> assign -> index -> flag.

1. **flag scans contigs one after another** (`src/missing_copy_flag.rs` `scan`). Each contig is an independent BAM
   sweep with its own piles and rows. The scan was ~90% of the post-assemble wall (OR6737 654 s after the index-once
   fix). Applied: contigs are scanned in parallel on a rayon pool of `--threads`, each with its own BAM reader, and
   joined in the original order, so `scan.tsv` and `consensus.fa` are unchanged. Check: on 8 KB3781 contigs, 32
   threads, both files byte-identical to the HEAD binary's (8,075 loci, 154 consensus sequences). The scan took 213 s
   on HEAD and 97 s with the change (`audit/r01`).
2. **families silently runs ~5x slower without minimap2-gate.** `mcl_families` falls back to aligning every pair and
   says so only in `PREFIX.families.log`. Products are the same, but the wall is ~5x; it happened to us on 2026-11-08.
   Applied: the driver prints a WARNING when the log carries the fallback note. Check: `bash -n`.
3. **The layered assembly is a second full `copy_assign` run** (`stage_assemble`, `PREFIX.final.gtf`). It doubles
   assemble CPU and BAM I/O to get one extra GTF. Not applied: a single run writing both GTFs needs the binary to
   branch after the shared passes. Candidate for a later round.
4. **`merged` and `all` are two genome-wide entry points.** `merged`'s assign sweeps every contig of the BAM header,
   without `--tied-prefilter`. `all` reads only the families' regions, with the prefilter. Not applied: switching
   `merged` to family regions is not proven byte-identical (the tied-outside-catalog registry sees different reads).
   Candidate: retire `merged` once nothing depends on it.

## Round 2: the index stage and the layered assembly

5. **Every run directory builds its own splice index** (`stage_index`, `PREFIX.splice.mmi`, ~13 GB for CHM13). The
   index depends only on the FASTA, the part size and minimap2. Applied: shared as `<FASTA>.splice.I<GB>G.mmi` in
   `RUSTLE_INDEX_DIR`, or next to the FASTA, with a key (path, size, mtime, part size, minimap2 version). Falls back to
   `PREFIX.splice.mmi` when the FASTA lies under `$HOME` (storage rule) or its directory takes no new file. Built under a
   temporary name and renamed. Check: `index_path` on the four cases (writable, symlinked into /anvme, explicit dir,
   read-only dir).
6. **The layered assembly re-runs the whole assembly** (round 1, item 3). The four layers act after assembly, inside
   the per-region polish loop (`copy_assign_bin.rs` ~10,860-11,110), on the region's transcript list. Emitting both
   GTFs from one run would save one full BAM pass and half of assemble CPU. Not applied this round: it needs a cloned
   per-region list and a second writer through the TPM and subchain steps. Kept for a later round, with a byte-identity
   check on both GTFs.

## Round 3: the driver's own tests

7. **Three driver tests failed since the share edge rule became the default** (`tests/pipeline_driver_units.rs`,
   since e08aab2b). The stub `mcl_families --help` did not name `--min-exon-share`, so the driver refused it. The
   expected command still carried the coverage rule. The tests also did not pass `--threads`, so after 9a511731 the
   command took the node's CPU count. Applied: the stub advertises `--min-exon-share --sd-damp-bed`; the helper passes
   `--threads 4` and clears the edge-rule and SEDEF variables; the expected command is the share rule.
8. **A flaky mtime test** (same file, `gtfs`). It slept 50 ms between two writes and expected different mtimes;
   `/anvme` stores coarser times, so the "older families.gtf" case passed the guard. Applied: the first file is
   back-dated with `touch -d 2000-01-01`. Check for 7 and 8: `cargo test --test pipeline_driver_units
   --test pipeline_driver_dispatch --test pipeline_driver_samples --test pipeline_driver_sd --test pipeline_seed_pool`
   all pass. The samples tests need `samtools` on PATH (pixi env).

## Round 4: per-genome defaults the run must not forget

9. **SD damping was a default only when `RUSTLE_SEDEF_BED` was exported.** A run that forgot it ran without damping:
   no error, lower precision (HG002 held-out .703 vs .545 at about equal recall). Applied: when the variable is unset,
   the families stage uses the sidecar `FASTA.sedef.bed`; `RUSTLE_SEDEF_BED=` (empty) turns it off. The sidecar's
   size and mtime enter the resume record. Sidecar placed for CHM13
   (`~/data/chm13/chm13v2.0.fa.sedef.bed` -> `sedef_dl/incoming/sedef_send/final_human.bed`). Gorilla not placed:
   `$WORK` (`/home/atuin`) refused the link with "Disk quota exceeded". Check: new test
   `a_sedef_sidecar_turns_on_sd_damping_unless_the_variable_is_set` (none, sidecar, empty variable).
10. **The censat exclusion (row 1559) is opt-in in the driver**, though memory records it as adopted for A119b. Not
    applied: it drops families from assign (it changes products), so whether it is a default is the user's decision.
    The same sidecar form (`FASTA.censat.bed`) would make it automatic.

## Round 5: `utilities` subcommands

11. Every subcommand of `utilities` was checked for callers and documentation. `sample-counts`, `parcn`,
    `window-depth`, `psv-discover`, `candidate-augment` and `sd-scan` are called by scripts or tests. `end-artifacts`,
    `gff-to-gtf`, `locus-bed` and `mcl-port` have no script caller but are documented user tools (README, REPRODUCE,
    DATA, register). Nothing deleted.

## Round 6: the whole Rust test suite (job `audit/r06`)

12. Round 3 found driver tests broken by a default flip, so every test target is run (`cargo test --release
    --no-fail-fast`). Results are recorded below when the job ends.

## Round 7: scoping the single-run layered assembly

13. The post-assembly block of `copy_assign --assemble-only` (`copy_assign_bin.rs` 10,749-11,200) works on one list
    of GTF lines: per-contig polish (the four layers act here), sub-chain, regroup, TPM, attribute passes, bridge
    regroup, then the write. The layered run differs only in four flags acting inside that block; its BAM passes are
    the same. One run could fork the list after assembly and run the block twice, writing `PREFIX.final.gtf` from the
    second pass. Gain to be priced from the e2e run (`audit/e2e0`: base and layered logs run concurrently).

## Round 8: `tools/rustle_samples.sh` (multi-sample: pool, then split)

14. **Per-sample assigns run one after another.** Each `assign` call already uses `--region-threads`; samples share no
    state but the catalog symlinks. On 2-3 samples a sequential loop is the safe default (memory: one A119b assign
    peaks at ~16 GB). Not applied: no measured need, and parallel runs would compete for the node's I/O.
15. Input checks reviewed (shared @SQ, `.bai` present, mapped reads, no shared read names, pool key, merge failure not
    keying an old pool): correct. No change.

Round 6 result: every test target passes (40 `test result: ok`, 0 failed; lib 1,203 passed, 13 ignored). This run
includes the round-3 fixes.

Round 7 result: in `audit/e2e0`, the base and layered assemblies ran concurrently. The base run took ~357 s and the
layered run ~394 s, both after a 190 s best-AS table rebuild. Merging them would save one run's CPU and BAM reads, but
little wall. **Deprioritised.**

## Round 9: a run must not read a repo file that changes under it

16. **Editing the driver during a run broke the run** (`audit/e2e0`, first attempt). Bash reads a script as it
    executes it. The assemble stage, a separate `bash "$0" assemble` call, read a file edited mid-run and died with
    "unexpected EOF" after 10 minutes. A `git pull` during a multi-hour run does the same. Applied: the driver copies
    itself (`cp -p`: same size and mtime, so the resume records' driver line is unchanged) to `$TMPDIR` and execs the
    copy.
    - Stage calls of `all` run that same copy (`RP_SNAPSHOT`).
    - A call of the repo file, from a user or `rustle_samples.sh`, makes its own copy.
    - `RP_DRIVER_DIR` keeps the repo's `tools/` for the sibling scripts and the default `--bin`. Neither variable has
      the `RUSTLE_` prefix, so neither enters the resume records.
    - The copying process removes the copy at exit (`RP_OWNER_PID`; exec keeps the pid).

    Check: the five driver test targets pass, and two calls leave no copy behind.

## Round 10: today's bench scripts (`bench/soto/rnakmer.py`, `bench/denovo_eval/exonmap.py`)

17. **`rnakmer.py arms` printed "VERDICT 1598" for every row** (row 1599's log carried it). Applied: the label is
    neutral ("VERDICT held-out"). Check: both self-tests pass.
18. **`run` and `arms` in `rnakmer.py` are near-copies** (read loop, buckets, scoring). Not merged: `run` is the frozen
    code of registered row 1596, and the counter is closed (rows 1596, 1598, 1599 FAIL), so neither will change again.

## Round 11: the families stage on the e2e run

19. `families` took 465 s on 72 threads: 401.5 s in the minimap2 all-vs-all, 64 s in everything else. The 6,445 loci
    form 1,380 families. The exact levers are spent: the gate, `--dp-memo`, the AVX2 kernel, the PAF cache for reruns.
    What remains changes products. Aligning exon sequence only, or masking sequence > 10 kb from an exon (1.14 of the
    1.86 Gbp aligned lies within 10 kb of an exon), needs a registered row priced against row 1593 on HG002. Not applied
    here.

## Round 12: where assign spends its 549 s (`audit/e2e0`)

20. **Cross-chromosome families read their copy windows one query at a time.** `~xchrom~MCL13` (150 copies on 21
    chromosomes, 147 windows, 14.8 Mb) took 149.5 s to read 39,965 reads, then 135.5 s to assign. The whole-chr1 region
    took 113.1 s to read. The window store (`RUSTLE_WINDOW_SWEEP=1`, r1507: one sweep per chromosome, 5x on the panel,
    identical tables) is opt-in since it was OOM-killed at 250 GB on the uncut chr13 hotspot (r1507, r1543 L2 withdrawn).
    The censat exclusion (row 1559) now removes that hotspot from assign, so the store may fit. Check queued (job
    `audit/r12`): assign alone with the store on the e2e0 catalog, every table compared byte for byte with e2e0's, peak
    RSS and wall recorded. A pass makes it a candidate default when an exclude BED is set; the user decides.

## Round 13: tools the default path needs

21. **A missing minimap2 surfaced only after the assembly.** `all` runs minimap2 in families (all-vs-all), assign
    (read-star), index and flag; none is checked before the 10-30 min assemble stage. Applied: `stage_all` checks
    `${RUSTLE_MINIMAP2:-minimap2}` before the plan runs and exits 2 with the variable to set. samtools and python3 are
    not needed on the default path (only `--legacy-catalog`, `--candidates`, `merged`, and the opt-in skip lists use
    them), so they are not checked. Check: the five driver test targets pass, with and without minimap2 on PATH (the
    dispatch tests' stub `all` now sets `RUSTLE_MINIMAP2=true`).
22. The index stage took 124 s (minimap2 -d, 19.4 GB) and is now built once per genome (round 2), so running it beside
    the earlier stages would save at most 124 s on a genome's first run. Not applied.

## Round 14: the multi-sample driver and the data docs

23. **`tools/rustle_samples.sh` read itself while running** (round 9's bug, second script). It runs for hours (pool,
    assemble, families, one assign per sample); an edit or `git pull` meanwhile changes the lines bash has yet to read.
    Applied: the same snapshot as the driver (`RS_SNAPSHOT`, `RS_DIR` for the sibling driver and the default `--bin`,
    copy removed at exit by the owning process). Check: `pipeline_driver_samples` passes (10 tests), and no copy is left
    in `$TMPDIR`.
24. **`docs/DATA.md` named only `RUSTLE_SEDEF_BED`** for SD damping. Applied: it names the round-4 sidecar
    (`<FASTA>.sedef.bed`, `RUSTLE_SEDEF_BED=` to turn it off) and where CHM13's is linked.

## Round 15: the assign tail

25. **The last region's read-star runs on one star slot while the node idles.** `~xchrom~MCL13` is the last region to
    finish in `audit/e2e0`; its read-star aligns 2,025 molecules to 150 copies in 123.1 s on 16 minimap2 threads
    (`RUSTLE_STAR_THREADS`), with the other star slots (72 / 16 = 4) free. minimap2's output does not depend on `-t`
    (tables identical at 8, 32 and 64 threads, `docs/PERFORMANCE.md`). Applied (`src/family.rs`, read-star call): the
    call borrows every star slot nobody holds (`try_acquire`, never waits) and runs minimap2 on `star x (1 + borrowed)`
    threads, returning them when it ends; only under the slot scheduler (`RUSTLE_FAMILY_SLOTS` not 0). A family that
    starts while the slots are lent waits for that one call. Check queued: assign alone on the e2e0 catalog with this
    binary (`audit/r15`) and with the HEAD binary (`audit/r15b`, the baseline for wall), tables compared byte for byte.

## Round 16: `merged`

26. `merged` (assemble -> families -> genome-wide assign without the prefilter) has no caller in tools, figures or
    tests besides the dispatch test, and two bench notes (`bench/MERGED_PIPELINE.md`, `bench/PERFORMANCE_AND_IO.md`).
    Deleting `stage_merged` and `stage_assign_fam` would remove ~60 lines and one of two genome-wide entry points, but
    `merged` products differ from `all`'s, so retiring it is the user's call. Not applied; `merged_env` stays either way
    (`all`'s resume records use it).

## Round 17: the flag stage's chr13 tail

27. **Pass 2 of the flag scan re-read the BAM once per firing locus** (`src/missing_copy_flag.rs`, `scan`). In
    `audit/e2e0` the other 24 contigs were done 30 min into the flag stage. chr13 then ran alone for over an hour, on
    one thread at 25% CPU: the process had read 4.4 TB from a 96 GB BAM. Every locus whose de values split (the
    mixture trigger) re-queried its own span and decoded every record in it. The chr13 hotspot's loci overlap and nest,
    so the same blocks were decompressed again and again. Applied: pass 2 is one more sweep of the contig.
    - pile()'s capped set equals pass 1's capped set: the same records are admitted (aligned_read_from_record returns
      None only for records pass 1 also skips), with the same overlap rule (pinned by the existing pile test), the same
      name-hash order and the same cap.
    - So the sweep decodes only those names, each record once, and sorts, dedups and caps as pile() did.
    - pile() is kept for its test only. The trigger is one function (`mixture_fires`) used by both passes.

    Check queued (`audit/r17`): scan.tsv and consensus.fa compared with the round-1 binary's on the 8 gorilla contigs
    and on A119b chr20-22 + chrX; chr13 alone with the new binary for its wall, and its rows compared with e2e0's once
    that run ends.

    **Result (job 4454164):** byte-identical on both panels. Gorilla: 8,076 rows, 154 consensus. A119b chr20-22 + X:
    5,742 rows, 864 consensus. A119b went 823 -> 217 s. Gorilla went 105 -> 173 s, though new ran second, on a warm
    cache: the whole-contig sweep costs a full read where few loci fire (round 18). The chr13 run did not start because
    the node has no `/usr/bin/time`; it is re-queued in round 18.

## Round 18: pass 2 reads only what fired

28. **The round-17 sweep read whole contigs** (`src/missing_copy_flag.rs`, `scan`). Where few loci fire (gorilla) that
    is a full extra pass over the BAM for little work. Applied: the fired loci, sorted by start, are merged into spans
    (gaps under 100 kb joined), and the sweep queries each span. Spans are disjoint and each holds its own fired loci, so
    every locus sees the same records as before. A read crossing two spans feeds each span's loci once.
    Check queued (`audit/r18`, job 4454266): outputs against round 17's old-binary outputs on both panels, then chr13
    alone (wall and peak RSS).

    **Result (job 4454266):** byte-identical to the round-1 binary on both panels. Gorilla 99 s (old 105, round 17 173),
    peak 2.6 GB. A119b chr20-22 + X 228 s (old 823), peak 1.7 GB. A119b chr13 alone: 452 s, peak 2.8 GB, 7,733 loci,
    1,387 fired. The old binary was still on chr13 in e2e0 after 1.5 h. The chr13 rows are compared with e2e0's when
    it ends.

## Round 19: chr13's 452 s

29. **Pass 2's span queries re-read about 8x chr13's share of the BAM** (`audit/r19`, job 4454299, sampled every 10 s).
    chr13 holds 9.4% of A119b's reads, about 9 GB of the 96 GB BAM.
    - Pass 1 (0-170 s) read about 7-9 GB at about 1.7 cores: record decoding, not I/O.
    - Pass 2 (170-410 s) read 77 GB at about 320 MB/s through 56 spans (100 kb merge gap).
    Each span query also reads the BAI's coarse bins. Those hold the hotspot's long and spliced records, so the queries
    re-read the same blocks. Experiment run (job 4454364): merge gap 100 kb vs 5 Mb vs one span from the first fired
    locus to the last. Measured on chr13, the gorilla panel and the A119b panel: wall, bytes read, and identity with
    the round-1 binary. The gap is an env var (`RUSTLE_FLAG_SPAN_GAP`) for this experiment only. The winner is then
    hardcoded and the variable removed.

    **Result (job 4454364), all byte-identical to the round-1 binary (chr13: all three equal each other):**

    | panel | gap 100 kb | gap 5 Mb | one span |
    |---|---|---|---|
    | A119b chr13 | 455 s, 84 GB read | 355 s, 16 GB | 370 s, 15 GB |
    | A119b chr20-22 + X | 229 s, 119 GB | 185 s, 29 GB | 225 s, 23 GB |
    | gorilla 8 contigs | 103 s, 26 GB | 324 s, 19 GB | 175 s, 23 GB |

    The gorilla walls swing threefold on nearly equal bytes: its BAM is on `$WORK`, and the timings there are noise.
    Round 17's "105 -> 173 s" was most likely the same noise. 5 Mb and one span read the same bytes. Applied: one span
    from the first fired locus to the last, the simplest form, and the round-18 merge loop and the env var are gone.
    Identical by construction to the measured `gap=1e15` arm (every fired locus merges into one span).

## Round 20: one thread per contig in pass 1

30. **The flag scan's long pole is one contig's single-threaded pass 1** (`src/missing_copy_flag.rs`, `scan`). The scan
    runs contigs in parallel, one thread each. A119b chr13 (9.4% of the reads) decodes for about 170 s at about 50 MB/s
    while the other threads are idle. Loci are independent: each sees the records overlapping it, and no state crosses
    loci. So a contig's position-sorted loci can be split into groups that run in parallel. Joined in order, the groups
    give the same rows and consensus in the same order. Each sweep now queries only its loci's range, not
    `chrom:1-len`; outside it, records feed no locus.
    Experiment run (job 4454467): the whole-genome A119b scan with groups off / 10 Mb / 30 Mb
    (`RUSTLE_FLAG_GROUP_MB`, for this experiment only), plus the gorilla panel at 10 Mb against the round-1 binary.

    **Result (job 4454467):** all byte-identical. Gorilla panel at 10 Mb equals the round-1 binary (8,076 rows).
    Whole-genome A119b (71,534 rows, 72 threads):

    | groups | wall | read |
    |---|---|---|
    | off (one per contig) | 509 s | 200 GB |
    | 10 Mb | 294 s | 511 GB |
    | 30 Mb | 283 s | 267 GB |

    Applied: 30 Mb groups, hardcoded, and the env var is gone. Whole-genome flag scan: hours (old binary, e2e0 still on
    chr13 after 4 h) -> 283 s. Unit tests 8/8, `tests/o2_origin_resolve` passes (with minimap2 on PATH).

    What was compared with what:
    - The three whole-genome runs were compared only with each other. All three carry the pass-1 range change, so
      that alone does not tie them to the round-1 binary.
    - Against the round-1 binary: the gorilla panel (`kb_10`). Also wg_30's chr20-22 + X rows, cut out of the
      whole-genome table: byte-identical to `r17/a_old.scan.tsv`, and all 864 of its consensus records match.
    - chr13: e2e0 hit its 6 h limit (it cannot be extended: permission denied), so there is no old-binary chr13
      table. Check queued (`audit/r20b`, job 4454622): two hotspot windows, chr13 14.0-14.25 and 15.0-15.2 Mb (1,691
      transcripts, about 436 fired), with the round-1 binary vs the final binary. Loci are independent, so their rows
      must match. The same job measures the whole-genome peak RSS (sacct is down): 72 group readers instead of 25 contig
      readers.

    e2e0 (job 4454128) ended at its 6 h limit inside the flag scan (round-1 binary, still on chr13). It started at
    07:58 and was killed at 13:39. Its earlier stages: assemble 595 s, families 465 s, assign 549 s, index 124 s. With
    rounds 17-20 the flag scan is 283 s genome-wide, so the same run would be about 34 min end to end plus the verdict
    step.

    r20b first try (job 4454622) failed twice over:
    - the new binary was named without `./`, so it never ran;
    - the round-1 binary ran 4 h on the 1,691-transcript hotspot windows and hit the job limit.
    Re-queued (job 4455435) on a smaller window, chr13 15.5-15.6 Mb (129 transcripts, about 23 fired).
    The assign checks r12 / r15b / r15 had all exited at once (127): their nodes have no `/usr/bin/time`.
    `audit/rss.py` stands in (wall and peak RSS, the same two lines). They are re-queued in a chain after r20b, so their
    walls do not share a BAM: r15b 4455436, r15 4455437, r12 4455438.

    **Result (job 4455435):**
    - chr13 15.5-15.6 Mb, round-1 binary vs final binary: byte-identical (115 rows, 23 fired, 23 consensus). Wall
      155 s vs 43 s.
    - Whole-genome A119b with the final binary: 314 s, peak RSS 16.6 GB (72 group readers). Equal to `r20/wg_30`.
    Rounds 17-20 are adopted: identity holds on every slice compared with the round-1 binary (gorilla 8 contigs,
    A119b chr20-22 + X, chr13 hotspot window).

## Round 21: PERFORMANCE.md showed the 2026-11-02 walls

31. **`docs/PERFORMANCE.md`, "Where the time goes", still gave assign as 10.75 h** and had no flag row. A reader would
    chase a stage that now takes 9 min. Applied: a current table from e2e0 and the flag scan (assemble 595, families
    465, assign 549, index 124, flag scan 283-314 s). The old table is kept under "History". Two switch rows record
    rounds 17-20.
    Held: overlapping flag with assign in `stage_all`. The gain is bounded by assign's single-region tail, which
    r15's slot borrowing targets. It would also blur the per-stage wall and RSS records and complicate `--to`/`--from`.
    Decide after r15.

## Results of the assign checks (rounds 12 and 15)

Run one after another on one node each (jobs 4455436-8). Each runs the assign stage alone on e2e0's catalog, and its 7
tables are compared with e2e0's.

| run | wall | peak RSS | MCL13 read-star | tables |
|---|---|---|---|---|
| r15b baseline (defaults) | 546 s | 35.8 GB | 116.8 s | 7/7 identical |
| r15 idle star slots borrowed | 536 s | 40.9 GB | 84.1 s | 7/7 identical |
| r12 window store (`RUSTLE_WINDOW_SWEEP=1`) | 2,914 s | 38.8 GB | — | 7/7 identical |

- **Round 15, slot borrowing:** the tail family's read-star drops 28%, but the stage gains only 10 s (inside run-to-run
  noise) for 5 GB more peak RSS. The tail is not read-star alone. Not adopted; the patch stays in `target_r15` only and
  is reverted from `src/family.rs`.
- **Round 12, window store:** the OOM is gone (38.8 GB, not 250 GB), but genome-wide it is 5.3x slower than the default
  windowed fetch now that region threads, the censat exclusion and the read budget are in. It stays opt-in.
  PERFORMANCE.md's row is updated.
- **Round 21's hold on the flag/assign overlap stands as closed:** the tail it would fill is about 2.5 min of a 9 min
  stage. Not worth the bookkeeping.

## Round 22: the families stage's minimap2

32. **Families is 465 s, and 402 s of it is one minimap2 all-vs-all** (`-x asm20 -X -t 72`, the loci FASTA against
    itself: 74,562 loci, 1.8 GB). The rest: MCL 46 s, everything else about 15 s. The question is whether that one
    process keeps 72 cores busy. Minimap2 builds the index on few threads, and the tail of the last query batch can
    leave cores idle. Measurement queued (job 4455561): the families stage alone, `--no-cache`, on e2e0's assembly,
    with minimap2's CPU seconds sampled every 5 s. Products are compared with e2e0's.

    **Measurement (job 4455561):** the families products are identical to e2e0's (clusters, copies, gated PAF), with
    the all-vs-all at 409 s. minimap2 averages 38.9 of 72 cores. Of its 405 s, 151 s run on fewer than 8 cores:
    - 0-91 s: the index build, which minimap2 runs on about 3 threads;
    - three gaps of 15-20 s (about 210, 285 and 355 s): the ends of query batches (`-K`, 500 Mb by default; 1.8 GB of
      queries makes 4 batches), each waiting for its slowest query;
    - the final ~30 s tail.
    Each query is mapped on its own, so one batch should give the same PAF without the gaps. Experiment queued
    (job 4455567): the exact command, default vs `-K 4G`, with each PAF compared with e2e0's.

    **Result (job 4455567):** the exact all-vs-all command. Default batch 399 s, 53.2 GB; `-K 4G` 326 s, 55.8 GB.
    Both PAFs are byte-identical to e2e0's. Applied in `src/family.rs`: `-K 4G` is appended after the PAF cache key is
    built. The output does not depend on it, so existing caches still hit. Check queued (job 4455577): the families
    stage end to end with the new `mcl_families` (other binaries unchanged), products against e2e0's.

    **Stage check (job 4455577):** families end to end with the new `mcl_families`: 385 s (475 s in the round-22
    measurement, 465 s in e2e0). The all-vs-all took 327 s. Clusters, copies and gated PAF are identical to e2e0's;
    `families.tsv` is not a families product, which the measurement run showed. Adopted.

## Round 23: the assemble stage

33. **Nothing to fix.**
    - e2e0's first job (4454126) failed after assemble with `line 1221: unexpected EOF while looking for matching '"'`.
      The round-13 edit of `tools/rustle_pipeline.sh` landed while that job ran from the repo file. bash reads a script
      as it goes, which is what the driver's self-snapshot now prevents (it was added after that job started). Already
      fixed.
    - Assemble's 595 s = the as_table pass, 194 s (one per BAM, cached next to it; the parallel version failed in
      r1567), plus the assembly, 398 s. The 25 per-contig `reads_in_region ... (0 reads)` lines are assemble-only
      copy_assign loading each contig, not wasted sweeps: the "0" counts assigned reads.
    - The assembly is I/O-bound and was tuned in r1564-1568. No lever found inside the audit's identity rule.

## Round 24: the flag verdict step

34. **The flag stage's second step (consensus alignment + verdicts) was never timed.** e2e0 hit its limit in the scan
    first. Measurement queued (job 4455583): `--from-scan` on the round-20 whole-genome scan, with the shared 12.3 GB
    splice index.

    **Result (job 4455583):** 41 s, peak 41.1 GB (the splice index in memory). 71,533 rows resumed, 6,302 consensus
    sequences aligned. Verdicts: scattered 3,777, reference_absent_candidate 1,207, unannotated_paralogue 975, no_home
    214, rna_editing 129. The whole flag stage is now about 6 min (scan about 300 s + 41 s). Nothing to fix.

## Round 25: what `all` needs besides its own binaries

35. **Nothing to fix.** The driver names `awk`, `sort`, `minimap2`/`minimap2-gate`, `samtools`, `python3` and `sedef`.
    On the default `all` path only minimap2 is needed; the round-13 preflight checks it. minimap2-gate is optional, a
    warning plus about 5x the families wall. The rest is opt-in:
    - samtools: the candidates stage (`--candidates`), the legacy catalog (`--legacy-catalog`) and `merged`. Default
      assign builds its regions from the copy table.
    - python3 (`annotated_families.py`): only with `RUSTLE_EXCLUDE_BED` or `RUSTLE_ANNOTATED_ONLY`.
    - sedef: only the `sd` stage.
    Queued meanwhile: `audit/e2e1` (job 4455585), `all` from scratch with every binary built from the current tree.
    Each product is compared with e2e0's, the scan with the round-20 tables, and the verdicts with round 24's.

## Check: `all` from scratch after rounds 1-25 (`audit/e2e1`, job 4455585)

Every binary built from the current tree; A119b, 72 threads, censat exclusion and SD damping as in e2e0.

| stage | e2e0 | e2e1 |
|---|---|---|
| assemble | 595 s (as_table pass included) | 397 s (as_table cached next to the BAM) |
| families | 465 s | 372 s |
| assign | 549 s | 475 s |
| index | 124 s (built) | 0 s (shared, present) |
| flag | > 5 h, killed at the 6 h limit | 327 s |
| **total** | not completed | **1,572 s (26 min)** |

All 16 compared products are byte-identical: `gtf`, `final.gtf`, `families.gtf`, the families clusters, copies and
gated PAF, all 7 assign tables (params included), the flag scan table and consensus (against `r20b/wg_new`), and
`flag.missing_copy.tsv` (against round 24). A fresh BAM adds the about-190 s as_table pass once.

## Round 26: the driver printed 0.0 GB peak memory

36. **The end-of-run table read "0.0 GB" for every stage on nodes without GNU time** (`tools/rustle_pipeline.sh`,
    `stage_all` summary). The `.time` record then has no `peak_rss_kb` line, and awk printed the empty value as 0. A
    reader would believe the stages use no memory. Applied: "peak n/a (no GNU time)" when the line is missing. Both
    cases were checked by hand; the driver dispatch tests pass (6).

## Round 27: MCL in families

37. **Nothing to fix.** MCL is 36-46 s of families. `mcl_pair_rule_counted` (`src/family.rs`) already computes the
    expansion and inflation by column in parallel, with sums in a fixed order (bit-identical).

## Round 28: where assign's 475 s go, and why the window store lost

38. **Assign (e2e1):** 416 regions on 8 region threads. Summed per region: BAM fetch 2,181 s, detect 961 s. That is about
    393 s per thread against a 475 s wall; the tail is `~xchrom~MCL13` (fetch 107 s + detect 134 s). The 231
    cross-chromosome families take 1,339 s of the fetch, re-reading blocks other regions also read (the r1505
    observation).
    **Why the window store (round 12) took 2,914 s:** since r1543's OOM fix it builds the store per batch of 64
    regions. Each batch sweeps every chromosome it touches twice (light, then full records) before its regions compute:
    7 batches, light 1,307 s + full 1,220 s, all serial ahead of the compute. One batch would sweep once (about 300 s on
    this BAM, I/O-bound) and then compute (about 140-200 s on 8 threads). That is no better than the default's 475 s,
    and it risks the r1507 OOM. Closed without a run. The store stays opt-in, and its cost model is now recorded here.

## Round 30: the round-20 query range at a contig's end

39. **Nothing to fix.** Pass 1 now queries `[first locus start, last locus end]` instead of `chrom:1-len`. A locus
    starting past the contig end (a GTF from another assembly version) gives a start beyond the contig. Tested on the
    `o2_origin_resolve` fixture (c1 is 5 kb, locus at 6,001-7,000): the e2e1 binary scans it with no reads and exits 0,
    as before. noodles accepts the range. A clamp was written and then removed: it guarded nothing.

## Round 29: assign region threads (measuring)

    Assign peaks at 35.8 GB now, not the 199 GB that set the default of 8 region threads (before the censat exclusion),
    and the BAM fetch dominates (2,181 s of fetch vs 961 s of detect). Job 4455603: 16 and 24 region threads, assign
    alone on e2e0's catalog, 7 tables against e2e0's.

    **Result (job 4455603):** 7/7 tables identical at both. 16 region threads 572 s (40.8 GB), 24 region threads 476 s
    (42.0 GB), against 8 at 475 s (e2e1) and 546 s (r15b). All inside the 475-572 s run-to-run band: the fetch is
    shared I/O on one BAM, so more readers buy nothing. Default of 8 stays. Closed.

## Round 31: `all` on a second species, with nothing extra

40. **Every check so far ran on A119b with its SEDEF and censat BEDs.** The objective is a pipeline that runs alone.
    Queued (`audit/e2e_kb`, job 4455622): `all` from scratch on gorilla KB3781 (23 GB BAM, GGO.fa) with no extra env.
    The BAM's and the FASTA's directories are on `$WORK`, whose file quota is at its hard limit, so the shared as_table
    and the shared splice index cannot be written there. The run shows whether both fall back cleanly.

    **Result (job 4455622): passes, 756 s end to end, exit 0.** assemble 299 s, families 133 s, assign 127 s, index
    134 s, flag 61 s.
    - The as_table went to the run directory (`kb.molecules.tsv`, one BAM pass) because the BAM's directory takes no
      file.
    - The splice index went to `kb.splice.mmi` (13.6 GB) because the FASTA's directory takes no file.
    - With no SEDEF sidecar, families ran undamped, as documented.
    - Verdicts: scattered 248, reference_absent_candidate 84, no_home 2, contamination 1, rna_editing 1.
    - The end-of-run table now reads "peak n/a (no GNU time)" (round 26).
    One cost to know: on such a layout every run rebuilds the 134 s index and the as_table.
    Both have a switch: `RUSTLE_INDEX_DIR=<writable dir>` shares the index, and `--as-table PATH` shares the table (its
    header names the BAM's path and size, and a table from another BAM is rebuilt, round 32).

## Round 32: sharing the as_table when the BAM's directory is read-only

41. **Nothing to fix; round 31's note corrected.** `--as-table PATH` already exists (`stage_assemble`). The path is used
    when its first line names this BAM's real path and size and it is not older than the BAM (`as_table_ok`).
    Otherwise it is rebuilt there under a temporary name and renamed. Together with `RUSTLE_INDEX_DIR` this covers
    read-only input directories.

## Round 33: `--low-mem` after the speed rounds

42. **`--low-mem` had not run since rounds 17-22.** That profile is 1 region thread, the FASTA indexed on the fly in
    1-Gbp parts and `--max-window-records 500000`. Queued (`audit/e2e_kb_lowmem`, job 4455642): `all --low-mem` on
    gorilla, products compared with e2e_kb's. By design, flag may differ (split-index hits) and assign may differ
    (windows above the record cap are skipped). Assemble and families must match.

    **Stopped by the user during assign (job 4455642 cancelled).** Partial result: assemble 637 s (1 region thread),
    families 135 s. `gtf`, `final.gtf`, `families.gtf`, the families clusters and copies are byte-identical to e2e_kb's.
    The assign and flag differences under `--low-mem` were not measured.

## Round 34: the README could not get a new user to a run

43. **`README.md` had no run command and no word on minimap2.** "Build, test, run" listed `cargo build` and the tests
    only. A new user does not learn that `all` needs `minimap2` on PATH, which the preflight refuses without. Nor that
    families is about 5x slower without `minimap2-gate` (`tools/minimap2_gate/build.sh`), or how the BAM must be
    aligned. Applied: a short block with the gate build, one `all` command, the inputs (indexed BAM aligned as in
    `docs/DATA.md`, FASTA with `.fai`), the two shared sidecars with their switches (`--as-table`,
    `RUSTLE_INDEX_DIR`), and the A119b wall. No test reads the README.

## Round 35: stale wall figures in the docs

44. README, REPRODUCE, DATA and MODULE_STATUS quote no stale stage walls. The open preregistrations carry their own
    dated estimates, as they should. `docs/PERFORMANCE.md`'s current table lacked the verdict step and the end-to-end
    total. Applied: rows for the verdict step (41 s, 41 GB) and `all` (1,572 s, e2e1; gorilla 756 s).

## Round 36: `pile()` is test-only

45. **Kept, with the reason written down.** Since round 17, `pile()` (`src/missing_copy_flag.rs`) runs only in its
    test, which could read as dead code. The test still pins what rounds 18-20 rely on: a BAI region query returns the
    records that start before the range but reach into it (its spliced read, 0-based [500,1100), against a locus at
    1,000). That is why both sweeps can query [first locus, last locus] instead of the whole contig. Its comment now
    says so. Unit tests 8/8.

## Round 37: `--threads` left out ran four stages' steps on 4 threads

46. **Without `--threads`, only families used the job's CPUs** (commit 9a511731). as_table, both flag steps and the
    index build ran on the default 4. The flag scan then keeps 4 of its 30 Mb groups going at once, about 4x its
    72-thread wall. A user who follows a README without `--threads` gets a slow run for no reason. Applied
    (`tools/rustle_pipeline.sh`): one `JOB_THREADS`, set after argument parsing: `--threads` if given, else
    `SLURM_CPUS_PER_TASK` or `nproc`. Families now uses it too, replacing its own copy of the rule. The settings
    records keep `$THREADS`, so a resumed run does not rerun. Left alone: catalog and candidates (opt-in, their thread
    dependence unchecked) and `merged`'s region threads (memory). Driver tests pass (units 4, dispatch 6, samples 10).
    Check queued (job 4455689): as_table, index, flag scan and verdict at 4 threads, against e2e_kb's 72-thread products.

    **Result (job 4455689): all identical at 4 threads to e2e_kb's 72-thread products**: as_table, splice index, flag
    scan table and consensus, verdicts. Walls at 4 threads: as_table 595 s (72 threads: inside assemble's 299 s), index
    158 s (134 s), scan 344 s (flag stage 61 s in all at 72), verdict 30 s. Adopted: a run without `--threads` on a
    72-core gorilla job saves about 10 min.

## Round 38: the index waits behind assign for no reason

47. **The index needs only the FASTA, yet `all` builds it after assign** (`stage_all`). On a genome's first run, or
    whenever the FASTA's directory takes no file, that is about 130 s of serial wall (gorilla 134 s of 756 s). Applied:
    when the plan has an index stage, no `--from`/`--to` is given, and the index is missing or its key is stale, `all`
    builds it in the background from the start. The build uses the same timing wrapper and records. The index step
    waits for it (before the resume check, so a half-built index is never skipped past) and replays its log. The build
    still writes a temporary name and renames, so a failed earlier stage leaves a finished index or none. Driver tests
    pass. They pass `--index`, so this path is checked by a run: job 4455730, gorilla `all` with a fresh
    `RUSTLE_INDEX_DIR` and `--as-table` on e2e_kb's table (which also exercises that switch). Products and the index
    are compared with e2e_kb's.

    **Result (job 4455730): passes.** 459 s end to end, exit 0. All 9 compared products plus the index are
    byte-identical to e2e_kb's. Stage walls: assemble 136 s (table reused through `--as-table`, so that switch works),
    families 134 s, assign 122 s, index 144 s (background, under assemble; the wait found it done), flag 67 s. Against
    e2e_kb's 756 s: about 160 s comes from the reused table and 134 s from the overlapped index. A new genome's first
    run now pays the index under assemble instead of after assign. Adopted.

## Round 39: the flag scan could start after assemble

48. **The scan reads only the BAM and `families.gtf`, an assemble product, yet it waits for families and assign.**
    Overlapping it with them (about 850 s on A119b) could take up to about 290 s off. But assign's fetch also reads
    the BAM, and families' minimap2 uses every core, so contention may eat the gain. Measured before any driver code
    (job 4455735): families (no cache) and then assign on e2e1's assembly, with the scan running alongside from the
    start. Walls are compared with e2e1's serial 372 + 475 + about 290 s, and products with e2e1's.

    **Measurement (job 4455735):** with the scan alongside from the start, the scan took 302 s (about 290 alone),
    families 419 s (372), assign 471 s (475). Together 890 s against about 1,137 s serial, so about 247 s saved.
    fam.copies, assignments, quant, scan table and consensus are identical to e2e1's.
    **Applied (`tools/rustle_pipeline.sh`):**
    - The scan is its own internal step, `stage_flag_scan` (dispatch `flag_scan`).
    - In `all`, once assemble has run (so flag will run) and with no `--from`/`--to`, it starts in the background.
    - The flag step waits for it, fails if it failed, then runs flag with `RP_FLAG_SCAN_DONE=1`, which skips the scan.
    - A direct `flag` call still scans first.
    `tests/pipeline_driver_dispatch.rs`: the plan's order is checked without the `flag_scan` line, which may land
    anywhere after assemble. A new test asserts it runs once, after assemble and before flag. 7 pass, 3 runs in a row.
    Check queued (`audit/e2e2`, job 4455739): `all` from scratch on A119b, all 16 products against e2e0's.

    **Result (`audit/e2e2`, job 4455739): passes.** `all` from scratch on A119b took 1,328 s (22.1 min), against
    1,572 s in e2e1. All 16 products are byte-identical. Stage walls: assemble 401 s, families 380 s (372 in e2e1),
    assign 502 s (475), index 0 s, flag 43 s (the verdict step; the scan ran under families and assign). Adopted.

## Round 40: `--low-mem` and the background scan

49. **The round-39 background scan ignored `--low-mem`.** That profile bounds memory: 1 region thread and the index in
    1-Gbp parts. Yet the scan (A119b 16.6 GB) would have run beside families' all-vs-all (about 55 GB). The index
    background build was already off there, since `--low-mem` plans no index stage. Applied: the scan launch requires
    `LOW_MEM=0`. Dispatch tests pass. Origin checked (`git fetch`): not moved, main is 5 commits ahead.
    Queued (job 4455755): `all --low-mem` on gorilla to completion (round 33 was cancelled at assign), with the node's
    memory in use sampled every 10 s.

## Round 41: a failed run left its background scan running

50. **When a stage failed, `all` exited but its background scan ran on.** The scan writes `PREFIX.flag_scan.*` in place,
    and a rerun that skips assemble runs flag's own scan into the same files. Under Slurm the job end kills it;
    interactively nothing does. Applied: on a stage failure, and on a failed background index, `all` stops the scan
    (`pkill -P` on its wrapper's PID, then the wrapper; no name matching). A background index build is left to
    finish: it writes a temporary name and renames, and the result is reusable. The first version broke an existing
    test: under `set -e`, `pkill` with nothing left to kill returned 1, and `all` exited 1 instead of the stage's
    status. Fixed with `|| true`. New test `a_failed_stage_stops_the_background_flag_scan`: the stub scan sleeps 3 s and
    then logs `scan_end`, assign fails, and 5 s later `scan_end` must be absent. Dispatch 8/8 (3 runs), units 4,
    samples 10.

## Round 42: "22 min" holds from scratch only

51. A resumed run that skips assemble, or a `--low-mem` run, scans after assign as before (about +290 s on A119b).
    That is correct, but a reader quoting 22 min would file it as a regression. README and PERFORMANCE.md now say
    "from scratch" and name both cases.

## Round 43: the round-40 low-mem run, and what bounds families' memory

52. **Result (`audit/r40`, job 4455755): `all --low-mem` on gorilla runs to the end.** 1,549 s, exit 0. Assemble,
    families, assign and the scan are byte-identical to the normal run (`audit/e2e_kb`). `flag.missing_copy.tsv`
    differs in 7 of 16,467 rows, with no verdict class changed: the split index (1-Gbp parts) places a few
    consensus sequences differently, by design. Node memory in use peaked at 43.7 GB at 651 s, inside families
    (611 s 37 GB, 641 s 41 GB, 671 s 44 GB). `--low-mem` promises a split index and one assign region; it says nothing
    about families, and families sets the peak. Queued (job 4455848): families' all-vs-all on the gorilla loci at
    threads 72/16 and `-K` 4G/500M, peak RSS each, PAF must stay identical. If fewer threads or a smaller batch
    cuts the peak without changing the PAF, `--low-mem` gets it; otherwise the doc names families as the peak.

    **Result (`audit/r43`, job 4455848):** families' all-vs-all on the gorilla loci, peak RSS by arm. Threads set the
    peak, `-K` barely does. The PAF is byte-identical in every arm:

    | threads | `-K` | wall | peak RSS |
    |---|---|---|---|
    | 72 | 4G | 127 s | 38.3 GB |
    | 16 | 4G | 219 s | 17.7 GB |
    | 72 | 500M | 143 s | 36.4 GB |
    | 8 | 500M | 391 s | 14.2 GB |

    Applied: under `--low-mem`, families runs on at most 16 threads, and the driver's `--low-mem` comment says so.
    That halves families' peak for +92 s on gorilla. Threads are part of the families PAF cache key, so a `--low-mem`
    run does not reuse a 72-thread run's PAF: it aligns once more. That is acceptable for a memory profile. New test
    `low_mem_caps_family_threads_at_16` (stub `mcl_families`, records `--threads`): it fails with the cap removed and
    passes with it. Driver units 5/5, dispatch 8/8. PERFORMANCE.md gains a `--low-mem` row.

## Round 44: the low-mem run with families capped

53. **Result (`audit/r44`, job 4455851): `all --low-mem` on gorilla with families at 16 threads.** 1,648 s, exit 0
    (r40: 1,549 s). Families took 224 s (134). All 9 compared products are byte-identical to r40, flag verdicts
    included. Node memory in use during families is now about 25 GB (r40: 44 GB). The run's peak moved to flag's
    align step: 39.7 GB at 1,563 s, with a 6.4 GB baseline. Queued (job 4455947): flag's minimap2 call (`-I 1G`,
    on-the-fly index parts) alone at 72/16/4 threads, peak RSS each, PAF compared.

## Round 45: flag's align under `--low-mem`

54. **Flag's minimap2 call (`-I 1G` on-the-fly parts) scales with threads too** (`audit/r45`, job 4455947): 72 threads
    171 s / 43.6 GB, 16 threads 203 s / 25.4 GB, 4 threads 317 s / 18.8 GB, PAF byte-identical in all. Applied: one
    driver variable, `MM_THREADS`, caps families' all-vs-all and flag's align at 16 threads under `--low-mem` (replacing
    round 43's families-only line). Verified (job 4455993): `flag --low-mem` on a copy of r44 wrote a verdict table
    byte-identical to r44's, peak RSS 24.9 GB. A first attempt (job 4455987) exited 2 because a named `flag` needs
    `--index`; its "same" compared the untouched copy and was discarded. Driver tests: dispatch 8, samples 10, units 5.

# Second goal (user, 2026-11-10): 15 more rounds (46-60), every finding applied

Goal: the pipeline much faster, self-contained, and as close as possible to O1 and O2. Unlike the first goal there is no
"last rounds need not be applied" exemption. Limits kept: byte-identical products unless a round fixes a bug (then a failing
test first and a before/after diff); no rule is adopted and no default flipped under the audit banner (registered questions
stay registered); deletions of opt-in hatches are collected for the user. Before round 46 the uncommitted tree was saved
(`audit/g2/pre/`: `tracked.patch`, `untracked.tgz`), and `all` was run from scratch on A119b with the tree's binaries frozen
in `audit/g2/bin0` (`audit/g2/e2e_base`, job 4459898) as the reference for the rounds below.

## Round 46: the guided assembly (today's port, uncommitted)

55. **`utilities guide` read the BAM on one thread** (`guide_run`): 1,469 s on gorilla OR6737, longer than the whole
    assemble stage. The pass only counts full-chain primary alignments per guide chain, per contig. Applied: one reader
    per guide contig through the `.bai` (`bam::seek_to_contig`), contigs in parallel on `--threads` (new option, the driver
    passes its own). Check (`audit/g2/r46`, job 4459906): 75 s at 72 threads, `guided.gtf` byte-identical to the
    single-threaded port's (47,369 multi-exon chains, 7,349 guide transcripts added).
56. **The command had no test.** New `tests/utilities_guide.rs`: a two-contig BAM; an assembled guide chain is not
    duplicated, a chain with two primary full-chain alignments is rescued with `reads "2"`, a chain with only a secondary
    alignment and one with no read are not written, and the output is the same at 1 and 4 threads. It failed first on
    the missing `--threads`.
57. **A resumed `all --guide` that skipped assemble wrote no guided file**, and neither did a guide added to a finished
    run: the step lived inside `stage_assemble` and the guide is not in assemble's settings record. Applied: `guided` is
    its own stage (`stage_guided`; `rustle_pipeline.sh guided --guide G ...` works alone). `all` calls it once assemble is
    behind it, run or skipped, when `PREFIX.guided.gtf` is missing or older than `PREFIX.gtf` or the guide. It sits outside
    the stage plan, so writing it never makes families or assign rerun. A direct `assemble --guide` still writes it.
    New test `all_with_a_guide_writes_the_guided_assembly_even_when_assemble_is_skipped` (failed first: no `guided` call;
    then: full run, rerun does nothing, product removed -> `guided` alone). Driver tests: dispatch 9, samples 10, sd 10, units 5.
58. **Mistake made and undone in this round:** the first edit of `utilities.rs` anchored on a line that also opens an
    earlier function and cut 1,640 lines; the file was restored from the round's snapshot and the edit redone inside
    `guide_run`. Edits to the large files are now anchored after the function's own `fn` line.

Looked at and left: **constant columns in the default tables** (e2e2, A119b). `as_margin` 0 (assign admits AS ties only),
`aligner_disagreement` 0 (opt-in gate), `depth_cn` NA (no DNA), `run_p` 1 / `run_top` empty (one sequencing run in the BAM,
`n_runs` = 1 at all 6,302 fired loci), the legacy-catalog columns of `fam.copies.tsv` (`core_hull`, `sd_depth`, ...): all as
designed. `junction_only` is 0 in all 752 families and `resolvable_j` = `resolvable_psv` = 70,055: no read is decided by a
copy-specific junction alone. That agrees with registers 679 and 681 (tied reads do not carry copy-specific junctions) and
no failing case was found, so it is recorded, not changed.

## Round 47: what the layered second assembly costs (measurement; not changed)

58b. **`PREFIX.final.gtf` is a whole second assemble run beside the first** (registers 1490/1491 keep it out of the
    families input). Measured twice (`audit/g2/r47`, job 4459907, frozen binaries): layers off 345 and 339 s, layers on
    415 and 414 s, so the second run costs about 72 s of wall although it repeats all the work. It is not thread
    oversubscription (`audit/g2/r47b`, job 4459999: 72, 36 and 24 rayon threads per run give 399, 397 and 413 s, both GTFs
    identical each time). The only way to save it is one process that forks each contig's transcripts after the shared
    polish and runs the end layers and the whole emit tail (sub-chain, TPM, tags, bridge regroup, side tables) twice.
    That tail is about 450 lines holding ten pieces of mutable state for two products that must both stay
    byte-identical. NOT APPLIED: the cost is 5% of an `all` run and the change is the riskiest in this goal. It is on
    the list for the user; `RUSTLE_ASSEMBLY_LAYERS=off` already skips the cost when `final.gtf` is not wanted.

## Round 48: where a BAM pass spends its time (measurement, one doc correction)

59. **The note "BAM passes are I/O-bound on /anvme (~610 MB/s one stream)" was wrong in its cause** (`audit/g2/r48`).
    Measured on the A119b BAM: the file system gives about 1.0 to 1.07 GB/s in total whatever the number of streams, a
    warm page cache 4 to 5 GB/s; yet the `as_table` pass takes 170 to 181 s on a warm BAM against 194 to 206 s cold. A
    cached file saves about 25 s of 200: the pass is bound by inflating and decoding records, not by the disk. Loading
    the `.bai` costs 0.05 s warm (2.2 s cold). Consequence for the rounds below: more parallel streams do not help,
    cheaper work per record does. A first version of the job measured nothing (the stream offset expanded to 0 in the
    outer shell, so every "parallel" stream re-read cached bytes); it was rerun from a quoted script and only the rerun
    is quoted here. Applied: the memory note and `docs/PERFORMANCE.md` wording.

## Round 49: compiler warnings

60. **36 warnings hid the ones that matter.** `cargo fix` applied 22 (unused imports, unneeded `mut`) in `family.rs`,
    `copy_assign_bin.rs`, `mcl_families.rs`, `sedef/seed.rs`, `sedef/decompose.rs`, `copy_phasing/tied_reads.rs`. It also
    removed three imports that only tests use (`HashSet` in `family.rs`; `gtf_loci`, `junction_count`,
    `fold_parts_into_loci` in `mcl_families.rs`), which broke the test build; they are back under `#[cfg(test)]`. 14 are
    left and were read one by one: a timing macro's last assignment, three functions only tests call, unread fields of
    `family_gfa`. None is a bug. No product can change (no executable statement was touched); the full suite and the
    byte-identity run of round 51 cover it.

## Round 52: what the driver needs from outside itself

61. **A missing `python3` stopped `all` after assemble and families, with no message.** Assign's two skip lists
    (`RUSTLE_EXCLUDE_BED`, the adopted censat exclusion, and `RUSTLE_ANNOTATED_ONLY`) are built by
    `tools/annotated_families.py`; its failure was `|| exit 2` with the reason only in `PREFIX.assign.depth.log`.
    Applied: `all` checks for the interpreter at the start, next to the minimap2 check, only when a skip list is set, and
    names the variable; a direct `assign` now says the skip list could not be built and where to read why. `RUSTLE_PYTHON`
    names another interpreter (as `RUSTLE_MINIMAP2` does). New test
    `all_refuses_at_the_start_when_a_skip_list_needs_a_python_that_is_not_there` (failed first: exit 0 and assemble ran).
    Looked at and left: `samtools` is called only by the opt-in `candidates`, `--legacy-catalog` and `merged` paths; the
    default `all` needs minimap2, the Rust binaries, and python3 only with a skip list.

## Round 50: inflate workers per region reader

62. **Each region reader inflated BGZF blocks on at most 4 threads, whatever the machine.** An old comment said workers
    do not matter; that was measured when every record was converted to `RecordBuf` (the consumer was the bound).
    Assemble's sweep reads lazy records now, and there inflation is the bound (`audit/g2/r50`, job 4459935, A119b, 72
    cores): assemble 397 s at 4 workers, 309 s at 8, 275 s at 16, `gtf` and `final.gtf` identical to the baseline each
    time. Assign, still on `RecordBuf` in those binaries, did not move (477, 471, 464 s; 6/6 tables identical): round 51.
    Applied: one helper, `bam::bgzf_workers` (the two copies of the default are gone): `RUSTLE_BAM_THREADS`, else a
    quarter of the cores, at least min(4, cores), at most 16. New test
    `region_reader_workers_scale_with_the_machine_and_obey_the_override` (failed first: no such function).

## Round 51: region readers build their reads from the lazy record

63. **Every region reader converted each BAM record to a `RecordBuf` first** (all data fields parsed, the sequence decoded
    base by base), then copied what it wanted; the tied prefilter's `Light` pass decoded sequences only to drop them.
    Measured on A119b chr20 (297k records, 16 inflate workers): 3.6 s with the conversion, 2.2 s reading the same fields
    from the lazy record. Applied: one routine, `lazy_reads`, gives the primary-pool read and the assignment read for a
    record in any `ReadMode`; the plain reader, the cached reader (`BamIndexCache`, the one assign uses) and the window
    sweep call it when no read-home table is set. The `RecordBuf` code stays for the home table and the secondary AS
    ratio. The test `streaming_and_buffered_readers_apply_the_same_home_rule` now compares the two paths field by field
    (name, CIGAR, sequence, qualities, MAPQ, AS, ms, de, ts, flags) on the plain reader and on the cached reader in
    `Full`, `Light` and `Only`; a wrong AS, a wrong ms and a wrong light CIGAR were each planted and each failed it.
    Full suite 1,540 passed, 0 failed.
64. **Mistake in this round, kept on record:** the first version changed only the plain reader. Its check (`audit/g2/r51`,
    job 4459998) gave 6/6 identical tables and no speed-up, because assign reads through the cached reader, a second
    copy of the same loop. Finding 63 is the corrected version (`audit/g2/r51b`, job 4460055).
65. **Result for assign: identical, not faster.** 6/6 tables identical at the default, at 4 workers and with
    `RUSTLE_WINDOW_SWEEP=1`; wall 488, 450 and 1,319 s against the baseline's 472 to 477 s. Summed fetch time is unchanged
    (2,123 s against 2,105 s over 416 regions). So assign's fetch is not record decoding: a region's copy windows are
    inflated on one thread, twice (light pass, then the tied molecules' full records), and 8 regions run at once on a
    72-core node. The sweep stays opt-in (it reads each batch's hull). Round 54 measures the region concurrency.

## Round 53: checkpoint run

`all` from scratch on A119b with rounds 46 to 52 applied (`audit/g2/r53`, job 4460019): **1,193 s against 1,315 s**,
all 16 products byte-identical to `e2e_base`. Assemble 283 s (was 422), families 387, assign 479, index 0, flag 42.

## Round 54: how many regions assign reads at once

66. **Assign ran 8 regions under a 4M-record budget on any machine.** Those numbers were set (r1543) when the satellite
    families were still assigned and memory peaked at 199 GB. With the adopted exclusion the regions are light, each is
    inflated on one thread, and a 72-core node sat mostly idle during the fetch. Measured on A119b, assign alone, six
    tables identical to the baseline in every arm (`audit/g2/r54`, `r54b`, jobs 4460095 and 4460128; peak RSS sampled from
    /proc): 8 / 4M 484 s, 34 GB; 16 / 4M 444 s; 24 / 4M 435 s; 16 / 8M 349 s, 40 GB; 24 / 16M 286 s, 40 GB; 36 / 32M
    281 s, 43 GB. The budget is the gate, and the gain stops at 24 / 16M. Applied: the driver picks 24 / 16M only on
    the measured configuration (>= 48 cores, >= 128 GB, `RUSTLE_EXCLUDE_BED` set, not `--low-mem`), else 8 / 4M as before;
    `RUSTLE_REGION_THREADS` and `RUSTLE_REGION_READ_BUDGET` still override. Assemble keeps its own 8 regions. New test
    `assign_reads_more_regions_at_once_only_on_the_measured_configuration` (failed first: no such function).
    Not done: reading each window once instead of twice. It needs each molecule's genome-wide placement count, which
    the binary side table does not carry, or full records held for a whole region; with the fetch now spread over 24
    regions the tail is one family (`~xchrom~MCL13`, about 230 s of fetch plus assignment), which that would not shorten much.

## Gorilla checkpoint (rounds 46 to 53)

`all` from scratch on gorilla KB3781, no extra environment (`audit/g2/kb_base` with the pre-goal binaries, `kb_r55` with
the current ones, jobs 4460129 and 4460130): 594 s against 561 s, all 16 products byte-identical. A new genome's index
build (142 s, in the background) is in both.

## Round 55: the flag stage no longer waits for assign

67. **`flag`'s align + verdict ran after assign although it reads nothing of families or assign** (the scan, the loci of
    assemble, the index). Applied: the background job `all` starts at the end of assemble is now the whole stage
    (`flag_bg` = scan, then align + verdict; its wall goes to `PREFIX.flag.time`), and the `flag` step of the plan only
    collects it. Two cases keep the old order, scan early and verdict last: `--candidates` (the verdict reads the
    candidates table) and a genome's first run (the index is still being built in the background). `--low-mem`,
    `--from`, `--to` and a resumed run that skips assemble are unchanged (no background job). Six stage-order tests were
    rewritten first and failed; `all_runs_the_flag_scan_once_in_the_background_after_assemble` also pins the
    `--candidates` case. Driver tests: dispatch 11, samples 10, sd 10, units 5.

## Round 56: checkpoint run

`all` from scratch on A119b, rounds 46 to 55 applied (`audit/g2/r56`, job 4460189): **964 s against 1,315 s** (-27%),
all 16 products byte-identical to `e2e_base`. Assemble 288 s (was 422), families 399 (378), assign 277 (472), flag 60 in
the background (43, after assign). The run before it, without round 55 (`audit/g2/r55`, job 4460133), took 1,126 s with
the same 16 identical products; its families stage took 517 s on that node, against 378 to 399 s in the four other
runs of this goal, so stage times between nodes differ by more than round 55's gain and only the sum of several runs
should be compared.

## Round 57: gorilla with the final driver

`all` from scratch on gorilla KB3781 (`audit/g2/kb_r57`, job 4460196): 546 s against 594 s with the pre-goal binaries,
16/16 products byte-identical. This run builds its index in the background, so it covers the branch of round 55 that
keeps the verdict last.

## Round 58: assemble's time after the sweep, and contigs at once

68. **Assemble swept 8 contigs under a 20M-record budget on any machine**, numbers from before the sweep read lazy
    records. Measured on A119b, assemble alone, `gtf`, `final.gtf` and `families.gtf` identical to the baseline in every
    arm (`audit/g2/r58`, job 4460197): 8 / 20M 291 s, 9 GB per process; 16 / 40M 237 s, 11 GB; 25 / 80M 212 s, 13 GB
    (two processes run side by side). Applied: 24 / 80M on a large node (>= 48 cores, >= 128 GB, not `--low-mem`), else
    8 / 20M; `RUSTLE_REGION_THREADS` still overrides. New test
    `assemble_sweeps_more_contigs_at_once_only_on_a_large_node` (failed first: no such function).
69. **Where the rest goes (timeline of the 8 / 20M arm, log sampled every 5 s):** the sweep ends near 130 s; the
    per-contig polish, one contig after another, runs to about 225 s; the base run has written everything by 242 s; the
    stage ends at 291 s, when the layered second run finishes. NOT APPLIED, for the user's list: the per-contig
    polish loop is serial (about 85 s) and could run contigs in parallel, and the layered run's tail (about 50 s) is the
    same cost as round 47. Both sit in the emit code of finding 58b.

## Round 59: final runs

`all` from scratch with rounds 46 to 58 applied, binaries equal to the tree's build (cmp), full suite 1,542 passed, 0
failed; driver tests dispatch 12, samples 10, sd 10, units 5.

| sample | before the goal | after | products |
|---|---|---|---|
| A119b on CHM13 (`audit/g2/r59`, job 4460207) | 1,315 s (assemble 422, families 378, assign 472, flag 43) | **918 s** (assemble 243, families 404, assign 270, flag 64 in the background) | 16/16 byte-identical to `e2e_base` |
| gorilla KB3781 (`audit/g2/kb_r59`, job 4460208) | 594 s (assemble 296, families 130, assign 141, index 145, flag 27) | **510 s** (assemble 244, families 127, assign 111, index 141, flag 27) | 16/16 byte-identical to `kb_base`; peak `copy_assign` RSS 10 GB |

Families is now the longest stage on A119b and was not changed: 98% of it is the minimap2 all-vs-all, whose
registered lever (exon-only alignment) is a rule question, not an audit one.

## Round 60: closing review of the goal's own changes

The diff of rounds 46 to 58 was read again against the code it replaced.
70. A comment in `stage_assemble` still gave the old 20M budget as the only value: corrected.
71. Nothing else found. Checked on purpose: `lazy_reads` against `alignment_read_from_record`, `aligned_read_from_record`
    and `bam_read_for` field by field (unmapped and start-less records, the `Only` name test, the light CIGAR, `ms`
    only on the cached readers, missing qualities are the same raw bytes on both paths); the background flag job on a
    failed stage (still one process deep, so the existing kill reaches the binary), on resume (a skipped assemble
    starts no background job and flag runs in the foreground as before) and with `--no-flag`, `--from`, `--to`.
    One limit kept as it is: the machine rules read the node's total memory, not a job's memory limit; on a node with
    48 or more cores and a job limit under about 48 GB set `RUSTLE_REGION_THREADS=8` (or use `--low-mem`).

72. **(advisor review) A failed stage left `minimap2` running under the background flag job.** `pkill -P` reached the
    job's children only; `missing_copy_flag` starts `minimap2` one level further down. `kill_tree` in the driver walks
    the whole process tree. The kill test got a grandchild first and failed on it (`grand_end`), then passed.
73. The driver exports `RUSTLE_BAM_THREADS` from the job's thread count when it is unset, so the inflate workers of
    round 50 follow `--threads` and not the node's core count.

## Round 61: the per-contig polish in parallel (finding 69, applied)

`copy_assign_bin.rs`: the loop that polished one contig after another is now a parallel map over contigs. Each
contig gets its own slice of the TSS evidence, buffers its log lines, and returns its models, floors, attributes and
counts; the merge runs in contig order, so the GTF and the log come out in the old order.
- Proof (`audit/g2/r62`, job 4460232, A119b assemble twice): `gtf`, `final.gtf`, `families.gtf`, `bridges.tsv`,
  `bridge_junctions.tsv`, `params.tsv` byte-identical to `e2e_base` both times; both logs hold the same lines in the
  same order (the only differing line is one that prints a load time). Full suite 1,542 passed, 0 failed.
- Assemble: 243 s -> 171 s and 140 s (two runs on one node, the second with a warm cache).
- The other half of finding 58b (one process writing both GTFs) stays not applied.

## Round 62: do the tables of a run agree with each other (O1 and O2 outputs)

`tools/check_run_tables.py PREFIX` reads the family and assignment tables of one run and tests 21 statements: cluster
sizes equal their member rows, a member locus is in one cluster, every family has two or more copies numbered from 0
without gaps, exons are sorted, disjoint and span the copy, `fam.copies.fa` has one record per copy row in the same
order, every assigned copy exists, a read has one row per family, assigned rows are AS ties with at least one
candidate and are never orphans, `n_reads_hard` equals the rows whose best copy it is (orphans out, as the code
comment above the quant write says), `families.tsv` read counts equal the assignment rows.
74. **All 21 hold on A119b (`r59`) and on gorilla (`kb_r59`). No bug found.** Two statements in my first version
    failed and were mine: I assumed the hard count used assigned rows only and that an assigned read always has two
    candidates. Both are documented otherwise in the code (`sole_candidate` column; the quant scope note).
75. **Seen, not changed (rule questions):** on A119b 946 of 3,698 assigned molecules are assigned in more than one
    family (gorilla: 0 of 722). The run reports this itself and `RUSTLE_XFAM_RECONCILE` is the registered switch.
76. **Copy pairs of different families that share sequence (A119b, 645 pairs among the 2,582 assigned-family
    copies):** 285 are on opposite strands; of the 360 on the same strand 141 share no exon base, 145 share exons over
    less than half of the shorter copy, and 74 share exons over half or more (example: `DN_chr1_16302609_2` in MCL142
    with 24 reads against `DN_chr1_16302612_2` in MCL114 with 2 reads). `catalog_overlaps` compares spans without
    strand or exons, and its comment called every such pair a defect. The comment is corrected; the rule is not
    changed, because the pair list feeds the `shared_locus` stratum and the 74 are the open locus-rule question.
77. **Missing test added:** `lazy_and_buffered_readers_agree_on_clipped_reverse_tagged_records`. The round-51 test
    used a fixture with no soft clips, no reverse strand, no sequence and no tags. The new one has all of them plus a
    secondary and a supplementary record, in the three read modes. Three deliberate breaks of `lazy_reads` (strand,
    qualities, `de`) were run against it and each one failed it.

## Round 63: does `all` need anything it does not say

`all` on `test_data/synthetic_family` under `env -i` with `PATH=/usr/bin:/bin`: stops at once with "no minimap2 on
PATH" (exit 2). With `minimap2` added and nothing else (no pixi environment, no python, empty `HOME`): exit 0, every
product written, with the known warning that `minimap2-gate` is absent (findings 2 and 35: optional, families
about 5x slower without it; the 403 s figure is without it too, it is not installed here). Limit: that fixture has no family, so assign's read work was not run in the empty environment; the
binaries start no program other than `minimap2` (and `mmseqs` when asked), read in round 52.

## Round 64: closing runs with everything applied

`all` from scratch, driver of this tree, binaries `audit/g2/r62/bin` (the tree's code; the tree's later build differs
only by a comment and a test). Single runs; the same configuration varies by about 100 s between nodes.

| sample | before the goal | after | products | memory |
|---|---|---|---|---|
| A119b on CHM13 (`r63`, job 4460234) | 1,315 s | **886 s** (assemble 169, families 403, assign 312, flag 78 in the background) | 16/16 byte-identical to `e2e_base`; 21/21 table checks | see round 65 |
| gorilla KB3781 (`kb_r63`, job 4460235) | 594 s | **488 s** (assemble 215, families 135, assign 111, index 143, flag 27) | 16/16 byte-identical to `kb_base`; 21/21 table checks | not measured (sampler bug, round 65) |

`r61` (job 4460223, same driver, binaries before the parallel polish): 911 s, 16/16. End to end the four A119b runs
of the final driver gave 964, 918, 911, 886 and 851 s: the last step is inside node noise; the stage figure is the claim
(assemble 242-243 s -> 140, 169, 171 s).

## Round 65: closing review (advisor), memory measured again, opt-in modes

78. **My memory sampler was wrong.** It added each process's resident size to the next process's owner (`Uid` comes
    before `VmRSS` in `/proc/PID/status`), so the 16, 38 and 20 GB of rounds 61 and 64 are void. Corrected sampler,
    A119b `all` again (`audit/g2/r65`, job 4460474, 851 s, 16/16 identical, 21/21 checks):
    - assemble stage, both `copy_assign` runs together: **25 GB** peak (was 2 x 13 GB before the parallel polish: no rise);
    - largest single `copy_assign` high-water mark: **50 GB** (assign; 40 GB in round 54, 34 GB before the goal);
    - all processes of the job summed (assign + the background flag job + their minimap2): **85 GB** peak.
    The fast defaults apply only on nodes with 128 GB or more and never under `--low-mem`, which also keeps flag in
    the foreground (driver line 1324), so nothing was changed. On a 128 GB node 85 GB leaves little room for other jobs.
79. **A table check that could not fail** (`non-assigned rows carry no copy`: its list was always empty, and the
    statement is false, tied rows name their best copy). Removed; 21 checks, all hold on `r63`, `kb_r63`, `r65`.
80. **Parallel polish on the opt-in paths** (`audit/g2/opt65`, job 4460475, gorilla KB3781, binaries before and after):
    `assemble --guide` gives the same `gtf`, `final.gtf`, `families.gtf`, `guided.gtf` and `bridges.tsv`;
    `RUSTLE_MONO_RULE=fallback,either RUSTLE_EXTEND_TRANSFRAGS=1` gives the same `gtf`, `final.gtf`, `families.gtf`,
    `bridges.tsv`.
- What rounds 62 and 63 were: checks of the written tables and a clean-environment run. The decision code of
  families (`mcl_families`) and assign was not read again line by line in this goal.

## Round 66: the open list applied, and the decision code of families and assign read

The list left "for the user" after round 65 was taken up: the deletion candidate, the rest of assemble's tail, and a
line-by-line read of the families and assign decision paths (two reviewers, read-only, each told what was already
checked). Every fix below had a failing test first.

81. **Deleted: `RUSTLE_WINDOW_SWEEP` / `RUSTLE_SWEEP_BATCH`** (r1507, 1,319 s against 480 s for assign): the store, its
    builder, `BamIndexCache::sweep_windows`, `sweep_supported`, `tests/window_sweep.rs` and its row in
    `PERFORMANCE.md`. About 300 lines. The prereg doc stays as the record.

**Families (`mcl_families`, the graph, SD damping, MCL, the copy table):**
82. **An SD file of another assembly weakened every edge and said nothing.** With no contig in common no pair is
    SD-linked, so `--sd-damp 0.25` multiplied every edge by 0.25 (exit 0, one log line "0 edge(s) SD-linked"). Now
    refused: "none of its contigs is a contig of the loci". Test `sd_pairs_report_whether_they_share_a_contig_with_the_loci`.
83. **`--min-exon-share` on a PAF without `cg:Z:` wrote zero families** (only with `--paf`; the pipeline's own
    all-vs-all always has CIGARs). Such records are now counted and the run refuses, as `--shared-exon-columns` did.
    Test `exon_share_counts_records_that_carry_no_cigar`.
84. **An exon ending before its start, or starting at 0, wrapped its length in release builds** and the locus lost its
    edges. The GTF parse now refuses it. Test `an_exon_with_impossible_coordinates_is_refused`. None in real input.
85. Dead branch removed (`Verdict::Below { contained }` could never be reached with `contained` true).
- Read and found clean: minus-strand PAF projection, 0/1-based conversions, edge symmetry, ordering (no hash order
  reaches an output), empty and single-locus input, NaN handling. Not read: annotation-mode units, the relations and
  hierarchy writers after the clusters, `absorb_loci`.
- Seen, a consequence of two registered rules, not changed: under `--sd-damp 0.25` at inflation 2.0 a lone pair with
  no SD link can never be a family (0.25 is below the pair boundary 0.2956).

**Assign (candidates, the tie gate, deciding sites, status, the tables):**
86. **The gate's log line was wrong under `--tied-prefilter`, which the pipeline always passes.** It counted only the
    tied molecules as "all" ("244802 of 244802 molecules ... 0 skipped"), so the warning for a BAM without `AS` tags
    could never fire. The light pass's counts are now used. Test `tied_prefilter_reports_the_same_gate_counts`.
87. **THIS ONE CHANGES O2 RESULTS. "A tied placement outside every unit" was tested by the placement's span, not its
    aligned blocks.** A tied placement whose intron merely crossed a copy counted as lying on it, so the molecule was
    not marked `tie_outside_catalog` and could end `assigned` although that placement was never scored. The code's
    own rule elsewhere is blocks (`has_block_in_any_copy_on`, `block_overlap`). Now `unit_under` (blocks) at the three
    places. Test `a_placement_lies_on_a_unit_only_where_its_blocks_do`. Effect, same input:
    - A119b: 6,773 of 200,817 rows change; 6,707 are ambiguous rows that only gain the flag; 53 rows go assigned ->
      tied, 8 go tied -> assigned. Assigned rows 4,860 -> 4,815, molecules 3,698 -> 3,662 (in more than one family
      946 -> 937); 103 quant rows move.
    - gorilla KB3781: 126 rows change, 1 assigned -> tied (722 -> 721), 3 quant rows move.
    To go back: restore the span test in `unit_under` (one line). Figures quoted from earlier runs use the old rule.
88. **`params.tsv` said `rna_editing_filter true` on the default path, where the filter is not applied** (pooled
    molecule observations pass no editing columns to the pair test). The record and the flag's help now say so
    (`false`). Applying the filter there would be a new rule: not done.
89. `--sibling-report` printed identity 1.0 and 0 columns for a read whose p-value was exactly 0 against every
    competitor (111 or more deciding sites): the nearest sibling is now named. Opt-in; no test (needs such a read).
90. `--read-star-junctions` decided nothing without `--read-star-unit` while `params.tsv` recorded it as on: refused.
    Test `read_star_junctions_is_refused_where_it_does_nothing`.
91. Dead work removed: a per-molecule CIGAR walk whose result was discarded (`genome_rate`).
    Also guarded: the per-record `fill_psv_obs` pass and the editing-column detection, which feed nothing on the
    pooled path with the read filter off (the binary's default). No new test: the end-to-end tables are the check.
- Read and found clean: the certificate arithmetic, the p-value tail, exact-tie LLR, determinism, subtraction
  guards, the skip list. Not read: the `--phase`, `--gtf`, `--posterior` and EM branches, the union certificate.

**Speed (products identical):**
92. **The layered assembly no longer holds the pipeline.** Families needs `PREFIX.families.gtf`; the layered run that
    writes `PREFIX.final.gtf` ends about 25 s later on A119b. Under `all` (not `--low-mem`) it finishes in the
    background, `all` waits for it before it ends and stops it when a stage fails. This replaces finding 58b: with
    the polish in parallel (round 61) the cost of two processes is those 25 s, not 70, and they are now off the path.
    Test `all_lets_the_layered_assembly_finish_under_families`.
93. **A large read-star alignment borrows idle thread slots.** One A119b family (2,025 molecules x 150 copies) took
    126 s on its slot's 16 threads. An alignment of 50,000 or more read-copy pairs now also takes the slots nobody
    holds (minimap2's output does not depend on `-t`). Test `a_large_alignment_borrows_the_free_slots_and_returns_them`.

## Round 67: closing runs of round 66

Single runs on one 72-core node each. Runs differ by 100 s and more between nodes: gorilla assemble took 206 s on
f0811 (`kb_r68`) and 381 s on f0575 (`kb_r69`) with binaries that differ only in assign code. Differences below
that are not findings.

| run | what | A119b `all` | assemble / families / assign / flag | gorilla KB3781 `all` |
|---|---|---|---|---|
| `r66`, `kb_r66` | fixes 81 to 90 and the dead work of 91 | 841 s | | |
| `r67`, `kb_r67` | + layered assembly in the background (92) | 845 s | 149 / 392 / 302 / 62 | 478 s |
| `r68`, `kb_r68` | + slot borrowing (93) | 834 s | 147 / 403 / 281 / 70 | 478 s |
| `r69`, `kb_r69` | + the guard of finding 91 (final binaries) | 819 s | 153 / 397 / 267 / 68 | 648 s (f0575), rerun `kb_r71` 500 s (f0332) |

- Products: every table and GTF of `r67`, `r68`, `r69` is identical to `r66`, and of `kb_r67` to `kb_r69` to `kb_r66`
  (16 of 16). Against the baseline before round 66 the four assign tables differ, by finding 87 and the `params.tsv`
  row of finding 88; the other 12 are identical.
- `tools/check_run_tables.py`: 21 of 21 checks hold on each run.
- Slot borrowing: the one large alignment went from 126.5 s to 120.1 and 120.6 s. Small; minimap2 does not scale there.
- Memory (A119b, sampled from /proc): the two assemble runs together 25 GB; one `copy_assign` up to 42 to 50 GB; all
  processes of the job together 85 to 88 GB. Gorilla 36 GB. Each fast setting has its own gate: assign's 24
  regions need 48 cores, 128 GB and the exclusion BED; assemble's region settings read the node's memory; the
  background layered assembly (92) has no memory check, only `--low-mem` turns it off. The peak while it overlaps
  families was not sampled on its own; it is inside the 85 to 88 GB job figure.
- Full suite on the final tree: 1,550 passed, 0 failed. Baseline before this goal: A119b 1,315 s, gorilla 594 s.

## Round 68: reading each assign window once (findings 65 and 66, measured, closed)

94. **One pass per window is not faster and doubles the memory.** The one-pass read already exists: assign without
    `--tied-prefilter` reads every window once with full records. Measured on one node, assign alone on the `r69`
    products, arms alternated (`audit/g2/r70`, job 4460601), six assign tables identical to `r69` in every arm:
    two passes 329 and 276 s, 40 and 42 GB; one pass 319 and 293 s, 85 and 89 GB. Only 244,802 of 6,557,717
    molecules in the windows are tied, so the second pass is small and holding every full record is what costs.
    What was measured is that existing one-pass path. The variant findings 65 and 66 named, one pass that keeps
    only tied molecules, was not built: it needs the placement counts in the side table (a format change), and my
    estimate of its gain, 30 to 40 s, is below the 53 s spread of the same arm on the same node. Not applied.

## Round 69: the opt-in branches read, and a family's windows read at once

A reviewer (read-only) read the branches round 66 had left: assign's `--phase`, `--posterior`, EM and part of the
union certificate; families' hierarchy and relations writers and `absorb_loci`. The tests below were written with
the fixes, not before them; each new function did not exist before. None of this touches a default product.

95. **`--emit-hierarchy` read every exon one base too far left.** 0-based half-open blocks were passed to code that
    expects GFF coordinates, so a duplicon or SD block ending where an exon starts counted as touching it. Now
    converted in `hierarchy_exons`. Test `hierarchy_exons_are_gff_coordinates`.
96. **`--posterior`'s zone columns were wrong:** the end was the largest copy start, the zone was empty for copies
    not named `DN_...`, and contigs could mix (read from the code, not run; the start now comes from the copy's
    span, not its name). Now the span of the zone's copies on the
    top copy's contig (`posterior_zone`). Test `posterior_zone_is_the_span_of_its_copies_on_one_contig`.
97. **`--absorb-loci` with `--emit-relations` could stop the run at the end** ("not a locus of the families GTF")
    after the clusters were written: absorbed genes were not kept resolvable as span-collided ones are. Now they
    are (`absorbed_under_keepers`). Not reproduced on data; test extends `dedupe` test in `family.rs`.
98. `--em-eps` and `--em-max-iter` are written to `params.tsv` when `--em` is on (two `--em` runs with different
    settings looked the same). Default `params.tsv` unchanged.
99. `--union-certificate`, `--phase`, `--posterior`, `--em` are refused with `--assemble-only`, where they did
    nothing and were recorded as on. Test `assignment_flags_are_refused_under_assemble_only`.
100. **Seen, not changed (a rule):** `absorb_loci` ignores `--absorb-ratio` when both loci have 0 reads (a GTF
    without read counts): any contained locus is absorbed. The code comment names the 0-read case, so it may be
    meant. `--absorb-loci` is opt-in and not adopted.
- Still not read: assign's `--gtf` branches and rescue block, the phase graph callees, the union helpers
  (`merge_same_locus`, `project_posterior`, `tie_row`, `recount_after_union`); `annotation_families`, the
  annotation-mode read-chain units, `write_locus_rep_copies`, the hierarchy external table.

101. **Tried and taken out: a family's copy windows read at once.** The slowest A119b family (147 windows) fetched
    in 86 to 123 s one window after another. On 8 plain threads it fetched in 14 to 16 s, but the stage did not get
    faster: assign alone, same node, old 332 and 261 s, new 275 and 273 s (`audit/g2/r73`, job 4460768), tables
    identical. On the rayon pool the stage got slower (459 and 482 s against 311 and 283 s, `audit/g2/r72`): a
    region waiting for its windows picked up other regions' work. Reverted; the code reads windows as before.
    **This closes findings 65 and 66 by measurement:** taking 70 to 100 s of fetch out of the slowest family did
    not move the stage, so removing only its second pass would not either (measured for families of 32 or more
    windows; the estimate in round 68 is replaced by this).
- Closing runs on the final binaries (`audit/g2/r74e`, `kb_r74`, jobs 4460933 and 4460934): A119b `all` 831 s
  (153 / 403 / 273 / 56), gorilla 479 s; 16 of 16 products identical to `r69` and `kb_r66`; 21 of 21 table checks on
  each; job peak 87 GB. Full suite on the final tree: 1,553 passed, 0 failed.

## Round 70: the rest of the opt-in code read

A second read-only pass over what round 69 listed as not read. Every finding is from reading the code; none was
reproduced on data, and none touches a default product.

Applied (no test: the flag has no fixture anywhere in the suite, which is itself a gap):
102. **`--no-cross-family-exon-overlap` removed units trimmed to nothing and so shifted every later unit's index**,
    while `merged_into` and the read-through evidence still pointed at the old positions: a panic past the end, or
    a unit written as merged into the wrong one. The emptied unit now keeps its place and is skipped when written.
103. **Same flag: `units.fa` kept the sequence of the exons before the trim.** Rebuilt after the trim (`unit_seq`).

Applied after the first report of this round, each with a test (108 and 109 failed first; 104's test came with it):
104. `--rescue-singletons`: every unspliced rescued read shared the empty intron chain, so reads of copies far apart
    became one transcript from the smallest start to the largest end. Unspliced reads now group by the copy they
    were assigned to (`rescue_key`). Test `unspliced_rescues_group_by_copy_and_spliced_ones_by_chain`.
108. `sd_blocks`: hulls are 1-based closed, SD pairs 0-based half-open, and the overlap test mixed them: a pair
    starting on the base after a hull linked it. Test `sd_blocks_read_hulls_and_pairs_in_their_own_coordinates`.
109. `--union-certificate` left `uniq` / `uniq_agree` in `families.tsv` unchanged for rows it tied, so `uniq` could
    exceed the assigned count. `recount_after_union` now gives the contribution back (the tie-outside gate did this
    on its own; it now goes through the same function). Test extended in
    `tie_outside_gate_after_all_regions_demotes_only_assigned_outside_rows`.

Found, not applied, with the reason:
105. `--phase`: the reference base of a column is not complemented for families whose frame is the minus strand,
    so the reference walk in `.phase.gfa` and the divergence counts are wrong there. The fix needs to know which
    copy supplied each column, which the family record does not carry: a data-structure change, not a small fix.
106. `--rescue-singletons` with `--assembly-polish mono|full`: rescued lines carry `support`, the polish reads
    `reads` and takes 0. Two fixes are possible and they give different outputs (keep rescues out of the polish, or
    let them count their support): a choice about the experiment.
107. `--core-majority-inclusive` on a family of two: the threshold `n/2 - 1` is 0, so gaps between SD hits count
    as core. That is the registered rule's own arithmetic at n = 2; changing it changes the rule.
- Read and clean: the union helpers, `copy_status`, the GFA writers, `write_locus_rep_copies`, `relations_tsv`,
  the hierarchy external table, the `--gtf-copy-set` grouping.
- Still not read: `annotation_families` graph, MCL and cluster folding (`family.rs` 61449-61790, 61950-62586,
  63115-63690), `family_container` readers, the polish callees (`tss_decide`, `tes_*`, `junction_snap`,
  `mono_models`, regroup, sub-chain), `copy_span_lifts`, `tie_set_of`, `kept_segment`.
- Closing runs on the final tree (`audit/g2/r76e`, `kb_r76`, jobs 4461254 and 4461253): A119b `all` 847 s, gorilla
  477 s; 16 of 16 products identical to `r69` and `kb_r66`; 21 of 21 table checks on each. Full suite: 1,555
  passed, 0 failed. (Before findings 104, 108, 109: `r75e` 813 s, `kb_r75` 482 s, same products.)

## For the user (not applied, or not the audit's to decide)

- **Finding 87 changes O2 results** (assigned rows on A119b 4,860 -> 4,815). It is a bug fix under the code's own
  block rule, with a failing test first, but it sits in the read-assignment layer: keep it or restore the span test
  (one line in `unit_under`). Figures quoted from earlier runs use the old rule.
- **Families is the longest stage (about 400 s of 820 on A119b):** 98% is the minimap2 all-vs-all, already one batch
  (`-K 4G`, round 22). Round 66 fixed three input faults there and found no further lever that keeps the products;
  the remaining lever is exon-only alignment, a rule question.
- **One process writing both GTFs (finding 58b)** is replaced by finding 92: the second run is off the path.
- **Registered questions untouched:** `RUSTLE_MONO_RULE` (r1626), `RUSTLE_EXTEND_TRANSFRAGS`, the tie-group locus rule,
  exon-only alignment for families, retiring `merged`, which GTF feeds guided mode.
- **Observed, already reported by the run itself:** on A119b 645 copy pairs of different families share genomic
  sequence, and those pairs double-claim 423 molecules (a molecule counted once per pair, `a1.assign.log`). Over all
  families, with finding 87: 4,815 assigned rows are 3,662 molecules, 937 of them assigned in more than one family.
- Found in opt-in code and not applied: findings 105, 106, 107 (round 70), each needing your choice. Not read: the list at the end of round 70.
- `absorb_loci` with 0-read loci (finding 100): a rule, yours to decide.
- Nothing is committed.

## Decisions by the user (2026-11-11)

- Finding 87: **keep** the block rule. Assign figures from runs before round 66 use the old span rule.
- Findings 105, 106, 107: leave recorded, not fixed.
- Finding 100 (`absorb_loci` at 0 reads): leave as is.
- The unread opt-in code listed in round 70: stop here, no further pass.
- Added after the decisions: test `unit_seq_joins_the_exons_and_follows_the_strand` pins the routine finding 103
  uses (test-only change; binaries unchanged). Finding 102 (unit positions kept after the trim) still has no test:
  it needs an annotation-mode fixture with two families sharing an exon, which the suite does not have.

## Audit of 2026-11-19 (on PSU Roar, one ponytail-audit pass over the whole repo)

20 findings were reported. The user asked to apply every must-fix and should-fix finding (1-15) except finding 2. On
finding 2 the user said: the PREFIX lock was meant for a home computer (i5, 32 GB) and is no longer needed on the HPC
clusters, so it is not extended and no new lock is added (finding 3 uses a unique temporary name instead). Findings 16-20
(nice to have) were not applied.

| # | what | done | check |
|---|---|---|---|
| 1 | `all` ran its child `assemble` on the primary pool (the assemble-alone default, register 1812) while `PREFIX.assemble.env` recorded good; `rustle_samples.sh` had the same gap | `stage_all` passes `--seed-pool "$SEED_POOL"` to its children; `rustle_samples.sh` passes `--seed-pool ${RUSTLE_SEED_POOL:-good}` unless a pool flag is given | test `all_passes_its_pool_to_the_child_assemble`; SIRV `all --to families` logs "seeding pool: good" |
| 3 | two `rustle_samples.sh` runs on one prefix wrote the same `POOLED.tmp` | `POOLED.tmp.$$` | `bash -n` |
| 4 | `PREFIX.final.pid` held a bare pid: a reused pid was waited on or killed; a killed layered run looked finished | the pid file holds pid, host and start time (`final_alive`); the background run writes `PREFIX.final.rc`; no rc file = the run stopped, `all` drops `assemble.env` and fails | dispatch test `all_lets_the_layered_assembly_finish_under_families` (stub updated); SIRV run wrote `s.final.rc` = 0 |
| 5 | `copy_assign` found a contig missing from the FASTA only when its region ran (hours in); a BAM of another build ran on the wrong bases; after a failed region the others ran on | startup check against `FASTA.fai` (refused without one; contig names, region ends, and under `--genome-wide` exact lengths); the parallel sweep stops taking regions after a failure | full `cargo test`; SIRV md5 unchanged |
| 6 | `RUSTLE_SEED_SPLICED=support` with `--read-junction-snap`, `--extend-transfrags` or the direct modes checked chains the support rule never saw | refused in that combination | build |
| 7 | the representatives cache keyed a file-valued setting (`RUSTLE_FOOTPRINT_WINDOWS=w.bed`) on its path alone | `env_fingerprint` adds the file's size and mtime when a value names a file; other keys unchanged | test `the_cache_key_changes_when_a_file_named_by_a_setting_changes` |
| 8 | a set but unparsable setting (`0,98`) silently became the default | `rustle::env_parse` / `env_opt` stop the run naming the variable, at the sites the audit named; `RUSTLE_SEED_SPLICED` / `RUSTLE_SEED_UNSPLICED` refuse unknown words | test `unset_or_empty_is_the_default_and_a_set_value_must_parse` |
| 9 | the resume record held the speed knobs the driver derives from the allocation, so a resume in a differently sized job redid every stage | `merged_env` leaves out `RUSTLE_{BAM,FAMILY,STAR,REGION}_THREADS` and `RUSTLE_REGION_READ_BUDGET` | test `the_resume_record_ignores_the_derived_speed_knobs`; the existing key test unchanged |
| 10 | `merged` ignored the genome's `.exclude.bed` sidecar; a product rewritten after the record was skipped as the merged run's | `merged` is refused when the sidecar exists and `RUSTLE_EXCLUDE_BED` is unset (as for the other skip lists); a phase is skipped only when its product is not newer than `PREFIX.merged.env` | driver tests |
| 11 | `eq_classes.py` ordered loci with equal start by set order (run-to-run differences) | the locus name breaks ties | pytest |
| 12 | `eq_classes.py`, `psi_bunches.py` took the read's alignment strand, not the transcript's (minimap2 `ts:A:-`) | `tx_strand()` flips on `ts:A:-`. The bench scripts behind rows 1826 and 1828 are NOT changed: their reads are oriented FLNC, where the two agree, and their registered numbers stay reproducible | pytest |
| 13 | no test ran the tools' selftests; `gtf_exclude_bed --selftest` exited 1 on success; nothing checked the copy table famreps writes | `tools/test_selftests.py` runs 9 selftests; `pixi run test-py` runs them with the bench tests (155 pass); new `tools/copy_table_check.py` (copy_idx 0..n-1 per family, unique, FASTA records = rows) runs inside `stage_famreps` before the table is replaced (passes on the OR6737 and KB3781 tables) | pytest |
| 14 | a failed index build left its ~13 GB temporary file beside the genome | removed on failure and by an EXIT trap added to the snapshot's trap, not replacing it | trap checked in isolation |
| 15 | unquoted `$RUSTLE_PYTHON` | **withdrawn**: quoting broke 14 dispatch tests, because `RUSTLE_PYTHON` is meant to be a command with arguments (`pixi run python`, the tests' `bash py.sh`); kept unquoted with a `shortcut:` comment (a python path with a space is not supported) | dispatch tests |

Default products: SIRV `final.gtf` md5 `5cb0e7e4e7eb4cccc8f9def75cce2f65` (assemble alone and `all`) unchanged. Full
`cargo test` passes (37 test binaries), the 90 driver tests pass, `pixi run test-py` 155 pass. Behaviour changes a user
can meet: `copy_assign` refuses a FASTA without `.fai` or of another build; a malformed `RUSTLE_*` number stops the run;
`all` on a pool-less command line really seeds with the good pool (before this fix, `all` runs seeded with primaries
whatever `assemble.env` said, see finding 1).
