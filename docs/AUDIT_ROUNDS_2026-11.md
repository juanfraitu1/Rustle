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
