# O1 layers, phases 1-2: sedef-rs as the SD and duplicon layers — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Vendor sedef-rs into Rustle as `utilities sd-scan` and a cached `sd` driver stage, let `mcl_families --emit-hierarchy --sd-cache` read its SD pairs and duplicons, and run the pre-registered gates V1/V1b (sedef-rs vs SEDEF on CHM13 v2.0 and gorilla) and V2 (sedef-rs duplicons vs DupMasker on the human F_O universe).

**Architecture:**
- **Vendoring.** The sandbox's seven modules are copied verbatim, with only their crate paths rewritten.
- **Three recorded changes:**
  - stages 1-2 become whole-genome-ready: chain merging stays within one chromosome pair (the sandbox fused blocks of different chromosome pairs), window buffers are released and pair counting is sharded, so a whole genome fits this machine's 25 GB;
  - duplication families link each segment's two copies (the sandbox's documented intent, which its code missed) and keep their copies;
  - a new `scan` module runs the four stages through files, on a local thread pool, and writes a duplicon table in the `--hierarchy-duplicons` format.
- **Driver:** a new `sd` stage, keyed by the FASTA's md5 and sd-scan's version key.
- **Hierarchy flag:** `--sd-cache PREFIX`, and `parse_sd` learns the sd-scan pairs format.
- **Scoring:** a new Python script scores V1. V2 reuses register row 1263's `score_units.py` and plan 1's `hier_match.py`, both unchanged.

**Tech Stack:** Rust 2021 (clap 4, rayon, flate2, rustc-hash 2), bash, Python 3 (pytest, scipy).

**Spec:** `docs/superpowers/specs/2026-10-13-o1-layered-hierarchy-design.md` (45098918). This plan covers:
- section 1: the `sd` stage, and the hierarchy's SD and duplicon layers;
- section 2: V1, V1b and V2;
- section 4: steps 1-2.

## Facts measured while writing this plan (2026-10-13)

1. **The sandbox reproduces its own reference.** Built from its sources (offline, from cached crates), it reproduces its chimp chr1 reference outputs (`/mnt/c/Users/jfris/Desktop/sd/sedef-rs/data/chr1.*`) byte for byte. Its `Cargo.lock` pins `either 1.19.0`, which is not cached here.

   | Output | Rows |
   |---|---|
   | candidates | 32,925 |
   | blocks | 10,663 |
   | refined | 10,663 (ok 10,472, self 39, weak 152) |
   | segments | 8,713 |
   | families | 3,316 |

   | Stage | Peak RSS (KB) | Time |
   |---|---|---|
   | seed | 3,145,772 | 10.5 s |
   | chain | — | 0.04 s |
   | refine | 954,900 | 15.0 s |
   | decompose | 560,896 | 30.7 s |

2. **Two sandbox properties block whole-genome runs; Task 2 fixes both.**
   - **Memory.** On chimp chr1-3 (`/mnt/linuxdisk/tmp/sedef_probe/chimp3.fa`, 620 Mbp; md5 7a1f8b2fb4ae4faf71b6aca2d836af10), the sandbox's seed peaked at 9,321,528 KB in 53.4 s (42.9 s of it counting 81,017 candidates). Its candidates are kept at `/mnt/linuxdisk/tmp/sedef_probe/cand.c3.tsv` (81,018 lines; md5 2ba16689fdf37b15e41656736917ba26).
     - Every window keeps its 2,000-hash k-mer buffer (16 KB) behind a 256-hash sketch: about 25 GB on a human genome.
     - The pair-count map holds every window pair sharing at least one sketch value, while only pairs sharing at least 40 survive.
     - Projected whole-genome peak: about 50 GB, against 25 GB of RAM.
   - **Chain merging ignores chromosomes.** `chain_and_merge` (pass 2) and `merge_blocks` union overlapping chains and blocks by window coordinates alone.
     - On a multi-chromosome genome, blocks of different chromosome pairs fuse, and the fused block takes the names of whichever entry a HashMap met first: wrong, and different from run to run.
     - On chimp chr1-3, two runs gave different segment counts (12,539 and 12,537) from identical candidates. A one-chromosome input cannot show it.
3. **The sandbox's "families" are not duplication families.**
   - `decompose::families` unions overlapping segment endpoints, but never a segment's own two endpoints.
   - Its documentation says otherwise: "Connected components are duplication families", "inverted copies included".
   - On a toy genome with two planted duplications, every one of its 10 families has a single copy.
   - Linking the two endpoints (Task 3) gives:
     - toy genome: 5 families of 2 copies;
     - chimp chr1: 1,033 families (766 with 2 copies, none with 1). The largest five hold 3,750 (21.5% of 17,426 copies), 3,626, 2,734, 1,690 and 1,178 copies. The span sum is 9,210,601 bp.
   - Segments are unchanged by the link.
4. **The toy genome (Task 4's fixture) through the sandbox:**
   - An 8 kb inverted copy is called as four 2 kb pairs: `chain` links windows only along forward diagonals.
   - A 5 kb direct copy is called as one pair.
   - Refined ends fall one base short: `14999` for a copy ending at 15,000; an inverted partner starts at `124001` for a window starting at 124,000.
5. **Masked bases.** Lowercase (soft-masked) bases count as k-mer bases (`kmer::encode_base`); `N` breaks k-mers.
6. **Genomes and SEDEF references.**

   | Input | Path | Notes |
   |---|---|---|
   | CHM13 v2.0 | `/mnt/linuxdisk/home/juanfraitu/winloci_data/chm13v2.0.fa` | 3,117,292,070 bp; 25 records (chr1-22, X, Y, M) |
   | Gorilla | `/mnt/linuxdisk/home/juanfraitu/_from_wsl/winloci_scratch/GGO.fasta` | 3,545,850,636 bp; 26 records; NC_ accession as first word, then a description |
   | CHM13 v1.0 | `/mnt/linuxdisk/home/juanfraitu/winloci_data/soto_replication/t2t-chm13-v1.0.fa.gz` | bgzip; 3,056,916,522 bp; chr1-22 and X; lengths equal the DupMasker/F_O `chm13.draft_v1.0` |
   | SEDEF, CHM13 | `/mnt/linuxdisk/tmp/o1_close/lift/sedef_pairs.v2.bed` | 150,833 rows, 34 columns, no header; md5 7bd521edff9bd1e27b8141027da37a01 |
   | SEDEF, gorilla | `/mnt/linuxdisk/home/juanfraitu/winloci_data/GGO_sedef_final.bed` | 253,029 rows plus one `#` line, 34 columns; md5 991ed34384832354c02df227e27f65b6 |
   | v1.0 → v1.1 chain | `/mnt/linuxdisk/home/juanfraitu/winloci_data/soto_replication/v1.0_to_v1.1_rdna_merged.chain` | |

   - The CHM13 SEDEF file is v1.0 SEDEF lifted to v2.0. Every sequence name it uses is in chm13v2.0, and no end lies beyond its chromosome.
   - The gorilla SEDEF file names all 26 GGO.fasta sequences, including the mitochondrion NC_011120.1, and no end lies beyond its sequence's length.
   - Identity in both is column 21 (fracMatch, a fraction).
7. **V2 inputs (the r1263 protocol):**
   - F_O: `/mnt/linuxdisk/tmp/ourhom/loci.paf` + `loci.gff3`;
   - DupMasker: `/mnt/linuxdisk/home/juanfraitu/winloci_data/duplicons/chm13.draft_v1.0_plus38Y_dupmasker_colors.bed`;
   - SEDEF v1.0: `/mnt/linuxdisk/home/juanfraitu/winloci_data/soto_replication/final_v1_clean.bed`;
   - Soto external families: `/mnt/linuxdisk/tmp/hier_flag/soto_external.tsv`;
   - registered comparator: `/mnt/linuxdisk/tmp/hier_flag/final/`. On it, `hier_match.py --map key` gives `exact_any` = any 245, expansion 215, **expansion_not_family 10**, family 231 (`/mnt/linuxdisk/tmp/o1_close/fo_check/key.json`).
8. **Figure code hashes today.** figures/samples.py `driver_stage_code` gives:

   | Stage | Hash |
   |---|---|
   | assemble | cf14d218d1cbe774 |
   | families | 58c0422b24a99e8d |
   | candidates | 635739635b54c07a |
   | catalog | 4c987eb6e5839a46 |
   | assign | 6016f8834bed63a8 |
   | assign_fam | e10889fd43f16d99 |
   | merged | f340f2721ccc3065 |
   | flag | 9cff53027b6ef944 |
   | index | 21eb18b2cc02dabe |

   The hash leaves out other stages' `stage_*()` bodies and every column-0 `case "$STAGE" in` … `esac` block, but keeps shared top-level lines. So the new stage must be a function plus a `case` block, nothing else.
9. **The module registry's parser.** In `docs/MODULE_STATUS.md`, the registry test reads the first backticked token of every line. A token ending in `.rs` counts as a listed `src/family.rs` module, so a note about `src/sedef/` must not open with such a token.
10. **The plan's Rust code was compiled and run before the plan was handed over.** Scratch crate: `/mnt/linuxdisk/tmp/sedef_probe/plancheck` (the sandbox modules, with Tasks 2-4 applied as written here).
    - All 25 tests pass, including the toy genome's exact pairs, families and duplicon tables.
    - Chimp chr1: candidates, blocks, pairs and segments identical to the sandbox's; families `1033 / 0 / 17426 / 3750 3626 2734 1690 1178 / 9210601`; peak 758,880 KB, 68 s.
    - Chimp chr1-3: candidates identical to the sandbox's (set and order); blocks 14,558 (14,544 with the merge bug), pairs 14,558 (ok 14,262), segments 12,554, families 1,568, duplicon rows 24,466; two runs byte-identical on every table; peak 1,804,760 and 1,900,468 KB in two builds, 212 s.
    - Power-law projection of peak memory (from chr1 and chr1-3): CHM13 v2.0 about 7.3-8.0 GB, gorilla about 8.2-9.0 GB.

## Rulings (decisions this plan takes; the spec does not settle them)

- **R1.** The sandbox's `main.rs` (its four-stage CLI) is not vendored; `scan` replaces it. Parameter sweeps stay in the sandbox.
- **R2.** Sharded counting (Task 2), and window b's sequence as the candidate order's last key, change no output on a one-sequence genome. On many sequences they fix the order of ties that the sandbox left to its hash map.
- **R3.** Families link the two endpoints of each segment (Task 3).
  - This is the sandbox's documented intent, and the only reading under which a family can serve as a duplicon.
  - `span_bp` is summed per chromosome.
  - `families.tsv` therefore differs from the sandbox's, while candidates, blocks, pairs and segments must not.
- **R4.** A sequence's name is the first word of its FASTA header, which is how the GFFs and SEDEF files name sequences.
- **R5. V1's universe.**
  - CHM13: the v2.0 bases that have a v1.0 ancestor (the chain's aligned blocks: chr1-22 and X; no Y, no M, no rDNA arrays). The reference is v1.0 SEDEF, lifted.
  - Gorilla: the whole FASTA.
  - The lift dropped 925 SEDEF rows (151,758 → 150,833). Where their regions fall inside the universe, they count against sedef-rs's precision. This is recorded, not corrected.
- **R6. V2 runs on the F_O universe (CHM13 v1.0).**
  - That is where r1263's protocol, comparator and baseline exist. It costs a third whole-genome sd-scan (v1.0).
  - Only the duplicon track is swapped; the SD layer stays SEDEF v1.0, so the comparison isolates the duplicons.
  - `--sd-cache` itself is exercised on A119b (Task 9), descriptively.
- **R7. V2's second bar.**
  - It uses `exact_any.expansion_not_family` from `hier_match.py --map key` (close-O1 amendment 7), recomputed on the re-run DupMasker arm.
  - Bar = ⌈0.9 × that count⌉: 9 if the count is 10 again.
  - It is computed once, at the held-out step.
- **R8. Whole-genome runs.**
  - They go through the driver stage into a persistent cache, `/mnt/linuxdisk/home/juanfraitu/sd_cache`, under a 23 GiB address-space cap.
  - Task 4 gates them on the memory projection.
- **R9.** `count_shards` defaults to 256, and window sketches and genome sequences are shrunk to their length (Tasks 2 and 4). Measured (Fact 10): chimp chr1 0.76 GB, chr1-3 1.8-1.9 GB; projected about 7-9 GB for a whole genome. Each extra shard costs one pass over the sorted entries, minutes in total.
- **R10.** Two more fixes against the sandbox, each output-neutral on chimp chr1 (every table byte-identical):
  - chain merging is restricted to one chromosome pair (Task 2): a bug fix;
  - `dedup_segments` orders tied representatives by (a1, b1, strand) (Task 3): its order no longer follows HashMap order.

## Global Constraints

- **Branch and repo:** work in `/mnt/linuxdisk/home/juanfraitu/rustle_m2_soto` on branch `o1-layers`, created from `main` in Task 1. Never commit to `main`.
- **Builds:**
  - only as `CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target TMPDIR=/mnt/linuxdisk/tmp cargo build --release ...`, with output redirected to a file, under `RLOCK_TIMEOUT=900 bash tools/rlock.sh heavy`;
  - tests run the same way (`cargo test --release ...`);
  - binaries land in `/mnt/linuxdisk/home/juanfraitu/rustle_target/release/`.
- **Jobs:**
  - one heavy job at a time (`tools/rlock.sh heavy`); light jobs (< 2 GB, < 3 min) via `tools/rlock.sh light`;
  - never `pkill -f`;
  - jobs longer than 10 minutes: Bash `run_in_background: true`, then wait for the completion notification (no sleep loops).
- **Commits:**
  - path-limited (`git add <paths>`), never `git add -A` or `git add .`;
  - never add `bench/hierarchy/*` or `bench/ideal_expression/*` (another session's untracked files);
  - every commit message ends with:
    ```
    Co-Authored-By: <your model name> <noreply@anthropic.com>
    Claude-Session: https://claude.ai/code/session_015Xtq7GodhTGPYBzHDp3o4z
    ```
- **The sandbox** `/mnt/c/Users/jfris/Desktop/sd/sedef-rs` is the user's: read it, never write it.
- **Files:** delete or overwrite nothing outside this plan's outputs. Outputs go under `/mnt/linuxdisk/tmp/o1_layers/` and `/mnt/linuxdisk/home/juanfraitu/sd_cache/`.
- **No existing product may change:**
  - `mcl_families` without `--sd-cache`, and `parse_sd` on SEDEF rows, behave as before;
  - every existing driver stage keeps its `driver_stage_code` (Fact 8).
- **Pre-registration:** the prereg (Task 7) is committed before any whole-genome sd-scan output is compared with SEDEF or DupMasker.
- **Python tests:** `python3 -m pytest -q <dir>`.

## Review Focus

1. **A gzipped or bgzipped genome.** CHM13 v1.0 is `.fa.gz`; it must give the same tables as the plain FASTA. Test: Task 4, `a_gzipped_genome_gives_the_same_tables`.
2. **FASTA headers with descriptions, and two records with the same first word.** Names are the first word; a repeated name is refused, never silently merged. Test: Task 4, `records_are_named_by_the_first_word_and_a_repeated_name_is_refused`.
3. **A genome without duplications.** Every table is written with its header only, and chain, refine and decompose do not panic on empty input. Test: Task 4, `a_genome_without_duplications_writes_header_only_tables`.
4. **A stale cache.** The same FASTA path with new content, a new sd-scan version, or a scan that died half-way must each lead to a new scan. Tests: Task 5, `a_changed_genome_or_version_key_is_scanned_again` and `a_failed_scan_is_not_cached`.
5. **Mixed hierarchy inputs.** `--sd-cache` together with `--hierarchy-duplicons` or `--hierarchy-sd` is refused, and sd-scan's `self`/`weak` rows never enter the SD layer. Tests: Task 9, `sd_cache_names_both_hierarchy_inputs_and_excludes_the_two_flags` and `parse_sd_reads_sd_scan_pairs_ok_rows_with_percent_identity`.

---

### Task 1: Vendor sedef-rs verbatim (`src/sedef/`, its 12 tests)

**Files:**
- Create: `src/sedef/mod.rs`; `src/sedef/{anchors,chain,decompose,fasta,kmer,refine,seed}.rs` (copies)
- Modify: `src/lib.rs` (after `pub mod family;`), `Cargo.toml` (`[dependencies]`), `docs/MODULE_STATUS.md` (a note at the end)

**Interfaces:**
- Consumes: nothing.
- Produces: `rustle::sedef::{anchors, chain, decompose, fasta, kmer, refine, seed}` with the sandbox's public API unchanged. Signatures later tasks use:
  - `fasta::read_fasta(path: &Path, chr_filter: Option<&str>) -> Result<Vec<(String, Vec<u8>)>, Box<dyn Error>>`
  - `seed::{SeedParams, WindowSketch {seq: u32, start: u64, end: u64, sketch: Vec<u64>}, SeqInfo {name, len}, Candidate {a: usize, b: usize, shared: u32}, build_windows(&[(String, Vec<u8>)], &SeedParams) -> (Vec<WindowSketch>, Vec<SeqInfo>), find_candidates, est_jaccard(u32, usize) -> f64}`
  - `chain::{shared_for_identity(f64, usize, usize) -> u32, load_pairs(&Path, u32), load_blocks(&Path), chain_and_merge(&[PairRec], &ChainParams) -> (Vec<Block>, ChainStats), ChainParams {window: u64, gap: u32, min_score: u32, merge_frac: f64, sketch: usize, k: usize, min_cell_id: f64, min_block_id: f64, min_shared_cell: u32}, Block {chr_a, a0, a1, chr_b, b0, b1, cells: u32, score: u32, est_id: f64}}`
  - `refine::{RefineParams {k, max_gap: u64, diag_bin: i64, min_span: u64, max_pairs: u64}, Refined {chr_a, ra0, ra1, chr_b, rb0, rb1, strand: char, anchors: usize, est_id: f64, ok: bool, self_pair: bool, stage2_id: f64}, refine_blocks(&[Block], &HashMap<String, Vec<u8>>, &RefineParams) -> Vec<(Refined, BlockAnchors)>, covered_len(&[(u64, u64)]) -> u64}`
  - `anchors::{BlockAnchors, write(&Path, usize, &[BlockAnchors]), read(&Path) -> Result<(usize, Vec<BlockAnchors>), _>}`
  - `decompose::{DecomposeParams {k, max_gap: u64, diag_bin: i64, min_seg: u64, max_pairs: u64}, Segment, Family, decompose(Vec<BlockAnchors>, &DecomposeParams) -> (Vec<Segment>, Vec<Family>, DecomposeStats)}`

This task ports tested code with its own tests, so there is no red step. The tests come with the files.

- [ ] **Step 1: Branch and output directory**

```bash
cd /mnt/linuxdisk/home/juanfraitu/rustle_m2_soto
git switch -c o1-layers main
mkdir -p /mnt/linuxdisk/tmp/o1_layers
```

- [ ] **Step 2: Copy the seven modules and rewrite their crate paths**

```bash
cd /mnt/linuxdisk/home/juanfraitu/rustle_m2_soto
SB=/mnt/c/Users/jfris/Desktop/sd/sedef-rs/src
md5sum $SB/{anchors,chain,decompose,fasta,kmer,refine,seed}.rs
mkdir -p src/sedef
for m in anchors chain decompose fasta kmer refine seed; do cp "$SB/$m.rs" src/sedef/; done
sed -i -E 's/crate::(anchors|chain|decompose|fasta|kmer|refine|seed)\b/crate::sedef::\1/g' src/sedef/*.rs
for m in anchors chain decompose fasta kmer refine seed; do
  diff <(sed -E 's/crate::sedef::/crate::/g' src/sedef/$m.rs) "$SB/$m.rs" > /dev/null && echo "$m ok" || echo "$m DIFFERS"
done
```

Expected md5s:

| File | md5 |
|---|---|
| anchors.rs | 4b31457182ff3b13547434ac20c2384d |
| chain.rs | eadd615dacc6da648ffc480114b8fc03 |
| decompose.rs | dcbf640590c583ae1d4ff6545230e0dd |
| fasta.rs | 7a364ba96b2c631ea224bd7a3d1e2a01 |
| kmer.rs | f56c0eedb15d36fe5455d4f19befb896 |
| refine.rs | 67e9c83a01f7fc4763b204221083ff90 |
| seed.rs | 9d79d48c254691514d0672dd35711a19 |

The loop prints seven `ok` lines. If an md5 differs, the sandbox changed since this plan: STOP and report.

- [ ] **Step 3: Write `src/sedef/mod.rs`**

```rust
//! **STATUS:** OTHER-BINARY — `utilities sd-scan` and the driver's `sd` stage; opt-in input of
//! `mcl_families --emit-hierarchy --sd-cache`. Adoption as the SD / duplicon source is decided by
//! `docs/PREREG_sedef_rs_layers_2026-10-13.md`.
//!
//! Segmental duplications (SDs) of a genome with sedef-rs, a pure-Rust re-implementation of the SEDEF/BISER cascade
//! (seed -> chain -> refine -> decompose), vendored on 2026-10-13 from the user's sandbox
//! `/mnt/c/Users/jfris/Desktop/sd/sedef-rs` (docs/superpowers/specs/2026-10-13-o1-layered-hierarchy-design.md). The
//! sandbox stays the place to experiment; changes are ported here deliberately and listed below.
//!
//! Vendored files (md5 of the sandbox copies): anchors.rs 4b31457182ff3b13547434ac20c2384d, chain.rs
//! eadd615dacc6da648ffc480114b8fc03, decompose.rs dcbf640590c583ae1d4ff6545230e0dd, fasta.rs
//! 7a364ba96b2c631ea224bd7a3d1e2a01, kmer.rs f56c0eedb15d36fe5455d4f19befb896, refine.rs
//! 67e9c83a01f7fc4763b204221083ff90, seed.rs 9d79d48c254691514d0672dd35711a19. The sandbox's main.rs (its four-stage
//! CLI) is not vendored.
//!
//! Changes against the sandbox, each with its test:
//! - every `crate::<module>` path reads `crate::sedef::<module>`.

pub mod anchors;
pub mod chain;
pub mod decompose;
pub mod fasta;
pub mod kmer;
pub mod refine;
pub mod seed;
```

- [ ] **Step 4: Register the module and its dependency; add the registry note**

In `src/lib.rs`, after the line `pub mod family;`, add:

```rust

/// Vendored sedef-rs (segmental duplications): `utilities sd-scan`, the driver's `sd` stage.
pub mod sedef;
```

In `Cargo.toml`, after the line `fxhash = "0.2"`, add:

```toml
rustc-hash = "2" # src/sedef (vendored sedef-rs): FxHashMap / FxHashSet, as in the sandbox
```

At the end of `docs/MODULE_STATUS.md`, append the following lines. No line may start its first backticked token with a name ending in `.rs` (Fact 9).

```markdown

**Vendored crate, not a family module:** `src/sedef/` (sedef-rs, segmental duplications; vendored 2026-10-13 from
the user's sandbox) — OTHER-BINARY: `utilities sd-scan` and the driver's `sd` stage; an opt-in input of
`mcl_families --emit-hierarchy --sd-cache`. Its adoption as the SD and duplicon source is decided by
`docs/PREREG_sedef_rs_layers_2026-10-13.md`.
```

- [ ] **Step 5: Build and run the ported tests, plus the module-status tests**

```bash
cd /mnt/linuxdisk/home/juanfraitu/rustle_m2_soto
RLOCK_TIMEOUT=900 bash tools/rlock.sh heavy env CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target TMPDIR=/mnt/linuxdisk/tmp \
  cargo test --release --lib sedef:: > /mnt/linuxdisk/tmp/o1_layers/t1_test.log 2>&1; echo "exit $?"
grep -E "^test |test result" /mnt/linuxdisk/tmp/o1_layers/t1_test.log
RLOCK_TIMEOUT=900 bash tools/rlock.sh heavy env CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target TMPDIR=/mnt/linuxdisk/tmp \
  cargo test --release --lib module_status > /mnt/linuxdisk/tmp/o1_layers/t1_status.log 2>&1; echo "exit $?"
grep "test result" /mnt/linuxdisk/tmp/o1_layers/t1_status.log
```

Expected:
- `exit 0` twice;
- 12 ported tests `ok`:
  - `sedef::anchors::tests::roundtrip`
  - `sedef::chain::tests::{diagonal_run_forms_one_chain, phase_drift_stays_in_one_chain, distinct_nested_duplications_do_not_merge, phase_shifted_fragments_merge_into_one_block}`
  - `sedef::kmer::tests::{canonical_is_strand_symmetric, skips_ambiguous_bases}`
  - `sedef::refine::tests::{refines_direct_duplication, detects_inverted_duplication, mosaic_hub_deflates_identity, stored_anchors_match_fresh_extraction_on_refined_coords, combinatorial_junk_doesnt_form_long_runs}`
- `test result: ok. 12 passed`;
- the module-status tests pass.

Dead-code warnings from `src/sedef` are expected (the sandbox had 7) and are left alone.

- [ ] **Step 6: Commit**

```bash
git add src/sedef/ src/lib.rs Cargo.toml Cargo.lock docs/MODULE_STATUS.md
git commit -m "sedef: vendor sedef-rs (seed, chain, refine, decompose) as src/sedef with its 12 tests

Verbatim copies of the user's sandbox modules (md5s in src/sedef/mod.rs), crate paths rewritten.
Plan docs/superpowers/plans/2026-10-13-o1-layers-phase-1-2.md, Task 1.

Co-Authored-By: <your model name> <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_015Xtq7GodhTGPYBzHDp3o4z"
```

---

### Task 2: Whole-genome readiness of stages 1-2 (`chain.rs` merges per chromosome pair; `seed.rs` memory)

**Files:**
- Modify: `src/sedef/chain.rs` (two merge conditions; one test in its `mod tests`), `src/sedef/seed.rs` (`SeedParams`, `build_windows`, `find_candidates`; new `count_shard`; new `#[cfg(test)] mod tests`), `src/sedef/mod.rs` (changes list)

**Interfaces:**
- Consumes: Task 1's `chain` and `seed` modules.
- Produces:
  - `chain::chain_and_merge` and `merge_blocks` union chains and blocks only within one (chr_a, chr_b) pair. Signatures are unchanged.
  - `SeedParams` gains `pub count_shards: usize`.
  - `seed::build_windows`: every `WindowSketch::sketch` has `capacity() == len()`.
  - `pub fn find_candidates(wins: &mut [WindowSketch], p: &SeedParams) -> Vec<Candidate>`. It empties every `wins[i].sketch` and returns the candidates unsorted, in shard order.

Two red-green cycles: chain first, then seed.

- [ ] **Step 1: Write the failing chain test**

In `src/sedef/chain.rs`'s `mod tests`, after the test `phase_shifted_fragments_merge_into_one_block` (before the module's closing `}`), add:

```rust
    #[test]
    fn chains_and_blocks_of_different_chromosome_pairs_never_merge() {
        // vendored change (2026-10-13): one run of three cells at the same window coordinates on (a, b) and on (x, y)
        // gives two blocks, each with its own names (the sandbox merged them into one block named after whichever
        // chain its HashMap met first)
        let on = |ca: &str, cb: &str, s_a: u64, s_b: u64| PairRec {
            chr_a: ca.into(),
            s_a,
            e_a: s_a + 2000,
            chr_b: cb.into(),
            s_b,
            e_b: s_b + 2000,
            shared: 30,
        };
        let mut pairs = vec![];
        for i in 0..3u64 {
            pairs.push(on("a", "b", i * 2000, 20_000 + i * 2000));
            pairs.push(on("x", "y", i * 2000, 20_000 + i * 2000));
        }
        let (blocks, stats) = chain_and_merge(&pairs, &params());
        assert_eq!(stats.n_chains, 2);
        let mut got: Vec<(String, u64, u64, String, u64, u64)> =
            blocks.iter().map(|b| (b.chr_a.clone(), b.a0, b.a1, b.chr_b.clone(), b.b0, b.b1)).collect();
        got.sort();
        assert_eq!(
            got,
            vec![
                ("a".to_string(), 0, 6000, "b".to_string(), 20_000, 26_000),
                ("x".to_string(), 0, 6000, "y".to_string(), 20_000, 26_000)
            ]
        );
    }
```

- [ ] **Step 2: Run it and watch it fail**

```bash
cd /mnt/linuxdisk/home/juanfraitu/rustle_m2_soto
RLOCK_TIMEOUT=900 bash tools/rlock.sh heavy env CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target TMPDIR=/mnt/linuxdisk/tmp \
  cargo test --release --lib chains_and_blocks > /mnt/linuxdisk/tmp/o1_layers/t2_red_chain.log 2>&1; echo "exit $?"
grep -E "^test |panicked|left|right" /mnt/linuxdisk/tmp/o1_layers/t2_red_chain.log | head -6
```

Expected: `exit 101` and `chains_and_blocks_of_different_chromosome_pairs_never_merge ... FAILED`. The sandbox code merges the two runs into one block (verified while writing this plan).

- [ ] **Step 3: Fix the two merge conditions**

In `chain_and_merge` (pass 2), replace:

```rust
        for &j in &active {
            if overlap_frac(chains[i].a0, chains[i].a1, chains[j].a0, chains[j].a1)
                >= p.merge_frac
```

with:

```rust
        for &j in &active {
            if chains[i].ca == chains[j].ca
                && chains[i].cb == chains[j].cb
                && overlap_frac(chains[i].a0, chains[i].a1, chains[j].a0, chains[j].a1)
                    >= p.merge_frac
```

In `merge_blocks`, replace:

```rust
        for &j in &active {
            if overlap_frac64(blocks[i].a0, blocks[i].a1, blocks[j].a0, blocks[j].a1)
                >= merge_frac
```

with:

```rust
        for &j in &active {
            if blocks[i].chr_a == blocks[j].chr_a
                && blocks[i].chr_b == blocks[j].chr_b
                && overlap_frac64(blocks[i].a0, blocks[i].a1, blocks[j].a0, blocks[j].a1)
                    >= merge_frac
```

Both sweeps stay correct. They visit every earlier entry whose span still reaches the current start, and the name test only removes unions between chromosome pairs.

- [ ] **Step 4: Run it and watch it pass**

```bash
RLOCK_TIMEOUT=900 bash tools/rlock.sh heavy env CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target TMPDIR=/mnt/linuxdisk/tmp \
  cargo test --release --lib sedef::chain > /mnt/linuxdisk/tmp/o1_layers/t2_green_chain.log 2>&1; echo "exit $?"
grep "test result" /mnt/linuxdisk/tmp/o1_layers/t2_green_chain.log
```

Expected: `exit 0`, `test result: ok. 5 passed`.

- [ ] **Step 5: Write the failing seed tests** (append to `src/sedef/seed.rs`)

```rust
#[cfg(test)]
mod tests {
    use super::*;

    fn params(shards: usize) -> SeedParams {
        SeedParams { k: 16, window: 2000, sketch: 8, min_kmers: 1, max_postings: 3, min_shared: 2, count_shards: shards }
    }

    fn win(id: u64, sketch: &[u64]) -> WindowSketch {
        WindowSketch { seq: 0, start: id * 2000, end: id * 2000 + 2000, sketch: sketch.to_vec() }
    }

    /// Hash 9 sits in 4 windows (> max_postings 3) and pairs nothing. Shared counts: (0,1) = {1,2,3} = 3;
    /// (0,2) = {3} = 1; (1,2) = {3,4} = 2; (3,4) = {5} = 1. With min_shared 2: (0,1) and (1,2).
    fn windows() -> Vec<WindowSketch> {
        vec![win(0, &[1, 2, 3, 9]), win(1, &[1, 2, 3, 4, 9]), win(2, &[3, 4, 9]), win(3, &[5, 9]), win(4, &[5, 6])]
    }

    fn sorted(c: Vec<Candidate>) -> Vec<(usize, usize, u32)> {
        let mut v: Vec<(usize, usize, u32)> = c.into_iter().map(|c| (c.a, c.b, c.shared)).collect();
        v.sort();
        v
    }

    #[test]
    fn sharded_counting_gives_the_one_map_candidates() {
        let one = sorted(find_candidates(&mut windows(), &params(1)));
        assert_eq!(one, vec![(0, 1, 3), (1, 2, 2)]);
        for s in [2, 3, 7, 64, 256] {
            assert_eq!(sorted(find_candidates(&mut windows(), &params(s))), one, "{s} shards");
        }
    }

    #[test]
    fn counting_releases_every_sketch() {
        let mut w = windows();
        find_candidates(&mut w, &params(4));
        assert!(w.iter().all(|x| x.sketch.is_empty()));
    }

    #[test]
    fn window_sketches_keep_only_their_hashes() {
        // vendored change (2026-10-13): the sandbox kept each window's 2,000-hash buffer (16 KB) behind a 256-hash
        // sketch, ~25 GB on a human genome
        let mut x = 7u64;
        let seq: Vec<u8> = (0..20_000)
            .map(|_| {
                x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                b"ACGT"[(x >> 62) as usize]
            })
            .collect();
        let p = SeedParams { k: 16, window: 2000, sketch: 256, min_kmers: 32, max_postings: 100, min_shared: 40, count_shards: 1 };
        let (wins, _) = build_windows(&[("s".to_string(), seq)], &p);
        assert_eq!(wins.len(), 10);
        assert!(wins.iter().all(|w| w.sketch.len() == 256 && w.sketch.capacity() == 256));
    }
}
```

- [ ] **Step 6: Run them and watch them fail**

```bash
RLOCK_TIMEOUT=900 bash tools/rlock.sh heavy env CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target TMPDIR=/mnt/linuxdisk/tmp \
  cargo test --release --lib sedef::seed > /mnt/linuxdisk/tmp/o1_layers/t2_red_seed.log 2>&1; echo "exit $?"
grep -E "^error" /mnt/linuxdisk/tmp/o1_layers/t2_red_seed.log | head -5
```

Expected: a non-zero exit, with `error[E0560]: struct `SeedParams` has no field named `count_shards`` (or E0063).

- [ ] **Step 7: Implement the seed changes**

In `SeedParams`, after `pub min_shared: usize,`, add:

```rust
    /// Vendored change (2026-10-13): pairs are counted in this many shards (a pair belongs to the shard of its lower
    /// window index modulo `count_shards`), so at most one shard map per thread is alive. Every value gives the same
    /// candidates; 1 is the sandbox's single map.
    pub count_shards: usize,
```

In `build_windows`, replace the line:

```rust
                hashes.truncate(p.sketch);
```

with:

```rust
                hashes.truncate(p.sketch);
                // vendored change: release the k-mer buffer's spare capacity (the sandbox kept 16 KB per window)
                hashes.shrink_to_fit();
```

Replace the whole `find_candidates` function (from its `/// Inverted index over sketch values + pair counting.` doc line to its closing brace) with:

```rust
/// Inverted index over sketch values + pair counting.
///
/// Implementation: flatten (hash, window) entries and parallel-sort -- group iteration then streams linearly
/// through memory instead of chasing one allocated Vec per hash, and `par_sort_unstable` uses all cores.
///
/// Vendored changes (2026-10-13): each window's sketch is released once its entries are written (nothing reads it
/// again), and the pairs are counted in `p.count_shards` independent shards, so a whole genome's pair map is never
/// held at once. Candidates come back in shard order; the caller sorts them.
pub fn find_candidates(wins: &mut [WindowSketch], p: &SeedParams) -> Vec<Candidate> {
    let mut entries: Vec<(u64, u32)> = Vec::with_capacity(wins.iter().map(|w| w.sketch.len()).sum());
    for (id, w) in wins.iter_mut().enumerate() {
        for &h in &w.sketch {
            entries.push((h, id as u32));
        }
        w.sketch = Vec::new();
    }
    entries.par_sort_unstable();
    let shards = p.count_shards.max(1) as u32;
    let per_shard: Vec<Vec<Candidate>> =
        (0..shards).into_par_iter().map(|s| count_shard(&entries, s, shards, p)).collect();
    per_shard.into_iter().flatten().collect()
}

/// The candidates whose lower window index is `shard` modulo `shards`. Inside one hash group the entries are sorted
/// by window index and a window holds a hash once (its sketch is deduplicated), so `a < b` whenever x < y: the
/// shard test can be made once per x.
fn count_shard(entries: &[(u64, u32)], shard: u32, shards: u32, p: &SeedParams) -> Vec<Candidate> {
    let mut counts: FxHashMap<(u32, u32), u32> =
        FxHashMap::with_capacity_and_hasher(entries.len() / 64 / shards as usize, Default::default());
    let n = entries.len();
    let mut i = 0;
    while i < n {
        let h = entries[i].0;
        let mut j = i + 1;
        while j < n && entries[j].0 == h {
            j += 1;
        }
        let m = j - i;
        if m >= 2 && m <= p.max_postings {
            for x in i..j {
                let a = entries[x].1;
                if a % shards != shard {
                    continue;
                }
                for y in (x + 1)..j {
                    let b = entries[y].1;
                    debug_assert!(a < b);
                    *counts.entry((a, b)).or_insert(0) += 1;
                }
            }
        }
        i = j;
    }
    counts
        .into_iter()
        .filter(|(_, c)| *c >= p.min_shared as u32)
        .map(|((a, b), c)| Candidate { a: a as usize, b: b as usize, shared: c })
        .collect()
}
```

`seed.rs`'s `use std::collections::HashMap;` may become unused. Leave the warning; keep the diff minimal.

In `src/sedef/mod.rs`, append to the changes list:

```rust
//! - `chain::chain_and_merge` (pass 2) and `merge_blocks` join chains / blocks only within one chromosome pair: the
//!   sandbox merged by window coordinates alone, fusing blocks of different chromosome pairs under the names its
//!   HashMap met first (test `chains_and_blocks_of_different_chromosome_pairs_never_merge`).
//! - `seed::build_windows` keeps only each sketch's hashes (the sandbox kept the 2,000-hash buffer, 16 KB per window),
//!   and `seed::find_candidates` counts pairs in `SeedParams::count_shards` shards, releasing each sketch as its entries
//!   are written: the same candidates, bounded memory (tests `window_sketches_keep_only_their_hashes`,
//!   `sharded_counting_gives_the_one_map_candidates`, `counting_releases_every_sketch`).
```

- [ ] **Step 8: Run the tests and watch them pass**

```bash
RLOCK_TIMEOUT=900 bash tools/rlock.sh heavy env CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target TMPDIR=/mnt/linuxdisk/tmp \
  cargo test --release --lib sedef:: > /mnt/linuxdisk/tmp/o1_layers/t2_green.log 2>&1; echo "exit $?"
grep -E "test result|FAILED|panicked" /mnt/linuxdisk/tmp/o1_layers/t2_green.log
```

Expected: `exit 0`, `test result: ok. 16 passed` (12 ported + 1 chain + 3 seed).

- [ ] **Step 9: Commit**

```bash
git add src/sedef/chain.rs src/sedef/seed.rs src/sedef/mod.rs
git commit -m "sedef: chain merges within a chromosome pair; seed keeps only sketches and counts pairs in shards

Sandbox chain fused blocks of different chromosome pairs (wrong, run-to-run different on multi-chromosome input);
sandbox seed kept a 16 KB buffer per window and one pair map (projected ~50 GB on a human genome).
Plan docs/superpowers/plans/2026-10-13-o1-layers-phase-1-2.md, Task 2.

Co-Authored-By: <your model name> <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_015Xtq7GodhTGPYBzHDp3o4z"
```

---

### Task 3: Duplication families link a segment's two copies; dedup order is total (`src/sedef/decompose.rs`)

**Files:**
- Modify: `src/sedef/decompose.rs` (`Family`, `families`, `dedup_segments`'s final sort; new `#[cfg(test)] mod tests`), `src/sedef/mod.rs` (changes list)

**Interfaces:**
- Consumes: Task 1's `decompose`.
- Produces:
  - `Family { id: usize, n_copies: usize, chrs: Vec<String>, span_bp: u64, copies: Vec<(String, u64, u64)> }`.
    - `copies` holds the family's distinct copies, sorted.
    - `n_copies` counts segment endpoints, as in the sandbox.
  - Ids `0..` follow the order: n_copies descending, span_bp descending, chrs ascending, then copies ascending.
  - `dedup_segments` returns segments ordered by anchors descending, then (chr_a, a0, chr_b, b0, a1, b1, strand).

- [ ] **Step 1: Write the failing tests** (append to `src/sedef/decompose.rs`)

```rust
#[cfg(test)]
mod tests {
    use super::*;

    fn seg(ca: &str, a0: u64, a1: u64, cb: &str, b0: u64, b1: u64, strand: char) -> Segment {
        Segment { chr_a: ca.into(), a0, a1, chr_b: cb.into(), b0, b1, strand, anchors: 100, est_id: 99.0 }
    }

    #[test]
    fn a_segment_joins_its_two_copies_and_overlap_joins_further_copies() {
        // S1 chrA:1000-3000 ~ chrB:5000-7000 (inverted); S2 chrB:5100-7000 ~ chrC:100-2000 (its chrB side overlaps
        // S1's by 1900 / 1900 >= 0.8); S3 chrA:50000-52000 ~ chrA:90000-92000 (a family of its own)
        let fams = families(&[
            seg("chrA", 1000, 3000, "chrB", 5000, 7000, '-'),
            seg("chrB", 5100, 7000, "chrC", 100, 2000, '+'),
            seg("chrA", 50000, 52000, "chrA", 90000, 92000, '+'),
        ]);
        assert_eq!(fams.len(), 2);
        assert_eq!((fams[0].id, fams[0].n_copies), (0, 4));
        assert_eq!(fams[0].chrs, vec!["chrA".to_string(), "chrB".to_string(), "chrC".to_string()]);
        let want: Vec<(String, u64, u64)> = vec![
            ("chrA".to_string(), 1000, 3000),
            ("chrB".to_string(), 5000, 7000),
            ("chrB".to_string(), 5100, 7000),
            ("chrC".to_string(), 100, 2000),
        ];
        assert_eq!(fams[0].copies, want);
        // span per chromosome: chrA 2000 + chrB 2000 (5000-7000 covers 5100-7000) + chrC 1900
        assert_eq!(fams[0].span_bp, 5900);
        assert_eq!((fams[1].id, fams[1].n_copies, fams[1].span_bp), (1, 2, 4000));
    }

    #[test]
    fn equal_families_are_ordered_by_their_copies() {
        // same copy count, span and chromosome: the family whose first copy comes first is F0, whatever the
        // HashMap's order (each call builds new maps with new random seeds)
        let segs = [seg("chrA", 9000, 11000, "chrA", 30000, 32000, '+'), seg("chrA", 1000, 3000, "chrA", 20000, 22000, '+')];
        for _ in 0..20 {
            let f = families(&segs);
            assert_eq!(f.len(), 2);
            assert_eq!(f[0].copies[0], ("chrA".to_string(), 1000, 3000));
            assert_eq!(f[1].copies[0], ("chrA".to_string(), 9000, 11000));
            assert_eq!(f.iter().map(|x| x.id).collect::<Vec<_>>(), vec![0, 1]);
        }
    }

    #[test]
    fn dedup_orders_tied_representatives_by_their_ends_and_strand() {
        // the same anchors and starts on opposite strands stay two segments (dedup joins one strand only); the
        // sandbox left their order to HashMap order
        for _ in 0..20 {
            let out = dedup_segments(vec![
                seg("chrA", 0, 2000, "chrB", 0, 2000, '-'),
                seg("chrA", 0, 2000, "chrB", 0, 2000, '+'),
            ]);
            assert_eq!(out.iter().map(|s| s.strand).collect::<Vec<_>>(), vec!['+', '-']);
        }
    }
}
```

- [ ] **Step 2: Run them and watch them fail**

```bash
RLOCK_TIMEOUT=900 bash tools/rlock.sh heavy env CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target TMPDIR=/mnt/linuxdisk/tmp \
  cargo test --release --lib sedef::decompose > /mnt/linuxdisk/tmp/o1_layers/t3_red.log 2>&1; echo "exit $?"
grep -E "^error" /mnt/linuxdisk/tmp/o1_layers/t3_red.log | head -5
```

Expected: a non-zero exit, with `error[E0609]: no field `copies` on type `Family``.

- [ ] **Step 3: Implement**

Replace `pub struct Family { ... }` with:

```rust
pub struct Family {
    pub id: usize,
    pub n_copies: usize,
    pub chrs: Vec<String>,
    pub span_bp: u64,
    /// Vendored change (2026-10-13): the family's distinct copies `(chrom, start, end)`, sorted.
    pub copies: Vec<(String, u64, u64)>,
}
```

Replace the whole `families` function (from `/// Group segment endpoints into duplication families` to its closing brace) with:

```rust
/// Group segment endpoints into duplication families: the two endpoints of a segment are one family, and nodes on
/// the same chromosome overlapping >= 80% (of the shorter; any strand) are one family. Connected components are
/// duplication families.
///
/// Vendored changes (2026-10-13): (1) the sandbox never joined a segment's two endpoints, so its families were
/// clusters of overlapping endpoints with every copy of a duplication in a family of its own; (2) `span_bp` sums the
/// covered bases per chromosome (the sandbox merged intervals of different chromosomes); (3) the distinct copies are
/// kept (`Family::copies`) and break the remaining ties of the order, so ids no longer follow HashMap order.
fn families(segs: &[Segment]) -> Vec<Family> {
    #[derive(Clone)]
    struct Node {
        chr: String,
        s: u64,
        e: u64,
    }
    let nodes: Vec<Node> = segs
        .iter()
        .flat_map(|s| {
            [
                Node { chr: s.chr_a.clone(), s: s.a0, e: s.a1 },
                Node { chr: s.chr_b.clone(), s: s.b0, e: s.b1 },
            ]
        })
        .collect();
    let n = nodes.len();
    let mut dsu = Dsu::new(n);
    for i in 0..segs.len() {
        dsu.union(2 * i, 2 * i + 1);
    }
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&i, &j| (&nodes[i].chr, nodes[i].s).cmp(&(&nodes[j].chr, nodes[j].s)));
    let mut active: Vec<usize> = Vec::new();
    for &i in &order {
        active.retain(|&j| nodes[j].e > nodes[i].s);
        for &j in &active {
            if nodes[i].chr == nodes[j].chr
                && overlap_frac(nodes[i].s, nodes[i].e, nodes[j].s, nodes[j].e) >= 0.8
            {
                dsu.union(i, j);
            }
        }
        active.push(i);
    }
    let mut agg: HashMap<usize, Vec<(String, u64, u64)>> = HashMap::new();
    for (i, nd) in nodes.iter().enumerate() {
        agg.entry(dsu.find(i)).or_default().push((nd.chr.clone(), nd.s, nd.e));
    }
    let mut out: Vec<Family> = agg
        .into_values()
        .map(|of_root| {
            let n_copies = of_root.len();
            let mut copies = of_root;
            copies.sort();
            copies.dedup();
            let mut chrs: Vec<String> = copies.iter().map(|c| c.0.clone()).collect();
            chrs.dedup();
            let span_bp = chrs
                .iter()
                .map(|c| {
                    let iv: Vec<(u64, u64)> = copies.iter().filter(|x| &x.0 == c).map(|x| (x.1, x.2)).collect();
                    crate::sedef::refine::covered_len(&iv)
                })
                .sum();
            Family { id: 0, n_copies, chrs, span_bp, copies }
        })
        .collect();
    out.sort_by(|x, y| {
        y.n_copies
            .cmp(&x.n_copies)
            .then_with(|| y.span_bp.cmp(&x.span_bp))
            .then_with(|| x.chrs.cmp(&y.chrs))
            .then_with(|| x.copies.cmp(&y.copies))
    });
    for (i, f) in out.iter_mut().enumerate() {
        f.id = i;
    }
    out
}
```

In `dedup_segments`, replace its final sort:

```rust
    out.sort_by(|x, y| {
        y.anchors
            .cmp(&x.anchors)
            .then_with(|| (&x.chr_a, x.a0, &x.chr_b, x.b0).cmp(&(&y.chr_a, y.a0, &y.chr_b, y.b0)))
    });
```

with:

```rust
    out.sort_by(|x, y| {
        y.anchors.cmp(&x.anchors).then_with(|| {
            (&x.chr_a, x.a0, &x.chr_b, x.b0, x.a1, x.b1, x.strand)
                .cmp(&(&y.chr_a, y.a0, &y.chr_b, y.b0, y.a1, y.b1, y.strand))
        })
    });
```

In `src/sedef/mod.rs`, append to the changes list:

```rust
//! - `decompose::dedup_segments` breaks ties of its order with (a1, b1, strand), so segments.tsv no longer follows
//!   HashMap order (test `dedup_orders_tied_representatives_by_their_ends_and_strand`).
//! - `decompose::families` joins the two endpoints of every segment (the module's documented intent; the sandbox
//!   joined only overlapping endpoints, so each copy was its own family), sums `span_bp` per chromosome, and keeps
//!   each family's distinct copies, which also break order ties
//!   (tests `a_segment_joins_its_two_copies_and_overlap_joins_further_copies`,
//!   `equal_families_are_ordered_by_their_copies`). families.tsv therefore differs from the sandbox's; chimp chr1:
//!   1,033 families (sandbox 3,316), largest 3,750 of 17,426 copies.
```

- [ ] **Step 4: Run the tests and watch them pass**

```bash
RLOCK_TIMEOUT=900 bash tools/rlock.sh heavy env CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target TMPDIR=/mnt/linuxdisk/tmp \
  cargo test --release --lib sedef:: > /mnt/linuxdisk/tmp/o1_layers/t3_green.log 2>&1; echo "exit $?"
grep -E "test result|FAILED|panicked" /mnt/linuxdisk/tmp/o1_layers/t3_green.log
```

Expected: `exit 0`, `test result: ok. 19 passed` (16 after Task 2, plus 3).

- [ ] **Step 5: Commit**

```bash
git add src/sedef/decompose.rs src/sedef/mod.rs
git commit -m "sedef: duplication families join each segment's two copies, keep their copies, order ties by copies; total dedup order

The sandbox's families never linked a segment's endpoints (each copy its own family), contrary to its own docs.
Plan docs/superpowers/plans/2026-10-13-o1-layers-phase-1-2.md, Task 3.

Co-Authored-By: <your model name> <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_015Xtq7GodhTGPYBzHDp3o4z"
```

---

### Task 4: `sedef::scan`, `utilities sd-scan`, the chimp chr1 equivalence and the memory gate

**Files:**
- Create: `src/sedef/scan.rs`, `tests/sd_scan.rs`, `bench/sedef/equivalence_chr1.sh`, `bench/sedef/memory_probe.sh`
- Modify: `src/sedef/mod.rs` (`pub mod scan;`, changes list), `src/utilities.rs` (module doc list, `Cmd::SdScan`, its match arm)

**Interfaces:**
- Consumes:
  - Task 2: `seed::find_candidates(&mut [WindowSketch], &SeedParams)` and `SeedParams::count_shards`;
  - Task 3: `decompose::Family::copies`;
  - Task 1: the signatures listed there.
- Produces (later tasks rely on these):
  - `rustle::sedef::scan::{ScanParams (Default), ScanParams::version_key(&self) -> String, ScanSummary, run(genome: &Path, out: &str, threads: usize, p: &ScanParams) -> anyhow::Result<ScanSummary>, read_genome(&Path) -> anyhow::Result<Vec<(String, Vec<u8>)>>, write_duplicons(&str, &[decompose::Family]) -> anyhow::Result<usize>, toy_genome_fasta(planted: bool) -> String}`.
  - Output files for `--out PREFIX`:
    - `PREFIX.sd.candidates.tsv`, `PREFIX.sd.blocks.tsv`;
    - `PREFIX.sd.pairs.tsv`: 11 columns, `#chr_a start_a end_a chr_b start_b end_b strand anchors est_id stage2_id status`;
    - `PREFIX.sd.anchors.bin`, `PREFIX.sd.segments.tsv`;
    - `PREFIX.sd.families.tsv`: `#family_id n_copies span_bp chromosomes`;
    - `PREFIX.sd.duplicons.tsv`: `#chrom start end family_id`, rows `chrom<TAB>start<TAB>end<TAB>F<id>`;
    - `PREFIX.sd.params.tsv`: `key<TAB>value`, including `version_key` and `bases`.
  - CLI: `utilities sd-scan --genome FA --out PREFIX [--threads 4]`, and `utilities sd-scan --version-key`.

- [ ] **Step 1: Write the failing unit tests, the fixture and a stub module**

Create `src/sedef/scan.rs` with only the fixture and the tests (the module's code comes in Step 4). Add `pub mod scan;` to `src/sedef/mod.rs` after `pub mod refine;`.

```rust
//! `utilities sd-scan` (Step 4 writes the module).

/// The test genome of this module's tests and of `tests/sd_scan.rs`. chrA and chrB, 200,000 bp each, from one LCG
/// (x = x * 6364136223846793005 + 1442695040888963407, base = "ACGT"[x >> 62], seed 7). With `planted`:
/// chrA:50,000-58,000, every 200th base changed (A->C->G->T->A), is copied inverted over chrB:120,000-128,000; then
/// chrB:10,000-15,000, every 100th base changed, is copied direct over chrB:170,000-175,000. Headers carry a
/// description after the name. Checked while writing the plan against the sandbox binary (plan Facts 3-4).
#[doc(hidden)]
pub fn toy_genome_fasta(planted: bool) -> String {
    fn seq(n: usize, x: &mut u64) -> Vec<u8> {
        (0..n)
            .map(|_| {
                *x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                b"ACGT"[(*x >> 62) as usize]
            })
            .collect()
    }
    fn mutate(s: &[u8], every: usize) -> Vec<u8> {
        s.iter()
            .enumerate()
            .map(|(i, &c)| {
                if i % every == every - 1 {
                    match c {
                        b'A' => b'C',
                        b'C' => b'G',
                        b'G' => b'T',
                        _ => b'A',
                    }
                } else {
                    c
                }
            })
            .collect()
    }
    fn rc(s: &[u8]) -> Vec<u8> {
        s.iter()
            .rev()
            .map(|&c| match c {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                _ => b'A',
            })
            .collect()
    }
    let mut x = 7u64;
    let a = seq(200_000, &mut x);
    let mut b = seq(200_000, &mut x);
    if planted {
        let inv = rc(&mutate(&a[50_000..58_000], 200));
        b.splice(120_000..128_000, inv);
        let dir = mutate(&b[10_000..15_000], 100);
        b.splice(170_000..175_000, dir);
    }
    let mut out = String::new();
    for (name, s) in [("chrA", &a), ("chrB", &b)] {
        out.push_str(&format!(">{name} toy genome\n"));
        for line in s.chunks(80) {
            out.push_str(std::str::from_utf8(line).unwrap());
            out.push('\n');
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use std::path::Path;

    const PLANTED_PAIRS: &str = "#chr_a\tstart_a\tend_a\tchr_b\tstart_b\tend_b\tstrand\tanchors\test_id\tstage2_id\tstatus\n\
chrB\t10000\t14999\tchrB\t170000\t174999\t+\t4200\t99.0\t97.8\tok\n\
chrA\t52000\t53999\tchrB\t124001\t126000\t-\t1840\t99.5\t99.5\tok\n\
chrA\t50000\t51999\tchrB\t126001\t128000\t-\t1840\t99.5\t99.4\tok\n\
chrA\t56000\t57999\tchrB\t120001\t122000\t-\t1840\t99.5\t99.4\tok\n\
chrA\t54000\t55999\tchrB\t122001\t124000\t-\t1840\t99.5\t99.2\tok\n";

    const PLANTED_FAMILIES: &str = "#family_id\tn_copies\tspan_bp\tchromosomes\n\
F0\t2\t9998\tchrB\n\
F1\t2\t3998\tchrA;chrB\n\
F2\t2\t3998\tchrA;chrB\n\
F3\t2\t3998\tchrA;chrB\n\
F4\t2\t3998\tchrA;chrB\n";

    const PLANTED_DUPLICONS: &str = "#chrom\tstart\tend\tfamily_id\n\
chrA\t50000\t51999\tF1\n\
chrA\t52000\t53999\tF2\n\
chrA\t54000\t55999\tF3\n\
chrA\t56000\t57999\tF4\n\
chrB\t10000\t14999\tF0\n\
chrB\t120001\t122000\tF4\n\
chrB\t122001\t124000\tF3\n\
chrB\t124001\t126000\tF2\n\
chrB\t126001\t128000\tF1\n\
chrB\t170000\t174999\tF0\n";

    fn tmp(name: &str) -> std::path::PathBuf {
        let d = std::env::temp_dir().join(format!("sd_scan_{name}_{}", std::process::id()));
        std::fs::create_dir_all(&d).unwrap();
        d
    }

    fn scan_toy(dir: &Path, name: &str, planted: bool) -> (ScanSummary, String) {
        let fa = dir.join(format!("{name}.fa"));
        std::fs::write(&fa, toy_genome_fasta(planted)).unwrap();
        let out = dir.join(name).to_str().unwrap().to_string();
        let s = run(&fa, &out, 2, &ScanParams::default()).expect("sd-scan");
        (s, out)
    }

    fn read(out: &str, ext: &str) -> String {
        std::fs::read_to_string(format!("{out}.sd.{ext}")).unwrap()
    }

    #[test]
    fn the_planted_duplications_come_out_as_pairs_and_two_copy_families() {
        let d = tmp("planted");
        let (s, out) = scan_toy(&d, "toy", true);
        assert_eq!(
            (s.sequences, s.bases, s.candidates, s.blocks, s.pairs, s.pairs_ok, s.segments, s.families, s.duplicon_rows),
            (2, 400_000, 7, 5, 5, 5, 5, 5, 10)
        );
        assert_eq!(read(&out, "pairs.tsv"), PLANTED_PAIRS);
        assert_eq!(read(&out, "families.tsv"), PLANTED_FAMILIES);
        assert_eq!(read(&out, "duplicons.tsv"), PLANTED_DUPLICONS);
        assert!(read(&out, "params.tsv").starts_with(&format!("version_key\t{}\n", ScanParams::default().version_key())));
    }

    #[test]
    fn two_runs_in_one_process_write_identical_tables() {
        let d = tmp("twice");
        let (_, a) = scan_toy(&d, "one", true);
        let (_, b) = scan_toy(&d, "two", true);
        for ext in ["candidates.tsv", "blocks.tsv", "pairs.tsv", "anchors.bin", "segments.tsv", "families.tsv", "duplicons.tsv"] {
            assert_eq!(
                std::fs::read(format!("{a}.sd.{ext}")).unwrap(),
                std::fs::read(format!("{b}.sd.{ext}")).unwrap(),
                "{ext}"
            );
        }
    }

    #[test]
    fn a_gzipped_genome_gives_the_same_tables() {
        let d = tmp("gz");
        let (_, plain) = scan_toy(&d, "plain", true);
        let gz = d.join("toy.fa.gz");
        let mut enc = flate2::write::GzEncoder::new(std::fs::File::create(&gz).unwrap(), flate2::Compression::default());
        enc.write_all(toy_genome_fasta(true).as_bytes()).unwrap();
        enc.finish().unwrap();
        let out = d.join("gz").to_str().unwrap().to_string();
        run(&gz, &out, 2, &ScanParams::default()).expect("sd-scan on .fa.gz");
        for ext in ["candidates.tsv", "pairs.tsv", "segments.tsv", "duplicons.tsv"] {
            assert_eq!(read(&plain, ext), read(&out, ext), "{ext}");
        }
    }

    #[test]
    fn a_genome_without_duplications_writes_header_only_tables() {
        let d = tmp("none");
        let (s, out) = scan_toy(&d, "none", false);
        assert_eq!((s.candidates, s.blocks, s.pairs, s.segments, s.families, s.duplicon_rows), (0, 0, 0, 0, 0, 0));
        for ext in ["candidates.tsv", "blocks.tsv", "pairs.tsv", "segments.tsv", "families.tsv", "duplicons.tsv"] {
            assert_eq!(read(&out, ext).lines().count(), 1, "{ext}: header only");
        }
    }

    #[test]
    fn records_are_named_by_the_first_word_and_a_repeated_name_is_refused() {
        let d = tmp("names");
        let fa = d.join("g.fa");
        std::fs::write(&fa, ">chr1 described here\nACGT\n>chr2\nAC\n").unwrap();
        let g = read_genome(&fa).unwrap();
        assert_eq!(g.iter().map(|(n, _)| n.as_str()).collect::<Vec<_>>(), vec!["chr1", "chr2"]);
        assert!(g.iter().all(|(_, s)| s.capacity() == s.len()), "sequences are shrunk to their length");
        std::fs::write(&fa, ">chr1 first\nACGT\n>chr1 second\nACGT\n").unwrap();
        let e = read_genome(&fa).unwrap_err().to_string();
        assert!(e.contains("two FASTA records are named chr1"), "{e}");
    }

    #[test]
    fn the_version_key_names_every_parameter_but_the_shard_count() {
        let p = ScanParams::default();
        let k = p.version_key();
        assert!(k.starts_with("sedef-rs 0.1.0 vendored 2026-10-13 scan v1;"), "{k}");
        assert_eq!(ScanParams { count_shards: 7, ..p.clone() }.version_key(), k);
        assert_ne!(ScanParams { min_seg: 999, ..p }.version_key(), k);
    }
}
```

- [ ] **Step 2: Write the failing CLI test** `tests/sd_scan.rs`

```rust
//! `utilities sd-scan` (src/sedef/scan.rs) from the command line: the toy genome's tables, `--version-key`, and the
//! arguments a scan needs (plan docs/superpowers/plans/2026-10-13-o1-layers-phase-1-2.md, Task 4).

use std::path::PathBuf;
use std::process::Command;

const BIN: &str = env!("CARGO_BIN_EXE_utilities");

fn scratch(name: &str) -> PathBuf {
    let d = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("sd_scan").join(name);
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

#[test]
fn sd_scan_writes_every_table_of_the_toy_genome() {
    let d = scratch("toy");
    let fa = d.join("toy.fa");
    std::fs::write(&fa, rustle::sedef::scan::toy_genome_fasta(true)).unwrap();
    let o = Command::new(BIN)
        .args(["sd-scan", "--genome"])
        .arg(&fa)
        .arg("--out")
        .arg(d.join("toy"))
        .args(["--threads", "2"])
        .output()
        .unwrap();
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    for ext in ["candidates.tsv", "blocks.tsv", "pairs.tsv", "anchors.bin", "segments.tsv", "families.tsv", "duplicons.tsv", "params.tsv"] {
        assert!(d.join(format!("toy.sd.{ext}")).exists(), "{ext}");
    }
    let params = std::fs::read_to_string(d.join("toy.sd.params.tsv")).unwrap();
    assert!(params.contains("\nfamilies\t5\n") && params.contains("\nduplicon_rows\t10\n"), "{params}");
}

#[test]
fn version_key_needs_no_genome_and_a_scan_needs_the_genome() {
    let o = Command::new(BIN).args(["sd-scan", "--version-key"]).output().unwrap();
    assert!(o.status.success());
    assert_eq!(
        String::from_utf8_lossy(&o.stdout).trim(),
        rustle::sedef::scan::ScanParams::default().version_key()
    );
    let o = Command::new(BIN).args(["sd-scan", "--out", "x"]).output().unwrap();
    assert!(!o.status.success());
    assert!(String::from_utf8_lossy(&o.stderr).contains("--genome"));
}
```

- [ ] **Step 3: Run both and watch them fail**

```bash
RLOCK_TIMEOUT=900 bash tools/rlock.sh heavy env CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target TMPDIR=/mnt/linuxdisk/tmp \
  cargo test --release --lib sedef::scan > /mnt/linuxdisk/tmp/o1_layers/t4_red.log 2>&1; echo "exit $?"
grep -E "^error" /mnt/linuxdisk/tmp/o1_layers/t4_red.log | head -5
```

Expected: a non-zero exit; `cannot find function `run`` / `cannot find struct ... `ScanParams``.

- [ ] **Step 4: Implement `src/sedef/scan.rs`**

Replace the module's first line (`//! `utilities sd-scan` (Step 4 writes the module).`) with the code below. Keep `toy_genome_fasta` and the tests module as written in Step 1.

```rust
//! `utilities sd-scan`: the four sedef-rs stages in one call, through files, on a thread pool of its own
//! (docs/superpowers/specs/2026-10-13-o1-layered-hierarchy-design.md, section 1). For `--out PREFIX`:
//!
//! - `PREFIX.sd.candidates.tsv`  stage 1 (seed): candidate window pairs
//! - `PREFIX.sd.blocks.tsv`      stage 2 (chain): SD blocks
//! - `PREFIX.sd.pairs.tsv`       stage 3 (refine): bp-exact SD pairs, strand, identity (percent), status ok/self/weak
//! - `PREFIX.sd.anchors.bin`     stage 3: the anchors stage 4 reads
//! - `PREFIX.sd.segments.tsv`    stage 4 (decompose): elementary segment pairs
//! - `PREFIX.sd.families.tsv`    stage 4: duplication families
//! - `PREFIX.sd.duplicons.tsv`   one row per distinct family copy, `chrom start end F<id>`: the `--hierarchy-duplicons`
//!   format (`family_hierarchy::DupIndex::from_text`)
//! - `PREFIX.sd.params.tsv`      the version key, the inputs, the counts and the seconds per stage
//!
//! With the default parameters (the sandbox CLI's), the candidates, blocks, pairs and segments equal the sandbox's
//! `seed`, `chain`, `refine --anchors-out` and `decompose --anchors` outputs byte for byte, except that a sequence is
//! named by the first word of its FASTA header (the sandbox kept the whole header line).

use anyhow::{anyhow, bail, Context, Result};
use std::collections::{HashMap, HashSet};
use std::io::Write;
use std::path::Path;
use std::time::Instant;

use crate::sedef::{anchors, chain, decompose, fasta, refine, seed};

/// Bumped whenever an output of `run` can change for the same genome and parameters; part of the version key.
pub const SCAN_VERSION: u32 = 1;

/// The sandbox CLI's default parameters (`sedef_rs seed|chain|refine|decompose --help`), plus the output-neutral
/// shard count of the seed stage's pair counting.
#[derive(Clone, Debug, PartialEq)]
pub struct ScanParams {
    pub k: usize,
    pub window: usize,
    pub sketch: usize,
    pub min_kmers: usize,
    pub max_postings: usize,
    pub min_shared: usize,
    pub chain_gap: u32,
    pub chain_min_score: u32,
    pub chain_merge_frac: f64,
    pub min_cell_id: f64,
    pub min_block_id: f64,
    pub refine_max_gap: u64,
    pub refine_diag_bin: i64,
    pub refine_min_span: u64,
    pub decompose_max_gap: u64,
    pub decompose_diag_bin: i64,
    pub min_seg: u64,
    pub max_pairs: u64,
    pub count_shards: usize,
}

impl Default for ScanParams {
    fn default() -> Self {
        ScanParams {
            k: 16,
            window: 2000,
            sketch: 256,
            min_kmers: 32,
            max_postings: 100,
            min_shared: 40,
            chain_gap: 2,
            chain_min_score: 20,
            chain_merge_frac: 0.5,
            min_cell_id: 88.0,
            min_block_id: 90.0,
            refine_max_gap: 400,
            refine_diag_bin: 50,
            refine_min_span: 500,
            decompose_max_gap: 200,
            decompose_diag_bin: 50,
            min_seg: 1000,
            max_pairs: 1024,
            count_shards: 256,
        }
    }
}

impl ScanParams {
    /// Everything that can change an output, on one line; the `sd` stage hashes it into its cache key.
    /// `count_shards` is left out: it changes no output.
    pub fn version_key(&self) -> String {
        format!(
            "sedef-rs 0.1.0 vendored 2026-10-13 scan v{SCAN_VERSION}; seed k={} window={} sketch={} min_kmers={} \
             max_postings={} min_shared={}; chain gap={} min_score={} merge_frac={} min_cell_id={} min_block_id={}; \
             refine max_gap={} diag_bin={} min_span={}; decompose max_gap={} diag_bin={} min_seg={}; max_pairs={}; \
             names=first-word",
            self.k,
            self.window,
            self.sketch,
            self.min_kmers,
            self.max_postings,
            self.min_shared,
            self.chain_gap,
            self.chain_min_score,
            self.chain_merge_frac,
            self.min_cell_id,
            self.min_block_id,
            self.refine_max_gap,
            self.refine_diag_bin,
            self.refine_min_span,
            self.decompose_max_gap,
            self.decompose_diag_bin,
            self.min_seg,
            self.max_pairs
        )
    }
}

/// What `run` wrote, counted.
#[derive(Clone, Debug, PartialEq)]
pub struct ScanSummary {
    pub sequences: usize,
    pub bases: u64,
    pub candidates: usize,
    pub blocks: usize,
    pub pairs: usize,
    pub pairs_ok: usize,
    pub segments: usize,
    pub families: usize,
    pub duplicon_rows: usize,
}

/// Every FASTA record (plain, gzip or bgzip), named by the first word of its header line, each sequence shrunk to its
/// length. A record without a name, two records with one name, or a file without records is refused.
pub fn read_genome(path: &Path) -> Result<Vec<(String, Vec<u8>)>> {
    let recs = fasta::read_fasta(path, None).map_err(|e| anyhow!("reading {}: {e}", path.display()))?;
    let mut seen = HashSet::new();
    let mut out = Vec::with_capacity(recs.len());
    for (header, seq) in recs {
        let name = header.split_whitespace().next().unwrap_or("").to_string();
        if name.is_empty() {
            bail!("{}: a FASTA record has no name", path.display());
        }
        if !seen.insert(name.clone()) {
            bail!("{}: two FASTA records are named {name}", path.display());
        }
        let mut seq = seq;
        seq.shrink_to_fit(); // the reader grows each sequence by doubling: up to 2x its length
        out.push((name, seq));
    }
    if out.is_empty() {
        bail!("{}: no FASTA records", path.display());
    }
    Ok(out)
}

/// The four stages on a pool of `threads` threads of its own (the sandbox's `build_global` can run once per process).
pub fn run(genome: &Path, out: &str, threads: usize, p: &ScanParams) -> Result<ScanSummary> {
    let threads = threads.max(1);
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .context("building the sd-scan thread pool")?;
    pool.install(|| run_stages(genome, out, threads, p))
}

fn run_stages(genome: &Path, out: &str, threads: usize, p: &ScanParams) -> Result<ScanSummary> {
    let f = |ext: &str| format!("{out}.sd.{ext}");
    let mut secs: Vec<(&str, f64)> = Vec::new();

    let t = Instant::now();
    let (sequences, bases, candidates) = seed_stage(genome, &f("candidates.tsv"), p)?;
    secs.push(("seed", t.elapsed().as_secs_f64()));
    eprintln!("[sd-scan] seed: {sequences} sequences, {bases} bp, {candidates} candidate window pairs ({:.1} s)", secs[0].1);

    let t = Instant::now();
    let blocks = chain_stage(&f("candidates.tsv"), &f("blocks.tsv"), p)?;
    secs.push(("chain", t.elapsed().as_secs_f64()));
    eprintln!("[sd-scan] chain: {blocks} blocks ({:.1} s)", secs[1].1);

    let t = Instant::now();
    let (pairs, pairs_ok) = refine_stage(genome, &f("blocks.tsv"), &f("pairs.tsv"), &f("anchors.bin"), p)?;
    secs.push(("refine", t.elapsed().as_secs_f64()));
    eprintln!("[sd-scan] refine: {pairs} pairs, {pairs_ok} ok ({:.1} s)", secs[2].1);

    let t = Instant::now();
    let (segments, families, duplicon_rows) =
        decompose_stage(&f("anchors.bin"), &f("segments.tsv"), &f("families.tsv"), &f("duplicons.tsv"), p)?;
    secs.push(("decompose", t.elapsed().as_secs_f64()));
    eprintln!(
        "[sd-scan] decompose: {segments} segments, {families} families, {duplicon_rows} duplicon rows ({:.1} s)",
        secs[3].1
    );

    let s = ScanSummary { sequences, bases, candidates, blocks, pairs, pairs_ok, segments, families, duplicon_rows };
    write_params(&f("params.tsv"), genome, threads, p, &s, &secs)?;
    Ok(s)
}

fn create(path: &str) -> Result<std::io::BufWriter<std::fs::File>> {
    Ok(std::io::BufWriter::new(std::fs::File::create(path).with_context(|| format!("creating {path}"))?))
}

/// Stage 1, as the sandbox's `run_seed`; the genome is released before the pairs are counted.
/// Returns (sequences, bases, candidates).
fn seed_stage(genome: &Path, out: &str, p: &ScanParams) -> Result<(usize, u64, usize)> {
    let seqs = read_genome(genome)?;
    let sequences = seqs.len();
    let bases: u64 = seqs.iter().map(|(_, s)| s.len() as u64).sum();
    let sp = seed::SeedParams {
        k: p.k,
        window: p.window,
        sketch: p.sketch,
        min_kmers: p.min_kmers,
        max_postings: p.max_postings,
        min_shared: p.min_shared,
        count_shards: p.count_shards,
    };
    let (mut wins, infos) = seed::build_windows(&seqs, &sp);
    drop(seqs);
    let mut cands = seed::find_candidates(&mut wins, &sp);
    sort_candidates(&mut cands, &wins);
    let mut w = create(out)?;
    writeln!(w, "#chr_a\tstart_a\tend_a\tchr_b\tstart_b\tend_b\tshared\test_jaccard")?;
    for c in &cands {
        let (wa, wb) = (&wins[c.a], &wins[c.b]);
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}",
            infos[wa.seq as usize].name,
            wa.start,
            wa.end,
            infos[wb.seq as usize].name,
            wb.start,
            wb.end,
            c.shared,
            seed::est_jaccard(c.shared, sp.sketch),
        )?;
    }
    w.flush()?;
    Ok((sequences, bases, cands.len()))
}

/// The sandbox's candidate order (most shared first, then window a's sequence and start, then window b's start),
/// with window b's sequence as the last key: the sandbox left ties between sequences to its single map's order,
/// which the shards change. On a one-sequence genome there are no such ties.
pub fn sort_candidates(cands: &mut [seed::Candidate], wins: &[seed::WindowSketch]) {
    cands.sort_by(|x, y| {
        y.shared.cmp(&x.shared).then_with(|| {
            (wins[x.a].seq, wins[x.a].start, wins[x.b].start, wins[x.b].seq)
                .cmp(&(wins[y.a].seq, wins[y.a].start, wins[y.b].start, wins[y.b].seq))
        })
    });
}

/// Stage 2, as the sandbox's `run_chain`. Returns the number of blocks.
fn chain_stage(cands: &str, out: &str, p: &ScanParams) -> Result<usize> {
    let min_shared_cell = chain::shared_for_identity(p.min_cell_id / 100.0, p.sketch, p.k);
    let pairs = chain::load_pairs(Path::new(cands), min_shared_cell).map_err(|e| anyhow!("reading {cands}: {e}"))?;
    let cp = chain::ChainParams {
        window: p.window as u64,
        gap: p.chain_gap,
        min_score: p.chain_min_score,
        merge_frac: p.chain_merge_frac,
        sketch: p.sketch,
        k: p.k,
        min_cell_id: p.min_cell_id,
        min_block_id: p.min_block_id,
        min_shared_cell,
    };
    let (mut blocks, _) = chain::chain_and_merge(&pairs, &cp);
    blocks.sort_by(|x, y| {
        y.score
            .cmp(&x.score)
            .then_with(|| (&x.chr_a, x.a0, &x.chr_b, x.b0).cmp(&(&y.chr_a, y.a0, &y.chr_b, y.b0)))
    });
    let mut w = create(out)?;
    writeln!(w, "#chr_a\tstart_a\tend_a\tchr_b\tstart_b\tend_b\tcells\tscore\test_id")?;
    for b in &blocks {
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.1}",
            b.chr_a, b.a0, b.a1, b.chr_b, b.b0, b.b1, b.cells, b.score, b.est_id
        )?;
    }
    w.flush()?;
    Ok(blocks.len())
}

/// Stage 3, as the sandbox's `run_refine --anchors-out`. Returns (pairs, ok pairs).
fn refine_stage(genome: &Path, blocks_path: &str, out: &str, anchors_out: &str, p: &ScanParams) -> Result<(usize, usize)> {
    let blocks = chain::load_blocks(Path::new(blocks_path)).map_err(|e| anyhow!("reading {blocks_path}: {e}"))?;
    let seqs: HashMap<String, Vec<u8>> = read_genome(genome)?.into_iter().collect();
    let rp = refine::RefineParams {
        k: p.k,
        max_gap: p.refine_max_gap,
        diag_bin: p.refine_diag_bin,
        min_span: p.refine_min_span,
        max_pairs: p.max_pairs,
    };
    let refined_pairs = refine::refine_blocks(&blocks, &seqs, &rp);
    drop(seqs);
    let (refined, anchored): (Vec<refine::Refined>, Vec<anchors::BlockAnchors>) = refined_pairs.into_iter().unzip();
    anchors::write(Path::new(anchors_out), p.k, &anchored).map_err(|e| anyhow!("writing {anchors_out}: {e}"))?;
    drop(anchored);
    let mut w = create(out)?;
    writeln!(w, "#chr_a\tstart_a\tend_a\tchr_b\tstart_b\tend_b\tstrand\tanchors\test_id\tstage2_id\tstatus")?;
    let mut ok = 0;
    for r in &refined {
        let status = if !r.ok {
            "weak"
        } else if r.self_pair {
            "self"
        } else {
            ok += 1;
            "ok"
        };
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.1}\t{:.1}\t{}",
            r.chr_a, r.ra0, r.ra1, r.chr_b, r.rb0, r.rb1, r.strand, r.anchors, r.est_id, r.stage2_id, status
        )?;
    }
    w.flush()?;
    Ok((refined.len(), ok))
}

/// Stage 4, as the sandbox's `decompose --anchors`, plus the duplicon table. Returns (segments, families, rows).
fn decompose_stage(anchors_path: &str, seg_out: &str, fam_out: &str, dup_out: &str, p: &ScanParams) -> Result<(usize, usize, usize)> {
    let (k, recs) = anchors::read(Path::new(anchors_path)).map_err(|e| anyhow!("reading {anchors_path}: {e}"))?;
    let dp = decompose::DecomposeParams {
        k,
        max_gap: p.decompose_max_gap,
        diag_bin: p.decompose_diag_bin,
        min_seg: p.min_seg,
        max_pairs: p.max_pairs,
    };
    let (segs, fams, _) = decompose::decompose(recs, &dp);
    let mut w = create(seg_out)?;
    writeln!(w, "#chr_a\tstart_a\tend_a\tchr_b\tstart_b\tend_b\tstrand\tanchors\test_id")?;
    for s in &segs {
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.1}",
            s.chr_a, s.a0, s.a1, s.chr_b, s.b0, s.b1, s.strand, s.anchors, s.est_id
        )?;
    }
    w.flush()?;
    let mut w = create(fam_out)?;
    writeln!(w, "#family_id\tn_copies\tspan_bp\tchromosomes")?;
    for f in &fams {
        writeln!(w, "F{}\t{}\t{}\t{}", f.id, f.n_copies, f.span_bp, f.chrs.join(";"))?;
    }
    w.flush()?;
    let rows = write_duplicons(dup_out, &fams)?;
    Ok((segs.len(), fams.len(), rows))
}

/// One row per distinct copy of every family, `chrom start end F<id>`, sorted by chromosome, start, end and family:
/// the `--hierarchy-duplicons` input (`family_hierarchy::DupIndex::from_text` skips the `#` header).
pub fn write_duplicons(path: &str, fams: &[decompose::Family]) -> Result<usize> {
    let mut rows: Vec<(&str, u64, u64, usize)> = fams
        .iter()
        .flat_map(|f| f.copies.iter().map(move |(c, s, e)| (c.as_str(), *s, *e, f.id)))
        .collect();
    rows.sort();
    let mut w = create(path)?;
    writeln!(w, "#chrom\tstart\tend\tfamily_id")?;
    for (c, s, e, id) in &rows {
        writeln!(w, "{c}\t{s}\t{e}\tF{id}")?;
    }
    w.flush()?;
    Ok(rows.len())
}

fn write_params(path: &str, genome: &Path, threads: usize, p: &ScanParams, s: &ScanSummary, secs: &[(&str, f64)]) -> Result<()> {
    let mut w = create(path)?;
    writeln!(w, "version_key\t{}", p.version_key())?;
    writeln!(w, "genome\t{}", genome.display())?;
    let bytes = std::fs::metadata(genome).with_context(|| format!("reading {}", genome.display()))?.len();
    writeln!(w, "genome_bytes\t{bytes}")?;
    writeln!(w, "threads\t{threads}")?;
    writeln!(w, "count_shards\t{}", p.count_shards)?;
    for (k, v) in [
        ("sequences", s.sequences as u64),
        ("bases", s.bases),
        ("candidates", s.candidates as u64),
        ("blocks", s.blocks as u64),
        ("pairs", s.pairs as u64),
        ("pairs_ok", s.pairs_ok as u64),
        ("segments", s.segments as u64),
        ("families", s.families as u64),
        ("duplicon_rows", s.duplicon_rows as u64),
    ] {
        writeln!(w, "{k}\t{v}")?;
    }
    for (k, v) in secs {
        writeln!(w, "seconds_{k}\t{v:.1}")?;
    }
    w.flush()?;
    Ok(())
}
```

In `src/sedef/mod.rs`, append to the changes list:

```rust
//! - `scan` (new) replaces the sandbox's main.rs: the four stages through files on a local thread pool,
//!   first-word sequence names, candidates ordered with window b's sequence as the last key, and the duplicon table
//!   `PREFIX.sd.duplicons.tsv` (tests in `scan` and tests/sd_scan.rs; chimp chr1 equality:
//!   bench/sedef/equivalence_chr1.sh).
```

In `src/utilities.rs`:

1. In the module doc list (the `//!   utilities ...` lines at the top), after the `candidate-augment` lines, add:

```rust
//!   utilities sd-scan --genome FA --out PREFIX [--threads 4]     (or --version-key)
```

2. In `enum Cmd`, after the `CandidateAugment { ... },` variant, add:

```rust
    /// Segmental duplications of a genome with the vendored sedef-rs (`src/sedef/`; docs/superpowers/specs/
    /// 2026-10-13-o1-layered-hierarchy-design.md): seed -> chain -> refine -> decompose, written as
    /// OUT.sd.{candidates,blocks,pairs,segments,families,duplicons,params}.tsv + OUT.sd.anchors.bin, with the sandbox
    /// CLI's default parameters. `--version-key` prints the line the driver's `sd` stage keys its cache on, and exits.
    SdScan {
        /// Genome FASTA (plain, gzip or bgzip); sequences are named by the first word of their header.
        #[arg(long, required_unless_present = "version_key")]
        genome: Option<String>,
        /// Output prefix.
        #[arg(long, required_unless_present = "version_key")]
        out: Option<String>,
        #[arg(long, short = 't', default_value_t = 4)]
        threads: usize,
        #[arg(long)]
        version_key: bool,
    },
```

3. In `main`'s `match cli.cmd`, after the `Cmd::CandidateAugment { ... } => { ... }` arm, add:

```rust
        Cmd::SdScan { genome, out, threads, version_key } => {
            let p = rustle::sedef::scan::ScanParams::default();
            if version_key {
                println!("{}", p.version_key());
                return Ok(());
            }
            let (genome, out) = (genome.expect("clap: --genome"), out.expect("clap: --out"));
            let s = rustle::sedef::scan::run(std::path::Path::new(&genome), &out, threads, &p)?;
            eprintln!(
                "[sd-scan] {} sequences, {} bp: {} candidates, {} blocks, {} pairs ({} ok), {} segments, {} families, \
                 {} duplicon rows -> {out}.sd.*",
                s.sequences, s.bases, s.candidates, s.blocks, s.pairs, s.pairs_ok, s.segments, s.families, s.duplicon_rows
            );
            Ok(())
        }
```

- [ ] **Step 5: Run the tests and watch them pass**

```bash
cd /mnt/linuxdisk/home/juanfraitu/rustle_m2_soto
RLOCK_TIMEOUT=900 bash tools/rlock.sh heavy env CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target TMPDIR=/mnt/linuxdisk/tmp \
  cargo test --release --lib sedef:: > /mnt/linuxdisk/tmp/o1_layers/t4_green.log 2>&1; echo "exit $?"
grep -E "test result|FAILED|panicked" /mnt/linuxdisk/tmp/o1_layers/t4_green.log
RLOCK_TIMEOUT=900 bash tools/rlock.sh heavy env CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target TMPDIR=/mnt/linuxdisk/tmp \
  cargo test --release --test sd_scan > /mnt/linuxdisk/tmp/o1_layers/t4_cli.log 2>&1; echo "exit $?"
grep -E "test result|FAILED|panicked" /mnt/linuxdisk/tmp/o1_layers/t4_cli.log
```

Expected: `exit 0` twice; `test result: ok. 25 passed` (lib `sedef::`), then `test result: ok. 2 passed`.

If `the_planted_duplications_...` fails on the pairs text, compare the toy FASTA with the plan's generator. The Python reference is `/mnt/linuxdisk/tmp/sedef_probe/toy/gen.py`; `python3 gen.py X.fa` writes the genome the expected values came from. `cmp` it against `toy_genome_fasta(true)` written to a file. Fix the fixture, never the expected text.

- [ ] **Step 6: Build the binary; check equality with the sandbox on chimp chr1**

Create `bench/sedef/equivalence_chr1.sh`:

```bash
#!/bin/bash
# `utilities sd-scan` vs the sandbox sedef-rs on chimp chr1 (plan docs/superpowers/plans/2026-10-13-o1-layers-phase-1-2.md,
# Task 4): candidates, blocks, pairs (= the sandbox's refined) and segments byte-identical after the header-name change
# ("chr1 mPanTro3-v2.1" -> "chr1"); families.tsv is compared by its size profile, since the vendored families join a
# segment's two copies (src/sedef/mod.rs).
#   equivalence_chr1.sh UTILITIES_BIN [SANDBOX_DATA_DIR] [OUT_DIR]
set -euo pipefail
BIN=$1; D=${2:-/mnt/c/Users/jfris/Desktop/sd/sedef-rs/data}; O=${3:-/mnt/linuxdisk/tmp/o1_layers/equiv_chr1}
mkdir -p "$O"
/usr/bin/time -f $'wall_s\t%e\npeak_rss_kb\t%M' -o "$O/time.tsv" "$BIN" sd-scan --genome "$D/chr1.fa" --out "$O/chr1" --threads 5 2> "$O/log"
rc=0
for pair in candidates:candidates blocks:blocks refined:pairs segments:segments; do
  ref=${pair%%:*}; ours=${pair#*:}
  if cmp -s <(sed 's/chr1 mPanTro3-v2.1/chr1/g' "$D/chr1.$ref.tsv") "$O/chr1.sd.$ours.tsv"; then echo "$ours identical"; else echo "$ours DIFFERS"; rc=1; fi
done
awk -F'\t' '!/^#/{n++; c[n]=$2; tot+=$2; sp+=$3; if($2==1)s++} END{asort(c); printf "families %d one-copy %d copies %d largest %d %d %d %d %d span_sum %d\n", n, s, tot, c[n], c[n-1], c[n-2], c[n-3], c[n-4], sp}' "$O/chr1.sd.families.tsv"
cat "$O/time.tsv"
exit $rc
```

Run:

```bash
cd /mnt/linuxdisk/home/juanfraitu/rustle_m2_soto
RLOCK_TIMEOUT=900 bash tools/rlock.sh heavy env CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target TMPDIR=/mnt/linuxdisk/tmp \
  cargo build --release --bin utilities > /mnt/linuxdisk/tmp/o1_layers/t4_build.log 2>&1; echo "build exit $?"
RLOCK_TIMEOUT=900 bash tools/rlock.sh heavy bash bench/sedef/equivalence_chr1.sh /mnt/linuxdisk/home/juanfraitu/rustle_target/release/utilities
```

Expected:
- `build exit 0`;
- `candidates identical`, `blocks identical`, `pairs identical`, `segments identical`;
- `families 1033 one-copy 0 copies 17426 largest 3750 3626 2734 1690 1178 span_sum 9210601`;
- `peak_rss_kb` about 760,000 (plan check: 758,880; the sandbox's seed alone took 3,145,772).

If any table DIFFERS, STOP and report the first differing line (`cmp` / `diff | head`). Do not adjust the code to match without a ruling.

- [ ] **Step 7: The memory gate (chimp chr1 and chr1-3, projection to whole genomes)**

Create `bench/sedef/memory_probe.sh`:

```bash
#!/bin/bash
# Peak memory of `utilities sd-scan` on chimp chr1 and chr1-3, and its power-law projection to the whole genomes
# (plan docs/superpowers/plans/2026-10-13-o1-layers-phase-1-2.md, Task 4 gate: <= 20 GB to go on).
#   memory_probe.sh UTILITIES_BIN [OUT_DIR]
set -euo pipefail
BIN=$1; O=${2:-/mnt/linuxdisk/tmp/o1_layers/memprobe}; D=/mnt/c/Users/jfris/Desktop/sd/sedef-rs/data
mkdir -p "$O"
[ -s "$O/chimp3.fa" ] || zcat "$D/mPanTro3-v2.1.fa.gz" | awk '/^>/{n++} n>3{exit} {print}' > "$O/chimp3.fa"
md5sum "$O/chimp3.fa"
for item in chr1:$D/chr1.fa chimp3:$O/chimp3.fa; do
  tag=${item%%:*}; fa=${item#*:}
  /usr/bin/time -f '%e %M' -o "$O/$tag.time" "$BIN" sd-scan --genome "$fa" --out "$O/$tag" --threads 5 2> "$O/$tag.log"
  echo "$tag wall_s peak_kb: $(cat "$O/$tag.time")"
done
# a second chr1-3 run must write the same tables (multi-sequence determinism, Fact 2)
"$BIN" sd-scan --genome "$O/chimp3.fa" --out "$O/chimp3_rerun" --threads 5 2> "$O/chimp3_rerun.log"
for ext in candidates.tsv blocks.tsv pairs.tsv segments.tsv families.tsv duplicons.tsv; do
  cmp -s "$O/chimp3.sd.$ext" "$O/chimp3_rerun.sd.$ext" && echo "chimp3 rerun $ext identical" || echo "chimp3 rerun $ext DIFFERS"
done
grep -E "^(blocks|pairs|pairs_ok|segments|families|duplicon_rows)\t" "$O/chimp3.sd.params.tsv" | tr '\n' ' '; echo
python3 - "$O" <<'EOF'
import math, sys
o = sys.argv[1]
pts = []
for tag in ("chr1", "chimp3"):
    bases = int(next(l.split("\t")[1] for l in open(f"{o}/{tag}.sd.params.tsv") if l.startswith("bases\t")))
    pts.append((bases, int(open(f"{o}/{tag}.time").read().split()[1])))
(x1, y1), (x2, y2) = pts
e = math.log(y2 / y1) / math.log(x2 / x1)
for g, size in (("CHM13 v2.0", 3_117_292_070), ("gorilla", 3_545_850_636)):
    print(f"projected peak {g}: {y2 * (size / x2) ** e / 1048576:.1f} GB (exponent {e:.2f})")
EOF
```

Run, and compare the chr1-3 candidates with the sandbox's (Fact 2) as sets:

```bash
cd /mnt/linuxdisk/home/juanfraitu/rustle_m2_soto
RLOCK_TIMEOUT=900 bash tools/rlock.sh heavy bash bench/sedef/memory_probe.sh /mnt/linuxdisk/home/juanfraitu/rustle_target/release/utilities
diff <(sed 's/ Pan troglodytes[^\t]*//g' /mnt/linuxdisk/tmp/sedef_probe/cand.c3.tsv | sort) <(sort /mnt/linuxdisk/tmp/o1_layers/memprobe/chimp3.sd.candidates.tsv) > /dev/null \
  && echo "chimp3 candidates: the sandbox's set" || echo "chimp3 candidates DIFFER"
cmp -s <(sed 's/ Pan troglodytes[^\t]*//g' /mnt/linuxdisk/tmp/sedef_probe/cand.c3.tsv) /mnt/linuxdisk/tmp/o1_layers/memprobe/chimp3.sd.candidates.tsv \
  && echo "chimp3 candidates: the sandbox's order too" || echo "chimp3 candidates: same set, order differs only in cross-sequence ties (ruling R2)"
```

Expected (plan check, Fact 10):
- chimp3.fa md5 7a1f8b2fb4ae4faf71b6aca2d836af10;
- six `chimp3 rerun ... identical` lines;
- blocks 14558, pairs 14558, pairs_ok 14262, segments 12554, families 1568, duplicon_rows 24466;
- chr1 peak about 760,000 KB and chimp3 about 1,800,000-1,900,000 KB;
- projected peaks about 7-8 GB (CHM13 v2.0) and 8-9 GB (gorilla);
- `chimp3 candidates: the sandbox's set` and `... the sandbox's order too`.

**Gate:** STOP and report BLOCKED with the numbers if a projection exceeds 20 GB, the candidate sets differ, or a rerun table differs. The controller rules (e.g. more shards) before any whole-genome run.

- [ ] **Step 8: Commit**

```bash
chmod +x bench/sedef/equivalence_chr1.sh bench/sedef/memory_probe.sh
git add src/sedef/scan.rs src/sedef/mod.rs src/utilities.rs tests/sd_scan.rs bench/sedef/equivalence_chr1.sh bench/sedef/memory_probe.sh
git commit -m "sedef: utilities sd-scan (four stages, local pool, duplicon table); equal to the sandbox on chimp chr1

<paste the four 'identical' lines, the families line, peak_rss_kb, and the two projections here>
Plan docs/superpowers/plans/2026-10-13-o1-layers-phase-1-2.md, Task 4.

Co-Authored-By: <your model name> <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_015Xtq7GodhTGPYBzHDp3o4z"
```

---

### Task 5: The driver's `sd` stage (cached per genome)

**Files:**
- Create: `tests/pipeline_driver_sd.rs`
- Modify: `tools/rustle_pipeline.sh`:
  - a header comment block;
  - a usage line;
  - `stage_sd()` and a `case "$STAGE"` block, both placed just above the line `if [ "$STAGE" = index ]; then`.

**Interfaces:**
- Consumes: Task 4's `utilities sd-scan --genome --out --threads` and `utilities sd-scan --version-key`.
- Produces:
  - `tools/rustle_pipeline.sh sd --fasta G --out PREFIX [--threads N] [--bin DIR]`;
  - the cache directory `${RUSTLE_SD_CACHE_DIR:-PREFIX.cache}/sd/<md5 of G>-<first 8 hex of md5(version key)>/`, holding `genome.sd.*`, `key.tsv` (`sd<TAB>md5<TAB>version key<TAB>G`) and `DONE`;
  - `PREFIX.sd.prefix`, holding the line `<cache dir>/genome` (what `--sd-cache` takes);
  - `PREFIX.sd.log`.

- [ ] **Step 1: Write the failing driver tests** `tests/pipeline_driver_sd.rs`

```rust
//! `tools/rustle_pipeline.sh sd` (plan docs/superpowers/plans/2026-10-13-o1-layers-phase-1-2.md, Task 5) against a
//! STUB `utilities` that records each sd-scan call and writes the two tables the stage summarises: what is tested is
//! what the driver decides (the cache key, reuse, the prefix file, the arguments the stage needs).

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const DRIVER: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tools/rustle_pipeline.sh");

fn scratch(name: &str) -> PathBuf {
    let d = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("pipeline_driver_sd").join(name);
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(d.join("bin")).unwrap();
    let stub = d.join("bin/utilities");
    std::fs::write(
        &stub,
        "#!/bin/bash\n\
         if [ \"$1 $2\" = 'sd-scan --version-key' ]; then echo \"${STUB_KEY:-stub key 1}\"; exit 0; fi\n\
         echo \"$*\" >> \"$(dirname \"$0\")/calls.txt\"\n\
         [ -n \"${STUB_FAIL:-}\" ] && exit 3\n\
         while [ $# -gt 0 ]; do case \"$1\" in --out) OUT=$2; shift 2;; *) shift;; esac; done\n\
         printf '#h\\nchr1\\t0\\t2000\\tchr2\\t0\\t2000\\t+\\t9\\t99.0\\t99.0\\tok\\n' > \"$OUT.sd.pairs.tsv\"\n\
         printf '#h\\nF0\\t2\\t4000\\tchr1;chr2\\n' > \"$OUT.sd.families.tsv\"\n",
    )
    .unwrap();
    Command::new("chmod").args(["+x", stub.to_str().unwrap()]).status().unwrap();
    std::fs::write(d.join("g.fa"), ">chr1\nACGT\n").unwrap();
    d
}

fn sd(dir: &Path, env: &[(&str, &str)]) -> Output {
    let mut c = Command::new("bash");
    c.arg(DRIVER)
        .args(["sd", "--fasta"])
        .arg(dir.join("g.fa"))
        .arg("--out")
        .arg(dir.join("run"))
        .arg("--bin")
        .arg(dir.join("bin"))
        .env_remove("RUSTLE_SD_CACHE_DIR")
        .env_remove("STUB_KEY")
        .env_remove("STUB_FAIL")
        .env("TMPDIR", dir);
    for (k, v) in env {
        c.env(k, v);
    }
    c.output().expect("bash failed to spawn")
}

fn calls(dir: &Path) -> Vec<String> {
    std::fs::read_to_string(dir.join("bin/calls.txt")).unwrap_or_default().lines().map(String::from).collect()
}

fn prefix(dir: &Path) -> String {
    std::fs::read_to_string(dir.join("run.sd.prefix")).unwrap().trim().to_string()
}

fn err(o: &Output) -> String {
    String::from_utf8_lossy(&o.stderr).to_string()
}

/// `<cache>/sd/<32 hex>-<8 hex>/genome`
fn assert_key_dir(p: &str, cache: &Path) {
    let rest = p
        .strip_prefix(&format!("{}/sd/", cache.display()))
        .unwrap_or_else(|| panic!("{p} is not under {}/sd/", cache.display()));
    let (key, tail) = rest.split_once('/').unwrap();
    assert_eq!(tail, "genome");
    let (md5, ver) = key.split_once('-').unwrap();
    assert!(md5.len() == 32 && ver.len() == 8, "{key}");
    assert!(md5.chars().chain(ver.chars()).all(|c| c.is_ascii_hexdigit()), "{key}");
}

#[test]
fn sd_runs_sd_scan_once_per_genome_without_a_bam() {
    let d = scratch("once");
    let o = sd(&d, &[]);
    assert!(o.status.success(), "{}", err(&o));
    let c = calls(&d);
    assert_eq!(c.len(), 1);
    let p = prefix(&d);
    assert_key_dir(&p, &d.join("run.cache"));
    assert_eq!(c[0], format!("sd-scan --genome {} --out {p} --threads 4", d.join("g.fa").display()));
    let dir = Path::new(&p).parent().unwrap();
    assert!(dir.join("DONE").exists());
    assert!(std::fs::read_to_string(dir.join("key.tsv")).unwrap().starts_with("sd\t"));
    assert!(err(&o).contains("1 SD pairs (ok), 1 duplication families"), "{}", err(&o));
    let o = sd(&d, &[]);
    assert!(o.status.success(), "{}", err(&o));
    assert_eq!(calls(&d).len(), 1, "a cached genome is not scanned again");
    assert!(err(&o).contains("sd: cached"), "{}", err(&o));
    assert_eq!(prefix(&d), p);
}

#[test]
fn a_changed_genome_or_version_key_is_scanned_again() {
    let d = scratch("changed");
    assert!(sd(&d, &[]).status.success());
    let first = prefix(&d);
    std::fs::write(d.join("g.fa"), ">chr1\nACGTA\n").unwrap();
    assert!(sd(&d, &[]).status.success());
    let second = prefix(&d);
    assert_ne!(first, second);
    assert!(sd(&d, &[("STUB_KEY", "stub key 2")]).status.success());
    assert_ne!(prefix(&d), second);
    assert_eq!(calls(&d).len(), 3);
}

#[test]
fn a_failed_scan_is_not_cached() {
    let d = scratch("failed");
    let o = sd(&d, &[("STUB_FAIL", "1")]);
    assert!(!o.status.success());
    assert!(err(&o).contains("sd-scan failed"), "{}", err(&o));
    assert!(sd(&d, &[]).status.success());
    assert_eq!(calls(&d).len(), 2, "the failed scan left no DONE, so the next call scans");
}

#[test]
fn sd_needs_an_existing_fasta_and_takes_the_cache_dir_from_the_environment() {
    let d = scratch("args");
    let o = Command::new("bash").arg(DRIVER).args(["sd", "--out"]).arg(d.join("run")).arg("--bin").arg(d.join("bin")).output().unwrap();
    assert_eq!(o.status.code(), Some(2));
    assert!(err(&o).contains("sd needs --fasta and --out"), "{}", err(&o));
    let o = Command::new("bash")
        .arg(DRIVER)
        .args(["sd", "--fasta"])
        .arg(d.join("missing.fa"))
        .arg("--out")
        .arg(d.join("run"))
        .arg("--bin")
        .arg(d.join("bin"))
        .output()
        .unwrap();
    assert_eq!(o.status.code(), Some(2));
    assert!(err(&o).contains("is missing or empty"), "{}", err(&o));
    let cache = d.join("shared_cache");
    assert!(sd(&d, &[("RUSTLE_SD_CACHE_DIR", cache.to_str().unwrap())]).status.success());
    assert_key_dir(&prefix(&d), &cache);
}
```

- [ ] **Step 2: Run them and watch them fail**

```bash
RLOCK_TIMEOUT=900 bash tools/rlock.sh heavy env CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target TMPDIR=/mnt/linuxdisk/tmp \
  cargo test --release --test pipeline_driver_sd > /mnt/linuxdisk/tmp/o1_layers/t5_red.log 2>&1; echo "exit $?"
grep -E "test result|panicked" /mnt/linuxdisk/tmp/o1_layers/t5_red.log | head -6
```

Expected: a non-zero exit. The tests fail because today the driver answers `sd` with `need --bam, --fasta, --out` or `unknown stage sd`.

- [ ] **Step 3: Implement the stage**

In `tools/rustle_pipeline.sh`, insert immediately ABOVE the line `if [ "$STAGE" = index ]; then` (the shared `--bam` check):

```bash
# sd (docs/superpowers/specs/2026-10-13-o1-layered-hierarchy-design.md): the genome's segmental duplications with the
# vendored sedef-rs (`utilities sd-scan`), once per genome: cached in RUSTLE_SD_CACHE_DIR (default PREFIX.cache) under
# sd/<md5 of the FASTA>-<md5 of sd-scan's version key, 8 hex>/genome.sd.*, DONE written last. PREFIX.sd.prefix names
# the cached prefix that `mcl_families --emit-hierarchy --sd-cache` takes. Needs --fasta and --out only (no BAM), so it
# is a stage function plus a `case "$STAGE"` block that exits before the shared --bam check: figures/samples.py's
# driver_stage_code leaves every other stage's code hash as it was.
stage_sd() {
  [ -n "$FASTA" ] && [ -n "$OUT" ] || { echo "[rustle_pipeline] sd needs --fasta and --out" >&2; exit 2; }
  [ -s "$FASTA" ] || { echo "[rustle_pipeline] sd: $FASTA is missing or empty" >&2; exit 2; }
  local key md5 dir
  key=$("$BIN/utilities" sd-scan --version-key) || { echo "[rustle_pipeline] sd: $BIN/utilities has no sd-scan; rebuild it" >&2; exit 2; }
  md5=$(md5sum < "$FASTA" | cut -c1-32)
  dir=${RUSTLE_SD_CACHE_DIR:-$OUT.cache}/sd/$md5-$(printf '%s' "$key" | md5sum | cut -c1-8)
  if [ -f "$dir/DONE" ]; then
    echo "[rustle_pipeline] $(date +%H:%M:%S) sd: cached ($dir)" >&2
  else
    mkdir -p "$dir"
    echo "[rustle_pipeline] $(date +%H:%M:%S) sd: utilities sd-scan on $FASTA -> $dir" >&2
    "$BIN/utilities" sd-scan --genome "$FASTA" --out "$dir/genome" --threads "$THREADS" > "$OUT.sd.log" 2>&1 \
      || { echo "[rustle_pipeline] sd: sd-scan failed (exit $?); see $OUT.sd.log" >&2; exit 1; }
    printf 'sd\t%s\t%s\t%s\n' "$md5" "$key" "$FASTA" > "$dir/key.tsv"
    touch "$dir/DONE"
  fi
  printf '%s\n' "$dir/genome" > "$OUT.sd.prefix"
  echo "[rustle_pipeline] $(date +%H:%M:%S) sd: $(awk -F'\t' '!/^#/ && $11=="ok"' "$dir/genome.sd.pairs.tsv" | wc -l) SD pairs (ok), $(awk '!/^#/' "$dir/genome.sd.families.tsv" | wc -l) duplication families; --sd-cache $dir/genome ($OUT.sd.prefix)" >&2
}
case "$STAGE" in
  sd) stage_sd; exit 0;;
esac
```

In the header comment block, add these lines immediately above the line `# (thesis record: families/catalog = O1, candidates = O3, ...`:

```bash
#   sd        the genome's segmental duplications and duplicons            utilities sd-scan (vendored sedef-rs)
#             once per genome, cached by the FASTA's md5 and sd-scan's version key in RUSTLE_SD_CACHE_DIR (default
#             PREFIX.cache)/sd/<key>/genome.sd.*; PREFIX.sd.prefix names that prefix for `mcl_families --emit-hierarchy
#             --sd-cache`. Needs --fasta and --out only; `all` does not run it.
```

Then, after the line `#   tools/rustle_pipeline.sh cache-ls --out PREFIX      list what PREFIX.cache holds (cache-clear: delete it)`, add:

```bash
#   tools/rustle_pipeline.sh sd --fasta G --out PREFIX [--threads N] [--bin DIR]   the genome's SDs, cached (stage sd)
```

- [ ] **Step 4: Run the tests and watch them pass; check the other stages' code hashes**

```bash
cd /mnt/linuxdisk/home/juanfraitu/rustle_m2_soto
RLOCK_TIMEOUT=900 bash tools/rlock.sh heavy env CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target TMPDIR=/mnt/linuxdisk/tmp \
  cargo test --release --test pipeline_driver_sd --test pipeline_driver_units > /mnt/linuxdisk/tmp/o1_layers/t5_green.log 2>&1; echo "exit $?"
grep -E "test result|FAILED|panicked" /mnt/linuxdisk/tmp/o1_layers/t5_green.log
python3 - <<'EOF'
import sys
sys.path.insert(0, "figures")
import samples
want = {"assemble": "cf14d218d1cbe774", "families": "58c0422b24a99e8d", "candidates": "635739635b54c07a",
        "catalog": "4c987eb6e5839a46", "assign": "6016f8834bed63a8", "assign_fam": "e10889fd43f16d99",
        "merged": "f340f2721ccc3065", "flag": "9cff53027b6ef944", "index": "21eb18b2cc02dabe"}
got = {s: samples.driver_stage_code("tools/rustle_pipeline.sh", s) for s in want}
bad = {s: (want[s], got[s]) for s in want if want[s] != got[s]}
print("driver_stage_code unchanged for all 9 stages" if not bad else f"CHANGED: {bad}")
EOF
```

Expected:
- `exit 0`;
- `test result: ok. 4 passed` (pipeline_driver_sd) and the existing pipeline_driver_units count, all passed;
- `driver_stage_code unchanged for all 9 stages`.

If a hash changed, the edit touched shared code: move it into the function or the `case` block.

- [ ] **Step 5: Commit**

```bash
git add tools/rustle_pipeline.sh tests/pipeline_driver_sd.rs
git commit -m "driver: sd stage (utilities sd-scan once per genome, cached by FASTA md5 + sd-scan version key)

Other stages' driver_stage_code unchanged (function + case block only).
Plan docs/superpowers/plans/2026-10-13-o1-layers-phase-1-2.md, Task 5.

Co-Authored-By: <your model name> <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_015Xtq7GodhTGPYBzHDp3o4z"
```

---

### Task 6: The V1 comparison script (`bench/sedef/v1_compare.py`)

**Files:**
- Create: `bench/sedef/v1_compare.py`, `bench/sedef/test_v1_compare.py`

**Interfaces:**
- Consumes: `bench/o1_hier/v1_to_v2.py::read_chain(lines) -> {chrom: [(v1_start, v2_start, size)]}` (plan 1, unchanged).
- Produces: the CLI `v1_compare.py --sedef REF.bed --sedefrs PAIRS.tsv (--chain CHAIN | --fai FAI) [--min-id 0.90] [--min-len 1000] --out OUT.json`. It prints one JSON line with `recall`, `precision`, `ref_bp`, `test_bp`, `ref_pairs` and `test_pairs`, and writes the full JSON (adding `per_chrom`, `recall_bins`, `universe_bp` and `inputs`).

- [ ] **Step 1: Write the failing tests** `bench/sedef/test_v1_compare.py`

```python
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import v1_compare as V  # noqa: E402


def sedef_row(c1, s1, e1, c2, s2, e2, ident):
    f = [c1, str(s1), str(e1), c2, str(s2), str(e2)] + ["x"] * 14 + [str(ident)] + ["y"] * 13
    assert len(f) == 34
    return "\t".join(f) + "\n"


def rs_row(c1, s1, e1, c2, s2, e2, est, status):
    return "\t".join([c1, str(s1), str(e1), c2, str(s2), str(e2), "+", "100", str(est), str(est), status]) + "\n"


def test_sedef_rows_are_kept_by_identity_and_the_shorter_side():
    lines = ["#header\n", sedef_row("chr1", 0, 2000, "chr2", 0, 2000, 0.95), sedef_row("chr1", 0, 2000, "chr2", 0, 2000, 0.89),
             sedef_row("chr1", 0, 2000, "chr2", 0, 999, 0.99)]
    assert V.read_pairs(lines, "sedef", 0.90, 1000) == [(("chr1", 0, 2000), ("chr2", 0, 2000), 0.95)]


def test_sedefrs_rows_keep_ok_only_with_percent_identity():
    lines = ["#h\n", rs_row("chr1", 0, 2000, "chr2", 0, 2000, 97.6, "ok"), rs_row("chr1", 0, 2000, "chr1", 0, 2000, 100.0, "self"),
             rs_row("chr3", 0, 2000, "chr4", 0, 2000, 95.0, "weak"), rs_row("chr5", 0, 2000, "chr6", 0, 2000, 89.9, "ok")]
    p = V.read_pairs(lines, "sedefrs", 0.90, 1000)
    assert [(a, b) for a, b, _ in p] == [(("chr1", 0, 2000), ("chr2", 0, 2000))]
    assert abs(p[0][2] - 0.976) < 1e-12


def test_recall_and_precision_count_only_universe_bases():
    ref = [(("chr1", 0, 1000), ("chr1", 5000, 7000), 0.95)]
    test = [(("chr1", 500, 6000), ("chr2", 0, 1000), 0.95)]
    res, _ = V.compare(ref, test, {"chr1": [(0, 10_000)]})   # chr2 is outside the universe
    assert (res["ref_bp"], res["ref_covered"]) == (3000, 1500)
    assert (res["test_bp"], res["test_covered"]) == (5500, 1500)
    assert res["recall"] == 0.5
    assert abs(res["precision"] - 1500 / 5500) < 1e-12


def test_universe_from_chain_is_the_v2_side_of_the_aligned_blocks():
    chain = ["chain 100 chr1 1000 + 0 1000 chr1 1010 + 0 1010 1\n", "300 0 10\n", "700\n"]
    # v1 [0,300) -> v2 [0,300); 10 bases inserted on v2; v1 [300,1000) -> v2 [310,1010)
    assert V.universe_from_chain(chain) == {"chr1": [(0, 300), (310, 1010)]}


def test_recall_bins_split_by_identity_and_length():
    ref = [(("chr1", 0, 2000), ("chr1", 10_000, 12_000), 0.96), (("chr1", 20_000, 26_000), ("chr1", 30_000, 36_000), 0.995)]
    b = V.recall_bins(ref, {"chr1": [(0, 2000), (20_000, 23_000)]}, {"chr1": [(0, 100_000)]})
    assert b["identity"]["0.95-0.98"] == {"bp": 4000, "recall": 0.5}
    assert b["identity"]["0.99-1.0000001"] == {"bp": 12000, "recall": 0.25}
    assert b["length"]["1000-5000"]["recall"] == 0.5
    assert b["length"]["5000-20000"]["recall"] == 0.25
```

- [ ] **Step 2: Run them and watch them fail**

```bash
cd /mnt/linuxdisk/home/juanfraitu/rustle_m2_soto && python3 -m pytest -q bench/sedef/test_v1_compare.py 2>&1 | tail -3
```

Expected: `ModuleNotFoundError: No module named 'v1_compare'`.

- [ ] **Step 3: Implement** `bench/sedef/v1_compare.py`

```python
#!/usr/bin/env python3
"""V1 / V1b of docs/PREREG_sedef_rs_layers_2026-10-13.md: sedef-rs SD calls (`utilities sd-scan` PREFIX.sd.pairs.tsv)
against SEDEF's (34-column BED), base-pair-wise.

    v1_compare.py --sedef REF.bed --sedefrs PREFIX.sd.pairs.tsv (--chain V1_TO_V2.chain | --fai GENOME.fai)
                  [--min-id 0.90] [--min-len 1000] --out OUT.json

Kept pairs: SEDEF identity = column 21 (fracMatch) >= --min-id; sedef-rs status `ok` and est_id / 100 >= --min-id; both
tools: min(side lengths) >= --min-len. Each tool's SD regions = both sides of its kept pairs, merged per chromosome and
intersected with the UNIVERSE: the v2.0 bases with a v1.0 ancestor in the chain (--chain: CHM13 v2.0 against lifted
v1.0 SEDEF) or every base of the genome (--fai). recall = SEDEF region bases covered by sedef-rs regions / SEDEF region
bases; precision = sedef-rs region bases covered by SEDEF regions / sedef-rs region bases. Descriptive: per chromosome,
and the recall of the SEDEF pairs binned by identity and by length (the bases of each pair's two sides inside the
universe that sedef-rs regions cover; overlapping pairs count their shared bases more than once)."""
import argparse
import bisect
import collections
import json
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "o1_hier"))
from v1_to_v2 import read_chain  # noqa: E402

ID_BINS = [(0.90, 0.95), (0.95, 0.98), (0.98, 0.99), (0.99, 1.0000001)]
LEN_BINS = [(1_000, 5_000), (5_000, 20_000), (20_000, 100_000), (100_000, 10**12)]


def read_pairs(lines, fmt, min_id, min_len):
    """-> [((chrom, start, end), (chrom, start, end), identity)] kept by the filters; fmt 'sedef' or 'sedefrs'"""
    out = []
    for ln in lines:
        if not ln.strip() or ln.startswith("#"):
            continue
        f = ln.rstrip("\n").split("\t")
        if fmt == "sedef":
            ident = float(f[20])
        else:
            if f[10] != "ok":
                continue
            ident = float(f[8]) / 100.0
        a, b = (f[0], int(f[1]), int(f[2])), (f[3], int(f[4]), int(f[5]))
        if ident >= min_id and min(a[2] - a[1], b[2] - b[1]) >= min_len:
            out.append((a, b, ident))
    return out


def merge(iv):
    out = []
    for s, e in sorted(iv):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [tuple(x) for x in out]


def regions(pairs):
    by = collections.defaultdict(list)
    for a, b, _ in pairs:
        for c, s, e in (a, b):
            by[c].append((s, e))
    return {c: merge(v) for c, v in by.items()}


def intersect(x, y):
    """the bases in both merged, sorted interval lists, as a merged list"""
    out, i, j = [], 0, 0
    while i < len(x) and j < len(y):
        s, e = max(x[i][0], y[j][0]), min(x[i][1], y[j][1])
        if s < e:
            out.append((s, e))
        if x[i][1] < y[j][1]:
            i += 1
        else:
            j += 1
    return out


def total(iv):
    return sum(e - s for s, e in iv)


def covered_one(s, e, by, starts):
    """bases of [s, e) covered by the merged, sorted intervals `by` (`starts` = their starts)"""
    j = max(0, bisect.bisect_right(starts, s) - 1)
    tot = 0
    while j < len(by) and by[j][0] < e:
        tot += max(0, min(e, by[j][1]) - max(s, by[j][0]))
        j += 1
    return tot


def universe_from_chain(lines):
    """the v2 (query) side of every aligned block of the chain: {chrom: merged intervals}"""
    return {c: merge([(q, q + size) for _, q, size in blocks]) for c, blocks in read_chain(lines).items()}


def universe_from_fai(lines):
    return {f[0]: [(0, int(f[1]))] for f in (ln.split("\t") for ln in lines if ln.strip())}


def compare(ref_pairs, test_pairs, universe):
    """-> (result dict, the sedef-rs regions clipped to the universe)"""
    clip = lambda reg: {c: intersect(v, universe.get(c, [])) for c, v in reg.items()}  # noqa: E731
    r, t = clip(regions(ref_pairs)), clip(regions(test_pairs))
    per = {}
    for c in sorted(set(r) | set(t)):
        rv, tv = r.get(c, []), t.get(c, [])
        per[c] = dict(ref_bp=total(rv), test_bp=total(tv), ref_covered=total(intersect(rv, tv)),
                      test_covered=total(intersect(tv, rv)))
    rb, tb = sum(p["ref_bp"] for p in per.values()), sum(p["test_bp"] for p in per.values())
    rc, tc = sum(p["ref_covered"] for p in per.values()), sum(p["test_covered"] for p in per.values())
    return dict(recall=rc / rb if rb else None, precision=tc / tb if tb else None, ref_bp=rb, test_bp=tb,
                ref_covered=rc, test_covered=tc, per_chrom=per), t


def recall_bins(ref_pairs, test_regions, universe):
    """recall of the SEDEF pairs per identity bin and per length bin; `test_regions` are already inside the universe"""
    us = {c: [s for s, _ in v] for c, v in universe.items()}
    ts = {c: [s for s, _ in v] for c, v in test_regions.items()}
    acc = {"identity": collections.defaultdict(lambda: [0, 0]), "length": collections.defaultdict(lambda: [0, 0])}
    for a, b, ident in ref_pairs:
        length = min(a[2] - a[1], b[2] - b[1])
        for kind, value, bins in (("identity", ident, ID_BINS), ("length", length, LEN_BINS)):
            lab = next((f"{lo}-{hi}" for lo, hi in bins if lo <= value < hi), None)
            if lab is None:
                continue
            for c, s, e in (a, b):
                acc[kind][lab][0] += covered_one(s, e, universe.get(c, []), us.get(c, []))
                acc[kind][lab][1] += covered_one(s, e, test_regions.get(c, []), ts.get(c, []))
    return {k: {lab: dict(bp=v[0], recall=v[1] / v[0] if v[0] else None) for lab, v in sorted(d.items())}
            for k, d in acc.items()}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sedef", required=True)
    ap.add_argument("--sedefrs", required=True)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--chain")
    g.add_argument("--fai")
    ap.add_argument("--min-id", type=float, default=0.90)
    ap.add_argument("--min-len", type=int, default=1000)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    with open(a.sedef) as fh:
        ref = read_pairs(fh, "sedef", a.min_id, a.min_len)
    with open(a.sedefrs) as fh:
        test = read_pairs(fh, "sedefrs", a.min_id, a.min_len)
    with open(a.chain or a.fai) as fh:
        universe = universe_from_chain(fh) if a.chain else universe_from_fai(fh)
    res, test_regions = compare(ref, test, universe)
    res.update(inputs=vars(a), ref_pairs=len(ref), test_pairs=len(test),
               universe_bp=sum(total(v) for v in universe.values()),
               recall_bins=recall_bins(ref, test_regions, universe))
    with open(a.out, "w") as fh:
        json.dump(res, fh, indent=1)
    print(json.dumps({k: res[k] for k in ("recall", "precision", "ref_bp", "test_bp", "ref_pairs", "test_pairs")}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 4: Run the tests and watch them pass**

```bash
cd /mnt/linuxdisk/home/juanfraitu/rustle_m2_soto && python3 -m pytest -q bench/sedef/test_v1_compare.py bench/o1_hier 2>&1 | tail -2
```

Expected: `57 passed`: 5 here, plus plan 1's 52 in `bench/o1_hier`, unchanged.

- [ ] **Step 5: Commit**

```bash
git add bench/sedef/v1_compare.py bench/sedef/test_v1_compare.py
git commit -m "bench/sedef: V1 comparison (sd-scan vs SEDEF, bp recall/precision inside a universe, binned recall)

Plan docs/superpowers/plans/2026-10-13-o1-layers-phase-1-2.md, Task 6.

Co-Authored-By: <your model name> <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_015Xtq7GodhTGPYBzHDp3o4z"
```

---

### Task 7: Pre-registration (V1, V1b, V2), before any whole-genome comparison

**Files:**
- Create: `docs/PREREG_sedef_rs_layers_2026-10-13.md`

**Interfaces:**
- Consumes: the commit hashes of Tasks 4 (`sd-scan`) and 6 (`v1_compare.py`), and the version key printed by `utilities sd-scan --version-key`.
- Produces: the prereg that Tasks 8, 10 and 11 append their results to, under `## Results`.

- [ ] **Step 1: Write the prereg**

Fill the three `<...>` slots:
- `git log --format=%h -1 -- src/sedef/scan.rs`;
- `git log --format=%h -1 -- bench/sedef/v1_compare.py`;
- the output of `/mnt/linuxdisk/home/juanfraitu/rustle_target/release/utilities sd-scan --version-key`.

Keep everything else verbatim.

```markdown
# sedef-rs as the SD and duplicon layers: V1, V1b, V2 (2026-10-13)

**Pre-registered before any whole-genome sd-scan output is compared with SEDEF or DupMasker.** Spec:
`docs/superpowers/specs/2026-10-13-o1-layered-hierarchy-design.md`, section 2. Plan:
`docs/superpowers/plans/2026-10-13-o1-layers-phase-1-2.md`. Tool: `utilities sd-scan` (vendored sedef-rs, `src/sedef/`,
commit <scan commit>), default parameters only; version key: `<version key>`. V1 scorer: `bench/sedef/v1_compare.py`
(commit <v1_compare commit>).

**Not re-proposed.**
- Register 1263 (our own duplicon units, NOT ADOPTED at 0.638 held-out) is not re-run; V2 reuses its scorer and its
  split unchanged.
- SD or duplicon structure deciding family membership (register 676, 677/997, 817, 818) is not tried: sd-scan only
  annotates.

## Known before running (recorded now)

- Inverted SDs come out as window-sized pieces: on the toy genome, an 8 kb inverted copy gives four 2 kb pairs (`chain`
  links windows only along forward diagonals).
  - The base-pair union is unaffected except at a copy's ends.
  - The duplicons of an inverted SD are split at 2 kb.
- SD ends are clipped by up to one 2 kb window (sandbox README caveat 6), and refined pairs stop one base short.
  Short SDs (1-5 kb) are the most exposed (recall).
- sd-scan identities are k-mer lower bounds; at V1b's 0.98 cut this costs recall.
- sd-scan families are connected components of segment and >= 80%-overlap links. On chimp chr1 the largest holds 21.5% of
  all copies (3,750 of 17,426). Coarse duplicons in mosaic hubs are a risk for V2.
- The lifted CHM13 reference lost 925 of 151,758 SEDEF rows (rDNA rows and unliftable ends). Where their bases fall inside
  the universe, they count against sd-scan's precision.
- The vendored copy fixes two sandbox bugs that only multi-chromosome genomes show (`src/sedef/mod.rs`): chain merging
  across chromosome pairs, and family ids that followed HashMap order. Its runs are reproducible table for table
  (chimp chr1-3, two runs).

## V1 — SD calls, sd-scan vs SEDEF (both genomes)

**Inputs.**

| Genome | sd-scan run | SEDEF reference |
|---|---|---|
| CHM13 v2.0 | `chm13v2.0.fa` | `/mnt/linuxdisk/tmp/o1_close/lift/sedef_pairs.v2.bed` (v1.0 `final_v1_clean.bed`, lifted by the T2T chain; md5 7bd521edff9bd1e27b8141027da37a01) |
| Gorilla | `GGO.fasta` (mGorGor1 v2.0) | `winloci_data/GGO_sedef_final.bed` (md5 991ed34384832354c02df227e27f65b6) |

**Kept pairs.**
- SEDEF: column 21 (fracMatch) >= 0.90.
- sd-scan: status `ok` and est_id / 100 >= 0.90.
- Both: the shorter side >= 1,000 bp.

**Regions.** Both sides of every kept pair, merged per chromosome, then intersected with the universe:
- CHM13: the v2.0 bases inside the aligned blocks of `v1.0_to_v1.1_rdna_merged.chain`, i.e. bases with a v1.0
  ancestor (chr1-22 and X; not Y, M or the rDNA arrays);
- gorilla: every base of GGO.fasta.

**Measures.**
- recall = SEDEF region bases covered by sd-scan regions / SEDEF region bases;
- precision = sd-scan region bases covered by SEDEF regions / sd-scan region bases.

**Bar:** recall >= 0.80 AND precision >= 0.80, on CHM13 v2.0 AND on gorilla.

**Reading, fixed now.**
- **Pass:** sd-scan is the SD layer's source for every genome (the `sd` stage).
- **Fail:** sd-scan is not adopted. The SD layer keeps external SEDEF where it exists (human, gorilla); genomes without
  SEDEF get no SD layer.

## V1b — SD98

The same comparison, with 0.98 in place of 0.90 for both tools.

**Bar:** recall >= 0.80 AND precision >= 0.80, on both genomes.

**Reading.**
- **Pass:** the Soto layer (phase 3) takes its SD98 pairs from sd-scan.
- **Fail:** the Soto layer takes SEDEF SD98 where SEDEF exists, and alignment-based identity refinement in sedef-rs
  becomes its own task.

## Descriptive (V1 and V1b, no bar)

- Per chromosome: bases and coverage.
- Recall of the SEDEF pairs by identity bin and by length bin.
- Pair counts.
- Wall time and peak RSS of each sd-scan run (also the CHM13 v1.0 run, V2's input).

## V2 — duplicons, sd-scan vs DupMasker (human, F_O universe, CHM13 v1.0)

**Arms.** Both run with the current binary; families are not recomputed from scratch, only re-clustered from the same PAF.
- **DupMasker:**

      mcl_families --paf /mnt/linuxdisk/tmp/ourhom/loci.paf --gff /mnt/linuxdisk/tmp/ourhom/loci.gff3
        --min-exonic-bp 1 --min-shared-exon-frac 0.60 --dump-pairs --emit-hierarchy
        --hierarchy-duplicons winloci_data/duplicons/chm13.draft_v1.0_plus38Y_dupmasker_colors.bed
        --hierarchy-sd winloci_data/soto_replication/final_v1_clean.bed
        --hierarchy-external /mnt/linuxdisk/tmp/hier_flag/soto_external.tsv

  This is the registered F_O command (10-09), re-run.
- **sd-scan:** the same command, with `--hierarchy-duplicons <CHM13 v1.0 sd cache prefix>.sd.duplicons.tsv`. The SD
  layer stays SEDEF, so only the duplicons differ.

**H0 (guard).** `clusters.tsv`, `loci.tsv` and `pairs.tsv` are byte-identical between the arms.

**Scorer.** `bench/paf_duplicons/score_units.py`, unchanged (register 1263):
- the cover of a family (>= 2 members) = its expansions plus its `unique` class;
- family score = maximum-weight bipartite matching on Jaccard / max(#blocks);
- family-hash split: sha256 even = DEV, odd = HELD-OUT.
- Its `share_exonic_with_paf_unit` output is the sd-scan duplicons here.

**Bars.**
- **V2a:** mean family score over the HELD-OUT families >= 0.80.
- **V2b:** Soto families equal to an expansion strictly below a family, with sd-scan duplicons, >= ⌈0.9 × the DupMasker
  arm's count⌉ (10 on the 10-09 run, so >= 9).
  - Measured by `bench/o1_hier/hier_match.py MEMBERS soto_external.tsv OUT --fold LOCI --map key --null 0`,
    `exact_any.expansion_not_family`.
  - All Soto families; computed once, at the held-out step.

**Procedure.**
- DEV is scored and written first.
- At most ONE amendment may follow, registered here before HELD-OUT is computed.
- HELD-OUT decides, and V2b is computed with it.

**Reading, fixed now.**
- **V2a and V2b hold:** sd-scan duplicons replace DupMasker as the expansion layer's source, for every genome.
  DupMasker stays available for comparison runs.
- **Otherwise:** sd-scan gives the SD level only; human expansions stay on DupMasker; apes have no expansion layer.

**Descriptive (no bar).**
- P2: sign test against the one-block null on structured families.
- P3: Soto EXPANSION placements per arm.
- P4: the two named fusions.
- Exact expansion sensitivity and precision; blocks per family.
- The NPIP family (MCL9) cover.
- Members' exonic bases carrying a duplicon, per arm.
- DEV and HELD-OUT side by side.
- Whether the re-run DupMasker arm equals the registered 10-09 files.

## Results
```

- [ ] **Step 2: Commit** (nothing whole-genome may have been compared yet)

```bash
git add docs/PREREG_sedef_rs_layers_2026-10-13.md
git commit -m "prereg: sedef-rs as the SD and duplicon layers (V1, V1b vs SEDEF; V2 vs DupMasker on F_O), before any run

Plan docs/superpowers/plans/2026-10-13-o1-layers-phase-1-2.md, Task 7.

Co-Authored-By: <your model name> <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_015Xtq7GodhTGPYBzHDp3o4z"
```

---

### Task 8: Whole-genome `sd` runs (CHM13 v2.0, gorilla, CHM13 v1.0); V1 and V1b scored

**Files:**
- Create: `bench/sedef/run_genome.sh`
- Modify: `docs/PREREG_sedef_rs_layers_2026-10-13.md` (`## Results`), `docs/NEGATIVE_RESULTS_REGISTER.md` (two rows)

**Interfaces:**
- Consumes:
  - Task 5's `sd` stage;
  - Task 6's `v1_compare.py`;
  - Task 7's prereg, which must be committed: `git log --oneline -1 -- docs/PREREG_sedef_rs_layers_2026-10-13.md`.
- Produces: `/mnt/linuxdisk/tmp/o1_layers/{chm13v2,ggo,chm13v1}.sd.prefix`, each holding one line, the cache prefix `/mnt/linuxdisk/home/juanfraitu/sd_cache/sd/<key>/genome`. Tasks 9-11 read them.

- [ ] **Step 1: The run script** `bench/sedef/run_genome.sh`

```bash
#!/bin/bash
# One whole-genome `sd` stage into the persistent cache, under a 23 GiB address-space cap, timed
# (plan docs/superpowers/plans/2026-10-13-o1-layers-phase-1-2.md, Task 8). An allocation over the cap aborts sd-scan
# instead of waking the kernel's OOM killer.
#   run_genome.sh NAME FASTA   ->  /mnt/linuxdisk/tmp/o1_layers/NAME.sd.{prefix,log,time}
set -euo pipefail
NAME=$1; FA=$2
O=/mnt/linuxdisk/tmp/o1_layers; mkdir -p "$O"
export RUSTLE_SD_CACHE_DIR=/mnt/linuxdisk/home/juanfraitu/sd_cache
ulimit -v 24117248
cd "$(dirname "$0")/../.."
/usr/bin/time -f $'wall_s\t%e\npeak_rss_kb\t%M' -o "$O/$NAME.sd.time" \
  bash tools/rustle_pipeline.sh sd --fasta "$FA" --out "$O/$NAME" --threads 4 --bin /mnt/linuxdisk/home/juanfraitu/rustle_target/release
cat "$O/$NAME.sd.time"
```

- [ ] **Step 2: Run the three genomes, one at a time**

Use Bash with `run_in_background: true`, one command per call, each after the previous one's completion notification. Expect roughly 20-40 minutes each.

```bash
cd /mnt/linuxdisk/home/juanfraitu/rustle_m2_soto && chmod +x bench/sedef/run_genome.sh
RLOCK_TIMEOUT=5400 bash tools/rlock.sh heavy bash bench/sedef/run_genome.sh chm13v2 /mnt/linuxdisk/home/juanfraitu/winloci_data/chm13v2.0.fa
```
```bash
cd /mnt/linuxdisk/home/juanfraitu/rustle_m2_soto
RLOCK_TIMEOUT=5400 bash tools/rlock.sh heavy bash bench/sedef/run_genome.sh ggo /mnt/linuxdisk/home/juanfraitu/_from_wsl/winloci_scratch/GGO.fasta
```
```bash
cd /mnt/linuxdisk/home/juanfraitu/rustle_m2_soto
RLOCK_TIMEOUT=5400 bash tools/rlock.sh heavy bash bench/sedef/run_genome.sh chm13v1 /mnt/linuxdisk/home/juanfraitu/winloci_data/soto_replication/t2t-chm13-v1.0.fa.gz
```

After each run, check it:

```bash
O=/mnt/linuxdisk/tmp/o1_layers; for n in chm13v2 ggo chm13v1; do [ -s $O/$n.sd.prefix ] || continue; P=$(cat $O/$n.sd.prefix); echo "== $n $P"; cat $O/$n.sd.time; grep -v '^genome\b' $P.sd.params.tsv | tr '\n' ' '; echo; done
```

Expected per run: the time file shows `wall_s` and `peak_rss_kb`, and `params.tsv` shows `sequences` 25 / 26 / 23.

**STOP** and report BLOCKED with the last 30 lines of `$O/<name>.sd.log` in any of these cases:
- the driver exits non-zero;
- a `memory allocation ... failed` abort;
- rlock timeout 124.

- [ ] **Step 3: Score V1 and V1b**

```bash
cd /mnt/linuxdisk/home/juanfraitu/rustle_m2_soto
L=/mnt/linuxdisk/tmp/o1_close/lift; W=/mnt/linuxdisk/home/juanfraitu/winloci_data; O=/mnt/linuxdisk/tmp/o1_layers
P2=$(cat $O/chm13v2.sd.prefix); PG=$(cat $O/ggo.sd.prefix)
for id in 0.90 0.98; do
  RLOCK_TIMEOUT=600 bash tools/rlock.sh light python3 bench/sedef/v1_compare.py --sedef $L/sedef_pairs.v2.bed --sedefrs $P2.sd.pairs.tsv \
    --chain $W/soto_replication/v1.0_to_v1.1_rdna_merged.chain --min-id $id --out $O/v1_chm13v2_$id.json
  RLOCK_TIMEOUT=600 bash tools/rlock.sh light python3 bench/sedef/v1_compare.py --sedef $W/GGO_sedef_final.bed --sedefrs $PG.sd.pairs.tsv \
    --fai /mnt/linuxdisk/home/juanfraitu/_from_wsl/winloci_scratch/GGO.fasta.fai --min-id $id --out $O/v1_ggo_$id.json
done
```

Each call prints one JSON line. A cut passes when `recall` >= 0.80 and `precision` >= 0.80 on both genomes. The cut 0.90 is V1; the cut 0.98 is V1b.

- [ ] **Step 4: Write the results into the prereg and the register**

Append to `docs/PREREG_sedef_rs_layers_2026-10-13.md`, under `## Results`, filling the tables from the JSON files and the `.sd.time` / `params.tsv` files:

```markdown
### V1 and V1b (YYYY-MM-DD)

| genome | cut | SEDEF pairs | sd-scan pairs | SEDEF bp | sd-scan bp | recall | precision |
|---|---|---|---|---|---|---|---|
| CHM13 v2.0 | 0.90 | ... | ... | ... | ... | ... | ... |
| gorilla | 0.90 | ... | ... | ... | ... | ... | ... |
| CHM13 v2.0 | 0.98 | ... | ... | ... | ... | ... | ... |
| gorilla | 0.98 | ... | ... | ... | ... | ... | ... |

**V1: PASS | FAIL. V1b: PASS | FAIL.** Reading (fixed in advance): ...

Recall by SEDEF identity / length bin (V1, both genomes): ...

| run | wall s | peak RSS GB | sequences | candidates | blocks | pairs (ok) | segments | families | duplicon rows |
|---|---|---|---|---|---|---|---|---|---|
| CHM13 v2.0 | ... |
| gorilla | ... |
| CHM13 v1.0 | ... |
```

Append two rows to `docs/NEGATIVE_RESULTS_REGISTER.md`, numbered after its last row (1278 and 1279 if still free), in the register's five-column form:

```markdown
| 1278 | YYYY-MM-DD | sedef-rs SD calls vs SEDEF, V1 (`docs/PREREG_sedef_rs_layers_2026-10-13.md`; `utilities sd-scan`, vendored `src/sedef/`) | sd-scan's SDs (>= 1 kb, >= 0.90) and SEDEF's cover each other at >= 0.80, bp-wise, on CHM13 v2.0 and gorilla | <✅ PASS / ⛔ FAIL>: CHM13 recall R / precision P; gorilla R / P; <one sentence on where the misses sit, from the bins> |
| 1279 | YYYY-MM-DD | sedef-rs SD98 vs SEDEF SD98, V1b (same prereg) | the same at >= 0.98 | <✅ PASS / ⛔ FAIL>: CHM13 R / P; gorilla R / P; <reading> |
```

- [ ] **Step 5: Commit**

```bash
git add bench/sedef/run_genome.sh docs/PREREG_sedef_rs_layers_2026-10-13.md docs/NEGATIVE_RESULTS_REGISTER.md
git commit -m "V1/V1b: sd-scan vs SEDEF on CHM13 v2.0 and gorilla (<PASS/FAIL>, <PASS/FAIL>); whole-genome sd caches

<the four recall/precision pairs; wall time and peak RSS per genome>
Plan docs/superpowers/plans/2026-10-13-o1-layers-phase-1-2.md, Task 8.

Co-Authored-By: <your model name> <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_015Xtq7GodhTGPYBzHDp3o4z"
```

---

### Task 9: The hierarchy reads the sd cache (`parse_sd` + `mcl_families --sd-cache`)

**Files:**
- Modify:
  - `src/family.rs`: `family_hierarchy::parse_sd` and its doc comment (around line 78735), plus a test in that module's `mod tests` (around line 79016);
  - `src/mcl_families.rs`: a new `sd_cache` field after `hierarchy_sd` (around line 267); a new fn `hierarchy_inputs`; the `ensure!` at about line 1145; the `write_hierarchy` call at about line 1414; a test in `mod tests`.

**Interfaces:**
- Consumes:
  - Task 4's `PREFIX.sd.pairs.tsv` (11 columns, status in the last) and `PREFIX.sd.duplicons.tsv`;
  - Task 8's `/mnt/linuxdisk/tmp/o1_layers/chm13v2.sd.prefix`.
- Produces:
  - `mcl_families --emit-hierarchy --sd-cache PREFIX`;
  - `fn hierarchy_inputs(args: &Args) -> (Option<String>, Option<String>)`, returning the duplicon track path and the SD pairs path.

- [ ] **Step 1: Write the failing tests**

In `src/family.rs`'s `family_hierarchy` `mod tests`, after `parse_sd_reads_identity_from_column_21`, add:

```rust
        #[test]
        fn parse_sd_reads_sd_scan_pairs_ok_rows_with_percent_identity() {
            let text = "#chr_a\tstart_a\tend_a\tchr_b\tstart_b\tend_b\tstrand\tanchors\test_id\tstage2_id\tstatus\n\
                        chr1\t100\t2100\tchr2\t500\t2500\t-\t1840\t97.6\t97.1\tok\n\
                        chr1\t0\t900\tchr1\t0\t900\t+\t50\t100.0\t96.6\tself\n\
                        chr3\t0\t400\tchr4\t0\t400\t+\t10\t91.0\t90.0\tweak\n";
            let p = parse_sd(text);
            assert_eq!(p.len(), 1);
            assert_eq!(p[0].a, ("chr1".to_string(), 100, 2100));
            assert_eq!(p[0].b, ("chr2".to_string(), 500, 2500));
            assert!((p[0].identity.unwrap() - 0.976).abs() < 1e-12);
        }
```

In `src/mcl_families.rs`'s `mod tests`, after `pair_rule_defaults_to_exact_mcl_and_accepts_exactly_the_two_rules`, add:

```rust
    #[test]
    fn sd_cache_names_both_hierarchy_inputs_and_excludes_the_two_flags() {
        let parse = |v: &[&str]| {
            Args::try_parse_from(["mcl_families", "--out", "o"].into_iter().chain(v.iter().copied()))
        };
        let a = parse(&["--emit-hierarchy", "--sd-cache", "/c/sd/k/genome"]).expect("parse");
        assert_eq!(
            hierarchy_inputs(&a),
            (Some("/c/sd/k/genome.sd.duplicons.tsv".to_string()), Some("/c/sd/k/genome.sd.pairs.tsv".to_string()))
        );
        let b = parse(&["--hierarchy-duplicons", "d.bed", "--hierarchy-sd", "s.bed"]).expect("parse");
        assert_eq!(hierarchy_inputs(&b), (Some("d.bed".to_string()), Some("s.bed".to_string())));
        assert_eq!(hierarchy_inputs(&parse(&[]).expect("parse")), (None, None));
        assert!(parse(&["--sd-cache", "p", "--hierarchy-duplicons", "d.bed"]).is_err());
        assert!(parse(&["--sd-cache", "p", "--hierarchy-sd", "s.bed"]).is_err());
    }
```

- [ ] **Step 2: Run them and watch them fail**

```bash
cd /mnt/linuxdisk/home/juanfraitu/rustle_m2_soto
RLOCK_TIMEOUT=900 bash tools/rlock.sh heavy env CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target TMPDIR=/mnt/linuxdisk/tmp \
  cargo test --release --lib parse_sd > /mnt/linuxdisk/tmp/o1_layers/t9_red_lib.log 2>&1; echo "exit $?"
grep -E "test result|panicked" /mnt/linuxdisk/tmp/o1_layers/t9_red_lib.log
RLOCK_TIMEOUT=900 bash tools/rlock.sh heavy env CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target TMPDIR=/mnt/linuxdisk/tmp \
  cargo test --release --bin mcl_families sd_cache > /mnt/linuxdisk/tmp/o1_layers/t9_red_bin.log 2>&1; echo "exit $?"
grep -E "^error" /mnt/linuxdisk/tmp/o1_layers/t9_red_bin.log | head -3
```

Expected:
- the lib run fails, with `assertion ... left: 3 ... right: 1` from `parse_sd_reads_sd_scan_pairs_...` (all three rows are kept today);
- the bin run fails to compile: `cannot find function `hierarchy_inputs``.

- [ ] **Step 3: Implement**

In `src/family.rs`, replace everything from the line `/// One SEDEF pair (`chrom1 s1 e1 chrom2 s2 e2 ...`, 0-based half-open), identity from column 21 when present.` through the closing brace of `parse_sd` (the `SdPair` struct and the function) with:

```rust
    /// One SD pair, both sides 0-based half-open, from either input format (told apart per row, see `parse_sd`).
    #[derive(Clone, Debug, PartialEq)]
    pub struct SdPair {
        pub a: (String, u64, u64),
        pub b: (String, u64, u64),
        pub identity: Option<f64>,
    }

    /// SD pairs from either format:
    /// - a SEDEF BED (`chrom1 s1 e1 chrom2 s2 e2 ...`): identity from column 21 when present;
    /// - `utilities sd-scan` pairs (`PREFIX.sd.pairs.tsv`: 11 columns, the last `ok` / `self` / `weak`): only `ok` rows
    ///   are kept, identity = column 9 (`est_id`, percent) / 100 (docs/PREREG_sedef_rs_layers_2026-10-13.md).
    pub fn parse_sd(text: &str) -> Vec<SdPair> {
        let mut out = vec![];
        for line in text.lines() {
            if line.starts_with('#') || line.trim().is_empty() {
                continue;
            }
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 6 {
                continue;
            }
            let p = |x: &str| x.parse::<u64>().ok();
            let (Some(s1), Some(e1), Some(s2), Some(e2)) = (p(f[1]), p(f[2]), p(f[4]), p(f[5])) else { continue };
            let sd_scan = f.len() == 11 && matches!(f[10], "ok" | "self" | "weak");
            if sd_scan && f[10] != "ok" {
                continue;
            }
            let identity = if sd_scan {
                f[8].parse::<f64>().ok().map(|x| x / 100.0)
            } else {
                f.get(20).and_then(|x| x.parse::<f64>().ok())
            };
            out.push(SdPair { a: (f[0].to_string(), s1, e1), b: (f[3].to_string(), s2, e2), identity });
        }
        out
    }
```

(`SdPair`'s fields are unchanged; only its doc comment and `parse_sd` change.)

In `src/mcl_families.rs`, after the `hierarchy_sd: Option<String>,` field, add:

```rust
    /// `utilities sd-scan` output prefix for `--emit-hierarchy` (the driver's `sd` stage writes it to `PREFIX.sd.prefix`):
    /// duplicons from `<PREFIX>.sd.duplicons.tsv` (sedef-rs duplication families) and the SD layer from
    /// `<PREFIX>.sd.pairs.tsv` (its `ok` rows). It takes the place of `--hierarchy-duplicons` and `--hierarchy-sd`, which it
    /// excludes (docs/superpowers/specs/2026-10-13-o1-layered-hierarchy-design.md).
    #[arg(long, conflicts_with_all = ["hierarchy_duplicons", "hierarchy_sd"])]
    sd_cache: Option<String>,
```

Add, next to `write_hierarchy` (above its doc comment):

```rust
/// `--emit-hierarchy`'s duplicon track and SD pairs: `--sd-cache PREFIX` names both, else the two flags.
fn hierarchy_inputs(args: &Args) -> (Option<String>, Option<String>) {
    match &args.sd_cache {
        Some(p) => (Some(format!("{p}.sd.duplicons.tsv")), Some(format!("{p}.sd.pairs.tsv"))),
        None => (args.hierarchy_duplicons.clone(), args.hierarchy_sd.clone()),
    }
}
```

Replace the `ensure!` that reads `"--emit-hierarchy needs --hierarchy-duplicons (the DupMasker colors BED)"` with:

```rust
    anyhow::ensure!(
        !args.emit_hierarchy || hierarchy_inputs(&args).0.is_some(),
        "--emit-hierarchy needs --hierarchy-duplicons (the DupMasker colors BED) or --sd-cache (a `utilities sd-scan` prefix)"
    );
```

In the `if args.emit_hierarchy { ... write_hierarchy(...) }` block, replace the two arguments `args.hierarchy_duplicons.as_deref().unwrap_or_default(),` and `args.hierarchy_sd.as_deref(),` with `dup.as_deref().unwrap_or_default(),` and `sd.as_deref(),`. Declare them first, as the block's first line:

```rust
        let (dup, sd) = hierarchy_inputs(&args);
```

Finally, in the `emit_hierarchy` field's doc comment, change `Needs `--gff` and `--hierarchy-duplicons`.` to `Needs `--gff` and `--hierarchy-duplicons` (or `--sd-cache`).`

- [ ] **Step 4: Run the tests and watch them pass; the full suite**

```bash
cd /mnt/linuxdisk/home/juanfraitu/rustle_m2_soto
RLOCK_TIMEOUT=900 bash tools/rlock.sh heavy env CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target TMPDIR=/mnt/linuxdisk/tmp \
  cargo test --release > /mnt/linuxdisk/tmp/o1_layers/t9_full.log 2>&1; echo "exit $?"
grep -E "test result" /mnt/linuxdisk/tmp/o1_layers/t9_full.log | awk '{p+=$4; f+=$6} END{print "passed", p, "failed", f}'
grep -E "^test .* FAILED|^---- " /mnt/linuxdisk/tmp/o1_layers/t9_full.log | head
```

Expected:
- every test passes except the 3 known `pipeline_driver_dispatch` failures (exit 127, present before this plan);
- the new `parse_sd_reads_sd_scan_...` and `sd_cache_names_...` pass, and so does `parse_sd_reads_identity_from_column_21`;
- any other failure is this task's to fix.

Run this step with `run_in_background: true` if it exceeds the 10-minute Bash limit.

- [ ] **Step 5: The real run on A119b (descriptive: does the cache drive the hierarchy end to end?)**

```bash
cd /mnt/linuxdisk/home/juanfraitu/rustle_m2_soto
RLOCK_TIMEOUT=900 bash tools/rlock.sh heavy env CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target TMPDIR=/mnt/linuxdisk/tmp \
  cargo build --release --bin mcl_families > /mnt/linuxdisk/tmp/o1_layers/t9_build.log 2>&1; echo "build exit $?"
P2=$(cat /mnt/linuxdisk/tmp/o1_layers/chm13v2.sd.prefix); S=human_A119b; R=/mnt/linuxdisk/tmp/rustle_figures/runs/$S
O=/mnt/linuxdisk/tmp/o1_layers/hier_sdcache/$S; mkdir -p $O
RLOCK_TIMEOUT=590 bash tools/rlock.sh heavy /mnt/linuxdisk/home/juanfraitu/rustle_target/release/mcl_families --paf $R/$S.fam.loci.paf --gff $R/$S.fam.loci.gff3 \
  --min-exonic-bp 1 --min-shared-exon-frac 0.60 --emit-hierarchy --sd-cache $P2 --out $O/ours > $O/log 2>&1; echo "exit $?"
grep "emit-hierarchy" $O/log /mnt/linuxdisk/tmp/o1_close/hier/$S/log
cmp $O/ours.clusters.tsv /mnt/linuxdisk/tmp/o1_close/hier/$S/ours.clusters.tsv && echo "families unchanged (annotation only)"
```

Expected:
- `exit 0`;
- two `--emit-hierarchy:` summary lines (sd-scan vs plan 1's lifted DupMasker + SEDEF): members, superfamilies and families equal, while expansions, duplicons owned, SD blocks and SD complexes may differ;
- `families unchanged (annotation only)`.

Record both lines in the commit message.

- [ ] **Step 6: Commit**

```bash
git add src/family.rs src/mcl_families.rs
git commit -m "hierarchy: --sd-cache PREFIX (sd-scan duplicons + ok SD pairs); parse_sd reads sd-scan pairs

<the two emit-hierarchy lines (sd-scan vs DupMasker/SEDEF) on human_A119b; families unchanged>
Plan docs/superpowers/plans/2026-10-13-o1-layers-phase-1-2.md, Task 9.

Co-Authored-By: <your model name> <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_015Xtq7GodhTGPYBzHDp3o4z"
```

---

### Task 10: V2 DEV (duplicons, F_O universe)

**Files:**
- Modify: `docs/PREREG_sedef_rs_layers_2026-10-13.md` (`## Results`: a V2 DEV subsection)

**Interfaces:**
- Consumes:
  - Task 8's `/mnt/linuxdisk/tmp/o1_layers/chm13v1.sd.prefix`;
  - `bench/paf_duplicons/score_units.py` (unchanged);
  - the `mcl_families` binary built in Task 9.
- Produces: `/mnt/linuxdisk/tmp/o1_layers/v2/{dupmasker,sdscan}/ours.*` and `score_dev.json`, which Task 11 reuses.

- [ ] **Step 1: Build both arms, check H0, compare the DupMasker arm with the registered run**

```bash
cd /mnt/linuxdisk/home/juanfraitu/rustle_m2_soto
W=/mnt/linuxdisk/home/juanfraitu/winloci_data; H=/mnt/linuxdisk/tmp/hier_flag; V=/mnt/linuxdisk/tmp/o1_layers/v2
B=/mnt/linuxdisk/home/juanfraitu/rustle_target/release; P1=$(cat /mnt/linuxdisk/tmp/o1_layers/chm13v1.sd.prefix)
DM=$W/duplicons/chm13.draft_v1.0_plus38Y_dupmasker_colors.bed
for arm in dupmasker sdscan; do
  if [ $arm = dupmasker ]; then DUP=$DM; else DUP=$P1.sd.duplicons.tsv; fi
  mkdir -p $V/$arm
  RLOCK_TIMEOUT=590 bash tools/rlock.sh heavy $B/mcl_families --paf /mnt/linuxdisk/tmp/ourhom/loci.paf --gff /mnt/linuxdisk/tmp/ourhom/loci.gff3 \
    --min-exonic-bp 1 --min-shared-exon-frac 0.60 --dump-pairs --emit-hierarchy --hierarchy-duplicons $DUP \
    --hierarchy-sd $W/soto_replication/final_v1_clean.bed --hierarchy-external $H/soto_external.tsv --out $V/$arm/ours > $V/$arm/log 2>&1
  echo "$arm exit $?"; grep "emit-hierarchy" $V/$arm/log
done
for f in clusters.tsv loci.tsv pairs.tsv; do cmp $V/dupmasker/ours.$f $V/sdscan/ours.$f && echo "H0 $f identical"; done
for f in clusters.tsv loci.tsv pairs.tsv hierarchy.members.tsv hierarchy.nodes.tsv hierarchy.edges.tsv hierarchy.external.tsv; do
  cmp -s $V/dupmasker/ours.$f $H/final/ours.$f && echo "dupmasker arm $f = the 10-09 run" || echo "dupmasker arm $f DIFFERS from the 10-09 run"
done
```

Expected:
- both arms `exit 0`;
- three `H0 ... identical` lines. If one is missing, STOP: H0 failed, so the arms are not comparable;
- the 10-09 comparison is descriptive: report each line. If a file differs, the re-run arm is the comparator (both arms share one binary). Give the reason in the results, e.g. `diff <(cut -f1-3 ...) | head`.

- [ ] **Step 2: Score DEV only**

```bash
cd /mnt/linuxdisk/home/juanfraitu/rustle_m2_soto
W=/mnt/linuxdisk/home/juanfraitu/winloci_data; V=/mnt/linuxdisk/tmp/o1_layers/v2; P1=$(cat /mnt/linuxdisk/tmp/o1_layers/chm13v1.sd.prefix)
RLOCK_TIMEOUT=600 bash tools/rlock.sh light python3 bench/paf_duplicons/score_units.py $V/dupmasker $V/sdscan dev /mnt/linuxdisk/tmp/ourhom/loci.gff3 \
  $W/duplicons/chm13.draft_v1.0_plus38Y_dupmasker_colors.bed $P1.sd.duplicons.tsv > $V/score_dev.json; echo "exit $?"
python3 -c 'import json,sys; d=json.load(open(sys.argv[1])); print({k: d[k] for k in ("families","P1_mean_family_score","member_weighted_score","families_identical","exact_sensitivity","exact_precision","P2_wins","P2_losses","P2_p")})' $V/score_dev.json
```

Do NOT run `heldout`, and do NOT run `hier_match.py` on these arms in this task.

- [ ] **Step 3: Write the DEV result into the prereg**

Append under `## Results`:

```markdown
### V2 DEV (YYYY-MM-DD; HELD-OUT not computed)

H0: clusters / loci / pairs identical between the arms (yes | no). DupMasker arm vs the registered 10-09 run: <per file>.

| measure (DEV, N families, S structured) | value |
|---|---|
| mean family score (P1 on DEV) | ... |
| member-weighted score | ... |
| families identical | ... / N |
| exact expansion sensitivity / precision | ... / ... |
| P2: wins / losses / ties vs the one-block null, sign test p | ... |
| P3: Soto families exactly one expansion, DupMasker / sd-scan (all families) | ... / ... |
| P4: PKD1P6-NPIPP1, PDXDC2P-NPIPB14P (sd-scan arm) | ... |
| members' exonic bases carrying a duplicon, DupMasker / sd-scan | ... / ... |
| NPIP (MCL9) cover, DupMasker / sd-scan | ... / ... |
| blocks per family, DupMasker / sd-scan | ... / ... |

Worst families (the 12 lowest scores): ...
```

- [ ] **Step 4: Commit**

```bash
git add docs/PREREG_sedef_rs_layers_2026-10-13.md
git commit -m "V2 DEV: sd-scan duplicons vs DupMasker on F_O (mean family score <x>; held-out not computed)

Plan docs/superpowers/plans/2026-10-13-o1-layers-phase-1-2.md, Task 10.

Co-Authored-By: <your model name> <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_015Xtq7GodhTGPYBzHDp3o4z"
```

**Controller step between Tasks 10 and 11 (not an implementer's).**
- Read the DEV result. Decide whether to register the prereg's one allowed amendment, e.g. a different duplicon construction from the same sd-scan segments.
- An amendment is written into the prereg as `### Amendment 1 (registered before HELD-OUT)` and committed. Its implementation is a task of its own, with tests; if it changes an output, also bump `SCAN_VERSION`.
- With no amendment, Task 11 runs on the arms as they are. In either case the decision is recorded in the ledger.

---

### Task 11: V2 HELD-OUT and V2b; verdicts; docs

**Files:**
- Modify:
  - `docs/PREREG_sedef_rs_layers_2026-10-13.md` (`## Results`: HELD-OUT, V2b, verdict);
  - `docs/NEGATIVE_RESULTS_REGISTER.md` (one row);
  - `docs/PENDING_2026-10-11.md` (a new section at the end);
  - `docs/superpowers/specs/2026-10-13-o1-layered-hierarchy-design.md` (a results pointer at the end);
  - `docs/MODULE_STATUS.md` (the note's last sentence);
  - `src/sedef/mod.rs` (the STATUS lines' last sentence).

**Interfaces:**
- Consumes: Task 10's arms and DEV result. If the controller registered an amendment, use the arms its implementation task rebuilt.
- Produces: the V2 verdict that phase 3's plan reads.

- [ ] **Step 1: HELD-OUT and V2b**

```bash
cd /mnt/linuxdisk/home/juanfraitu/rustle_m2_soto
git log --oneline -3 -- docs/PREREG_sedef_rs_layers_2026-10-13.md
W=/mnt/linuxdisk/home/juanfraitu/winloci_data; H=/mnt/linuxdisk/tmp/hier_flag; V=/mnt/linuxdisk/tmp/o1_layers/v2; P1=$(cat /mnt/linuxdisk/tmp/o1_layers/chm13v1.sd.prefix)
RLOCK_TIMEOUT=600 bash tools/rlock.sh light python3 bench/paf_duplicons/score_units.py $V/dupmasker $V/sdscan heldout /mnt/linuxdisk/tmp/ourhom/loci.gff3 \
  $W/duplicons/chm13.draft_v1.0_plus38Y_dupmasker_colors.bed $P1.sd.duplicons.tsv > $V/score_heldout.json; echo "exit $?"
python3 -c 'import json,sys; d=json.load(open(sys.argv[1])); print({k: d[k] for k in ("families","P1_mean_family_score","member_weighted_score","families_identical","exact_sensitivity","exact_precision","P2_wins","P2_losses","P2_p")})' $V/score_heldout.json
for arm in dupmasker sdscan; do
  RLOCK_TIMEOUT=600 bash tools/rlock.sh light python3 bench/o1_hier/hier_match.py $V/$arm/ours.hierarchy.members.tsv $H/soto_external.tsv $V/$arm/key.json \
    --fold $V/$arm/ours.loci.tsv --map key --null 0 > /dev/null
  python3 -c 'import json,sys; print(sys.argv[1], json.load(open(sys.argv[2]))["exact_any"])' $arm $V/$arm/key.json
done
```

The `git log` line must show the DEV commit, and the amendment's if there was one, before this run.

**Verdict:**
- **V2a** holds when held-out `P1_mean_family_score` >= 0.80.
- **V2b** holds when the sd-scan arm's `expansion_not_family` >= ⌈0.9 × the DupMasker arm's⌉.
- **V2 PASS** = both hold.

- [ ] **Step 2: Write the verdicts**

Append to the prereg's `## Results`:

```markdown
### V2 HELD-OUT and V2b (YYYY-MM-DD)

<the same table as DEV, for HELD-OUT; DEV beside it>

V2b: Soto families equal to an expansion strictly below a family (`exact_any.expansion_not_family`, key mapping):
DupMasker <n>, sd-scan <m>; bar ⌈0.9 × n⌉ = <b>.

**V2: PASS | FAIL (V2a <PASS/FAIL>, V2b <PASS/FAIL>).** Reading (fixed in advance): <the matching reading, verbatim>.
```

Append one register row (after the last, 1280 if free):

```markdown
| 1280 | YYYY-MM-DD | sedef-rs duplicons as the expansion layer's source vs DupMasker, V2 (`docs/PREREG_sedef_rs_layers_2026-10-13.md`; F_O, CHM13 v1.0; scorer of row 1263) | expansions built from sd-scan duplicons agree with DupMasker's (held-out mean family score >= 0.80) and keep Soto's sub-family expansions (>= 90% of DupMasker's) | <✅ ADOPTED / ⛔ NOT ADOPTED>: held-out P1 <x> (DEV <y>), V2b <m> vs <n>; <one sentence on where it disagrees> |
```

Append to `docs/PENDING_2026-10-11.md`:

```markdown
## O1 layers, phases 1-2: sedef-rs SD and duplicon layers (YYYY-MM-DD)

- **sd-scan.** `utilities sd-scan` (vendored sedef-rs, `src/sedef/`) and the driver's `sd` stage. Caches:
  `/mnt/linuxdisk/home/juanfraitu/sd_cache` (CHM13 v2.0, gorilla, CHM13 v1.0; prefixes in
  `/mnt/linuxdisk/tmp/o1_layers/*.sd.prefix`). `mcl_families --emit-hierarchy --sd-cache PREFIX` reads it.
- **Verdicts:** V1 <PASS/FAIL> (rows 1278); V1b <PASS/FAIL> (1279); V2 <PASS/FAIL> (1280). Prereg:
  `docs/PREREG_sedef_rs_layers_2026-10-13.md`.
- **Next:** phase 3 of `docs/superpowers/specs/2026-10-13-o1-layered-hierarchy-design.md` (the Rust Soto core, S1; assembly
  copy number, CN1; S2). Its SD98 source follows V1b's reading. Its plan is written from these results.
```

Append to the spec:

```markdown

## Results of phases 1-2 (YYYY-MM-DD)

V1 <PASS/FAIL>, V1b <PASS/FAIL>, V2 <PASS/FAIL>: `docs/PREREG_sedef_rs_layers_2026-10-13.md` (register 1278-1280). The
readings fixed there decide what phases 3-4 build on.
```

In `docs/MODULE_STATUS.md`'s note and in `src/sedef/mod.rs`'s STATUS lines, replace "Adoption as the SD / duplicon source is decided by" (mod.rs) and "Its adoption as the SD and duplicon source is decided by" (MODULE_STATUS.md) with the outcome, for example: `Adopted as the SD source (V1 PASS); duplicons not adopted (V2 FAIL): see`. Keep the prereg path after it.

- [ ] **Step 3: Rebuild and run the sedef tests once more (mod.rs changed)**

```bash
cd /mnt/linuxdisk/home/juanfraitu/rustle_m2_soto
RLOCK_TIMEOUT=900 bash tools/rlock.sh heavy env CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target TMPDIR=/mnt/linuxdisk/tmp \
  cargo test --release --lib -- sedef:: module_status > /mnt/linuxdisk/tmp/o1_layers/t11_test.log 2>&1; echo "exit $?"
grep "test result" /mnt/linuxdisk/tmp/o1_layers/t11_test.log
```

Expected: `exit 0`; all pass (libtest takes several filters after `--`).

- [ ] **Step 4: Commit**

```bash
git add docs/PREREG_sedef_rs_layers_2026-10-13.md docs/NEGATIVE_RESULTS_REGISTER.md docs/PENDING_2026-10-11.md \
  docs/superpowers/specs/2026-10-13-o1-layered-hierarchy-design.md docs/MODULE_STATUS.md src/sedef/mod.rs
git commit -m "V2: sd-scan duplicons vs DupMasker on F_O, held-out (<PASS/FAIL>); phase 1-2 verdicts recorded

<held-out P1, V2b counts; V1/V1b/V2 one line each>
Plan docs/superpowers/plans/2026-10-13-o1-layers-phase-1-2.md, Task 11.

Co-Authored-By: <your model name> <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_015Xtq7GodhTGPYBzHDp3o4z"
```
