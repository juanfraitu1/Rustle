# Closing O1, phases 0 + A: MCL two-node fix and the hierarchy-aware evaluation

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or
> superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:**
1. Fix MCL's floating-point decision on two-node components.
2. Re-cluster the four default catalogs.
3. Score human families against Soto with a hierarchy-aware Hungarian matcher: dev A119b, held-out testis.

**Architecture:**
- **Rust (one private helper in `src/family.rs` `annotation_families`):** a connected two-node component is one
  cluster after MCL.
- **Re-clustering:** `mcl_families --paf/--gff` on each catalog's cached all-vs-all.
- **Python bench (`bench/o1_hier/`):**
  - lifts the v1.0 DupMasker / SEDEF / CAT inputs to CHM13 v2.0 through the T2T v1.0 -> v1.1 chain;
  - builds the Soto truth with exons;
  - scores hierarchy nodes against it with scipy's Hungarian solver and two label-shuffle nulls.
- Phases V (`sedef-rs`) and B (apes) get their own plan after A's results.

**Tech stack:** Rust (rustle crate), Python 3 (numpy 2.4, scipy 1.17), samtools.

**Spec:** `docs/superpowers/specs/2026-10-13-close-o1-design.md`

## Global Constraints

- **Build and run rules:**
  - Build only with `CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target TMPDIR=/mnt/linuxdisk/tmp`, always
    `--release`, cargo output to a log file.
  - Every cargo / `mcl_families` / famsim run goes through `bash tools/rlock.sh heavy` (`RLOCK_TIMEOUT=590`,
    600000 ms tool timeout). Python scoring goes through `bash tools/rlock.sh light` unless a step says heavy.
  - Foreground only; never `pkill`. Another Claude session (PID 2532845, cwd the Desktop) was alive when this plan was
    written; the locks keep heavy jobs apart, so never bypass them.
- **Commits:**
  - path-limited only: never `git add -A`; another session's untracked `bench/hierarchy/*` and
    `bench/ideal_expression/*` must not be added;
  - every commit ends with
    `Co-Authored-By: <your model> <noreply@anthropic.com>` and
    `Claude-Session: https://claude.ai/code/session_015Xtq7GodhTGPYBzHDp3o4z`.
- **Coordinates:**
  - The default catalogs are CHM13 **v2.0**. DupMasker colors, SEDEF calls and CAT v4 are **v1.0**, and are lifted by
    the T2T chain `winloci_data/soto_replication/v1.0_to_v1.1_rdna_merged.chain` (chr1-22, X, M of v2.0 = v1.1).
  - Truth genes join by Gene ID (S1C `Gene ID` = CAT v4 column 19), never by name.
- **Pre-registration:** every scored run is pre-registered (`docs/PREREG_*.md` committed) before it runs. Nothing above
  a Results section is edited after a run.
- **Fix 0 scope:** it changes ONLY connected two-node components. Every other cluster must stay identical as a member
  set with identical values. Cluster IDs are ranks, so they renumber: compare member sets, never IDs.
- **Bars** (spec, plus the rulings below):
  - **E1:** beat the global null at p < 0.001.
  - **E1b:** beat the within-family null at p < 0.001.
  - **E2:** |held-out exact rate - dev exact rate| <= 0.10, rate = exact / evaluable.
- **Outputs:**
  - all run outputs go to `/mnt/linuxdisk/tmp/o1_close/` (~55 GB free);
  - the shipped catalogs in `/mnt/linuxdisk/tmp/rustle_figures/runs/` are read, never written.

## Rulings made while writing this plan (deviations from the spec text; review these)

1. **SEDEF input and lift.**
   - The 10-09 hierarchy used `winloci_data/soto_replication/final_v1_clean.bed` (151,758 SEDEF pairs, identity in
     column 21), not `HSA_sedef_pairs.bed`; its column 21 is not an identity.
   - Nothing was lifted on 10-09: the F_O universe was v1.0.
   - Here the inputs are lifted by the T2T chain. Liftoff is not used: it can swap near-identical paralogs.
   - Probe: DupMasker 199,154 rows lifted, SEDEF lifted, 1,000 / 1,000 sampled intervals sequence-identical for each.
   - Dropped rows: GRCh38 chrY 6,306 (no chain); acrocentric rDNA-region rows 355.
2. **Truth mapping is exon + strand, not span overlap.** A Soto gene joins a member when an exon shares >= 1 bp on the
   same strand with an exon of the member's locus (or of a locus folded into it). On the 10-09 F_O data, raw span overlap
   cut the exact-expansion count from 215 to 169 by joining antisense and intronic neighbours; exon + strand gives 206.
3. **Evaluable truth family = >= 2 members on our loci** (the 10-09 convention: "the 400 Soto families with >= 2
   members").
   - E2's rate is exact / evaluable, so a tissue's expression coverage does not decide it.
   - Unconditioned numbers (exact / 491, gene recall) are reported next to it.
4. **E1b, a within-family null, is added.**
   - The spec's global null is passed trivially: on F_O it averages 0.03 exact matches.
   - The within-family null keeps families, superfamilies and SD complexes, and moves only the expansions. It is the
     test that controls depth.
   - On F_O it came within 7 exact matches of the observed value.
5. **The pair rule runs after MCL, in `annotation_families::mcl` only.**
   - This is equivalent to "before MCL": flow never crosses components.
   - The separate CSC MCL in `src/utilities.rs:633` (a `utilities` subcommand, not the families stage) is untouched.
6. **Re-clustering drops `--fasta` / `--emit-units`.** Neither changes `clusters.tsv` or `loci.tsv`.
   - Probe: the HEAD binary through `--paf --gff` reproduced `human_testis.fam.clusters.tsv` and `.fam.loci.tsv`
     byte-identically in 5.7 s.
   - The shipped copy tables are not rebuilt (out of scope).

## Facts established while writing this plan (do not re-derive)

- **MCL defaults** (`src/mcl_families.rs`): inflation 2.8, prune **1e-9**, max_iter 100.
  - The two-node probe (Python port, same operation order) at 2.8 / 1e-9 splits 160 of the 500 weights 0.500..0.999,
    including 0.99, 0.985, 0.90 and 0.75.
  - With a triangle and an isolated node in the same graph, 0.99 JOINS and 0.90 splits. Convergence is tested on all
    columns at once, so a pair's fate also depends on unrelated components.
- **`HomologyGraph.edges`** keys are `(i, j)` with `i < j`. `attractors` returns clusters sorted by (size desc, first
  member).
- **The hierarchy flag:**
  - It writes `<out>.hierarchy.members.tsv` with columns `sd_complex superfamily family member strand ... expansions
    ... sd_blocks sd_components sd_max_identity`. `expansions` is a comma list that includes `<family>.unique`
    pseudo-expansions.
  - Node levels: `sd_complex superfamily family expansion` (scored), plus `duplicon sd_block external` (not scored:
    duplicons have no members, and SD blocks are not in the spec's levels).
  - Member keys are GFF 1-based `chrom:start-end`. The fold map `<out>.loci.tsv` maps record keys to representatives,
    and also lists representatives that are not members (single-locus clusters).
- **Dry run on v2.0** (testis, HEAD binary, lifted inputs, no truth):
  - 1,115 members, 342 families, 398 expansion nodes, 335 superfamilies, 149 SD complexes (largest 547 members);
  - families byte-identical with the flag;
  - every member key is in `.fam.loci.gff3`.
- **Soto S1C** (`bench/soto/soto_famCN_S1C.tsv`): 2,572 rows.
  - 2,458 rows in 491 families proper; 114 `Unassigned_*`.
  - All 2,458 Gene IDs are found in CAT v4 (`winloci_data/soto_replication/cat_v4.bed`, BED12, v1.0, 232,792
    transcripts).
  - 1,961 rows have spans identical to `hier_flag/soto_external.tsv`. The other 497 differ because the hull over all
    transcripts is wider.
- **Scorer regression on the 10-09 F_O run (v1.0, guided universe):**
  - `--map key` with `hier_flag/soto_external.tsv`: 491 labels, 400 evaluable, 12 one-member, 79 none; `exact_any`
    expansion **215** (= the Rust placement); hierarchy exact 245 vs flat 228.
  - `--map exon` with the v1.0 truth table: 402 evaluable; `exact_any` expansion **206**; hierarchy exact 231 vs
    flat 217.
- **Run times:** testis re-clustering 5.7 s / 350 MB; hierarchy 2.6 s / 440 MB. A119b's PAF is 11x larger; budget
  ~2 min, heavy lock.
- **Probe outputs** are in `/mnt/linuxdisk/tmp/o1_close/probe/` (scratch; ignore).

## Review Focus

1. **A Soto gene that joins no member.** This covers a gene that is intronic, antisense, on an unclustered locus or
   absent. It must be counted in coverage and against the unconditioned rates, never dropped silently (Task 4:
   `test_coverage_separates_...`, `test_a_shared_exon_on_the_same_strand_joins_...`).
2. **A family and an expansion with the identical member set.** The pair must be matched once, preferring the family
   (the named tie policy) (Task 4: `test_identical_sets_at_two_levels_...`, `test_one_to_one_...`).
   As built, the tie policy is canonical placement (8f0515ef; spec amendment 5): see the prereg's Tie policy.
3. **An interval with an end in a chain gap, or on a chromosome without a chain.**
   - Such an interval is unliftable and counted.
   - An interval spanning an indel lifts with it and is flagged, never shifted silently (Task 3:
     `test_an_end_inside_a_gap_is_unliftable`, `test_an_interval_spanning_the_indel_...`).
4. **The pair rule.** It must not touch a triangle, a longer path or an edge-free node, and must keep the shipped
   order (Task 1: `the_pair_rule_touches_only_two_node_components`, `a_pair_joins_inside_a_mixed_graph_...`).
5. **The nulls.**
   - They keep node sizes and nesting.
   - `within_family` must leave every family, superfamily and SD complex unchanged.
   - Self-match is perfect (Task 4: `test_permutations_keep_sizes_...`, `test_self_match_is_perfect_...`).

---

### Task 1: pre-register fix 0, then the two-node rule in `mcl` (TDD), build, famsim

**Files:**
- Create: `docs/archive/2026-10/merged/PREREG_mcl_two_node_2026-10-13.md`
- Modify: `src/family.rs`, inside `pub mod annotation_families` (opens near line 59803):
  - `pub fn mcl(` (~60875): its last line `attractors(&col, n)`;
  - add a private helper right after `fn attractors` (which ends ~60982);
  - add tests to the module's `mod tests { use super::*; ... }` (opens ~61106), right after the test
    `mcl_splits_two_cliques_joined_by_a_single_bridge` (~61318).

**Interfaces:**
- Produces:
  - `/mnt/linuxdisk/tmp/o1_close/mcl_families.before` (HEAD, pre-fix) and `.after` (fixed). Task 2 uses both,
    Task 6 uses `.after`.
  - The signature `mcl(g: &HomologyGraph, inflation: f64, max_iter: usize, prune: f64) -> Vec<Vec<usize>>` and its output
    order (size desc, then first member) are unchanged.

- [ ] **Step 1: Write and commit the prereg**

Create `docs/archive/2026-10/merged/PREREG_mcl_two_node_2026-10-13.md` with exactly:

```markdown
# PREREG — MCL two-node components: the fix (2026-10-13)

**Bug** (`project_mcl_two_node_fragility`, found 10-02 by the famsim ladder; `bench/FAMSIM.md`, "What the two bold rows
are", item 1).
- `annotation_families::mcl` adds self-loops of 1.0, while an edge weighs identity x cov_longer < 1.
- On a connected two-node component the column converges toward a tie between "self" and "mate". The convergence test
  (`|u - v| < 1e-7`, applied to all columns at once) stops at a residual whose SIGN picks the attractor.
- At the shipped inflation 2.8 / prune 1e-9, 160 of the 500 weights 0.500..0.999 split a lone pair, among them 0.99,
  0.985, 0.90 and 0.75.
- With other components in the graph the set changes: next to a triangle, 0.99 joins and 0.90 splits. So a pair's
  existence also depended on unrelated families.
- Three or more nodes join. 57% of genome-wide families have two members.

**Fix.**
- After the attractors are read, a connected component of exactly two nodes (an edge whose two endpoints have no other
  edge) is one cluster.
- Every other cluster passes through unchanged, in the same order (size desc, first member).
- Flow never crosses components, so this equals deciding pairs before MCL.
- Scope: `annotation_families::mcl`, i.e. the `mcl_families` families stage. The separate CSC MCL in
  `src/utilities.rs` is not changed.

**Bars.**
- **F1** (unit, red first):
  - every probe weight, and all 500 weights 0.500..0.999, give one cluster `[0, 1]`;
  - a graph with a triangle, a pair at 0.90 and an edge-free node gives `[[0,1,2],[3,4],[5]]`;
  - the helper leaves a path, a triangle and an edge-free node alone.
- **F2** (the four default catalogs, re-clustered from their cached `.fam.loci.paf` / `.fam.loci.gff3`):
  - H0: the pre-fix binary reproduces each shipped `.fam.clusters.tsv` and `.fam.loci.tsv` byte-identically. A failure
    stops the task.
  - `bench/o1_hier/pair_diff.py` (member sets; IDs are ranks and are ignored) must give `pairs_split` 0 and
    `other_changes` 0.
  - Reported per catalog: pairs joined, clusters and members before and after.
- **F3** (famsim): after the fix, `snp_0.01` and `snp_0.1` (human SNRPB and orangutan SPAG7 templates) and human
  `shuffle` have `denovo_intact` and `guided_intact` = 1. The before state is recorded.
- **F4:** `cargo test --release` fails only the 3 known `tests/pipeline_driver_dispatch.rs` tests.

**Reading.** If F1-F4 pass, the fix is the default.
- It changes the output of the default families stage.
- The shipped catalogs are not rewritten here. The next families-stage run, and the paused figure rescoring, pick it up.

## Results

(pending)
```

```bash
git add docs/archive/2026-10/merged/PREREG_mcl_two_node_2026-10-13.md
git commit -m "PREREG MCL two-node fix (F1-F4)"
```

- [ ] **Step 2: Keep the pre-fix binary and record the famsim before-state** (before ANY cargo command: a full
  `cargo test` rebuilds the binaries)

```bash
mkdir -p /mnt/linuxdisk/tmp/o1_close/famsim
T=/mnt/linuxdisk/home/juanfraitu/rustle_target/release
echo "binary $(stat -c %Y $T/mcl_families), last src commit $(git log -1 --format=%ct -- src/ Cargo.toml)"
```

- The binary's mtime must be newer than the last src commit. When this plan was written: binary 12:59, commit
  12:42 on 10-13.
- If it is older, first run
  `RLOCK_TIMEOUT=590 CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target TMPDIR=/mnt/linuxdisk/tmp bash tools/rlock.sh heavy bash -c "cargo build --release > /mnt/linuxdisk/tmp/o1_close/build_head.log 2>&1"`.

Then:

```bash
cp /mnt/linuxdisk/home/juanfraitu/rustle_target/release/mcl_families /mnt/linuxdisk/tmp/o1_close/mcl_families.before
F=/mnt/linuxdisk/tmp/o1_close/famsim
RLOCK_TIMEOUT=590 bash tools/rlock.sh heavy python3 bench/famsim ladder --template-spec bench/famsim/examples/ladder_human_SNRPB.json --out $F/before_human --only snp_0.01,snp_0.1,shuffle > $F/before_human.log 2>&1
RLOCK_TIMEOUT=590 bash tools/rlock.sh heavy python3 bench/famsim ladder --template-spec bench/famsim/examples/ladder_ppy_SPAG7.json --out $F/before_ppy --only snp_0.01,snp_0.1 > $F/before_ppy.log 2>&1
cut -f1,4-7 $F/before_human/ladder.tsv $F/before_ppy/ladder.tsv
```

- Expected (the 10-02 tables, `bench/famsim/examples/ladder_*_2026-10-02.tsv`): `snp_0.01` and `snp_0.1` show
  `denovo_intact 0 guided_intact 0` for both templates, and human `shuffle` shows `denovo_intact 0`.
- Record whatever it shows. A rung that is already intact before the fix is noted, not "fixed".

- [ ] **Step 3: Write the failing tests** (in `annotation_families`' test module, right after
  `mcl_splits_two_cliques_joined_by_a_single_bridge`)

```rust
            fn pair_graph(w: f64) -> HomologyGraph {
                let mut g = HomologyGraph::default();
                g.genes.push(("c".to_string(), 1, 100));
                g.genes.push(("c".to_string(), 201, 300));
                g.edges.insert((0, 1), w);
                g
            }

            /// PREREG_mcl_two_node_2026-10-13 F1: a connected pair is ONE family at every weight. Before the fix the
            /// sign of a 1e-7 convergence residual decided it: the probe split 0.99, 0.985, 0.90 and 0.75, and 160 of
            /// the 500 weights 0.500..0.999 split at the shipped inflation 2.8 / prune 1e-9.
            #[test]
            fn a_connected_pair_is_one_family_at_every_weight() {
                for w in [0.999, 0.995, 0.99, 0.985, 0.98, 0.97, 0.95, 0.93, 0.90, 0.85, 0.80, 0.75] {
                    assert_eq!(mcl(&pair_graph(w), 2.8, 100, 1e-9), vec![vec![0, 1]], "w = {w}");
                }
                let split: Vec<f64> = (0..500)
                    .map(|i| 0.5 + 0.001 * i as f64)
                    .filter(|&w| mcl(&pair_graph(w), 2.8, 100, 1e-9).len() != 1)
                    .collect();
                assert!(split.is_empty(), "{} weights still split, e.g. {:?}", split.len(), &split[..split.len().min(8)]);
            }

            /// F1: inside a mixed graph the pair joins, while the triangle and the edge-free node stay as MCL leaves
            /// them, in the shipped order (size desc, then first member). A pair's fate also depended on the REST of the
            /// graph, because convergence is tested on all columns at once: here 0.90 splits and 0.99 joins.
            #[test]
            fn a_pair_joins_inside_a_mixed_graph_and_the_rest_is_untouched() {
                let mut g = HomologyGraph::default();
                for i in 0..6u64 {
                    g.genes.push(("c".to_string(), i * 100 + 1, i * 100 + 50));
                }
                for (a, b) in [(0, 1), (0, 2), (1, 2)] {
                    g.edges.insert((a, b), 0.95);
                }
                g.edges.insert((3, 4), 0.90); // splits in this graph before the fix
                assert_eq!(mcl(&g, 2.8, 100, 1e-9), vec![vec![0, 1, 2], vec![3, 4], vec![5]]);
            }
```

- [ ] **Step 4: Run them, see the red**

Run:
`RLOCK_TIMEOUT=590 CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target TMPDIR=/mnt/linuxdisk/tmp bash tools/rlock.sh heavy bash -c "cargo test --release --lib -- a_connected_pair a_pair_joins > /mnt/linuxdisk/tmp/o1_close/t1_red.log 2>&1"; grep -E "^test |left:|test result" /mnt/linuxdisk/tmp/o1_close/t1_red.log`

Expected:
- both tests FAIL;
- `a_connected_pair_is_one_family_at_every_weight` fails at `w = 0.99` with `left: [[0], [1]]`;
- `a_pair_joins_inside_a_mixed_graph_and_the_rest_is_untouched` fails with `left: [[0, 1, 2], [3], [4], [5]]`.

This was checked on a standalone copy of `mcl` / `attractors`.

- [ ] **Step 5: Implement**

In `pub fn mcl(`, replace its last line `            attractors(&col, n)` with:

```rust
            join_two_node_components(&g.edges, attractors(&col, n))
```

Add right after `fn attractors` (its closing brace, ~60982):

```rust

        /// PREREG_mcl_two_node_2026-10-13: a connected component of exactly two nodes is ONE cluster. On two nodes MCL's
        /// column converges toward a tie between "self" and "mate", and the 1e-7 convergence test stopped at a residual
        /// whose SIGN chose join or split (0.99 split, 0.98 joined). Flow never crosses components, so a pair's fate is
        /// settled here and every other cluster passes through unchanged, in the same order.
        fn join_two_node_components(
            edges: &BTreeMap<(usize, usize), f64>,
            clusters: Vec<Vec<usize>>,
        ) -> Vec<Vec<usize>> {
            let mut deg: BTreeMap<usize, usize> = BTreeMap::new();
            for &(a, b) in edges.keys() {
                if a != b {
                    *deg.entry(a).or_default() += 1;
                    *deg.entry(b).or_default() += 1;
                }
            }
            let mut mate: BTreeMap<usize, usize> = BTreeMap::new();
            for &(a, b) in edges.keys() {
                if a != b && deg[&a] == 1 && deg[&b] == 1 {
                    mate.insert(a, b);
                    mate.insert(b, a);
                }
            }
            if mate.is_empty() {
                return clusters;
            }
            let mut out: Vec<Vec<usize>> = Vec::with_capacity(clusters.len());
            for c in clusters {
                match (c.len(), c.first().and_then(|x| mate.get(x))) {
                    (1, Some(&m)) if m > c[0] => out.push(vec![c[0], m]),
                    (1, Some(_)) => {}
                    _ => out.push(c),
                }
            }
            out.sort_by(|a, b| b.len().cmp(&a.len()).then(a[0].cmp(&b[0])));
            out
        }
```

Add this test after the two from Step 3:

```rust

            /// F1: the rule reads the GRAPH, not MCL's output: nodes of a path or a triangle that MCL left as singletons
            /// are not pairs (each has two edges), an edge-free node stays alone, and a joined pair takes its sorted place.
            #[test]
            fn the_pair_rule_touches_only_two_node_components() {
                let mut edges: BTreeMap<(usize, usize), f64> = BTreeMap::new();
                for (a, b) in [(0, 1), (1, 2), (3, 4), (3, 5), (4, 5)] {
                    edges.insert((a, b), 0.9);
                }
                let singles: Vec<Vec<usize>> = (0..7).map(|i| vec![i]).collect();
                assert_eq!(join_two_node_components(&edges, singles.clone()), singles);
                edges.insert((7, 8), 0.9);
                let mut want: Vec<Vec<usize>> = vec![vec![7, 8]];
                want.extend((0..7).map(|i| vec![i]));
                assert_eq!(join_two_node_components(&edges, (0..9).map(|i| vec![i]).collect()), want);
            }
```

- [ ] **Step 6: Run the lib suite, then the full suite**

```bash
RLOCK_TIMEOUT=590 CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target TMPDIR=/mnt/linuxdisk/tmp bash tools/rlock.sh heavy bash -c "cargo test --release --lib > /mnt/linuxdisk/tmp/o1_close/t1_lib.log 2>&1"; grep -E "test result|FAILED|panicked" /mnt/linuxdisk/tmp/o1_close/t1_lib.log
RLOCK_TIMEOUT=590 CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target TMPDIR=/mnt/linuxdisk/tmp bash tools/rlock.sh heavy bash -c "cargo test --release --no-fail-fast > /mnt/linuxdisk/tmp/o1_close/t1_full.log 2>&1"; grep -E "test result|FAILED" /mnt/linuxdisk/tmp/o1_close/t1_full.log
```

Expected:
- lib: `test result: ok. 1051 passed` (1,048 before + 3).
- full: every result line `ok` except `tests/pipeline_driver_dispatch.rs` with its 3 known failures (1,287 passed /
  3 failed before this task; now 1,290 / 3).
- If the full run exceeds 590 s, rerun it with `RLOCK_TIMEOUT=2400` in the background
  (`run_in_background`).
- Any OTHER failing test pinned a split pair. Report its name and assertion. Do not edit it before the controller
  rules.

- [ ] **Step 7: Build, keep the fixed binary, famsim after**

```bash
RLOCK_TIMEOUT=590 CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target TMPDIR=/mnt/linuxdisk/tmp bash tools/rlock.sh heavy bash -c "cargo build --release > /mnt/linuxdisk/tmp/o1_close/t1_build.log 2>&1"; echo "exit $?"
cp /mnt/linuxdisk/home/juanfraitu/rustle_target/release/mcl_families /mnt/linuxdisk/tmp/o1_close/mcl_families.after
F=/mnt/linuxdisk/tmp/o1_close/famsim
RLOCK_TIMEOUT=590 bash tools/rlock.sh heavy python3 bench/famsim ladder --template-spec bench/famsim/examples/ladder_human_SNRPB.json --out $F/after_human --only snp_0.01,snp_0.1,shuffle > $F/after_human.log 2>&1
RLOCK_TIMEOUT=590 bash tools/rlock.sh heavy python3 bench/famsim ladder --template-spec bench/famsim/examples/ladder_ppy_SPAG7.json --out $F/after_ppy --only snp_0.01,snp_0.1 > $F/after_ppy.log 2>&1
cut -f1,4-7 $F/after_human/ladder.tsv $F/after_ppy/ladder.tsv
```

Expected:
- every listed rung shows `denovo_intact 1 guided_intact 1` and `denovo_F` / `guided_F` = 1.0000;
- `cmp` of `mcl_families.before` and `.after` differs (the binary changed).

- [ ] **Step 8: Record F1, F3, F4 and commit**

Replace `(pending)` in the prereg's Results with three lines:
- F1: the red messages from Step 4, then green.
- F3: the before and after ladder rows.
- F4: the full-suite line.

Write "F2: Task 2" as a fourth line.

```bash
git add src/family.rs docs/archive/2026-10/merged/PREREG_mcl_two_node_2026-10-13.md
git commit -m "mcl: a connected two-node component is one cluster (PREREG_mcl_two_node_2026-10-13 F1/F3/F4; red-first tests)"
```

---

### Task 2: re-cluster the four default catalogs (F2)

> **Note (2026-10-13, after the final review):** the code embedded in this task predates the review fix rounds. The
> repo code is authoritative.

**Files:**
- Create: `bench/o1_hier/pair_diff.py`, `bench/o1_hier/test_pair_diff.py`, `bench/o1_hier/recluster.sh`
- Modify (results only): `docs/archive/2026-10/merged/PREREG_mcl_two_node_2026-10-13.md`; `docs/NEGATIVE_RESULTS_REGISTER.md` (append a
  row); `docs/THESIS_OBJECTIVES.md` (the O1 row of the STATUS 2026-10-11 table); `docs/archive/2026-10/merged/PENDING_2026-10-11.md`

**Interfaces:**
- Consumes: `/mnt/linuxdisk/tmp/o1_close/mcl_families.{before,after}` (Task 1).
- Produces:
  - `/mnt/linuxdisk/tmp/o1_close/recluster/{before,after}/<S>.fam.{clusters,loci,params}.tsv` for S in
    `human_A119b human_testis chimp_PTR gorilla_KB3781`;
  - `recluster/<S>.pair_diff.json`.
  - Task 6 compares its families with `recluster/after/<S>.fam.clusters.tsv`.

- [ ] **Step 1: Write the failing tests** — `bench/o1_hier/test_pair_diff.py`:

```python
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
import pair_diff as P  # noqa: E402

V3, V2 = ("3", "1.0000", "0.1000", "NA"), ("2", "1.0000", "0.0000", "NA")


def test_two_joined_singletons_are_the_only_allowed_change():
    before = {frozenset("abc"): V3, frozenset("gh"): V2}
    after = {frozenset("abc"): V3, frozenset("de"): V2, frozenset("gh"): V2}
    r = P.diff(before, after)
    assert (r["pairs_joined"], r["pairs_split"], r["other_changes"], r["unchanged"]) == (1, 0, 0, 2)
    assert (r["clusters_before"], r["clusters_after"], r["members_before"], r["members_after"]) == (2, 3, 5, 7)


def test_a_changed_member_set_is_an_other_change():
    r = P.diff({frozenset("abc"): V3}, {frozenset("ab"): V2, frozenset("cd"): V2})
    assert r["other_changes"] == 3 and r["pairs_joined"] == 0


def test_the_same_set_with_other_values_is_an_other_change():
    assert P.diff({frozenset("ab"): V2}, {frozenset("ab"): ("2", "0.5000", "0.0000", "NA")})["other_changes"] == 1


def test_a_pair_that_disappears_is_a_split():
    r = P.diff({frozenset("ab"): V2, frozenset("abc"[2:] + "xy"): V3}, {frozenset("cxy"): V3})
    assert (r["pairs_split"], r["other_changes"]) == (1, 0)


def test_load_groups_rows_by_cluster_and_ignores_ids(tmp_path):
    h = "cluster_id\tsize\tdensity\tfrac_in\tcorroborated\tchrom\tstart\tend\n"
    a, b = tmp_path / "a.tsv", tmp_path / "b.tsv"
    a.write_text(h + "MCL0\t2\t1.0000\t0.0000\tNA\tchr1\t1\t10\nMCL0\t2\t1.0000\t0.0000\tNA\tchr1\t21\t30\n")
    b.write_text(h + "MCL7\t2\t1.0000\t0.0000\tNA\tchr1\t21\t30\nMCL7\t2\t1.0000\t0.0000\tNA\tchr1\t1\t10\n")
    assert P.load(str(a)) == P.load(str(b)) == {frozenset({"chr1:1-10", "chr1:21-30"}): V2}
```

Run `python3 -m pytest -q bench/o1_hier/test_pair_diff.py`. Expected: `ModuleNotFoundError: No module named
'pair_diff'`.

- [ ] **Step 2: Implement** — `bench/o1_hier/pair_diff.py`:

```python
#!/usr/bin/env python3
"""PREREG_mcl_two_node_2026-10-13 F2: two clusters.tsv compared as member SETS (cluster ids are ranks, so they renumber).

    pair_diff.py BEFORE.tsv AFTER.tsv OUT.json

clusters.tsv lists clusters with >= min-size (2) members, one row per member. Allowed: a 2-member cluster in AFTER
whose members are absent from BEFORE (two singletons joined). Every other difference — a member set present on one side
only, or the same set with other size/density/frac_in/corroborated values — is an `other_change`; a 2-member set of
BEFORE whose members are absent from AFTER is a `pair_split`.
"""
import collections
import json
import sys

VALUES = ("size", "density", "frac_in", "corroborated")


def load(path):
    """-> {frozenset(member 'chrom:start-end'): (size, density, frac_in, corroborated)}."""
    members, vals = collections.defaultdict(list), {}
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ix = {k: hdr.index(k) for k in ("cluster_id", "chrom", "start", "end") + VALUES}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            cid = f[ix["cluster_id"]]
            members[cid].append(f"{f[ix['chrom']]}:{f[ix['start']]}-{f[ix['end']]}")
            vals[cid] = tuple(f[ix[k]] for k in VALUES)
    return {frozenset(m): vals[c] for c, m in members.items()}


def diff(before, after):
    in_before = frozenset().union(*before) if before else frozenset()
    in_after = frozenset().union(*after) if after else frozenset()
    r = collections.Counter(unchanged=0, pairs_joined=0, pairs_split=0, other_changes=0)
    for s, v in after.items():
        if s in before:
            r["unchanged" if before[s] == v else "other_changes"] += 1
        elif len(s) == 2 and not s & in_before:
            r["pairs_joined"] += 1
        else:
            r["other_changes"] += 1
    for s in before:
        if s not in after:
            r["pairs_split" if len(s) == 2 and not s & in_after else "other_changes"] += 1
    r["clusters_before"], r["clusters_after"] = len(before), len(after)
    r["members_before"], r["members_after"] = len(in_before), len(in_after)
    return dict(r)


def main(argv=None):
    a = argv if argv is not None else sys.argv[1:]
    res = diff(load(a[0]), load(a[1]))
    res["verdict"] = "PASS" if res["pairs_split"] == 0 and res["other_changes"] == 0 else "FAIL"
    with open(a[2], "w") as fh:
        json.dump(res, fh, indent=1)
    print(json.dumps(res))


if __name__ == "__main__":
    sys.exit(main())
```

Run `python3 -m pytest -q bench/o1_hier/test_pair_diff.py`. Expected: `5 passed`.

- [ ] **Step 3: `bench/o1_hier/recluster.sh`**

```bash
#!/bin/bash
# Re-cluster one default catalog from its cached all-vs-all (docs/archive/2026-10/merged/PREREG_mcl_two_node_2026-10-13.md F2).
#   recluster.sh BIN SAMPLE OUTDIR
# --fasta / --emit-units are not passed: neither changes clusters.tsv or loci.tsv (human_testis: byte-identical, 5.7 s).
set -euo pipefail
BIN=$1; S=$2; O=$3
R=/mnt/linuxdisk/tmp/rustle_figures/runs/$S
mkdir -p "$O"
"$BIN" --paf "$R/$S.fam.loci.paf" --gff "$R/$S.fam.loci.gff3" \
  --min-exonic-bp 1 --min-shared-exon-frac 0.60 --out "$O/$S.fam" > "$O/$S.log" 2>&1
```

- [ ] **Step 4: H0 — the pre-fix binary reproduces every shipped catalog**

```bash
C=/mnt/linuxdisk/tmp/o1_close
for S in human_testis chimp_PTR gorilla_KB3781 human_A119b; do
  RLOCK_TIMEOUT=590 bash tools/rlock.sh heavy bash bench/o1_hier/recluster.sh $C/mcl_families.before $S $C/recluster/before
  R=/mnt/linuxdisk/tmp/rustle_figures/runs/$S
  cmp $C/recluster/before/$S.fam.clusters.tsv $R/$S.fam.clusters.tsv && cmp $C/recluster/before/$S.fam.loci.tsv $R/$S.fam.loci.tsv && echo "$S H0 ok"
done
```

Expected: `H0 ok` for all four. A119b's PAF is 2.0 GB: if its run exceeds 590 s, rerun that one in the background
with `RLOCK_TIMEOUT=3000`. If any `cmp` differs, STOP and report the first differing lines; do not continue.

- [ ] **Step 5: The fixed binary, and the diff**

```bash
C=/mnt/linuxdisk/tmp/o1_close
for S in human_testis chimp_PTR gorilla_KB3781 human_A119b; do
  RLOCK_TIMEOUT=590 bash tools/rlock.sh heavy bash bench/o1_hier/recluster.sh $C/mcl_families.after $S $C/recluster/after
  python3 bench/o1_hier/pair_diff.py $C/recluster/before/$S.fam.clusters.tsv $C/recluster/after/$S.fam.clusters.tsv $C/recluster/$S.pair_diff.json
done
```

Expected: `"verdict": "PASS"` (`pairs_split` 0, `other_changes` 0) for every catalog, with `pairs_joined` > 0 in at
least the two human catalogs.
- If a catalog FAILS, STOP and report its JSON plus three example member sets from each side.
- Do not loosen the checker.

- [ ] **Step 6: Record F2, the new counts, and commit**

1. **Prereg Results:** replace the "F2: Task 2" line with:
   - one line per catalog: `pairs_joined`, clusters before -> after, members before -> after;
   - the H0 verdict.
2. **Register:** append to `docs/NEGATIVE_RESULTS_REGISTER.md` the row
   `| 1276 | 2026-10-13 | MCL two-node components (docs/archive/2026-10/merged/PREREG_mcl_two_node_2026-10-13.md) | A connected pair's existence was decided by the sign of a 1e-7 convergence residual (and by unrelated components) | <verdict: F1-F4 with the per-catalog joins> |`.
   Check the last row number first; use the next free one.
3. **THESIS_OBJECTIVES:** in the O1 row of the STATUS 2026-10-11 table, after the sentence listing the four catalogs'
   families / loci, insert:
   ` After the MCL two-node fix (2026-10-13, docs/archive/2026-10/merged/PREREG_mcl_two_node_2026-10-13.md; re-clustered from the cached all-vs-all, shipped files unchanged): A119b <a> / <b>, testis <c> / <d>, chimp <e> / <f>, gorilla <g> / <h>.`
   The numbers are `clusters_after` / `members_after`.
4. **PENDING:** append to `docs/archive/2026-10/merged/PENDING_2026-10-11.md` under a new heading `## O1 two-node fix (2026-10-13)` one
   bullet:
   - the shipped `runs/*` catalogs predate the fix;
   - the next families-stage run picks it up;
   - figure rescoring (paused) must use post-fix catalogs.

```bash
git add bench/o1_hier/pair_diff.py bench/o1_hier/test_pair_diff.py bench/o1_hier/recluster.sh docs/archive/2026-10/merged/PREREG_mcl_two_node_2026-10-13.md docs/NEGATIVE_RESULTS_REGISTER.md docs/THESIS_OBJECTIVES.md docs/archive/2026-10/merged/PENDING_2026-10-11.md
git commit -m "MCL two-node fix: four default catalogs re-clustered (F2: only split pairs change); register + O1 counts"
```

---

### Task 3: the v1.0 -> v2.0 lift and the Soto truth table

> **Note (2026-10-13, after the final review):** the code embedded in this task predates the review fix rounds. The
> repo code is authoritative.

**Files:**
- Create: `bench/o1_hier/v1_to_v2.py`, `bench/o1_hier/test_v1_to_v2.py`, `bench/o1_hier/soto_truth.py`,
  `bench/o1_hier/test_soto_truth.py`

**Interfaces:**
- Produces in `/mnt/linuxdisk/tmp/o1_close/lift/`:
  - `dupmasker_colors.v2.bed`, `sedef_pairs.v2.bed` (hierarchy inputs, Task 6), with `dupmasker.json` and
    `sedef.json` (lift reports with sequence verification);
  - `soto_truth.v1.tsv` (Task 4's regression check) and `soto_truth.v2.tsv` (Task 6's truth), with
    `soto_truth.v{1,2}.json`.
- Functions:
  - `v1_to_v2.read_chain(lines) -> {chrom: [(v1_start, v2_start, size)]}`;
  - `v1_to_v2.Lifter(chain).lift(chrom, s, e) -> (s2, e2, crosses_gap) | None` (0-based half-open);
  - `soto_truth.cat_genes`, `soto_truth.build`, `soto_truth.s1c_rows`.

- [ ] **Step 1: Write the failing lift tests** — `bench/o1_hier/test_v1_to_v2.py`:

```python
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
import v1_to_v2 as L  # noqa: E402

# v1 chr16 1000 bp, v2 995 bp: 500 aligned, 5 bp only in v1, 495 aligned
CHAIN = """chain 995 chr16 1000 + 0 1000 chr16 995 + 0 995 1
500 5 0
495

chain 1000 chr1 1000 + 0 1000 chr1 1000 + 0 1000 2
1000
"""


def lifter():
    return L.Lifter(L.read_chain(CHAIN.splitlines()))


def test_chain_blocks():
    assert L.read_chain(CHAIN.splitlines()) == {"chr16": [(0, 0, 500), (505, 500, 495)], "chr1": [(0, 0, 1000)]}


def test_identity_chromosome_lifts_unchanged():
    assert lifter().lift("chr1", 100, 200) == (100, 200, False)


def test_shift_after_a_deletion():
    assert lifter().lift("chr16", 600, 700) == (595, 695, False)


def test_an_interval_spanning_the_indel_lifts_with_it_and_is_flagged():
    assert lifter().lift("chr16", 400, 600) == (400, 595, True)


def test_an_end_inside_a_gap_is_unliftable():
    assert lifter().lift("chr16", 502, 510) is None
    assert lifter().lift("chr16", 490, 503) is None


def test_unknown_chromosome_is_unliftable():
    assert lifter().lift("chrY", 1, 5) is None


def test_rows_pairs_and_one_based():
    lines = ["chr16\t601\t700\tID_1\n", "chrY\t1\t5\tID_2\n"]
    out, c, _ = L.lift_rows(lifter(), lines, one_based=True)
    assert out == ["chr16\t596\t695\tID_1\n"] and c["rows_dropped"] == 1 and c["unliftable|chrY"] == 1
    out, c, _ = L.lift_rows(lifter(), ["chr1\t10\t20\tchr16\t600\t700\tS\n", "chr1\t10\t20\tchr16\t502\t510\tS\n"],
                            pairs=True)
    assert out == ["chr1\t10\t20\tchr16\t595\t695\tS\n"] and c["rows_dropped"] == 1
    out, _, _ = L.lift_rows(lifter(), ["chr16\t600\t700\tSD1\t0\t+\t600\t700\t1,1,1\n"])
    assert out == ["chr16\t595\t695\tSD1\t0\t+\t595\t695\t1,1,1\n"]
```

Run `python3 -m pytest -q bench/o1_hier/test_v1_to_v2.py`. Expected: `ModuleNotFoundError: No module named
'v1_to_v2'`.

- [ ] **Step 2: Implement** — `bench/o1_hier/v1_to_v2.py`:

```python
#!/usr/bin/env python3
"""CHM13 v1.0 -> v2.0 coordinates through the T2T v1.0 -> v1.1 chain (v2.0 = v1.1 + the HG002 chrY).

    v1_to_v2.py lift CHAIN IN OUT [--pairs] [--one-based] [--report R.json] [--verify V1_FA V2_FA [--n 1000] [--seed 7]]

An interval lifts only if both of its ends fall inside aligned blocks of its chromosome's chain. An end inside a chain
gap (an indel, the rDNA arrays) or on a chromosome without a chain (GRCh38 chrY in the plus38Y files) makes it
unliftable: the row is dropped and counted, never shifted. An interval whose ends sit in different blocks spans an
indel; it lifts with the indel and is counted (`crosses_gap`). `--pairs` lifts columns 1-3 and 4-6 (SEDEF pairs; the row
is dropped if either side fails); `--one-based` reads and writes GFF-style closed starts.
"""
import argparse
import bisect
import collections
import json
import random
import subprocess
import sys


def read_chain(lines):
    """{chrom: [(v1_start, v2_start, size)]} — one + strand chain per chromosome, same name on both sides (asserted)."""
    out, cur, t, q = {}, None, 0, 0
    for line in lines:
        f = line.split()
        if not f:
            continue
        if f[0] == "chain":
            assert f[2] == f[7] and f[4] == f[9] == "+", line
            assert f[2] not in out, f"more than one chain for {f[2]}"
            cur, t, q = out.setdefault(f[2], []), int(f[5]), int(f[10])
            continue
        size = int(f[0])
        cur.append((t, q, size))
        if len(f) == 3:
            t += size + int(f[1])
            q += size + int(f[2])
    return out


class Lifter:
    def __init__(self, chain):
        self.blocks = chain
        self.starts = {c: [b[0] for b in bs] for c, bs in chain.items()}

    def _block(self, c, x):
        i = bisect.bisect_right(self.starts[c], x) - 1
        b = self.blocks[c]
        return i if i >= 0 and x < b[i][0] + b[i][2] else None

    def lift(self, chrom, s, e):
        """0-based half-open [s, e) on v1.0 -> (s2, e2, crosses_gap) on v2.0, or None."""
        if chrom not in self.blocks or e <= s:
            return None
        i, j = self._block(chrom, s), self._block(chrom, e - 1)
        if i is None or j is None:
            return None
        b = self.blocks[chrom]
        s2 = s - b[i][0] + b[i][1]
        e2 = e - 1 - b[j][0] + b[j][1] + 1
        return (s2, e2, i != j) if e2 > s2 else None


def lift_rows(lifter, lines, pairs=False, one_based=False):
    """-> (lifted lines, counts, [(v1 region, v2 region)] per lifted interval) for TSV/BED lines."""
    off = 1 if one_based else 0
    out, c, regions = [], collections.Counter(), []
    for line in lines:
        if line.startswith(("#", "track", "browser")) or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        new, got = list(f), []
        for ci, si, ei in [(0, 1, 2)] + ([(3, 4, 5)] if pairs else []):
            r = lifter.lift(f[ci], int(f[si]) - off, int(f[ei]))
            if r is None:
                c[f"unliftable|{f[ci]}"] += 1
                got = None
                break
            new[si], new[ei] = str(r[0] + off), str(r[1])
            c["crosses_gap"] += r[2]
            got.append(((f[ci], int(f[si]) - off, int(f[ei])), (f[ci], r[0], r[1])))
        if got is None:
            c["rows_dropped"] += 1
            continue
        if not pairs and len(f) >= 8 and f[6] == f[1] and f[7] == f[2]:  # BED thickStart/thickEnd follow the span
            new[6], new[7] = new[1], new[2]
        out.append("\t".join(new) + "\n")
        c["rows_lifted"] += 1
        regions.extend(got)
    return out, dict(sorted(c.items())), regions


def fetch(fa, regions):
    """samtools faidx -> {region string: upper-case sequence}; regions are 0-based half-open (chrom, s, e)."""
    names = [f"{c}:{s + 1}-{e}" for c, s, e in regions]
    seqs, cur = {}, None
    for i in range(0, len(names), 500):
        res = subprocess.run(["samtools", "faidx", fa] + names[i:i + 500], capture_output=True, text=True, check=True)
        for line in res.stdout.splitlines():
            if line.startswith(">"):
                cur = line[1:]
                seqs[cur] = []
            else:
                seqs[cur].append(line.upper())
    return {k: "".join(v) for k, v in seqs.items()}


def verify(v1_fa, v2_fa, pairs_of_regions, n, seed):
    """Sequence identity of n sampled lifted intervals (those not spanning an indel)."""
    rng = random.Random(seed)
    sample = [p for p in pairs_of_regions if (p[0][2] - p[0][1]) == (p[1][2] - p[1][1])]
    sample = rng.sample(sample, min(n, len(sample)))
    a, b = fetch(v1_fa, [p[0] for p in sample]), fetch(v2_fa, [p[1] for p in sample])
    same = sum(a[f"{p[0][0]}:{p[0][1] + 1}-{p[0][2]}"] == b[f"{p[1][0]}:{p[1][1] + 1}-{p[1][2]}"] for p in sample)
    return {"checked": len(sample), "identical": same, "fraction": same / len(sample) if sample else None}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("lift")
    p.add_argument("chain"); p.add_argument("src"); p.add_argument("dst"); p.add_argument("--report")
    p.add_argument("--pairs", action="store_true"); p.add_argument("--one-based", action="store_true")
    p.add_argument("--verify", nargs=2, metavar=("V1_FA", "V2_FA"))
    p.add_argument("--n", type=int, default=1000); p.add_argument("--seed", type=int, default=7)
    a = ap.parse_args(argv)
    lifter = Lifter(read_chain(open(a.chain)))
    out, counts, regions = lift_rows(lifter, open(a.src), a.pairs, a.one_based)
    with open(a.dst, "w") as fh:
        fh.writelines(out)
    rep = {"src": a.src, "dst": a.dst, "counts": counts}
    if a.verify:
        rep["verify"] = verify(a.verify[0], a.verify[1], regions, a.n, a.seed)
    if a.report:
        with open(a.report, "w") as fh:
            json.dump(rep, fh, indent=1)
    print(json.dumps(rep, indent=1))


if __name__ == "__main__":
    sys.exit(main())
```

Run the tests. Expected: `7 passed`.

- [ ] **Step 3: Write the failing truth tests** — `bench/o1_hier/test_soto_truth.py`:

```python
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
import soto_truth as S  # noqa: E402
import v1_to_v2 as L  # noqa: E402


def bed(chrom, s, e, strand, sizes, offs, gid):
    f = [chrom, str(s), str(e), "x", "0", strand, str(s), str(e), "0", str(len(sizes)),
         ",".join(map(str, sizes)) + ",", ",".join(map(str, offs)) + ","] + ["."] * 6 + [gid]
    return "\t".join(f) + "\n"


def test_exons_are_the_union_of_the_transcripts():
    lines = [bed("chr1", 100, 400, "+", [50, 100], [0, 200], "G1"), bed("chr1", 120, 500, "+", [60, 80], [0, 300], "G1")]
    assert S.cat_genes(lines) == {"G1": ("chr1", "+", [(100, 180), (300, 400), (420, 500)])}


def test_rows_join_by_gene_id_and_write_gff_coordinates():
    genes = {"G1": ("chr1", "-", [(100, 180), (300, 500)])}
    out, c = S.build([("ID_1", "G1"), ("ID_1", "G9")], genes)
    assert out == [("ID_1", "G1", "chr1", 101, 500, "-", "101-180,301-500")]
    assert c["genes_not_in_cat"] == 1 and c["labels"] == 1


def test_lift_drops_an_exon_in_a_gap_and_shifts_the_rest():
    chain = L.read_chain("chain 1 chr1 1000 + 0 1000 chr1 995 + 0 995 1\n500 5 0\n495\n".splitlines())
    genes = {"G1": ("chr1", "+", [(498, 503), (600, 700)])}
    out, c = S.build([("ID_1", "G1")], genes, L.Lifter(chain))
    assert out == [("ID_1", "G1", "chr1", 596, 695, "+", "596-695")] and c["exons_unliftable"] == 1


def test_s1c_drops_unassigned(tmp_path):
    p = tmp_path / "s1c.tsv"
    p.write_text("Gene Name\tGene ID\tFamily ID\nA\tG1\tID_1\nB\tG2\tUnassigned_3\n")
    assert S.s1c_rows(str(p)) == [("ID_1", "G1")]
```

Run `python3 -m pytest -q bench/o1_hier/test_soto_truth.py`. Expected: `ModuleNotFoundError: No module named
'soto_truth'`.

- [ ] **Step 4: Implement** — `bench/o1_hier/soto_truth.py`:

```python
#!/usr/bin/env python3
"""Soto 2025 families as a truth table with exons: S1C (Gene ID -> Family ID) joined to CAT v4 (CHM13 v1.0 BED12) by
Gene ID (column 19), never by name; optionally lifted to CHM13 v2.0 through the T2T v1.0 -> v1.1 chain.

    soto_truth.py S1C.tsv CAT_V4.bed OUT.tsv [--chain CHAIN] [--report R.json]

OUT columns: label gene_id chrom start end strand exons — GFF 1-based closed; `exons` = the union of the exon blocks of
every transcript of the gene, `s-e` comma-joined; start/end = their hull. `Unassigned_*` labels are not families and
are dropped. With --chain an exon with an end in a chain gap is dropped and counted, and a gene left without exons is
dropped and counted.
"""
import argparse
import collections
import csv
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from v1_to_v2 import Lifter, read_chain  # noqa: E402


def merge(blocks):
    """0-based half-open blocks -> their union, sorted."""
    out = []
    for s, e in sorted(blocks):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [tuple(b) for b in out]


def cat_genes(lines):
    """CAT v4 BED12 -> {gene_id: (chrom, strand, [0-based exon blocks, unioned])}."""
    acc = {}
    for line in lines:
        f = line.rstrip("\n").split("\t")
        if len(f) < 19:
            continue
        chrom, s, strand, gid = f[0], int(f[1]), f[5], f[18]
        sizes = [int(x) for x in f[10].rstrip(",").split(",")]
        offs = [int(x) for x in f[11].rstrip(",").split(",")]
        g = acc.setdefault(gid, (chrom, strand, []))
        assert g[0] == chrom and g[1] == strand, f"{gid}: transcripts on different chromosomes or strands"
        g[2].extend((s + o, s + o + z) for o, z in zip(offs, sizes))
    return {gid: (c, st, merge(b)) for gid, (c, st, b) in acc.items()}


def s1c_rows(path):
    """-> [(Family ID, Gene ID)] for the families proper (no `Unassigned_*`)."""
    with open(path) as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    h = rows[0]
    gi, fi = h.index("Gene ID"), h.index("Family ID")
    return [(r[fi], r[gi]) for r in rows[1:] if len(r) > fi and r[fi] and not r[fi].startswith("Unassigned_")]


def build(rows, genes, lifter=None):
    out, c = [], collections.Counter()
    for label, gid in sorted(set(rows)):
        if gid not in genes:
            c["genes_not_in_cat"] += 1
            continue
        chrom, strand, ex = genes[gid]
        if lifter is not None:
            lifted = [lifter.lift(chrom, s, e) for s, e in ex]
            c["exons_unliftable"] += sum(r is None for r in lifted)
            ex = merge((r[0], r[1]) for r in lifted if r is not None)
            if not ex:
                c["genes_unliftable"] += 1
                continue
        out.append((label, gid, chrom, ex[0][0] + 1, ex[-1][1], strand, ",".join(f"{s + 1}-{e}" for s, e in ex)))
        c["rows"] += 1
    c["labels"] = len({r[0] for r in out})
    return out, dict(sorted(c.items()))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("s1c"); ap.add_argument("cat"); ap.add_argument("out")
    ap.add_argument("--chain"); ap.add_argument("--report")
    a = ap.parse_args(argv)
    with open(a.cat) as fh:
        genes = cat_genes(fh)
    lifter = Lifter(read_chain(open(a.chain))) if a.chain else None
    out, counts = build(s1c_rows(a.s1c), genes, lifter)
    with open(a.out, "w") as fh:
        fh.write("label\tgene_id\tchrom\tstart\tend\tstrand\texons\n")
        fh.writelines("\t".join(map(str, r)) + "\n" for r in out)
    rep = {"inputs": vars(a), "counts": counts}
    if a.report:
        with open(a.report, "w") as fh:
            json.dump(rep, fh, indent=1)
    print(json.dumps(rep, indent=1))


if __name__ == "__main__":
    sys.exit(main())
```

Run `python3 -m pytest -q bench/o1_hier`. Expected: all pass (`4 passed` for this file).

- [ ] **Step 5: Lift the hierarchy inputs with sequence verification, build the truth tables**

```bash
W=/mnt/linuxdisk/home/juanfraitu/winloci_data; C=$W/soto_replication/v1.0_to_v1.1_rdna_merged.chain; L=/mnt/linuxdisk/tmp/o1_close/lift
V1=$W/soto_replication/t2t-chm13-v1.0.fa.gz; V2=$W/chm13v2.0.fa; mkdir -p $L
RLOCK_TIMEOUT=300 bash tools/rlock.sh light python3 bench/o1_hier/v1_to_v2.py lift $C $W/duplicons/chm13.draft_v1.0_plus38Y_dupmasker_colors.bed $L/dupmasker_colors.v2.bed --verify $V1 $V2 --report $L/dupmasker.json
RLOCK_TIMEOUT=300 bash tools/rlock.sh light python3 bench/o1_hier/v1_to_v2.py lift $C $W/soto_replication/final_v1_clean.bed $L/sedef_pairs.v2.bed --pairs --verify $V1 $V2 --report $L/sedef.json
RLOCK_TIMEOUT=300 bash tools/rlock.sh light python3 bench/o1_hier/soto_truth.py bench/soto/soto_famCN_S1C.tsv $W/soto_replication/cat_v4.bed $L/soto_truth.v1.tsv --report $L/soto_truth.v1.json
RLOCK_TIMEOUT=300 bash tools/rlock.sh light python3 bench/o1_hier/soto_truth.py bench/soto/soto_famCN_S1C.tsv $W/soto_replication/cat_v4.bed $L/soto_truth.v2.tsv --chain $C --report $L/soto_truth.v2.json
```

Expected (measured while writing the plan):
- **DupMasker:** `rows_lifted` 199154; `rows_dropped` 6661 (`unliftable|chrY` 6306, the GRCh38 Y; the rest on
  chr13/14/15/21/22 rDNA regions); `verify` 1000 / 1000 identical.
- **SEDEF:** `rows_lifted` 150833, `rows_dropped` 925, `rows_crossing_gap` 156, `sides_crossing_gap` 247 (corrected in Task 3's fix round: the first expected value, 305, also counted sides of dropped rows); `verify` 1000 / 1000.
- **Truth v1 and v2:** each `rows` 2458, `labels` 491; v2 `exons_unliftable` 0.

A verify fraction below 1.0 or a count off by more than 1% means STOP and report.

- [ ] **Step 6: Commit**

```bash
git add bench/o1_hier/v1_to_v2.py bench/o1_hier/test_v1_to_v2.py bench/o1_hier/soto_truth.py bench/o1_hier/test_soto_truth.py
git commit -m "o1_hier: CHM13 v1.0 -> v2.0 chain lift (sequence-verified) and the Soto truth with exons (S1C x CAT v4 by Gene ID)"
```

---

### Task 4: the hierarchy-aware scorer (TDD) and its regression on the 10-09 run

> **Note (2026-10-13, after the final review):** the code embedded in this task predates the review fix rounds. The
> repo code is authoritative.

**Files:**
- Create: `bench/o1_hier/hier_match.py`, `bench/o1_hier/test_hier_match.py`, `bench/o1_hier/advisor_table.py`,
  `bench/o1_hier/test_advisor_table.py`

**Interfaces:**
- Consumes: `lift/soto_truth.v1.tsv` (Task 3).
- Produces:
  - the CLI `hier_match.py MEMBERS TRUTH OUT.json --loci-gff GFF [--fold LOCI_TSV] [--map exon|key] [--null N]
    [--seed S]`;
  - its JSON keys: `coverage`, `hierarchy`, `flat`, `exact_any`, `null.{global,within_family}.{p_exact,
    p_mean_jaccard,...}`, `per_truth`;
  - `advisor_table.py OUT.tsv NAME=JSON ...`.

- [ ] **Step 1: Write the failing tests** — `bench/o1_hier/test_hier_match.py`:

```python
import os
import random
import sys

sys.path.insert(0, os.path.dirname(__file__))
import hier_match as H  # noqa: E402

NODES = {
    ("family", "F1"): frozenset({"a", "b", "c", "d"}),
    ("expansion", "F1.x1"): frozenset({"a", "b"}),
    ("expansion", "F1.x2"): frozenset({"c", "d"}),
    ("family", "F2"): frozenset({"e", "f"}),
    ("expansion", "F2.x1"): frozenset({"e", "f"}),  # identical to its family: the tie policy picks the family
    ("superfamily", "SF1"): frozenset({"a", "b", "c", "d", "e", "f"}),
}
MEMBERS = sorted("abcdef")
FAMILY_OF = {"a": "F1", "b": "F1", "c": "F1", "d": "F1", "e": "F2", "f": "F2"}


def test_truth_one_level_down_is_found_at_the_expansion_level():
    truth = {"T1": frozenset({"a", "b"}), "T2": frozenset({"c", "d"})}
    m = {x.label: x for x in H.Scorer(NODES, MEMBERS).run(truth)}
    assert (m["T1"].level, m["T1"].node_id, m["T1"].jaccard) == ("expansion", "F1.x1", 1.0)
    assert (m["T2"].level, m["T2"].node_id) == ("expansion", "F1.x2")


def test_the_flat_score_misses_what_the_hierarchy_finds():
    truth = {"T1": frozenset({"a", "b"}), "T2": frozenset({"c", "d"})}
    fam = {k: v for k, v in NODES.items() if k[0] == "family"}
    assert H.summary(H.Scorer(fam, MEMBERS).run(truth))["exact"] == 0
    assert H.summary(H.Scorer(NODES, MEMBERS).run(truth))["exact"] == 2


def test_identical_sets_at_two_levels_are_matched_once_with_the_named_tie_policy():
    (m,) = H.Scorer(NODES, MEMBERS).run({"T": frozenset({"e", "f"})})
    assert (m.level, m.node_id) == ("family", "F2")


def test_one_to_one_gives_a_contested_node_to_one_truth_family_only():
    truth = {"T1": frozenset({"e", "f"}), "T2": frozenset({"e", "f"})}
    ms = H.Scorer({("family", "F2"): NODES[("family", "F2")]}, MEMBERS).run(truth)
    assert sorted(m.jaccard for m in ms) == [0.0, 1.0]
    assert H.exact_any(truth, NODES) == {"any": 2, "expansion": 2, "family": 2}


def R(chrom, s, e, strand, exons):
    return H.Record(chrom, s, e, strand, exons)


LOCI = {  # key: (chrom, strands, exons) as load_loci_gff returns them
    "chr1:100-400": ("chr1", frozenset("+"), [(100, 150), (350, 400)]),
    "chr1:1000-1400": ("chr1", frozenset("+"), [(1000, 1100), (1300, 1400)]),
    "chr1:2000-2100": ("chr1", frozenset("-"), [(2000, 2100)]),
    "chr1:3000-3100": ("chr1", frozenset("+"), [(3000, 3100)]),   # folded into chr1:1000-1400
    "chr1:5000-5100": ("chr1", frozenset("+"), [(5000, 5100)]),   # a catalog locus outside every family
}
MEM = ["chr1:100-400", "chr1:1000-1400", "chr1:2000-2100"]
FOLD = {"chr1:3000-3100": "chr1:1000-1400"}


def test_a_shared_exon_on_the_same_strand_joins_and_an_intron_or_the_other_strand_does_not():
    truth = {"T": [R("chr1", 120, 130, "+", [(120, 130)]),      # exon on exon, same strand: joins chr1:100-400
                   R("chr1", 200, 300, "+", [(200, 300)]),      # inside the intron only: no member
                   R("chr1", 2050, 2060, "+", [(2050, 2060)]),  # exon on exon, other strand: no member
                   R("chr1", 3050, 3060, "+", [(3050, 3060)])]}  # a folded record's exon: its representative
    sets, per, cov = H.truth_sets(truth, MEM, fold=FOLD, loci=LOCI)
    assert sets["T"] == frozenset({"chr1:100-400", "chr1:1000-1400"})
    assert [len(x) for x in per["T"]] == [1, 0, 0, 1]
    assert (cov["records_on_members"], cov["records_without_locus"]) == (2, 2)


def test_coverage_separates_an_unclustered_locus_from_no_locus_and_counts_empty_families():
    truth = {"A": [R("chr1", 120, 130, "+", [(120, 130)]), R("chr1", 5010, 5020, "+", [(5010, 5020)])],
             "B": [R("chr1", 9000, 9100, "+", [(9000, 9100)])]}
    sets, per, cov = H.truth_sets(truth, MEM, fold=FOLD, loci=LOCI)
    assert (cov["records_on_unclustered_locus"], cov["records_without_locus"]) == (1, 1)
    assert cov["labels_without_member"] == ["B"] and cov["labels_one_member"] == 1 and cov["labels_no_member"] == 1
    node = {("family", "F"): frozenset(MEM)}
    (m,) = H.Scorer(node, MEM).run({"A": sets["A"]})
    assert H.gene_recall(m, node, per) == 0.5  # the record on the unclustered locus counts against the family


def test_key_mode_needs_identical_coordinates(tmp_path):
    t = tmp_path / "t.tsv"
    t.write_text("chr1\t100\t400\tT1\nchr1\t3000\t3100\tT1\nchr1\t101\t400\tT2\nchr3\t1\t9\tUnassigned_4\n")
    truth = H.load_truth(str(t))
    assert sorted(truth) == ["T1", "T2"] and truth["T2"][0].exons == [(101, 400)]
    sets, _, _ = H.truth_sets(truth, MEM, fold=FOLD, how="key")
    assert sets["T1"] == frozenset({"chr1:100-400", "chr1:1000-1400"}) and sets["T2"] == frozenset()


def test_the_soto_truth_table_is_read_with_its_strand_and_exons(tmp_path):
    t = tmp_path / "t.tsv"
    t.write_text("label\tgene_id\tchrom\tstart\tend\tstrand\texons\nID_1\tG1\tchr1\t101\t500\t-\t101-180,301-500\n")
    assert H.load_truth(str(t)) == {"ID_1": [H.Record("chr1", 101, 500, "-", [(101, 180), (301, 500)])]}


def test_loci_gff_pools_exons_and_uses_the_span_without_exon_lines(tmp_path):
    g = tmp_path / "l.gff3"
    g.write_text("##gff-version 3\n"
                 "chr1\t.\tgene\t10\t90\t.\t+\t.\tID=gene-A;Name=A\n"
                 "chr1\t.\texon\t10\t20\t.\t+\t.\tParent=gene-A;gene=A\n"
                 "chr1\t.\texon\t80\t90\t.\t+\t.\tParent=gene-A;gene=A\n"
                 "chr2\t.\tgene\t5\t50\t.\t-\t.\tID=gene-B;Name=B\n")
    assert H.load_loci_gff(str(g)) == {"chr1:10-90": ("chr1", frozenset("+"), [(10, 20), (80, 90)]),
                                        "chr2:5-50": ("chr2", frozenset("-"), [(5, 50)])}


def test_permutations_keep_sizes_and_within_family_keeps_families():
    rng = random.Random(7)
    for kind in ("global", "within_family"):
        p = H.permutation(MEMBERS, FAMILY_OF, kind, rng)
        assert sorted(p) == MEMBERS and sorted(p.values()) == MEMBERS
        moved = H.relabel(NODES, p)
        assert all(len(moved[k]) == len(v) for k, v in NODES.items())
    p = H.permutation(MEMBERS, FAMILY_OF, "within_family", rng)
    assert all(FAMILY_OF[p[m]] == FAMILY_OF[m] for m in MEMBERS)
    moved = H.relabel(NODES, p)
    assert all(moved[k] == v for k, v in NODES.items() if k[0] in ("family", "superfamily"))


def test_self_match_is_perfect_and_the_null_is_not():
    truth = {k[1]: v for k, v in NODES.items() if k[0] == "expansion"}
    sc = H.Scorer(NODES, MEMBERS)
    assert H.summary(sc.run(truth))["exact"] == 3
    ex, mj = H.null(sc, truth, MEMBERS, FAMILY_OF, "global", 200, 7)
    assert len(ex) == 200 and max(ex) <= 3 and sum(ex) / len(ex) < 3
    assert H.p_value(ex, 3) > 1 / 201  # six members: chance hits are common, as they should be


def test_cli_end_to_end(tmp_path):
    mem = tmp_path / "m.tsv"
    mem.write_text("sd_complex\tsuperfamily\tfamily\tmember\texpansions\n"
                   "CX0\tSF0\tMCL0\tchr1:1-10\tMCL0.x1\nCX0\tSF0\tMCL0\tchr1:21-30\tMCL0.x1\n"
                   "CX0\tSF0\tMCL0\tchr1:41-50\tMCL0.x2\nCX0\tSF0\tMCL0\tchr1:61-70\tMCL0.x2\n")
    tr = tmp_path / "t.tsv"
    tr.write_text("chr1\t2\t5\tA\nchr1\t25\t26\tA\nchr1\t45\t46\tB\nchr1\t65\t69\tB\n")
    gff = tmp_path / "l.gff3"
    gff.write_text("".join(f"chr1\t.\tgene\t{s}\t{s + 9}\t.\t+\t.\tID=gene-{s}\n" for s in (1, 21, 41, 61)))
    out = tmp_path / "o.json"
    H.main([str(mem), str(tr), str(out), "--loci-gff", str(gff), "--null", "20"])
    import json
    r = json.load(open(out))
    assert r["hierarchy"]["exact"] == 2 and r["flat"]["exact"] == 0
    assert r["hierarchy"]["exact_by_level"] == {"expansion": 2}
    assert set(r["null"]) == {"global", "within_family"}
```

Run `python3 -m pytest -q bench/o1_hier/test_hier_match.py`. Expected: `ModuleNotFoundError: No module named
'hier_match'`.

- [ ] **Step 2: Implement** — `bench/o1_hier/hier_match.py`:

```python
#!/usr/bin/env python3
"""Hierarchy-aware O1 scorer (docs/superpowers/specs/2026-10-13-close-o1-design.md, section A).

    hier_match.py MEMBERS_TSV TRUTH_TSV OUT.json [--loci-gff LOCI_GFF3] [--fold LOCI_TSV] [--map exon|key]
                  [--null 1000] [--seed 7]        (--loci-gff is required by --map exon, the default)

Nodes: every family, expansion, superfamily and SD complex of `<out>.hierarchy.members.tsv`, each as its set of member
loci (`chrom:start-end`, GFF 1-based). Duplicons, SD blocks and external families are not nodes.
Truth: the `soto_truth.py` table (`label gene_id chrom start end strand exons`, GFF 1-based) or the 4-column
`chrom start end label` of `--hierarchy-external` (no strand; the span is its one exon); `Unassigned_*` labels are not
families. `--map exon` (default, de novo catalogs): a truth gene joins a member when one of its exons shares >= 1 bp,
on the same strand, with an exon of the member's locus or of a locus folded into it (`--fold`, the run's
`<out>.loci.tsv`; exons from `--loci-gff`, the catalog's `.fam.loci.gff3`). `--map key`: identical coordinates only,
directly or through the fold map (the 10-09 F_O universe, where the members ARE the annotation records).
A truth family is EVALUABLE with >= 2 members. Matching: one-to-one Hungarian (scipy linear_sum_assignment) on
-Jaccard; rows in label order, columns in (level rank, node id) order, rank family < expansion < superfamily <
sd_complex: the named tie policy. `exact_any` asks only whether an equal node exists (tie-free, not one-to-one).
Nulls keep every node's size and the nesting: `global` relabels all members by one random permutation per round;
`within_family` permutes members inside their own family, so families, superfamilies and SD complexes are unchanged
and only the expansions move.
"""
import argparse
import bisect
import collections
import json
import random
import sys

import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.sparse import csr_matrix

LEVELS = ("family", "expansion", "superfamily", "sd_complex")
RANK = {level: i for i, level in enumerate(LEVELS)}
Match = collections.namedtuple("Match", "label level node_id jaccard recall precision inter t_size n_size")
Record = collections.namedtuple("Record", "chrom start end strand exons")


def parse_key(k):
    c, r = k.rsplit(":", 1)
    s, e = r.split("-")
    return c, int(s), int(e)


def load_nodes(members_tsv):
    """-> ({(level, node_id): frozenset(member)}, {member: family})."""
    nodes, family_of = collections.defaultdict(set), {}
    with open(members_tsv) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ix = {k: hdr.index(k) for k in ("sd_complex", "superfamily", "family", "member", "expansions")}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            m = f[ix["member"]]
            family_of[m] = f[ix["family"]]
            for level in ("sd_complex", "superfamily", "family"):
                if f[ix[level]] not in ("", "-"):
                    nodes[(level, f[ix[level]])].add(m)
            for x in f[ix["expansions"]].split(","):
                if x not in ("", "-"):
                    nodes[("expansion", x)].add(m)
    return {k: frozenset(v) for k, v in nodes.items()}, family_of


def load_truth(path):
    """-> {label: [Record]} (GFF 1-based closed); `Unassigned_*` labels dropped."""
    out = collections.defaultdict(list)
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0] == "label":
                continue
            if len(f) >= 7:
                label, chrom, s, e, strand = f[0], f[2], int(f[3]), int(f[4]), f[5]
                exons = [tuple(int(x) for x in b.split("-")) for b in f[6].split(",")]
            elif len(f) >= 4:
                label, chrom, s, e, strand = f[3], f[0], int(f[1]), int(f[2]), "."
                exons = [(s, e)]
            else:
                continue
            if not label.startswith("Unassigned_"):
                out[label].append(Record(chrom, s, e, strand, exons))
    return dict(out)


def load_fold(path):
    """`<out>.loci.tsv` (annotation -> representative) -> {record key: member key}."""
    with open(path) as fh:
        fh.readline()
        return dict(line.rstrip("\n").split("\t")[:2] for line in fh if line.strip())


def load_loci_gff(path):
    """Every locus of a catalog (`.fam.loci.gff3`: gene lines + exon lines by Parent) -> {key: (chrom, strands, exons)};
    a locus without exon lines is its span (as `exonless_span`); loci sharing a key pool their exons and strands."""
    genes, exons = {}, collections.defaultdict(list)
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            at = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            if f[2] == "gene":
                genes[at["ID"]] = (f"{f[0]}:{f[3]}-{f[4]}", f[0], f[6], int(f[3]), int(f[4]))
            elif f[2] == "exon":
                exons[at["Parent"]].append((int(f[3]), int(f[4])))
    out = {}
    for gid, (key, chrom, strand, s, e) in genes.items():
        _, strands, ex = out.get(key, (chrom, frozenset(), []))
        out[key] = (chrom, strands | {strand}, ex + (exons.get(gid) or [(s, e)]))
    return out


class SpanIndex:
    """Closed intervals per chromosome with a payload; overlap queries by binary search."""

    def __init__(self, items):
        by = collections.defaultdict(list)
        for c, s, e, p in items:
            by[c].append((s, e, p))
        self.by, self.starts, self.maxlen = by, {}, {}
        for c, v in by.items():
            v.sort(key=lambda x: (x[0], x[1]))
            self.starts[c] = [x[0] for x in v]
            self.maxlen[c] = max(x[1] - x[0] for x in v)

    def query(self, c, s, e):
        v = self.by.get(c)
        if not v:
            return []
        lo = bisect.bisect_left(self.starts[c], s - self.maxlen[c])
        hi = bisect.bisect_right(self.starts[c], e)
        return [p for a, b, p in v[lo:hi] if a <= e and s <= b]


def _strand_ok(a, b):
    return "." in a or "." in b or bool(a & b)


def truth_sets(truth, members, fold=None, loci=None, how="exon"):
    """truth {label: [Record]}; members = member keys; fold {record key: member key}; loci = load_loci_gff (exon mode).
    -> (sets {label: frozenset(member)}, per_record {label: [frozenset(member)]}, coverage)."""
    mset = set(members)
    fold = {r: m for r, m in (fold or {}).items() if m in mset}
    cov, all_idx = collections.Counter(), None
    if how == "exon":
        items = []
        for key, m in [(m, m) for m in sorted(mset)] + sorted(fold.items()):
            if key not in loci:
                cov["member_loci_missing_from_gff"] += 1
                c, s, e = parse_key(key)
                items.append((c, s, e, (frozenset("."), m)))
                continue
            chrom, strands, ex = loci[key]
            items.extend((chrom, s, e, (strands, m)) for s, e in ex)
        idx = SpanIndex(items)
        all_idx = SpanIndex([(chrom, s, e, (strands, 1)) for chrom, strands, ex in loci.values() for s, e in ex])

        def hit(r, index=idx):
            return {p for s, e in r.exons for st, p in index.query(r.chrom, s, e) if _strand_ok(st, {r.strand})}
    else:
        def hit(r, index=None):
            k = f"{r.chrom}:{r.start}-{r.end}"
            return {k} if k in mset else ({fold[k]} if k in fold else set())
    sets, per_record, without = {}, {}, []
    for label in sorted(truth):
        recs = [frozenset(hit(r)) for r in truth[label]]
        per_record[label] = recs
        sets[label] = frozenset().union(*recs)
        cov["records"] += len(recs)
        for r, h in zip(truth[label], recs):
            if h:
                cov["records_on_members"] += 1
                cov["records_on_several_members"] += len(h) > 1
            elif all_idx is not None and hit(r, all_idx):
                cov["records_on_unclustered_locus"] += 1
            else:
                cov["records_without_locus" if all_idx is not None else "records_without_member"] += 1
        n = len(sets[label])
        cov["labels"] += 1
        cov["labels_evaluable" if n >= 2 else ("labels_one_member" if n == 1 else "labels_no_member")] += 1
        if n == 0:
            without.append(label)
    shared = collections.Counter(m for s in sets.values() for m in s)
    out = dict(sorted(cov.items()))
    out["members_in_several_truth_families"] = sum(1 for v in shared.values() if v > 1)
    out["labels_without_member"] = without
    return sets, per_record, out


def _incidence(sets, index):
    rows = [index[m] for s in sets for m in s]
    cols = [j for j, s in enumerate(sets) for _ in s]
    return csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(index), len(sets)))


class Scorer:
    """One-to-one Hungarian between truth sets and nodes on -Jaccard, with the named tie policy (rows in label
    order, columns in (level rank, node id) order)."""

    def __init__(self, nodes, members):
        self.keys = sorted(nodes, key=lambda k: (RANK[k[0]], k[1]))
        self.index = {m: i for i, m in enumerate(sorted(members))}
        self.sizes = np.array([len(nodes[k]) for k in self.keys], float)
        self.N = _incidence([nodes[k] for k in self.keys], self.index)

    def run(self, truth):
        labels = sorted(truth)
        out = {}
        if labels and self.keys:
            T = _incidence([truth[x] for x in labels], self.index)
            inter = (T.T @ self.N).tocsc()
            used = np.unique(inter.nonzero()[1])
            if used.size:
                I = inter[:, used].toarray()
                ts = np.array([len(truth[x]) for x in labels], float)
                J = I / (ts[:, None] + self.sizes[used][None, :] - I)
                for i, j in zip(*linear_sum_assignment(-J)):
                    if I[i, j] > 0:
                        level, nid = self.keys[used[j]]
                        n = self.sizes[used[j]]
                        out[labels[i]] = Match(labels[i], level, nid, float(J[i, j]), float(I[i, j] / ts[i]),
                                               float(I[i, j] / n), int(I[i, j]), int(ts[i]), int(n))
        return [out.get(x) or Match(x, None, None, 0.0, 0.0, 0.0, 0, len(truth[x]), 0) for x in labels]


def summary(matches, nodes=None, per_record=None):
    n = len(matches)
    hit = [m for m in matches if m.level]
    exact = sum(m.jaccard == 1.0 for m in matches)
    out = {
        "n_truth": n,
        "n_matched": len(hit),
        "exact": exact,
        "exact_rate": exact / n if n else None,
        "mean_jaccard": sum(m.jaccard for m in matches) / n if n else None,
        "macro_recall": sum(m.recall for m in matches) / n if n else None,
        "macro_precision": sum(m.precision for m in hit) / len(hit) if hit else None,
        "micro_recall": sum(m.inter for m in matches) / sum(m.t_size for m in matches) if n else None,
        "micro_precision": sum(m.inter for m in hit) / sum(m.n_size for m in hit) if hit else None,
        "by_level": dict(sorted(collections.Counter(m.level for m in hit).items())),
        "exact_by_level": dict(sorted(collections.Counter(m.level for m in matches if m.jaccard == 1.0).items())),
    }
    if nodes is not None and per_record is not None and n:
        out["macro_gene_recall"] = sum(gene_recall(m, nodes, per_record) for m in matches) / n
    return out


def gene_recall(m, nodes, per_record):
    """Share of the truth family's RECORDS (with or without a member) that have a member in the matched node."""
    recs = per_record[m.label]
    if not m.level or not recs:
        return 0.0
    node = nodes[(m.level, m.node_id)]
    return sum(1 for r in recs if r & node) / len(recs)


def exact_any(truth, nodes):
    """Tie-free: per level, how many truth sets equal SOME node of that level (no one-to-one)."""
    levels_of = collections.defaultdict(set)
    for (level, _), s in nodes.items():
        levels_of[s].add(level)
    c = collections.Counter()
    for s in truth.values():
        ls = levels_of.get(s, set())
        for level in ls:
            c[level] += 1
        c["any"] += bool(ls)
    return dict(sorted(c.items()))


def permutation(members, family_of, kind, rng):
    """A bijection member -> member: `global` over all members, `within_family` inside each family."""
    groups = collections.defaultdict(list)
    for m in sorted(members):
        groups["*" if kind == "global" else family_of[m]].append(m)
    out = {}
    for g in sorted(groups):
        dst = groups[g][:]
        rng.shuffle(dst)
        out.update(zip(groups[g], dst))
    return out


def relabel(sets, mapping):
    return {k: frozenset(mapping[m] for m in v) for k, v in sets.items()}


def null(scorer, truth, members, family_of, kind, n, seed):
    """Relabelling the truth by sigma^-1 scores exactly as relabelling every node by sigma, so the nodes stay fixed.
    -> (exact counts, mean Jaccards), one per round."""
    rng = random.Random(f"{seed}:{kind}")
    ex, mj = [], []
    for _ in range(n):
        ms = scorer.run(relabel(truth, permutation(members, family_of, kind, rng)))
        ex.append(sum(m.jaccard == 1.0 for m in ms))
        mj.append(sum(m.jaccard for m in ms) / len(ms) if ms else 0.0)
    return ex, mj


def p_value(values, observed):
    return (1 + sum(v >= observed for v in values)) / (1 + len(values))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("members")
    ap.add_argument("truth")
    ap.add_argument("out")
    ap.add_argument("--fold")
    ap.add_argument("--loci-gff")
    ap.add_argument("--map", choices=("exon", "key"), default="exon")
    ap.add_argument("--null", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=7)
    a = ap.parse_args(argv)
    nodes, family_of = load_nodes(a.members)
    members = sorted(family_of)
    if a.map == "exon" and not a.loci_gff:
        ap.error("--map exon needs --loci-gff (the catalog's .fam.loci.gff3)")
    sets, per_record, cov = truth_sets(load_truth(a.truth), members, load_fold(a.fold) if a.fold else None,
                                       load_loci_gff(a.loci_gff) if a.loci_gff else None, a.map)
    ev = {k: v for k, v in sets.items() if len(v) >= 2}
    hier_scorer = Scorer(nodes, members)
    fam_nodes = {k: v for k, v in nodes.items() if k[0] == "family"}
    hier, flat = hier_scorer.run(ev), Scorer(fam_nodes, members).run(ev)
    res = {"inputs": vars(a), "n_members": len(members),
           "n_nodes": dict(sorted(collections.Counter(k[0] for k in nodes).items())), "coverage": cov,
           "hierarchy": summary(hier, nodes, per_record), "flat": summary(flat, nodes, per_record),
           "exact_any": exact_any(ev, nodes), "null": {}}
    for kind in ("global", "within_family"):
        if a.null > 0:
            ex, mj = null(hier_scorer, ev, members, family_of, kind, a.null, a.seed)
            res["null"][kind] = {"rounds": a.null, "exact_mean": sum(ex) / len(ex), "exact_max": max(ex),
                                 "mean_jaccard_mean": sum(mj) / len(mj), "mean_jaccard_max": max(mj),
                                 "p_exact": p_value(ex, res["hierarchy"]["exact"]),
                                 "p_mean_jaccard": p_value(mj, res["hierarchy"]["mean_jaccard"])}
    fm = {m.label: m for m in flat}
    res["per_truth"] = [{"label": m.label, "n_records": len(per_record[m.label]), "n_members": m.t_size,
                         "level": m.level, "node": m.node_id, "jaccard": round(m.jaccard, 6),
                         "recall": round(m.recall, 6), "precision": round(m.precision, 6),
                         "gene_recall": round(gene_recall(m, nodes, per_record), 6),
                         "flat_node": fm[m.label].node_id, "flat_jaccard": round(fm[m.label].jaccard, 6)}
                        for m in hier]
    with open(a.out, "w") as fh:
        json.dump(res, fh, indent=1)
    print(json.dumps({k: v for k, v in res.items() if k not in ("per_truth", "inputs")}, indent=1))


if __name__ == "__main__":
    sys.exit(main())
```

Run `python3 -m pytest -q bench/o1_hier/test_hier_match.py`. Expected: `12 passed`.

- [ ] **Step 3: The advisor table** — test first, `bench/o1_hier/test_advisor_table.py`:

```python
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
import advisor_table as A  # noqa: E402


def test_row_reads_a_hier_match_json(tmp_path):
    d = {"coverage": {"labels": 10, "labels_evaluable": 6, "labels_one_member": 1, "labels_no_member": 3,
                      "records_on_members": 20, "records_on_unclustered_locus": 4, "records_without_locus": 6},
         "hierarchy": {"exact": 3, "exact_rate": 0.5, "mean_jaccard": 0.8, "macro_recall": 0.9, "macro_precision": 0.85,
                       "macro_gene_recall": 0.7, "by_level": {"family": 4, "expansion": 2}},
         "flat": {"exact": 2, "exact_rate": 1 / 3, "mean_jaccard": 0.7, "macro_recall": 0.75, "macro_precision": 0.9},
         "exact_any": {"expansion": 3}, "null": {"global": {"p_exact": 0.000999}, "within_family": {"p_exact": 0.02}}}
    r = A.row("s", d)
    assert (r["evaluable"], r["hier_exact"], r["hier_exact_of_all"], r["lands_expansion"], r["lands_superfamily"]) == (6, 3, 0.3, 2, 0)
    assert (r["flat_exact_rate"], r["p_global"], r["p_within_family"]) == (0.333, 0.001, 0.02)
```

(red: `ModuleNotFoundError`), then `bench/o1_hier/advisor_table.py`:

```python
#!/usr/bin/env python3
"""The advisor-facing O1 table: flat vs hierarchy-aware matching, where each truth family lands, and coverage.

    advisor_table.py OUT.tsv NAME=hier_match.json [NAME=hier_match.json ...]

One row per sample; also printed as a markdown table for the results section.
"""
import json
import sys

COLS = ["sample", "truth_families", "evaluable", "one_member", "no_member", "records_on_members",
        "records_on_unclustered_locus", "records_without_locus", "flat_exact", "flat_exact_rate", "hier_exact",
        "hier_exact_rate", "hier_exact_of_all", "flat_mean_jaccard", "hier_mean_jaccard", "flat_macro_recall",
        "hier_macro_recall", "flat_macro_precision", "hier_macro_precision", "hier_gene_recall", "lands_family",
        "lands_expansion", "lands_superfamily", "lands_sd_complex", "exact_any_expansion", "p_global", "p_within_family"]


def row(name, d):
    c, h, f = d["coverage"], d["hierarchy"], d["flat"]
    lv, nl = h["by_level"], d.get("null", {})
    r3 = lambda x: "NA" if x is None else round(x, 3)  # noqa: E731
    return {"sample": name, "truth_families": c["labels"], "evaluable": c.get("labels_evaluable", 0),
            "one_member": c.get("labels_one_member", 0), "no_member": c.get("labels_no_member", 0),
            "records_on_members": c.get("records_on_members", 0),
            "records_on_unclustered_locus": c.get("records_on_unclustered_locus", "NA"),
            "records_without_locus": c.get("records_without_locus", c.get("records_without_member", 0)),
            "flat_exact": f["exact"], "flat_exact_rate": r3(f["exact_rate"]), "hier_exact": h["exact"],
            "hier_exact_rate": r3(h["exact_rate"]), "hier_exact_of_all": r3(h["exact"] / c["labels"]),
            "flat_mean_jaccard": r3(f["mean_jaccard"]), "hier_mean_jaccard": r3(h["mean_jaccard"]),
            "flat_macro_recall": r3(f["macro_recall"]), "hier_macro_recall": r3(h["macro_recall"]),
            "flat_macro_precision": r3(f["macro_precision"]), "hier_macro_precision": r3(h["macro_precision"]),
            "hier_gene_recall": r3(h.get("macro_gene_recall")),
            "lands_family": lv.get("family", 0), "lands_expansion": lv.get("expansion", 0),
            "lands_superfamily": lv.get("superfamily", 0), "lands_sd_complex": lv.get("sd_complex", 0),
            "exact_any_expansion": d["exact_any"].get("expansion", 0),
            "p_global": r3(nl.get("global", {}).get("p_exact")),
            "p_within_family": r3(nl.get("within_family", {}).get("p_exact"))}


def main(argv=None):
    a = argv if argv is not None else sys.argv[1:]
    rows = [row(n, json.load(open(p))) for n, p in (x.split("=", 1) for x in a[1:])]
    with open(a[0], "w") as fh:
        fh.write("\t".join(COLS) + "\n")
        fh.writelines("\t".join(str(r[k]) for k in COLS) + "\n" for r in rows)
    print("| " + " | ".join(COLS) + " |")
    print("|" + "---|" * len(COLS))
    for r in rows:
        print("| " + " | ".join(str(r[k]) for k in COLS) + " |")


if __name__ == "__main__":
    sys.exit(main())
```

Run `python3 -m pytest -q bench/o1_hier`. Expected: `29 passed` (5 + 7 + 4 + 12 + 1).

- [ ] **Step 4: Regression on the 10-09 F_O run** (v1.0, guided universe; the members ARE the annotation records)

```bash
H=/mnt/linuxdisk/tmp/hier_flag; L=/mnt/linuxdisk/tmp/o1_close/lift; O=/mnt/linuxdisk/tmp/o1_close/fo_check; mkdir -p $O
RLOCK_TIMEOUT=300 bash tools/rlock.sh light python3 bench/o1_hier/hier_match.py $H/final/ours.hierarchy.members.tsv $H/soto_external.tsv $O/key.json --fold $H/final/ours.loci.tsv --map key --null 200 > /dev/null
RLOCK_TIMEOUT=300 bash tools/rlock.sh light python3 bench/o1_hier/hier_match.py $H/final/ours.hierarchy.members.tsv $L/soto_truth.v1.tsv $O/exon.json --fold $H/final/ours.loci.tsv --loci-gff /mnt/linuxdisk/tmp/ourhom/loci.gff3 --null 200 > /dev/null
python3 bench/o1_hier/advisor_table.py $O/advisor.tsv FO_key=$O/key.json FO_exon=$O/exon.json
```

Expected exactly (measured while writing the plan; any difference is a transcription error):
- **`key.json`** (this reproduces the Rust placement of 10-09):
  - coverage `labels` 491, `labels_evaluable` 400, `labels_one_member` 12, `labels_no_member` 79;
  - `exact_any.expansion` 215; `hierarchy.exact` 245; `flat.exact` 228; `hierarchy.mean_jaccard` 0.8245 (4 dp);
  - `null.global.exact_mean` 0.045; `null.within_family.exact_mean` 236.93.
- **`exon.json`:**
  - `labels_evaluable` 402; `exact_any.expansion` 206; `hierarchy.exact` 231; `flat.exact` 217;
  - `hierarchy.mean_jaccard` 0.8056; `null.within_family.exact_mean` 223.9.

- [ ] **Step 5: Commit**

```bash
git add bench/o1_hier/hier_match.py bench/o1_hier/test_hier_match.py bench/o1_hier/advisor_table.py bench/o1_hier/test_advisor_table.py
git commit -m "o1_hier: hierarchy-aware Hungarian scorer (exon+strand truth mapping, named tie policy, global + within-family nulls, flat continuity); reproduces the 10-09 placement (215)"
```

---

### Task 5: pre-register the human evaluation (E1, E1b, E2)

**Files:**
- Create: `docs/archive/2026-10/merged/PREREG_o1_hierarchy_eval_2026-10-13.md`

- [ ] **Step 1: Write the prereg** with exactly:

```markdown
# PREREG — O1 hierarchy-aware evaluation, human (2026-10-13)

Spec: `docs/superpowers/specs/2026-10-13-close-o1-design.md` section A; plan
`docs/superpowers/plans/2026-10-13-close-o1-phase-0A.md`.

**Question.** Are Soto 2025's families nodes of our family hierarchy? The flat one-to-one match against our MCL
families counts every Soto family that sits one level down as a miss.

**Catalogs.**
- Dev: human A119b. Held-out: human testis.
- Both are default catalogs, re-clustered with the MCL two-node fix (`docs/archive/2026-10/merged/PREREG_mcl_two_node_2026-10-13.md`).

**Hierarchy.**
- Command: `mcl_families.after --paf <S>.fam.loci.paf --gff <S>.fam.loci.gff3 --min-exonic-bp 1
  --min-shared-exon-frac 0.60 --emit-hierarchy --hierarchy-duplicons lift/dupmasker_colors.v2.bed --hierarchy-sd
  lift/sedef_pairs.v2.bed`.
- Inputs: DupMasker colors and SEDEF `final_v1_clean.bed`, lifted from CHM13 v1.0 by the T2T chain (sequence-verified).
- No `--hierarchy-external`: the truth never enters the build.
- Its families must equal the re-clustered catalog (`cmp`).

**Truth.** `lift/soto_truth.v2.tsv`:
- S1C families proper (491; `Unassigned_*` excluded), genes joined to CAT v4 by Gene ID;
- exons = the union over the gene's transcripts, chain-lifted to v2.0.

**Mapping.**
- A Soto gene joins a member when one of its exons shares >= 1 bp, on the same strand, with an exon of the member's
  locus or of a locus folded into it.
- Exons come from the catalog's `.fam.loci.gff3`; the fold map is the run's `ours.loci.tsv`.
- A family is EVALUABLE with >= 2 members. The others are reported as coverage.

**Scorer.** `bench/o1_hier/hier_match.py --map exon --null 1000 --seed 7`.
- One-to-one Hungarian (scipy) on -Jaccard between evaluable truth families and every family, expansion,
  superfamily and SD-complex node.
- Tie policy: rows in label order; columns in (level rank family < expansion < superfamily < sd_complex, node id)
  order.
  As built (8f0515ef), canonical placement follows: see the prereg's Tie policy (spec amendment 5).
- Nulls, both keeping node sizes and nesting:
  - **global**: one random permutation of all members per round;
  - **within_family**: members permuted inside their own family. Families, superfamilies and SD complexes are
    unchanged; only expansions move.

**Bars** (`exact` = Hungarian matches with Jaccard 1; `rate` = exact / evaluable):
- **E1:** on each catalog, `null.global.p_exact` < 0.001 AND `null.global.p_mean_jaccard` < 0.001.
- **E1b:** on each catalog, `null.within_family.p_exact` < 0.001 AND `null.within_family.p_mean_jaccard` < 0.001.
- **E2:** |rate(testis) - rate(A119b)| <= 0.10.

**Reading, fixed now:**
- **E1 + E1b + E2 pass:** "Soto's families are nodes of our hierarchy". The O1 headline for Soto becomes, per catalog,
  the exact rate and mean Jaccard, next to the flat score and the coverage.
- **E1 + E2 pass, E1b fails:** the claim stops at the family level and above. Our expansions recover Soto no better
  than random within-family splits of the same sizes.
- **E2 fails:** sample-specific; only the per-catalog numbers are reported.
- **E1 fails:** no claim.

**Reported regardless:**
- coverage (records on members / on unclustered loci / without a locus; labels evaluable / one member / none; members
  in several truth families);
- the hierarchy and flat summaries (exact, rate, mean Jaccard, macro and micro recall / precision, gene recall);
- `exact_any` per level, and where matches land (`by_level`);
- the unconditioned exact / 491;
- the per-truth table;
- the advisor table (`bench/o1_hier/advisor_table.py`).

**Disclosures.**
- E1b, the exon + strand mapping and the >= 2-member evaluable rule were fixed after the scorer's regression run on the
  10-09 F_O universe (v1.0, guided). There the global null averaged 0.03 exact matches; the within-family null averaged
  224 of the observed 231; span-overlap mapping cut exact expansions from 215 to 169.
- No de novo catalog has been scored by this scorer before this prereg.
- A non-Soto member inside a matched node counts against precision. Soto's truth covers only its SD98 gene set, so
  precision is a lower bound.

## Results

(pending)
```

```bash
git add docs/archive/2026-10/merged/PREREG_o1_hierarchy_eval_2026-10-13.md
git commit -m "PREREG O1 hierarchy-aware evaluation, human (E1 global null, E1b within-family null, E2 dev/held-out)"
```

---

### Task 6: build the human hierarchies, score them, record

**Files:**
- Create: `docs/hierarchy_results/o1_close_human_2026-10-13.tsv` (the advisor table)
- Modify: `docs/archive/2026-10/merged/PREREG_o1_hierarchy_eval_2026-10-13.md` (Results only)
- Modify: `docs/NEGATIVE_RESULTS_REGISTER.md` (append a row)
- Modify: `docs/THESIS_OBJECTIVES.md` (the O1 row of the STATUS 2026-10-11 table, as the reading allows)

**Interfaces:**
- Consumes:
  - `mcl_families.after` and `recluster/after/*` (Tasks 1-2);
  - `lift/*.v2.*` (Task 3);
  - `hier_match.py` and `advisor_table.py` (Task 4);
  - the committed prereg (Task 5).

- [ ] **Step 1: Check that the prereg is committed, then build both hierarchies**

```bash
git log --oneline -1 -- docs/archive/2026-10/merged/PREREG_o1_hierarchy_eval_2026-10-13.md   # must print the Task 5 commit
L=/mnt/linuxdisk/tmp/o1_close/lift; B=/mnt/linuxdisk/tmp/o1_close/mcl_families.after
for S in human_testis human_A119b; do
  R=/mnt/linuxdisk/tmp/rustle_figures/runs/$S; O=/mnt/linuxdisk/tmp/o1_close/hier/$S; mkdir -p $O
  RLOCK_TIMEOUT=590 bash tools/rlock.sh heavy $B --paf $R/$S.fam.loci.paf --gff $R/$S.fam.loci.gff3 --min-exonic-bp 1 --min-shared-exon-frac 0.60 \
    --emit-hierarchy --hierarchy-duplicons $L/dupmasker_colors.v2.bed --hierarchy-sd $L/sedef_pairs.v2.bed --out $O/ours > $O/log 2>&1
  grep "emit-hierarchy" $O/log
  cmp $O/ours.clusters.tsv /mnt/linuxdisk/tmp/o1_close/recluster/after/$S.fam.clusters.tsv && echo "$S: families = the re-clustered catalog"
done
```

Expected: both `cmp` lines print. If one differs, STOP.

- [ ] **Step 2: Score** (heavy lock: 1,000 rounds x 2 nulls on ~10k members)

```bash
L=/mnt/linuxdisk/tmp/o1_close/lift
for S in human_testis human_A119b; do
  R=/mnt/linuxdisk/tmp/rustle_figures/runs/$S; O=/mnt/linuxdisk/tmp/o1_close/hier/$S
  RLOCK_TIMEOUT=1800 bash tools/rlock.sh heavy python3 bench/o1_hier/hier_match.py $O/ours.hierarchy.members.tsv $L/soto_truth.v2.tsv $O/hier_match.json \
    --loci-gff $R/$S.fam.loci.gff3 --fold $O/ours.loci.tsv --null 1000 --seed 7 > $O/hier_match.log 2>&1; echo "$S exit $?"
done
python3 bench/o1_hier/advisor_table.py docs/hierarchy_results/o1_close_human_2026-10-13.tsv A119b=/mnt/linuxdisk/tmp/o1_close/hier/human_A119b/hier_match.json testis=/mnt/linuxdisk/tmp/o1_close/hier/human_testis/hier_match.json
```

- If a scoring call would exceed the 600000 ms tool timeout, run it with `run_in_background` and
  `RLOCK_TIMEOUT=3000`.
- Then compute the verdicts from the two JSONs exactly as the prereg defines them: E1, E1b, E2 with
  rate = `hierarchy.exact_rate`.

- [ ] **Step 3: Record**

1. **Prereg Results** (replace `(pending)`):
   - the three verdicts with their numbers (p-values, rates and their difference);
   - the advisor table's markdown, as printed;
   - per catalog: `exact_any`, `by_level`, the null summaries, and coverage;
   - the reading that applies, quoted from the prereg.
2. **Register row** (next free number; 1277 if Task 2 used 1276):
   `| 1277 | 2026-10-13 | O1 hierarchy-aware evaluation vs Soto (docs/archive/2026-10/merged/PREREG_o1_hierarchy_eval_2026-10-13.md) | Soto's families are nodes of our hierarchy (dev A119b, held-out testis) | <verdicts + headline numbers> |`.
3. **THESIS_OBJECTIVES:** edit the O1 row only as the reading allows.
   - All three pass: add the hierarchy claim with both catalogs' exact rates and mean Jaccard next to the flat ones.
   - E1b fails: add the family-level claim only.
   - Otherwise: add one clause, "hierarchy-aware Soto evaluation: not claimed (<which bar failed>)".

```bash
git add docs/archive/2026-10/merged/PREREG_o1_hierarchy_eval_2026-10-13.md docs/NEGATIVE_RESULTS_REGISTER.md docs/THESIS_OBJECTIVES.md docs/hierarchy_results/o1_close_human_2026-10-13.tsv
git commit -m "O1 hierarchy-aware evaluation, human: E1/E1b/E2 results + advisor table"
```

---

## After this plan

- Phase V (`sedef-rs` vendored as `src/sedef/`, V1/V2) and phase B (gorilla, chimp vs Liftoff-Soto + Yoo 2025) get
  the next plan, written from this plan's results.
- V1's CHM13 SEDEF reference is `lift/sedef_pairs.v2.bed` from Task 3, not `HSA_sedef_pairs.bed`.
- The whole-branch review runs at the end. Ask before pushing.
